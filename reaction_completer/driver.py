import logging
import re
from collections import defaultdict
from functools import reduce
from operator import or_
from tokenize import TokenError
from nltk.metrics.distance import edit_distance
from reaction_completer import ReactionCompleter
from reaction_completer.errors import (
    TooFewPrecursors, TooManyPrecursors,
    CannotBalance, FormulaException)
from reaction_completer.formatting import render_reaction
from reaction_completer.material import MaterialInformation

__author__ = 'Haoyan Huo'
__maintainer__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'


def substitute_element_vars(targets):
    # Generate all possible combinations of (target, element_substitution)
    targets_to_balance = []
    for target in targets:
        has_element_vars = len(target['elements_vars']) != 0
        target_elements = reduce(or_, [set(x['elements'].keys()) for x in target['composition']])
        element_vars_used = set(target['elements_vars'].keys()) & target_elements
        if has_element_vars and len(element_vars_used) > 0:
            for sub in element_vars_used:
                for subs in target['elements_vars'][sub]:
                    targets_to_balance.append((target, {sub: subs}))
        else:
            targets_to_balance.append((target, None))

    return targets_to_balance


def material_dict_to_info(material_dict, sub_dict=None):
    compositions = []
    for comp in material_dict['composition']:
        compositions.append({
            'amount': comp['amount'],
            'elements': comp['elements'],
        })
    return MaterialInformation(
        material_dict['material_string'],
        material_dict['material_formula'],
        compositions, sub_dict)


def screen_good_precursors(precursors):
    precursor_objects = []
    for precursor in precursors:
        try:
            precursor_objects.append(material_dict_to_info(precursor))
        except FormulaException as e:
            logging.debug('Dropping precursors %s because %r', precursor['material_string'], e)
            continue

    return precursor_objects


def try_balance(precursors_to_balance, target, substitution, all_precursors):
    target_to_balance = material_dict_to_info(target, substitution)
    completer = ReactionCompleter(precursors_to_balance, target_to_balance)
    solution = completer.compute_reactions()

    return (
        target_to_balance.material_formula,
        solution,
        substitution,
        render_reaction(all_precursors, target, solution, substitution)
    )


def find_precursors_in_same_sentence(precursors_to_balance, sentences):
    """
    Try to eliminate the precursors not in the same sentence.
    Also finds sets that don't come from a material name (such as manganese nitrates).
    """
    precursor_candidates = []

    def is_word_material(m):
        return bool(re.match(r'^[\w\s()]+$', m.material_string))

    for sentence in sentences:
        candidates = []
        for precursor in precursors_to_balance:
            if precursor.material_formula in sentence:
                candidates.append(precursor)
        if candidates:
            precursor_candidates.append(candidates)

        candidates_no_conversion = [
            x for x in candidates
            if edit_distance(x.material_formula, x.material_string) < len(x.material_string) * 0.5]
        if candidates_no_conversion:
            precursor_candidates.append(candidates_no_conversion)

    # Find the list of precursors that are in the same sentence
    for sentence in sentences:
        candidates = []
        for precursor in precursors_to_balance:
            if precursor.material_string in sentence:
                candidates.append(precursor)

        candidates_no_conversion = [
            x for x in candidates
            if edit_distance(x.material_formula, x.material_string) < len(x.material_string) * 0.5]
        if candidates_no_conversion:
            precursor_candidates.append(candidates_no_conversion)

        if candidates:
            precursor_candidates.append(candidates)

            # Make a copy of materials that don't come from English words
            # if a similar material has been found.
            materials_by_chemistry = defaultdict(list)
            for p in candidates:
                materials_by_chemistry[tuple(sorted(p.all_elements))].append(p)

            for chemistry in set(materials_by_chemistry):
                materials = materials_by_chemistry[chemistry]
                if len(materials) > 1:
                    materials_by_chemistry[chemistry] = [x for x in materials if not is_word_material(x)]
                    if not materials_by_chemistry[chemistry]:
                        materials_by_chemistry[chemistry] = materials

            candidates_no_words = []
            for i in materials_by_chemistry.values():
                candidates_no_words.extend(i)
            precursor_candidates.append(candidates_no_words)

    return precursor_candidates


def balance_recipe(precursors, targets, sentences=None):
    """
    Balance a recipe extracted using synthesis project pipeline.

    If argument "sentences" is a list of sentences in the paragraph, when
    too many materials are given and we are unable to determine the right
    set of precursors, this function will try to gather precursors in the
    same sentence, and use the subset to balance the reaction.

    Example usage:
    >>>
    >>> precursors = [
    >>>     {
    >>>         "material_formula": "SrCO3",
    >>>         "material_string": "SrCO3",
    >>>         "composition": [
    >>>             {
    >>>                 "formula": "SrCO3",
    >>>                 "elements": {"Sr": "1.0", "C": "1.0", "O": "3.0"},
    >>>                 "amount": "1.0"
    >>>             }
    >>>         ],
    >>>     },
    >>>     {
    >>>         "material_formula": "Al2O3",
    >>>         "material_string": "Al2O3",
    >>>         "composition": [
    >>>             {
    >>>                 "formula": "Al2O3",
    >>>                 "elements": {"Al": "2.0", "O": "3.0"},
    >>>                 "amount": "1.0"
    >>>             }
    >>>         ],
    >>>     },
    >>>     {
    >>>         "material_formula": "MnO",
    >>>         "material_string": "MnO",
    >>>         "composition": [
    >>>             {
    >>>                 "formula": "MnO",
    >>>                 "elements": {"Mn": "1.0", "O": "1.0"},
    >>>                 "amount": "1.0"
    >>>             }
    >>>         ],
    >>>     },
    >>>     {
    >>>         "material_formula": "Fe2O3",
    >>>         "material_string": "Fe2O3",
    >>>         "composition": [
    >>>             {
    >>>                 "formula": "Fe2O3",
    >>>                 "elements": {"Fe": "2.0", "O": "3.0"},
    >>>                 "amount": "1.0"
    >>>             }
    >>>         ],
    >>>     },
    >>>     {
    >>>         "material_formula": "ZrO2",
    >>>         "material_string": "ZrO2",
    >>>         "composition": [
    >>>             {
    >>>                 "formula": "ZrO2",
    >>>                 "elements": {"Zr": "1.0", "O": "2.0"},
    >>>                 "amount": "1.0"
    >>>             }
    >>>         ]
    >>>     },
    >>>     {
    >>>         "material_formula": "H2O",
    >>>         "material_string": "H2O",
    >>>         "composition": [
    >>>             {
    >>>                 "formula": "H2O",
    >>>                 "elements": {"O": "1.0", "H": "2.0"},
    >>>                 "amount": "1.0"
    >>>             }
    >>>         ]
    >>>     },
    >>> ]
    >>> targets = [
    >>>     {
    >>>         "material_formula": "Sr6(A2O4)6",
    >>>         "material_string": "Sr6(A2O4)6",
    >>>         "composition": [
    >>>             {
    >>>                 "formula": "Sr6(Fe2O4)6",
    >>>                 "elements": {"A": "12.0", "O": "24.0", "Sr": "6.0"},
    >>>                 "amount": "1.0"
    >>>             }
    >>>         ],
    >>>         "elements_vars": {
    >>>             "A": ["Fe", "Al"]
    >>>         },
    >>>         "additives": ["Mn2+"]
    >>>     },
    >>> ]
    >>> text = [
    >>>     "SrCO3, Al2O3, MnO and Fe2O3 are used to synthesize Mn2+doped-Sr6(A2O4)6, A=Fe, Al.",
    >>>     "Milling media is ZrO2",
    >>>     "There is some H2O found in the final product."
    >>> ]
    >>>
    >>> reactions = balance_recipe(precursors, targets, text)
    >>> print('Found', len(reactions), 'reactions')
    >>> for reaction in reactions:
    >>>     print(reaction)

    :param precursors: List of precursors
    :param targets: List of targets
    :param sentences: List of sentences
    :return:
    """
    sentences = sentences or []

    targets_to_balance = substitute_element_vars(targets)
    precursors_to_balance = screen_good_precursors(precursors)

    solutions = []
    for target, substitution in targets_to_balance:
        try:
            # Skip if substitution is bad
            target_object = material_dict_to_info(target, substitution)
        except FormulaException as e:
            logging.debug('Failed to convert target! %r', e)
            continue

        try:
            solution = try_balance(precursors_to_balance, target, substitution, precursors)
            solutions.append(solution)
        except TooFewPrecursors:
            precursor_candidates = list(filter(lambda x: not x.is_hco, precursors_to_balance))
            try:
                solution = try_balance(precursor_candidates, target, substitution, precursors)
                solutions.append(solution)
            except (CannotBalance, TokenError) as e_subset:
                logging.debug(
                    'Failed trying inorganic precursor subset for target: %s, '
                    'precursors: %r: %r',
                    target_object.material_formula,
                    [x.material_formula for x in precursor_candidates], e_subset)
        except TooManyPrecursors as e:
            precursor_candidates = find_precursors_in_same_sentence(precursors_to_balance, sentences)

            if not precursor_candidates:
                logging.debug('No possible precursor subsets for target: %s, precursors: %r: %r',
                              target_object.material_formula,
                              [x.material_formula for x in precursors_to_balance], e)
            else:
                success = False
                # Iterate over all candidate precursors, and find the first success
                # reaction that can be completed.
                for candidates in precursor_candidates:
                    try:
                        solution = try_balance(candidates, target, substitution, precursors)
                        solutions.append(solution)
                        success = True
                        break
                    except (CannotBalance, TokenError) as e_subset:
                        logging.debug(
                            'Failed trying precursor subset for target: %s, '
                            'precursors: %r: %r',
                            target_object.material_formula,
                            [x.material_formula for x in candidates], e_subset)

                if not success:
                    logging.debug('Cannot find a subset of precursors for '
                                  'target: %s, precursors: %r: %r',
                                  target_object.material_formula,
                                  [x.material_formula for x in precursors_to_balance], e)

        except (CannotBalance, TokenError) as e:
            logging.debug('Cannot balance reaction for target: %s, precursors: %r: %r',
                          target_object.material_formula,
                          [x.material_formula for x in precursors_to_balance], e)
        except Exception as e:
            logging.warning('Unexpected error for target: %s, precursors: %r: %r',
                            target_object.material_formula,
                            [x.material_formula for x in precursors_to_balance], e)
    return solutions
