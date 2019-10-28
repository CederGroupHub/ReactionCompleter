import logging
from functools import reduce
from operator import or_

import sympy
from sympy import Matrix, symbols

from reaction_completer.errors import (
    StupidRecipe, TooManyPrecursors, TooFewPrecursors)
from reaction_completer.formatting import simplify_print
from reaction_completer.material import MaterialInformation

__author__ = 'Haoyan Huo'
__maintainer__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'


class ReactionCompleter(object):
    def __init__(self, precursors: [MaterialInformation],
                 target: MaterialInformation,
                 target_min_nv=2):
        """
        A reaction completer that takes a set of precursors and a target,
        then calculates the possible reactions, using sympy for symbolic
        derivation.

        :param precursors: List of precursors.
        :type precursors: list(MaterialInformation)
        :param target: The target material.
        :type target: MaterialInformation
        :param target_min_nv:
        """
        self.precursors = precursors
        self.target = target
        self.target_min_nv = target_min_nv

        self._precursor_candidates = []
        self._decomposition_chemicals = {}
        self._exchange_chemicals = {}
        self._linear_eq = {}

        self._inspect_target()
        self._prepare_precursors()
        self._setup_linear_equation()

    def _inspect_target(self):
        """
        Prepare the target material into a ready-to-use structure.
        """
        if len(self.target.nv_elements) < self.target_min_nv:
            raise StupidRecipe(
                'Target must have more than 1 non volatile elements, got %r: %s' %
                (self.target.nv_elements, self.target.material_formula))

    def _prepare_precursors(self):
        # find the set of precursors

        seen_precursors = set()
        for precursor in self.precursors:
            # Skip precursors that are seen
            if precursor.material_formula in seen_precursors:
                continue

            seen_precursors.add(precursor.material_formula)

            if precursor.all_elements_dict == self.target.all_elements_dict:
                # TODO: we need a smarter comparison
                raise StupidRecipe('Precursor list contains target')

            if len(precursor.all_elements) == 0:
                logging.debug(
                    'Skipping empty precursor %s: %s',
                    precursor.material_formula)
                continue

            if len(precursor.nv_elements - self.target.nv_elements) > 0:
                logging.debug(
                    'Skipping precursor %s because it '
                    'has excessive chemical elements',
                    precursor.material_formula)
                continue

            self._precursor_candidates.append(precursor)
            self._decomposition_chemicals.update(precursor.decompose_chemicals)

        self._exchange_chemicals.update(self.target.exchange_chemicals)

        if len(self._precursor_candidates) == 0:
            raise StupidRecipe('Precursor candidates is empty')

        # check for equality
        precursors_nv_elements = reduce(
            or_, [x.nv_elements for x in self._precursor_candidates])
        missing_elements = self.target.nv_elements - precursors_nv_elements
        if len(missing_elements) > 0:
            raise StupidRecipe(
                'Precursor candidates do not '
                'provide non volatile elements: %r' % missing_elements)

    def _setup_linear_equation(self):
        all_elements = reduce(
            or_,
            [x.all_elements for x in self._precursor_candidates] +
            [self.target.all_elements] +
            [set(x) for x in self._exchange_chemicals.values()] +
            [set(x) for x in self._decomposition_chemicals.values()])
        all_elements = sorted(list(all_elements))

        # Create the symbols that will be used for linear eq.
        chemical_symbols = ''
        for i in range(len(self._precursor_candidates)):
            chemical_symbols += 'p%d, ' % i
        for i in range(len(self._decomposition_chemicals)):
            chemical_symbols += 'r%d, ' % i
        for i in range(len(self._exchange_chemicals)):
            chemical_symbols += 'e%d, ' % i
        chemical_symbols += 't'
        chemical_symbols = symbols(chemical_symbols)

        coefficient_matrix = []
        which_side = []

        def fill_row(material_elements):
            row = [material_elements.get(element, 0) for element in all_elements]
            coefficient_matrix.append(row)

        for precursor in self._precursor_candidates:
            fill_row(precursor.all_elements_dict)
            which_side.append('fl')

        for chemical in sorted(self._decomposition_chemicals):
            fill_row(self._decomposition_chemicals[chemical])
            which_side.append('dr')

        for chemical in sorted(self._exchange_chemicals):
            fill_row(self._exchange_chemicals[chemical])
            which_side.append('dl')

        target_elements = self.target.all_elements_dict
        target_vector = [target_elements.get(i, 0) for i in all_elements]

        coefficient_matrix = Matrix(coefficient_matrix).T

        target_vector = Matrix(target_vector)

        self._linear_eq.update({
            'chemical_symbols': chemical_symbols,
            'coefficient_matrix': coefficient_matrix,
            'target_vector': target_vector,
            'which_side': which_side,
            'all_elements': all_elements,
        })

    def _render_reaction(self, solution: tuple):
        balanced = {
            'left': {},
            'right': {self.target.material_formula: '1'}
        }

        solution = list(solution)
        which_side = self._linear_eq['which_side']

        precursor_solutions = solution[:len(self._precursor_candidates)]
        precursor_side = which_side[:len(self._precursor_candidates)]
        del solution[:len(self._precursor_candidates)]
        del which_side[:len(self._precursor_candidates)]

        decomposition_solutions = solution[:len(self._decomposition_chemicals)]
        decomposition_side = which_side[:len(self._decomposition_chemicals)]
        del solution[:len(self._decomposition_chemicals)]
        del which_side[:len(self._decomposition_chemicals)]

        exchange_solutions = solution.copy()
        exchange_side = which_side.copy()

        def decide_side_value(s, val):
            if s[0] == 'f':
                if s[1] == 'l':
                    return 'left', val
                elif s[1] == 'r':
                    return 'right', -val
            elif s[0] == 'd':
                if not isinstance(val, sympy.Float):
                    value_zero = val.evalf(
                        subs={x: 0.001 for x in val.free_symbols})
                    value_negative = float(value_zero) < 0
                else:
                    value_negative = float(val) < 0

                if s[1] == 'l':
                    return ('left', val) if not value_negative else ('right', -val)
                elif s[1] == 'r':
                    return ('right', -val) if value_negative else ('left', val)

        for precursor, amount, side in zip(
                self._precursor_candidates, precursor_solutions, precursor_side):
            material_formula = precursor.material_formula
            side, value = decide_side_value(side, amount)

            value_s = simplify_print(value)
            if value_s != '0':
                balanced[side][material_formula] = value_s

        for chemical, amount, side in zip(
                sorted(self._decomposition_chemicals), decomposition_solutions, decomposition_side):
            side, value = decide_side_value(side, amount)

            value_s = simplify_print(value)
            if value_s != '0':
                balanced[side][chemical] = value_s

        for chemical, amount, side in zip(
                sorted(self._exchange_chemicals), exchange_solutions, exchange_side):
            side, value = decide_side_value(side, amount)

            value_s = simplify_print(value)
            if value_s != '0':
                balanced[side][chemical] = value_s

        return balanced

    def compute_reactions(self):
        try:
            a = self._linear_eq['coefficient_matrix']
            b = self._linear_eq['target_vector']
            solution, params = a.gauss_jordan_solve(b)

            if len(params) > 0:
                raise TooManyPrecursors(
                    'Too many precursors to balance %r ==> %r' % (
                        [x.material_formula for x in self.precursors],
                        self.target.material_formula))

            solution = solution.T[:1, :]
        except ValueError:
            raise TooFewPrecursors('Too few precursors to balance')

        return self._render_reaction(solution)
