from tokenize import TokenError

from sympy.parsing.sympy_parser import parse_expr

from reaction_completer.errors import FormulaException
from reaction_completer.periodic_table import NON_VOLATILE_ELEMENTS, ELEMENTS

__author__ = 'Haoyan Huo'
__maintainer__ = 'Haoyan Huo'
__email__ = 'haoyan.huo@lbl.gov'


class MaterialInformation(object):
    def __init__(self, material_string, material_formula,
                 material_composition, substitution_dict=None):
        """
        Represents information about a material.

        material_string is the original text excerpt from the synthesis
        paragraph. It must match exactly with the original sequence of
        characters. It is used in some heuristics to find list of
        precursors mentioned in same sentences.

        material_formula is a string representation of the material.
        This is usually the human-readable representation, and often
        includes some conventions. For example, Fe2O3⋅H2O. This is
        useful for some heuristics discovery, such as the detection of
        H2O release.

        material_composition is either one of the following:

        1. A dictionary with keys 'amount' and 'elements'. The 'amount'
        suggests the molar amount of this formula. The 'elements' is a
        dictionary, whose keys are elements, values are amounts.
        2. A list of dictionaries, the dictionaries follow the same
            structure as above.

        For example, material composition can take the form of:
            [
                {'amount': '1.0', 'elements': {'O': 1, 'H': 2}},
                ...
            ]

        substitution_dict is either None, indicating no substitution is
        to be made; or a dictionary containing the substitution of
        elements in the material_composition dictionary.
        """

        # Type checks
        if not isinstance(material_string, str):
            raise TypeError('material_string must be a string, got %r' % type(material_string))
        if not isinstance(material_formula, str):
            raise TypeError('material_formula must be a string, got %r' % type(material_formula))

        if not isinstance(material_composition, (list, tuple)):
            material_composition = [material_composition]

        # Ensure the composition has right data types
        all_elements = set()
        for composition in material_composition:
            if set(composition.keys()) != {'amount', 'elements'}:
                raise ValueError('Illegal composition dictionary %r. '
                                 'You should only put keys "amount" '
                                 'and "elements"' % set(composition.keys()))
            if not isinstance(composition['amount'], str):
                composition['amount'] = str(composition['amount'])

            for element, amount in composition['elements'].items():
                all_elements.add(element)
                if not isinstance(element, str):
                    raise TypeError('composition.elements keys must be str, got %r' % type(element))
                if not isinstance(amount, str):
                    composition['elements'][element] = str(amount)

        substituted_elements = set()
        if substitution_dict is not None:
            if not isinstance(substitution_dict, dict):
                raise TypeError('substitution_dict must be a dict, got %r' % type(substitution_dict))
            for from_element, to_element in substitution_dict.items():
                if not isinstance(from_element, str):
                    raise TypeError('substitution_dict keys must be strings, got %r' % type(from_element))
                if not isinstance(to_element, str):
                    raise TypeError('substitution_dict values must be strings, got %r' % type(to_element))
                if from_element not in all_elements:
                    raise ValueError('substitution_dict keys must be in composition.elements, got %s' % from_element)
                substituted_elements.add(from_element)
                if to_element not in ELEMENTS:
                    raise ValueError('substitution_dict values must be chemical elements, got %s' % to_element)
        for element in all_elements - substituted_elements:
            if element not in ELEMENTS:
                raise FormulaException(
                    'composition.elements has non-chemical, non-substitutional element: %s' % element)

        self.material_string = material_string
        self.material_formula = material_formula
        self.material_composition = material_composition
        self.substitution_dict = substitution_dict or {}

        self.non_volatile_elements = {}
        self.other_elements = {}

        self._parse()

    def _parse(self):
        for component in self.material_composition:
            try:
                fraction = parse_expr(component['amount'])
            except (SyntaxError, TokenError):
                raise FormulaException(
                    'Sympy cannot parse component molar fraction: %s'
                    % component['amount'])

            for element, amount_s in component['elements'].items():
                element = self.substitution_dict.get(element, element)

                try:
                    amount = parse_expr(amount_s)
                except (SyntaxError, TokenError):
                    raise FormulaException(
                        'Sympy cannot parse element amount: %s'
                        % amount_s)

                if element not in ELEMENTS:
                    raise FormulaException(
                        '%s is not a valid chemical element' % element)

                if element in NON_VOLATILE_ELEMENTS:
                    if element not in self.non_volatile_elements:
                        self.non_volatile_elements[element] = fraction * amount
                    else:
                        self.non_volatile_elements[element] += fraction * amount
                else:
                    if element not in self.other_elements:
                        self.other_elements[element] = fraction * amount
                    else:
                        self.other_elements[element] += fraction * amount

    def __str__(self):
        return '<MaterialInformation for %s>' % self.material_formula

    def __repr__(self):
        return self.__str__()

    def __unicode__(self):
        return self.__str__()

    @property
    def nv_elements_dict(self):
        return self.non_volatile_elements

    @property
    def nv_elements(self):
        return set(self.non_volatile_elements.keys())

    @property
    def v_elements_dict(self):
        return self.other_elements

    @property
    def v_elements(self):
        return set(self.other_elements.keys())

    @property
    def all_elements_dict(self):
        a = self.non_volatile_elements.copy()
        a.update(self.other_elements)
        return a

    @property
    def all_elements(self):
        return self.nv_elements | self.v_elements

    @staticmethod
    def _compare_composition(comp1, comp2):
        if comp1.keys() != comp2.keys() or len(comp1) == 0:
            return False
        scaling = set()
        for k in comp1:
            try:
                value1, value2 = float(comp1[k]), float(comp2[k])
            except ValueError:
                return False
            scaling.add(value1 / value2)
        return len(scaling) == 1

    def _component_search_helper(self, target_comp, exclude_nv=False):
        for comp in self.material_composition:
            if exclude_nv:
                subset = {x: y for x, y in comp['elements'].items() if x not in NON_VOLATILE_ELEMENTS}
            else:
                subset = comp['elements']
            if self._compare_composition(subset, target_comp):
                return True
        return False

    @property
    def is_hco(self):
        return self.all_elements == {'C', 'H', 'O'}

    _COMP_WATER = {'H': 2, 'O': 1}

    @property
    def has_water(self):
        return self._component_search_helper(self._COMP_WATER)

    _COMP_ACETATE = {'C': 2, 'H': 3, 'O': 2}
    _COMP_ACETATE_CHARGED = {'C': 2, 'H': 3, 'O': 2, 'e-': 1}

    @property
    def has_acetate(self):
        return self._component_search_helper(
            self._COMP_ACETATE, exclude_nv=True)

    _COMP_NITRATE = {'N': 1, 'O': 3}
    _COMP_NITRATE_CHARGED = {'N': 1, 'O': 3, 'e-': 1}

    @property
    def has_nitrate(self):
        return self._component_search_helper(
            self._COMP_NITRATE, exclude_nv=True)

    _COMP_HYDROXIDE = {'H': 1, 'O': 1}
    _COMP_HYDROXIDE_CHARGED = {'H': 1, 'O': 1, 'e-': 1}

    @property
    def has_hydroxide(self):
        return self._component_search_helper(
            self._COMP_HYDROXIDE, exclude_nv=True)

    _COMP_CARBONATE = {'C': 1, 'O': 3}

    @property
    def has_carbonate(self):
        return self._component_search_helper(
            self._COMP_CARBONATE, exclude_nv=True)

    _COMP_AMMONIUM = {'H': 4, 'N': 1}
    _COMP_AMMONIUM_CHARGED = {'H': 4, 'N': 1, 'e-': -1}

    @property
    def has_ammonium(self):
        return self._component_search_helper(
            self._COMP_AMMONIUM, exclude_nv=True)

    @property
    def decompose_chemicals(self):
        decompose = {}

        def add_solution_chemicals():
            decompose['[OH-]'] = self._COMP_HYDROXIDE_CHARGED
            decompose['H2O'] = self._COMP_WATER

        # Find whether there are water or not
        if self.has_water:
            decompose['H2O'] = self._COMP_WATER

        if self.has_acetate:
            decompose['[CH3COO-]'] = self._COMP_ACETATE_CHARGED
            add_solution_chemicals()

        if self.has_nitrate:
            decompose['[NO3-]'] = self._COMP_NITRATE_CHARGED
            add_solution_chemicals()

        if self.has_hydroxide:
            decompose['[OH-]'] = self._COMP_HYDROXIDE_CHARGED
            add_solution_chemicals()

        if self.has_carbonate:
            decompose['CO2'] = {'C': 1, 'O': 2}

        # FIXME: material_string is different from material_formula! How to better determine decompose chemicals?
        if 'NH4' in self.material_formula:
            decompose['NH3'] = {'H': 3, 'N': 1}

        return decompose

    @property
    def exchange_chemicals(self):
        absorption = {}

        # This justifies the usage of O2 for oxide synthesis
        # For any oxides, elements may be oxidized or reduced by
        # O2 or H2, or any other gases during synthesis. However
        # Almost at all times we may write the reducing/oxidizing
        # using only O2: if O2 appears at the LHS, the element is
        # oxidized, if O2 appears at the RHS, the element is reduced.
        if 'O' in self.v_elements:
            absorption['O2'] = {'O': 2}

        return absorption
