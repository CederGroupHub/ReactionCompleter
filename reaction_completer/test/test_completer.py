import re
from unittest import TestCase

from reaction_completer import balance_recipe


def simple_parse(composition_string):
    """
    Parse compact composition strings such as:
    1.0-D:1.0+C:1.0+O:3.0;
    2.0-Fe:2.0+O:3.0?
    D:K,Na;E:O,H<
    Na+&K+
    """
    composition_string = re.sub(r'\s', '', composition_string)

    additives = []
    if '<' in composition_string:
        composition_string, additives_string = composition_string.split('<')
        additives = additives_string.split('&')

    elements_vars = {}
    if '?' in composition_string:
        composition_string, elements_vars_string = composition_string.split('?')
        for element_var in elements_vars_string.split(';'):
            element, values = element_var.split(':')
            elements_vars[element] = values.split(',')

    compositions = composition_string.split(';')
    material_composition = []
    for composition in compositions:
        comp_amount, elements = composition.split('-')
        elements_dict = {}
        for pair in elements.split('+'):
            element, amount = pair.split(':')
            elements_dict[element] = amount
        material_composition.append({
            'amount': comp_amount,
            'elements': elements_dict
        })
    return material_composition, elements_vars, additives


def make_material(data):
    material_formula, material_string, composition_string = data
    composition, elements_vars, additives = simple_parse(composition_string)
    return {
        'material_formula': material_formula,
        'material_string': material_string,
        'composition': composition,
        'elements_vars': elements_vars,
        'additives': additives,
    }


def make_materials(data):
    return [make_material(x) for x in data]


class TestSimple(TestCase):
    def test_basic_completer(self):
        """
        BaCO3 + TiO2 ==== BaTiO3 + CO2
        """
        precursors = make_materials([
            ("BaCO3", "BaCO3", "1.0-Ba:1.0+C:1.0+O:3.0"),
            ("TiO2", "TiO2", "1.0-Ti:1.0+O:2.0"),
        ])
        targets = make_materials([
            ("BaTiO3", "BaTiO3", "1.0-Ba:1.0+Ti:1.0+O:3.0"),
        ])

        reactions = balance_recipe(precursors, targets)
        self.assertListEqual(reactions, [(
            'BaTiO3',
            {
                'left': {'BaCO3': '1', 'TiO2': '1'},
                'right': {'BaTiO3': '1', 'CO2': '1'}
            },
            None,
            '1 BaCO3 + 1 TiO2 == 1 BaTiO3 + 1 CO2'
        )])


class TestElementSubstitution(TestCase):
    def test_simple(self):
        """
        SrCO3 + Al2O3 + Fe2O3 ==== Sr6(A2O4)6, A=Al, Fe
        """
        precursors = make_materials([
            ('SrCO3', 'SrCO3', '1.0-Sr:1.0+C:1.0+O:3.0'),
            ('Al2O3', 'Al2O3', '1.0-Al:2.0+O:3.0'),
            ('MnO', 'MnO', '1.0-Mn:1.0+O:1.0'),
            ('Fe2O3', 'Fe2O3', '1.0-Fe:2.0+O:3.0'),
            ('ZrO2', 'ZrO2', '1.0-Zr:1.0+O:2.0'),
            ('H2O', 'H2O', '1.0-H:2.0+O:2.0'),
        ])
        targets = make_materials([
            ('Sr6(A2O4)6', 'Sr6(A2O4)6', '1.0-A:12.0+O:24.0+Sr:6.0?A:Fe,Al<Mn2+')
        ])
        text = [
            "SrCO3, Al2O3, MnO and Fe2O3 are used to synthesize Mn2+doped-Sr6(A2O4)6, A=Fe, Al.",
            "Milling media is ZrO2",
            "There is some H2O found in the final product."
        ]

        reactions = balance_recipe(precursors, targets, text)
        self.assertListEqual(reactions, [
            (
                'Sr6(A2O4)6', {
                    'left': {'SrCO3': '6', 'Fe2O3': '6'},
                    'right': {'Sr6(A2O4)6': '1', 'CO2': '6'}
                },
                {'A': 'Fe'},
                '6 Fe2O3 + 6 SrCO3 == 1 Sr6(A2O4)6 + 6 CO2; A = Fe ; target Sr6(A2O4)6 with additives Mn2+ via MnO'
            ),
            (
                'Sr6(A2O4)6', {
                    'left': {'SrCO3': '6', 'Al2O3': '6'},
                    'right': {'Sr6(A2O4)6': '1', 'CO2': '6'}
                },
                {'A': 'Al'},
                '6 Al2O3 + 6 SrCO3 == 1 Sr6(A2O4)6 + 6 CO2; A = Al ; target Sr6(A2O4)6 with additives Mn2+ via MnO'
            )
        ])
