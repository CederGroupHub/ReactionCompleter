from reaction_completer.test.reaction_tester import TestReaction


class TestSimple(TestReaction):
    def test_basic_completer(self):
        reactions = self.balance_equation([
            ("BaCO3", "BaCO3", "1.0-Ba:1.0+C:1.0+O:3.0"),
            ("TiO2", "TiO2", "1.0-Ti:1.0+O:2.0"),
        ], [
            ("BaTiO3", "BaTiO3", "1.0-Ba:1.0+Ti:1.0+O:3.0"),
        ])
        self.assertReactionsEqual(reactions, [
            '1 BaCO3 + 1 TiO2 == 1 BaTiO3 + 1 CO2'
        ])

    def test_h2o(self):
        reactions = self.balance_equation([
            ("C4H6CuO4·H2O", "C4H6CuO4·H2O", "1.0-C:4+H:6+Cu:1+O:4;1.0-H:2+O:1"),
            ("Cr2O3", "Cr2O3", "1-Cr:2+O:3"),
        ], [
            ("CuCrO2", "CuCrO2", "1-Cr:1+Cu:1+O:2"),
        ])
        self.assertReactionsEqual(reactions, [
            '1 C4H6CuO4·H2O + 0.5 Cr2O3 + 3.75 O2 == 1 CuCrO2 + 4 CO2 + 4 H2O'
        ])


class TestElementSubstitution(TestReaction):
    def test_simple(self):
        """
        SrCO3 + Al2O3 + Fe2O3 ==== Sr6(A2O4)6, A=Al, Fe
        """
        reactions = self.balance_equation([
            ('SrCO3', 'SrCO3', '1.0-Sr:1.0+C:1.0+O:3.0'),
            ('Al2O3', 'Al2O3', '1.0-Al:2.0+O:3.0'),
            ('MnO', 'MnO', '1.0-Mn:1.0+O:1.0'),
            ('Fe2O3', 'Fe2O3', '1.0-Fe:2.0+O:3.0'),
            ('ZrO2', 'ZrO2', '1.0-Zr:1.0+O:2.0'),
            ('H2O', 'H2O', '1.0-H:2.0+O:2.0'),
        ], [
            ('Sr6(A2O4)6', 'Sr6(A2O4)6', '1.0-A:12.0+O:24.0+Sr:6.0?A:Fe,Al<Mn2+')
        ], [
            "SrCO3, Al2O3, MnO and Fe2O3 are used to synthesize Mn2+doped-Sr6(A2O4)6, A=Fe, Al.",
            "Milling media is ZrO2",
            "There is some H2O found in the final product."
        ])

        self.assertReactionsEqual(reactions, [
            '6 Fe2O3 + 6 SrCO3 == 1 Sr6(A2O4)6 + 6 CO2; A = Fe ; target Sr6(A2O4)6 with additives Mn2+ via MnO',
            '6 Al2O3 + 6 SrCO3 == 1 Sr6(A2O4)6 + 6 CO2; A = Al ; target Sr6(A2O4)6 with additives Mn2+ via MnO'
        ])
