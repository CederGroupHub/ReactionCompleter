from reaction_completer.test.reaction_tester import TestReaction


class TestSolutionBased(TestReaction):
    def test_acetate(self):
        reactions = self.balance_equation([
            ('Co(CH3COO)2·4H2O', 'Co(CH3COO)2·4H2O', 'Co+C:4+H:6+O:4;4--H:2+O'),
            ('Li(CH3COO)·2H2O', 'Li(CH3COO)·2H2O', 'Li+C:2+H:3+O:2;2--H:2+O'),
        ], [
            ('LiCoO2', 'LiCoO2', 'Li+Co+O:2')
        ])

        self.assertReactionsEqual(reactions, [
            '1 Co(CH3COO)2·4H2O + 1 Li(CH3COO)·2H2O + 0.25 O2 + 3 [OH-] == '
            '1 LiCoO2 + 7.5 H2O + 3 [CH3COO-]'
        ])

    def test_nitrate(self):
        reactions = self.balance_equation([
            ('Sm(NO3)3', 'Sm(NO3)3', 'Sm+N:3+O:9'),
            ('Co(NO3)3', 'Co(NO3)3', 'Co+N:3+O:9'),
            ('Sr(NO3)2', 'Sr(NO3)2', 'Sr+N:2+O:6')
        ], [
            ('Sm1-xSrxCoO3', 'Sm1-xSrxCoO3', 'Sm:1-x+Sr:x+Co+O:3')
        ])

        self.assertReactionsEqual(reactions, [
            '1 Co(NO3)3 + 0.25*x O2 + 1-x Sm(NO3)3 + x Sr(NO3)2 + 6-x [OH-] == '
            '1 Sm1-xSrxCoO3 + 3-0.5*x H2O + 6-x [NO3-]'
        ])

    def test_hydroxide(self):
        reactions = self.balance_equation([
            ('LiOH·H2O', 'LiOH·H2O', 'Li+O+H;H:2+O'),
            ('Mn(NO3)2·6H2O', 'Mn(NO3)2·6H2O', 'Mn+N:2+O:6;6--H:2+O'),
            ('Fe(NO3)3·9H2O', 'Fe(NO3)3·9H2O', 'Fe+N:3+O:9;9--H:2+O')
        ], [
            ('LiFeO2-Li2MnO3', 'LiFeO2-Li2MnO3', 'Li+Fe+O:2;Li:2+Mn+O:3')
        ])

        self.assertReactionsEqual(reactions, [
            '1 Fe(NO3)3·9H2O + 3 LiOH·H2O + 1 Mn(NO3)2·6H2O + 0.5 O2 + 5 [OH-] == '
            '1 LiFeO2-Li2MnO3 + 22 H2O + 5 [NO3-]'
        ])

    def test_ammonia(self):
        reactions = self.balance_equation([
            ('Li(COOCH3)', 'Li(COOCH3)', 'Li+C:2+O:2+H:3'),
            ('Mn(COOCH3)2·H2O', 'Mn(COOCH3)2·H2O', 'Mn+C:4+O:4+H:6;H:2+O'),
            ('NH4H2PO4', 'NH4H2PO4', 'N:1+H:6+P+O:4')
        ], [
            ('LiMnPO4', 'LiMnPO4', 'Li+Mn+P+O:4')
        ])

        self.assertReactionsEqual(reactions, [
            '1 Li(COOCH3) + 1 Mn(COOCH3)2·H2O + 1 NH4H2PO4 + 3 [OH-] == '
            '1 LiMnPO4 + 4 H2O + 1 NH3 + 3 [CH3COO-]'
        ])


class TestSimple(TestReaction):
    def test_carbonate(self):
        reactions = self.balance_equation([
            ("BaCO3", "BaCO3", "Ba:1.0+C:1.0+O:3.0"),
            ("TiO2", "TiO2", "Ti:1.0+O:2.0"),
        ], [
            ("BaTiO3", "BaTiO3", "Ba:1.0+Ti:1.0+O:3.0"),
        ])
        self.assertReactionsEqual(reactions, [
            '1 BaCO3 + 1 TiO2 == 1 BaTiO3 + 1 CO2'
        ])

    def test_h2o(self):
        reactions = self.balance_equation([
            ("CuO·H2O", "CuO·H2O", "Cu+O:2;1.0--H:2+O:1"),
            ("Cr2O3", "Cr2O3", "Cr:2+O:3"),
        ], [
            ("CuCrO2", "CuCrO2", "Cr:1+Cu:1+O:2"),
        ])
        self.assertReactionsEqual(reactions, [
            '1 CuO·H2O + 0.5 Cr2O3 == 1 CuCrO2 + 0.75 O2 + 1 H2O'
        ])


class TestElementSubstitution(TestReaction):
    def test_simple(self):
        """
        SrCO3 + Al2O3 + Fe2O3 ==== Sr6(A2O4)6, A=Al, Fe
        """
        reactions = self.balance_equation([
            ('SrCO3', 'SrCO3', 'Sr:1.0+C:1.0+O:3.0'),
            ('Al2O3', 'Al2O3', 'Al:2.0+O:3.0'),
            ('MnO', 'MnO', 'Mn:1.0+O:1.0'),
            ('Fe2O3', 'Fe2O3', 'Fe:2.0+O:3.0'),
            ('ZrO2', 'ZrO2', 'Zr:1.0+O:2.0'),
            ('H2O', 'H2O', 'H:2.0+O:2.0'),
        ], [
            ('Sr6(A2O4)6', 'Sr6(A2O4)6', 'A:12.0+O:24.0+Sr:6.0?A:Fe,Al<Mn2+')
        ], [
            "SrCO3, Al2O3, MnO and Fe2O3 are used to synthesize Mn2+doped-Sr6(A2O4)6, A=Fe, Al.",
            "Milling media is ZrO2",
            "There is some H2O found in the final product."
        ])

        self.assertReactionsEqual(reactions, [
            '6 Fe2O3 + 6 SrCO3 == 1 Sr6(A2O4)6 + 6 CO2; A = Fe ; target Sr6(A2O4)6 with additives Mn2+ via MnO',
            '6 Al2O3 + 6 SrCO3 == 1 Sr6(A2O4)6 + 6 CO2; A = Al ; target Sr6(A2O4)6 with additives Mn2+ via MnO'
        ])
