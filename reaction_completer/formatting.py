import re

import sympy
from sympy.core.numbers import NegativeOne, One, Zero
from sympy.printing.precedence import precedence

from reaction_completer.errors import ExpressionPrintException, FormulaException
from reaction_completer.material import MaterialInformation
from reaction_completer.periodic_table import ELEMENTS

FLOAT_ROUND = 3  # 3 decimal places 0.001
_FLOAT_RE = re.compile(r"""
        (?P<sign>[-+])?             # Sign of the float
        (?=\d|\.\d)                 # Make sure there is some number following 
        (?P<int>\d*)                # Integer part
        (\.(?P<frac>\d*))?          # Fraction part
        ([eE](?P<exp>[-+]?\d+))?    # Exponential
        """, re.VERBOSE)


def nicely_print_float(f_s):
    """
    Print a float number nicely.
    :param f_s: string of the float number.
    :return:
    """
    m = _FLOAT_RE.match(f_s)
    if not m:
        raise ValueError('This is not a float!')

    integer = m.group('int') or '0'
    fraction = m.group('frac') or ''
    exp = int(m.group('exp') or 0)
    sign = m.group('sign') or ''

    while exp > 0:
        if len(fraction):
            integer += fraction[0]
            fraction = fraction[1:]
        else:
            integer += '0'
        exp -= 1
    fraction = fraction.rstrip('0')
    sign = '-' if sign == '-' else ''
    floating_number = sign + integer + ('.' + fraction if len(fraction) else '')

    return floating_number


def simplify_print(expr: sympy.Expr):
    if isinstance(expr, sympy.Float):
        # Just a float number.
        return nicely_print_float(
            str(expr.round(FLOAT_ROUND)))
    elif isinstance(expr, sympy.Add):
        expression = ''
        for i, ele in enumerate(expr.args):
            ele_str = simplify_print(ele)

            if ele_str == '0':
                continue

            if ele_str[0] == '-' or expression == '':
                expression += ele_str
            else:
                expression += '+' + ele_str

        if expression == '':
            return '0'
        else:
            return expression
    elif isinstance(expr, sympy.Mul):
        coefficient, _ = expr.as_coeff_Mul()
        if coefficient < 0:
            expr = -expr
            sign = '-'
        else:
            sign = ''

        exps = []
        for arg in expr.as_ordered_factors():
            exp = simplify_print(arg)
            if exp == '0':
                return '0'
            if exp != '1':
                if precedence(arg) < precedence(expr):
                    exp = '(%s)' % exp
                exps.append(exp)

        expression = sign + '*'.join(exps)
        return expression
    elif isinstance(expr, NegativeOne):
        return '-1'
    elif isinstance(expr, One):
        return '1'
    elif isinstance(expr, Zero):
        return '0'
    elif isinstance(expr, sympy.Symbol):
        return expr.name
    elif isinstance(expr, sympy.Integer):
        return str(expr.p)
    else:
        raise ExpressionPrintException(
            'Do not know how to print %r: %r' % (type(expr), expr))


ions_regex = re.compile('|'.join(sorted(ELEMENTS, key=lambda x: (-len(x), x))))
OMIT_IONS = {'O', 'H', 'N'}


def find_ions(string):
    return ions_regex.findall(string)


def render_reaction(precursors, target, reaction, element_substitution=None):
    element_substitution = element_substitution or {}

    left_strings = []
    # Alphabetic ordering for left part
    for material, amount in sorted(reaction['left'].items()):
        left_strings.append('%s %s' % (amount, material))
    right_strings = []
    # Alphabetic ordering for right part, except
    # for target material it's the first element
    for material, amount in sorted(
            reaction['right'].items(),
            key=lambda x: (x[0] != target['material_formula'], x[0], x[1])):
        right_strings.append('%s %s' % (amount, material))

    reaction_string = ' + '.join(left_strings) + ' == ' + ' + '.join(right_strings)

    if len(element_substitution) > 0:
        subs_string = ', '.join(['%s = %s' % (sub, subs)
                                 for (sub, subs) in element_substitution.items()])
        reaction_string += '; ' + subs_string

    # Populate additives
    if target['additives']:
        additive_ions = set(find_ions(' '.join(target['additives'])))
        added_anions = set()
        additive_precursors = []

        for precursor in precursors:
            compositions = []
            for comp in precursor['composition']:
                compositions.append({
                    'amount': comp['amount'],
                    'elements': comp['elements'],
                })

            try:
                mat_info = MaterialInformation(
                    precursor['material_string'],
                    precursor['material_formula'],
                    compositions)

                if mat_info.all_elements and \
                        any(x in additive_ions for x in mat_info.all_elements):
                    added_anions.update(mat_info.all_elements)
                    additive_precursors.append(precursor['material_formula'])
            except FormulaException:
                pass

        reaction_string += ' ; target %s with additives %s via %s' % (
            target['material_formula'],
            ', '.join(target['additives']),
            ', '.join(sorted(additive_precursors))
        )

    return reaction_string
