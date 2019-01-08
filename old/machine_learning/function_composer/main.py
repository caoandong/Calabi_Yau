import numpy as np
import sympy as sp
from sympy import simplify


def expr_in_bank(expr, expr_bank):
    for expr_tmp in expr_bank:
        if simplify(expr-expr_tmp) == 0:
            return True
    return False

def parse_instruction(op, expr_1, expr_2, expr_bank):
    op = op % 4
    expr = expr_1
    if op == 0:
        expr_new = expr_1 + expr_2
    elif op == 1:
        expr_new = expr_1 - expr_2
    elif op == 2:
        expr_new = expr_1*expr_2
    elif op == 3:
        # careful about infty
        expr_new = expr_1/expr_2
    elif op == 4:
        # careful about indef
        expr_new = sp.sqrt(expr_1)
    # if expr_new != 0 and (not expr_in_bank(expr_new, expr_bank)):
    if expr_new != 0:
        expr_bank.append(expr_new)
        expr = expr_new
    return expr, expr_bank

def composer(instruction, expr_bank):
    expr = 0
    for line in instruction:
        op, expr_2_idx = line
        expr_new = 0
        expr_1 = expr
        expr_2 = expr_bank[expr_2_idx]
        expr, expr_bank = parse_instruction(op, expr_1, expr_2, expr_bank)
    return expr, expr_bank

def init_expr_bank(var):
    expr_bank = []
    for var1 in var:
        sqrt = sp.sqrt(var1)
        expr_bank.append(sqrt)
        for var2 in var:
            expr_tmp = var1 + var2
            if expr_tmp != 0 and (not expr_in_bank(expr_tmp, expr_bank)):
                expr_bank.append(expr_tmp)
            expr_tmp = var1 - var2
            if expr_tmp != 0 and (not expr_in_bank(expr_tmp, expr_bank)):
                expr_bank.append(expr_tmp)
            expr_tmp = var1*var2
            if expr_tmp != 0 and (not expr_in_bank(expr_tmp, expr_bank)):
                expr_bank.append(expr_tmp)
            expr_tmp = var1/var2
            if expr_tmp != 0 and (not expr_in_bank(expr_tmp, expr_bank)):
                expr_bank.append(expr_tmp)
            expr_tmp = var2/var1
            if expr_tmp != 0 and (not expr_in_bank(expr_tmp, expr_bank)):
                expr_bank.append(expr_tmp)
    return expr_bank

# Try a random function composer
def random_composer(num_var, num_ops, n_step):
    # ops = {0:"plus", 1:"minus", 2:"mult", 3:"div", 4:"sqrt"}
    ops = ["plus", "minus", "mult", "div", "sqrt"]
    var = []
    for i in range(num_var):
        var.append(sp.symbols('x'+str(i)))
    expr_bank = init_expr_bank(var)
    print 'obtained expressions:'
    print expr_bank
    instruction = []
    expr = 0
    for i in range(n_step):
        op = np.random.randint(num_ops)
        expr_2_idx = np.random.randint(len(expr_bank))
        expr_2 = expr_bank[expr_2_idx]
        line = [op, expr_2]
        expr, expr_bank =  parse_instruction(op, expr, expr_2, expr_bank)
    print 'output:', expr
    return simplify(expr)

random_composer(2, 5, 10)
