import sympy as sym
from sympy import sin, cos
from sympy.utilities.lambdify import implemented_function
import interval as ival
import numpy as np
import sys


def derive_matrix(g, v):
    """
    Function for calculating partial derivative of matrix g
    :param g : array to be derived
    :param v : variables for derivative
    :return gv: derived matrix
    """
    g_v_all = []
    for i in range(len(v)):
        g_v_all.append(sym.diff(g, v[i]))  # Calculate derivative of G with respect to v
    gv = sym.Matrix()
    for i in range(len(g_v_all)):
        gv = sym.Matrix([gv, g_v_all[i]])
    gv = gv.reshape(g.shape[0], len(v)).T
    return gv


def derived_recurrent_form(f, v, u, l):
    """
    Produce derived recurrent function
    :param f: old right-hand side
    :param v: list of checking intervals
    :param u: list of fixed intervals
    :param l: lamda matrix
    :return: function of derived recurrent form
    """
    param = [u] + [l]
    g = sym.Matrix()
    f_rec = recurrent_form(f, v, l)
    for j in range(len(v)):
        g_v = derive_matrix(sym.Matrix([f_rec[j]]), v)
        g = sym.Matrix([g, g_v])
    g = function_replacer(g)
    return sym.lambdify([v, param], g)


def mysin(x):
    """
    Interval sin
    :param x: interval
    :return: interval sin(x)
    """
    return ival.sin(x)


def mycos(x):
    """
    Interval cos
    :param x: interval
    :return: interval cos(x)
    """
    return ival.cos(x)


def recurrent_form(f, V, lam):
    """
    Produces the right-hand side of the equivalent recurrent form for the equation f(v) = 0
    :param f: old right-hand side
    :param V: variables
    :param lam: lamda matrix
    :return: the recurrent right-hand side v - lam * f(v)
    """
    v = sym.Matrix()
    for i in range(len(V)):
        v = v.row_insert(i, sym.Matrix([V[i]]))
    lam = sym.Matrix([lam]).reshape(len(V), len(V))
    rec_form = v - lam * f
    print("#####")
    print("Recurrent form ".upper())
    print("#####")
    print(rec_form)
    return rec_form  # Equivalent recurrent transformation


def function_replacer(fv):
    mysin1 = implemented_function(sym.Function('mysin1'), lambda x: mysin(x))
    mycos1 = implemented_function(sym.Function('mycos1'), lambda x: mycos(x))
    fv = fv.replace(sin, mysin1)
    fv = fv.replace(cos, mycos1)
    return fv

def derived_f(f, v, u):
    """
    Produce numerical derived function
    :param f: old right-hand side
    :param v: list of checking intervals
    :param u: list of fixed intervals
    :return: function of numerical derived form
    """
    fv = derive_matrix(f, v)
    fv = function_replacer(fv)
    return sym.lambdify([v, u], fv)


def centered_form(f, V, C, param):
    """
    Centered interval form
    :param f: old right-hand side
    :param V: list of checking intervals
    :param C: matrix of points from interval V
    :param param: parameters
    :return: function for calculating centered interval form
    """
    g_fin = sym.Matrix()
    g_fin_elementwise = []
    C = sym.Matrix([C]).reshape(len(V), len(V))
    # print("#####")
    # print("Derived F".upper())
    # print("#####")
    for i in range(len(V)):
        v = sym.Matrix()
        for j in range(len(V)):
            v = v.row_insert(j, sym.Matrix([V[j]]))
        g_v = derive_matrix(sym.Matrix([f[i]]), v)
        # print(g_v)
        c = C[i, ::].T
        v_c = v - c
        subsv = []
        for j in range(len(V)):
            subsv.append((V[j], c[j]))
        f_s = f.subs(subsv)
        g_eval = sym.Matrix([f_s[i]]) + g_v.T * v_c  # Classical central form
        g_fin = sym.Matrix([g_fin, g_eval])
        g_fin_elementwise.append(sym.lambdify([V, c, param], function_replacer(g_eval)))
    g_fin = function_replacer(g_fin)
    return sym.lambdify([V, C, param], g_fin), g_fin_elementwise


def classical_krawczyk_extension(f, u, v, l, c):

    """
    Krawczyk_evalutation function (centered form of the recurrent form of the system of the equations)
    :param f: old right-hand side
    :param u: list of fixed intervals
    :param v: list of checking intervals
    :param l: lamda matrix
    :param c: list of points in v
    :return: function of centered form from recurrent form
    """
    param = [u] + [l]
    return centered_form(recurrent_form(f, v, l), v, c, param)

def hansen_sengupta_extension(f, u, v, lam, c):
    n = len(v)
    param = [u] + [lam]
    lam = np.array(lam).reshape(n, n)
    c_matrix = np.array(c).reshape(n, n)
    c = np.array(c[0:2]).reshape(n, 1)
    c = sym.Matrix(c)
    v = sym.Matrix(v)
    subsv = []
    for j in range(n):
        subsv.append((v[j], c[j]))
    f_c = f.subs(subsv)
    g_v = derive_matrix(f, v)
    lam = sym.Matrix(lam)
    A = lam * f_c
    B = lam * g_v
    diag_elements = []
    for i in range(n):
        for j in range(n):
            if i == j:
                diag_elements.append(B[i, j])
    D = sym.diag(*diag_elements)
    M = B - D
    g = c - D**(-1)*A + D**(-1)*M*(c - v)
    g_fin = function_replacer(g)
    g = function_replacer(g)
    print("HS extension formula")
    for eq in g:
        print(eq)
    print("HS extension formula elementwise")
    g_fin_elementwise = []
    subsv = []
    for k in range(1, n):
        for j in range(n):
            subsv.append((c[j], c_matrix[k, j]))
            g[k] = g[k].subs(subsv)
    for i in range(n):
        print(g[i])
        g_fin_elementwise.append(sym.lambdify([v, c_matrix[i, ::], param], g[i]))
    return sym.lambdify([v, c, param], g_fin), g_fin_elementwise