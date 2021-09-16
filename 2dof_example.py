import sympy as sym
import numpy as np
import interval as ival
import sys
from check_box import check_box
from ExtensionClass import ClassicalKrawczykExtension, BicenteredKrawczykExtension, HansenSenguptaExtension
from plotter_support_functions import uni_plotter
import matplotlib.pyplot as plt
from LoggerClass import Logger


def plot_circles(r1 = 3, r2 = 18, a = 5, l=6):
    circle = plt.Circle((-(r2 - r1)/2, -a/2), radius=l, fc='y', fill=False)
    plt.gca().add_patch(circle)
    circle = plt.Circle(((r2 - r1)/2, -a/2), radius=l, fc='y', fill=False)
    plt.gca().add_patch(circle)
    circle = plt.Circle((-(r2 - r1)/2, a/2), radius=l, fc='y', fill=False)
    plt.gca().add_patch(circle)
    circle = plt.Circle(((r2 - r1)/2, a/2), radius=l, fc='y', fill=False)
    plt.gca().add_patch(circle)


def symbolic_2dof_func(l1=6, l2=6, a=5, b=10.5):
    """
    Creating symbol variables for circle eq. system
    :return: symbolic eq. system,
            symbolic U (fixed boxes),
            symbolic V (checking boxes),
            symbolic Vmid,
            symbolic C
    """
    v = sym.symbols('v1, v2')
    u = sym.symbols('u1, u2')
    f = sym.Matrix([[-(u[0] - v[0] + b)**2 - (u[1] - a/2)**2 + l1**2],
                    [-(u[0] - v[1] + b)**2 - (u[1] + a/2)**2 + l2**2]])
    return f, u, v

N = 41  # The number of boxes on uniform grid
##### 2-DOF
left = 3
right = 12
f_sym, u_sym, v_sym = symbolic_2dof_func(b=(right - left)/2 + left)
v1 = ival.Interval([left, right])
v2 = ival.Interval([left, right])
v_ival = [v1, v2]
u_upper = 15  # the width of the of the 2-dimensional square
grid = np.linspace(-20, 20, N + 1)  # The vector to build size-dim. grid
size = 2  # The dimension of uniform grid
eps = 1e-6
coef = 2
classical_krawczyk_extension = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=False)
area_boxes_classical_krawczyk, border_boxes_classical_krawczyk = check_box(grid, size, v_ival,\
                                                                           classical_krawczyk_extension, eps)
uni_plotter(area_boxes_classical_krawczyk, border_boxes_classical_krawczyk, u_upper, "Classical Krawczyk", size=2)
plot_circles(left, right)
# classical_krawczyk_extension_elementwise = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=True)
# area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise = check_box(grid, size, v_ival,\
#                                                                            classical_krawczyk_extension_elementwise, eps)
# uni_plotter(area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise, u_upper,
#             "Classical Krawczyk Elementwise", size=2)
# plot_circles()
bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
bicentered_krawczyk_loger = Logger(grid, size, v_ival, eps, bicentered_krawczyk_extension)
area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk = check_box(grid, size, v_ival,\
                                                                           bicentered_krawczyk_extension, eps)
uni_plotter(area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk, u_upper, "Bicentered Krawczyk", size=2,
            logger=bicentered_krawczyk_loger)
plot_circles(left, right)
hansen_sengupta_extension = HansenSenguptaExtension(f_sym, v_sym, u_sym, coef=1, is_elementwise=False)
hansen_sengupta_loger = Logger(grid, size, v_ival, eps, hansen_sengupta_extension, decomp=False)
area_boxes_hansen_sengupta, border_boxes_hansen_sengupta = check_box(grid, size, v_ival,\
                                                                           hansen_sengupta_extension, eps, log=False, decomposition=False)
uni_plotter(area_boxes_hansen_sengupta, border_boxes_hansen_sengupta, u_upper, "Hansen-Sengupta", size=2, logger=hansen_sengupta_loger)
plot_circles(left, right)
plt.show()