import sympy as sym
import numpy as np
import interval as ival
import sys
from check_box import check_box
from ExtensionClass import ClassicalKrawczykExtension, BicenteredKrawczykExtension, HansenSenguptaExtension
from plotter_support_functions import uni_plotter
import matplotlib.pyplot as plt
from LoggerClass import Logger
from sympy import sin, cos


def plot_circles(r1 = 3, r2 = 18, ang1 = None, ang2 = None):
    circle = plt.Circle((0, 0), radius=r1, fc='y', fill=False)
    plt.gca().add_patch(circle)
    circle = plt.Circle((0, 0), radius=r2, fc='y', fill=False)
    plt.gca().add_patch(circle)
    if ang1!=None and ang2!=None:
        plt.plot([r1*np.cos(ang2), r2*np.cos(ang2)], [r1*np.sin(ang2), r2*np.sin(ang2)], c="black")
        plt.plot([r1 * np.cos(ang1), r2 * np.cos(ang1)], [r1 * np.sin(ang1), r2 * np.sin(ang1)], c="black")


def symbolic_posypkin_robot_func(l1=6, l2=6, a=5, b=10.5):
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
    f = sym.Matrix([[v[0]*cos(v[1]) - u[0]],
                    [v[0]*sin(v[1]) - u[1]]])
    return f, u, v

N = 40  # The number of boxes on uniform grid
##### 2-DOF
left_v1 = 3
right_v1 = 12
left_v2 = np.pi/6
right_v2 = (2*np.pi)/3
f_sym, u_sym, v_sym = symbolic_posypkin_robot_func()
v1 = ival.Interval([left_v1, right_v1])
v2 = ival.Interval([left_v2, right_v2])
v_ival = [v1, v2]
u_upper = 20  # the width of the of the 2-dimensional square
grid = np.linspace(-u_upper, u_upper, N + 1)  # The vector to build size-dim. grid
size = 2  # The dimension of uniform grid
eps = 1e-6
coef = 2

classical_krawczyk_extension = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=False)
classical_krawczyk_loger = Logger(grid, size, v_ival, eps, classical_krawczyk_extension)
area_boxes_classical_krawczyk, border_boxes_classical_krawczyk = check_box(grid, size, v_ival,\
                                                                           classical_krawczyk_extension, eps)
uni_plotter(area_boxes_classical_krawczyk, border_boxes_classical_krawczyk, u_upper, "Classical Krawczyk", size=2, logger=classical_krawczyk_loger)
plot_circles(left_v1, right_v1, left_v2, right_v2)
classical_krawczyk_extension_elementwise = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=True)
classical_krawczyk_elementwise_loger = Logger(grid, size, v_ival, eps, classical_krawczyk_extension_elementwise)
area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise = check_box(grid, size, v_ival,\
                                                                           classical_krawczyk_extension_elementwise, eps)
uni_plotter(area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise, u_upper,
            "Classical Krawczyk Elementwise", size=2, logger=classical_krawczyk_elementwise_loger)
plot_circles(left_v1, right_v1, left_v2, right_v2)
bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
bicentered_krawczyk_loger = Logger(grid, size, v_ival, eps, bicentered_krawczyk_extension)
area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk = check_box(grid, size, v_ival,\
                                                                           bicentered_krawczyk_extension, eps)
uni_plotter(area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk, u_upper, "Bicentered Krawczyk", size=2,
            logger=bicentered_krawczyk_loger)
plot_circles(left_v1, right_v1, left_v2, right_v2)
hansen_sengupta_extension = HansenSenguptaExtension(f_sym, v_sym, u_sym, coef=1, is_elementwise=False)
hansen_sengupta_loger = Logger(grid, size, v_ival, eps, hansen_sengupta_extension, decomp=False)
area_boxes_hansen_sengupta, border_boxes_hansen_sengupta = check_box(grid, size, v_ival,\
                                                                           hansen_sengupta_extension, eps, log=False, decomposition=False)
uni_plotter(area_boxes_hansen_sengupta, border_boxes_hansen_sengupta, u_upper, "Hansen-Sengupta", size=2, logger=hansen_sengupta_loger)
plot_circles(left_v1, right_v1, left_v2, right_v2)
hansen_sengupta_extension_elementwise = HansenSenguptaExtension(f_sym, v_sym, u_sym, coef=1, is_elementwise=True)
hansen_sengupta_loger_elementwise = Logger(grid, size, v_ival, eps, hansen_sengupta_extension_elementwise, decomp=False)
area_boxes_hansen_sengupta_elementwise, border_boxes_hansen_sengupta_elementwise = check_box(grid, size, v_ival,\
                                                                           hansen_sengupta_extension_elementwise, eps, log=False, decomposition=False)
uni_plotter(area_boxes_hansen_sengupta_elementwise, border_boxes_hansen_sengupta_elementwise, u_upper, "Hansen-Sengupta Elementwise", size=2, logger=hansen_sengupta_loger_elementwise)
plot_circles(left_v1, right_v1, left_v2, right_v2)
plt.show()