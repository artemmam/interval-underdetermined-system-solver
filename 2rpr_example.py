import sympy as sym
import numpy as np
import interval as ival
import sys
from check_box import check_box, check_one_box
from ExtensionClass import ClassicalKrawczykExtension, BicenteredKrawczykExtension, HansenSenguptaExtension
from plotter_support_functions import uni_plotter, plot_one_box
import matplotlib.pyplot as plt
from LoggerClass import Logger


def plot_circles(r1=3, r2=15, d=8):
    circle = plt.Circle((-0.5*d, 0), radius=r1, fc='y', fill=False)
    plt.gca().add_patch(circle)
    circle = plt.Circle((0.5*d, 0), radius=r1, fc='y', fill=False)
    plt.gca().add_patch(circle)
    circle = plt.Circle((-0.5*d, 0), radius=r2, fc='y', fill=False)
    plt.gca().add_patch(circle)
    circle = plt.Circle((0.5*d, 0), radius=r2, fc='y', fill=False)
    plt.gca().add_patch(circle)


def symbolic_2rpr_func(d=8):
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
    f = sym.Matrix([[v[0] ** 2 - (u[0] + 0.5 * d) ** 2 - u[1] ** 2],
                    [v[1] ** 2 - (u[0] - 0.5 * d) ** 2 - u[1] ** 2]])
    return f, u, v

N = 30  # The number of boxes on uniform grid
##### 2-rpr
f_sym, u_sym, v_sym = symbolic_2rpr_func()
v1 = ival.Interval([3., 15.])
v2 = ival.Interval([3., 15.])
v_ival = [v1, v2]
u_upper = 15  # the width of the of the 2-dimensional square
grid = np.linspace(-u_upper, u_upper, N + 1)  # The vector to build size-dim. grid
size = 2  # The dimension of uniform grid
eps = 1e-3  # accuracy
coef = 2  # Coefficient
grid_v = np.linspace(v1[0], v1[1], 10 + 1)
# box = [ival.Interval([-11, -10]), ival.Interval([0, 1])]
# bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
# check = check_one_box(box, v_ival, bicentered_krawczyk_extension, eps, log=True)
# print("$$$$$")
#
# check_one_box(box, v_ival, bicentered_krawczyk_extension, eps, log=True, strategy="Enlarge", grid=grid_v, dim=2)
# plot_one_box(box, check, u_upper)
# plot_circles()
# plt.show()
classical_krawczyk_extension = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=False)
area_boxes_classical_krawczyk, border_boxes_classical_krawczyk = check_box(grid, size, v_ival,\
                                                                           classical_krawczyk_extension, eps)
uni_plotter(area_boxes_classical_krawczyk, border_boxes_classical_krawczyk, u_upper, "Classical Krawczyk", size=2)
plot_circles()
#
# area_boxes_classical_krawczyk_infl, border_boxes_classical_krawczyk_infl = check_box(grid, size, v_ival,\
#                                                                            classical_krawczyk_extension, eps, strategy="Inflation", grid_v=grid_v, dim_v=2)
# uni_plotter(area_boxes_classical_krawczyk_infl, border_boxes_classical_krawczyk_infl, u_upper, "Classical Krawczyk Inflation", size=2)
# plot_circles()

# classical_krawczyk_extension_elementwise = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=True)
# area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise = check_box(grid, size, v_ival,\
#                                                                            classical_krawczyk_extension_elementwise, eps)
# uni_plotter(area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise, u_upper,
#             "Classical Krawczyk Elementwise", size=2)
bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk = check_box(grid, size, v_ival,\
                                                                           bicentered_krawczyk_extension, eps)
uni_plotter(area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk, u_upper, "Bicentered Krawczyk", size=2)
plot_circles()

# area_boxes_bicentered_krawczyk_infl, border_boxes_bicentered_krawczyk_infl = check_box(grid, size, v_ival,\
#                                                                            bicentered_krawczyk_extension, eps, strategy="Inflation", grid_v=grid_v, dim_v=2)
# uni_plotter(area_boxes_bicentered_krawczyk_infl, border_boxes_bicentered_krawczyk_infl, u_upper, "Bicentered Krawczyk Inflation", size=2)
# plot_circles()
#
hansen_sengupta_extension = HansenSenguptaExtension(f_sym, v_sym, u_sym, coef=1.0, is_elementwise=False)
hansen_sengupta_loger = Logger(grid, size, v_ival, eps, hansen_sengupta_extension, decomp=False)
area_boxes_hansen_sengupta, border_boxes_hansen_sengupta = check_box(grid, size, v_ival,\
                                                                           hansen_sengupta_extension, eps, log=False, decomposition=False)
uni_plotter(area_boxes_hansen_sengupta, border_boxes_hansen_sengupta, u_upper, "Hansen-Sengupta", size=2, logger=hansen_sengupta_loger)
plot_circles()
#
# area_boxes_hansen_sengupta_infl, border_boxes_hansen_sengupta_infl = check_box(grid, size, v_ival,\
#                                                                            hansen_sengupta_extension, eps, log=False, decomposition=False,
#                                                                      strategy="Inflation", grid_v=grid_v, dim_v=2)
# uni_plotter(area_boxes_hansen_sengupta_infl, border_boxes_hansen_sengupta_infl, u_upper, "Hansen-Sengupta Inflation", size=2, logger=hansen_sengupta_loger)
# plot_circles()

plt.show()