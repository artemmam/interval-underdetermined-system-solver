import sympy as sym
import numpy as np
import interval as ival
import sys
from check_box import check_box, check_one_box, reccur_func_enlarge
from LoggerClass import Logger
from ExtensionClass import ClassicalKrawczykExtension, BicenteredKrawczykExtension, HansenSenguptaExtension
from plotter_support_functions import uni_plotter, plot_one_box
import matplotlib.pyplot as plt

def symbolic_sphere_func():
    """
    Creating symbol variables for circle eq. system
    :return: symbolic eq. system,
            symbolic U (fixed boxes),
            symbolic V (checking boxes),
            symbolic Vmid,
            symbolic C
    """
    v = [sym.symbols('v1')]
    u = sym.symbols('u1, u2')
    f = sym.Matrix([(u[0]) ** 2 +
                    (u[1]) ** 2 + v[0] ** 2 - 1])
    return f, u, v

N = 30  # The number of boxes on uniform grid
##### 1d circle
f_sym, u_sym, v_sym = symbolic_sphere_func()
v1 = ival.Interval([0, 1.0])
v_ival = [v1]
u_upper = 1.2  # the width of the of the 2-dimensional square
#derived_reccurent_form = derived_reccurent_form(f, V, U, Vmid)
grid = np.linspace(-u_upper, u_upper, N + 1)  # The vector to build size-dim. grid
#print(grid)
size = 2  # The dimension of uniform grid
eps = 1e-6  # accuracy
coef = 2  # Coefficient

# box = [ival.Interval([0.8, 0.9]), ival.Interval([0.5, 0.6])]
# bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
# check = check_one_box(box, v_ival, bicentered_krawczyk_extension, eps, log=True)
# print("$$$$$")
# Inflation
grid_v = np.linspace(v1[0], v1[1], 10 + 1)
# check_one_box(box, v_ival, bicentered_krawczyk_extension, eps, log=True, strategy="Enlarge", grid=grid_v, dim=1)
# plot_one_box(box, check, u_upper)
# circle2 = plt.Circle((0, 0), 1, fc="y", fill = False)
# plt.gca().add_patch(circle2)
# plt.show()
classical_krawczyk_extension = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=False)
classical_krawczyk_loger = Logger(grid, size, v_ival, eps, classical_krawczyk_extension, decomp=False)
area_boxes_classical_krawczyk_infl, border_boxes_classical_krawczyk_infl = check_box(grid, size, v_ival,\
                                                                           classical_krawczyk_extension, eps, log=False,
                                                                           decomposition=False, strategy="Inflation", grid_v=grid_v, dim_v=1)
uni_plotter(area_boxes_classical_krawczyk_infl, border_boxes_classical_krawczyk_infl, u_upper,
            "Classical Krawczyk Infaltion", size=2, logger=classical_krawczyk_loger)
circle2 = plt.Circle((0, 0), 1, fc="y", fill = False)
plt.gca().add_patch(circle2)
area_boxes_classical_krawczyk, border_boxes_classical_krawczyk = check_box(grid, size, v_ival,\
                                                                           classical_krawczyk_extension, eps, log=False,
                                                                           decomposition=False)
uni_plotter(area_boxes_classical_krawczyk, border_boxes_classical_krawczyk, u_upper,
            "Classical Krawczyk", size=2, logger=classical_krawczyk_loger)
circle2 = plt.Circle((0, 0), 1, fc="y", fill = False)
plt.gca().add_patch(circle2)
# classical_krawczyk_extension_elementwise = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=True)
# area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise = check_box(grid, size, v_ival,\
#                                                                            classical_krawczyk_extension_elementwise, eps)
# uni_plotter(area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise, u_upper,
#             "Classical Krawczyk Elementwise", size=2)
# circle3 = plt.Circle((0, 0), 1, fc="y", fill = False)
# plt.gca().add_patch(circle3)
bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
area_boxes_bicentered_krawczyk_infl, border_boxes_bicentered_krawczyk_infl = check_box(grid, size, v_ival,\
                                                                           bicentered_krawczyk_extension, eps,
                                                                             strategy="Inflation", grid_v=grid_v, dim_v=1)
uni_plotter(area_boxes_bicentered_krawczyk_infl, border_boxes_bicentered_krawczyk_infl, u_upper, "Bicentered Krawczyk Inflation", size=2)
circle2 = plt.Circle((0, 0), 1, fc="y", fill = False)
plt.gca().add_patch(circle2)

area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk = check_box(grid, size, v_ival,\
                                                                           bicentered_krawczyk_extension, eps)
uni_plotter(area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk, u_upper, "Bicentered Krawczyk", size=2)
circle2 = plt.Circle((0, 0), 1, fc="y", fill = False)
plt.gca().add_patch(circle2)
#
# hansen_sengupta_extension = HansenSenguptaExtension(f_sym, v_sym, u_sym, coef=1, is_elementwise=False)
# hansen_sengupta_loger = Logger(grid, size, v_ival, eps, hansen_sengupta_extension, decomp=False)
# area_boxes_hansen_sengupta, border_boxes_hansen_sengupta = check_box(grid, size, v_ival,\
#                                                                      hansen_sengupta_extension, eps,
#                                                                      log=False, decomposition=False)
# uni_plotter(area_boxes_hansen_sengupta, border_boxes_hansen_sengupta, u_upper,
#             "Hansen-Sengupta", size=2, logger=hansen_sengupta_loger)
# circle2 = plt.Circle((0, 0), 1, fc="y", fill = False)
# plt.gca().add_patch(circle2)
plt.show()