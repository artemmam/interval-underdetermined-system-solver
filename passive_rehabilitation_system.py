import sympy as sym
import numpy as np
import interval as ival
import sys
from check_box import check_box
from ExtensionClass import ClassicalKrawczykExtension, BicenteredKrawczykExtension, HansenSenguptaExtension
from plotter_support_functions import uni_plotter
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from LoggerClass import Logger
from sympy import sin, cos
import math
import itertools as it


def distance(x1, y1, x2, y2):
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def angle(vector1, vector2):
    x1, y1 = vector1
    x2, y2 = vector2
    inner_product = x1*x2 + y1*y2
    len1 = math.hypot(x1, y1)
    len2 = math.hypot(x2, y2)
    return math.degrees(math.acos(inner_product/(len1*len2)))


def get_coordinates(al, beta, a=4, b=2):
    x = a*np.cos(al) + b*np.cos(al + beta)
    y = a*np.sin(al) + b*np.sin(al + beta)
    return x, y


def plot_area(ax1, a = 4, b = 2, ang1_0 = None, ang1_1 = None,  ang2_0 = None, ang2_1 = None):
    circle = plt.Circle((0, 0), radius=abs(a - b), fc='y', fill=False)
    plt.gca().add_patch(circle)
    circle = plt.Circle((0, 0), radius=a + b, fc='y', fill=False)
    plt.gca().add_patch(circle)
    A = list(it.product([ang1_0, ang1_1], [ang2_0, ang2_1]))
    x = []
    y = []
    # print(A)
    dots = []
    for el in A:
        print(math.degrees(el[0]), math.degrees(el[1]))
        XY = get_coordinates(*el, a, b)
        print(XY)
        dots.append(XY)
    dots = np.array(dots)
    # print(dots)
    # print(dots[::, 0])
    # print(dots[::, 1])
    plt.scatter(dots[::, 0], dots[::, 1])
    # ang1 = angle(dots[0], dots[3])
    # ang2 = angle(dots[1], dots[2])
    # print(ang1, ang2)
    # r1 = distance(0, 0, *dots[0])
    # r2 = distance(0, 0, *dots[1])
    # print(r1, r2)
    # ang3 = angle(dots[0], (r1, 0))
    # ang4 = angle(dots[1], (r2, 0))
    # print(ang3, ang4)
    # a = Arc((0, 0), 2 * r1, 2 * r1, 0, ang3, ang3 + ang1, color='red', lw=1)
    # ax1.add_patch(a)
    # a1 = Arc((0, 0), 2 * r2, 2 * r2, 0, ang4, ang4 + ang2, color='red', lw=1)
    # ax1.add_patch(a1)
    # ax1.plot([dots[0][0], dots[1][0]], [dots[0][1], dots[1][1]], color='red')
    # ax1.plot([dots[2][0], dots[3][0]], [dots[2][1], dots[3][1]], color='red')


def symbolic_pasive_rehabilitation_system_func(l_a=4, l_b=2):
    """
    Creating symbol variables for circle eq. system
    :return: symbolic eq. system,
            symbolic U (fixed boxes),
            symbolic V (checking boxes),
            symbolic C
    """
    v = sym.symbols('v1, v2')
    u = sym.symbols('u1, u2')
    f = sym.Matrix([[u[0] - l_a * cos(v[0]) - l_b * cos(v[0] + v[1])],
                    [u[1] - l_a * sin(v[0]) - l_b * sin(v[0] + v[1])]]
                   )
    return f, u, v

N = 5  # The number of boxes on uniform grid
##### 2-DOF
a = 4
b = 2
left_v1 = np.pi/6
right_v1 = np.pi
left_v2 = np.pi/6
right_v2 = np.pi
f_sym, u_sym, v_sym = symbolic_pasive_rehabilitation_system_func(a, b)
v1 = ival.Interval([left_v1, right_v1])
v2 = ival.Interval([left_v2, right_v2])
v_ival = [v1, v2]
u_upper = 10  # the width of the of the 2-dimensional square
grid = np.linspace(2, 3, N + 1)  # The vector to build size-dim. grid
size = 2  # The dimension of uniform grid
eps = 1e-6
coef = 2

fig1 = plt.figure(figsize=(8, 8))
ax1 = fig1.add_subplot(1, 1, 1)
# classical_krawczyk_extension = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=False)
# area_boxes_classical_krawczyk, border_boxes_classical_krawczyk = check_box(grid, size, v_ival,\
#                                                                            classical_krawczyk_extension, eps)
# uni_plotter(area_boxes_classical_krawczyk, border_boxes_classical_krawczyk, u_upper, "Classical Krawczyk", size=2)
# plot_area(ax1, a, b, left_v1, right_v1, left_v2, right_v2)
# classical_krawczyk_extension_elementwise = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=True)
# area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise = check_box(grid, size, v_ival,\
#                                                                            classical_krawczyk_extension_elementwise, eps)
# uni_plotter(area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise, u_upper,
#             "Classical Krawczyk Elementwise", size=2)
# plot_circles()
bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
bicentered_krawczyk_loger = Logger(grid, size, v_ival, eps, bicentered_krawczyk_extension, decomp=False)
area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk = check_box(grid, size, v_ival,\
                                                                           bicentered_krawczyk_extension, eps)
uni_plotter(area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk, u_upper, "Bicentered Krawczyk", size=2,
            logger=bicentered_krawczyk_loger, ax=ax1, fig=fig1)
# plot_circles(left_v1, right_v1, left_v2, right_v2)
# hansen_sengupta_extension = HansenSenguptaExtension(f_sym, v_sym, u_sym, coef=1, is_elementwise=False)
# hansen_sengupta_loger = Logger(grid, size, v_ival, eps, hansen_sengupta_extension, decomp=False)
# area_boxes_hansen_sengupta, border_boxes_hansen_sengupta = check_box(grid, size, v_ival,\
#                                                                            hansen_sengupta_extension, eps, log=False, decomposition=False)
# uni_plotter(area_boxes_hansen_sengupta, border_boxes_hansen_sengupta, u_upper, "Hansen-Sengupta", size=2, logger=hansen_sengupta_loger)
# plot_circles()
plot_area(ax1, a, b, left_v1, right_v1, left_v2, right_v2)
plt.show()