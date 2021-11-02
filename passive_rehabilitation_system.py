import sympy as sym
import numpy as np
import interval as ival
import sys
from check_box import check_box, check_one_box, diam
from ExtensionClass import ClassicalKrawczykExtension, BicenteredKrawczykExtension, HansenSenguptaExtension
from plotter_support_functions import uni_plotter, plot_one_box
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from LoggerClass import Logger
from sympy import sin, cos
import math
import itertools as it
import argparse


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


def get_coordinates_first_bar(al, a=4):
    x = a*np.cos(al)
    y = a*np.sin(al)
    return x, y


def plot_area(ax1, a = 4, b = 2, ang1_0 = None, ang1_1 = None,  ang2_0 = None, ang2_1 = None):
    circle = plt.Circle((0, 0), radius=abs(a - b), fc='y', fill=False)
    plt.gca().add_patch(circle)
    circle = plt.Circle((0, 0), radius=a + b, fc='y', fill=False)
    plt.gca().add_patch(circle)
    A = list(it.product([ang1_0, ang1_1], [ang2_0, ang2_1]))
    dots = []
    for el in A:
        XY = get_coordinates(*el, a, b)
        dots.append(XY)
    dots = np.array(dots)
    # plt.scatter(dots[::, 0], dots[::, 1])
    # print(dots)
    ang1 = angle(dots[0], dots[2])
    ang2 = angle(dots[1], dots[3])
    r1 = distance(0, 0, *dots[0])
    r2 = distance(0, 0, *dots[1])
    ang3 = angle(dots[0], (r1, 0))
    ang4 = angle(dots[1], (r2, 0))
    arc1 = Arc((0, 0), 2 * r1, 2 * r1, 0, ang3, ang3 + ang1, color='red', lw=1)
    ax1.add_patch(arc1)
    arc2 = Arc((0, 0), 2 * r2, 2 * r2, 0, ang4, ang4 + ang2, color='red', lw=1)
    ax1.add_patch(arc2)
    # ax1.plot([dots[0][0], dots[1][0]], [dots[0][1], dots[1][1]], color='red')
    # ax1.plot([dots[2][0], dots[3][0]], [dots[2][1], dots[3][1]], color='red')
    ax1.scatter(dots[::, 0], dots[::, 1], color="red")
    dots_bar1 = []
    for el in [ang1_0, ang1_1]:
        dots_bar1.append(get_coordinates_first_bar(el, a))
    dots_bar1 = np.array(dots_bar1)
    for j in range(2):
        for i in range(2):
            X1 = [0, dots_bar1[j, 0], dots[i + j*2, 0]]
            Y1 = [0, dots_bar1[j, 1], dots[i + j*2, 1]]
            ax1.plot(X1, Y1, linestyle='--', color="blue")
    ax1.scatter(dots_bar1[::, 0], dots_bar1[::, 1], color="blue")
    # print(dots_bar1)
    r1_bar1 = distance(*dots_bar1[0], *dots[0])
    arc3 = Arc(dots_bar1[0], 2 * r1_bar1, 2 * r1_bar1, 0, math.degrees(ang1_0 + ang2_0), math.degrees(ang2_1 + ang1_0), color='red', lw=1)
    ax1.add_patch(arc3)
    r2_bar1 = distance(*dots_bar1[1], *dots[2])
    arc4 = Arc(dots_bar1[1], 2 * r2_bar1, 2 * r2_bar1, 0, math.degrees(ang1_1 + ang2_0), math.degrees(ang1_1 + ang2_1),
               color='red', lw=1)
    ax1.add_patch(arc4)
    ax1.set_xlabel("u1")
    ax1.set_ylabel("u2")
    ax1.arrow(0, 0, 0, 6.1, head_width=0.15, head_length=0.1, fc='k', ec='k')
    ax1.arrow(0, 0, 6.1, 0, head_width=0.15, head_length=0.1, fc='k', ec='k')



def get_minmax_xy(a = 4, b = 2, ang1_0 = None, ang1_1 = None,  ang2_0 = None, ang2_1 = None):
    A = list(it.product([ang1_0, ang1_1], [ang2_0, ang2_1]))
    dots = []
    for el in A:
        XY = get_coordinates(*el, a, b)
        dots.append(XY)
    dots = np.array(dots)
    return min(dots[::, 0]), max(dots[::, 0]), min(dots[::, 1]), max(dots[::, 1])


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

N = 15  # The number of boxes on uniform grid
##### 2-DOF
a = 4
b = 2
parser = argparse.ArgumentParser(description="Angles in radians")
parser.add_argument('-v1_0', dest="v1_0", type=float)
parser.add_argument('-v1_1', dest="v1_1", type=float)
parser.add_argument('-v2_0', dest="v2_0", type=float)
parser.add_argument('-v2_1', dest="v2_1", type=float)
args = parser.parse_args()
print(args)
left_v1 = math.radians(args.v1_0)
right_v1 = math.radians(args.v1_1)
left_v2 = math.radians(args.v2_0)
right_v2 = math.radians(args.v2_1)
f_sym, u_sym, v_sym = symbolic_pasive_rehabilitation_system_func(a, b)
v1 = ival.Interval([left_v1, right_v1])
v2 = ival.Interval([left_v2, right_v2])
v_ival = [v1, v2]
# u_lims = 6  # the width of the of the 2-dimensional square
ux_lower, ux_upper, uy_lower, uy_upper = get_minmax_xy(a, b, left_v1, right_v1, left_v2, right_v2)
u_x = [ux_lower - 1, ux_upper + 1]
u_y = [uy_lower - 1, uy_upper + 1]
u_lims = [u_x, u_y]
grid_u1 = np.linspace(u_x[0], u_x[1], N + 1)  # The vector to build size-dim. grid
grid_u2 = np.linspace(u_y[0], u_y[1], N + 1)  # The vector to build size-dim. grid
grid = [grid_u1, grid_u2]
size = 2  # The dimension of uniform grid
eps = 1e-6
coef = 2
Nv = 10
grid_v1 = np.linspace(v1[0], v1[1], Nv + 1)
grid_v2 = np.linspace(v2[0], v2[1], Nv + 1)
grid_v = [grid_v1, grid_v2]
# print(grid_v)
# fig1 = plt.figure(figsize=(8, 8))
# ax1 = fig1.add_subplot(1, 1, 1)
# ^^^^^
# box = [ival.Interval([2.008, 2.463]), ival.Interval([4.366, 4.593])] # if diam <0.82, then border
# 45 90 45 90
# [1.552, 2.008] [4.593, 4.821] border (inside)
# [2.008, 2.463], [4.366, 4.593] - outside (border)
# [-0.269, 0.187], [5.276, 5.504] border (inside)
# print(diam(box))
# bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False, log=True)
# check_one_box(box, v_ival, bicentered_krawczyk_extension, eps, log=True)
# print("$$$$$")

# check = check_one_box(box, v_ival, bicentered_krawczyk_extension, eps, log=True, strategy="Enlarge", grid=grid_v, dim=2, uniform=False)
# plot_one_box(box, check, u_lims, ax=ax1)
# plot_area(ax1, a, b, left_v1, right_v1, left_v2, right_v2)
# ^^^^^
# fig1 = plt.figure(figsize=(8, 8))
# ax1 = fig1.add_subplot(1, 1, 1)
# classical_krawczyk_extension = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=False)
# area_boxes_classical_krawczyk, border_boxes_classical_krawczyk = check_box(grid, size, v_ival,\
#                                                                            classical_krawczyk_extension, eps, uniform_u=False)
# uni_plotter(area_boxes_classical_krawczyk, border_boxes_classical_krawczyk, u_lims, "Classical Krawczyk", size=2)
# plot_area(ax1, a, b, left_v1, right_v1, left_v2, right_v2)
# classical_krawczyk_extension_elementwise = ClassicalKrawczykExtension(f_sym, v_sym, u_sym, is_elementwise=True)
# area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise = check_box(grid, size, v_ival,\
#                                                                            classical_krawczyk_extension_elementwise, eps)
# uni_plotter(area_boxes_classical_krawczyk_elementwise, border_boxes_classical_krawczyk_elementwise, u_upper,
#             "Classical Krawczyk Elementwise", size=2)
# plot_circles()
###
fig1 = plt.figure(figsize=(8, 8))
ax1 = fig1.add_subplot(1, 1, 1)
# bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
# bicentered_krawczyk_loger = Logger(grid, size, v_ival, eps, bicentered_krawczyk_extension, decomp=False, uniform_u=False)
# area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk = check_box(grid, size, v_ival,\
#                                                                            bicentered_krawczyk_extension, eps, uniform_u=False)
# uni_plotter(area_boxes_bicentered_krawczyk, border_boxes_bicentered_krawczyk, u_lims, "Bicentered Krawczyk", size=2,
#             logger=bicentered_krawczyk_loger, ax=ax1, fig=fig1)
plot_area(ax1, a, b, left_v1, right_v1, left_v2, right_v2)
# plt.savefig('./fig/passive-rehabilitation-system-Bicentered _Krawczyk_'+str(N) + "_" + str(Nv) + "_" +
#             str(args.v1_0)  + "_" + str(args.v1_1) + "_" + str(args.v2_0) + "_" + str(args.v2_1) + '.png')
###
# fig1 = plt.figure(figsize=(8, 8))
# ax1 = fig1.add_subplot(1, 1, 1)
# bicentered_krawczyk_loger_infl = Logger(grid, size, v_ival, eps, bicentered_krawczyk_extension, decomp=False,
#                                         uniform_u=False, strategy="Inflation", dim=2, grid_v=grid_v)
# area_boxes_bicentered_krawczyk_infl, border_boxes_bicentered_krawczyk_infl = check_box(grid, size, v_ival,\
#                                                                            bicentered_krawczyk_extension, eps,
#                                                                              strategy="Inflaction", dim_v=2,
#                                                                              grid_v=grid_v, uniform_u=False)
# uni_plotter(area_boxes_bicentered_krawczyk_infl, border_boxes_bicentered_krawczyk_infl, u_lims, "Bicentered Krawczyk Inflaction", size=2,
#             logger=bicentered_krawczyk_loger_infl, ax=ax1, fig=fig1)
# plot_area(ax1, a, b, left_v1, right_v1, left_v2, right_v2)
# plt.savefig('./fig/passive-rehabilitation-system-Bicentered _Krawczyk_inflation_'+str(N) + "_" + str(Nv) + "_" +
#             str(args.v1_0)  + "_" + str(args.v1_1) + "_" + str(args.v2_0) + "_" + str(args.v2_1) + '.png')
###
# hansen_sengupta_extension = HansenSenguptaExtension(f_sym, v_sym, u_sym, coef=1, is_elementwise=False)
# hansen_sengupta_loger = Logger(grid, size, v_ival, eps, hansen_sengupta_extension, decomp=False)
# area_boxes_hansen_sengupta, border_boxes_hansen_sengupta = check_box(grid, size, v_ival,\
#                                                                            hansen_sengupta_extension, eps, log=False, decomposition=False)
# uni_plotter(area_boxes_hansen_sengupta, border_boxes_hansen_sengupta, u_upper, "Hansen-Sengupta", size=2, logger=hansen_sengupta_loger)
# plot_circles()
# plot_area(ax1, a, b, left_v1, right_v1, left_v2, right_v2)
plt.show()