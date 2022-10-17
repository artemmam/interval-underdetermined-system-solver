import sympy as sym
import numpy as np
import interval as ival
from sympy import sin, cos, nonlinsolve
import math
import sys
from check_box import check_box, check_one_box
from ExtensionClass import ClassicalKrawczykExtension, BicenteredKrawczykExtension, HansenSenguptaExtension
from plotter_support_functions import uni_plotter, plot_one_box
import matplotlib.pyplot as plt
# from LoggerClass import Logger
import argparse
from TestingExampleClass import Example
import itertools as it
# from LoggerClass import Logger


def save_boxes(name, boxes):
    with open("./boxes/" + name, "w+") as outfile:
        for box in boxes:
            for dim in box:
                for a in dim:
                    outfile.write(str(a) + " ")
            outfile.write("\n")

def distance(x1, y1, x2, y2):
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def angle(vector1, vector2):
    x1, y1 = vector1
    x2, y2 = vector2
    inner_product = x1*x2 + y1*y2
    len1 = math.hypot(x1, y1)
    len2 = math.hypot(x2, y2)
    return math.degrees(math.acos(inner_product/(len1*len2)))


def get_coordinates(al, beta, a=8, b=5, d=9):
    x = a*np.cos(al) + b*np.cos(al + beta)
    y = a*np.sin(al) + b*np.sin(al + beta)
    # [(u[0] - a * cos(v[0])) ** 2 + (u[1] - a * sin(v[0])) ** 2 - b ** 2],
    # [(u[0] - d - b * cos(v[1])) ** 2 + (u[1] - b * sin(v[1])) ** 2 - a ** 2]])
    return x, y

# def get_coordinates_first_bar(al, a=4):
#     x = a*np.cos(al)
#     y = a*np.sin(al)
#     return x, y


def plot_area(ax, params=[]):
    circle = plt.Circle((-4.5, 0), radius=3, fc='y', fill=False)
    circle1 = plt.Circle((4.5, 0), radius=3, fc='y', fill=False)
    ax.add_patch(circle)
    ax.add_patch(circle1)
    circle = plt.Circle((-4.5, 0), radius=13, fc='y', fill=False)
    circle1 = plt.Circle((4.5, 0), radius=13, fc='y', fill=False)
    ax.add_patch(circle)
    ax.add_patch(circle1)


def get_minmax_xy(f_sym, v_sym, u_sym, ang1_0 = None, ang1_1 = None,  ang2_0 = None, ang2_1 = None):
    A = list(it.product([ang1_0, ang1_1], [ang2_0, ang2_1]))
    dots = []
    print(A)
    for el in A:
        print("---")
        print(el)
        subs = ((v_sym[0], el[0]),
                (v_sym[1], el[1]))
        f_subs = f_sym.subs(subs)
        # print(f_subs)
        sol = nonlinsolve(f_subs, u_sym)
        for x in sol:
            dots.append(x)
    dots = np.array(dots)
    print(dots)
    return min(dots[::, 0]), max(dots[::, 0]), min(dots[::, 1]), max(dots[::, 1])


def symbolic_transformed_dextar_func(a, b, d=8):
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
    f = sym.Matrix([[(u[0] + d/2 - a*cos(v[0]))**2 + (u[1] - a*sin(v[0]))**2 - b**2],
                    [(u[0] - d/2 - a*cos(v[1]))**2 + (u[1] - a*sin(v[1]))**2 - b**2]])
    return f, u, v

# N = 30  # The number of boxes on uniform grid
# ##### 2-rpr
# f_sym, u_sym, v_sym = symbolic_2rpr_func()
# v1 = ival.Interval([3., 15.])
# v2 = ival.Interval([3., 15.])
# v_ival = [v1, v2]
# u_upper = 15  # the width of the of the 2-dimensional square
# grid = np.linspace(-u_upper, u_upper, N + 1)  # The vector to build size-dim. grid
# size = 2  # The dimension of uniform grid
# eps = 1e-3  # accuracy
# coef = 2  # Coefficient
# grid_v = np.linspace(v1[0], v1[1], 10 + 1)


parser = argparse.ArgumentParser(description="Angles in radians")
parser.add_argument('-Nu', dest="Nu", type=int)
parser.add_argument('-Nv', dest="Nv", type=int)
parser.add_argument('--parallel', dest="parallel", action='store_true')
parser.add_argument('--record_time', dest="record_time", action='store_true')
parser.add_argument('--plotting', dest="plotting", action='store_true')
parser.add_argument('-v1', dest="v1", type=str)
parser.add_argument('-v2', dest="v2", type=str)

args = parser.parse_args()
# print(args)
if args.parallel:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    world_size = comm.Get_size()
    rank = comm.Get_rank()
else:
    rank = 0
    world_size = 1
N = args.Nu  # The number of boxes on uniform grid
v1_0, v1_1 = np.fromstring(args.v1, dtype=int, sep=',')
v2_0, v2_1 = np.fromstring(args.v2, dtype=int, sep=',')
left_v1 = math.radians(v1_0)
right_v1 = math.radians(v1_1)
left_v2 = math.radians(v2_0)
right_v2 = math.radians(v2_1)
a, b, d = 8, 5, 9
f_sym, u_sym, v_sym = symbolic_transformed_dextar_func(a, b, d)
v1 = ival.Interval([left_v1, right_v1])
v2 = ival.Interval([left_v2, right_v2])
v_ival = [v1, v2]
# borders = get_minmax_xy(f_sym, v_sym, u_sym, ang1_0 = left_v1, ang1_1 = right_v1,  ang2_0 = left_v2, ang2_1 = right_v2)
# print(borders)
# sys.exit(1)
# u_lims = 6  # the width of the of the 2-dimensional square
# ux_lower, ux_upper, uy_lower, uy_upper = get_minmax_xy(a, b, left_v1, right_v1, left_v2, right_v2)
u_l = -20
u_u = 20
u_x = [u_l, u_u]
u_y = [u_l, u_u]
u_lims = [u_x, u_y]
u1 = ival.Interval([u_l, u_u])
u2 = ival.Interval([u_l, u_u])
u_ini = [u1, u2]
grid_u1 = np.linspace(u_x[0], u_x[1], N + 1)  # The vector to build size-dim. grid
grid_u2 = np.linspace(u_y[0], u_y[1], N + 1)  # The vector to build size-dim. grid
grid_u = [grid_u1, grid_u2]
u_dim = 2  # The dimension of uniform grid
eps = 1e-6
side = u1.width()/N
eps_bnb = np.sqrt(2*side**2) - 0.1
coef = 2
Nv = args.Nv
grid_v1 = np.linspace(v1[0], v1[1], Nv + 1)
grid_v2 = np.linspace(v2[0], v2[1], Nv + 1)
grid_v = [grid_v1, grid_v2]
v_dim = 2
# area_params = [r1, r2, d]
# save_fig_params = [N, Nv, r1, r2, d, args.parallel]
save_fig_params = [N, Nv, left_v1, right_v1, left_v2, right_v2, a, b, d, args.parallel]

bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
#**********
# box = [ival.Interval([8, 9]), ival.Interval([8, 9])]
# temp = check_one_box(box, v_ival, bicentered_krawczyk_extension, eps, log=True, max_iter=9, decomposition=False, strategy = "Enlargement", grid = grid_v, dim = 2, uniform=False)
# print(temp)
# sys.exit(1)
#**********
#####
# Bicentered_Krawczyk_Enlargment_V = Example(bicentered_krawczyk_extension, parallel=args.parallel, record_time=False, strategy="Enlargment")
#
# area_boxes, border_boxes = Bicentered_Krawczyk_Enlargment_V.check_box_branch(u_ini, v_ival, eps1=eps, eps2=8, grid_v=grid_v, v_dim=v_dim,
#                                            uniform_v=False)
# print("Enlargment V time bisection, ", Bicentered_Krawczyk_Enlargment_V.time)
# Bicentered_Krawczyk_Enlargment_V.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area,
#                                           area_params=area_params, save_fig=args.plotting, title = "Bicentered_Krawczyk_Enlargment_V_2RPR_branch")
######
# Bicentered_Krawczyk_Enlargment_V_bnb = Example(bicentered_krawczyk_extension, parallel=args.parallel, record_time=False, strategy="Enlargment")
#
# area_boxes, border_boxes = Bicentered_Krawczyk_Enlargment_V_bnb.check_box_branch(u_ini, v_ival, eps1=eps, eps2=8, eps3=1.7, mod="BNB")
# print("Enlargment BNB V time bisection, ", Bicentered_Krawczyk_Enlargment_V_bnb.time)
# Bicentered_Krawczyk_Enlargment_V_bnb.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area,
#                                           area_params=area_params, save_fig=args.plotting, title = "Bicentered_Krawczyk_Enlargment_V BNB_2RPR_branch")
######
Bicentered_Krawczyk_Default = Example(bicentered_krawczyk_extension, parallel=args.parallel, record_time=False, strategy="Default")
# BicLoggerPlot = Logger(grid_u, u_dim, v_ival, eps, bicentered_krawczyk_extension, uniform_u=False,
#                        strategy="Default", grid_v=grid_v, dim=v_dim, decomp=True
#                        )
#
area_boxes, border_boxes = Bicentered_Krawczyk_Default.check_box_branch(u_ini, v_ival, eps_krawczyk=eps, eps_bnb=eps_bnb)
# print("Default V time bisection, ", Bicentered_Krawczyk_Default.time)
Bicentered_Krawczyk_Default.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area, area_params=save_fig_params,
                                     save_fig=args.plotting, title="Bicentered_Krawczyk_DexTar branch")

# ####

# area_boxes, border_boxes = Bicentered_Krawczyk_Default.check_box_branch(u_ini, v_ival, eps_krawczyk=eps, eps_bnb=eps_bnb ,decomposition=True)
# print("BnB decomposition time, ", Bicentered_Krawczyk_Default.time)
# Bicentered_Krawczyk_Default.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area, area_params=[],
#                                      save_fig=args.plotting, title="Bicentered_Krawczyk_Decomposition_DexTar branch")
# #####
# Bicentered_Krawczyk_Enlargment_V = Example(bicentered_krawczyk_extension, parallel=args.parallel, record_time=False, strategy="Enlargement")
# #
# area_boxes, border_boxes = Bicentered_Krawczyk_Enlargment_V.check_box_branch(u_ini, v_ival, eps_krawczyk=eps, eps_bnb=eps_bnb, grid_v=grid_v, v_dim=v_dim,
#                                            uniform_v=False)
# print("BnB enlargement time, ", Bicentered_Krawczyk_Enlargment_V.time)
# Bicentered_Krawczyk_Enlargment_V.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area, area_params=[], save_fig=args.plotting, title = "Bicentered_Krawczyk_Enlargment DexTar branch")
if rank == 0:
    print("Num procs", world_size)
    print("BnB time, ", Bicentered_Krawczyk_Default.time)
    # print("Enlargment V time, ", Bicentered_Krawczyk_Enlargment_V.time)
    # Bicentered_Krawczyk_Enlargment_V.plotting(area_boxes, border_boxes, u_lims, save_fig=args.plotting, title = "Bicentered_Krawczyk_Enlargment_simple_DexTar")
#     save_boxes("dextar_simple_inside_" + str(N) + "_" + str(Nv) + ".txt", area_boxes)
#     save_boxes("dextar_simple_border_" + str(N) + "_" + str(Nv) + ".txt", border_boxes)
if not args.parallel:
    plt.show()