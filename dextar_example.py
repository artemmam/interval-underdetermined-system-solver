import sympy as sym
import numpy as np
import interval as ival
import sys
from check_box import check_box, check_one_box, check_box_parallel
from ExtensionClass import ClassicalKrawczykExtension, BicenteredKrawczykExtension, HansenSenguptaExtension
from plotter_support_functions import uni_plotter, plot_one_box
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
# from LoggerClass import Logger
from sympy import sin, cos
import math
import itertools as it
import argparse
import matplotlib
from TestingExampleClass import Example

from timeit import default_timer as timer


# matplotlib.use('Agg')
def write_time(file, size, n, time):
    f = open(file + ".txt", "a+")
    f.write(str(size) + ", " + str(n) + ", " + str(time) + "\n")
    f.close()


def write_time_per_proc(file, rank, time):
    f = open(file + ".txt", "a+")
    f.write(str(rank) + ": " + str(time) + "\n")
    f.close()


def distance(x1, y1, x2, y2):
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def angle(vector1, vector2):
    x1, y1 = vector1
    x2, y2 = vector2
    inner_product = x1 * x2 + y1 * y2
    len1 = math.hypot(x1, y1)
    len2 = math.hypot(x2, y2)
    return math.degrees(math.acos(inner_product / (len1 * len2)))


def get_coordinates(al, beta, a=4, b=2):
    x = a * np.cos(al) + b * np.cos(al + beta)
    y = a * np.sin(al) + b * np.sin(al + beta)
    return x, y


def get_coordinates_first_bar(al, a=4):
    x = a * np.cos(al)
    y = a * np.sin(al)
    return x, y


def plot_area(ax1, a=4, b=2, ang1_0=None, ang1_1=None, ang2_0=None, ang2_1=None):
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
    print(dots)
    ang1 = angle(dots[0], dots[2])
    ang2 = angle(dots[1], dots[3])
    r1 = distance(0, 0, *dots[0])
    r2 = distance(0, 0, *dots[1])
    ang3 = angle(dots[0], (r1, 0))
    ang4 = angle(dots[1], (r2, 0))
    print(ang3, ang3 + ang1)
    arc1 = Arc((0, 0), 2 * r1, 2 * r1, 0, ang3, ang3 + ang1, color='red', lw=1)
    ax1.add_patch(arc1)
    print(ang4, ang4 + ang2)
    arc2 = Arc((0, 0), 2 * r2, 2 * r2, 0, ang4, ang4 + ang2, color='red', lw=1)
    ax1.add_patch(arc2)
    ax1.scatter(dots[::, 0], dots[::, 1], color="red")
    dots_bar1 = []
    for el in [ang1_0, ang1_1]:
        dots_bar1.append(get_coordinates_first_bar(el, a))
    r1_bar1 = distance(*dots_bar1[0], *dots[0])
    arc3 = Arc(dots_bar1[0], 2 * r1_bar1, 2 * r1_bar1, 0, math.degrees(ang1_0 + ang2_0), math.degrees(ang2_1 + ang1_0),
               color='red', lw=1)
    ax1.add_patch(arc3)
    r2_bar1 = distance(*dots_bar1[1], *dots[2])
    arc4 = Arc(dots_bar1[1], 2 * r2_bar1, 2 * r2_bar1, 0, math.degrees(ang1_1 + ang2_0), math.degrees(ang1_1 + ang2_1),
               color='red', lw=1)
    ax1.add_patch(arc4)


def get_minmax_xy(a=4, b=2, ang1_0=None, ang1_1=None, ang2_0=None, ang2_1=None):
    A = list(it.product([ang1_0, ang1_1], [ang2_0, ang2_1]))
    dots = []
    for el in A:
        XY = get_coordinates(*el, a, b)
        dots.append(XY)
    dots = np.array(dots)
    return min(dots[::, 0]), max(dots[::, 0]), min(dots[::, 1]), max(dots[::, 1])


def symbolic_dextar_system_func(l_a=7.2, l_b=2, l_c = 2, l_d = 7.2, d=6):
    """
    Creating symbol variables for circle eq. system
    :return: symbolic eq. system,
            symbolic U (fixed boxes),
            symbolic V (checking boxes),
            symbolic C
    """
    v = sym.symbols('v1, v2, v3, v4')
    u = sym.symbols('u1, u2')
    f = sym.Matrix([
                    [u[0] - l_b * cos(v[2]) - l_a * cos(v[0]) - d/2],
                    [u[1] - l_d * cos(v[1]) - l_c * cos(v[3]) + d/2],
                    [u[0] - l_a * sin(v[0]) - l_b * sin(v[2])],
                    [u[1] - l_d * sin(v[1]) - l_c * sin(v[3])]
                    ])
    return f, u, v


# a = 4.5
# b = 4.5
parser = argparse.ArgumentParser(description="Angles in radians")
parser.add_argument('-Nu', dest="Nu", type=int)
parser.add_argument('-Nv', dest="Nv", type=int)
parser.add_argument('--parallel', dest="parallel", action='store_true')
parser.add_argument('--record_time', dest="record_time", action='store_true')
parser.add_argument('--plotting', dest="plotting", action='store_true')
parser.add_argument('-v1', dest="v1", type=str)
parser.add_argument('-v2', dest="v2", type=str)
parser.add_argument('-v3', dest="v3", type=str)
parser.add_argument('-v4', dest="v4", type=str)

args = parser.parse_args()
print(args)
# print(args)
if args.parallel:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    world_size = comm.Get_size()
    rank = comm.Get_rank()
else:
    rank = 0
    world_size = 1
path = "./log/rehub_sys_bicentered_krawczyk_enlargement_time_procs" + str(rank)
# path = 'rehub_sys_bicentered_krawczyk_enlargement_time_procs'
# if rank == 0:
# open('rehub_sys_bicentered_krawczyk_enlargement_time_procs.txt', 'w').close()
if args.parallel:
    comm.Barrier()
N = args.Nu  # The number of boxes on uniform grid
v1_0, v1_1 = np.fromstring(args.v1, dtype=int, sep=',')
v2_0, v2_1 = np.fromstring(args.v2, dtype=int, sep=',')
v3_0, v3_1 = np.fromstring(args.v3, dtype=int, sep=',')
v4_0, v4_1 = np.fromstring(args.v4, dtype=int, sep=',')
left_v1 = math.radians(v1_0)
right_v1 = math.radians(v1_1)
left_v2 = math.radians(v2_0)
right_v2 = math.radians(v2_1)
left_v3 = math.radians(v3_0)
right_v3 = math.radians(v3_1)
left_v4 = math.radians(v4_0)
right_v4 = math.radians(v4_1)
f_sym, u_sym, v_sym = symbolic_dextar_system_func()
v1 = ival.Interval([left_v1, right_v1])
v2 = ival.Interval([left_v2, right_v2])
v3 = ival.Interval([left_v1, right_v1])
v4 = ival.Interval([left_v2, right_v2])
v_ival = [v1, v2, v3, v4]
# u_lims = 6  # the width of the of the 2-dimensional square
ux_lower, ux_upper, uy_lower, uy_upper = -12.2, 12.2, -9.2, 9.2#get_minmax_xy(a, b, left_v1, right_v1, left_v2, right_v2)
u_x = [ux_lower - 1, ux_upper + 1]
u_y = [uy_lower - 1, uy_upper + 1]
u_lims = [u_x, u_y]
grid_u1 = np.linspace(u_x[0], u_x[1], N + 1)  # The vector to build size-dim. grid
grid_u2 = np.linspace(u_y[0], u_y[1], N + 1)  # The vector to build size-dim. grid
grid_u = [grid_u1, grid_u2]
u_dim = 2  # The dimension of uniform grid
eps = 1e-6
coef = 2
Nv = args.Nv
grid_v1 = np.linspace(v1[0], v1[1], Nv + 1)
grid_v2 = np.linspace(v2[0], v2[1], Nv + 1)
grid_v3 = np.linspace(v3[0], v3[1], Nv + 1)
grid_v4 = np.linspace(v4[0], v4[1], Nv + 1)
grid_v = [grid_v1, grid_v2, grid_v3, grid_v4]
v_dim = 4
# if rank == 0:
# start = timer()
area_params = [left_v1, right_v1, left_v2, right_v2]

save_fig_params = [N, Nv, left_v1, right_v1, left_v2, right_v2, args.parallel]
bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
Bicentered_Krawczyk_Enlargment_V = Example(bicentered_krawczyk_extension, path=path, parallel=args.parallel,
                                           record_time=False, strategy="Enlargment", log=False)
area_boxes, border_boxes = Bicentered_Krawczyk_Enlargment_V.check_box(grid_u, u_dim, v_ival, eps, grid_v, v_dim,
                                                                      uniform_u=False, uniform_v=False)
if args.record_time:
    write_time_per_proc(path, rank, Bicentered_Krawczyk_Enlargment_V.time)
Bicentered_Krawczyk_Enlargment_V.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area,
                                          area_params=area_params, save_fig=args.plotting,
                                          title="Bicentered_Krawczyk_Enlargment_V_DexTAR",
                                          save_fig_params=save_fig_params)

Bicentered_Krawczyk_Default = Example(bicentered_krawczyk_extension, parallel=args.parallel, record_time=False, strategy="Default")
area_boxes, border_boxes = Bicentered_Krawczyk_Default.check_box(grid_u, u_dim, v_ival, eps, uniform_u=False)
# if args.record_time:
# write_time_per_proc("./rehub_sys_bicentered_krawczyk_time_procs", rank, Bicentered_Krawczyk_Default.time)
Bicentered_Krawczyk_Default.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area,
area_params=area_params, save_fig=args.plotting, title = "Bicentered_Krawczyk_Default_rehab_system", save_fig_params=save_fig_params)
# Bicentered_Krawczyk_Default = Example(bicentered_krawczyk_extension, parallel=args.parallel, record_time=False, strategy="Default")
#
# area_boxes, border_boxes = Bicentered_Krawczyk_Default.check_box_branch(u_ini, v_ival, eps1=eps, eps2=2)
# print("Default V time bisection, ", Bicentered_Krawczyk_Default.time)
# Bicentered_Krawczyk_Default.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area,
#                                           area_params=area_params, save_fig=args.plotting, title = "Bicentered_Krawczyk_Default_2RPR branch")
# area_boxes, border_boxes = Bicentered_Krawczyk_Default.check_box(grid_u, u_dim, v_ival, eps, uniform_u=False)
# print("Default V time basic, ", Bicentered_Krawczyk_Default.time)
# Bicentered_Krawczyk_Default.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area,
#                                           area_params=area_params, save_fig=args.plotting, title = "Bicentered_Krawczyk_Default_2RPR basic")
if not args.parallel:
    plt.show()