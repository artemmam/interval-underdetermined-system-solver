import sympy as sym
import numpy as np
import interval as ival
import sys
from check_box import check_box, check_one_box
from ExtensionClass import ClassicalKrawczykExtension, BicenteredKrawczykExtension, HansenSenguptaExtension
from plotter_support_functions import uni_plotter, plot_one_box
import matplotlib.pyplot as plt
from LoggerClass import Logger
import argparse
from TestingExampleClass import Example


def plot_area(ax1, r1=3, r2=15, d=8):
    circle = plt.Circle((-0.5*d, 0), radius=r1, fc='y', fill=False)
    ax1.add_patch(circle)
    circle = plt.Circle((0.5*d, 0), radius=r1, fc='y', fill=False)
    ax1.add_patch(circle)
    circle = plt.Circle((-0.5*d, 0), radius=r2, fc='y', fill=False)
    ax1.add_patch(circle)
    circle = plt.Circle((0.5*d, 0), radius=r2, fc='y', fill=False)
    ax1.add_patch(circle)


def symbolic_2rpr_func(d = 8):
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
parser.add_argument('-r1', dest="r1", type=int)
parser.add_argument('-r2', dest="r2", type=int)
parser.add_argument('-d', dest="d", type=int)

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
r1 = args.r1
r2 = args.r2
d = args.d
f_sym, u_sym, v_sym = symbolic_2rpr_func(d)
v1 = ival.Interval([r1, r2])
v2 = ival.Interval([r1, r2])
v_ival = [v1, v2]
# u_lims = 6  # the width of the of the 2-dimensional square
# ux_lower, ux_upper, uy_lower, uy_upper = get_minmax_xy(a, b, left_v1, right_v1, left_v2, right_v2)
u_x = [-r2 - 1, r2 + 1]
u_y = [-r2 - 1, r2 + 1]
u_lims = [u_x, u_y]
u1 = ival.Interval([-r2 - 1, r2 + 1])
u2 = ival.Interval([-r2 - 1, r2 + 1])
u_ini = [u1, u2]
grid_u1 = np.linspace(u_x[0], u_x[1], N + 1)  # The vector to build size-dim. grid
grid_u2 = np.linspace(u_y[0], u_y[1], N + 1)  # The vector to build size-dim. grid
grid_u = [grid_u1, grid_u2]
u_dim = 2  # The dimension of uniform grid
eps = 1e-6
coef = 2
Nv = args.Nv
if args.Nv:
    grid_v1 = np.linspace(v1[0], v1[1], Nv + 1)
    grid_v2 = np.linspace(v2[0], v2[1], Nv + 1)
    grid_v = [grid_v1, grid_v2]
v_dim = 2
area_params = [r1, r2, d]
save_fig_params = [N, Nv, r1, r2, d, args.parallel]
bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
# f_num = bicentered_krawczyk_extension.lambdify_f()
# u_box = [ival.Interval([-9.6, -6.4]), ival.Interval([6.4, 9.6])]
# print(f_num([3, 3], u_box))
# sys.exit(1)
#####
# Bicentered_Krawczyk_Enlargment_V = Example(bicentered_krawczyk_extension, parallel=args.parallel, record_time=False, strategy="Enlargment")
#
# area_boxes, border_boxes = Bicentered_Krawczyk_Enlargment_V.check_box_branch(u_ini, v_ival, eps1=eps, eps2=8, grid_v=grid_v, v_dim=v_dim,
#                                            uniform_v=False)
# print("Enlargment V time bisection, ", Bicentered_Krawczyk_Enlargment_V.time)
# Bicentered_Krawczyk_Enlargment_V.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area,
#                                           area_params=area_params, save_fig=args.plotting, title = "Bicentered_Krawczyk_Enlargment_V_2RPR_branch")
# ######
# Bicentered_Krawczyk_Enlargment_V_bnb = Example(bicentered_krawczyk_extension, parallel=args.parallel, record_time=False, strategy="Enlargment")
#
# area_boxes, border_boxes = Bicentered_Krawczyk_Enlargment_V_bnb.check_box_branch(u_ini, v_ival, eps1=eps, eps2=8, eps3=1.7, mod="BNB")
# print("Enlargment BNB V time bisection, ", Bicentered_Krawczyk_Enlargment_V_bnb.time)
# Bicentered_Krawczyk_Enlargment_V_bnb.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area,
#                                           area_params=area_params, save_fig=args.plotting, title = "Bicentered_Krawczyk_Enlargment_V BNB_2RPR_branch")
######
# BicLoggerPlot = Logger(grid_u, u_dim, v_ival, eps, bicentered_krawczyk_extension, uniform_u=False)
Bicentered_Krawczyk_Default = Example(bicentered_krawczyk_extension, parallel=args.parallel, record_time=False, strategy="Default")
#
area_boxes, border_boxes = Bicentered_Krawczyk_Default.check_box_branch(u_ini, v_ival, eps1=eps, eps2=2)
print("Default V time bisection, ", Bicentered_Krawczyk_Default.time)
Bicentered_Krawczyk_Default.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area,
                                          area_params=area_params, save_fig=args.plotting, title = "Bicentered_Krawczyk_Default_2RPR branch")
#####
area_boxes, border_boxes = Bicentered_Krawczyk_Default.check_box(grid_u, u_dim, v_ival, eps, uniform_u=False)
print("Default V time basic, ", Bicentered_Krawczyk_Default.time)
Bicentered_Krawczyk_Default.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_area,
                                          area_params=area_params, save_fig=args.plotting, title = "Bicentered_Krawczyk_Default_2RPR basic")
if not args.parallel:
    plt.show()