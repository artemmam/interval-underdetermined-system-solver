import sympy as sym
import numpy as np
import interval as ival
import sys
from check_box import check_box, check_one_box
from ExtensionClass import ClassicalKrawczykExtension, BicenteredKrawczykExtension, HansenSenguptaExtension
from plotter_support_functions import uni_plotter
import matplotlib.pyplot as plt
from LoggerClass import Logger
import argparse
from TestingExampleClass import Example
import itertools as it

def plot_circles(ax, left = 3, right = 15, a = 5, l=6):
    circle = plt.Circle((left, -a/2), radius=l, fc='y', fill=False)
    ax.add_patch(circle)
    circle = plt.Circle((right, -a/2), radius=l, fc='y', fill=False)
    ax.add_patch(circle)
    circle = plt.Circle((left, a/2), radius=l, fc='y', fill=False)
    ax.add_patch(circle)
    circle = plt.Circle((right, a/2), radius=l, fc='y', fill=False)
    ax.add_patch(circle)


def symbolic_2dof_func(l1=6, l2=6, a=5):
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
    f = sym.Matrix([[-(u[0] - v[0])**2 - (u[1] - a/2)**2 + l1**2],
                    [-(u[0] - v[1])**2 - (u[1] + a/2)**2 + l2**2]])
    return f, u, v


parser = argparse.ArgumentParser(description="Angles in radians")
parser.add_argument('-Nu', dest="Nu", type=int)
parser.add_argument('-Nv', dest="Nv", type=int)
parser.add_argument('--parallel', dest="parallel", action='store_true')
parser.add_argument('--record_time', dest="record_time", action='store_true')
parser.add_argument('--plotting', dest="plotting", action='store_true')
parser.add_argument('-left', dest="left", type=int)
parser.add_argument('-right', dest="right", type=int)
parser.add_argument('-l', dest="l", type=int)

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

left = args.left
right = args.right
l = args.l
f_sym, u_sym, v_sym = symbolic_2dof_func(l1=l, l2=l)
v1 = ival.Interval([left, right])
v2 = ival.Interval([left, right])
v_ival = [v1, v2]
# borders = get_minmax_xy(f_sym, v_sym, u_sym, ang1_0 = left_v1, ang1_1 = right_v1,  ang2_0 = left_v2, ang2_1 = right_v2)
# print(borders)
# plot_circles(left, right)
# plt.show()
# sys.exit(1)
# u_lims = 6  # the width of the of the 2-dimensional square
# ux_lower, ux_upper, uy_lower, uy_upper = get_minmax_xy(a, b, left_v1, right_v1, left_v2, right_v2)
# u_l = -15
# u_u = 15
u_x = [v1.mid() - 2*l, v1.mid() + 2*l]
u_y = [-10, 10]
u_lims = [u_x, u_y]
# u1 = ival.Interval([u_l, u_u])
# u2 = ival.Interval([u_l, u_u])
# u_ini = [u1, u2]
grid_u1 = np.linspace(u_x[0], u_x[1], N + 1)  # The vector to build size-dim. grid
grid_u2 = np.linspace(u_y[0], u_y[1], N + 1)  # The vector to build size-dim. grid
grid_u = [grid_u1, grid_u2]
u_dim = 2  # The dimension of uniform grid
eps = 1e-6
# side = u1.width()/N
# eps_bnb = np.sqrt(2*side**2) - 0.1
coef = 2
Nv = args.Nv
grid_v1 = np.linspace(v1[0], v1[1], Nv + 1)
grid_v2 = np.linspace(v2[0], v2[1], Nv + 1)
grid_v = [grid_v1, grid_v2]
v_dim = 2
# print(grid_v1)
# area_params = [r1, r2, d]
# save_fig_params = [N, Nv, r1, r2, d, args.parallel]
save_fig_params = [N, Nv, left, right, args.parallel]
# print(Nv, N, v_ival)
area_params = [left, right, 5, l]
bicentered_krawczyk_extension = BicenteredKrawczykExtension(f_sym, v_sym, u_sym, coef=coef, is_elementwise=False)
#**********
# box = [ival.Interval([1.5, 2.25]), ival.Interval([0.0, 0.75])]
# temp = check_one_box(box, v_ival, bicentered_krawczyk_extension, eps, log=True, max_iter=9, decomposition=False, strategy = "Enlargement", grid_v = grid_v, dim_v= 2, uniform_v=False)
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
BicLoggerPlot = Logger(grid_u, u_dim, v_ival, eps, bicentered_krawczyk_extension, uniform_u=False,
                       strategy="Default", grid_v=grid_v, dim=v_dim,
                       )
#
# area_boxes, border_boxes = Bicentered_Krawczyk_Default.check_box_branch(u_ini, v_ival, eps1=eps, eps2=eps_bnb)
# print("Default V time bisection, ", Bicentered_Krawczyk_Default.time)
# Bicentered_Krawczyk_Default.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_circles, area_params=area_params, save_fig=args.plotting, title = "Bicentered_Krawczyk_Default_2-DOF_branch")
# ####
area_boxes, border_boxes = Bicentered_Krawczyk_Default.check_box(grid_u, u_dim, v_ival, eps, uniform_u=False)
# print("Default V time basic, ", Bicentered_Krawczyk_Default.time)
Bicentered_Krawczyk_Default.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_circles, area_params=area_params,
                                     save_fig=args.plotting, title = "Bicentered_Krawczyk_Default_2-DOF", logger=BicLoggerPlot)
#####
# Bicentered_Krawczyk_Enlargment_V = Example(bicentered_krawczyk_extension, parallel=args.parallel, record_time=False, strategy="Enlargement")
#
# area_boxes, border_boxes = Bicentered_Krawczyk_Enlargment_V.check_box(grid_u, u_dim, v_ival, eps=eps, grid_v=grid_v, v_dim=v_dim,
#                                            uniform_v=False)
# Bicentered_Krawczyk_Enlargment_V.plotting(area_boxes, border_boxes, u_lims, plot_area=plot_circles, area_params=area_params, save_fig=args.plotting, title = "Bicentered_Krawczyk_Enlargment_2-DOF")
# if rank == 0:
#     print("Enlargment V time, ", Bicentered_Krawczyk_Enlargment_V.time)
#     Bicentered_Krawczyk_Enlargment_V.plotting(area_boxes, border_boxes, u_lims, save_fig=args.plotting, title = "Bicentered_Krawczyk_Enlargment_simple_DexTar")
#     save_boxes("dextar_simple_inside_" + str(N) + "_" + str(Nv) + ".txt", area_boxes)
#     save_boxes("dextar_simple_border_" + str(N) + "_" + str(Nv) + ".txt", border_boxes)
if not args.parallel:
    plt.show()