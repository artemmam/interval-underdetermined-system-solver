import interval as ival
import numpy as np
import itertools as it
from timeit import default_timer as timer
from symbolic_support_functions import lambdify_f
# from pathos.multiprocessing import ProcessingPool as Pool
# import dill  # the code below will fail without this line
# import pickle
import sympy as sym
import copy


def write_time_per_proc(file, rank, time, message = ""):
    f = open(file + ".txt", "a+")
    f.write(str(rank) + ": " + str(time) + " | " + message +   "\n")
    f.close()


def diam(A):
    s = 0
    for v in A:
        s += v.width()**2
    return np.sqrt(s)

def unpack_array(a):
    b = []
    for proc in a:
        if len(proc)!=0:
            for box in proc:
                b.append(box)
    return b

def make_boxes_list(grid, dim, uniform=True):
    """
    Make list of boxes in dim dimension from vector grid
    :param grid: vector on which grid is constructed
    :param dim:  the dimensional of grid
    :return: the list of boxes in dim
    """
    if uniform == True:
        grid_size = len(grid) - 1
        grid_intervals = []
        grid_variants = []
        for i in range(grid_size):
            grid_intervals.append(ival.Interval([grid[i], grid[i + 1]]))
        for i in range(dim):
            grid_variants.append(grid_intervals)
        grid_n_size = list(it.product(*grid_variants))
    else:
        grid_variants = []
        grid_numbers = np.shape(grid)[0]
        for i in range(grid_numbers):
            grid_intervals = []
            one_grid = grid[i]
            grid_size = len(one_grid) - 1
            for j in range(grid_size):
                grid_intervals.append(ival.Interval([one_grid[j], one_grid[j + 1]]))
            grid_variants.append(grid_intervals)
        grid_n_size = list(it.product(*grid_variants))
    return grid_n_size


def reccur_func_mi(box, v_init, eps, extension, max_iter=10, log=False, decomposition=False, eps_decomp=np.pi, depth=12):
    v_iter = v_init.copy()
    n = len(v_init)
    k = 0
    f_num = extension.lambdify_f()
    if depth == 0:
        return "border"
    if log:
        print("*****")
        print("Natural :", f_num(v_iter, box).reshape(-1))
    check_1d = []
    for i in range(len(v_iter)):
        f_num_i = extension.lambdify_f(idx=i)
        a = f_num_i(v_iter[i][0], box)
        b = f_num_i(v_iter[i][1], box)
        if a[1] <= 0 <= b[0] or a[0] >= 0 >= b[1]:
            check_1d.append("R")
        elif f_num_i(v_iter[i], box)[0] > 0 or f_num_i(v_iter[i], box)[1] < 0:
            check_1d.append('N')
        else:
            check_1d.append('U')
    rvl_R = [x == 'R' for x in check_1d]
    rvl_N = [x == 'N' for x in check_1d]
    # print(rvl_R)
    if np.all(rvl_R):
        return "inside"
    elif any(rvl_N):
        return 'outside'
    v_ext = extension.calculate_extension(box, v_iter, log=log).reshape(-1)
    if log:
        print("Number of iteration =", k)
        print("box", box)
        print("Old V = ", v_iter)
        print("New V = ", v_ext)
    check_ext_intervals = [type(x) == ival.interval.ExtendedInterval for x in v_ext]
    if not(np.all(check_ext_intervals)):
        qil = []
        idx = np.arange(2)[check_ext_intervals]
        for j in idx:
            for i in range(n):
                c = v_iter[i].mid()
                q = c - v_ext[j][i]
                print(q)
                if q[0] <= v_iter[i][0] and q[1] >= v_iter[i][1]:
                    qil.append([ival.Interval([v_iter[i][0], c]), ival.Interval([c, v_iter[i][1]])])
                else:
                    if log:
                        print("xi changed from ", v_iter[i], " to ", q)
                    qil.append(q)
    else:
        qil = []
        for i in range(n):
            qil.append([v_ext[0][i], v_ext[1][i]])
    if log:
        print("qil", qil)
    rv = "outside"
    for z in qil:
        new_v = []
        for i in range(n):
            new_v.append(z[i].intersec(v_iter[i]))
        if log:
            print("Go ", new_v)
        res = reccur_func_mi(box, new_v, eps, extension, max_iter=10, log=log, decomposition=False, eps_decomp=np.pi, depth=depth-1)
        if res == "inside":
            rv = "inside"
            break
        elif res == "border":
            rv = "border"
    return rv

def reccur_func(box, v_init, eps, extension, max_iter=10, log=False, decomposition=False, eps_decomp=np.pi):
    v_iter = v_init.copy()
    n = len(v_init)
    k = 0
    f_num = extension.lambdify_f()
    while True:
        if log:
            print("*****")
            print("Natural :", f_num(v_iter, box).reshape(-1))
        # check_1d = []
        # for i in range(len(v_iter)):
        #     f_num_i = extension.lambdify_f(idx=i)
        #     a = f_num_i(v_iter[i][0], box)
        #     b = f_num_i(v_iter[i][1], box)
        #     if a[1] <= 0 <= b[0] or a[0] >= 0 >= b[1]:
        #         check_1d.append("R")
        #     elif f_num_i(v_iter[i], box)[0] > 0 or f_num_i(v_iter[i], box)[1] < 0:
        #         check_1d.append('N')
        #     else:
        #         check_1d.append('U')
        # rvl_R = [x == 'R' for x in check_1d]
        # rvl_N = [x == 'N' for x in check_1d]
        # # print(rvl_R)
        # if np.all(rvl_R):
        #     return "inside"
        # elif any(rvl_N):
        #     return 'outside'
        if k>5:
            for i in range(n):
                v_iter[i].scale(1.1)
                v_iter[i] = v_iter[i].intersec(v_init[i])
        check = True
        v_ext = extension.calculate_extension(box, v_iter, log=log).reshape(-1)
        if log:
            print("Number of iteration =", k)
            print("box", box)
            print("Old V = ", v_iter)
            print("New V = ", v_ext)
        for i in range(n):
            if not (v_ext[i].isIn(v_iter[i])):
                check = False
                break
        if check:
            return "inside"
        for i in range(n):
            if v_iter[i].isNoIntersec(v_ext[i]):
                return "outside"
            else:
                v_iter[i] = v_iter[i].intersec(v_ext[i])
        if k > max_iter:
            if decomposition:
                if log:
                    print("Bisection")
                    print(diam(v_iter))
                if diam(v_iter) > eps_decomp:
                    v_l, v_r = separate_box(v_iter)
                    if log:
                        print("Left", v_l)
                    res_l = reccur_func(box, v_l, eps, extension, max_iter=10, log=log, decomposition=decomposition, eps_decomp=eps_decomp)
                    if log:
                        print("Right", v_r)
                    res_r = reccur_func(box, v_r, eps, extension, max_iter=10, log=log, decomposition=decomposition, eps_decomp=eps_decomp)
                    if res_l == "inside" or res_r == "inside":
                        return "inside"
            return "border"
            # v_prev = v_iter.copy()
        k += 1


def reccur_func_enlarge(box, v_init, v_ival, eps, extension, max_iter=10, log=False, decomposition=False, eps_decomp=np.pi):
    v_iter = v_init.copy()
    n = len(v_init)
    # print("box", box)
    k = 0
    f_num = extension.lambdify_f()
    # if log:
    #     f_centered = extension.centered_form()

    while True:
        if log:
            print("*****")
            # print(v_iter, box)
            print("box", box)
            print("V_ival", v_ival)
            print("Old V = ", v_iter)
            print("Natural :", f_num(v_iter, box).reshape(-1))
        # for nat_ext in f_num(v_iter, box).reshape(-1):
        #     if not ival.Interval([0, 0]).isInIns(nat_ext):
        #         return "outside"
        # check_1d = []
        # for i in range(len(v_iter)):
        #     left = f_num[i](v_iter[i][0])
        #     right = f_num[i](v_iter[i][1])
        #     if right <= 0 <= left or left >= 0 >= right:
        #         check_1d.append("R")
        #     elif left > 0 or right < 0:
        #         check_1d.append('N')
        # rvl_R = [x == 'R' for x in check_1d]
        # rvl_N = [x == 'N' for x in check_1d]
        # if np.all(rvl_R):
        #     return "inside"
        # elif any(rvl_N):
        #     return 'outside'
        if k > 5:
            for i in range(n):
                v_iter[i].scale(1.1)
                v_iter[i] = v_iter[i].intersec(v_ival[i])
        check = True
        v_ext = extension.calculate_extension(box, v_iter, log=log).reshape(-1)
        if log:
            print("Number of iteration =", k)
            print("New V = ", v_ext)
        for i in range(n):
            if not (v_ext[i].isIn(v_iter[i])):
                check = False
                break
        if check:
            return "inside"
        for i in range(n):
            if v_ext[i].isNoIntersec(v_ival[i]):
                return "outside"
            else:
                v_iter[i] = v_ext[i].intersec(v_ival[i])
        if k > max_iter:
            if decomposition:
                if log:
                    print("Bisection")
                if diam(v_init) > eps_decomp:
                    v_l, v_r = separate_box(v_init)
                    if log:
                        print("Left", v_l)
                    res_l = reccur_func_enlarge(box, v_l, v_ival, eps, extension, max_iter=10, log=log,
                                        decomposition=decomposition, eps_decomp=eps_decomp)
                    if log:
                        print("Right", v_r)
                    res_r = reccur_func_enlarge(box, v_r, v_ival, eps, extension, max_iter=10, log=log,
                                        decomposition=decomposition, eps_decomp=eps_decomp)
                    if res_l == "inside" or res_r == "inside":
                        return "inside"
                    elif res_l == "outside" and res_r == "outside":
                        return "outside"
            return "border"
        k += 1


def reccur_func_elementwise(box, v_init, eps, extension, max_iter=10, log=False, decomposition=False, level=0):
    v_iter = v_init.copy()
    n = len(v_init)
    v_prev = v_iter.copy()
    k = 0
    f_num = extension.lambdify_f()
    while True:
        if log:
            print("Natural :", f_num(v_iter, box).reshape(-1))
        for nat_ext in f_num(v_iter, box).reshape(-1):
            if not ival.Interval([0, 0]).isInIns(nat_ext):
                return "outside"
        check = True
        v_ext_funcs, param, c = extension.calculate_extension(box, v_iter)
        c = np.array(c).reshape(n, n)
        if log:
            print("*****")
            print("Number of iteration =", k)
            print("box", box)
            print("Old V = ", v_iter)
            print(param, c)
            print("New V = ")
        for i in range(n):
            try:
                v_ext = v_ext_funcs[i](v_iter, c[i], param).reshape(-1)[0]
            except:
                v_ext = v_ext_funcs[i](v_iter, c[i], param)
            if log:
                print(v_ext)
            if not (v_ext.isIn(v_iter[i])):
                check = False
            if v_iter[i].isNoIntersec(v_ext):
                return "outside"
            else:
                v_iter[i] = v_iter[i].intersec(v_ext)
        if check:
            return "inside"
        if abs(diam(v_iter) - diam(v_prev)) / (0.5 * abs(diam(v_iter) + diam(v_prev))) < eps or k > max_iter:
            return "border"
        v_prev = v_iter.copy()
        k += 1


def check_box_parallel(grid, dim, v_ival, extension, eps, log=False, max_iter=10, decomposition=False, eps_decomp=None, strategy = "Default",
              grid_v = None, dim_v=None, uniform_v = True, uniform_u = True, path = "", enlargement=False, boxes=[]):
    """
    Function for checking boxes on dim-dimensional uniform grid with checker method
    :param grid: 1-d grid
    :param dim: number of the dimensions
    :param v_ival: vector of not fixed interval variables
    :param extension: Extension calculator-class object, contains param info and calculate interval extension
    :param eps: error
    :param log: turn on log info printing
    :return: list of inside boxes, list of border boxes
    """
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    start = timer()
    area_boxes = []
    border_boxes = []
    if len(boxes) == 0:
        grid = np.array(grid)
        grid_size = len(grid) - 1
        all_boxes = make_boxes_list(grid, dim, uniform_u)
    else:
        all_boxes = boxes
    for i, box in enumerate(all_boxes):
        if i % size == rank:
            # print("rank = ", rank, "; ", i, "/", len(all_boxes) - 1)
            # print(i, "/", len(all_boxes) - 1)
            if extension.is_elementwise:
                temp = reccur_func_elementwise(all_boxes[i], v_ival, eps, extension, max_iter, log=log,
                                               decomposition=decomposition)
            else:
                temp = check_1d(box, v_ival, extension, False)
                if temp != "inside":
                    temp = check_box_posypkin(extension.f, extension.u, extension.v, v_ival, box)
                    if temp == "border":
                        if strategy == "Default":
                            if enlargement:
                                v = v_ival
                                temp = reccur_func_enlarge(all_boxes[i], v, v_ival, eps, extension, max_iter, log=log,
                                                           decomposition=decomposition, eps_decomp=eps_decomp)
                            else:
                                temp = reccur_func(all_boxes[i], v_ival, eps, extension, max_iter, log=log,
                                                   decomposition=decomposition, eps_decomp=eps_decomp)
                        else:
                            temp = "outside"
                            grid_v = np.array(grid_v)
                            v_boxes = make_boxes_list(grid_v, dim_v, uniform_v)
                            temp_list = []
                            for v in v_boxes:
                                if len(v) == 1:
                                    v = [v[0]]
                                else:
                                    v = np.array(v)
                                if enlargement:
                                    temp_loc = reccur_func_enlarge(all_boxes[i], v, v_ival, eps, extension, max_iter, log=log,
                                                                   decomposition=decomposition)
                                else:
                                    temp_loc = reccur_func(all_boxes[i], v, eps, extension, max_iter, log=log,
                                                           decomposition=decomposition, eps_decomp=eps_decomp)
                                if temp_loc == "inside":
                                    temp = "inside"
                                    break
                                else:
                                    temp_list.append(temp_loc)
                            if temp != "inside":
                                check = [True if temp_list[i] == "border" else False for i in range(len(temp_list))]
                                if np.any(check):
                                    temp = "border"
            if temp == 'inside':
                area_boxes.append(all_boxes[i])
            elif temp == 'border':
                border_boxes.append(all_boxes[i])
    # end = timer()
    # write_time_per_proc(path, rank, end - start, "total = " + str(len(all_boxes)/size) + "; area = " + str(len(area_boxes)) +  "; border = " + str(len(border_boxes)))
    # start = timer()
    area_boxes_g = comm.gather(area_boxes, root=0)
    border_boxes_g = comm.gather(border_boxes, root=0)
    # end = timer()
    # write_time_per_proc(path, rank, end - start, "GATHER")
    area_boxes_g_unpacking = []
    border_boxes_g_unpacking = []
    if rank == 0:
        area_boxes_g_unpacking = unpack_array(area_boxes_g)
        border_boxes_g_unpacking = unpack_array(border_boxes_g)
    return area_boxes_g_unpacking, border_boxes_g_unpacking

def check_1d(box, v, extension, log, eps_1d = np.pi/8):
    v_ival = copy.deepcopy(v)
    check = np.zeros(len(v_ival))
    for i in range(len(v_ival)):
        if log:
            print(i)
        f_num = lambdify_f(extension.f[i], extension.u, extension.v[i])
        if log:
            print("Natural 1d", f_num(box, v_ival[i]))
        # print(v_ival[i])
        # print(f_num(box, v_ival[i]))
        if ival.Interval([0, 0]).isIn(f_num(box, v_ival[i])):
            while v_ival[i][0] != v_ival[i][1]:
                if log:
                    # print(f_num(box, v_ival[i]))
                    print("v_ival", v_ival[i])
                    print("1 var", f_num(box, v_ival[i][0])[1], f_num(box, v_ival[i][1])[0])
                    print((f_num(box, v_ival[i][0])[1]<=0 and f_num(box, v_ival[i][1])[0]>=0))
                    print("2 var", f_num(box, v_ival[i][0])[0], f_num(box, v_ival[i][1])[1])
                    print((f_num(box, v_ival[i][0])[0]>=0 and f_num(box, v_ival[i][1])[1]<=0))
                if ((f_num(box, v_ival[i][0])[1] <= 0 and f_num(box, v_ival[i][1])[0] >= 0) or (
                        f_num(box, v_ival[i][0])[0] >= 0 and f_num(box, v_ival[i][1])[1] <= 0)):
                    check[i] = True
                    break
                elif f_num(box, v_ival[i][0])[1] <= 0 or f_num(box, v_ival[i][0])[0] >= 0:
                    v_ival[i][1] = v_ival[i][1] - eps_1d
                elif f_num(box, v_ival[i][1])[0] >= 0 or f_num(box, v_ival[i][1])[1] <= 0:
                    v_ival[i][0] = v_ival[i][0] + eps_1d
                else:
                    v_ival[i][0] = v_ival[i][0] + eps_1d
        else:
            return "outside"
    if np.all(check):
        # print(box)
        #         print("inside")
        return "inside"
    else:
        return "not inside"


#     else:
#         print("not inside")
def check_roots(xi, f, df, depth, log = False):
    if log:
        print("%%%iter = ", 12 - depth, "%%%")
    if log:
        print("In: xi = ", xi)
    if depth == 0:
#         print("RECURSION failed for ", xi)
        rv = 'U'
    else:
        fi = f(xi)
        fia = f(ival.Interval([xi[0], xi[0]]))
        fib = f(ival.Interval([xi[1], xi[1]]))
        dfi = df(xi)
        if log:
            print("fi: ", fi, fia, fib)
        if fia[1] <= 0 <= fib[0] or fia[0] >= 0 >= fib[1]:
            rv = 'R'
        elif fi[0] > 0 or fi[1] < 0:
            rv = 'N'
        else:
            c = xi.mid()
            ci = ival.Interval([c,c])
            fci = f(ci)
            qi = fci / dfi
            if log:
                print("qi", qi)
            if type(qi) is ival.ExtendedInterval:
                qil = [ci - qi[0], ci - qi[1]]
            else:
                q = ci - qi
                if q[0] <= xi[0] and q[1] >= xi[1]:
                    qil = [ival.Interval([xi[0], c]), ival.Interval([c, xi[1]])]
                else:
                    if log:
                        print("xi changed from ", xi, " to ", q)
                    qil = [q]
            if log:
                    print("qil", qil)
            rv = 'N'
            for z in qil:
                if log:
                    print("z", z)
                ni = z.intersec(xi)
                if ni[0] <= ni[1]:
                    if log:
                        print("Go ", ni)
                    res = check_roots(ni, f, df, depth - 1, log)
                    if res == 'R':
                        rv = 'R'
                        break
                    elif res == 'U':
                        rv = 'U'
    if log:
        print("rv = ", rv)
    return rv


def check_box_posypkin(f, u, v, v_ival, bb):
    area_boxes = []
    border_boxes = []
    # grid = np.array(grid)
    # grid_size = len(grid) - 1
    # all_boxes = make_boxes_list(grid, dim, uniform_u)
    df_sym = []
    for i in range(len(f)):
        df_sym.append(sym.diff(f[i], v[i]))
    f_num = []
    df_num = []
    for i in range(len(f)):
        f_num.append(lambdify_f(f[i], u, v[i]))
        df_num.append(lambdify_f(df_sym[i], u, v[i]))
        # print(bb)
    rvl = []
    for i in range(len(f)):
        rvl.append(check_roots(v_ival[i], lambda v: f_num[i](bb, v), lambda v: df_num[i](bb, v), 12))
    rvl_R = [x == 'R' for x in rvl]
    rvl_N = [x == 'N' for x in rvl]
    if all(rvl_R):
        return "inside"
    elif any(rvl_N):
        return "outside"
    else:
        return "border"



def check_box(grid, dim, v_ival, extension, eps, log=False, max_iter=10, decomposition=False, eps_decomp=None, strategy = "Default",
              grid_v = None, dim_v=None, uniform_v = True, uniform_u = True, enlargement=False, boxes=[]):
    """
    Function for checking boxes on dim-dimensional uniform grid with checker method
    :param grid: 1-d grid
    :param dim: number of the dimensions
    :param v_ival: vector of not fixed interval variables
    :param extension: Extension calculator-class object, contains param info and calculate interval extension
    :param eps: error
    :param log: turn on log info printing
    :return: list of inside boxes, list of border boxes
    """
    area_boxes = []
    border_boxes = []
    if len(boxes) == 0:
        grid = np.array(grid)
        grid_size = len(grid) - 1
        all_boxes = make_boxes_list(grid, dim, uniform_u)
    else:
        all_boxes = boxes
    # print(all_boxes[:5])
    # print(diam(all_boxes[0]))
    # all_boxes = [[ival.Interval([0.0, 0.625]), ival.Interval([12.5, 13.125])]]
    # TODO:check the intersection in enlargement bisection function for box [ival.Interval([0.0, 0.625]), ival.Interval([12.5, 13.125])]]
    # all_boxes = [[ival.Interval([2.5, 3.125]), ival.Interval([0.625, 1.25])]]
    # all_boxes = [[ival.Interval([-5.625, -5.0]), ival.Interval([-10.0, -9.375])]]
    # all_boxes = [[ival.Interval([0.0, 0.625]), ival.Interval([11.25, 11.875])]]
    # all_boxes = [[ival.Interval([0.0, 0.625]), ival.Interval([11.875, 12.5])]]
    # 2 dof 32 nodes, bad box ([-2.812, -1.875], [2.812, 3.75])
    # all_boxes = [[ival.Interval([-2.812, -1.875]), ival.Interval([2.812, 3.75])]]
    # bad_boxes = [[ival.Interval([0.0, 0.625]), ival.Interval([12.5, 13.125])],
    #              [ival.Interval([-5.625, -5.0]), ival.Interval([-10.0, -9.375])],
    #              [ival.Interval([2.5, 3.125]), ival.Interval([0.625, 1.25])],
    #              [ival.Interval([-6.875, -6.25]), ival.Interval([0.625, 1.25])],
    #              [ival.Interval([8.125, 8.75]), ival.Interval([-4.375, -3.75])],
    #              [ival.Interval([5.0, 5.625]), ival.Interval([-10.0, -9.375])],
    #              [ival.Interval([-6.875, -6.25]), ival.Interval([-8.125, -7.5])],
    #              [ival.Interval([-4.375, -3.75]), ival.Interval([10.625, 11.25])]
    #              ]
    # good_inside_boxes = [[ival.Interval([-2.5, -1.875]), ival.Interval([10.0, 10.625])],
    #                      [ival.Interval([0.0, 0.625]), ival.Interval([0.0, 0.625])],
    #                      [ival.Interval([7.5, 8.125]), ival.Interval([1.25, 1.875])],
    #                      [ival.Interval([5.625, 6.25]), ival.Interval([-6.875, -6.25])],
    #                      [ival.Interval([-2.5, -1.875]), ival.Interval([-3.125, -2.5])]
    #                      ]
    # all_boxes = [good_inside_boxes[0]]
    for i, box in enumerate(all_boxes):
        # if box[0][0] == 0 and box[0][1] == 2 and box[1][0] == 2 and box[1][1] == 4:
        print(box)
        # print(i, "/", len(all_boxes) - 1)
        # print(i)
        if extension.is_elementwise:
            temp = reccur_func_elementwise(all_boxes[i], v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
        else:
            print("1D")
            temp = "not inside"
            temp = check_1d(box, v_ival, extension, False)
            print(temp)
            if temp == "not inside":
                print("Newton")
                temp = check_box_posypkin(extension.f, extension.u, extension.v, v_ival, box)
                print(temp)
                if temp == "border":
                    if strategy == "Default":
                        if enlargement:
                            v = v_ival
                            print("Krawczyk")
                            temp = reccur_func_enlarge(all_boxes[i], v, v_ival, eps, extension, max_iter, log=log,
                                                       decomposition=decomposition, eps_decomp=eps_decomp)
                            print(temp)
                        else:
                            temp = reccur_func(all_boxes[i], v_ival, eps, extension, max_iter, log=log,
                                               decomposition=decomposition, eps_decomp=eps_decomp)
                            print(temp)
                    else:
                        temp = "outside"
                        grid_v = np.array(grid_v)
                        v_boxes = make_boxes_list(grid_v, dim_v, uniform_v)
                        temp_list = []
                        for v in v_boxes:
                            if len(v) == 1:
                                v = [v[0]]
                            else:
                                v = np.array(v)
                            if enlargement:
                                temp_loc = reccur_func_enlarge(all_boxes[i], v, v_ival, eps, extension, max_iter, log=log,
                                                               decomposition=decomposition)
                            else:
                                temp_loc = reccur_func(all_boxes[i], v, eps, extension, max_iter, log=log,
                                                       decomposition=decomposition, eps_decomp=eps_decomp)
                            if temp_loc == "inside":
                                temp = "inside"
                                break
                            else:
                                temp_list.append(temp_loc)
                        if temp != "inside":
                            check = [True if temp_list[i] == "border" else False for i in range(len(temp_list))]
                            if np.any(check):
                                temp = "border"
        if temp == 'inside':
            area_boxes.append(all_boxes[i])
        elif temp == 'border':
            border_boxes.append(all_boxes[i])
    return area_boxes, border_boxes


def check_one_box(box, v_ival, extension, eps, log=False, max_iter=9, decomposition=False, strategy="Default", grid_v=None, dim_v=None, uniform_v=False):
    if strategy == "Default":
        if extension.is_elementwise:
            temp = reccur_func_elementwise(box, v_ival, eps, extension, max_iter, log=log,
                                           decomposition=decomposition)
        else:
            temp = reccur_func(box, v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
        print(temp)
    else:
        # print(grid_v)
        grid_size = len(grid_v) - 1
        temp_list = []
        temp = "outside"
        grid_v = np.array(grid_v)
        v_boxes = make_boxes_list(grid_v, dim_v, uniform_v)
        # print(v_boxes)
        for v in v_boxes:
            # print("v=", v)
            if len(v) == 1:
                v = [v[0]]
            else:
                v = np.array(v)
            print(v)
            temp_infl = reccur_func_enlarge(box, v, v_ival, eps, extension, max_iter, log=log,
                                            decomposition=decomposition)
            print(temp_infl)
            if temp_infl == "inside":
                temp = "inside"
                break
            else:
                temp_list.append(temp_infl)
        if temp != "inside":
            check = [True if temp_list[i] == "border" else False for i in range(len(temp_list))]
            if np.any(check):
                temp = "border"
        print(temp)
    return temp


def check_box_branch(ini_box, v_ival, extension, eps, eps_bnb, eps_decomp, log=False, max_iter=10, decomposition=False, strategy="Default",
              grid_v = None, dim_v=None, uniform_v = True):
    """
    Function for checking boxes on dim-dimensional uniform grid with checker method
    :param grid: 1-d grid
    :param dim: number of the dimensions
    :param v_ival: vector of not fixed interval variables
    :param extension: Extension calculator-class object, contains param info and calculate interval extension
    :param eps: error
    :param log: turn on log info printing
    :return: list of inside boxes, list of border boxes
    """
    area_boxes = []
    border_boxes = []
    queueu = []
    queueu.append(ini_box)
    s = 0
    # eps_decomp = 2*np.pi
    while len(queueu)>0:
        # print(s)
        s += 1
        # eps_decomp = eps_decomp/2
        box = queueu.pop(0)
        if extension.is_elementwise:
            temp = reccur_func_elementwise(box, v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
        else:
            if strategy == "Default":
                temp = reccur_func(box, v_ival, eps, extension, max_iter, log=log, decomposition=decomposition, eps_decomp=eps_decomp)
            else:
                temp = "outside"
                grid_v = np.array(grid_v)
                v_boxes = make_boxes_list(grid_v, dim_v, uniform_v)
                temp_list = []
                for v in v_boxes:
                    if len(v) == 1:
                        v = [v[0]]
                    else:
                        v = np.array(v)
                    temp_infl = reccur_func_enlarge(box, v, v_ival, eps, extension, max_iter, log=log,
                                               decomposition=decomposition)
                    if temp_infl == "inside":
                        temp = "inside"
                        break
                    else:
                        temp_list.append(temp_infl)
                if temp != "inside":
                    check = [True if temp_list[i] == "border" else False for i in range(len(temp_list))]
                    # print(check)
                    if np.any(check):
                        temp = "border"
        # print(temp)
        if temp == 'inside':
            area_boxes.append(box)
        elif temp == 'border':
            if diam(box) <= eps_bnb:
                print(box, diam(box), eps_bnb)
                border_boxes.append(box)
            else:
                box_l, box_r = separate_box(box)
                queueu.append(box_l)
                queueu.append(box_r)
    print("Number of boxes = ", s)
    return area_boxes, border_boxes


def check_box_branch_parallel(ini_box, v_ival, extension, eps, eps_bnb, eps_decomp, log=False, max_iter=10, decomposition=False, strategy = "Default",
              grid_v = None, dim_v=None, uniform_v = True):
    """
    Function for checking boxes on dim-dimensional uniform grid with checker method
    :param grid: 1-d grid
    :param dim: number of the dimensions
    :param v_ival: vector of not fixed interval variables
    :param extension: Extension calculator-class object, contains param info and calculate interval extension
    :param eps: error
    :param log: turn on log info printing
    :return: list of inside boxes, list of border boxes
    """
    from MPI_queue import mpi_queue
    area_boxes = []
    border_boxes = []
    queue = []
    queue.append(ini_box)
    if strategy == "Default":
        args = {"v_init": v_ival, "eps": eps, "max_iter": max_iter, "log": log,
                "decomposition": decomposition, "eps_decomp": eps_decomp}
        v_boxes = []
    else:
        args = {"v_ival": v_ival, "eps": eps, "max_iter": max_iter, "log": log}
        grid_v = np.array(grid_v)
        v_boxes = make_boxes_list(grid_v, dim_v, uniform_v)
    s = 0
    mq = mpi_queue()
    if strategy == "Default":
        if mq.flag_main:
            mq.set_function("")
            mq.set_args(queue=queue, args=args, eps_bnb=eps_bnb, type=strategy, v_boxes=v_boxes)
            mq.execute()
            area_boxes = mq.results[0]
            border_boxes = mq.results[1]
            if log:
                print ('ERRORS:', mq.errors)
        else:
            mq.execute()
    return area_boxes, border_boxes
    

def separate_box(box):
    max_width = -np.inf
    new_box_left = np.empty_like(box)
    new_box_right = np.empty_like(box)
    separ_i = 0
    for i in range(len(box)):
        if box[i].width() > max_width:
            max_width = box[i].width()
            separ_i = i
    for i in range(len(box)):
        if i == separ_i:
            left_box = ival.Interval([box[i][0], box[i].mid()])
            right_box = ival.Interval([box[i].mid(), box[i][1]])
            new_box_left[i] = left_box
            new_box_right[i] = right_box
        else:
            new_box_left[i] = box[i]
            new_box_right[i] = box[i]
    return new_box_left, new_box_right