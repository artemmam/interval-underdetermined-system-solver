import interval as ival
import numpy as np
import itertools as it

from timeit import default_timer as timer
# from pathos.multiprocessing import ProcessingPool as Pool
# import dill  # the code below will fail without this line
# import pickle


def write_time_per_proc(file, rank, time, message = ""):
    f = open(file + ".txt", "a+")
    f.write(message + " | " +  str(rank) + ": " + str(time) + "\n")
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


def reccur_func(box, v_init, eps, extension, max_iter=10, log=False, decomposition=False, level=0):
    v_iter = v_init.copy()
    n = len(v_init)
    v_prev = v_iter.copy()
    v_ext = v_iter.copy()
    new_v_left = np.empty_like(v_init)
    new_v_right = np.empty_like(v_init)
    new_v_middle = np.empty_like(v_init)
    # print("box", box)
    k = 0
    f_num = extension.lambdify_f()
    # if log:
    #     f_centered = extension.centered_form()
    while True:
        if log:
            print("*****")
            print("Natural :", f_num(v_iter, box).reshape(-1))
        # for nat_ext in f_num(v_iter, box).reshape(-1):
        #     if not ival.Interval([0, 0]).isInIns(nat_ext):
        #         return "outside"
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
        # w_x = ival.n_width(v_iter)
        # kr_minux_m_x = []
        for i in range(n):
            # kr_minux_m_x.append(v_ext[i] - v_iter[i].mid())
            if not (v_ext[i].isIn(v_iter[i])):
                check = False
                break
        # norm_buf = ival.norm(kr_minux_m_x)
        #
        # if not norm_buf < w_x / 2:
        #     check = False
        if check:
            return "inside"
        for i in range(n):
            if v_iter[i].isNoIntersec(v_ext[i]):
                return "outside"
            else:
                v_iter[i] = v_iter[i].intersec(v_ext[i])
        if k > max_iter:
            if decomposition:
                if level < 3:
                    if log:
                        print("Decomposition")
                        print("Level: ", level)
                    for j in range(n):
                        separ_i = j
                        for i in range(n):
                            if i == separ_i:
                                d = v_init[i][1] - v_init[i][0]
                                v_left = ival.Interval([v_init[i][0], v_init[i][0] + d / 3])
                                v_middle = ival.Interval([v_init[i][0] + d / 3, v_init[i][1] - d / 3])
                                v_right = ival.Interval([v_init[i][1] - d / 3, v_init[i][1]])
                                new_v_left[i] = v_left
                                new_v_middle[i] = v_middle
                                new_v_right[i] = v_right
                            else:
                                new_v_left[i] = v_init[i]
                                new_v_middle[i] = v_init[i]
                                new_v_right[i] = v_init[i]
                        if log:
                            print("Left: ", new_v_left)
                            print("Middle: ", new_v_middle)
                            print("Right:", new_v_right)
                            print("GOING LEFT")
                        left_check = reccur_func(box, new_v_left, eps, extension, max_iter, log, decomposition, level=level+1)
                        if left_check == "inside":
                            return "inside"
                        if log:
                            print("GOING MIDDLE")
                        middle_check = reccur_func(box, new_v_middle, eps, extension, max_iter, log, decomposition, level=level+1)
                        if middle_check == "inside":
                            if log:
                                print(middle_check)
                            return "inside"
                        if log:
                            print("GOING RIGHT")
                        right_check = reccur_func(box, new_v_right, eps, extension, max_iter, log, decomposition, level=level+1)
                        if right_check == "inside":
                            return "inside"
            return "border"
        v_prev = v_iter.copy()
        k += 1


def reccur_func_enlarge(box, v_init, v_ival, eps, extension, max_iter=10, log=False, decomposition=False, level=0):
    max_iter = 10
    v_iter = v_init.copy()
    n = len(v_init)
    v_prev = v_iter.copy()
    v_ext = v_iter.copy()
    new_v_left = np.empty_like(v_init)
    new_v_right = np.empty_like(v_init)
    new_v_middle = np.empty_like(v_init)
    # print("box", box)
    k = 0
    f_num = extension.lambdify_f()
    # if log:
    #     f_centered = extension.centered_form()
    while True:
        if log:
            print("*****")
            print("Natural :", f_num(v_iter, box).reshape(-1))
        for nat_ext in f_num(v_iter, box).reshape(-1):
            if not ival.Interval([0, 0]).isInIns(nat_ext):
                return "outside"
        if k>5:
            for i in range(n):
                v_iter[i].scale(1.1)
                v_iter[i] = v_iter[i].intersec(v_ival[i])
        check = True
        v_ext = extension.calculate_extension(box, v_iter, log=log).reshape(-1)
        if log:
            print("Number of iteration =", k)
            print("box", box)
            print("V_ival", v_ival)
            print("Old V = ", v_iter)
            print("New V = ", v_ext)
        # w_x = ival.n_width(v_iter)
        # kr_minux_m_x = []
        # for i in range(n):
        #     kr_minux_m_x.append(v_ext[i] - v_iter[i].mid())
        for i in range(n):
            if not (v_ext[i].isIn(v_iter[i])):
                check = False
                break
        # norm_buf = ival.norm(kr_minux_m_x)
        #
        # if log:
        #     print("Norm = ", norm_buf, 3)
        #     print("N-width/2 = ", w_x/2, 3)
        # for i in range(n):
        #     if not norm_buf[i] < w_x[i]/2:
        #         check = False
        # if not norm_buf < w_x / 2:
        #     check = False
        if check:
            return "inside"
        for i in range(n):
            if v_ext[i].isNoIntersec(v_ival[i]):
                return "outside"
            else:
                v_iter[i] = v_ext[i].intersec(v_ival[i])
        if abs(diam(v_iter) - diam(v_ival)) / (0.5 * abs(diam(v_iter) + diam(v_ival))) < eps or k > max_iter:
            if decomposition:
                if level < 3:
                    if log:
                        print("Decomposition")
                        print("Level: ", level)
                    for j in range(n):
                        separ_i = j
                        for i in range(n):
                            if i == separ_i:
                                d = v_init[i][1] - v_init[i][0]
                                v_left = ival.Interval([v_init[i][0], v_init[i][0] + d / 3])
                                v_middle = ival.Interval([v_init[i][0] + d / 3, v_init[i][1] - d / 3])
                                v_right = ival.Interval([v_init[i][1] - d / 3, v_init[i][1]])
                                new_v_left[i] = v_left
                                new_v_middle[i] = v_middle
                                new_v_right[i] = v_right
                            else:
                                new_v_left[i] = v_init[i]
                                new_v_middle[i] = v_init[i]
                                new_v_right[i] = v_init[i]
                        if log:
                            print("Left: ", new_v_left)
                            print("Middle: ", new_v_middle)
                            print("Right:", new_v_right)
                            print("GOING LEFT")
                        left_check = reccur_func(box, new_v_left, eps, extension, max_iter, log, decomposition, level=level+1)
                        if left_check == "inside":
                            return "inside"
                        if log:
                            print("GOING MIDDLE")
                        middle_check = reccur_func(box, new_v_middle, eps, extension, max_iter, log, decomposition, level=level+1)
                        if middle_check == "inside":
                            if log:
                                print(middle_check)
                            return "inside"
                        if log:
                            print("GOING RIGHT")
                        right_check = reccur_func(box, new_v_right, eps, extension, max_iter, log, decomposition, level=level+1)
                        if right_check == "inside":
                            return "inside"
            return "border"
        v_prev = v_iter.copy()
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


def check_box_parallel(grid, dim, v_ival, extension, eps, log=False, max_iter=10, decomposition=False, strategy = "Default",
              grid_v = None, dim_v=None, uniform_v = True, uniform_u = True):
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
    grid = np.array(grid)
    grid_size = len(grid) - 1
    all_boxes = make_boxes_list(grid, dim, uniform_u)
    for i, box in enumerate(all_boxes):
        if i % size == rank:
            # print("rank = ", rank, "; ", i, "/", len(all_boxes) - 1)
            # print(i, "/", len(all_boxes) - 1)
            if extension.is_elementwise:
                temp = reccur_func_elementwise(all_boxes[i], v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
            else:
                if strategy == "Default":
                    temp = reccur_func(all_boxes[i], v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
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
                        temp_infl = reccur_func_enlarge(all_boxes[i], v, v_ival, eps, extension, max_iter, log=log,
                                                   decomposition=decomposition)
                        if temp_infl == "inside":
                            temp = "inside"
                            break
                        else:
                            temp_list.append(temp_infl)
                    if temp!="inside":
                        check = [True if temp_list[i] == "border" else False for i in range(len(temp_list))]
                        if np.any(check):
                            temp = "border"
            if temp == 'inside':
                area_boxes.append(all_boxes[i])
            elif temp == 'border':
                border_boxes.append(all_boxes[i])
    end = timer()
    write_time_per_proc("bicentered_krawczyk_enlarge_time_procs", rank, end - start, "total = " + str(len(all_boxes)/size) + "; area = " + str(len(area_boxes)) +  "; border = " + str(len(border_boxes)))
    start = timer()
    area_boxes_g = comm.gather(area_boxes, root=0)
    border_boxes_g = comm.gather(border_boxes, root=0)
    end = timer()
    write_time_per_proc("bicentered_krawczyk_enlarge_time_procs", rank, end - start, "GATHER")
    area_boxes_g_unpacking = []
    border_boxes_g_unpacking = []
    if rank == 0:
        area_boxes_g_unpacking = unpack_array(area_boxes_g)
        border_boxes_g_unpacking = unpack_array(border_boxes_g)
    return area_boxes_g_unpacking, border_boxes_g_unpacking
    

def check_box(grid, dim, v_ival, extension, eps, log=False, max_iter=10, decomposition=False, strategy = "Default",
              grid_v = None, dim_v=None, uniform_v = True, uniform_u = True):
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
    grid = np.array(grid)
    grid_size = len(grid) - 1
    all_boxes = make_boxes_list(grid, dim, uniform_u)
    # print(all_boxes[81])
    for i, box in enumerate(all_boxes):
        print(i, "/", len(all_boxes) - 1)
        # print(i)
        if extension.is_elementwise:
            temp = reccur_func_elementwise(all_boxes[i], v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
        else:
            if strategy == "Default":
                temp = reccur_func(all_boxes[i], v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
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
                    # print(v)
                    temp_infl = reccur_func_enlarge(all_boxes[i], v, v_ival, eps, extension, max_iter, log=log,
                                               decomposition=decomposition)
                    if temp_infl == "inside":
                        temp = "inside"
                        break
                    else:
                        temp_list.append(temp_infl)
                if temp!="inside":
                    check = [True if temp_list[i] == "border" else False for i in range(len(temp_list))]
                    # print(check)
                    if np.any(check):
                        temp = "border"
        if temp == 'inside':
            area_boxes.append(all_boxes[i])
        elif temp == 'border':
            border_boxes.append(all_boxes[i])
    return area_boxes, border_boxes


def check_one_box(box, v_ival, extension, eps, log=False, max_iter=9, decomposition=False, strategy = "Default", grid = None, dim = None, uniform=True):
    if strategy == "Default":
        if extension.is_elementwise:
            temp = reccur_func_elementwise(box, v_ival, eps, extension, max_iter, log=log,
                                           decomposition=decomposition)
        else:
            temp = reccur_func(box, v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
        print(temp)
    else:
        grid_v = np.array(grid)
        grid_size = len(grid_v) - 1
        temp_list = []
        temp = "outside"
        v_boxes = make_boxes_list(grid, dim, uniform=uniform)
        for v in v_boxes:
            print("v=", v)
            if len(v) == 1:
                v = [v[0]]
            else:
                v = np.array(v)
            # print(v)
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
            # print(check)
            if np.any(check):
                temp = "border"
        print(temp)
    return temp
