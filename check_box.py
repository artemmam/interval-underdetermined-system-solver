import interval as ival
import numpy as np
import itertools as it
from timeit import default_timer as timer
# from pathos.multiprocessing import ProcessingPool as Pool
# import dill  # the code below will fail without this line
# import pickle


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


def reccur_func(box, v_init, eps, extension, max_iter=10, log=False, decomposition=False, eps_decomp=np.pi):
    v_iter = v_init.copy()
    n = len(v_init)
    v_prev = v_iter.copy()
    v_ext = v_iter.copy()
    new_v_left = np.empty_like(v_init)
    new_v_right = np.empty_like(v_init)
    new_v_middle = np.empty_like(v_init)
    k = 0
    f_num = extension.lambdify_f()
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
            print(v_iter, box)
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
                # if v_iter[i].isIn(v_ext[i]):
                #     v_iter[i] = v_ext[i].intersec(v_iter[i])
                # else:
                v_iter[i] = v_ext[i].intersec(v_ival[i])
        if abs(diam(v_iter) - diam(v_ival)) / (0.5 * abs(diam(v_iter) + diam(v_ival))) < eps or k > max_iter:
            if decomposition:
                if log:
                    print("Bisection")
                    print(v_iter, diam(v_iter), eps_decomp)
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


def check_box_parallel(grid, dim, v_ival, extension, eps, log=False, max_iter=10, decomposition=False, eps_decomp=None, strategy = "Default",
              grid_v = None, dim_v=None, uniform_v = True, uniform_u = True, path = "", enlargement=False):
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
                temp = reccur_func_elementwise(all_boxes[i], v_ival, eps, extension, max_iter, log=log,
                                               decomposition=decomposition)
            else:
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
    

def check_box(grid, dim, v_ival, extension, eps, log=False, max_iter=10, decomposition=False, eps_decomp=None, strategy = "Default",
              grid_v = None, dim_v=None, uniform_v = True, uniform_u = True, enlargement=False):
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
    # print(diam(all_boxes[0]))
    all_boxes = [[ival.Interval([0.0, 0.625]), ival.Interval([0.0, 0.625])]]
    for i, box in enumerate(all_boxes):
        # print(box)
        # print(i, "/", len(all_boxes) - 1)
        # print(i)
        if extension.is_elementwise:
            temp = reccur_func_elementwise(all_boxes[i], v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
        else:
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
            # print(check)
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
        # elif temp == 'outside' and diam(box)>eps1:
        #     box_l, box_r = separate_box(box)
        #     queueu.append(box_l)
        #     queueu.append(box_r)
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
    # from mpi4py import MPI
    # from mpi4py.futures import MPIPoolExecutor, MPICommExecutor
    # comm = MPI.COMM_WORLD
    from MPI_queue import mpi_queue
    area_boxes = []
    border_boxes = []
    # with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
    # with MPIPoolExecutor() as executor:
    #     if executor is not None:
    # print("Executor rank: ", comm.Get_rank())
    # print("World size: ", comm.Get_size())
    queue = []
    queue.append(ini_box)
    # args = [v_ival, eps]
    if strategy == "Default":
        args = {"v_init": v_ival, "eps": eps, "max_iter": max_iter, "log": log,
                "decomposition": decomposition, "eps_decomp": eps_decomp}
        v_boxes = []
    else:
        args = {"v_ival": v_ival, "eps": eps, "max_iter": max_iter, "log": log}
        grid_v = np.array(grid_v)
        v_boxes = make_boxes_list(grid_v, dim_v, uniform_v)
    # print(all_boxes[81])
    s = 0
    mq = mpi_queue()
    if strategy == "Default":
        if mq.flag_main:
            mq.set_function("")
            mq.set_args(queue=queue, args=args, eps_bnb=eps_bnb, type=strategy, v_boxes=v_boxes)
            mq.execute()
            area_boxes = mq.results[0]
            border_boxes = mq.results[1]
            # time.sleep(1)
            # print ('RESULTS:', mq.results)
            if log:
                print ('ERRORS:', mq.errors)
        else:
            mq.execute()
    #
    # while len(queueu)>0:
    #     # print(queueu)
    #     s+=1
    #     box = queueu.pop(0)
    #     # print(box)
    #     # print(i, "/", len(all_boxes) - 1)
    #     # print(i)
    #     if extension.is_elementwise:
    #         temp = reccur_func_elementwise(box, v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
    #     else:
    #         if strategy == "Default":
    #             temp = reccur_func(box, v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
    #         else:
    #             temp = "outside"
    #             grid_v = np.array(grid_v)
    #             v_boxes = make_boxes_list(grid_v, dim_v, uniform_v)
    #             temp_list = []
    #             for v in v_boxes:
    #                 if len(v) == 1:
    #                     v = [v[0]]
    #                 else:
    #                     v = np.array(v)
    #                 # print(v)
    #                 temp_infl = reccur_func_enlarge(box, v, v_ival, eps, extension, max_iter, log=log,
    #                                            decomposition=decomposition)
    #                 if temp_infl == "inside":
    #                     temp = "inside"
    #                     break
    #                 else:
    #                     temp_list.append(temp_infl)
    #             if temp!="inside":
    #                 check = [True if temp_list[i] == "border" else False for i in range(len(temp_list))]
    #                 # print(check)
    #                 if np.any(check):
    #                     temp = "border"
    #     # print(temp)
    #     if temp == 'inside':
    #         area_boxes.append(box)
    #     elif temp == 'border':
    #         if diam(box)<eps1:
    #             border_boxes.append(box)
    #         else:
    #             box_l, box_r = separate_box(box)
    #             queueu.append(box_l)
    #             queueu.append(box_r)
    #     # elif temp == 'outside' and diam(box)>eps1:
    #     #     box_l, box_r = separate_box(box)
    #     #     queueu.append(box_l)
    #     #     queueu.append(box_r)
    # print("Number of boxes = ", s)
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