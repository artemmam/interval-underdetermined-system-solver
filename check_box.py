import interval as ival
import numpy as np
import itertools as it
from pathos.multiprocessing import ProcessingPool as Pool
import dill  # the code below will fail without this line
import pickle


def diam(A):
    s = 0
    for v in A:
        s += v.width()**2
    return np.sqrt(s)

def make_boxes_list(grid, dim):
    """
    Make list of boxes in dim dimension from vector grid
    :param grid: vector on which grid is constructed
    :param dim:  the dimensional of grid
    :return: the list of boxes in dim
    """
    grid_size = len(grid) - 1
    grid_intervals = []
    grid_variants = []
    for i in range(grid_size):
        grid_intervals.append(ival.Interval([grid[i], grid[i + 1]]))
    for i in range(dim):
        grid_variants.append(grid_intervals)
    grid_n_size = list(it.product(*grid_variants))
    return grid_n_size


def reccur_func(box, v_init, eps, extension, max_iter=20, log=False, decomposition=False, level=0):
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
    if log:
        f_centered = extension.centered_form()
    while True:
        if log:
            print("Natural :", f_num(v_iter, box).reshape(-1))
            c = []
            for i in range(n):
                for j in range(n):
                    c.append(v_iter[i].mid())
            c = np.array(c).reshape(n, n).T.reshape(n * n)
            print("Centered :", f_centered(v_iter, c, box))
        for nat_ext in f_num(v_iter, box).reshape(-1):
            if not ival.Interval([0, 0]).isInIns(nat_ext):
                return "outside"

        check = True
        v_ext = extension.calculate_extension(box, v_iter, log=log).reshape(-1)
        if log:
            print("*****")
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
        if abs(diam(v_iter) - diam(v_prev)) / (0.5 * abs(diam(v_iter) + diam(v_prev))) < eps or k > max_iter:
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


def reccur_func_elementwise(box, v_init, eps, max_iter, extension, log=False):
    v_iter = v_init.copy()
    n = len(v_init)
    v_prev = v_iter.copy()
    # v_ext = v_iter[0].copy()
    # print("box", box)
    k = 0
    while True:
        check = True
        v_ext_funcs, param, c = extension.calculate_extension(box, v_iter)
        c = np.array(c).reshape(n, n)
        # print("*****")
        # print("Number of iteration =", k)
        # print("Old V = ", v_iter)
        for i in range(n):
            v_ext = v_ext_funcs[i](v_iter, c[i], param).reshape(-1)[0]
            # print("New V = ", v_ext)
            if not (v_ext.isIn(v_iter[i])):
                v_iter[i] = v_iter[i].intersec(v_ext)
                check = False
            if v_iter[i].isNoIntersec(v_ext):
                return "outside"
        if check:
            return "inside"
        if abs(diam(v_iter) - diam(v_prev)) / (0.5 * abs(diam(v_iter) + diam(v_prev))) < eps or k > max_iter:
            return "border"
        v_prev = v_iter.copy()
        k += 1


def check_box(grid, dim, v_ival, extension, eps, log=False, max_iter=20, decomposition=False):
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
    all_boxes = make_boxes_list(grid, dim)
    for i in range(grid_size**dim):
        # print(i)
        if extension.is_elementwise:
            temp = reccur_func_elementwise(all_boxes[i], v_ival, eps, max_iter, extension, log=log)
        else:
            temp = reccur_func(all_boxes[i], v_ival, eps, extension, max_iter, log=log, decomposition=decomposition)
        if temp == 'inside':
            area_boxes.append(all_boxes[i])
        elif temp == 'border':
            border_boxes.append(all_boxes[i])
    return area_boxes, border_boxes
