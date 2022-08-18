import numpy as np
from check_box import make_boxes_list
from check_box import reccur_func, reccur_func_elementwise, check_one_box
import os
import pyautogui


class Logger:

    def __init__(self, grid, size, v_init, eps, ext_calcul, uniform_u=True, decomp=False, strategy="Default", grid_v=None, dim=None):
        """
        :param grid: grid
        :param size: the size of the grid
        :param v_init: initial value for V
        :param eps: error
        :param ext_calcul: extension calculator class
        """
        self.__grid = grid
        self.__size = size
        self.__v_init = v_init
        self.__eps = eps
        self.__ext_calcul = ext_calcul
        self.__decomp = decomp
        self.__uniform_u = uniform_u
        self.__strategy = strategy
        self.__grid_v = grid_v
        self.__dim = dim
    def find_box(self, x, y):
        L = len(self.__v_init)
        print("\n\n\n")
        print("LOGGING")
        all_boxes = make_boxes_list(self.__grid, self.__size, self.__uniform_u)
        print(all_boxes)
        if self.__size == 2:
            for box in all_boxes:
                if (x > box[0][0] and x<box[0][1] and y > box[1][0] and y<box[1][1]):
                    print("Box:", box)
                    break
        if self.__size == 1:
            for box in all_boxes:
                if x > box[0][0] and x<box[0][1]:
                    print("Box:", box)
                    break
        print("v_init:", self.__v_init)
        # print("Lambda", self.__ext_calcul.calculate_lam(self.__v_init, box).reshape((L, L)))
        if self.__strategy == "Inflation":
            check = check_one_box(box, self.__v_init, self.__ext_calcul, self.__eps, log=True, strategy=self.__strategy, grid = self.__grid_v, dim=self.__dim)
        else:
            if self.__ext_calcul.is_elementwise:
                check = reccur_func_elementwise(box, self.__v_init, self.__eps, self.__ext_calcul, log=True, decomposition=self.__decomp)
            else:
                check = reccur_func(box, self.__v_init, self.__eps, self.__ext_calcul, log=True,
                                    decomposition=self.__decomp)
        print(check)
        print("-----\n\n\n")