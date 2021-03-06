import numpy as np
from check_box import make_boxes_list
from check_box import reccur_func
import os
import pyautogui


class Logger:

    def __init__(self, grid, size, v_init, eps, ext_calcul, decomp=False):
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

    def find_box(self, x, y):
        L = len(self.__v_init)
        print("\n\n\n")
        print("LOGGING")
        all_boxes = make_boxes_list(self.__grid, self.__size)
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
        check = reccur_func(box, self.__v_init, self.__eps, self.__ext_calcul, log=True, decomposition=self.__decomp)
        print(check)
        print("-----\n\n\n")