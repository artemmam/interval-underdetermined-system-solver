from symbolic_support_functions import derived_f, classical_krawczyk_extension, derived_recurrent_form, \
    hansen_sengupta_extension, centered_form, function_replacer
import numpy as np
import interval as ival
import sympy as sym



class BaseExtension:
    def __init__(self, f, v, u, coef=1, is_elementwise=False):
        """
        :param f: system of equations
        :param u: box to check
        :param v: variables for checking
        :param coef: coefficient for varying lambda matrix
        """
        self.__f = f
        self.__u = u
        self.__v = v
        self.__coef = coef
        self.__is_elementwise = is_elementwise
        self.__numeric_extension = self.get_numeric_extension()
        self.__derived_f = derived_f(self.f, self.v, self.u)

    @property
    def f(self):
        return self.__f

    @property
    def v(self):
        return self.__v

    @property
    def u(self):
        return self.__u

    @property
    def coef(self):
        return self.__coef

    @property
    def numeric_extension(self):
        return self.__numeric_extension

    @property
    def is_elementwise(self):
        return self.__is_elementwise

    def get_numeric_extension(self):
        ...

    def calculate_lam(self, v, u, coef=1):
        """
        Function for calculation matrix lambda (L = (mid(F'))**-1)
        :param v: variables
        :param u: box to check
        :return: matrix lambda
        """
        f_derived_num = self.__derived_f
        f_derived_num = f_derived_num(v, u)
        n = len(v)
        m = np.zeros_like(f_derived_num)
        for i in range(len(v)):
            for j in range(len(v)):
                m[i, j] = coef*ival.valueToInterval(f_derived_num[i, j]).mid()
        m = m.astype(np.float64)
        if np.linalg.det(m) == 0:
            return np.linalg.inv(m + np.eye(n)).reshape(n*n)
        else:
            return np.linalg.inv(m).reshape(n*n)

    def lambdify_f(self):
        f = function_replacer(self.f)
        return sym.lambdify([self.v, self.u], f)

    def centered_form(self):
        param = self.u
        c = []
        n = len(self.v)
        for i in range(len(self.v)):
            for j in range(len(self.v)):
                c.append(sym.symbols("c" + str(i) + str(j)))
        c = np.array(c).reshape(n, n).T.reshape(n * n)
        # centered_form_num = centered_form(self.f, self.v, c, param)
        return centered_form(self.f, self.v, c, param)[0]


class ClassicalKrawczykExtension(BaseExtension):
    def get_numeric_extension(self):
        """
                Function for transformation symbolic Krawczyk operator function into numeric
                :return: numeric Krawczyk operator calculation function
                """
        c = []
        lam = []
        for i in range(len(self.v)):
            for j in range(len(self.v)):
                c.append(sym.symbols("c" + str(i) + str(j)))
                lam.append(sym.symbols("lam" + str(i) + str(j)))
        if self.is_elementwise:
            return classical_krawczyk_extension(self.f, self.u, self.v, lam, c)[1]
        else:
            return classical_krawczyk_extension(self.f, self.u, self.v, lam, c)[0]

    def calculate_extension(self, box, v, log = False):
        """
        Function for calculation interval Krawczyk extension
        :param box: box to check
        :param V: variables for checking
        :return: interval vector
        """
        n = len(v)
        lamda_matrix = self.calculate_lam(v, box, self.coef)
        param = [box] + [lamda_matrix]
        c = []
        for i in range(n):
            for j in range(n):
                c.append(v[i].mid())
        c = np.array(c).reshape(n, n).T.reshape(n*n)
        if self.is_elementwise:
            return self.numeric_extension, param, c
        else:
            return np.array(self.numeric_extension(v, c, param))


class BicenteredKrawczykExtension(BaseExtension):
    def __init__(self, f, u, v, coef=1, is_elementwise=False):
        """
        :param f: system of equations
        :param u: box to check
        :param v: variables for checking
        :param coef: int, coefficient for Bicentered Krawczyk extension
        """
        super().__init__(f, u, v, coef, is_elementwise=False)
        self.__numeric_derived_recurrent_form = self.get_numeric_derived_recurrent_form()

    @property
    def numeric_derived_recurrent_form(self):
        return self.__numeric_derived_recurrent_form

    def get_numeric_extension(self):
        """
                Function for transformation symbolic Krawczyk operator function into numeric
                :return: numeric Krawczyk operator calculation function
                """
        c = []
        lam = []
        for i in range(len(self.v)):
            for j in range(len(self.v)):
                c.append(sym.symbols("c" + str(i) + str(j)))
                lam.append(sym.symbols("lam" + str(i) + str(j)))
        if self.is_elementwise:
            return classical_krawczyk_extension(self.f, self.u, self.v, lam, c)[1]
        else:
            return classical_krawczyk_extension(self.f, self.u, self.v, lam, c)[0]

    def get_numeric_derived_recurrent_form(self):
        lam = []
        for i in range(len(self.v)):
            for j in range(len(self.v)):
                lam.append(sym.symbols("lam" + str(i) + str(j)))
        return derived_recurrent_form(self.f, self.v, self.u, lam)


    def calcul_new_c(self, V, L, box):
        """
        Function for calculation cmin and cmax for bicentered Krawczyk
        :param V: variables for checking
        :param L: lambda-matrix
        :param box: box to check
        :return: intervals c_min and c_max
        """
        param = [box] + [L]
        new_v = self.numeric_derived_recurrent_form(V, param)
        new_v = new_v.reshape(len(V), len(V))
        n = len(new_v)
        c_max = np.zeros_like(new_v)
        c_min = np.zeros_like(new_v)
        for i in range(n):
            for j in range(n):
                if new_v[i][j][1] <= 0:
                    c_min[i][j] = V[j][1]
                elif new_v[i][j][0] >= 0:
                    c_min[i][j] = V[j][0]
                else:
                    c_min[i][j] = (new_v[i][j][1] * V[j][0] - new_v[i][j][0] * V[j][1]) / (new_v[i][j][1] - new_v[i][j][0])
        for i in range(n):
            for j in range(n):
                if new_v[i][j][1] <= 0:
                    c_max[i][j] = V[j][0]
                elif new_v[i][j][0] >= 0:
                    c_max[i][j] = V[j][1]
                else:
                    c_max[i][j] = (new_v[i][j][0] * V[j][0] - new_v[i][j][1] * V[j][1]) / (new_v[i][j][0] - new_v[i][j][1])
        return c_min, c_max

    def calculate_extension(self, box, v, log = False):
        """
        Function for calculation interval Bicentered Krawczyk extension
        :param box: box to check
        :param V: variables for checking
        :return: interval vector
        """
        lambda_matrix = self.calculate_lam(v, box, self.coef)
        param = [box] + [lambda_matrix]
        c_min, c_max = self.calcul_new_c(v, lambda_matrix, box)
        c_min = c_min.reshape(len(v) * len(v))
        c_max = c_max.reshape(len(v) * len(v))
        v_ext_min, v_ext_max = np.array(self.numeric_extension(v, c_min, param)), \
                               np.array(self.numeric_extension(v, c_max, param))
        v_bic = []
        for i in range(len(v)):
            v_bic.append(v_ext_min[i][0].intersec(v_ext_max[i][0]))
        return np.array(v_bic).T


class HansenSenguptaExtension(BaseExtension):
    def get_numeric_extension(self):
        """
                Function for transformation symbolic extension function into numeric
                :return: numeric extension function
                """
        c = []
        lam = []
        for i in range(len(self.v)):
            for j in range(len(self.v)):
                c.append(sym.symbols("c" + str(i) + str(j)))
                lam.append(sym.symbols("lam" + str(i) + str(j)))
        if self.is_elementwise:
            return hansen_sengupta_extension(self.f, self.u, self.v, lam, c)[1]
        else:
            return hansen_sengupta_extension(self.f, self.u, self.v, lam, c)[0]

    def calculate_extension(self, box, v, log = False):
        """
        Function for calculation interval extension
        :param box: box to check
        :param v: variables for checking
        :return: interval vector
        """
        n = len(v)
        lamda_matrix = self.calculate_lam(v, box, self.coef)
        param = [box] + [lamda_matrix]
        c = []
        for i in range(n):
            for j in range(n):
                c.append(v[i].mid())
        c = np.array(c).reshape(n, n).T.reshape(n*n)
        if self.is_elementwise:
            return self.numeric_extension, param, c
        else:
            if log:
                print("///")
                print("HS input")
                print("v", v)
                print("c", c)
                print("lambda", lamda_matrix )
                print("///")
            return np.array(self.numeric_extension(v, c[0:2], param))