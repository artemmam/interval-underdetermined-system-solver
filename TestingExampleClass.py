from check_box import check_box_parallel, check_box, check_box_branch
import matplotlib.pyplot as plt
from plotter_support_functions import uni_plotter
from timeit import default_timer as timer

class Example:
    def __init__(self, extension, parallel=False, record_time=False, strategy="Default", name="", path ="1", log=False):
        """
        :param f: system of equations
        :param u: box to check
        :param v: variables for checking
        :param coef: coefficient for varying lambda matrix
        """
        self.extension = extension
        self.parallel = parallel
        self.record_time = record_time
        self.strategy = strategy
        self.time = 0
        self.name = name
        self.path = path
        self.log = log

    def write_time(self, file, rank, size, n, time):
        f = open(file + ".txt", "a+")
        f.write(str(rank) + str(size) + ", " + str(n) + ", " + str(time) + "\n")
        f.close()

    def save_boxes(self, title, boxes):
        with open("./boxes/" + title, "w+") as outfile:
            outfile.write("\n".join(boxes))

    def check_box(self, grid_u, dim_u, v_ival, eps, grid_v=None, v_dim=None, uniform_u=None, uniform_v=None):
        start = timer()
        if self.parallel:
            area_boxes, border_boxes = check_box_parallel(grid_u, dim_u, v_ival, \
                               self.extension, eps,
                               strategy=self.strategy, dim_v=v_dim,
                               grid_v=grid_v, uniform_u=uniform_u, uniform_v=uniform_v, path = self.path)
        else:
            area_boxes, border_boxes = check_box(grid_u, dim_u, v_ival, \
                               self.extension, eps,
                               strategy=self.strategy, dim_v=v_dim,
                               grid_v=grid_v, uniform_u=uniform_u, uniform_v=uniform_v, log=self.log)
        end = timer()
        self.time = end - start
        return area_boxes, border_boxes

    def check_box_branch(self, box, v_ival, eps1, eps2, eps3=None, mod="Default", decomposition=False, grid_v=None, v_dim=None, uniform_v=None):
        start = timer()
        if self.parallel:
            area_boxes, border_boxes = check_box_branch(box, v_ival, \
                               self.extension, eps1, eps2,
                               strategy=self.strategy, dim_v=v_dim,
                               grid_v=grid_v, uniform_v=uniform_v)
        else:
            area_boxes, border_boxes = check_box_branch(box, v_ival, \
                               self.extension, eps1, eps2, eps2=eps3,
                               strategy=self.strategy, decomposition=decomposition, mod=mod, dim_v=v_dim,
                               grid_v=grid_v, uniform_v=uniform_v)
        end = timer()
        self.time = end - start
        return area_boxes, border_boxes


    def plotting(self, area_boxes, border_boxes, u_lims, plot_area=None, area_params=None, save_fig=False, title = "", save_fig_params = "", logger=None):
        size = len(u_lims)
        if self.parallel:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            world_size = comm.Get_size()
            rank = comm.Get_rank()
        else:
            rank = 0
        if rank == 0:
            if size!=3:
                fig1 = plt.figure(figsize=(8, 8))
                ax1 = fig1.add_subplot(1, 1, 1)
            else:
                fig1 = plt.figure(figsize=(8, 8))
                ax1 = fig1.add_subplot(111, projection='3d')
            uni_plotter(area_boxes, border_boxes, u_lims, title, size=size,
            ax=ax1, fig=fig1, logger=logger)
            if plot_area:
                plot_area(ax1, *area_params)
            if save_fig:
                print("Save figure")
                str_save_fig_params = ""
                for param in save_fig_params:
                    str_save_fig_params+=str(param) + "_"
                plt.savefig('./fig/' + title + "_" + str_save_fig_params + '.png')