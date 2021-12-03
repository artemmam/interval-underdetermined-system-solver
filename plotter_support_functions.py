import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib import colors
import numpy as np
import sys

def plot_circles(r1, r2, d, ax):
    circle = Circle((-0.5*d, 0), radius=r1, fc='y', fill=False)
    ax.add_patch(circle)
    circle = Circle((0.5*d, 0), radius=r1, fc='y', fill=False)
    ax.add_patch(circle)
    circle = Circle((-0.5*d, 0), radius=r2, fc='y', fill=False)
    ax.add_patch(circle)
    circle = Circle((0.5*d, 0), radius=r2, fc='y', fill=False)
    ax.add_patch(circle)


def plot_linear_cube(ax, x, y, z, dx, dy, dz, color='green'):
    xx = [x, x, x+dx, x+dx, x]
    yy = [y, y+dy, y+dy, y, y]
    kwargs = {'alpha': 0.5, 'color': color, 'linewidth': 1}
    ax.plot3D(xx, yy, [z]*5, **kwargs)
    ax.plot3D(xx, yy, [z+dz]*5, **kwargs)
    ax.plot3D([x, x], [y, y], [z, z+dz], **kwargs)
    ax.plot3D([x, x], [y+dy, y+dy], [z, z+dz], **kwargs)
    ax.plot3D([x+dx, x+dx], [y+dy, y+dy], [z, z+dz], **kwargs)
    ax.plot3D([x+dx, x+dx], [y, y], [z, z+dz], **kwargs)


def plot_opaque_cube(ax, x=10, y=20, z=30, dx=40, dy=50, dz=60):
    xx = np.linspace(x, x+dx, 2)
    yy = np.linspace(y, y+dy, 2)
    zz = np.linspace(z, z+dz, 2)

    xx, yy = np.meshgrid(xx, yy)

    ax.plot_surface(xx, yy, z)
    ax.plot_surface(xx, yy, z+dz)

    yy, zz = np.meshgrid(yy, zz)
    ax.plot_surface(x, yy, zz)
    ax.plot_surface(x+dx, yy, zz)

    xx, zz = np.meshgrid(xx, zz)
    ax.plot_surface(xx, y, zz)
    ax.plot_surface(xx, y+dy, zz)

class Clicker:
    def __init__(self, logger, fig):
        self.__logger = logger
        self.fig = fig
        self.cid = fig.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        ix = event.xdata
        iy = event.ydata
        self.__logger.find_box(ix, iy)




def uni_plotter(area_points, border_points, L2, title, logger = 0, ax = 0, fig = 0, size = 0):
    """
    Plotting the set of inside and border boxes
    :param area_points: area boxes set
    :param border_points: border boxes set
    :param L2: the size of boundary box
    :param title: the name of the method
    """
    #global gl_logger
    #gl_logger = logger

    plt.rcParams.update({'font.size': 18})

    if size == 1:
        if ax == 0:
            fig, ax = plt.subplots(figsize=(8, 2))
        x_min, y_min, x_max, y_max = x_lim[0] - L2 / 10, y_lim[0] - L2 / 10, x_lim[1] + L2 / 10, y_lim[1] + L2 / 10
        ax.set_xlim([x_min, x_max])
        ax.set_ylim([-1, 1])
        ax.plot([x_lim[0], x_lim[1]], [0, 0], linewidth = 12, alpha = 0.4)
        for i in range(len(area_points)):
            ax.plot([area_points[i][0][0], area_points[i][0][1]],[0, 0], marker = "|", color = "green")
        for i in range(len(border_points)):
            ax.plot([border_points[i][0][0], border_points[i][0][1]],[0, 0], marker = "|", color = "yellow")
    elif size == 2:
        if ax == 0:
            fig, ax = plt.subplots(figsize=(8, 8))
        x_lim = L2[0]  # Left border of rectangle which we use to build uniform grid
        y_lim = L2[1]
        x_min, y_min, x_max, y_max = x_lim[0] - abs(x_lim[0]) / 10, y_lim[0] - abs(y_lim[0]) / 10, \
                                     x_lim[1] + abs(x_lim[1]) / 10, y_lim[1] + abs(y_lim[1]) / 10
        ax.set_ylim([y_min, y_max])
        ax.set_xlim([x_min, x_max])
        rect1 = Rectangle([x_lim[0], y_lim[0]], x_lim[1] - x_lim[0], y_lim[1] - y_lim[0], fill=False, color='g',
                          linewidth=0.5)
        ax.add_patch(rect1)
        ax.axes.set_aspect('equal')
        for i in range(len(area_points)):  # Plot rectangles, which compose workspace area
            rect1 = Rectangle([area_points[i][0][0], area_points[i][1][0]],
                              area_points[i][0][1] - area_points[i][0][0],
                              area_points[i][1][1] - area_points[i][1][0],
                              fill=True, fc='green', color='black', linewidth=0.5, alpha=0.8)
            ax.add_patch(rect1)
        for i in range(len(border_points)):  # Plot rectangles, which compose the border of workspace area
            rect2 = Rectangle([border_points[i][0][0], border_points[i][1][0]],
                              border_points[i][0][1] - border_points[i][0][0],
                              border_points[i][1][1] - border_points[i][1][0],
                              fill=True, fc='yellow', color='black', linewidth=0.5, alpha=0.8)
            ax.add_patch(rect2)
    else:
        if ax == 0:
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(111, projection='3d')
        x_lim = L2[0]  # Left border of rectangle which we use to build uniform grid
        y_lim = L2[1]
        z_lim = L2[2]
        x_min, y_min, z_min, x_max, y_max, z_max = x_lim[0] - abs(x_lim[0]) / 10, y_lim[0] - abs(y_lim[0]) / 10, z_lim[0] - abs(z_lim[0]) / 10,\
                                     x_lim[1] + abs(x_lim[1]) / 10, y_lim[1] + abs(y_lim[1]) / 10,  z_lim[1] + abs(z_lim[1]) / 10
        ax.set_ylim([y_min, y_max])
        ax.set_xlim([x_min, x_max])
        ax.set_zlim([z_min, z_max])
        # rect1 = Rectangle([left_border, left_border], 2*right_border, 2*right_border, fill=False, color='g',
        #                   linewidth=2.0)
        plot_linear_cube(ax, x_min, y_min, z_min, abs(x_max - x_min), abs(y_max - y_min), abs(z_max - z_min), color="black")
        #ax.axes.set_aspect('equal')
        # plt.show()
        # sys.exit(1)
        for i in range(len(area_points)):  # Plot rectangles, which compose workspace area
            # rect1 = Rectangle([area_points[i][0][0], area_points[i][1][0]],
            #                   area_points[i][0][1] - area_points[i][0][0],
            #                   area_points[i][1][1] - area_points[i][1][0],
            #                   fill=True, fc='green', color='black', linewidth=1.0, alpha=1)
            # ax.add_patch(rect1)
            plot_linear_cube(ax, area_points[i][0][0], area_points[i][1][0], area_points[i][2][0],
                             area_points[i][0][1] - area_points[i][0][0], area_points[i][1][1] - area_points[i][1][0],
                             area_points[i][2][1] - area_points[i][2][0], color="green")
        for i in range(len(border_points)):  # Plot rectangles, which compose the border of workspace area
            # rect2 = Rectangle([border_points[i][0][0], border_points[i][1][0]],
            #                   border_points[i][0][1] - border_points[i][0][0],
            #                   border_points[i][1][1] - border_points[i][1][0],
            #                   fill=True, fc='yellow', color='black', linewidth=1.0, alpha=1)
            # ax.add_patch(rect2)
            plot_linear_cube(ax, border_points[i][0][0], border_points[i][1][0], border_points[i][2][0],
                             border_points[i][0][1] - border_points[i][0][0], border_points[i][1][1] - border_points[i][1][0],
                             border_points[i][2][1] - border_points[i][2][0], color="yellow")
    ax.set_title(title, fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    if logger != 0:
        Click_Obj = Clicker(logger, fig)

def plot_one_box(box, type, L2, ax = 0):
    plt.rcParams.update({'font.size': 18})
    x_lim = L2[0]  # Left border of rectangle which we use to build uniform grid
    y_lim = L2[1]
    if ax == 0:
        fig, ax = plt.subplots(figsize=(8, 8))
    x_min, y_min, x_max, y_max = x_lim[0] - abs(x_lim[0]) / 10, y_lim[0] - abs(y_lim[0]) / 10, \
                                 x_lim[1] + abs(x_lim[1]) / 10, y_lim[1] + abs(y_lim[1]) / 10
    ax.set_ylim([y_min, y_max])
    ax.set_xlim([x_min, x_max])
    rect1 = Rectangle([x_lim[0], y_lim[0]], x_lim[1] - x_lim[0], y_lim[1] - y_lim[0], fill=False, color='g',
                      linewidth=0.5)
    ax.add_patch(rect1)
    ax.axes.set_aspect('equal')
    if type == "inside":
        color = "green"
    elif type == "border":
        color = "yellow"
    else:
        color = "black"
    rect1 = Rectangle([box[0][0], box[1][0]],
                      box[0][1] - box[0][0],
                      box[1][1] - box[1][0],
                      fill=True, fc=color, color='black', linewidth=0.5, alpha=0.8)
    ax.add_patch(rect1)


def iter_plot(s_all, N):
    s_class = s_all[:N**2].reshape(N, N)
    s_bic = s_all[N**2:].reshape(N, N)
    fig = plt.figure(figsize=(16, 8))
    grid = AxesGrid(fig, 111,
                    nrows_ncols=(1, 2),
                    axes_pad=0.5,
                    cbar_mode='single',
                    cbar_location='right',
                    cbar_pad=0.1
                    )
    ax1 = grid[0].imshow(s_class, cmap = plt.cm.get_cmap('viridis_r', 12), interpolation='none')
    ax2 = grid[1].imshow(s_bic, cmap = plt.cm.get_cmap('viridis_r', 12), interpolation='none')
    grid[0].set_title("Classical")
    grid[1].set_title("Bicentered")
    ax = [ax1, ax2]
    vmin = min(image.get_array().min() for image in ax)
    vmax = max(image.get_array().max() for image in ax)
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    for im in ax:
        im.set_norm(norm)
    cbar1 = grid.cbar_axes[0].colorbar(ax1)
    plt.show()


def plot_all_methods(methods, points, L2, example, size, L1v = 0, L2v = 0, d = 0, r = 0):
    n = len(methods)
    if size == 1:
        fig, ax = plt.subplots(n, 1, figsize=(4, 12))
    else:
        fig, ax = plt.subplots(1, n, figsize=(12, 4))
    for i, method in enumerate(methods):
        uni_plotter(points[method].area_points, points[method].border_points, L2, method, logger=0, ax=ax[i], size = size)
        if example == "2rpr":
            plot_circles(L1v, L2v, d, ax[i])
        elif example == "circle":
            circle = plt.Circle((0, 0), radius=1, fc='y', fill=False)
            ax[i].add_patch(circle)
        elif example =="2links1angle":
            ax[i].plot([-r, r], [0, 0], color = "r", lw = 6, alpha = 0.3)







"""
def iter_plot(s_class, N, title):
    s_class = s_class.reshape(N, N)
    fig, ax = plt.subplots(figsize=(8, 8))
    plt.imshow(s_class, cmap = plt.cm.get_cmap('viridis_r', 10), interpolation='none')
    #plt.pcolor(s_class, cmap=plt.get_cmap('viridis_r', 50))

    plt.colorbar()
    ax.set_title(title + " method")
    """