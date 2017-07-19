#!/bin/python
import sgems
import math
import numpy as np
import matplotlib.pyplot as pltp
import pylab as plt
import scipy
import copy
from utils import UnionFind, Painter, Map3D
from pylab import plt
import numpy as np

"""
The base python sgems plugin.
"""
def read_params(a,j=''):
  for i in a:
    if (type(a[i])!=type({'a':1})):
      print j+"['"+str(i)+"']="+str(a[i])
    else:
      read_params(a[i],j+"['"+str(i)+"']")

class pdfsequential12:
    def __init__(self):
        pass

    def initialize(self, params):
        self.params=params
        return True

    def execute(self):
        head_property = self.params['head_prop']['property']
        head_grid = self.params['head_prop']['grid']

        tail_property = self.params['tail_prop']['property']
        tail_grid = self.params['tail_prop']['grid']

        #less_iqual = None
        #more_iqual = None
        cut_off = float(self.params['cut_off']['value'])
        new_prop = self.params['new_prop']['value']
        new_prop1 = self.params['new_prop1']['value']
        new_prop2 = self.params['new_prop2']['value']
        #mean = float(self.params['mean']['value'])
        #median = float(self.params['median']['value'])
        #variance = float(self.params['variance']['value'])
        fases = int(self.params['fases']['value'])
        simulation = int(self.params['simulation']['value'])
        runs = int(self.params['runs']['value'])
        #x = float(self.params['x']['value'])
        #y = float(self.params['y']['value'])
        #z = float(self.params['z']['value'])
        less_equal = self.params['rsless']['value']
        more_equal = self.params['rbgreater']['value']
        seeds = self.params['seeds']['value']
        seeds_ = map(int,seeds.strip().split(","))

        g = None
        ind = None
        w = None
        m = None
        datahist_global = None

        xgrid = sgems.get_property(head_grid,"_X_")
        ygrid = sgems.get_property(head_grid,"_Y_")
        zgrid = sgems.get_property(head_grid,"_Z_")
        grid_value = sgems.get_property(head_grid,head_property)
        r = sgems.get_dims(head_grid)

        xgrid1 = sgems.get_property(tail_grid,"_X_")
        ygrid1 = sgems.get_property(tail_grid,"_Y_")
        zgrid1 = sgems.get_property(tail_grid,"_Z_")
        n_ = sgems.get_property(tail_grid,tail_property)

        v_d = [] #value id
        l_c = [] #list cut off
        l_2 = [] #second list
        l_g = [] #global list reference_value
        h = 1

        n_rows = r[0]#30#len(ygrid)
        n_cols = r[1]#25#len(xgrid)
        n_levels = r[2]
        print "berthin", r, "n_", n_rows, n_cols, n_levels
        K = len(seeds_)
        matrix_dim = n_rows, n_cols

        def cut_off_(more_equal, less_equal):
            if (more_equal == '1' and less_equal == '0'):
                def w_o(cut_off): #define waste ore in 0 and 1
                    for i in grid_value:
                        if i >= cut_off:
                            l_c.append(i)
                        else:
                            l_c.append(0)
                    return sgems.set_property(head_grid,new_prop,l_c),l_c
                    #return l_c
                w_o(cut_off)
            elif (more_equal == '0' and less_equal == '1'):
                def w_o_(cut_off): #define waste as 0 and ore as grade
                    for i in grid_value:
                        if i <= cut_off:
                            l_c.append(i)
                        else:
                            l_c.append(0)
                    return sgems.set_property(head_grid,new_prop,l_c),l_c
                    #return l_c
                w_o_(cut_off)
            elif (more_equal == '0' and less_equal == '0'):
                print "Select option of cut off"
                self.dict_gen_params['execution_status'] = "ERROR"
            else:
                print "Error in execution type parameters"
                self.dict_gen_params['execution_status'] = "ERROR"

        grid1 = cut_off_(more_equal, less_equal)
        gridnt = l_c
        grid = np.array(gridnt).reshape(n_rows, n_cols, n_levels)

        plt_imshowN([grid, grid == 0, grid, grid == 0],
                    [0, 0, 1, 1], [dict(cmap='jet', interpolation='nearest')] * 4, dim=(2,2))

        """
        def cut_off_(grid_value, cut_off, more_equal, less_equal):
            more_equal, less_equal = map(int, [more_equal, less_equal])
            if more_equal + less_equal != 1:
                print 'Error, wrong cut_off values'
            else:
                if more_equal:
                    threshold = grid_value >= cut_off
                else:
                    threshold = grid_value <= cut_off
                return threshold * grid_value
        grid2 = cut_off_(grid_value, cut_off, more_equal, less_equal)
        print type(grid2), len(grid2)
        sgems.set_property(head_grid,new_prop,grid2)#print grid#plt_imshowN([head_grid, grid], [dict(cmap='jet', interpolation='nearest'), dict(cmap='jet', interpolation='nearest')], (1,2))
        """

        def quantil(p):
            return np.percentile(p, np.arange(0,100,5))

        def compare(s,h):
            c=[]
            for i in s:
                for j in h:
                    if i+j > 0:
                        c.append(((i-j)**2)/(i+j))
            d= sum(c)
            return d

        def run_union_find(painter, K, iter_label, max_steps):
            fake_painter = copy.deepcopy(painter)
            fake_painter.paint_Kregions(K, iter_label, max_steps)
            y0 = np.ndarray.tolist((fake_painter.connected_processed_map * fake_painter.grid).ravel())
            c_global= quantil(gridnt)
            c_list=quantil(y0)
            return compare(c_global,c_list), fake_painter

        def search_best_run(painter, K, n_iterations, max_steps):
            best_cost = 10000000000
            best_painter = None
            for iter_label in xrange(n_iterations):
                cost, fake_painter = run_union_find(painter, K, iter_label, max_steps)
                print 'iter', iter_label, 'cost', cost
                if cost < best_cost:
                    best_cost = cost
                    best_painter = fake_painter
            assert(best_painter != None)
            return (best_cost, best_painter)

        def repeat_all(grid,K, painter, n_iterations, n_steps_list, n_runs):
            for i, max_steps in zip(xrange(n_runs), n_steps_list):
                print i, "run"
                cost, painter = search_best_run(painter, K, n_iterations, max_steps)
                #plt_imshowN([painter.connected_map, painter.connected_map_ordem], [0, 0], [dict(cmap='jet', interpolation='nearest'), dict(cmap='jet')], dim=(1,2))
                plt_imshowN([painter.connected_map * (painter.grid != 0), painter.connected_map_ordem, painter.connected_map * (painter.grid != 0), painter.connected_map_ordem],
                            [0, 0, 1, 1], [dict(cmap='jet', interpolation='nearest')] * 4, dim=(2,2))
                #np.reshape(grid, (1,np.product(grid.shape)))
                maps_seed = (painter.connected_map * painter.connected_map)
                maps_seed_ = np.ndarray.tolist(maps_seed.ravel()) #maps__ = np.array(maps_).T
                sgems.set_property(head_grid,new_prop1+'-'+str(i),maps_seed_)
                maps_ = painter.connected_map_ordem
                maps__ = np.ndarray.tolist(maps_.ravel())
                #
                sgems.set_property(head_grid,new_prop2+'-'+str(i),maps__)
                ee = np.ma.masked_array(grid,mask=maps_==0)
                ee_ = ee.filled(0)
                ee__ = np.ndarray.tolist(ee_.ravel())
                sgems.set_property(head_grid,"sect_maps"+'-'+str(i),ee__)


        uf = UnionFind(n_rows * n_cols * n_levels)
        guto = Painter(seeds_, grid, uf, K)
        repeat_all(grid, K, painter = guto, n_iterations = simulation, n_steps_list = [fases] * runs, n_runs = runs)#cambiando

        return True

    def finalize(self):
        return True

    def name(self):
        return "pdfsequential12"

################################################################################
def get_plugins():
    return ["pdfsequential12"]

################################################################################
from pylab import plt
import numpy as np

def plt_imshowN(imgs, levels, params = None, dim = (1,2)):
    fig = plt.figure()
    for i in xrange(1, np.multiply(*dim) + 1):
        p = fig.add_subplot(dim[0], dim[1], i)
        if params is not None and params[i-1] is not None:
            p.imshow(imgs[i-1][..., levels[i - 1]], **params[i-1])
        else:
            p.imshow(imgs[i-1][..., levels[i - 1]])
    plt.show()

def plotSurface(M, level):
    n_elts_per_levels = self.n_rows * self.n_cols
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    xx, yy = np.meshgrid(range(self.n_cols), range(self.n_rows))
    surf = ax.plot_surface(xx, yy,  M[..., level], cmap = 'terrain', cstride = 1, rstride = 1, linewidth=0, antialiased=False)
    plt.show()

# plt_imshow([im1, im2], [None, {'cmap':'gray'}], dim=(1,2))"""
