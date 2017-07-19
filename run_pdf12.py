import pdfsequential12
import numpy as np
import utils
reload(utils)
reload(pdfsequential12)
from utils import Map3D
from pylab import plt

#n_rows, n_cols, n_levels = 100, 200, 3
#grid = np.array([np.random.randint(0, 10000) for i in xrange(n_rows * n_cols * n_levels)]).reshape(n_rows, n_cols, n_levels)

mapa = Map3D()
mapa.readFromFile('3dproob3.txt')
mapa.convertCoordinates()
n_rows, n_cols, n_levels = mapa.n_rows, mapa.n_cols, mapa.n_depth
grid = mapa.toMatrix()
print grid.shape
plt.imshow(grid[..., 0])
plt.show()
plt.imshow(grid[..., 1])
plt.show()
seeds = [(0,5,0), (10,0,0), (0,15,1), (15,10,1)]
limits = [n_rows, n_cols, n_levels]
assert (all([grid[(ith_seed)] for ith_seed in seeds]))
seeds = [5, 10*n_cols, n_rows*n_cols+15, n_rows*n_cols+15*n_cols+10]
pdf = pdfsequential12.pdfsequential12()
params = dict(
  cut_off = dict(value = 1),
  fases = dict(value = 100),
  simulation = dict(value = 2),
  runs = dict(value = 10),
  seeds = dict(value = ','.join(map(str, seeds))),
  rsless = dict(value = '0'),
  rbgreater = dict(value = '1'),
  #head_prop = dict(grid = grid[..., 0].reshape(29, 25, 1), property = None),
  head_prop = dict(grid = grid, property = None),
  tail_prop = dict(grid = None, property = None),
  new_prop = dict(value = 'prop'),
  new_prop1 = dict(value = 'prop1'),
  new_prop2 = dict(value = 'prop2'),
)

pdf.initialize(params=params)
pdf.execute()
