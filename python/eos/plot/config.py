import matplotlib
from matplotlib import rcParams

matplotlib.use('pgf')

# set some default values for plotting
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.serif'] = 'Computer Modern Sans serif'
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['font.weight'] = 400

matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['axes.linewidth'] = 1

matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.pad_inches'] = 0.1

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

matplotlib.rcParams['text.usetex'] = True

matplotlib.rcParams['pgf.preamble'] = [
    r'\usepackage[hidelinks]{hyperref}',
    r'\usepackage{amsmath}',
    r'\usepackage{xcolor}'
]
