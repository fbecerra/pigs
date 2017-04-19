import matplotlib
matplotlib.use('TkAgg')
from yt.config import ytcfg; ytcfg["yt","serialize"] = "True"
from yt.mods import *
from yt.analysis_modules.star_analysis.api import *
from yt.analysis_modules.level_sets.api import *
from fields import *
from constants import *
import numpy
import math
import re
import random
import scipy as sc
from scipy import interpolate
from scipy import ndimage
import pylab as pl
pl.ion()
from matplotlib import rc
fontsize=14
rc('text', usetex=True)
rc('font', **{'family':'serif','serif':'Computer Modern Roman', 'size':fontsize})
rc('axes', labelsize=fontsize)
rc('legend', fontsize=fontsize, numpoints=1, frameon=False)
rc('xtick', labelsize=fontsize)
rc('ytick', labelsize=fontsize)
rc('lines', lw=1.5, mew=0.3)
rc('grid', linewidth=0.5)
