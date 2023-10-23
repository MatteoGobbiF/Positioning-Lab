import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import sys

sys.path.append("..")
from utils import *

# Ionospheric correction parameters:
alpha = [7.4506*10**(-9), 1.4901*10**(-8), -5.9605*10**(-8), -1.1921*10**(-7)]
beta = [9.2160*10**(4), 1.3107*10**(5), -6.5536*10**(4), -5.2429*10**(5)]

