import matplotlib.pyplot as plt
import math
import sys
sys.path.append("..")

#Import data
from utils import *

dt0 = -7.661711424589e-05
dt1 = -3.183231456205e-12
dt2 =  0.000000000000e+00
sqrt_a = 5.153650835037e+03 # (sqrt of meters)
a = sqrt_a**2
e = 3.841053112410e-03
M0 = 1.295004883409e+00 #(radians)
Omega_0 = -2.241692424630e-01 #(radians)
Omega_dot = -8.386063598924e-09 #(radians/sec)
i_0 = 9.634782624741e-01 #(radians)
i_dot = -7.286017777600e-11 #(radians/sec)
w_0 = 9.419793734505e-01 #(radians)
w_dot = 0.0 #(radians/sec)
GMe = 3.986005e+14 #(m3/s2)
Omega_E_dot = 7.2921151467e-05 #(radians)

t_start = 0
t_end = 86399
t_step = 30
t = list(range(t_start, t_end, t_step))

clock_offsets = []

for epoch in t:
    clock_offsets.append(compute_clock_offset(epoch, dt0, dt1, dt2))

# plot
fig, ax = plt.subplots(figsize=(10,6))
ax.set(xlabel='Seconds in a day', ylabel='Clock-offset', title = 'Clock offset vs time')
ax.plot(t, clock_offsets, '-', color='blue')
plt.show()

coord_ORS = []
coord_ITRF_cart = []
coord_ITRF_geo = []

for Dt in t:
    n = math.sqrt(GMe/a**3)
    Mt = M0 + n*(Dt-t_start)
    eta = estimate_eta(0, 10e-10, Mt, e)
    psi = 1/math.tan((math.sqrt(1-e**2)*math.sin(eta))/(math.cos(eta)-e))
    r = (a*(1-e**2))/(1+e*math.cos(psi))
    w = w_0 + w_dot*(Dt-t_start)
    i = i_0 + i_dot*(Dt-t_start)
    Omega = Omega_0 + (Omega_dot - Omega_E_dot) * (Dt - t_start)

    x = r * math.cos(psi)
    y = r * math.sin(psi)
    z = 0

    coord_ORS.append([x,y,z])
    coord_ITRF_cart.append(ORS2ITRF(np.array([x,y,z]),Omega,i,w))
    coord_ITRF_geo.append(gc2gg(coord_ITRF_cart[-1]))


import pandas as pd
import geopandas as gpd
coord_ITRF_cart = np.array(coord_ITRF_cart).T
print(coord_ITRF_cart)
# Realize a dataframe containing satellite coordinates
df = pd.DataFrame()
df['time'] = t
df['Latitude'] = coord_ITRF_cart[0,:]
df['Longitude'] = coord_ITRF_cart[1,:]

# Transform the DataFrame in a GeoDataFrame
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude), crs = 3857)

# Load the basemap
world = gpd.read_file('world\world.shp')

# Plot the trajectory with world basemap
fig, ax = plt.subplots (figsize = (15,15))
world.plot(ax=ax)
ax.set(xlabel='Longitude', ylabel='Latitude', title='Satellite daily trajectory')
gdf.plot(ax = ax, marker='o', color='red')

plt.show()




