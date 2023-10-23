import math
import numpy as np

def deg2rad(deg):
    return deg/180*math.pi

def sex2deg(sex):
    return sex[0] + sex[1]/60 + sex[2]/3600

def deg2sex(deg):
    degrees = np.fix(deg)
    primes = np.fix((deg-degrees)*60)
    seconds = ((deg-degrees)*60-primes)*60
    return [degrees, primes, seconds]

def rad2deg(rad):
    return rad*180/math.pi

def sex2rad(sex):
    return deg2rad(sex2deg(sex))

def rad2sex(rad):
    return deg2sex(rad2deg(rad))

def compute_Rn (latitude, a = 6378137, e = 0.081819191042832):
    return a/math.sqrt(1- (e**2)*math.sin(latitude)**2)

def gg2gc (gg, a = 6378137, e = 0.081819191042832):
    Rn = compute_Rn(gg[0])

    x = (Rn + gg[2])*math.cos(gg[0])*math.cos(gg[1])
    y = (Rn + gg[2])*math.cos(gg[0])*math.sin(gg[1])
    z = (Rn*(1-e**2) + gg[2])*math.sin(gg[0])

    return [x, y, z]

def compute_R_0(lat, long):
    return np.array([[-math.sin(long), math.cos(long), 0], 
         [-math.sin(lat)*math.cos(long), -math.sin(lat)*math.sin(long), math.cos(lat)],
         [math.cos(lat)*math.cos(long), math.cos(lat)*math.sin(long), math.sin(lat)]])

def lc2gc (lc, Xogc, lat, long):
    
    return Xogc + compute_R_0(lat, long).T.dot(lc)

def compute_Rx(a):
    return np.array([[1, 0, 0],
                    [0, math.cos(a), -math.sin(a)],
                    [0, math.sin(a), math.cos(a)]])

def compute_Ry(a):
    return  np.array([[math.cos(a), 0, -math.sin(a)],
                     [0, 1, 0],
                     [math.sin(a), 0, math.cos(a)]])

def compute_Rz(a):
    return np.array([[math.cos(a), math.sin(a), 0],
                     [-math.sin(a), math.cos(a), 0],
                     [0, 0, 1]])

def compute_rotation_matrix(xi, eta, alpha):
       return (compute_Rz(alpha)@compute_Ry(eta)@compute_Rx(xi))

def gc2gg(gc, a = 6378137, b = 6356752.3141, e = 0.081819191042832): #If not specified use GRS80
    e2_b = (a**2-b**2)/b**2
    r = math.sqrt(gc[0]**2 + gc[1]**2)
    psi = math.atan2(gc[2], r*math.sqrt(1-e**2))

    phi = math.atan2(gc[2]+e2_b*b*math.sin(psi)**3, r-e**2*a*math.cos(psi)**3)
    lamb = math.atan2(gc[1], gc[0])

    Rn = compute_Rn(phi, a, e)

    h = r/math.cos(phi)-Rn

    return [phi, lamb, h]

def propagate_cov(C, R):
    return (R @ C @ R.T)

def compute_clock_offset(epoch, dt0, a, b):
    return (dt0 + a*epoch + b*epoch**2)

def estimate_eta(prev, threshhold, Mt, e):
    eta = Mt+e*math.sin(prev)
    if(abs(eta-prev)<threshhold):
        return eta
    else:
        return estimate_eta(eta, threshhold, Mt, e)

def ORS2ITRF(ors, Omega, i, w):
    return compute_Rz(-Omega)@compute_Rx(-i)@compute_Rz(-w)@ors.T

def iono_correction(lat, lon, az, el, time, ionoparams):
    # Initialization
    c = 299792458

    # Ionospheric parameters
    a0, a1, a2, a3, b0, b1, b2, b3 = ionoparams

    # Elevation from 0 to 90 degrees
    el = np.abs(el)

    # Conversion to semicircles
    lat = lat / 180
    lon = lon / 180
    az = az / 180
    el = el / 180

    psi = (0.0137 / (el + 0.11)) - 0.022

    phi = lat + psi * np.cos(az * np.pi)

    if phi > 0.416:
        phi = 0.416
    elif phi < -0.416:
        phi = -0.416

    # Geodetic longitude of the earth projection of the ionospheric intersection point
    lambda_ = lon + (psi * np.sin(az * np.pi)) / np.cos(phi * np.pi)

    # Geomagnetic latitude of the earth projection of the ionospheric intersection point
    ro = phi + 0.064 * np.cos((lambda_ - 1.617) * np.pi)

    # Local time in seconds
    t = lambda_ * 43200 + time

    if t >= 86400:
        t -= 86400
    elif t < 0:
        t += 86400

    # Obliquity factor
    f = 1 + 16 * (0.53 - el) ** 3

    a = a0 + a1 * ro + a2 * ro ** 2 + a3 * ro ** 3
    a = np.where(a < 0, 0, a)

    p = b0 + b1 * ro + b2 * ro ** 2 + b3 * ro ** 3
    # p = np.where(p < 72000, 72000, p)

    x = (2 * np.pi * (t - 50400)) / p

    # Ionospheric delay
    if abs(x) < 1.57:
        delay = c * f * (5e-9 + a * (1 - (x ** 2) / 2 + (x ** 4) / 24))

    else: 
        delay = c * f * 5e-9

    return delay
