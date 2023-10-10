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

def compute_rotation_matrix(xi, eta, alpha):
    Rx = np.array([[1, 0, 0],
                   [0, math.cos(xi), -math.sin(xi)],
                   [0, math.sin(xi), math.cos(xi)]])
    Ry = np.array([[math.cos(eta), 0, -math.sin(eta)],
                   [0, 1, 0],
                   [math.sin(eta), 0, math.cos(eta)]])
    Rz = np.array([[math.cos(alpha), math.sin(alpha), 0],
                   [-math.sin(alpha), math.cos(alpha), 0],
                   [0, 0, 1]])
    
    return (Rz@Ry@Rx)

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