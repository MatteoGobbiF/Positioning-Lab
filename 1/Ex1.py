import numpy as np
import math

import sys
sys.path.append("..")

from utils import *

xO_gg = [sex2rad([44,23,24.000]), sex2rad([8,56,20.000]), 70]

xO_gc = gg2gc(xO_gg)

xA_ll = [0, 30, 0]
xB_ll = [0, -30, 0]
xC_ll = [200, 0, 0]

xi = sex2rad([0, 0, 10.23])
eta = sex2rad([0, 0, 9.5])
alpha = sex2rad([30, 27, 18])

R_lc2ll = compute_rotation_matrix(xi, eta, alpha)
R_ll2lc = R_lc2ll.T

#Local cartesian

xA_lc = R_ll2lc.dot(xA_ll)
xB_lc = R_ll2lc.dot(xB_ll)
xC_lc = R_ll2lc.dot(xC_ll)
#print(round(xA_lc[0], 3))

#Global cartesian

xA_gc = lc2gc(xA_lc, xO_gc, xO_gg[0], xO_gg[1])
xB_gc = lc2gc(xB_lc, xO_gc, xO_gg[0], xO_gg[1])
xC_gc = lc2gc(xC_lc, xO_gc, xO_gg[0], xO_gg[1])

#Conversion to ETRF on the website
#Insert ETRF values manually

xA_ETRF = [4509854.8133, 709344.7329, 4439228.7611]
xB_ETRF = [4509885.8305, 709380.3979, 4439191.8018]
xC_ETRF = [4509773.4717, 709521.8492, 4439282.7125]

#Global geodetic

xA_gg_rad = gc2gg(xA_ETRF)
xB_gg_rad = gc2gg(xB_ETRF)
xC_gg_rad = gc2gg(xC_ETRF)

#Conversion to sex

xA_gg = [rad2sex(xA_gg_rad[0]), rad2sex(xA_gg_rad[1]), xA_gg_rad[2]]
xB_gg = [rad2sex(xB_gg_rad[0]), rad2sex(xB_gg_rad[1]), xB_gg_rad[2]]
xC_gg = [rad2sex(xC_gg_rad[0]), rad2sex(xC_gg_rad[1]), xC_gg_rad[2]]

#Covariance propagation

sigmaA = 0.02
sigmaB = 0.02
sigmaC = 0.1
sigmaO = 0.1

C_ll_A = np.diag([sigmaA**2, sigmaA**2, sigmaA**2])
C_ll_B = np.diag([sigmaB**2, sigmaB**2, sigmaB**2])
C_ll_C = np.diag([sigmaC**2, sigmaC**2, sigmaC**2])

C_gg_O = np.diag([sigmaO**2, sigmaO**2, sigmaO**2])

#This get the wrong result for some reason
C_lc_A = propagate_cov(C_ll_A, R_ll2lc)
C_lc_B = propagate_cov(C_ll_B, R_ll2lc)
C_lc_C = propagate_cov(C_ll_C, R_ll2lc)

print(C_lc_A)

R_0 = compute_R_0(xO_gg[0], xO_gg[1])

C_A = propagate_cov((C_gg_O + propagate_cov(C_lc_A, R_0.T)), R_0)
C_B = propagate_cov((C_gg_O + propagate_cov(C_lc_B, R_0.T)), R_0)
C_C = propagate_cov((C_gg_O + propagate_cov(C_lc_C, R_0.T)), R_0)

#Get the standard deviation in East, North, Up in cm for each point

deviationA = 100*(np.sqrt(np.diag(C_A)))
deviationB = 100*(np.sqrt(np.diag(C_B)))
deviationC = 100*(np.sqrt(np.diag(C_C)))




