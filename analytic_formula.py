#A reminder: the camera is a long way away in the z direction

import pylab, sys, os
import numpy as np
from scipy.special import *






if dictionary_SI['analytic'] ==0:
    None
    #Do NOT USE analytic = 0. This is reserved for no analytic model
    
elif dictionary_SI['analytic'] == 1:
    print "Analytic Sphere"
    def theory(QRadius):
        QR = QRadius*dictionary_SI['radius_1']
        def function(x):
            return ((3/(QR**3))*(np.sin(QR)-(QR*np.cos(QR))))**2
        b= np.asarray([function(QR[x]) if QR[x]!=0 else function(QR[x+1]) for x in range(len(QR))])
        return b


elif dictionary_SI['analytic'] == 2:
    print "Analytic Cyliner"
    def theory(distance_from_origin):
        QR = distance_from_origin*dictionary_SI['radius_1']
        def function(x):
            return (jv(1.,x)/(x))**2
        cylinder_intensity = np.asarray([function(QR[x]) if QR[x]!=0 else function(QR[x+1]) for x in range(len(QR))])
        return cylinder_intensity


elif dictionary_SI['analytic'] == 3:
    print "Core Shell analytic model"
    def theory(QRadius):
        def function(x):
            return ((dictionary_SI['rho_1']-dictionary_SI['rho_2'])*dictionary_SI['radius_2']*jv(1.,x*dictionary_SI['radius_2'])/x +dictionary_SI['rho_2']*dictionary_SI['radius_1']*jv(1., x*dictionary_SI['radius_1'])/x)
        b= np.asarray([function(QRadius[x]) if QRadius[x]!=0 else function(QRadius[x+1]) for x in range(len(QRadius))])
        return b**2
#
#__________  radius_2
#rho_1     |
#          |
#          |
#          |
#          |         radius_1
#          |         ________________rho = 0
#          |        |
#          |        |
#          |        |
#          |________|rho_2 (note: rho_2 should be entered as a negative number if desired)
#
#




elif dictionary_SI['analytic'] == 4:
    print "analytic gaussian cylinder"
    print "The standard deviation is given from the radius"
    def theory(QRadius):
        QR = QRadius*dictionary_SI['radius_2']
        def function(x):
            return np.exp(-(x/4)**2)
        b= np.asarray([function(QR[x]) if QR[x]!=0 else function(QR[x+1]) for x in range(len(QR))])
        return b**2



def theory_csv():
        global dictionary_SI
        QSize = dictionary_SI['QSize']
        pixels = dictionary_SI['pixels']
        EHC = dictionary_SI['EHC']
        z_dim = dictionary_SI['z_dim']
        x_theta = dictionary_SI['x_theta']

        Q = np.array([[[row*QSize/pixels-0.5*QSize, col*QSize/pixels-0.5*QSize] for col in range(int(pixels))] for row in range(int(pixels))])

        QRArray = np.array(Q[:,:,0]**2+Q[:,:,1]**2)**0.5
        alpha = np.angle( Q[:,:,0]+1j*Q[:,:,1] )
        x_theta = dictionary_SI['x_theta'] #gamma
        theta = np.arcsin(np.sum(Q**2, axis = 2)**0.5/(2*EHC))

        QZArray = np.sin(theta)*EHC*(   np.cos(alpha)*np.sin(x_theta)*np.cos(theta)-np.cos(x_theta)*np.sin(theta)   )
        Intensity = np.asarray([theory(QRadius*np.cos(x_theta)) for QRadius in QRArray])
        L = np.array([[2*np.sin(a*z_dim/2.)/a if a!=0 else z_dim for a in b] for b in QZArray])
        
        Intensity = np.asarray([theory(np.cos(x_theta)*np.sum(QRadius**2, axis = 1)**0.5) for QRadius in Q])
        Intensity = Intensity*(L**2)
        return Intensity/np.sum(Intensity)

