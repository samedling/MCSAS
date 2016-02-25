#A reminder: the camera is a long way away in the z direction

import pylab, sys, os
import numpy as np
from scipy.special import *

import global_vars as g

#    def theory(QRadius):
#        QR = QRadius*g.dictionary_SI['radius_1']
#        def function(x):
#            return ((3/(QR**3))*(np.sin(QR)-(QR*np.cos(QR))))**2
#        b= np.asarray([function(QR[x]) if QR[x]!=0 else function(QR[x+1]) for x in range(len(QR))])
#        return b
def t1sphere(QRadius):     #TODO: check this; it no longer returns a 2D array instead of a 1D one, but it might not be correct
    QR = QRadius*g.dictionary_SI['radius_1']
    def function(x):
        return ((np.sin(x)-x*np.cos(x))/x**3)**2
    return np.asarray([function(QR[x]) if QR[x]!=0 else function(QR[x+1]) for x in range(len(QR))])

def t2cylinder(distance_from_origin):
    QR = distance_from_origin*g.dictionary_SI['radius_1']
    def function(x):
        return (jv(1.,x)/(x))**2
    return np.asarray([function(QR[x]) if QR[x]!=0 else function(QR[x+1]) for x in range(len(QR))])

def t3coreshell(QRadius):
    def function(x):
        return ((g.dictionary_SI['rho_1']-g.dictionary_SI['rho_2'])*g.dictionary_SI['radius_2']*jv(1.,x*g.dictionary_SI['radius_2'])/x +g.dictionary_SI['rho_2']*g.dictionary_SI['radius_1']*jv(1., x*g.dictionary_SI['radius_1'])/x)
    b= np.asarray([function(QRadius[x]) if QRadius[x]!=0 else function(QRadius[x+1]) for x in range(len(QRadius))])
    return b**2

def t4gaussian(QRadius):
    QR = QRadius*g.dictionary_SI['radius_2']
    def function(x):
        return np.exp(-(x**2)/2.)
    b= np.asarray([function(QR[x]) if QR[x]!=0 else function(QR[x+1]) for x in range(len(QR))])
    return b



def theory(Q_radius):
   if g.dictionary_SI['analytic'] == 1:
      return t1sphere(Q_radius)
   elif g.dictionary_SI['analytic'] == 2:
      return t2cylinder(Q_radius)
   elif g.dictionary_SI['analytic'] == 3:
      return t3coreshell(Q_radius)
   elif g.dictionary_SI['analytic'] == 4:
      return t4gaussian(Q_radius)


#if g.dictionary_SI['analytic'] ==0:
#    None
#    #Do NOT USE analytic = 0. This is reserved for no analytic model
#    
#elif g.dictionary_SI['analytic'] == 1:
#    print "Analytic Sphere"
#    def theory(QRadius):
#        QR = QRadius*g.dictionary_SI['radius_1']
#        def function(x):
#            return ((3/(QR**3))*(np.sin(QR)-(QR*np.cos(QR))))**2
#        b= np.asarray([function(QR[x]) if QR[x]!=0 else function(QR[x+1]) for x in range(len(QR))])
#        return b
#
#
#elif g.dictionary_SI['analytic'] == 2:
#    print "Analytic Cylinder"
#    def theory(distance_from_origin):
#        QR = distance_from_origin*g.dictionary_SI['radius_1']
#        def function(x):
#            return (jv(1.,x)/(x))**2
#        cylinder_intensity = np.asarray([function(QR[x]) if QR[x]!=0 else function(QR[x+1]) for x in range(len(QR))])
#        return cylinder_intensity
#
#
#elif g.dictionary_SI['analytic'] == 3:
#    print "Core Shell analytic model"
#    def theory(QRadius):
#        def function(x):
#            return ((g.dictionary_SI['rho_1']-g.dictionary_SI['rho_2'])*g.dictionary_SI['radius_2']*jv(1.,x*g.dictionary_SI['radius_2'])/x +g.dictionary_SI['rho_2']*g.dictionary_SI['radius_1']*jv(1., x*g.dictionary_SI['radius_1'])/x)
#        b= np.asarray([function(QRadius[x]) if QRadius[x]!=0 else function(QRadius[x+1]) for x in range(len(QRadius))])
#        return b**2
##
##__________  radius_2
##rho_1     |
##          |
##          |
##          |
##          |         radius_1
##          |         ________________rho = 0
##          |        |
##          |        |
##          |        |
##          |________|rho_2 (note: rho_2 should be entered as a negative number if desired)
##
##
#
#
#
#
#elif g.dictionary_SI['analytic'] == 4:
#    print "analytic gaussian cylinder"
#    print "The standard deviation is given from the radius"
#    def theory(QRadius):
#        QR = QRadius*g.dictionary_SI['radius_2']
#        def function(x):
#            return np.exp(-(x/4)**2)
#        b= np.asarray([function(QR[x]) if QR[x]!=0 else function(QR[x+1]) for x in range(len(QR))])
#        return b**2






def theory_csv():
    g.dprint('Calling theory_csv.')
    QSize = g.dictionary_SI['QSize']
    xpixels,ypixels = g.dictionary_SI['pixels'].split()
    EHC = g.dictionary_SI['EHC']
    z_dim = g.dictionary_SI['z_dim']

    Q = np.array([[[row*QSize/float(xpixels)-0.5*QSize, col*QSize/float(ypixels)-0.5*QSize, 2*EHC*np.sin((((row-0.5*float(ypixels))**2 + (col-0.5*float(xpixels))**2)**0.5)*QSize/float(xpixels)/2/EHC)**2] for col in range(int(ypixels))] for row in range(int(xpixels))])

    QRArray = np.array(Q[:,:,0]**2+Q[:,:,1]**2+Q[:,:,2]**2)**0.5#absolute value of the q vector

    if g.dictionary_SI['analytic']!=1:#all cylindrical objects run the following.
        alpha = np.angle( Q[:,:,0]+1j*Q[:,:,1] )
        x_theta = g.dictionary_SI['x_theta'] #gamma
        theta = 0.5*np.arctan(QRArray/EHC)
        QZArray = QRArray*(np.cos(alpha)*np.sin(x_theta)*np.cos(theta)-np.cos(x_theta)*np.sin(theta))#The modified Qz, because of the rotations

        #This is here as it took me a long time to figure out what was going on with this.
        #For small angle approx/etc, the constant out the front of the form factor is proportional to L.
        #If you don't make this approximation, instead you get (2/qz)*sin(qz*L/2) which approaches L as qz goes to zero. 
        #This is why we multiply by L**2 (squaring is because we have already squared the normal form factor)
        #Then we want to multiply by the normal QRArray, as we are now interested in the total Z for the form factor, not just the z component.

        if g.dictionary_SI['Qz']==0:#No small angle approximation
           L= np.array([[2/(QZArray[x,y])*np.sin(QZArray[x,y]*z_dim/2.) if (QZArray[x,y])!=0 else z_dim for x in range(int(xpixels))] for y in range(int(ypixels))])
        else:
           L= np.array([[1/(QRArray[x,y]*np.cos(alpha[x,y])*np.sin(x_theta))*np.sin(QRArray[x,y]*np.cos(alpha[x,y])*np.sin(x_theta)*z_dim/2.) if (QRArray[x,y]*np.cos(alpha[x,y])*np.sin(x_theta))!=0 else z_dim for x in range(int(xpixels))] for y in range(int(ypixels))])
        Intensity = np.asarray([theory(QRadius) for QRadius in QRArray]*(L**2))
        return Intensity/np.sum(Intensity)
    else:
           Intensity = np.asarray([theory(QRadius) for QRadius in QRArray])
           return Intensity/np.sum(Intensity)






