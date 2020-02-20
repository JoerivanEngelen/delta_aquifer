# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 15:14:18 2020

Utility functions to calculate distances between two points on the geoid.
Adapted from the functions written by Rens van Beek, Utrecht University.

@author: engelen


"""

import numpy as np

deg2Rad= np.pi/180.

def getArcDistance(latA, lonA, latB, lonB, radius= 6371221.3, testVerbose= False):
  '''Computes the distance between two points, positioned by \
their geographic coordinates along the surface of a perfect sphere \
in units given by the radius used. Input variables include:\n
- latA, latB: the latitude of the points considered in decimal degrees,\n
- lonA, lonB: the longitude of the points considered in decimal degrees,\n
- radius: the radius of the sphere in metres, set by default to \
that of Earth (6371221.3 m).'''
        #-make arrays if needed
  if isinstance(latA,float):
                latA= np.array(latA)
  if isinstance(lonA,float):
                lonA= np.array(lonA)
  if isinstance(latB,float):
                latB= np.array(latB)
  if isinstance(lonB,float):
                lonB= np.array(lonB)
  #-pad latitudes, longitudes
  if latA.size != latB.size:
    latA= np.ones(latB.shape)*latA
  if lonA.size != lonB.size:
    lonA= np.ones(lonB.shape)*lonA
  #-convert all coordinates to radians
  pA= latA*deg2Rad; pB= latB*deg2Rad
  lA= lonA*deg2Rad; lB= lonB*deg2Rad
  dL= lB-lA; dP= pB-pA
  #-set arcDist default to zero and create mask of cells to be processed
  a= np.sin(0.5*dP)*np.sin(0.5*dP)+\
    np.cos(pA)*np.cos(pB)*\
    np.sin(0.5*dL)*np.sin(0.5*dL)
  arcDist= 2*np.arctan2(a**0.5,(1-a)**0.5)
  arcDist*= radius
  if testVerbose:
    print(' * along an ideal sphere of radius %f, the distance between the points at lat, lon' \
  '%f, %f and %f, %f respectively amounts to %f' %\
      (radius, latA[arcDist == arcDist.max()][0].astype(float), lonA[arcDist == arcDist.max()][0].astype(float),\
       latB[arcDist == arcDist.max()][0].astype(float),\
       lonB[arcDist == arcDist.max()][0].astype(float), arcDist[arcDist == arcDist.max()][0].astype(float)))
  #-return arcDist
  return arcDist

def getAzimuth(latA, lonA, latB, lonB, radius= 6371221.3, testVerbose= False):
  '''Returns the array of the azimuth between two points, positioned by \
their geographic coordinates along the surface of a perfect sphere \
in units given by the radius used. Input variables include:\n
- latA, latB: the latitude of the points considered in decimal degrees,\n
- lonA, lonB: the longitude of the points considered in decimal degrees,\n
- radius: the radius of the sphere in metres, set by default to \
that of Earth (6371221.3 m).
NOTE: azimuth is computed relative to the point specified by A and positive CW from N'''
        #-make arrays if needed
  if isinstance(latA,float):
                latA= np.array(latA)
  if isinstance(lonA,float):
                lonA= np.array(lonA)
  if isinstance(latB,float):
                latB= np.array(latB)
  if isinstance(lonB,float):
                lonB= np.array(lonB)
  #-pad latitudes, longitudes
  if latA.size != latB.size:
    latA= np.ones(latB.shape)*latA
  if lonA.size != lonB.size:
    lonA= np.ones(lonB.shape)*lonA
  pA= latA*deg2Rad; pB= latB*deg2Rad
  lA= lonA*deg2Rad; lB= lonB*deg2Rad
  dL= lB-lA; dP= pB-pA
  #-azimuth
  x= np.sin(dL)*np.cos(pB)
  y= np.cos(pA)*np.sin(pB)-np.sin(pA)*np.cos(pB)*np.cos(dL)
  azimuth= deg2Rad**-1*np.arctan2(x, y)+360.
  azimuth= azimuth % 360.

  #-return azimuth
  return azimuth

