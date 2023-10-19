#! /usr/bin/env python
import sys, os
import shutil
import subprocess
from datetime import datetime
import time as timeit
import pickle
import glob
import copy
import numpy as np
import scipy.stats as stats
import scipy.optimize as optimize
from scipy.signal import argrelmin
from scipy.integrate import simps
from scipy.spatial import Voronoi, Delaunay, ConvexHull
from sklearn.decomposition import PCA
import pytraj as pt
import parmed as pmd
from pymbar import mbar
import water_properties as wp

import surface_library as sl

# testing sortlib and waterlib versions in "ProteinDev" project 
sys.path.append('/home/drobins/ProteinDev/fortran')
import sortlib
import waterlib as wl

# only need this for python2, update when I get this framework functional on python3 
import matplotlib
try:
  os.environ["DISPLAY"]
except KeyError:
  shsowPlots = False
  matplotlib.use('Agg')

import matplotlib.pyplot as plt


Usage = """A library of classes and functions to compute the system-average and local order parameters for water structure.
"""

### I think this is all I can use the TrajObject class for without making things overly complex, but it may be useful to implement additional classes to deal with the finer details...
class TrajObject:
  """A class to manage trajectory information more compactly than implementing parmed and pytraj each time.
     Attributes:
                 topFile - topology file for this system
                 trajFile - trajectory file for this system
                 stride - skip every "stride" steps of the trajectory (int). Default=1
                 solResName - string defining the residue name for the non-water cosolvent. Default='(!:WAT)'
                 watResName - string defining the residue name for the water. Default='(:WAT)'
  """
  def __init__(self, topFile, trajFile=None, stride=1, solResName='(!:WAT)', watResName='(:WAT)'):
    self.topFile = topFile
    self.trajFile = trajFile
    self.stride = stride
    self.solResName = solResName
    self.watResName = watResName
    self.top = pmd.load_file(topFile)
    if trajFile is not None:
      self.traj = pt.iterload(trajFile, pt.load_parmed(self.top, traj=False) ,stride=stride)

  def getWatInds(self):
    """Given the above functions and __init__, we can extract the oxygen and hydrogen indice of water molecules
       Inputs:                                                                                                      
               self - contains all necessary objects, in this case just the resnames and loading functions
     
       Outputs:
               watInds - contains the water oxygen indices (numpy array of ints)
               watHInds - contains the water hydrogen indices (numpy array of ints)
               lenWat - # of atoms in a water atom
    """
    #traj = self.loadTraj()
    #top = self.loadTop()

    traj = self.traj
    top = self.top
    
    nWatAtoms = len(traj.top.select(self.watResName))
    watInds = traj.top.select(self.watResName+'&(!@H=)&(!@EP=)')
    watHInds = traj.top.select(self.watResName+'&(@H=)')
    if len(watInds)!=0:
      lenWat = int(nWatAtoms/len(watInds))
    else:
      lenWat = 0
    return watInds, watHInds, lenWat

  def getHeavyInds(self):
    """Given the above functions and __init__, we can extract the heavy atom indices
       Inputs:                  
               self - contains all necessary objects, in this case just the resnames and loading functions
     
       Outputs:
               heavyInds - contains heavy indices (numpy array of ints)
    """
    traj = self.traj
    top = self.top
    # exclude virtual atoms and hydrogens
    heavyInds = traj.top.select('(!@H=)&(!@EP=)')
    return heavyInds

  def getPhobicInds(self):
    """Given the above functions and __init__, we can extract the carbon and sulfur indices of all molecules
       Inputs:                                                                                                      
               self - contains all necessary objects, in this case just the resnames and loading functions
     
       Outputs:
               phobicInds - contains carbon and sulfur indices (numpy array of ints)
    """
    traj = self.traj
    top = self.top
    phobicInds = traj.top.select('(@C=)|(@S=)')
    return phobicInds

  def getPhilicInds(self):
    """Given the above functions and __init__, we can extract the nitrogen and oxygen indices of all molecules
       Inputs:                                                                                                      
               self - contains all necessary objects, in this case just the resnames and loading functions
     
       Outputs:
               philiicInds - contains oxygen and nitrogen indices (numpy array of ints)
    """
    traj = self.traj
    top = self.top
    philicInds = traj.top.select('(@O=)|(@N=)')
    return philicInds

  def getSolInds(self):
    """Given the above functions and __init__, we can extract the heavy atom, hydrogen, carbon, nitrogen, and oxygen indice of non-water cosolvent molecules
       Inputs:                                                                                                      
               self - contains all necessary objects, in this case just the resnames and loading functions
     
       Outputs:
               solInds - contains the cosolvent heavy atom indices (numpy array of ints)
               solHInds - contains the solute hydrogen indices (numpy array of ints)
               solCInds - contains the solute carbon indices (numpy array of ints)
               solNInds - contains the solute nitrogen indices (numpy array of ints)
               solOInds - contains the solute oxygen indices (numpy array of ints)

               solSInds - contains the solute sulfur indices (numpy array of ints)

    """
    traj = self.traj
    top = self.top
    solInds = traj.top.select(self.solResName+'&(!@H=)')
    solHInds = traj.top.select(self.solResName+'&(@H=)')
    solCInds = traj.top.select(self.solResName+'&(@C=)')
    solNInds = traj.top.select(self.solResName+'&(@N=)')
    solOInds = traj.top.select(self.solResName+'&(@O=)')
    solSInds = traj.top.select(self.solResName+'&(@S=)')
    return solInds, solHInds, solCInds, solNInds, solOInds, solSInds

