#! /usr/bin/env python
import sys, os
import argparse
from datetime import datetime
import time as timeit
import pickle
import glob
import copy
import numpy as np
import pytraj as pt
import parmed as pmd

Usage = "This class imports a simulation trajectory and uses parmed and pytraj to extract common atomic sub-indices."

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
    self.topFile = topFile # define topology file string
    self.trajFile = trajFile # define trajectory file string
    self.stride = stride # define trajectory sampling frequency
    self.solResName = solResName # define cpptraj mask for solute atoms
    self.watResName = watResName # define cpptraj mask for water atoms
    self.top = pmd.load_file(topFile) # parmed load topology
    # pytraj load trajectory
    if trajFile is not None:
      self.traj = pt.iterload(trajFile, pt.load_parmed(self.top, traj=False) ,stride=stride)

  def getWatInds(self):
    """This function extracts the oxygen and hydrogen indices of water molecules as defined by the residue masks (per cpptraj/pytraj).
       Outputs:
               watInds - contains the water oxygen indices (numpy array of ints)
               watHInds - contains the water hydrogen indices (numpy array of ints)
               lenWat - # of atoms in a water atom
    """
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
    """This function extracts the heavy indices of the system as defined by the residue masks (per cpptraj/pytraj).
       Outputs:
               heavyInds - contains heavy indices (numpy array of ints)
    """
    traj = self.traj
    top = self.top
    # exclude virtual atoms and hydrogens
    heavyInds = traj.top.select('(!@H=)&(!@EP=)')
    return heavyInds

  def getPhobicInds(self):
    """This function extracts the hydrophobic indices (C and S) of the system as defined by the residue masks (per cpptraj/pytraj).       
       Outputs:
               phobicInds - contains carbon and sulfur indices (numpy array of ints)
    """
    traj = self.traj
    top = self.top
    phobicInds = traj.top.select('(@C=)|(@S=)')
    return phobicInds

  def getPhilicInds(self):
    """This function extracts the hydrophilic indices (O and N) of the system as defined by the residue masks (per cpptraj/pytraj).       
       Outputs:
               philiicInds - contains oxygen and nitrogen indices (numpy array of ints)
    """
    traj = self.traj
    top = self.top
    philicInds = traj.top.select('(@O=)|(@N=)')
    return philicInds

  def getSolInds(self):
    """This function extracts co-solvent (non-water) indices of the system as defined by the residue masks (per cpptraj/pytraj).       
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

