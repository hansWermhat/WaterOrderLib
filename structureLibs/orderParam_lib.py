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
from TrajObject import TrajObject

from matplotlib.image import NonUniformImage

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

Usage = """A library of functions to compute the system-average and local order parameters for water structure.
"""

### Additional functionalities outside of the class
def getHBInds(top, frame, solInds, solHInds, solNInds, solOInds):
  """Here, we utilize parmed to more easily extract the indices of hydrogen bond partners for use with waterlib functions.
     Inputs:                                                                                                     
            top - parmed topology 
            frame - single pytraj trajectory frame
            solInds - list of solute (or water) heavy indices
            solHInds - list of solute (or water) hydrogen indices
            solNInds - list of solute (or water) nitrogen indices
            solOInds - list of solute (or water) oxygen indices
            
     Outputs:
            hbOInds - contains list of lists of acceptorO, donorOInds, and donorHOInds atoms
            hbNInds - contains list of lists of acceptorN, donorNInds, and donorHNInds atoms
    """
  
  acceptorOInds = []
  donorHOInds = []
  donorOInds = []

  acceptorNInds = []
  donorHNInds = []
  donorNInds = []
  
  ### Manually compute and order the H-B donors and aceptors for use in calculation
  # loop thru all atoms in the topology
  for i,atom in enumerate(top.atoms):

    # if oxgygen then proceed 
    if i in solOInds:
      iatoms = atom.bond_partners # get list of i-th oxygen atom bonding partners 
      count = 0 # initialize to zero to track # of hydrogens

      # loop thru all bonded atoms
      for j,jatom in enumerate(iatoms):
        name = jatom.name
        idx = jatom.idx
        if 'H' in name:
          donorHOInds.append(idx)
          count += 1

      acceptorOInds.append(i) # add oxygen atom to acceptor list
      # loop thru all bound H atoms
      for k in range(count):
        donorOInds.append(i) # add oxygen atom for each bound hydrogen atom

    # if Nitrogen then proceed 
    elif i in solNInds:
      iatoms = atom.bond_partners # get list of i-th nitrogen atom bonding partners 
      count = 0 # initialize to zero to track # of hydrogens

      # loop thru all bonded atoms
      for j,jatom in enumerate(iatoms):
        name = jatom.name
        idx = jatom.idx
        if 'H' in name:
          donorHNInds.append(idx)
          count += 1

      acceptorNInds.append(i) # add nitrogen atom to acceptor list
      # loop thru all bound H atoms
      for k in range(count):
        donorNInds.append(i) # add nitrogen atom for each bound hydrogen atom
    else:
      continue

  acceptorOInds = np.array(acceptorOInds, dtype=int)
  acceptorNInds = np.array(acceptorNInds, dtype=int)
  donorOInds = np.array(donorOInds, dtype=int)
  donorNInds = np.array(donorNInds, dtype=int)
  donorHOInds = np.array(donorHOInds, dtype=int)
  donorHNInds = np.array(donorHNInds, dtype=int)

  hbOInds = [acceptorOInds, donorOInds, donorHOInds]
  hbNInds = [acceptorNInds, donorNInds, donorHNInds]
  return hbOInds, hbNInds

### Cluster analysis using sorting library
def getClusters(hbMat):
  """We feed matrix of residue connectivity to determine the graphs.
     Inputs:                      
            hbMat - (nRes, nRes)-numpy array
            
     Outputs:
            clusterSize - contains list of cluster sizes
  """
  clusters = []
  resTrack = np.array([], dtype=int) # numpy array to track residues already accounted for (those participating in a cluster)
  for i in range(hbMat.shape[0]):
    # check to see if residue i is already in a cluster. if so, skip this iteration
    if i in resTrack:
      continue

    iCluster = np.zeros((1, hbMat.shape[0]), dtype=int)
    iCluster = np.zeros(hbMat.shape[0], dtype=int)
    M = np.int64(np.sum(hbMat[i,:])) # need to grab the number of edges directly connecting ith residue to hb partners
    N = hbMat.shape[0] # need to grab the number of residues
    vertex = i+1
    test = sortlib.depthfirstsort(vertex, hbMat, iCluster, M, N)
    iCluster = np.where(test==1)[0]
    resTrack = np.concatenate((resTrack, iCluster))

    if len(iCluster)==hbMat.shape[0]: # if all molecules are in single cluster, no need to continue
      clusters.append(iCluster)
      break
    elif len(iCluster)==0: # define 0-length cluster as equal to 1 to account for # of unbound molecules
      clusters.append(np.array([i]))
      continue
    else: # in all other cases, just add the set of molecular indices contributing to the cluster
      clusters.append(iCluster)
      
  return clusters

def getHBClusterStats(topFile, trajFile, acceptorInds, donorInds, donorHInds,
                      stride=1, distCut = 3.0, angCut = 150.0):
  """Here, we extract a matrix of residue-residue HBs. This should work fairly generally,
  e.g., watInds=boundWatInds and solInds=(nonBoundWatInds & solInds) should produce an
  hbMat with dimension (nRes_bound, nRes-nRes_bound)
     Inputs:                                                               
            topFile - topology name 
            trajFile - trajectory name
            acceptorInds - list of acceptor heavy atom indices
            donorInds - list of donor heavy atom indices
            donorHInds - list of donor hydrogen indices            
            stride - skip every "stride" steps of the trajectory (int). 
                     Default=1

     Outputs:
            hbMat - (nRes,nRes)-array of HB contacts between residues
    """
  obj = TrajObject(topFile, trajFile=trajFile, stride=stride, solResName=None,
                    watResName=None)
  top = obj.top #loadTop()     
  traj = obj.traj

  clusters = []
  for t, frame in enumerate(traj):
    # get box vectors and positions
    thisbox = np.reshape(np.array(frame.box.values[:3]),(1,3))
    thispos = np.array(frame.xyz)

    # get positions
    acceptorPos = thispos[acceptorInds]
    donorPos = thispos[donorInds]
    donorHPos = thispos[donorHInds]

    start = timeit.time()
    # get all hbs
    allHB = wl.generalhbonds(acceptorPos,
                             donorPos, donorHPos,
                             thisbox, distCut, angCut)
    
    # create mapping from acceptor Ind to atom index
    acceptMap = [acceptorInds[i] for i in range(len(acceptorInds))]
    donorHMap = [donorHInds[i] for i in range(len(donorHInds))]

    resAccept = np.zeros(allHB.shape[0], dtype=int)
    resDonorH = np.zeros(allHB.shape[1], dtype=int)

    # get acceptor heavy atom and donor hydrogen residue indices
    for i in range(resAccept.shape[0]):
      resAccept[i] = top.atoms[acceptMap[i]].residue.idx

    for i in range(resDonorH.shape[0]):
      resDonorH[i] = top.atoms[donorHMap[i]].residue.idx
  
    ### Manually compute and order the H-B donors and aceptors for use in calculation
    hbMat = np.zeros((len(top.residues), len(top.residues)))
    # loop thru all atoms in the topology
    for i in range(len(top.residues)):
      # convert residue acceptor and donor inds from atom basis to allHB basis
      res_acceptConv = np.where(resAccept==i)[0]
      res_donHConv = np.where(resDonorH==i)[0]

      acceptHB = allHB[res_acceptConv, :]
      donorHB = allHB[:, res_donHConv]

      donInds = np.unique(np.where(acceptHB==1)[1])
      acceptInds = np.unique(np.where(donorHB==1)[0])
      resPairs = np.concatenate((resAccept[acceptInds], resDonorH[donInds]))

      # update hbMat
      hbMat[i, resPairs] = 1
    
    iClusters = getClusters(hbMat)
    clusterSize = np.array([len(clust) for clust in iClusters 
                            if len(clust)!=1])
    
    clusters.append(clusterSize)
  
  clusters = np.concatenate(clusters)
  meanCluster = np.mean(clusters)
  return meanCluster

def getIonClusterStats(topFile, trajFile, Inds, chargeAssign,
                        stride=1, distCut = 3.4):
  """Here, we extract a matrix of contacts for non-HB cluster evaluation.
     Inputs:                                                               
            topFile - topology name 
            trajFile - trajectory name
            Inds - list of atom indices            
            chargeAssign - len(Inds) long list of ion charges
            stride - skip every "stride" steps of the trajectory (int). 
                     Default=1

     Outputs:
            hbMat - (nRes,nRes)-array of HB contacts between residues
    """
  obj = TrajObject(topFile, trajFile=trajFile, stride=stride, solResName=None,
                    watResName=None)
  top = obj.top #loadTop()     
  traj = obj.traj

  clusters = []
  cations = []
  
  cationInds = [i for i in range(len(Inds)) if chargeAssign[i]>0]
  anionInds = [i for i in range(len(Inds)) if chargeAssign[i]<0]

  for t, frame in enumerate(traj):
    # get box vectors and positions
    thisbox = np.reshape(np.array(frame.box.values[:3]),(1,3))
    thispos = np.array(frame.xyz)

    # get positions
    subPos = thispos[Inds]
    
    start = timeit.time()
    # get all pairs
    pairMat = wl.allnearneighbors(subPos, thisbox, 0.0, distCut)
          
    tClusters = getClusters(pairMat)
    
    # get charges of these clusters
    tCharges = [ chargeAssign[cluster] for cluster in tClusters ]

    # compute cluster size and charge
    clusterSize = np.array( [ len(clust) for clust in tClusters ] )
    clusterCharge = np.array( [ np.sum(charge) for charge in tCharges ] )  

    # loop through clusters and assign effective charges for each cation based on which cluster they exist in
    cationCharge = []
    for i, zEff in enumerate(clusterCharge):
      cluster = tClusters[i]
      if any(x in cluster for x in cationInds):
        cationCharge.append(zEff)
      else:
        continue

    cations.append(np.array(cationCharge))
    clusters.append(clusterSize)
  
  clusters = np.concatenate(clusters)
  cations = np.concatenate(cations)
  meanCluster = np.mean(clusters)
  meanCharge = np.mean(cations)
  
  clusterDist, bins = np.histogram(clusters,
                                   bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                   density=False)

  # save the file!                  
  np.savetxt('clusterDistribution.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
                                                  clusterDist],axis=1),
             header='# clusters    frequency',fmt="%.3e")

  return meanCluster

def getNeighborStats(topFile, trajFile, Inds1, Inds2, nAtoms1, nAtoms2, 
                     stride=1, distCut = 3.4, switch=False):
  """Here, we extract a matrix of contacts for non-HB cluster evaluation.
     Inputs:                                                               
            topFile - topology name 
            trajFile - trajectory name
            Inds1 - list of atom indices 1 
            Inds2 - list of atom indices 2
            nAtoms1 - number of atoms in each molecule of type 1
            nAtoms2 - number of atoms in each molecule of type 2
            switch - logical expression that's True if Inds1 and Inds2 are the same
            stride - skip every "stride" steps of the trajectory (int). 
                     Default=1

     Outputs:
            meanNeighbors - average number of neighbors near the minority atom set
    """
  obj = TrajObject(topFile, trajFile=trajFile, stride=stride, solResName=None,
                    watResName=None)
  top = obj.top #loadTop()     
  traj = obj.traj

  numberCoord = [] # list to store the coordination environment in each timestep
  for t, frame in enumerate(traj):
    # get box vectors and positions
    thisbox = np.reshape(np.array(frame.box.values[:3]),(1,3))
    thispos = np.array(frame.xyz)

    # get positions
    subPos1 = thispos[Inds1]
    subPos2 = thispos[Inds2]

    if switch:
      neighbors = wl.allnearneighbors(subPos1, thisbox, 0.0, distCut)
      nRes = int(len(Inds1)/nAtoms1)      
      resNumbers = np.zeros(nRes, dtype=int)
      
      for n in range(nRes):
        nNeighbors = neighbors[int(n*nAtoms1):int((n+1)*nAtoms1), :]
        nNeighbors[int(n*nAtoms1):int((n+1)*nAtoms1), 
                   int(n*nAtoms1):int((n+1)*nAtoms1)] = 0
        resNumbers[n] = len(np.unique(np.where(nNeighbors==1)[1]))

      numberCoord.append(resNumbers)

    else:
      neighbors = wl.nearneighbors(subPos1, subPos2, thisbox, 0.0, distCut)
      nRes1 = int(len(Inds1)/nAtoms1)
      nRes2 = int(len(Inds2)/nAtoms2)
      resNumbers = np.zeros(nRes1, dtype=int)
      
      for n in range(nRes1):
        nNeighbors = neighbors[int(n*nAtoms1):int((n+1)*nAtoms1), :]
        resNumbers[n] = len(np.unique(np.where(nNeighbors==1)[1]))

      numberCoord.append(resNumbers)          

      
  # concatenate numberCoord
  numberCoord = np.concatenate(numberCoord)
  meanCoord = np.mean(numberCoord)
  
  coordDist, bins = np.histogram(numberCoord,
                                 bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                                 density=False)

  # save the file!                  
  np.savetxt('coordDistribution.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
                                                coordDist],axis=1),
             header='# coords    frequency',fmt="%.3e")

  return meanCoord

### computing CI from block averaging samples                                  
def getCI(means):
  meanCI = means[int(0.5*len(means))]
  upperCI = means[int(0.975*len(means))] - meanCI
  lowerCI = meanCI - means[int(0.025*len(means))]
  return max(upperCI, lowerCI)

### block averaging functionality, takes in a list of properties and outputs error bars  
def blockAverage(vals, nBlocks=20):
  # compute the number of observations                                        
  nObs = len(vals)

  # set number of blocks to split traj into                                    
  obsBlocks = np.zeros(nBlocks)
  lenBlock = len(vals)/nBlocks

  for i in range(nBlocks):
    obsBlocks[i] = np.mean(vals[int(i*lenBlock):int((i+1)*lenBlock)])

  nSamp = nBlocks
  nResamp = 10000
  obsMeans = np.zeros(nResamp)

  # randomly select nResamp blocks to reconstruct the nSamp-long dataset      
  for n in range(nResamp):
    obsMeans[n] = np.mean( np.random.choice(obsBlocks, nSamp) )

  obsMeans = np.sort(obsMeans)

  # determine confidence intervals                                             
  CI = getCI(obsMeans)
  return CI

def getBoundWrap(topFile, frame, watInds, watHInds, 
                 solInds, solHInds, solCInds, solOInds, solNInds, solSInds,
                 cutoff=4.0, hbDist=3.0, hbAng=150.0):
  """Given a topology file and a trajectory frame, computes the bound, wrap, and non-shell water indices 
     Inputs:                                                                                                        
             topFile - topology file (path from cwd)                                              
             trajFile - trajectory file name needed for object (TrajObject class... maybe)
                  
             frame - frame fed into nc                                           
  
             watInds - water oxygen indices

             watHInds - water hydrogen indices

             solInds - all cosolvent (solute) heavy atom indices

             solHInds - cosolvent (solute) hydrogen atom indices

             solCInds - cosolvent (solute) carbon atom indices

             solOInds - cosolvent (solute) oxygen atom indices

             solNInds - cosolvent (solute) nitrogen atom indices

             solSInds - cosolvent (solute) sulfur atom indices


             cutoff - (Default 4-A) the hydration layer cutoff, in angstroms

             hbDist - (Default 3-A) the hb-disance cutoff, in angstroms

             hbAng - (Default 150-deg) the hb-angle cutoff, in degrees

     Outputs:                         
             boundInds - list of bound water indices                             
             wrapInds - list of wrap water indices
             shellInds - list of shell water indices (bound and wrap)
             nonShellInds - list of "non-shell" indices (neither bound nor wrap)
  """
  # load topologies and trajectories using the above class AnalysisObject        
  obj = TrajObject(topFile, trajFile=None, stride=1, solResName=None, 
                    watResName=None)

  top = obj.top #loadTop()

  # get hydrogen bonding atoms
  hbOInds, hbNInds = getHBInds(top, frame, solInds, solHInds, solNInds, solOInds)

  # delineate solute oxygen acceptor and donor indices
  sol_acceptorOInds = hbOInds[0]
  sol_donorOInds = hbOInds[1]
  sol_donorHOInds = hbOInds[2]

  # delineate solute nitrogen acceptor and donor indices
  sol_acceptorNInds = hbNInds[0]
  sol_donorNInds = hbNInds[1]
  sol_donorHNInds = hbNInds[2]

  pos = np.array(frame.xyz)
  thisbox = np.array(frame.box.values[:3])

  # extract water and solute positions
  watPos = pos[watInds]
  watHPos = pos[watHInds]
  solPos = pos[solInds]

  # extract acceptro and donor positions
  sol_acceptorOPos = pos[sol_acceptorOInds]
  sol_donorOPos = pos[sol_donorOInds]
  sol_donorHOPos = pos[sol_donorHOInds]
  
  sol_acceptorOPos = pos[sol_acceptorOInds]
  sol_donorOPos = pos[sol_donorOInds]
  sol_donorHOPos = pos[sol_donorHOInds]
      
  # compute water neighboring solutes
  neighbors = wl.nearneighbors(solPos, watPos, thisbox, 0.0, cutoff)
  mask = np.unique(np.where( neighbors==1 )[1])
  shellInds = watInds[mask]
  nonShellInds = np.delete(watInds, mask)
      
  ### Currently disabled: Search for waters unique to a single cosolvent molecule
  if False:
    shellPos = []
    shellInds = []
    bannedInds = []
    for i in range(len(solInds)/6):
      isolPos = solPos[i*6:(i+1)*6-1]
      neighbors = wl.nearneighbors(isolPos, watPos, thisbox, 0.0, cutoff)
      ishellInds = watInds[ np.unique( np.where( neighbors==1 )[1] ) ]
      for ind in ishellInds:
        if ind not in bannedInds:
          if ind not in shellInds:
            shellInds.append(ind)
          elif ind in shellInds:
            shellInds.remove(ind)
            bannedInds.append(ind)
          else:
            continue

    shellInds = np.array(shellInds)
  ### Currently Disabled

  # extract the acceptor and donor indices for the shell waters
  hbOInds, _ = getHBInds(top, frame, shellInds, watHInds, solNInds, shellInds)

  wat_acceptorOInds = hbOInds[0]
  wat_donorOInds = hbOInds[1]
  wat_donorHOInds = hbOInds[2]

  wat_acceptorOPos = pos[wat_acceptorOInds]
  wat_donorOPos = pos[wat_donorOInds]
  wat_donorHOPos = pos[wat_donorHOInds]

  # compute wat-wat HBs
  watWatHBs = wl.generalhbonds(wat_acceptorOPos, 
                               wat_donorOPos, wat_donorHOPos,
                               thisbox, hbDist, hbAng)

  # compute glyc-glyc HBs
  solSolHBs = wl.generalhbonds(sol_acceptorOPos,
                               sol_donorOPos, sol_donorHOPos,
                               thisbox, hbDist, hbAng)

  # compute HBs where water is the acceptor
  watSolHBs = wl.generalhbonds(wat_acceptorOPos,
                               sol_donorOPos, sol_donorHOPos,
                               thisbox, hbDist, hbAng)

  boundMask_wat = np.unique(np.where(watSolHBs==1)[0])
  
  # compute HBs where the cosolvent is the acceptor
  solWatHBs = wl.generalhbonds(sol_acceptorOPos,
                               wat_donorOPos, wat_donorHOPos,
                               thisbox, hbDist, hbAng)

  # set up dummy variable to quickly compute unique sol-wat HBs
  dummy = np.zeros(len(wat_donorOPos))
  dummy[np.unique(np.where(solWatHBs==1)[1])] = 1
      
  # compute shell water indices corresponding to bound wats
  boundMask_sol = np.where(np.ceil(0.5*(dummy[0::2]+dummy[1::2])))[0]

  # combine bound indices
  boundMask = np.sort(np.unique(np.concatenate([boundMask_wat, 
                                                boundMask_sol]))) 

  # grab wrap inds using mask
  mask = np.ones(len(shellInds), dtype=bool)
  mask[boundMask] = False
  wrapInds = shellInds[mask]
  boundInds = shellInds[boundMask]
            
  return boundInds, wrapInds, shellInds, nonShellInds

### For now this will do. It gives us a lot of useful information that we can used not only in further analyses, but for other functions in this library! Maybe I can add a functionality that gives error bars on the rdf from block averaging, but want to be a bit more confident in my implementation first...
def rdfCalc(topFile, trajFile, solResName='(!:WAT)', watResName='(:WAT)', binwidth=0.1,totbins=150, stride=1):
  """Given a topology file and trajectory (or a list of) file, computes the Ow-Ow, solute-Ow, and solute-solute 2D radial distribution funtions. Though our primary interests lie in water-species coordination, we may arbitrarily state the residue names according to pytraj/parmed masking schemes to get the coordination between species #1 (solResName) and species #2 (watResName). Here, solute atoms will be the set of heavy atoms.
     Inputs:
             topFile - topology file (path from cwd)
             trajFile - trajectory file, or list of them (path from cwd)
             solResName - string defining the residue name for the non-water cosolvent. Default='(!:WAT)'
             watResName - string defining the residue name for the water. Default='(:WAT)'
             binwidth - width of RDF histogram bins in angstroms (float). Default=0.1
             totbins - total # of RDF histogram bins (int). Default=200
             stride - skip every "stride" steps of the trajectory (int). Default=1
     Outputs:
             cutoff - estimate of the location of the 2nd sol-Ow min. (Angstroms)

  """
  # load topologies and trajectories using the above class AnalysisObject
  obj = TrajObject(topFile, trajFile, stride, solResName, watResName)
  top = obj.top #loadTop()
  traj = obj.traj #.loadTraj()

  # select water oxygen and hydrogen indices
  watInds, watHInds, lenWat = obj.getWatInds()

  # select cosolvent heavy
  solInds, solHInds, solCInds, solNInds, solOInds, solSInds = obj.getSolInds()

  # set up rdf and coord storage lists
  tot_rdf_OwOw = []; tot_rdf_SolOw = []; tot_rdf_SolSol = []
  tot_coord_OwOw = []; tot_coord_SolOw = []; tot_coord_SolSol = []
  tot_n1_OwOw = []; tot_n1_SolOw = []; tot_tParam = []

  # set up number of trajectory chunks to sample
  nChunks = 5
  chunkSize = int(len(traj)/nChunks)

  for c in range(nChunks):

    # define the rdf vectors to store data in with zeros. For now, only considering these three rdfs for simplicity, probably will only need OwOw and SolOw to probe wat. structure and find cutoff distance, respectively.
    rdf_OwOw = np.zeros(totbins)
    rdf_SolOw = np.zeros(totbins)
    rdf_SolSol = np.zeros(totbins)

    # loop through trajectory to compute average RDFs
    for t, frame in enumerate(traj[int(c*chunkSize):int((c+1)*chunkSize)]):
      # grab ith positions and box vectors
      thispos = np.array(frame.xyz)
      thisbox = np.reshape(np.array(frame.box.values[:3]),(1,3))
    
      thisWat = np.array(thispos[watInds])
      thisSol = np.array(thispos[solInds])
       
      dist = np.linspace(0,(totbins-1)*binwidth,totbins) + binwidth

      bulkdens=1.0 # set bulk density to 1.0 so we have local density plots
      rdf_OwOw = rdf_OwOw + wl.radialdistsame(thisWat,binwidth,totbins,bulkdens,thisbox)
    
      # include conditional to prevent calculation of sol rdfs for water
      if solInds==[]:
        continue
      else:
        rdf_SolSol = rdf_SolSol + wl.radialdistsame(thisSol,binwidth,totbins,bulkdens,thisbox)
        rdf_SolOw = rdf_SolOw + wl.radialdist(thisSol,thisWat,binwidth,totbins,bulkdens,thisbox)
  
    # calculate the average rdf
    rdf_OwOw = rdf_OwOw/(1.0+t)
    rdf_SolSol = rdf_SolSol/(1.0+t)
    rdf_SolOw = rdf_SolOw/(1.0+t)

    # append rdfs to sample array
    tot_rdf_OwOw.append(rdf_OwOw)
    tot_rdf_SolSol.append(rdf_SolSol)
    tot_rdf_SolOw.append(rdf_SolOw)

    # compute coordination numbers
    coord_OwOw = np.zeros(len(dist)-2)
    coord_SolOw = np.zeros(len(dist)-2)
    coord_SolSol = np.zeros(len(dist)-2)
    for j in range(2,len(dist)):
      coord_OwOw[j-2] = 8.0*np.pi*simps(rdf_OwOw[:j]*(dist[:j])**2.0, 
                                        dist[:j])
      if solInds!=[]:
        coord_SolOw[j-2] = 4.0*np.pi*simps(rdf_SolOw[:j]*(dist[:j])**2.0, 
                                           dist[:j])
        coord_SolSol[j-2] = 8.0*np.pi*simps(rdf_SolSol[:j]*(dist[:j])**2.0, 
                                            dist[:j])
  
    tot_coord_OwOw.append(coord_OwOw)
    tot_coord_SolSol.append(coord_SolSol)
    tot_coord_SolOw.append(coord_SolOw)

    # once again only calculate this if solInds is not empty
    if solInds!=[]:
      # estimate the second minimum in the Sol-Ow rdf
      relMin_SolOW = rdf_SolOw[argrelmin(rdf_SolOw)]
      cutoff = dist[argrelmin(rdf_SolOw)][:3]

      # estimate n1_SolOw @ "first" min.
      n1_SolOw = coord_SolOw[argrelmin(rdf_SolOw)[0][0]-2]
      tot_n1_SolOw.append(n1_SolOw)

    relMin_OwOw = rdf_OwOw[argrelmin(rdf_OwOw)]
    n1_OwOw = coord_OwOw[argrelmin(rdf_OwOw)[0][0]-2]

    rdf = rdf_OwOw[:argrelmin(rdf_OwOw)[0][0]]/rdf_OwOw[-1]
    rdf_dist = dist[:argrelmin(rdf_OwOw)[0][0]]
    rc = dist[argrelmin(rdf_OwOw)[0][0]]
          
    tParam = simps(rdf, rdf_dist)/rc
    tot_n1_OwOw.append(n1_OwOw)
    tot_tParam.append(tParam)

  # convert rdf and coord nums to arrays
  tot_rdf_OwOw = np.array(tot_rdf_OwOw)
  tot_rdf_SolOw = np.array(tot_rdf_SolOw)
  tot_rdf_SolSol = np.array(tot_rdf_SolSol)

  tot_coord_OwOw = np.array(tot_coord_OwOw)
  tot_coord_SolOw = np.array(tot_coord_SolOw)
  tot_coord_SolSol = np.array(tot_coord_SolSol)

  # compute mean and SE of mean for each quant
  rdf_OwOw_se = np.std(tot_rdf_OwOw, axis=0, ddof=1)/np.sqrt(nChunks-1)
  rdf_SolOw_se = np.std(tot_rdf_SolOw, axis=0, ddof=1)/np.sqrt(nChunks-1)
  rdf_SolSol_se = np.std(tot_rdf_SolSol, axis=0, ddof=1)/np.sqrt(nChunks-1)

  coord_OwOw_se = np.std(tot_coord_OwOw, axis=0, ddof=1)/np.sqrt(nChunks-1)
  coord_SolOw_se = np.std(tot_coord_SolOw, axis=0, ddof=1)/np.sqrt(nChunks-1)
  coord_SolSol_se = np.std(tot_coord_SolSol, axis=0, ddof=1)/np.sqrt(nChunks-1)

  n1_OwOw_se = np.std(np.array(tot_n1_OwOw), ddof=1)/np.sqrt(nChunks-1)
  n1_SolOw_se = np.std(np.array(tot_n1_SolOw), ddof=1)/np.sqrt(nChunks-1)
  tParam_se = np.std(np.array(tot_tParam), ddof=1)/np.sqrt(nChunks-1)

  n1_OwOw = np.mean(tot_n1_OwOw)
  n1_SolOw = np.mean(tot_n1_SolOw)
  tParam = np.mean(tot_tParam)

  np.savetxt('rdf.txt', np.stack([dist,
                                  rdf_OwOw, rdf_OwOw_se,
                                  rdf_SolSol, rdf_SolSol_se,
                                  rdf_SolOw, rdf_SolOw_se], 
                                 axis=1),
             header='pair distance (A)     Ow-Ow rdf     err     Sol-Sol rdf     err     Sol-Ow rdf     err',fmt="%.3e")
  np.savetxt('coord.txt', np.stack([dist[2:], 
                                    coord_OwOw, coord_OwOw_se,
                                    coord_SolSol, coord_SolSol_se,
                                    coord_SolOw, coord_SolOw_se],
                                   axis=1),
             header='pair distance (A)     Ow-Ow n1     err     Sol-Sol n1     err     Sol-Ow n1     err',fmt="%.3e")

  if solInds!=[]:
    return [n1_OwOw, n1_OwOw_se], [n1_SolOw, n1_SolOw_se], [tParam, tParam_se]
  else:
    return n1_OwOw, t

def hbCalc(topFile, trajFile, solResName='(!:WAT)', watResName='(:WAT)', stride=1):
  """Given a topology file and trajectory (or a list of) file, computes the average number of water and cosolvent hydrogen bonds per molecule
     Inputs:                                                                                                        
             topFile - topology file (path from cwd)                                                                
             trajFile - trajectory file, or list of them (path from cwd)                                            
             solResName - string defining the residue name for the non-water cosolvent. Default='(!:WAT)'               
             watResName - string defining the residue name for the water. Default='(:WAT)'                             
             stride - skip every "stride" steps of the trajectory (int). Default=1                                  
     Outputs:                     
             avgHBs - avg # species HBs                              
  """
  # load topologies and trajectories using the above class AnalysisObject  
  obj = TrajObject(topFile, trajFile, stride, solResName, watResName)
  top = obj.top #loadTop()
  traj = obj.traj #loadTraj()

  # select water oxygen and hydrogen indices                               
  watInds, watHInds, lenWat = obj.getWatInds()

  # select cosolvent heavy atoms
  solInds, solHInds, solCInds, solNInds, solOInds, solSInds = obj.getSolInds()
  hbOInds, hbNInds = getHBInds(top, traj[0], solInds, solHInds, solNInds, solOInds)

  # delineate solute oxygen acceptor and donor indices
  sol_acceptorOInds = hbOInds[0]
  sol_donorOInds = hbOInds[1]
  sol_donorHOInds = hbOInds[2]

  # delineate solute nitrogen acceptor and donor indices
  sol_acceptorNInds = hbNInds[0]
  sol_donorNInds = hbNInds[1]
  sol_donorHNInds = hbNInds[2]

  # extract the acceptor and donor indices for the shell waters
  hbOInds, _ = getHBInds(top, traj[0], watInds, watHInds, [], watInds)
  wat_acceptorOInds = hbOInds[0]
  wat_donorOInds = hbOInds[1]
  wat_donorHOInds = hbOInds[2]

  # compute number of per cosolute accepting heavy atoms and donor hydrogens
  nSol = traj[:1, solResName].topology.n_residues

  nAccO = int(len(sol_acceptorOInds)/nSol)
  nAccN = int(len(sol_acceptorNInds)/nSol)
  nDonO = int(len(sol_donorOInds)/nSol)
  nDonN = int(len(sol_donorNInds)/nSol)
  nDonHO = int(len(sol_donorHOInds)/nSol)
  nDonHN = int(len(sol_donorHNInds)/nSol)

  # set up empty list to hold number of wat-wat hbs
  numWatHBs = []
  numSolHBs = []

  for t,frame in enumerate(traj):
    pos = np.array(frame.xyz)
    thisbox = np.array(frame.box.values[:3])

    # extract water and solute positions
    watPos = pos[watInds]
    watHPos = pos[watHInds]
    solPos = pos[solInds]

    # extract acceptor and donor positions
    sol_acceptorOPos = pos[sol_acceptorOInds]
    sol_acceptorNPos = pos[sol_acceptorNInds]

    sol_donorOPos = pos[sol_donorOInds]
    sol_donorNPos = pos[sol_donorNInds]

    sol_donorHOPos = pos[sol_donorHOInds]
    sol_donorHNPos = pos[sol_donorHNInds]

    wat_acceptorOPos = pos[wat_acceptorOInds]
    wat_donorOPos = pos[wat_donorOInds]
    wat_donorHOPos = pos[wat_donorHOInds]

    watwatHBs = wl.generalhbonds(wat_acceptorOPos,
                                 wat_donorOPos,wat_donorHOPos,
                                 thisbox, 3.5, 120.0)

    watsolOHBs = wl.generalhbonds(wat_acceptorOPos,
                                  sol_donorOPos,sol_donorHOPos,
                                  thisbox, 3.5, 120.0)
    solwatOHBs = wl.generalhbonds(sol_acceptorOPos,
                                  wat_donorOPos,wat_donorHOPos,
                                  thisbox, 3.5, 120.0)

    watsolNHBs = wl.generalhbonds(wat_acceptorOPos,
                                  sol_donorNPos,sol_donorHNPos,
                                  thisbox, 3.5, 120.0)
    solwatNHBs = wl.generalhbonds(sol_acceptorNPos,
                                  wat_donorOPos,wat_donorHOPos,
                                  thisbox, 3.5, 120.0)

    solOsolOHBs = wl.generalhbonds(sol_acceptorOPos,
                                   sol_donorOPos,sol_donorHOPos,
                                   thisbox, 3.5, 120.0)
    solOsolNHBs = wl.generalhbonds(sol_acceptorOPos,
                                   sol_donorNPos,sol_donorHNPos,
                                   thisbox, 3.5, 120.0)
    solNsolOHBs = wl.generalhbonds(sol_acceptorNPos,
                                   sol_donorOPos,sol_donorHOPos,
                                   thisbox, 3.5, 120.0)
    solNsolNHBs = wl.generalhbonds(sol_acceptorNPos,
                                   sol_donorNPos,sol_donorHNPos,
                                   thisbox, 3.5, 120.0)


#    if watsolOHBs.shape[1]==0:
#      watsolOHBs = np.zeros((len(watInds),1))
#    if watsolNHBs.shape[1]==0:
#      watsolNHBs = np.zeros((len(watInds),1))

    ### HARD CODED just to get the job done here
    # compute per cosolute molecule water-cosolute bonds
    solOAcc = np.sum(solwatOHBs, axis=1) + np.sum(solOsolOHBs, axis=1)
    solOAcc += np.sum(solOsolNHBs, axis=1)

    solODon = np.sum(watsolOHBs, axis=0) + np.sum(solOsolOHBs, axis=0)
    solODon += np.sum(solNsolOHBs, axis=0)

    solOAcc = sum( [ solOAcc[i::nAccO] for i in range(nAccO) ] )
    solODon = sum( [ solODon[i::nDonO] for i in range(nDonO) ] )

    solNAcc = np.sum(solwatNHBs, axis=1) + np.sum(solNsolNHBs, axis=1)
    solNAcc += np.sum(solNsolOHBs, axis=1)

    solNDon = np.sum(watsolNHBs, axis=0) + np.sum(solNsolNHBs, axis=0)
    solNDon += np.sum(solOsolNHBs, axis=0)

    solNAcc = sum( [ solNAcc[i::nAccN] for i in range(nAccN) ] )
    solNDon = sum( [ solNDon[i::nDonN] for i in range(nDonN) ] )

    solTot = solNAcc + solNDon + solOAcc + solODon

    numSolHBs.append(solTot)

    # compute per water molecule hydrogen bonds
    watwatAcc = np.sum(watwatHBs, axis=1)
    watwatDon = np.sum(watwatHBs, axis=0)
    watwatDon = watwatDon[::2] + watwatDon[1::2]

    watsolOAcc = np.sum(watsolOHBs, axis=1)
    solwatODon = np.sum(solwatOHBs, axis=0)
    solwatODon = solwatODon[::2] + solwatODon[1::2]

    watsolNAcc = np.sum(watsolNHBs, axis=1)
    solwatNDon = np.sum(solwatNHBs, axis=0)
    solwatNDon = solwatNDon[::2] + solwatNDon[1::2]

    watwatTot = watwatAcc + watwatDon
    watsolOTot = watsolOAcc + solwatODon
    watsolNTot = watsolNAcc + solwatNDon
        
    numTot = watwatTot + watsolOTot + watsolNTot
    numWatHBs.append(numTot) 

  # check if we need to concatenate
  if type(numWatHBs[0]) is not np.ndarray:
    pass
  else:
    numWatHBs = np.concatenate(numWatHBs)
  
  if type(numSolHBs[0]) is not np.ndarray:
    pass
  else:
    numSolHBs = np.concatenate(numSolHBs)
  
  avgWatHBs = np.mean(numWatHBs)
  avgSolHBs = np.mean(numSolHBs)

  hbDist, bins = np.histogram(numWatHBs, 
                              bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                              density=False)

  # save the file!
  np.savetxt('hbDistribution_water.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
                                                   hbDist],axis=1),
             header='# hbs    frequency',fmt="%.3e")

  hbDist, bins = np.histogram(numSolHBs, 
                              bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                              density=False)
    
  # save the file!
  np.savetxt('hbDistribution_cosolv.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
                                                    hbDist],axis=1),
             header='# hbs    frequency',fmt="%.3e")
  return avgWatHBs, avgSolHBs

# use this to calculate voronoi volumes in voronoiCalc
def voronoi_volumes(points, boxL, numWats):
  # define arrays to store voronoi vol and area
  vol = np.zeros(len(points))
  area = np.zeros(len(points))

  # add points across periodic boundaries of the simulation box
  new_points = []
  new_points.append(points)
  for i, point in enumerate(points):
    if point[0]<0.5*boxL:
      new_points.append(np.array([0.-points[i,0]+0, 
                                  points[i, 1], points[i, 2]]))
    if point[0]>0.5*boxL:
      new_points.append(np.array([boxL-points[i,0]+boxL, 
                                  points[i, 1], points[i, 2]]))
    
    if point[1]<0.5*boxL:
      new_points.append(np.array([points[i, 0], 0.-points[i,1]+0, 
                                  points[i, 2]]))
    if point[1]>0.5*boxL:
      new_points.append(np.array([points[i, 0], boxL-points[i,1]+boxL, 
                                  points[i, 2]]))

    if point[2]<0.5*boxL:
      new_points.append(np.array([points[i, 0], points[i, 1],
                                  0.-points[i,2]+0]))
    if point[2]>0.5*boxL:
      new_points.append(np.array([points[i, 0], points[i, 1],
                                  boxL-points[i,2]+boxL]))

  points = np.vstack(new_points)
  v = Voronoi(points)
  
  for i, reg_num in enumerate(v.point_region[:numWats]):
    indices = v.regions[reg_num]
    if -1 in indices: # some regions can be opened
      vol[i] = np.inf
      area[i] = np.inf
    else:
      vol[i] = ConvexHull(v.vertices[indices], qhull_options='QJ').volume
      area[i] = ConvexHull(v.vertices[indices], qhull_options='QJ').area

  return vol, area

def voronoiCalc(topFile, trajFile, subInds=None, nPops=0, solResName='(!:WAT)', watResName='(:WAT)', stride=1, cutoff=4.0, hbDist=3.0, hbAng=150.0):
  """Given a topology file and trajectory (or a list of) file, computes the Ow-Ow, solute-Ow, and solute-solute 2D radial distribution funtions. Here, solute atoms will be the set of heavy atoms.
     Inputs:
             topFile - topology file (path from cwd)
             trajFile - trajectory file, or list of them (path from cwd)
             solResName - string defining the residue name for the non-water cosolvent. Default='(!:WAT)'
             watResName - string defining the residue name for the water. Default='(:WAT)'
             subInds - Options:
                       None - only get bulk quantity (Default)

                       'bound' - perform getBoundWrap() to extract 
                                 bound and wrap indices for:
                                 cutoff, hbDist and hbAng defined in voronoiCalc() 

                       list containing water indices for n-populations at each 
                       timestep 
                       (e.g., [[boundInds,wrapInds]_{0},...,[boundInds,wrapInds]_{t}])

             nPops - the number of populations in subInds (Default=1)

             stride - skip every "stride" steps of the trajectory (int). Default=1

             cutoff - solute hydration shell cutoff radius in Angstroms (Default=4.0)
             hbDist - hydrogen bond cutoff distance in Angstroms (Default=3.0)
             hbAng - hydrogen bond cutoff angle in degrees (Default=150.0)

     Outputs:
             avgLSI - estimate of the mean LSI (Angstrom^2)
             stdLSI - estimate of the standard dev. LSI (Angstrom^2)

  """
  # load topologies and trajectories using the above class AnalysisObject
  obj = TrajObject(topFile, trajFile, stride, solResName, watResName)
  top = obj.top #loadTop()
  traj = obj.traj #.loadTraj()

  # select water oxygen and hydrogen indices
  watInds, watHInds, lenWat = obj.getWatInds()

  # select cosolvent heavy
  solInds, solHInds, solCInds, solNInds, solOInds, solSInds = obj.getSolInds()

  # get all heavy indices
  heavyInds = np.concatenate((watInds,solInds))

  # create mapping function from water index to heavyInds equivalent
  mapHeavy = {watInds[i]:i for i in range(len(watInds))}

  # Create lists to store volume and area values.
  watVol = [ [] for i in range((nPops+1)) ]
  watArea = [ [] for i in range((nPops+1)) ]
  watEta = [ [] for i in range((nPops+1)) ]

  # Lists to store mean and var. volume, area, and asphericity
  avgArea = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  varArea = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  avgVol = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  varVol = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  avgEta = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  varEta = [ np.zeros(len(traj)) for i in range((nPops+1)) ]

  # loop through trajectory to compute average RDFs
  for t, frame in enumerate(traj):

    if subInds is None:
      inds = [ [ mapHeavy[ ind ] for ind in watInds ] ]
    elif subInds=='bound':
      boundInds, wrapInds, shellInds, nonShellInds = getBoundWrap(topFile, frame, 
                                                                  watInds, watHInds,
                                                                  solInds, solHInds, 
                                                                  solCInds, solOInds, 
                                                                  solNInds, solSInds,
                                                                  cutoff=cutoff, 
                                                                  hbDist=hbDist, 
                                                                  hbAng=hbAng)
      # use boundWrap to extract bound and wrap indices
      inds = [ mapHeavy[boundInds], mapHeavy[wrapInds] ] 
      nPops = 2 # set number of populations to account for bound/wrap
    else:
      # reorder and get subIndices for this frame
      inds = [ [ mapHeavy[ subInds[t][i][j] ] 
                 for j in range(len(subInds[t][i])) ] 
               for i in range(nPops) ]

    # grab ith positions and box vectors
    thispos = np.array(frame.xyz)
    thisbox = np.reshape(np.array(frame.box.values[:3]),(1,3))
 
    thisHeavy = np.array(thispos[heavyInds])

    # calculate the voronoi cell volumes
    VV = voronoi_volumes(thisHeavy, thisbox[0,0], len(watInds))
    Vol = VV[0][:len(watInds)]
    Area = VV[1][:len(watInds)]

    # get volume and area properties of all waters
    watVol[0].append( Vol[~np.isinf(Vol)] ); watArea[0].append( Area[~np.isinf(Area)] )
    watEta[0].append( Area[~np.isinf(Area)]**3.0/36.0/np.pi / 
                      Vol[~np.isinf(Vol)]**2.0 )

    avgVol[0][t] = np.mean( Vol[~np.isinf(Vol)] ); varVol[0][t] = np.var( Vol[~np.isinf(Vol)] )
    avgArea[0][t] = np.mean( Area[~np.isinf(Area)] ); varArea[0][t] = np.var( Area[~np.isinf(Area)] )
    avgEta[0][t] = np.mean( Area[~np.isinf(Area)]**3.0/36.0/np.pi/Vol[~np.isinf(Vol)]**2.0 )
    varEta[0][t] = np.var( Area[~np.isinf(Area)]**3.0/36.0/np.pi/Vol[~np.isinf(Vol)]**2.0 )

    # loop through subPositions
    for j in range(1, (nPops+1)):
      watVol[j].append( Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])] )
      watArea[j].append( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])] )
      watEta[j].append( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])]**3.0/36.0/np.pi / 
                        Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])]**2.0 )

      avgVol[j][t] = np.mean( Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])] )
      varVol[j][t] = np.var( Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])] )

      avgArea[j][t] = np.mean( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])] )
      varArea[j][t] = np.var( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])] )
      
      avgEta[j][t] = np.mean( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])]**3.0/36.0/np.pi / 
                              Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])]**2.0 )
      varEta[j][t] = np.var( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])]**3.0/36.0/np.pi /
                             Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])]**2.0 ) 

      
  # estimate simulation average and CIs using block averaging
  avgVol_mean = np.zeros((nPops+1)); avgVol_CI = np.zeros((nPops+1))
  varVol_mean = np.zeros((nPops+1)); varVol_CI = np.zeros((nPops+1))
  avgArea_mean = np.zeros((nPops+1)); avgArea_CI = np.zeros((nPops+1))
  varArea_mean = np.zeros((nPops+1)); varArea_CI = np.zeros((nPops+1))
  avgEta_mean = np.zeros((nPops+1)); avgEta_CI = np.zeros((nPops+1))
  varEta_mean = np.zeros((nPops+1)); varEta_CI = np.zeros((nPops+1))
 
  for j in range((nPops+1)):
    avgVol_CI[j] = blockAverage(avgVol[j][:]); avgVol_mean[j] = np.mean(avgVol[j][:])
    varVol_CI[j] = blockAverage(varVol[j][:]); varVol_mean[j] = np.mean(varVol[j][:])
    avgArea_CI[j] = blockAverage(avgArea[j][:]); avgArea_mean[j] = np.mean(avgArea[j][:])
    varArea_CI[j] = blockAverage(varArea[j][:]); varArea_mean[j] = np.mean(varArea[j][:])
    avgEta_CI[j] = blockAverage(avgEta[j][:]); avgEta_mean[j] = np.mean(avgEta[j][:])
    varEta_CI[j] = blockAverage(varEta[j][:]); varEta_mean[j] = np.mean(varEta[j][:])

    VolDist, bins = np.histogram(np.concatenate(watVol[j]), bins=500,
                                 range=[10.0, 60.0],
                                 density=False)

    # save the file!                                                            
    np.savetxt('VolDistribution_'+str(j)+'.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
                                                           VolDist],axis=1),
               header='water volume (A^3)    frequency',fmt="%.3e")
    
    AreaDist, bins = np.histogram(np.concatenate(watArea[j]), bins=500,
                                  range=[10.0, 100.0],
                                  density=False)
                                                 
    # save the file!                                                            
    np.savetxt('AreaDistribution_'+str(j)+'.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
                                                            AreaDist],axis=1),
               header='water area (A^2)    frequency',fmt="%.3e")
    
    EtaDist, bins = np.histogram(np.concatenate(watEta[j]), bins=500,
                                 range=[1.00, 2.5],
                                 density=False)

    # save the file!                                                            
    np.savetxt('EtaDistribution_'+str(j)+'.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
                                                           EtaDist],axis=1),
               header='asphericity    frequency',fmt="%.3e")


  avgVol = [avgVol_mean, avgVol_CI]; varVol = [varVol_mean, varVol_CI]
  avgArea = [avgArea_mean, avgArea_CI]; varArea = [varArea_mean, varArea_CI]
  avgEta = [avgEta_mean, avgEta_CI]; varEta = [varEta_mean, varEta_CI]

  return avgVol, varVol, avgArea, varArea, avgEta, varEta

def hydratedVolumeCalc(topFile, trajFile, subInds=None, nPops=0, solResName='(!:WAT)', watResName='(:WAT)', stride=1):
  """Given a topology file and trajectory (or a list of) file, computes the effective hydrated molecular volume of the non-water molecules. Here, solute atoms will be the set of heavy atoms.
     Inputs:
             topFile - topology file (path from cwd)
             trajFile - trajectory file, or list of them (path from cwd)
             solResName - string defining the residue name for the non-water cosolvent. Default='(!:WAT)'
             watResName - string defining the residue name for the water. Default='(:WAT)'
             subInds - list containing water indices for n-populations at each timestep
                       (e.g., [ [boundInds, wrapInds]_{0}, ..., [boundInds, wrapInds]_{t} ])

             nPops - the number of populations in subInds (Default=1)
             stride - skip every "stride" steps of the trajectory (int). Default=1
     Outputs:
             avgVol - estimate of the mean hydrated volume (Angstrom^3)
             varVol - estimate of the var. of the hydrated volume (Angstrom^3)

  """
  # load topologies and trajectories using the above class AnalysisObject
  obj = TrajObject(topFile, trajFile, stride, solResName, watResName)
  top = obj.top #loadTop()
  traj = obj.traj #.loadTraj()

  # select water oxygen and hydrogen indices
  watInds, watHInds, lenWat = obj.getWatInds()

  # select cosolvent heavy
  #solInds, solHInds, solCInds, solNInds, solOInds, solSInds = obj.getSolInds()
  solInds, solHInds, _, _, _, _ = obj.getSolInds()

  # get all heavy indices
  heavyInds = np.concatenate((watInds,solInds))

  # define opposite ordering
  orderedInds = np.concatenate((solInds, watInds))

  # create mapping function from water index to heavyInds equivalent
  mapHeavy = {watInds[i]:i for i in range(len(watInds))}

  # Create lists to store volume and area values.
  molVol = [ [] for i in range((nPops+1)) ]

  # Lists to store mean and var. volume, area, and asphericity
  avgVol = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  varVol = [ np.zeros(len(traj)) for i in range((nPops+1)) ]

  # loop through trajectory to compute average RDFs
  for t, frame in enumerate(traj):

    if subInds is None:
      inds = [ [ mapHeavy[ ind ] for ind in watInds ] ]
    else:
      # reorder and get subIndices for this frame
      inds = [ [ mapHeavy[ subInds[t][i][j] ] 
                 for j in range(len(subInds[t][i])) ] 
               for i in range(nPops) ]

    # grab ith positions and box vectors
    thispos = np.array(frame.xyz)
    thisbox = np.reshape(np.array(frame.box.values[:3]),(1,3))
 
    thisHeavy = np.array(thispos[heavyInds])
    thisOrdered = np.array(thispos[orderedInds])

    # get contact list with the non-water cosolvent
    contacts, _, _, _ = sl.voronoi_contacts(thisHeavy, 
                                            thisbox[0][0], len(solInds))

    print(contacts)
    print(contacts.shape)
    print(np.sum(contacts))
    stop

    # calculate the voronoi cell volumes
    VV = voronoi_volumes(thisHeavy, thisbox[0,0], len(watInds))
    Vol = VV[0][:len(watInds)]
    Area = VV[1][:len(watInds)]

    # get volume and area properties of all waters
#    watVol[0].append( Vol[~np.isinf(Vol)] ); watArea[0].append( Area[~np.isinf(Area)] )
#    watEta[0].append( Area[~np.isinf(Area)]**3.0/36.0/np.pi / 
#                      Vol[~np.isinf(Vol)]**2.0 )

#    avgVol[0][t] = np.mean( Vol[~np.isinf(Vol)] ); varVol[0][t] = np.var( Vol[~np.isinf(Vol)] )
#    avgArea[0][t] = np.mean( Area[~np.isinf(Area)] ); varArea[0][t] = np.var( Area[~np.isinf(Area)] )
#    avgEta[0][t] = np.mean( Area[~np.isinf(Area)]**3.0/36.0/np.pi/Vol[~np.isinf(Vol)]**2.0 )
#    varEta[0][t] = np.var( Area[~np.isinf(Area)]**3.0/36.0/np.pi/Vol[~np.isinf(Vol)]**2.0 )

    # loop through subPositions
#    for j in range(1, (nPops+1)):
#      watVol[j].append( Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])] )
#      watArea[j].append( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])] )
#      watEta[j].append( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])]**3.0/36.0/np.pi / 
#                        Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])]**2.0 )

#      avgVol[j][t] = np.mean( Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])] )
#      varVol[j][t] = np.var( Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])] )

#      avgArea[j][t] = np.mean( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])] )
#      varArea[j][t] = np.var( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])] )
      
#      avgEta[j][t] = np.mean( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])]**3.0/36.0/np.pi / 
#                              Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])]**2.0 )
#      varEta[j][t] = np.var( Area[inds[j-1]][~np.isinf(Area[inds[j-1]])]**3.0/36.0/np.pi /
#                             Vol[inds[j-1]][~np.isinf(Vol[inds[j-1]])]**2.0 ) 

      
  # estimate simulation average and CIs using block averaging
#  avgVol_mean = np.zeros((nPops+1)); avgVol_CI = np.zeros((nPops+1))
#  varVol_mean = np.zeros((nPops+1)); varVol_CI = np.zeros((nPops+1))
#  avgArea_mean = np.zeros((nPops+1)); avgArea_CI = np.zeros((nPops+1))
#  varArea_mean = np.zeros((nPops+1)); varArea_CI = np.zeros((nPops+1))
#  avgEta_mean = np.zeros((nPops+1)); avgEta_CI = np.zeros((nPops+1))
#  varEta_mean = np.zeros((nPops+1)); varEta_CI = np.zeros((nPops+1))
 
#  for j in range((nPops+1)):
#    avgVol_CI[j] = blockAverage(avgVol[j][:]); avgVol_mean[j] = np.mean(avgVol[j][:])
#    varVol_CI[j] = blockAverage(varVol[j][:]); varVol_mean[j] = np.mean(varVol[j][:])
#    avgArea_CI[j] = blockAverage(avgArea[j][:]); avgArea_mean[j] = np.mean(avgArea[j][:])
#    varArea_CI[j] = blockAverage(varArea[j][:]); varArea_mean[j] = np.mean(varArea[j][:])
#    avgEta_CI[j] = blockAverage(avgEta[j][:]); avgEta_mean[j] = np.mean(avgEta[j][:])
#    varEta_CI[j] = blockAverage(varEta[j][:]); varEta_mean[j] = np.mean(varEta[j][:])

#    VolDist, bins = np.histogram(np.concatenate(watVol[j]), bins=500,
#                                 range=[10.0, 60.0],
#                                 density=False)

    # save the file!                                                            
#    np.savetxt('VolDistribution_'+str(j)+'.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
#                                                           VolDist],axis=1),
#               header='water volume (A^3)    frequency',fmt="%.3e")
    
#    AreaDist, bins = np.histogram(np.concatenate(watArea[j]), bins=500,
#                                  range=[10.0, 100.0],
#                                  density=False)
                                                 
    # save the file!                                                            
#    np.savetxt('AreaDistribution_'+str(j)+'.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
#                                                            AreaDist],axis=1),
#               header='water area (A^2)    frequency',fmt="%.3e")
    
#    EtaDist, bins = np.histogram(np.concatenate(watEta[j]), bins=500,
#                                 range=[1.00, 2.5],
#                                 density=False)

    # save the file!                                                           
#    np.savetxt('EtaDistribution_'+str(j)+'.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
#                                                           EtaDist],axis=1),
#               header='asphericity    frequency',fmt="%.3e")


#  avgVol = [avgVol_mean, avgVol_CI]; varVol = [varVol_mean, varVol_CI]
#  avgArea = [avgArea_mean, avgArea_CI]; varArea = [varArea_mean, varArea_CI]
#  avgEta = [avgEta_mean, avgEta_CI]; varEta = [varEta_mean, varEta_CI]

  return 0

def threeBodyCalc(topFile, trajFile, subInds=None, nPops=0, solResName='(!:WAT)', watResName='(:WAT)', nBins=500, stride=1, output2D=False):
  """Given a topology file and trajectory (or a list of) file, computes the three-body angle distribution for the system-aveage and for as specific population (i.e. hydration shell, bound, or wrap)
     Inputs:                                                                                                        
             topFile - topology file (path from cwd)                                                                
             trajFile - trajectory file, or list of them (path from cwd)                                            
             solResName - string defining the residue name for the non-water cosolvent. Default='(!:WAT)'               
             watResName - string defining the residue name for the water. Default='(:WAT)'                             
             subInds - list containing water indices for n-populations at each timestep
                       (e.g., [ [boundInds, wrapInds]_{0}, ..., [boundInds, wrapInds]_{t} ])
             stride - skip every "stride" steps of the trajectory (int). Default=1                                  
             output2D - switch to save off wat-wat coords and generate 2D histogram of (theta, Nc) space (default=False)
     Outputs lists of distribution properties for each population(e.g., [i_bound, i_wrap]):   
             fracTet - pop. of tetrahedral waters                             
             avgCos - avg. cos from tet. portion of the distribution 
             varCos - var. cos from tet. portion of the distribution
             entropy - distribution entropy
             nWats - number of waters in each population
  """
  # load topologies and trajectories using the above class AnalysisObject                                           
  obj = TrajObject(topFile, trajFile, stride, solResName, watResName)
  top = obj.top #loadTop()
  traj = obj.traj #loadTraj()

  # select water oxygen and hydrogen indices                        
  watInds, watHInds, lenWat = obj.getWatInds()

  # select cosolvent heavy atoms
  solInds, solHInds, solCInds, solNInds, solOInds, solSInds = obj.getSolInds()

  # set up empty lists to hold 3body Dist. for water populations
  angles = [ [] for i in range((nPops+1)) ]
  
  # set up empty list to hold water-water coordination numbers for 3-body angles
  numbers = []
  angles_trial = []
  
  # lists to store distribution statistics for each water population
  nWats = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  pTet = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  avgCos = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  varCos = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  entropy = [ np.zeros(len(traj)) for i in range((nPops+1)) ]

  for t,frame in enumerate(traj):
    # get water positions and box info
    pos = np.array(frame.xyz)
    thisbox = np.array(frame.box.values[:3])
    watPos = pos[watInds]

    if subInds is None:
      subPos = watPos
    else:
      # grab sub-water population positions
      subPos = [ pos[subInds[t][i]] for i in range(nPops) ]

    # compute 3body-dists for all waters and compute statistics
    [jAngles, jNumbers] = wp.getCosAngs(watPos, watPos, thisbox)

    angles[0].append(jAngles)
    nWats[0][t] = len(watInds)
    if output2D:
      count = 0
      for n in jNumbers:
        count=int(n-1)
        while count>0:
          numbers.append( [int(n-1) for k in range(count)] )
          count = count - 1        
   
    if len(angles[0])!=0:
      angDist, bins, a, b, c, d = wp.tetrahedralMetrics(jAngles, nBins=nBins)
      pTet[0][t] = a; avgCos[0][t] = b; varCos[0][t] = c; entropy[0][t] = d
    else:
      pTet[0][t] = 0; avgCos[0][t] = 0; varCos[0][t] = 0; entropy[0][t] = 0

    # compute 3body-dists for each population and compute statistics
    for j in range(1, (nPops+1)):
      [jAngles, jNumbers] = wp.getCosAngs(subPos[j-1], watPos, thisbox)
      angles[j].append(jAngles)
      nWats[j][t] = len(subInds[t][j-1]) 

      if len(angles[j])!=0:
        angDist, bins, a, b, c, d = wp.tetrahedralMetrics(jAngles, nBins=nBins)
        pTet[j][t] = a; avgCos[j][t] = b; varCos[j][t] = c; entropy[j][t] = d
      else:
        pTet[j][t] = 0; avgCos[j][t] = 0; varCos[j][t] = 0; entropy[j][t] = 0

  nWats_mean = np.zeros((nPops+1)); nWats_CI = np.zeros((nPops+1))
  pTet_mean = np.zeros((nPops+1)); pTet_CI = np.zeros((nPops+1))
  avgCos_mean = np.zeros((nPops+1)); avgCos_CI = np.zeros((nPops+1))
  varCos_mean = np.zeros((nPops+1)); varCos_CI = np.zeros((nPops+1))
  entropy_mean = np.zeros((nPops+1)); entropy_CI = np.zeros((nPops+1))
  for j in range((nPops+1)):
    nWats_CI[j] = blockAverage(nWats[j][:]); nWats_mean[j] = np.mean(nWats[j][:])
    pTet_CI[j] = blockAverage(pTet[j][:]); pTet_mean[j] = np.mean(pTet[j][:])
    avgCos_CI[j] = blockAverage(avgCos[j][:]); avgCos_mean[j] = np.mean(avgCos[j][:])
    varCos_CI[j] = blockAverage(varCos[j][:]); varCos_mean[j] = np.mean(varCos[j][:])
    entropy_CI[j] = blockAverage(entropy[j][:]); entropy_mean[j] = np.mean(entropy[j][:])

  # combine mean and CI for each statistic
  pTet = [pTet_mean, pTet_CI]
  avgCos = [avgCos_mean, avgCos_CI]
  varCos = [varCos_mean, varCos_CI]
  entropy = [entropy_mean, entropy_CI]
  nWats = [nWats_mean, nWats_CI]

  # compute (avg.) sub-poulation 3-body distributions
  for j in range((nPops+1)):
    angles[j] = np.concatenate(angles[j])
    if len(angles[j])!=0:
      angDist, bins, a, b, c, d = wp.tetrahedralMetrics(angles[j], 
                                                        nBins=nBins)
      np.savetxt('3bDistribution_'+str(j)+'.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
                                                            angDist],axis=1),
                 header='3-body angle (deg)    frequency',fmt="%.3e")

  # if saving (theta,Nc) 2D-histogram, compute and save
  if output2D:
    numbers = np.concatenate(numbers).astype(float)
    
    # pre-define x and y edges
    xedges = np.arange(-1.5, 13.5, 1)
    yedges = np.linspace(0, 180, 500)
    H, xedges, yedges = np.histogram2d(numbers, angles[0],
                                       #angles[0], numbers, 
                                       bins=(xedges, yedges))

    H = H/np.sum(H)
        
    fig, ax = plt.subplots(figsize=(4,4))
    xcenters = (xedges[:-1] + xedges[1:]) / 2
    ycenters = (yedges[:-1] + yedges[1:]) / 2

    ax.imshow(H, interpolation='gaussian', cmap='viridis', aspect='auto',
              origin='lower',
              extent=(ycenters[0], ycenters[-1],
                      xcenters[1], xcenters[-1])
    )

    print('the mean three-body angle is {}'.format(np.mean(angles[0])) )
    print('the mean water-water coord. is {}'.format(np.mean(numbers)) )

    ax.axvline(x=50, linewidth=1, linestyle='dashed', color='white')
    ax.axvline(x=np.mean(angles[0]), 
               linewidth=1, linestyle='dashed', color='grey')
    ax.axvline(x=130, linewidth=1, linestyle='dashed', color='darkorange')
    ax.axhline(y=4, linewidth=1, linestyle='dashed', color='magenta')

    ax.set_xticks(np.arange(0, 200, 20))
    ax.set_yticks(np.arange(0, 13, 2))
    ax.set_xlim([0, 180])
    ax.set_ylim([0, 12])
    ax.set_xlabel(r'$\theta [^{\circ}]$')
    ax.set_ylabel(r'$N_{c}$')
    plt.savefig('3bDistribution_2D.png')

  return pTet, avgCos, varCos, entropy, nWats

def tetOrderCalc(topFile, trajFile, subInds=None, nPops=0, solResName='(!:WAT)', watResName='(:WAT)', stride=1):
  """Given a topology file and trajectory (or a list of) file, computes the tet. order parameter distribution for the system-aveage and for as specific population (i.e. hydration shell, bound, or wrap)
     Inputs:                                                                                                        
             topFile - topology file (path from cwd)                                                                
             trajFile - trajectory file, or list of them (path from cwd)                                            
             solResName - string defining the residue name for the non-water cosolvent. Default='(!:WAT)'               
             watResName - string defining the residue name for the water. Default='(:WAT)'                             
             subInds - list containing water indices for n-populations at each timestep
                       (e.g., [ [boundInds, wrapInds]_{0}, ..., [boundInds, wrapInds]_{t} ])
             stride - skip every "stride" steps of the trajectory (int). Default=1                                  
     Outputs lists of distribution properties for each population(e.g., [i_bound, i_wrap]):   
             avgQ - avg. tet. order param 
             varQ - var. tet. order param
  """
  # load topologies and trajectories using the above class AnalysisObject                   
  obj = TrajObject(topFile, trajFile, stride, solResName, watResName)
  top = obj.top #loadTop()
  traj = obj.traj #loadTraj()

  # select water oxygen and hydrogen indices                        
  watInds, watHInds, lenWat = obj.getWatInds()

  # select cosolvent heavy atoms
  solInds, solHInds, solCInds, solNInds, solOInds, solSInds = obj.getSolInds()

  # set up empty lists to hold q-Dist. for water populations
  qVals = [ [] for i in range((nPops+1)) ]
  
  # lists to store distribution statistics for each water population
  avgQ = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  varQ = [ np.zeros(len(traj)) for i in range((nPops+1)) ]

  for t,frame in enumerate(traj):
    # get water positions and box info
    pos = np.array(frame.xyz)
    thisbox = np.array(frame.box.values[:3])
    watPos = pos[watInds]

    if subInds is None:
      subPos = watPos
    else:
      # grab sub-water population positions
      subPos = [ pos[subInds[t][i]] for i in range(nPops) ]

    # compute q-dists for all waters and compute statistics
    jQ = wp.getOrderParamq(watPos, watPos, thisbox)
    qVals[0].append(jQ)
    avgQ[0][t] = np.mean(jQ); varQ[0][t] = np.var(jQ)

    # compute q-dists for each population and compute statistics
    for j in range(1, (nPops+1)):
      jQ = wp.getOrderParamq(subPos[j-1], watPos, thisbox)
      qVals[j].append(jQ)

      avgQ[j][t] = np.mean(jQ); varQ[j][t] = np.var(jQ)

  avgQ_mean = np.zeros((nPops+1)); avgQ_CI = np.zeros((nPops+1))
  varQ_mean = np.zeros((nPops+1)); varQ_CI = np.zeros((nPops+1))
  for j in range((nPops+1)):
    avgQ_CI[j] = blockAverage(avgQ[j][:]); avgQ_mean[j] = np.mean(avgQ[j][:])
    varQ_CI[j] = blockAverage(varQ[j][:]); varQ_mean[j] = np.mean(varQ[j][:])

  # combine mean and CI for each statistic
  avgQ = [avgQ_mean, avgQ_CI]
  varQ = [varQ_mean, varQ_CI]

  # compute (avg.) sub-poulation q-distributions
  for j in range((nPops+1)):
    qVals[j] = np.concatenate(qVals[j])
    qDist, bins = np.histogram(qVals[j], bins=500,
                               range=[0.0, 1.0],
                               density=False)

    np.savetxt('qDistribution_'+str(j)+'.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
                                                         qDist],axis=1),
               header='qVal    frequency',fmt="%.3e")

  return avgQ, varQ

def hexOrderCalc(topFile, trajFile, subInds=None, nPops=0, solResName='(!:WAT)', endResName='(:WAT)', stride=1, lowCut=0.0, highCut=7.0):
  """Given a topology file and trajectory (or a list of) file, computes the hexagonal. order parameter distribution for the system-aveage and for surface atoms
     Inputs:                                                                                                        
             topFile - topology file (path from cwd)                                                                
             trajFile - trajectory file, or list of them (path from cwd)                                            
             solResName - string defining the residue name for the non-water cosolvent. Default='(!:WAT)'               
             endResName - string defining the residue name for the chain ends. Default='(:WAT)'                             
             subInds - list containing water indices for n-populations at each timestep
                       (e.g., [ [boundInds, wrapInds]_{0}, ..., [boundInds, wrapInds]_{t} ])
             stride - skip every "stride" steps of the trajectory (int). Default=1                                  
     Outputs lists of distribution properties for each population(e.g., [i_bound, i_wrap]):   
             avgPsi - avg. psi-6 order param 
             varPsi - var. psi-6 order param
  """
  # load topologies and trajectories using the above class AnalysisObject                   
  obj = TrajObject(topFile, trajFile, stride, solResName, endResName)
  top = obj.top #loadTop()
  traj = obj.traj #loadTraj()

  # select water oxygen and hydrogen indices                        
  endInds, endHInds, lenEnd = obj.getWatInds()
  endInds = endInds[1::2]

  # select cosolvent heavy atoms
  solInds, solHInds, solCInds, solNInds, solOInds, solSInds = obj.getSolInds()

  # set up empty lists to hold psi-Dist. for water populations
  psiVals = [ [] for i in range((nPops+1)) ]
  
  # lists to store distribution statistics for each water population
  avgPsi = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  varPsi = [ np.zeros(len(traj)) for i in range((nPops+1)) ]

  for t,frame in enumerate(traj):
    # get water positions and box info
    pos = np.array(frame.xyz)
    thisbox = np.array(frame.box.values[:3])
    endPos = pos[endInds]

    if subInds is None:
      subPos = endPos
    else:
      # grab sub-end population positions
      subPos = [ pos[subInds[t][i]] for i in range(nPops) ]

    # compute psi6-dists for all end and compute statistics
    jPsi = wp.getOrderParamPsi(endPos, endPos, thisbox, 
                               lowCut=0.0, highCut=7.0)
    psiVals[0].append(jPsi)
    avgPsi[0][t] = np.mean(jPsi); varPsi[0][t] = np.var(jPsi)

    # compute psi-dists for each population and compute statistics
    for j in range(1, (nPops+1)):
      jPsi = wp.getOrderParamPsi(subPos[j-1], endPos, thisbox)
      psiVals[j].append(jPsi)

      avgPsi[j][t] = np.mean(jPsi); varPsi[j][t] = np.var(jPsi)

  avgPsi_mean = np.zeros((nPops+1)); avgPsi_CI = np.zeros((nPops+1))
  varPsi_mean = np.zeros((nPops+1)); varPsi_CI = np.zeros((nPops+1))
  for j in range((nPops+1)):
    avgPsi_CI[j] = blockAverage(avgPsi[j][:]); avgPsi_mean[j] = np.mean(avgPsi[j][:])
    varPsi_CI[j] = blockAverage(varPsi[j][:]); varPsi_mean[j] = np.mean(varPsi[j][:])

  # combine mean and CI for each statistic
  avgPsi = [avgPsi_mean, avgPsi_CI]
  varPsi = [varPsi_mean, varPsi_CI]

  # compute (avg.) sub-poulation psi6-distributions
  for j in range((nPops+1)):
    psiVals[j] = np.concatenate(psiVals[j])
    psiDist, bins = np.histogram(psiVals[j], bins=500,
                                 range=[0.0, 1.0],
                                 density=False)

    np.savetxt('psiDistribution_'+str(j)+'.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
                                                           psiDist],axis=1),
               header='psiVal    frequency',fmt="%.3e")

  return avgPsi, varPsi

def lsiCalc(topFile, trajFile, subInds=None, nPops=0, solResName='(!:WAT)', watResName='(:WAT)', stride=1):
  """Given a topology file and trajectory (or a list of) file, computes the local structure index (LSI) distribution for the system-aveage and for as specific population (i.e. hydration shell, bound, or wrap)
     Inputs:                                                                                                        
             topFile - topology file (path from cwd)                                                                
             trajFile - trajectory file, or list of them (path from cwd)                                            
             solResName - string defining the residue name for the non-water cosolvent. Default='(!:WAT)'               
             watResName - string defining the residue name for the water. Default='(:WAT)'                             
             subInds - list containing water indices for n-populations at each timestep
                       (e.g., [ [boundInds, wrapInds]_{0}, ..., [boundInds, wrapInds]_{t} ])
             stride - skip every "stride" steps of the trajectory (int). Default=1                                  
     Outputs lists of distribution properties for each population(e.g., [i_bound, i_wrap]):   
             avgLSI - avg. LSI 
             varLSI - var. LSI
  """
  # load topologies and trajectories using the above class AnalysisObject                   
  obj = TrajObject(topFile, trajFile, stride, solResName, watResName)
  top = obj.top #loadTop()
  traj = obj.traj #loadTraj()

  # select water oxygen and hydrogen indices                        
  watInds, watHInds, lenWat = obj.getWatInds()

  # select cosolvent heavy atoms
  solInds, solHInds, solCInds, solNInds, solOInds, solSInds = obj.getSolInds()

  # set up empty lists to hold LSI Dist. for water populations
  lsiVals = [ [] for i in range((nPops+1)) ]
  
  # lists to store distribution statistics for each water population
  avgLSI = [ np.zeros(len(traj)) for i in range((nPops+1)) ]
  varLSI = [ np.zeros(len(traj)) for i in range((nPops+1)) ]

  for t,frame in enumerate(traj):
    # get water positions and box info
    pos = np.array(frame.xyz)
    thisbox = np.array(frame.box.values[:3])
    watPos = pos[watInds]

    if subInds is None:
      subPos = watPos
    else:
      # grab sub-water population positions
      subPos = [ pos[subInds[t][i]] for i in range(nPops) ]

    # compute LSI-dists for all waters and compute statistics
    [jLSI, numLSI] = wp.getLSI(watPos, watPos, thisbox)
    lsiVals[0].append(jLSI)
    avgLSI[0][t] = np.mean(jLSI); varLSI[0][t] = np.var(jLSI)

    # compute LSI-dists for each population and compute statistics
    for j in range(1, (nPops+1)):
      [jLSI, numLSI] = wp.getLSI(subPos[j-1], watPos, thisbox, lowCut=0.0, highCut=3.7)
      lsiVals[j].append(jLSI)

      avgLSI[j][t] = np.mean(jLSI); varLSI[j][t] = np.var(jLSI)

  avgLSI_mean = np.zeros((nPops+1)); avgLSI_CI = np.zeros((nPops+1))
  varLSI_mean = np.zeros((nPops+1)); varLSI_CI = np.zeros((nPops+1))
  for j in range((nPops+1)):
    avgLSI_CI[j] = blockAverage(avgLSI[j][:]); avgLSI_mean[j] = np.mean(avgLSI[j][:])
    varLSI_CI[j] = blockAverage(varLSI[j][:]); varLSI_mean[j] = np.mean(varLSI[j][:])

  # combine mean and CI for each statistic
  avgLSI = [avgLSI_mean, avgLSI_CI]
  varLSI = [varLSI_mean, varLSI_CI]

  # compute (avg.) sub-poulation LSI distributions
  for j in range((nPops+1)):
    lsiVals[j] = np.concatenate(lsiVals[j])
    lsiDist, bins = np.histogram(lsiVals[j], bins=500,
                               range=[0.0, 0.3],
                               density=False)

    np.savetxt('lsiDistribution_'+str(j)+'.txt', np.stack([0.5*(bins[:-1]+bins[1:]),
                                                          lsiDist],axis=1),
               header='lsiVal [A^2]    frequency',fmt="%.3e")

  return avgLSI, varLSI

### Temporarily, a fixed volume HS is inserted. Should implement a function that calculates context-dependent HS radii as done in Jacob's solute_insertions.py script!
def chemPotCalc(topFile, trajFile, solResName='(!:WAT)', watResName='(:WAT)', probeRadius=3.3, keyword=False, stride=1):
  """Given a topology file and trajectory (or a list of) file, computes the tet. order param. distribution for the system-aveage and for as specific population (i.e. hydration shell, bound, or wrap)
     Inputs:                         
             topFile - topology file (path from cwd)
             trajFile - trajectory file, or list of them (path from cwd)
             solResName - string defining the residue name for the non-water cosolvent. Default='(!:WAT)'               
             watResName - string defining the residue name for the water. Default='(:WAT)' 
             probeRadius - HS radius in A. Default=3.3                         
             keyword - string to pass into funciton to determine wheter to calculate shell populations. Default=False
             stride - skip every "stride" steps of the trajectory (int). Default=1                                  
     Outputs:
             mu - -log(Pv[0]) [in kT]               
             avgN - avg. N from Pv[N]                
             avgN2 - avg. N^2 from Pv[N]                
  """
  # load topologies and trajectories using the above class AnalysisObject
  obj = TrajObject(topFile, trajFile, stride, solResName, watResName)
  top = obj.top #loadTop()
  traj = obj.traj #loadTraj()

  # select water oxygen and hydrogen indices
  watInds, watHInds, lenWat = obj.getWatInds()

  # select cosolvent heavy atoms
  solInds, solHInds, solCInds, solNInds, solOInds, solSInds = obj.getSolInds()
  
  # select all heavy atomss
  heavyInds = traj.top.select('(!@H=)&(!@EPW)')

  # cutoff for the second hydration shell of the cosolvent molecules on average
  cutoff = 4.2

  # if the keyword is true, we will set up the calculation for subsets of water molecules
  if keyword:
    # set up empty lists to hold 3body Dist. for waters in cosolvent shell (defined by rdfCalc
    numOverlap = np.arange(100)
    countOverlap = np.zeros(len(numOverlap))

    for t,frame in enumerate(traj):
      pos = np.array(frame.xyz)
      thisbox = np.array(frame.box.values[:3])

      # extract water and solute positions
      watPos = pos[watInds]
      solPos = pos[solInds]
      heavyPos = pos[heavyInds]

      # generate 1000 random insertions within shell of solutes
      count = 0
      numIns = 100000
      hsPos = np.zeros((numIns,3))
      thisTotOverlap = np.zeros(numIns, dtype=int)
      overlapBool = np.zeros((numIns, len(heavyPos)), dtype=int)

      while count<numIns:
        randX = 2.0 * (np.random.random(1)-0.5) * cutoff
        randY = 2.0 * (np.random.random(1)-0.5) * cutoff
        randZ = 2.0 * (np.random.random(1)-0.5) * cutoff
        randSq = np.sqrt(randX[0]**2.0+randY[0]**2.0+randZ[0]**2.0)
        if randSq>cutoff:
          continue
        else:
          randSolPos = pos[np.random.choice(solInds)]
          hsPos[count, 0] = randSolPos[0] + randX
          hsPos[count, 1] = randSolPos[1] + randY
          hsPos[count, 2] = randSolPos[2] + randZ
          count += 1

      overlapBool += wl.nearneighbors(hsPos, heavyPos, thisbox, 0.0, probeRadius)
      thisTotOverlap += np.sum(np.array(overlapBool, dtype=bool), axis=1)
      thisBins = np.arange(np.max(thisTotOverlap) + 1)
      countOverlap[thisBins] += np.bincount(thisTotOverlap)
     
#    print('Hard-sphere solute insertion probability: %f'%(-np.log(countOverlap[0]/np.sum(countOverlap))))

    #Save the distribution to file                                            
    np.savetxt('HS-solute_overlap_hist_Shell.txt', np.vstack((numOverlap, 
                                                              countOverlap)).T,
               header='Number of non-solute atoms overlapping           Histogram count')
    muHS = -np.log(countOverlap[0]/np.sum(countOverlap))
    avgN = np.dot(numOverlap,countOverlap)/np.sum(countOverlap)
    avgN2 = np.dot(numOverlap**2.0,countOverlap)/np.sum(countOverlap)

  # otherwise, we will just get the system-average
  else:

    # set up empty lists to hold 3body Dist. for waters in cosolvent shell (defined by rdfCalc
    numOverlap = np.arange(100)
    countOverlap = np.zeros(len(numOverlap))

    for t,frame in enumerate(traj):
      pos = np.array(frame.xyz)
      thisbox = np.array(frame.box.values[:3])

      # extract water and solute positions
      watPos = pos[watInds]
      solPos = pos[solInds]
      heavyPos = pos[heavyInds]

      # generate 1000 random insertions within shell of solutes
      count = 0
      numIns = 10000
      hsPos = np.zeros((numIns,3))
      thisTotOverlap = np.zeros(numIns, dtype=int)
      overlapBool = np.zeros((numIns, len(heavyPos)), dtype=int)

      hsPos[:,0] = np.random.random(numIns) * thisbox[0]
      hsPos[:,1] = np.random.random(numIns) * thisbox[1]
      hsPos[:,2] = np.random.random(numIns) * thisbox[2]

      overlapBool += wl.nearneighbors(hsPos, heavyPos, thisbox, 0.0, probeRadius)
      thisTotOverlap += np.sum(np.array(overlapBool, dtype=bool), axis=1)
      thisBins = np.arange(np.max(thisTotOverlap) + 1)
      countOverlap[thisBins] += np.bincount(thisTotOverlap)
     
#    print('Hard-sphere solute insertion probability: %f'%(-np.log(countOverlap[0]/np.sum(countOverlap))))

    #Save the distribution to file                                            
    np.savetxt('HS-solute_overlap_hist.txt', np.vstack((numOverlap, 
                                                              countOverlap)).T,
               header='Number of non-solute atoms overlapping           Histogram count')
    muHS = -np.log(countOverlap[0]/np.sum(countOverlap))
    avgN = np.dot(numOverlap,countOverlap)/np.sum(countOverlap)
    avgN2 = np.dot(numOverlap**2.0,countOverlap)/np.sum(countOverlap)

  return muHS, avgN, avgN2

### Let's test the functionality of the protein surface characterization library here. Specifically, looking at the contact area calculation
def contactAreaCalc(topFile, trajFile, solResName='(!:WAT)', watResName='(:WAT)', stride=1):
  """Given a topology file and trajectory (or a list of) file, fraction of cosolvent (solute) surface area incontact with bound, wrap, hydrophobic, and hydrophilic atoms
     Inputs:                                                                                 topFile - topology file (path from cwd)                                                                
            trajFile - trajectory file, or list of them (path from cwd)                                            
            solResName - string defining the residue name for the non-water cosolvent. Default='(!:WAT)'               
            watResName - string defining the residue name for the water. Default='(:WAT)'                             
             stride - skip every "stride" steps of the trajectory (int). Default=1                                  
     Outputs [bound, wrap, phobic, philic]:               
             totalArea - total area occupied by each population     
             fracArea - fraction of the area occupied by each population
  """
  # load topologies and trajectories using the above class AnalysisObject       
  obj = TrajObject(topFile, trajFile, stride, solResName, watResName)  

  top = obj.top #loadTop()
  traj = obj.traj #loadTraj()

  # convert given indices to heavy index basis
  def convertHeavyInds(heavyInds, targetInds):
    newTargetInds = np.array([i for i, k in enumerate(heavyInds)
                              if k in targetInds])
    return newTargetInds

  # select heavy atom indices
  heavyInds = obj.getHeavyInds()

  # select water oxygen and hydrogen indices                        
  watInds, watHInds, lenWat = obj.getWatInds()

  # select cosolvent heavy atoms
  solInds, solHInds, solCInds, solNInds, solOInds, solSInds = obj.getSolInds()
  hbOInds, hbNInds = getHBInds(top, traj[0], solInds, solHInds, solNInds, solOInds)

  # get solute bond partners by solute index
  solRes = []
  for i, iatom in enumerate(top.atoms):
    if i in solInds:
      iSolRes = []
      ires = iatom.residue.atoms
      for j, jatom in enumerate(ires):
        if 'H' not in jatom.name:
          iSolRes.append(jatom.idx)
      iSolRes = convertHeavyInds(heavyInds, iSolRes)
      solRes.append(iSolRes)

  # select hydrophobic and hydrophilic heavy atom indices
  phobicInds = obj.getPhobicInds()
  philicInds = obj.getPhilicInds()

  # cutoff for the second hydration shell of the cosolvent molecules on average
  cutoff = 4.2

  # set up empty lists to hold number of waters in cosolvent shell (defined by rdfCalc for Ow-Ow)
  nBound = np.zeros(len(traj))
  nWrap = np.zeros(len(traj))
  nPhobic = np.zeros(len(traj))
  nPhilic = np.zeros(len(traj))

  # set up empty lists to hold 3body Dist. for waters in cosolvent shell (defined by rdfCalc
  tot = np.zeros(len(traj))
  totBound = np.zeros(len(traj))
  totWrap = np.zeros(len(traj))
  totPhobic = np.zeros(len(traj))
  totPhilic = np.zeros(len(traj))

  fracBound = np.zeros(len(traj))
  fracWrap = np.zeros(len(traj))
  fracPhobic = np.zeros(len(traj))
  fracPhilic = np.zeros(len(traj))

  # get water, phobic, and philic indices in heavy index basis
  watHeavyInds = convertHeavyInds(heavyInds, watInds)
  solHeavyInds = convertHeavyInds(heavyInds, solInds)
  phobicHeavyInds = convertHeavyInds(heavyInds, phobicInds)
  philicHeavyInds = convertHeavyInds(heavyInds, philicInds)
  
  for t,frame in enumerate(traj):
    start = timeit.time()

    # get water positions and box info
    pos = np.array(frame.xyz)
    thisbox = np.array(frame.box.values[:3])
    watPos = pos[watInds]
    heavyPos = pos[heavyInds]
    
    start = timeit.time()
    contacts, _, _, _ = sl.voronoi_contacts(heavyPos, thisbox[0], len(heavyInds))

    # use getBoundWrap() to get bound, wrap, shell, and non-shell inds
    boundInds,wrapInds,shellInds,nonShellInds=getBoundWrap(topFile, frame, 
                                                           watInds, watHInds,
                                                           solInds, solHInds,
                                                           solCInds, solOInds,
                                                           solNInds, solSInds)
      
    # convert bound, wrap, shell, non-shell to heavy inds
    boundHeavyInds = np.int64(boundInds/lenWat)
    wrapHeavyInds = np.int64(wrapInds/lenWat) 
    shellHeavyInds = np.int64(shellInds/lenWat)
    nonShellHeavyInds = np.int64(nonShellInds/lenWat) 

    # grap wrap, bound, shell, and non-shell positions
    boundPos = pos[boundInds]
    wrapPos = pos[wrapInds]
    
    def getTotArea(contacts, solInds, targetInds):
      totTarget = 0
      tot = 0
      for i, sInd in enumerate(solInds):
        iContact = contacts[sInd, :]
        tot += np.sum(iContact)/2.0
        for j, targetInd in enumerate(targetInds):
          if targetInd in solRes[i]:
            continue
          else:
            totTarget += iContact[targetInd]/2.0        
      return totTarget, tot

    # compute total area occupied by each population
    totPhobic[t], tot[t] = getTotArea(contacts, solHeavyInds, phobicHeavyInds)
    totPhilic[t], _ = getTotArea(contacts, solHeavyInds, philicHeavyInds)
    #tot = totPhobic[t] + totPhilic[t]
    totBound[t], _ = getTotArea(contacts, solHeavyInds, boundHeavyInds)
    totWrap[t], _ = getTotArea(contacts, solHeavyInds, wrapHeavyInds)

    # compute fraction of area occupied by each population
    fracPhobic[t] = totPhobic[t]/tot[t]
    fracPhilic[t] = totPhilic[t]/tot[t]
    fracBound[t] = totBound[t]/tot[t]
    fracWrap[t] = totWrap[t]/tot[t]

  tot_CI=blockAverage(tot);tot=np.mean(tot)
  totPhobic_CI=blockAverage(totPhobic);totPhobic=np.mean(totPhobic)
  totPhilic_CI=blockAverage(totPhilic);totPhilic=np.mean(totPhilic)
  totBound_CI=blockAverage(totBound);totBound=np.mean(totBound)
  totWrap_CI=blockAverage(totWrap);totWrap=np.mean(totWrap)

  fracPhobic_CI=blockAverage(fracPhobic);fracPhobic=np.mean(fracPhobic)
  fracPhilic_CI=blockAverage(fracPhilic);fracPhilic=np.mean(fracPhilic)
  fracBound_CI=blockAverage(fracBound);fracBound=np.mean(fracBound)
  fracWrap_CI=blockAverage(fracWrap);fracWrap=np.mean(fracWrap)

  totArea = [tot, totPhobic, totPhilic, totBound, totWrap]
  totArea_CI = [tot_CI, totPhobic_CI, totPhilic_CI, totBound_CI, totWrap_CI]

  fracArea = [fracPhobic, fracPhilic, fracBound, fracWrap]
  fracArea_CI = [fracPhobic_CI, fracPhilic_CI, fracBound_CI, fracWrap_CI]

  return totArea, totArea_CI, fracArea, fracArea_CI

# I'm implementing a switch here. If TRUE, run the test script! 
switch = False

if switch:
  # getting top and struct names from command line
  topFile = sys.argv[1]
  trajFile = sys.argv[2]

  stride = 100
  obj = TrajObject(topFile, trajFile, stride=stride, 
                   solResName=(':MOL'), watResName=(':WAT'))
  top = obj.top                                                      
  traj = obj.traj                                               

  avgVol, varVol, avgArea, varArea, avgEta, varEta = voronoiCalc(topFile, 
                                                                 trajFile, 
                                                                 subInds=None,
                                                                 nPops=0, 
                                                                 solResName='(!:WAT)', 
                                                                 watResName='(:WAT)', 
                                                                 stride=100)

  print(avgVol)
  print(avgArea)
  print(avgEta)
  stop

  # select water oxygen and hydrogen indices                                
  watInds, watHInds, lenWat = obj.getWatInds()

  # select cosolvent heavy                                                  
  solInds, solHInds, solCInds, solNInds, solOInds, solSInds = obj.getSolInds()

  hbOInds, hbNInds = getHBInds(top, traj[0], solInds, solHInds, solNInds, solOInds)

  # delineate solute oxygen acceptor and donor indices
  sol_acceptorOInds = hbOInds[0]
  sol_donorOInds = hbOInds[1]
  sol_donorHOInds = hbOInds[2]

  # delineate solute nitrogen acceptor and donor indices
  sol_acceptorNInds = hbNInds[0]
  sol_donorNInds = hbNInds[1]
  sol_donorHNInds = hbNInds[2]

  # get water hbInds
  hbOInds, _= getHBInds(top, traj[0], 
                        watInds, watHInds, [], watInds)
  
  # delineate solute oxygen acceptor and donor indices
  wat_acceptorOInds = hbOInds[0]
  wat_donorOInds = hbOInds[1]
  wat_donorHOInds = hbOInds[2]

  # combine acceptor and donor atoms
  acceptorInds = np.concatenate(np.array([wat_acceptorOInds, 
                                          sol_acceptorOInds, 
                                          sol_acceptorNInds]))

  donorInds = np.concatenate(np.array([wat_donorOInds,
                                       sol_donorOInds, 
                                       sol_donorNInds]))

  donorHInds = np.concatenate(np.array([wat_donorHOInds, 
                                        sol_donorHOInds, 
                                        sol_donorHNInds]))

  # cutoff for the second hydration shell of the cosolvent molecules on average
  cutoff = 4.0

  # define number of populations that we want to study
  nPops = 4

  if 'boundFile.npy' in os.listdir():
    trialInds = np.load('boundFile.npy', allow_pickle=True)
    if ( len(trialInds[0])!=nPops )|( len(trialInds)!=len(traj) ):
      os.remove('boundFile.npy')
    else:
      subInds = trialInds

  if 'boundFile.npy' not in os.listdir():
    subInds = []
    for t, frame in enumerate(traj):
      # use getBoundWrap() to get bound, wrap, shell, and non-shell inds
      boundInds, wrapInds, shellInds, nonShellInds = getBoundWrap(topFile, frame, 
                                                                  watInds, watHInds,
                                                                  solInds, solHInds,
                                                                  solCInds, solOInds,
                                                                  solNInds, solSInds, cutoff = 4.6)
    
      subInds.append([ boundInds, wrapInds, shellInds, nonShellInds ])

    np.save('boundFile.npy', np.array(subInds, dtype=object), allow_pickle=True)
  
  stop

  ### getHBMat() test
  # get indices for each residue
  resInds = [traj.top.select(":"+str(i+1)) for i in range(len(top.residues))]

  shellInds = subInds[0][2]
  
  shellResInds = [ top.atoms[shellInds[i]].residue.idx for i in range(len(shellInds)) ]

  start = timeit.time()

  meanCluster = getClusterStats(topFile, trajFile, 
                                acceptorInds, 
                                donorInds, 
                                donorHInds,
                                stride=100)
  print(meanCluster)
  stop

  # testing getClusters() on shell water subset
  hbMat = hbMat[shellResInds, :]
  hbMat = hbMat[:, shellResInds]

  clusters = getClusters(hbMat)

  nClusters = len(clusters)
  
  print(nClusters)

  meanCluster = np.mean( np.array( [len(clusters[i]) for i in range(nClusters)] ) )
  print(meanCluster)
  stop


#  avgVol, varVol, avgArea, varArea, avgEta, varEta  = voronoiCalc(topFile,
#                                                                  trajFile,
#                                                                  subInds,
#                                                                  nPops=4,
#                                                                  solResName='(:MOL)', 
#                                                                  watResName='(:WAT)', 
#                                                                  stride=stride)


#  avgLSI, varLSI = lsiCalc(topFile, 
#                           trajFile, 
#                           subInds, 
#                           nPops=4, 
#                           solResName='(:MOL)', 
#                           watResName='(:WAT)', 
#                           stride=stride)

#  avgQ, varQ = tetOrderCalc(topFile, 
#                            trajFile, 
#                            subInds, 
#                            nPops=4, 
#                            solResName='(:MOL)', 
#                            watResName='(:WAT)', 
#                            stride=stride)

#  pTet, avgCos, varCos, entropy, numbers = threeBodyCalc(topFile, 
#                                                         trajFile, 
#                                                         subInds,
#                                                         nPops = nPops,
#                                                         solResName='(:MOL)', 
#                                                         watResName='(:WAT)', 
#                                                         stride=stride)


#  test = rdfCalc(topFile, 
#                 trajFile, 
#                 solResName='(!:WAT)', 
#                 watResName='(:WAT)',
#                 stride=1)

#  print(test)

#  avgHB = hbCalc(topFile,
#                 trajFile,
#                 solResName='(!:WAT)',
#                 watResName='(:WAT)',
#                 stride=100)

#  print(avgHB)
#  stop

  totArea,totArea_CI,fracArea,fracArea_CI = contactAreaCalc(topFile, 
                                                            trajFile, 
                                                            solResName='(:MOL)&(@O=)', 
                                                            watResName='(:WAT)', 
                                                            stride=500)

  print('The total solute surface area is {} +/- {} angstroms^2 \n'.format(totArea[0], totArea_CI[0]))
  print('The total solute surface area occupied by hydrophobic contacts is {} +/- {} angstroms^2 \n'.format(totArea[1], totArea_CI[1]))
  print('The total solute surface area occupied by hydrophilic contacts is {} +/- {} angstroms^2 \n'.format(totArea[2], totArea_CI[2]))
  print('The total solute surface area occupied by bound contacts is {} +/- {} angstroms^2 \n'.format(totArea[3], totArea_CI[3]))
  print('The total solute surface area occupied by wrap contacts is {} +/- {} angstroms^2 \n'.format(totArea[4], totArea_CI[4]))

  print('The fraction of solute surface area occupied by hydrophobic contacts is {} +/- {} angstroms^2 \n'.format(fracArea[0], fracArea_CI[0]))
  print('The fraction of solute surface area occupied by hydrophilic contacts is {} +/- {} angstroms^2 \n'.format(fracArea[1], fracArea_CI[1]))
  print('The fraction of solute surface area occupied by bound contacts is {} +/- {} angstroms^2 \n'.format(fracArea[2], fracArea_CI[2]))
  print('The fraction of solute surface area occupied by wrap contacts is {} +/- {} angstroms^2 \n'.format(fracArea[3], fracArea_CI[3]))

  stop

  pTet, avgCos, varCos, entropy, numbers = threeBodyCalc(topFile, 
                                                         trajFile, 
                                                         solResName='(:MOL)', 
                                                         watResName='(:WAT)', 
                                                         stride=200)

  pTet, avgCos, varCos, entropy, numbers = threeBodyCalc(topFile, 
                                                         trajFile, 
                                                         solResName='(:MOL)', 
                                                         watResName='(:WAT)',
                                                         stride=200)

  pTet_CI = pTet[1]
  pTet = pTet[0]
  pTet_bound = pTet[0]; pTet_wrap = pTet[1]; pTet_shell = pTet[2]; 
  pTet_nonShell = pTet[3]
  
  numbers_CI = numbers[1]
  numbers = numbers[0]
  nBound = numbers[0]; nWrap = numbers[1]; nShell = numbers[2]; 
  nNonShell = numbers[3]

  print('The tetrahedral population of the bound waters is {} +/- {}\n'.format(pTet_bound, pTet_CI[0]))
  print('The tetrahedral population of the wrap waters is {} +/- {}\n'.format(pTet_wrap, pTet_CI[1]))
  print('The tetrahedral population of the shell waters is {} +/- {}\n'.format(pTet_shell, pTet_CI[2]))
  print('The tetrahedral population of the non-shell waters is {} +/- {}\n'.format(pTet_nonShell, pTet_CI[3]))

  print('The number of bound waters is {} +/- {}\n'.format(nBound, numbers_CI[0]))
  print('The number of wrap waters is {} +/- {}\n'.format(nWrap, numbers_CI[1]))
  print('The number of shell waters is {} +/- {}\n'.format(nShell, numbers_CI[2]))
  print('The number of non-shell waters is {} +/- {}\n'.format(nNonShell, numbers_CI[3]))

  stop
