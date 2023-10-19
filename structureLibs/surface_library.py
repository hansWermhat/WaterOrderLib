import sys, os
import time as timeit

import numpy as np
import pyvista as pv
import trimesh as tm
from trimesh.curvature import discrete_gaussian_curvature_measure
from trimesh.curvature import discrete_mean_curvature_measure
from stl import mesh
import parmed as pmd
import pytraj as pt
#import waterlib as wl

# testing sortlib and waterlib versions in "ProteinDev" project                                           
sys.path.append('/home/drobins/scripts/ProteinDev/fortran')
import sortlib
import waterlib as wl
import imagelib as il

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource

from scipy.spatial import Voronoi, Delaunay, ConvexHull

def genSphere():
    u = np.linspace(0, np.pi, 30)
    v = np.linspace(0, 2*np.pi, 30)
    x = np.outer(np.sin(u), np.sin(v))
    y = np.outer(np.sin(u), np.cos(v))
    z = np.outer(np.cos(u), np.ones_like(v))
    return x,y,z

def goldenSpiral(n=100):
    """Golden spiral on a unit sphere to generate points for sasa estimation
    Inputs: n (integer)- default=100
    Outputs: pos (n, 3)-array of sphere points
    """
    inds = np.arange(0,n)
    goldenRatio = (1.0+5.0**0.5)/2.0
    theta = 2.0*np.pi*inds/goldenRatio
    phi = np.arccos(1.0-2.0*(inds+0.5)/n)
    pos =  np.array([np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi),
                     np.cos(phi)])
    pos = np.reshape(pos,(n,3))
    return pos

### Generate vdwRadii                                                         
def vdwAssign(top, nonSolName=['SOL','NA', 'CL'],
              vdwC=1.70, vdwN=1.55, vdwO=1.52, vdwS=1.80):
    vdw = []
    atomNames = []
    for i, res in enumerate(top.residues):
        if res.name not in nonSolName:
            for atom in res.atoms:
                if 'C' in atom.name[0]:
                    vdw.append(vdwC)
                    atomNames.append('C')
                elif 'O' in atom.name[0]:
                    vdw.append(vdwO)
                    atomNames.append('O')
                elif 'N' in atom.name[0]:
                    vdw.append(vdwN)
                    atomNames.append('N')
                elif 'S' in atom.name[0]:
                    vdw.append(vdwS)
                    atomNames.append('S')
    return vdw, atomNames

### Additional functionalities outside of the class                           
def getBonds(top, protInds):
  """Here, we utilize parmed to more elegantly extract the indices of hydroge\
n bond partners for use with waterlib functions.                              
     Inputs:                                                                 \
                                                                              
            top - parmed topology                                             
            frame - single pytraj trajectory frame                            
            solInds - list of solute (or water) heavy indices                  
                                                                              
     Outputs:                                                                 
            hbOInds - contains list of lists of acceptorO, donorOInds, and do\
norHOInds atoms                                                               
            hbNInds - contains list of lists of acceptorN, donorNInds, and do\
norHNInds atoms                                                               
    """

  numC = np.zeros((1, len(protInds)))
  numO = np.zeros((1, len(protInds)))
  numN = np.zeros((1, len(protInds)))
  numS = np.zeros((1, len(protInds)))
  ### Manually compute and order the H-B donors and aceptors for use in calculation                                                                        
  count = 0
  # loop thru all atoms in the topology                                       
  for i,atom in enumerate(top.atoms):
      if i in protInds:
          iatoms = atom.bond_partners                                          
          # loop thru all bonded atoms                                        
          for j,jatom in enumerate(iatoms):
              if jatom.name[0]=='C':
                  numC[:, count] += 1
              elif jatom.name[0]=='O':
                  numO[:, count] += 1
              elif jatom.name[0]=='N':
                  numN[:, count] += 1
              elif jatom.name[0]=='S':
                  numS[:, count] += 1
              else:
                  continue
          count += 1
  return numC, numO, numN, numS

### Create grid enclosing protein                                             
def sasaGrid(heavyPos, thisbox, cutoff):                                  
    # get minimum and maximum heavy atom positions                            
    minX=np.min(heavyPos[:,0]);minY=np.min(heavyPos[:,1])
    minZ=np.min(heavyPos[:,2])
    maxX=np.max(heavyPos[:,0]);maxY=np.max(heavyPos[:,1])
    maxZ=np.max(heavyPos[:,2])

    # set cutoff low array and make sure cutoff high array is formatted as    
    # (1, NPos)                                                               
    cutoff_low = np.zeros((1, heavyPos.shape[0]))
    cutoff_high = np.reshape(cutoff, (1, heavyPos.shape[0]))

    # generate cubic grid around protein                                      
    nBins = 50
    xSpan = np.linspace(0.80*minX, 1.20*maxX, nBins).reshape(1,nBins)
    ySpan = np.linspace(0.80*minY, 1.20*maxY, nBins).reshape(1,nBins)
    zSpan = np.linspace(0.80*minZ, 1.20*maxZ, nBins).reshape(1,nBins)
    xSpace = xSpan[:, 1] - xSpan[:, 0]
    ySpace = ySpan[:, 1] - ySpan[:, 0]
    zSpace = zSpan[:, 1] - zSpan[:, 0]

    X, Y, Z = np.meshgrid(xSpan, ySpan, zSpan, indexing='ij')
    pos = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])

    posBool = np.zeros((pos.shape[0], pos.shape[1]))#, bool)                  
    for i in range(pos.shape[1]):
        iPos = pos[:, i]
        iPos = iPos.reshape((1,3))
        overlap = wl.nearneighbors3(iPos, heavyPos, thisbox, cutoff_low,
                                    cutoff_high)
        minInd = np.argmin(np.abs(overlap))

        #if np.sum(overlap)!=0:                                               
        posBool[:, i] = np.array([1.0, 1.0, 1.0])*overlap[:, minInd]
            #[True, True, True]) #np.ones((3, i), bool)                       

    posBool = posBool[:, :, np.newaxis]
    posBool = np.reshape(posBool, (3, nBins, nBins, nBins))
    posBool = posBool[0, :, :, :]
    verts, faces, norms, vals = measure.marching_cubes(posBool, 0,
                                                       spacing=(xSpace,
                                                                ySpace,
                                                                zSpace))

    # shift the mesh over to the actual protein position                      
    mean = np.array([minX,minY,minZ])
    verts = verts + mean*0.9
    return verts, faces

### Create an instantaneous density interface mesh around protein        
def densityGrid(heavyPos, watPos, thisbox, level=0.016,
                minFrac = 0.7):
    maxFrac = 2.0-minFrac
    # get minimum and maximum heavy atom positions                       
    minX=np.min(heavyPos[:,0]);minY=np.min(heavyPos[:,1])
    minZ=np.min(heavyPos[:,2])
    maxX=np.max(heavyPos[:,0]);maxY=np.max(heavyPos[:,1])
    maxZ=np.max(heavyPos[:,2])
    allMin = np.min(heavyPos); allMax = np.max(heavyPos)

    # generate cubic grid around protein                                      
    nBins = 81

    xSpan = np.linspace(allMin-thisbox[0,0]/2.0, 
                        allMax+thisbox[0,0]/2.0, nBins).reshape(1,nBins)
    ySpan = np.linspace(allMin-thisbox[0,0]/2.0, 
                        allMax+thisbox[0,0]/2.0, nBins).reshape(1,nBins)
    zSpan = np.linspace(allMin-thisbox[0,0]/2.0, 
                        allMax+thisbox[0,0]/2.0, nBins).reshape(1,nBins)
    xSpace = xSpan[:, 1] - xSpan[:, 0]
    ySpace = ySpan[:, 1] - ySpan[:, 0]
    zSpace = zSpan[:, 1] - zSpan[:, 0]
    xSpan = xSpan[:, :-1] + xSpace
    ySpan = ySpan[:, :-1] + ySpace
    zSpan = zSpan[:, :-1] + zSpace

    # Per W-C OG paper, use smoothing length of 2.4-A (maybe 2.8 to be consistent with SASA calc
    dens, dens_norm = wl.willarddensityfield(watPos, xSpan, ySpan, zSpan, 
                                             thisbox,
                                             smoothlen=2.4)

    # use marching cubes to find half bulk density isosurface
    verts, faces, norms, vals = measure.marching_cubes(dens, level,
                                                       spacing=(xSpace,
                                                                ySpace,
                                                                zSpace))

    # shift the mesh over to the actual protein position                     
    verts = verts - allMin
    verts = verts-0.5*np.max(verts)
    return verts, faces

### Create an instantaneous density interface mesh around protein        
def densityVoxel(heavyPos, watPos, thisbox):                                
    # get minimum and maximum heavy atom positions                            
    minX=np.min(heavyPos[:,0]);minY=np.min(heavyPos[:,1])
    minZ=np.min(heavyPos[:,2])
    maxX=np.max(heavyPos[:,0]);maxY=np.max(heavyPos[:,1])
    maxZ=np.max(heavyPos[:,2])

    # generate cubic grid around protein                                      
    nBins = 11

    # generate x,y,z around molecule
    xSpan = np.linspace(0.8*minX, 1.2*maxX, nBins).reshape(1,nBins)
    ySpan = np.linspace(0.8*minY, 1.2*maxY, nBins).reshape(1,nBins)
    zSpan = np.linspace(0.8*minZ, 1.2*maxZ, nBins).reshape(1,nBins)
    
    # get x,y,z binwidth
    xWidth = xSpan[:, 1] - xSpan[:, 0]
    yWidth = ySpan[:, 1] - ySpan[:, 0]
    zWidth = zSpan[:, 1] - zSpan[:, 0]

    # shift x,y,z span
    xSpan = xSpan[:, :-1] + xWidth
    ySpan = ySpan[:, :-1] + yWidth
    zSpan = zSpan[:, :-1] + zWidth

    # Just a plain old density grid
    dens = wl.densityfield(watPos, xSpan, ySpan, zSpan, 
                           thisbox)
    return dens

### voronoi calculation                                                       
# use this to calculate voronoi volumes in voronoiCalc                        
def voronoi_contacts(points, boxL, numPos):
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

  contacts = np.zeros((numPos, numPos))
  proteinArea = np.zeros((1, numPos))
  proteinVol = np.zeros((1, numPos))
  watArea = np.zeros((1, numPos))
  for i in range(0, numPos):
      iRegNum = v.point_region[:numPos][i]
      iInd = np.array(v.regions[iRegNum])
      proteinArea[:, i] = ConvexHull(v.vertices[iInd],
                                     qhull_options='QJ').area
      proteinVol[:, i] = ConvexHull(v.vertices[iInd],
                                    qhull_options='QJ').volume

      for j in range(i+1, numPos):
          jRegNum = v.point_region[:numPos][j]
          jInd = np.array(v.regions[jRegNum])
          coInds = np.intersect1d(iInd, jInd)
          # Qhull will fail for nPos<4...                                     
          if len(coInds)>=4:
              contacts[i, j] = ConvexHull(v.vertices[coInds],
                                          qhull_options='QJ').area
              contacts[j, i] = contacts[i, j]
          # Hence, directly calculate the triangle area for a face with       
          # 3 contacts                                                        
          elif len(coInds)==3:
              contacts[i, j] = il.trianglearea(v.vertices[coInds])
              contacts[j, i] = contacts[i, j]
          else:
              continue
      watArea[:, i] = 2.0*proteinArea[:, i] - np.sum(contacts[i, :])
  return contacts, proteinArea, watArea, proteinVol

def localConnections(heavyPos, connMat, atomNames):
    """Generate connectivity graph to indicate heavy atom connections based on voronoi cell contact areas
    Inputs: 
             heavyPos (Npos,3)- heavy positions
             connMat (Npos,NPos)-array of contact area between atom pairs
             atomNames (NPos)- list of atom names ('C', 'N', 'O' or 'S') 
    Outputs: 
             connNum (1,NPos)-array of protein-protein "connections"
             concPhobic (1,NPos)-array of fraction of hydrophobic ('C' or 'S') atoms at each heavy atom positions 
    """
    connNum = [len(np.where(connMat[i, :]!=0)[0])
               for i in range(connMat.shape[0])]
    connNum = np.array(connNum).reshape(1, connMat.shape[0])

    connNumC = np.zeros((1, len(atomNames)))
    connNumO = np.zeros((1, len(atomNames)))
    connNumN = np.zeros((1, len(atomNames)))
    connNumS = np.zeros((1, len(atomNames)))
    for i in range(len(atomNames)):
        iInds = np.where(connMat[i, :]!=0)[0]
        iNames = [atomNames[k] for k in iInds]
        iNames.append(atomNames[i])
        for j in range(len(iNames)):
            if iNames[j]=='C':
                connNumC[:, i] += 1
            elif iNames[j]=='O':
                connNumO[:, i] += 1
            elif iNames[j]=='N':
                connNumN[:, i] +=1
            elif iNames[j]=='S':
                connNumS[:, i] +=1
            else:
                continue

    concC = connNumC/(1.0+connNum)
    concO = connNumO/(1.0+connNum)
    concN = connNumN/(1.0+connNum)
    concS = connNumS/(1.0+connNum)

    # (For now) hydrophobic atoms are C and S and hydrophilic are N and O     
    concPhobic = concC+concS
    return connNum, connNumC, connNumO, connNumN, connNumS, concPhobic

def connectPlot(heavyPos, connMat, atomProp, propName='figure'):
    """Generate connectivity graph to indicate heavy atom connections based on voronoi cell contact areas
    Inputs: heavyPos (Npos,3)- heavy positions
            connMat (Npos,NPos)- array of contact area between atom pairs
            atomProp (1,NPos)- array of atom properties  
            propName (string)- property names [Default='figure']
    Outputs: None
    """
    matplotlib.rcParams.update({'figure.figsize':[10,6]})
    matplotlib.rcParams.update({'font.size':8.0})
    matplotlib.rcParams.update({'axes.labelsize':8.0})
    matplotlib.rcParams.update({'legend.fontsize':8.0})
    matplotlib.rcParams.update({'xtick.labelsize':8.0})
    matplotlib.rcParams.update({'ytick.labelsize':8.0})
    matplotlib.rcParams.update({'lines.markersize':5.0})
    matplotlib.rcParams.update({'lines.linewidth':0.25})

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    cmap = plt.cm.get_cmap('RdBu_r')

    # color-scaled scatterplot for a given local atomic property
    p = ax.scatter(heavyPos[:, 0], heavyPos[:, 1], heavyPos[:, 2],
                   c=atomProp,
                   vmin=np.min(atomProp), vmax=np.max(atomProp),
                   cmap=cmap)

    # generate connectivity graph edges
    for i in range(connMat.shape[0]):
        iPos = heavyPos[i, :]
        for j in range(i, connMat.shape[0]):
            jPos = heavyPos[j, :]
            if connMat[i, j]!=0:
                tempPos = np.array([iPos, jPos]).reshape((2, 3))
                ax.plot3D(tempPos[:, 0], tempPos[:, 1], tempPos[:, 2],
                          color='black', marker=None)
                
    fig.colorbar(p, ax=ax, fraction=0.025, pad=0.0, location='left')
    plt.savefig(propName+'.png')
    return

### Check for overlaps between surface insertions and other solute heavy atoms
def sasaCalc(heavyPos, thisbox, vdwRadii, solRadius=1.4):
    n = 100
    unitPos = goldenSpiral(n)
    sasaPos = []
    # temporarily just look at the first protein heavy atom                   
    #subHeavyPos = np.reshape(heavyPos[510, :], (1, 3))                       
    sasa = np.zeros(heavyPos.shape[0])
    for i, iPos in enumerate(heavyPos):
        insPos = (vdwRadii[i]+solRadius)*np.array([unitPos])
        insPos += np.reshape(iPos, (1,3))
        insPos = insPos[0, :, :]
        overlapBool = np.zeros((1,insPos.shape[0]))
        for j in range(heavyPos.shape[0]):
            jPos = np.reshape(heavyPos[j, :], (1,3))
            if j==i:
                continue
            else:
                overlapBool += wl.nearneighbors(jPos, insPos,
                                                thisbox, 0,
                                                vdwRadii[j])

        inds = np.where(overlapBool==0)[1]
        sasaPos.append(insPos[inds, :])
        sasa[i] = (len(inds)/n)*4.0*np.pi*(solRadius+vdwRadii[i])

    inds = []
    for i, atom in enumerate(sasaPos):
        inds.append(i*np.ones(atom.shape[0]))

    return sasaPos, sasa, inds

### Generate a 3D representation of the SASA for a given frame
def sasaPlot(heavyPos, thisbox, vdwRadii, watRadius=1.4):
    """ Given a set of protein heavy atom positiions and vdw radii, generate a sasa mesh
    Input: heavyPos - (NPos, 3)-array of protein heavy atom positons
           thisbox - (1, 3)-array of the box vectors
           vdwRaii - (NPos)-list of vdw Radii
    Outputs: None
    """
    figure = plt.figure()
    ax = figure.add_subplot(111, projection='3d')
    
    ### Using marching_cubes, compute sasa faces and vertices                
    verts, faces = sasaGrid(heavyPos, thisbox,
                            cutoff=np.array(vdwRadii)+watRadius)

    ### Using trimesh, compute the gaussian curvature at each vertex
    tmesh = tm.Trimesh(vertices=verts, faces=faces)
    gauss = discrete_gaussian_curvature_measure(mesh=tmesh,
                                                points=tmesh.vertices,
                                                radius=0.0)

    ### for now just do the dumb implementation of taking the mean curvature at the face center (need to implement 3d barycentric, but there shouldn't be too much variation so it'll be okay to test)                                      
    gaussColors = []
    gaussFace = il.propertybarycentric(verts[faces], gauss[faces])

    ### let's quicly generate a custom color sequence for the face colors     
    cvals = [np.min(gauss), np.mean(gauss), np.max(gauss)]
    colors = ['blue', 'white', 'red']
    norm = plt.Normalize(min(cvals), max(cvals))
    tuples = list(zip(map(norm, cvals), colors))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('', tuples)
    cmappable = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap('RdBu_r'))

    for i in range(len(faces)):
        gaussColors.append(cmap(gaussFace[i]))

    # instead plot the marching cubes result                             
    ax.add_collection3d(Poly3DCollection(verts[faces],
                                         alpha=0.5,
                                         edgecolors=['0.9'],
                                         facecolors=gaussColors,
                                         linewidth=0.0))
    
    figure.colorbar(mappable=cmappable, ax=ax,
                    fraction=0.025, pad=0.0, location='left')

    dimMin = np.min(verts[faces]); dimMax = np.max(verts[faces])
    ax.set_xlim([0.90*dimMin,
                 1.10*dimMax])
    ax.set_ylim([0.90*dimMin,
                 1.10*dimMax])
    ax.set_zlim([0.90*dimMin,
                 1.10*dimMax])
    plt.savefig('sasaSurf.png')
    plt.close()
    return
###

### Generate a 3D representation of the SASA for a given frame
def densityPlot(heavyPos, watPos, thisbox, level=0.016, figure=None, ax=None, 
                colors=['blue','white','red'],
                check=False):
    """ Given a set of protein heavy atom and water oxygen positions, generate an inst. interface mesh
    Input: heavyPos - (NPos, 3)-array of protein heavy atom positons
           watPos - (Nwat, 3)-array of water oxygen positions
           thisbox - (1, 3)-array of the box vectors
    Outputs: None
    """
    if figure is None:
        figure = plt.figure()
        ax = figure.add_subplot(111, projection='3d')

    ### Using marching_cubes, compute density faces and vertices    
    verts, faces = densityGrid(heavyPos[:4,:], watPos, thisbox, level=level)

    ### Using trimesh, compute the gaussian curvature at each vertex
    tmesh = tm.Trimesh(vertices=verts, faces=faces)
    gauss = discrete_gaussian_curvature_measure(mesh=tmesh,
                                                points=tmesh.vertices,
                                                radius=0.0)

    ### for now just do the dumb implementation of taking the mean curvature at the face center (need to implement 3d barycentric, but there shouldn't be too much variation so it'll be okay to test)                                      
    gaussColors = []
    gaussFace = il.propertybarycentric(verts[faces], gauss[faces])

    ### let's quicly generate a custom color sequence for the face colors
    cvals = [np.min(gauss), np.mean(gauss), np.max(gauss)]
    norm = plt.Normalize(min(cvals), max(cvals))
    tuples = list(zip(map(norm, cvals), colors))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('', tuples)
  
    for i in range(len(faces)):
        gaussColors.append(cmap(gaussFace[i]))

    # instead plot the marching cubes result                             
    ax.add_collection3d(Poly3DCollection(verts[faces],
                                         alpha=0.1,
                                         edgecolors=['0.9'],
                                         facecolors=gaussColors,
                                         linewidth=0.1))

    ax.set_xlim([-0.15*thisbox[0,0], 0.15*thisbox[0,0]])
    ax.set_ylim([-0.15*thisbox[0,0], 0.15*thisbox[0,0]])
    ax.set_zlim([-0.15*thisbox[0,0], 0.15*thisbox[0,0]])

#    ax.view_init(azim=30, elev=100)

    # check if last plot, then plot atoms    
    if check:
        # get sphere coords
        xs, ys, zs = genSphere()
        ax.plot_surface(heavyPos[0,0]+xs*0.5, 
                        heavyPos[0,1]+ys*0.5, 
                        heavyPos[0,2]+zs*0.5,
                        color='r')
        ax.plot_surface(heavyPos[1,0]+xs*0.5, 
                        heavyPos[1,1]+ys*0.5, 
                        heavyPos[1,2]+zs*0.5,
                        color='y')
        ax.plot_surface(heavyPos[2,0]+xs*0.5, 
                        heavyPos[2,1]+ys*0.5, 
                        heavyPos[2,2]+zs*0.5,
                        color='gray')
        ax.plot_surface(heavyPos[3,0]+xs*0.5, 
                        heavyPos[3,1]+ys*0.5, 
                        heavyPos[3,2]+zs*0.5,
                        color='gray')

        figure.savefig('densitySurf.png')
    else:
        figure.savefig('densitySurf.png')

    return figure, ax
###
