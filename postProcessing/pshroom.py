#!/usr/bin/env python3
#[ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ]#
#[
#[ Perform some post-processing operations on
#[ MuSHrooM data.
#[
#[ This python script uses a number of functions (utilities)
#[ defined in pshroomUtil.py.
#[
#[ Manaure Francisquez.
#[
#[ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ]#

from pylab import *
import argparse
import numpy as np
import pshroomUtil as psu
import adios as ad
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from shutil import copyfile
from scipy.optimize import curve_fit
import sys

#[ Handle command-line arguments
parser = argparse.ArgumentParser(description='Compare two MuSHrooM runs')
parser.add_argument('dataDir', default=None,
                    help='Location of MuSHrooM data')
parser.add_argument('outDir', default=None,
                    help='Directory where new files are placed')
parser.add_argument('-d0', '--diagnostic0', default=None, nargs='+', type=int,
                    help='First component of diagnostic to plot.')
parser.add_argument('-d1', '--diagnostic1', default=None, nargs='+', type=int,
                    help='Second component of diagnostic to plot.')
parser.add_argument('-d2', '--diagnostic2', default=None, nargs='+', type=int,
                    help='Third component of diagnostic to plot.')
parser.add_argument('-fr', '--frame-range', default=None, nargs=2, type=int,
                    help='Indicate the first and last frame to plot (comment out to plot all)')
parser.add_argument('--nocopy', default=False, action='store_true',
                    help="Don't copy the shroom files for record keeping")

args = parser.parse_args()

#[ Get dataDir and outDir from command line (if passed in).
#[ Otherwise use hard-coded values.
if args.dataDir and args.outDir:
  dataDir = args.dataDir
  outDir  = args.outDir
else:
  #dataDir = '/dartfs-hpc/scratch/mana/shroom/parallelTest2/'        #[ Location of MuSHrooM data.
  #outDir  = '/dartfs-hpc/scratch/mana/shroom/parallelTest2/post/'    #[ Directory where new files are placed.
  #dataDir = '/ihome/mana/documents/multiscale/code/HM-2D/src/production/'        #[ Location of MuSHrooM data.
  #outDir  = '/ihome/mana/documents/multiscale/code/HM-2D/postProcessing/post/'    #[ Directory where new files are placed.
  #dataDir = '/Users/mana/Documents/Research/MIT/code/HM-2D/src/production/'        #[ Location of MuSHrooM data.
  #outDir  = '/Users/mana/Documents/Research/MIT/code/HM-2D/src/production/post/'    #[ Directory where new files are placed.
  dataDir = '/nobackup1/francisquez//shroom/35943806/'        #[ Location of MuSHrooM data.
  outDir  = '/nobackup1/francisquez//shroom/35943806/post/'    #[ Directory where new files are placed.
  #dataDir = '/scratch/gpfs/manaurer//shroom/36725/'        #[ Location of MuSHrooM data.
  #outDir  = '/scratch/gpfs/manaurer//shroom/36725/post/'    #[ Directory where new files are placed.

adiosIO = True         #[ =True ADIOS, =False binary (input simulation parameters below).

outDataFile      = True    #[ If true save the post-processed data to a file.
outFigureFile    = True    #[ If false, show figure on screen (default=True).
figureFileFormat = 'png'   #[ Can be png, pdf, ps, eps, svg.

#[ Variables we wish to process.
varNames = ['phik']

#[ Choose diagnostic to plot.
#[   =[0,0] 2D plane (real).
#[   =[0,1] 2D plane (square amplitude of complex).
#[   =[1,.] fk^2 v. time for each ky (=[1,0]), kx (=[1,1]) or |k| (=[1,2]).
#[   =[2,.] fk^2-based growth rate v. ky (=[2,0]), kx (=[2,1]), or |k| ([=2,2]).
#[   =[3,.] fk^2 v. ky (=[3,0]), kx (=[3,1]), or |k| ([=3,2]), averaged over sampleTime.
#[   =[4,.] fk(ky=0)^2 as a function of x and time.
#[   =[5,.] particle flux, electrons (=[5,0]) or ions (=[5,1]).
#[   =[6,.] heat flux, electrons (=[6,0]) or ions (=[6,1]).
#[   =[7,.,.] particle flux spectrum, electrons (=[7,0,.]) or ions (=[7,1,.]),
#[            summed over x (=[7,.,0]) or y (=[7,.,1]) and averaged over time.
#[   =[8,.,.] heat flux spectrum, electrons (=[8,0,.]) or ions (=[8,1,.]),
#[            summed over x (=[8,.,0]) or y (=[8,.,1]) and averaged over time.
d2 = 0
d1 = 0
d0 = 0
if args.diagnostic0:
  d0 = args.diagnostic0[0]
if args.diagnostic1:
  d1 = args.diagnostic1[0]
if args.diagnostic2:
  d2 = args.diagnostic2[0]
diagnostic = [d0, d1, d2]

#[ Indicate whether temperatures were evolved (default = true).
#evolveTemp = False

#[ Indicate the first and last frame to plot (comment out to plot all):
#frameInitEnd = [1000, 4000]
if args.frame_range:
  frameInitEnd = args.frame_range

#[ Number of frames in previous run (to be added to the file name of figures from this run).
#prevRunFrames = 4000

#[ Perform operation (e.g. measure growth rate) in the 'sampleTime' time window.
#sampleTime = [500.0, 1000.0]

#[ For 2D planes, choose whether to pcolormesh or contour plot.
#[   =0 pcolor
#[   =1 contour
op2Dplot   = 0

#[ Plot limits (x, y, colorbar). Comment out or make empty for automatically adjusted limits.
#plotLimits = [list(),[0.,1.4],list()]

#[ For data that changes by orders of magnitude one can request log scale along
#[ x (isAxisLog[0]), y (isAxisLog[1]) or in the colorbar (isAxisLog[2]). The
#[ default=False. If isAxisLog[2]=True one must also indicate colorbar limits in plotLimits.
isAxisLog = [True, True, False]

#[ Some diagnostics require that the user indicate the y-axis or colorbar label.
cLabel = r'$[eL_n/(\rho_sT_e)]\phi(x,y)$'        #[ potential.
#cLabel = r'$[L_n/(\rho_sn_0)]n_e(x,y)$'        #[ Density.
#cLabel = r'$[L_n/(\rho_sT_{e0})]T_e(x,y)$'        #[ Temperature.
#cLabel = r'$(e\rho_sL_n/T_e)\nabla^2\phi(x,y)$'  #[ vorticity.

yLabel = r'$|\phi_{k}|^2$'
#yLabel = r'$2\pi k\left\langle|\phi_{k}|^2\right\rangle/(2\pi)^2$'

#cLabel = r'$[eL_n/(\rho_sT_e)]^2|\hat{\phi}|^2$'        #[ potential.
#cLabel = r'$|\hat{n}_e|^2$'        #[ potential.

#[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END OF USER INPUTS (maybe) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ]#

print("")
print("    <... MuSHrooM post-processing cabin ...>   ")
print("")

#[ Define various font sizes for plots.
xyLabelFontSize = 16
titleFontSize   = 16
cLabelFontSize  = 16
legendFontSize  = 14
xyTickFontSize  = 12

#[ Check if user requested colorbar limits. If not, adjust dynamically.
try:
  plotLimits
except NameError:
  adjustPlotLim = [True, True, True]
  plotLimits    = [list(),list(),list()]
else:
  adjustPlotLim = [True, True, True]
  for d in range(3):
    if len(plotLimits[d])==2:
      adjustPlotLim[d] = False

#[ By default output figure to file.
try:
  outFigureFile
except NameError:
  outFigureFile=True
#[ By default set the number of frames from previous runs to zero.
try:
  prevRunFrames
except NameError:
  prevRunFrames=0
#[ By default do not use log scale for x, y or color axis.
try:
  isAxisLog
except NameError:
  isAxisLog=[False, False, False]

try:
  evolveTemp
except NameError:
  evolveTemp = True
else:
  if any((e=="tempElck" or e=="tempIonk") for e in varNames):
    sys.exit("  evolveTemp=False so can't process temperature files.")

#[ Define IO functions (different for ADIOS and binary).
if adiosIO:
  psIO = psu.psIObp()
  fTypeExtension = '.bp'
else:
  psIO = psu.psIObin()
  fTypeExtension = '.bin'

if outFigureFile:
  #.To have a virtual display use the line below.
  plt.switch_backend('agg')

#[ Add an extra slash to the directories in case the user didn't.
dataDir = dataDir+'/'
outDir  = outDir+'/'

#[ Check if post-processing directory exists. Create it if not.
psu.checkMkdir(outDir)

#[ Copy source and input files into post-processing directory for record keeping.
if not args.nocopy:
  copyfile(dataDir+'shroom.F90',    outDir+'shroom.F90')
  copyfile(dataDir+'shroom.in',     outDir+'shroom.in')
  copyfile(dataDir+'shroomFLAGS.h', outDir+'shroomFLAGS.h')
  copyfile(dataDir+'makefile',      outDir+'makefile')
  copyfile(dataDir+'setup.sh',      outDir+'setup.sh')

nVars = len(varNames)    #[ Number of variables to post-process.

#[ Read variables saved as attributes.
fileAttributes = psIO.getAttributes(dataDir+varNames[0]+fTypeExtension)
#[ If the attributes of this file were messed up, one can load those from an earlie run.
#dataDirAlt = '/nobackup1/francisquez//shroom/30906245/'        #[ Location of MuSHrooM data.
#fileAttributes = psIO.getAttributes(dataDirAlt+varNames[0]+fTypeExtension)
Nkx          = fileAttributes['Nkx']            #[ Number of distinct (absolute magnitude) kx modes.
Nky          = fileAttributes['Nky']            #[ Number of distinct (absolute magnitude) kx modes.
kxMin        = fileAttributes['kxMin']          #[ Minimum finite amplitude of kx modes.
kyMin        = fileAttributes['kyMin']          #[ Minimum finite amplitude of ky modes.
omSte        = fileAttributes['omSte']          #[ Parameter in the curvature drift frequency.
omde         = fileAttributes['omde']           #[ Parameter in the diamagnetic drift frequency.
tau          = fileAttributes['tau']            #[ Ratio of ion to electron temperature.
muMass       = fileAttributes['muMass']         #[ Square root of ion to electron mass ratio.
deltae       = fileAttributes['deltae']         #[ Electron (T_perp+T_par)/T.
deltaPerpe   = fileAttributes['deltaPerpe']     #[ Electron T_perp/T.
eta_e        = fileAttributes['eta_e']          #[ Electron L_n/L_T.
deltai       = fileAttributes['deltai']         #[ Ion (T_perp+T_par)/T.
deltaPerpi   = fileAttributes['deltaPerpi']     #[ Ion T_perp/T.
eta_i        = fileAttributes['eta_i']          #[ Ion L_n/L_T.
lambdaD      = fileAttributes['lambdaD']        #[ Normalized Debye length.
adiabaticElc = fileAttributes['adiabaticElc']   #[ Adiabatic electrons? =0 no, =else yes.
adiabaticIon = fileAttributes['adiabaticIon']   #[ Adiabatic ions? =0 no, =else yes.
HDmodel      = fileAttributes['HDmodel']        #[ Hyperdiffusion model.

speciesName = ['Elc','Ion']  #[ Auxiliary list.
speciesSubs = ['e','i']      #[ Species subscript letter.

Nk   = [Nkx, Nky]
kMin = [kxMin, kyMin]

#[ Define k-space arrays.
Nekx, Neky, kx, ky, kxSq, kySq, kSq  = psu.kGrid(Nkx,kxMin,Nky,kyMin)

kMag = np.sqrt(kSq)    #[ Magnitude of wavenumber vector k.

Nek  = np.array([Nekx, Neky], dtype='int')          #[ Number of elements in k-space array.
Nx   = np.array([Nekx, (Neky-1)*2], dtype='int')    #[ Number of cells in de-aliased real space.

dualDim = [1, 0]

#[ After establishing the number of elements in each direction
#[ assign value to byte-size variables (used if employing binary IO).
psIO.setByteSizes(Nek)

#[ Determine the number of frames saved in file.
sFrames = psIO.getFrameN(dataDir+varNames[0]+fTypeExtension)

#[ If user did not input total desired number of frames to plot, set to plot all frames.
try:
  frameInitEnd
except NameError:
  frameInitEnd = [0, sFrames-1]

#[ Total number of frames to plot.
pFrames = frameInitEnd[1]-frameInitEnd[0]+1

if sFrames < pFrames:
  sys.exit("  Error: cannot plot more frames than file contains. Terminating... ")


print("  Processing Shroom simulation with:")
print("    Nkx = ",Nkx,"  |  Nky = ",Nky)
print("  Located in ",dataDir)
print("")
print("  Will output diagnostic ", diagnostic, " for variable : ", varNames)
print("  in ",outDir)
print("")

time = 0.0

if diagnostic[0] == 0:

  if diagnostic[1] == 0:
    #[ Will plot in real space.
    Lx = np.array([2.0*np.pi/kMin[0], 2.0*np.pi/kMin[1]],dtype='float')    #[ Simulation box size.
    dx = np.array([Lx[0]/float(Nx[0]+np.mod(Nx[0],2)), Lx[1]/float(Nx[1]+np.mod(Nx[1],2))],dtype='float')    #[ Cell size.
    x  = [np.arange(Nx[0],dtype='float')*dx[0]-Lx[0]/2.0+(1.0-np.mod(Nx[0],2))*0.5*dx[0], \
          np.arange(Nx[1],dtype='float')*dx[1]-Lx[1]/2.0+(1.0-np.mod(Nx[1],2))*0.5*dx[1]]
    xLabel, yLabel = r'$x/\rho_s$', r'$y/\rho_s$'
  elif diagnostic[1] == 1:
    #[ Will plot in k-space.
    Nx = np.array([Nekx, Neky], dtype='int')
    dx = np.array([kxMin, kyMin], dtype='float')    #[ Cell size.
    x  = [kx, ky]
    #[ The above provide the grid as used in shroom. But for plotting
    #[ we want to shift the grid so kx=0 is in the center.
    x = [np.append(kx[Nkx:],kx[:Nkx]), ky]
    xLabel, yLabel = r'$k_x\rho_s$', r'$k_y\rho_s$'

  pData = np.zeros((Nx[0],Nx[1]))    #[ Array with data to be plotted.

  #[ Create enlarged array with zero value for last ky.
  inDataC = np.empty((Nekx,Neky), dtype='complex')

  #[ Mesh grid for plotting. It needs one more point along each direction
  #[ per matplotlib standards, i.e. need to provide nodal grid points, while
  #[ the data is cell-centered.
  xNodal = [np.append(np.insert(0.5*(x[0][:-1]+x[0][1:]),0,x[0][0]-dx[0]/2.),x[0][-1]+dx[0]/2.), \
            np.append(np.insert(0.5*(x[1][:-1]+x[1][1:]),0,x[1][0]-dx[1]/2.),x[1][-1]+dx[1]/2.)]
  X = [np.outer(xNodal[0],np.ones(xNodal[1].shape)), np.outer(np.ones(xNodal[0].shape),xNodal[1])]

  #[ Prepare figure.
  figProp1a = [5.6,4.8]
  ax1aPos   = [0.12, 0.135, 0.66, 0.77]
  cax1aPos  = [0.8, 0.135, 0.02, 0.77]
  fig1      = plt.figure(figsize=(figProp1a[0],figProp1a[1]))
  ax1a      = fig1.add_axes(ax1aPos)
  ax1a.set_xlabel(xLabel, fontsize=xyLabelFontSize)
  ax1a.set_ylabel(yLabel, fontsize=xyLabelFontSize, labelpad=-1)
  #[ Plot something to get the colorbar initialized.
  if op2Dplot == 0:
    cbar_ax1a = fig1.add_axes(cax1aPos)
    if not isAxisLog[2]:
      hpl1a = ax1a.pcolormesh(X[0], X[1], pData)
      cbara = plt.colorbar(hpl1a,ax=ax1a,cax=cbar_ax1a)
    else:
      hpl1a = ax1a.pcolormesh(X[0], X[1], pData,
                              norm=colors.SymLogNorm(linthresh=0.03,linscale=0.03,
                              vmin=plotLimits[2][0],vmax=plotLimits[2][1]),cmap='RdBu_r')
      cbara = plt.colorbar(hpl1a,ax=ax1a,cax=cbar_ax1a,extend='both')
    cbara.set_label(cLabel, rotation=270, labelpad=14, fontsize=cLabelFontSize)
    #[ If user specifies colorbar limits do not update colorbar.
    if not adjustPlotLim[2]:
      hpl1a.set_clim(plotLimits[2][0],plotLimits[2][1])
      cbara.set_clim(plotLimits[2][0],plotLimits[2][1])
  elif op2Dplot == 1:
    hpl1a = ax1a.contour(X[0], X[1], pData, colors='black')
  ax1a.set_title(r'Time, $c_s t/L_n = $'+('{:06.4f}'.format(float(time))))
  #[ Set the x-y limits.
  psu.setXYlim(ax1a, adjustPlotLim, [[xNodal[0][0],xNodal[0][-1]],[xNodal[1][0],xNodal[1][-1]]], plotLimits)

elif (diagnostic[0] > 0):

  #[ Extract the array of time stamps for desired frames.
  pTimes = psIO.getTimes(dataDir+varNames[0]+fTypeExtension,frameInitEnd[0],frameInitEnd[1])

  if (diagnostic[0] >= 5) and (diagnostic[0] <= 8):
    if diagnostic[0] <= 6:
      flux = np.zeros(pFrames)    #[ Array with particle/heat flux averaged over space in time.
    else:
      flux = np.zeros(Nk[dualDim[diagnostic[2]]])    #[ Array with particle/heat flux spectrum.
    #.Get additional arrays needed to compute fluxes.
    iky       = 1.0j*ky
    varNames  = ['phik']
    if evolveTemp:
      inDataAdd = [0.,0.]
    else:
      inDataAdd = [0.]
    if (diagnostic[0] == 5) or (diagnostic[0] == 7):  #[ Particle flux.
      flrFuncs = psu.calcFLR(kx,ky,diagnostic[1],tau,muMass,only=["Sb","Db","hatLap","avgJ0","Gamma0"])
      phiFac   = 1.0-flrFuncs["Gamma0"]
      #.These are the functions in the Poisson equation.
      #denFac   = flrFuncs["Sb"]
      #tempFac  = (1.0/flrFuncs["Db"])*0.5*flrFuncs["hatLap"]*flrFuncs["avgJ0"]
      #.These are the functions in the evolution equation.
      denFac   = flrFuncs["avgJ0"]
      tempFac  = 0.5*flrFuncs["hatLap"]*flrFuncs["avgJ0"]
    elif (diagnostic[0] == 6) or (diagnostic[0] == 8): #[ Heat flux.
      flrFuncs = psu.calcFLR(kx,ky,diagnostic[1],tau,muMass,only=["hatLap","avgJ0","Gamma0"])
      phiFac   = 0.0
      if diagnostic[1] == 0:  #[ Electrons.
        denFac  = (3./2.+deltaPerpe*0.5*flrFuncs["hatLap"])*flrFuncs["avgJ0"]
        tempFac = (1.0+0.5*flrFuncs["hatLap"])*flrFuncs["avgJ0"]
      elif diagnostic[1] == 1:  #[ Ions.
        denFac  = (3./2.+deltaPerpi*0.5*flrFuncs["hatLap"])*flrFuncs["avgJ0"]
        tempFac = (1./tau)*(1.0+0.5*flrFuncs["hatLap"])*flrFuncs["avgJ0"]

  elif diagnostic[0] == 4:

    #[ Will plot in real space along x.
    Lx = 2.0*np.pi/kMin[0]    #[ Simulation box size along x.
    dx = Lx/float(Nx[0]+np.mod(Nx[0],2))    #[ Cell size along x.
    x  = np.arange(Nx[0])*dx-Lx/2.0+(1.0-np.mod(Nx[0],2))*0.5*dx

    pData = np.zeros((pTimes.size,Nx[0]))    #[ Array with data to be plotted.

    #[ Create enlarged array with zero value for last ky.
    inDataC = np.empty((Nekx,1), dtype='complex')

    #[ Mesh grid for plotting.
    X = [np.outer(pTimes,np.ones(x.shape)), np.outer(np.ones(pTimes.shape),x)]

    #[ Prepare figure.
    figProp1a = [7,3.2]
    ax1aPos   = [0.125, 0.17, 0.73, 0.78]
    cax1aPos  = [0.88, 0.17, 0.03, 0.78]
    fig1      = plt.figure(figsize=(figProp1a[0],figProp1a[1]))
    ax1a      = fig1.add_axes(ax1aPos)

  elif (diagnostic[0] > 0):
    if (diagnostic[0] == 2) or (diagnostic[0] == 3):
      #[ If user does not indicate sampleTime, set it to the whole period.
      try:
        sampleTime
      except NameError:
        sampleTime = [pTimes[0],pTimes[-1]]
      #[ Terminate if sampleTime is outside of plotting range.
      if (sampleTime[0] > pTimes[-1]) or (sampleTime[1] > pTimes[-1]):
        sys.exit("  Error: sampleTime outside of plotting window. Terminating... ")

    #[ Array of square amplitudes as a function of kx or ky and time.
    if diagnostic[1] == 0:
      dNk  = Neky
      kVar = ky
    elif (diagnostic[1] == 1) or (diagnostic[1] == 2):
      dNk  = Nekx
      #[ For plotting we want to shift the grid so kx=0 is in the center.
      kVar = np.append(kx[Nkx:],kx[:Nkx])

    fkp  = np.zeros((Nekx,Neky))     #[ New Fourier space quantity derived from fk.
    fkr  = np.zeros((dNk,pFrames))   #[ Reduced fkr in time.

    if diagnostic[1] == 2:
      #[ Count the number of k's used in averaged over bands.
      kBandCount = np.zeros(dNk)
      for kk in range(dNk):
        for jk in range(Neky):
          for ik in range(Nekx):
            if ((kMag[jk,ik] > ky[kk]-0.5*kyMin) and ((kMag[jk,ik] < ky[kk]+0.5*kyMin))):
              kBandCount[kk] = kBandCount[kk]+1.0

  #[ Prepare figure.
  if (diagnostic[0] == 2) or (diagnostic[0] == 3) or ((diagnostic[0] >= 5) and (diagnostic[0] <= 8)):
    figProp1a = [7,3.2]
    ax1aPos   = [0.16, 0.17, 0.8, 0.8]
  elif (diagnostic[0] == 1):
    #[ Need more room for a k-spectrum colorbar.
    figProp1a = [7,3.2]
    ax1aPos   = [0.125, 0.17, 0.75, 0.8]
    cax1aPos  = [0.89, 0.17, 0.03, 0.8]
  fig1      = plt.figure(figsize=(figProp1a[0],figProp1a[1]))
  ax1a      = fig1.add_axes(ax1aPos)
  if (diagnostic[0] == 1):
    #[ Add colorbar for k spectrum.
    ax1a.set_prop_cycle(plt.cycler('color', plt.cm.viridis(np.linspace(0, 1, dNk))))
    cbar_ax1a = fig1.add_axes(cax1aPos)

#[ Set the tick font size.
for tick in ax1a.xaxis.get_major_ticks():
  tick.label.set_fontsize(xyTickFontSize)
for tick in ax1a.yaxis.get_major_ticks():
  tick.label.set_fontsize(xyTickFontSize)

#[ Loop over variables, if multiple.
for iV in varNames:

  #[ Open data file.
  fH = psIO.fOpen(dataDir+iV+fTypeExtension)
  #[ Read attributes (moves the position within the file for binary IO).
  psIO.skipAttributes(fH)

  if diagnostic[0] == 0:
    if diagnostic[1] == 0:
      outDir = outDir+iV.split('k')[0]+'/'    #[ Create a folder with variable name.
      psu.checkMkdir(outDir)
      outFigFileNameFmt = outDir+iV.split('k')[0]+"_xy-f%05d."+figureFileFormat
    elif diagnostic[1] == 1:
      outDir = outDir+iV+'Sq/'    #[ Create a folder with variable name.
      psu.checkMkdir(outDir)
      outFigFileNameFmt = outDir+iV+"Sq_kxky-f%05d."+figureFileFormat
  elif ((diagnostic[0] >= 5) and (diagnostic[0] <= 8)):
    #[ Will also need the particle density and temperature.
    fHadd = [ psIO.fOpen(dataDir+'den'+speciesName[diagnostic[1]]+'k'+fTypeExtension) ]
    psIO.skipAttributes(fHadd[0])
    if evolveTemp:
      fHadd.append( psIO.fOpen(dataDir+'temp'+speciesName[diagnostic[1]]+'k'+fTypeExtension) )
      psIO.skipAttributes(fHadd[1])

  fCntr = 0   #[ Frame counter: number of frames processed.
  for iF in range(frameInitEnd[0],frameInitEnd[1]+1):

    #[ Read current frame data and time stamp. readFrame can also
    #[ return other time-varying quantities. See pshroomUtil.py.
    time, inData, _, _, _, _, _, _, _ = psIO.readFrame(fH,iF)

    if diagnostic[0] == 0:
      if diagnostic[1] == 0:
        #[ Transform to real space and store data.
        inDataC = psIO.reorderComplexField(inData,Nk,Nek)    #[ Rearrange into order FFT expects.
        pData   = np.fft.irfftn(inDataC, axes=(0,1))  #[ Assumes nonunitary normalization, w/ 1/N factor in FFT_r2c.
      elif diagnostic[1] == 1:
        ##[ Plot square amplitude of complex data.
        #for ik in range(Nekx):
        #  for jk in range(Neky):
        #    pData[ik,jk] = psIO.getSqNormIJ(inData,ik,jk)

        #.Need to reorganize the data so the k's are in ascending order.
        for jk in range(Nx[1]):
          pData[:Nkx-1,jk] = inData[Nkx:,jk].real**2+inData[Nkx:,jk].imag**2
          pData[Nkx-1:,jk] = inData[:Nkx,jk].real**2+inData[:Nkx,jk].imag**2

      if op2Dplot == 0:
        hpl1a.set_array(pData.ravel())
        if adjustPlotLim[2]:
          hpl1a.autoscale()
        else:
          cbara.set_clim(plotLimits[2][0],plotLimits[2][1])
      elif op2Dplot == 1:
        plt.cla()
        hpl1a = ax1a.contour(X[0], X[1], pData, colors='black')
      ax1a.set_title(r'Time, $c_s t/L_n = $'+('{:06.4f}'.format(float(time))))
      plt.savefig(outFigFileNameFmt % (prevRunFrames+iF), dpi=200, format=figureFileFormat)

    elif ((diagnostic[0] >= 5) and (diagnostic[0] <= 8)):
      #[ Will plot particle flux in time. Read density and temperature.
      _, inDataAdd[0], _, _, _, _, _, _, _ = psIO.readFrame(fHadd[0],iF)
      if evolveTemp:
        _, inDataAdd[1], _, _, _, _, _, _, _ = psIO.readFrame(fHadd[1],iF)
      #[ Compute the velocity in real space.
      vExk = inData
      for jk in range(Neky):
        vExk[:,jk] = -iky[jk]*vExk[:,jk]
      #[ Compute the particle density in particle coordinates.
      if evolveTemp:
        denPk = denFac*inDataAdd[0]+tempFac*inDataAdd[1]+phiFac*inData
      else:
        denPk = denFac*inDataAdd[0]+phiFac*inData
      fluxDensityk = vExk*np.conjugate(denPk)

      if diagnostic[0] <= 6:
        #[ Volume average.
        for j in range(Nek[1]):
          for i in range(Nek[0]):
            flux[fCntr] = flux[fCntr]+np.real(fluxDensityk[i,j])
      else:
        #[ Sum over one dimension.
        flux = flux+np.abs(np.sum(fluxDensityk, axis=diagnostic[2])[:Nk[dualDim[diagnostic[2]]]])

    elif (diagnostic[0] == 4):

      #[ Transform to real space and store data.
      inDataC = psIO.reorderComplexField(inData[:,0],[Nk[0],1],[Nek[0],1])    #[ Rearrange into order FFT expects.
      pData[fCntr,:] = np.real(np.fft.ifftn(inDataC))

    elif (diagnostic[0] > 0):

      if (diagnostic[0] > 0):
        #[ Obtain square amplitude vs. time for each k.
        fkp = psIO.getSqNorm(inData)

      if diagnostic[1] == 0:
        #[ Average over kx to plot vs. ky.
        for jk in range(dNk):
          fkr[jk,fCntr] = np.mean(fkp[:,jk])
          #fkr[jk,fCntr] = fkp[0,jk]  #[ Use this to select kx=0 only. Change output filenames below accordingly.
      elif diagnostic[1] == 1:
        #[ Reorganize so negative kx's are first. Average over ky to plot vs. kx.
        fkr[:Nkx-1,fCntr] = np.mean(fkp[Nkx:,:],axis=1)
        fkr[Nkx-1:,fCntr] = np.mean(fkp[:Nkx,:],axis=1)
      elif diagnostic[1] == 2:
        #[ Average quantity over a band of k's.
        for kk in range(dNk):
          for jk in range(Neky):
            for ik in range(Nekx):
              if ((kMag[jk,ik] > ky[kk]-0.5*kyMin) and ((kMag[jk,ik] < ky[kk]+0.5*kyMin))):
                fkr[kk,fCntr] = fkr[kk,fCntr]+fkp[ik,jk]
          fkr[kk,fCntr] = (2.0*np.pi*ky[kk]*(fkr[kk,fCntr]/kBandCount[kk])/((2.0*np.pi)**2))*float(np.prod(Nx))  #[ Python defaults to norm="backward" in irfft, and shroom uses nonunitary FFT.

    fCntr = fCntr + 1    #[ Increment frame counter.

  #[ END OF LOOP OVER TIME FRAMES........ ]#

  fH.close()

  if ((diagnostic[0] >= 5) and (diagnostic[0] <= 8)):
    #[ Close density and temperature files.
    fHadd[0].close()
    if evolveTemp:
      fHadd[1].close()

    if (diagnostic[0] == 5) or (diagnostic[0] == 7):
      fluxLabel, fluxLetter = 'n', '\Gamma'
    elif (diagnostic[0] == 6) or (diagnostic[0] == 8):
      fluxLabel, fluxLetter = 'heat', 'Q'

    if diagnostic[0] <= 6:  #[ Volume average flux in time.
      print(" Mean flux in last 1/4 of simulation = ",np.mean(flux[3*(np.size(flux)//4):]))
      print("")

      if outDataFile:  #[ Save the data to a file.
        psIO.dataToFile(fileAttributes,{'time':pTimes,'flux':flux},outDir+fluxLabel+'Flux'+speciesName[diagnostic[1]]+'-xyAv.bp')

      hpl1a = ax1a.plot(pTimes, flux)  #[ Plot flux vs. time.

      xLabel = r'$c_s t/L_n$'
      if (diagnostic[0] == 5):
        yLabel = r'$\Gamma_{%s} \thinspace (c_sn_0\rho_s^2/L_n^2)$' % speciesSubs[diagnostic[1]]
      elif (diagnostic[0] == 6):
        yLabel = r'$Q_{%s}\thinspace (c_sn_0T_{i0}\rho_s^2/L_n^2)$' % speciesSubs[diagnostic[1]]

      ##[ Load data from Smith's 1997 PoP figure 2 (particle flux).
      #smithFluxData = psu.getSmithFlux()
      #hpl1b = ax1a.plot(smithFluxData[:,0], smithFluxData[:,1], linestyle='--')
      #ax1a.legend([r'MuSHrooM',r'SH 1997'], fontsize=legendFontSize, frameon=False)

      #[ Set the x-y limits.
      psu.setXYlim(ax1a, adjustPlotLim, [[pTimes[0],pTimes[-1]],[None,None]], plotLimits, useDefaultY=False)

      figLabel = '-xyAv'
    else:    #[ Flux spectrum.
      flux = flux/fCntr  #[ Complete the time average.

      if diagnostic[2] == 0:
        kVar = ky
        xLabel, yLabel = r'$k_y\rho_s$', r'$\langle\sum_{k_x}\,\hat{%s}_{%s}\rangle_t$' % (fluxLetter,speciesSubs[diagnostic[1]])
      elif diagnostic[2] == 1:
        kVar = kx[:Nkx]
        xLabel, yLabel = r'$k_x\rho_s$', r'$\langle\sum_{k_y}\,\hat{%s}_{%s}\rangle_t$' % (fluxLetter,speciesSubs[diagnostic[1]])

      if outDataFile:  #[ Save the data to a file.
        if diagnostic[2] == 0:
          figLabel = '-kxSumTimeAv-f'+str(frameInitEnd[0])+'-'+str(frameInitEnd[1])
          psIO.dataToFile(fileAttributes,{'ky':kVar,'flux':flux},outDir+fluxLabel+'Flux'+speciesName[diagnostic[1]]+figLabel+'.bp')
        elif diagnostic[2] == 1:
          figLabel = '-kySumTimeAv-f'+str(frameInitEnd[0])+'-'+str(frameInitEnd[1])
          psIO.dataToFile(fileAttributes,{'kx':kVar,'flux':flux},outDir+fluxLabel+'Flux'+speciesName[diagnostic[1]]+figLabel+'.bp')

     #[ Plot flux spectrum.
      if isAxisLog[0] and isAxisLog[1]:
        hpl1a = ax1a.loglog(kVar, flux)
      elif isAxisLog[0]:
        hpl1a = ax1a.semilogx(kVar, flux)
      elif isAxisLog[1]:
        hpl1a = ax1a.semilogy(kVar, flux)
      else:
        hpl1a = ax1a.plot(kVar, flux)

    ax1a.set_xlabel(xLabel, fontsize=xyLabelFontSize, labelpad=-2)
    ax1a.set_ylabel(yLabel, fontsize=xyLabelFontSize, labelpad=-1)

    if outFigureFile:
      plt.savefig(outDir+fluxLabel+'Flux'+speciesName[diagnostic[1]]+figLabel+'.'+figureFileFormat, format=figureFileFormat)
    else:
      plt.show()
  elif (diagnostic[0] == 1):
    #[ Create scalar mappable for colorbar indicating wavenumber k.
    sm    = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=np.amin(kVar), vmax=np.amax(kVar)))
    sm._A = []
    cbar  = plt.colorbar(sm,ax=ax1a,cax=cbar_ax1a)

    #[ Plot square amplitude vs. k.
    for ik in range(1,dNk):    #[ If starting at ik=1, don't plot k=0 mode.
      hpl1a = ax1a.semilogy(pTimes, fkr[ik,:])

    #[ Set the x-y limits.
    psu.setXYlim(ax1a, adjustPlotLim, [[pTimes[0],pTimes[-1]],[None,None]], plotLimits, useDefaultY=False)

    xLabel = r'$c_s t/L_n$'
    if diagnostic[1] == 0:
      cLabel, fNameSuffix = r'$k_y\rho_s$', '-kxAv.'
    elif diagnostic[1] == 1:
      cLabel, fNameSuffix = r'$k_x\rho_s$', '-kyAv.'
    elif diagnostic[1] == 2:
      cLabel, fNameSuffix = r'$|k|\rho_s$', '-kBandAv.'

    ax1a.set_xlabel(xLabel, fontsize=xyLabelFontSize, labelpad=-2)
    ax1a.set_ylabel(yLabel, fontsize=xyLabelFontSize, labelpad=-1)
    cbar.set_label(cLabel, rotation=270, labelpad=+20, fontsize=cLabelFontSize)

    if outFigureFile:
      if (diagnostic[0] == 1):
        fNamePrefix = iV+'Sq'

      plt.savefig(outDir+fNamePrefix+fNameSuffix+figureFileFormat, format=figureFileFormat)
    else:
      plt.show()
  elif (diagnostic[0] == 2):
    #[ Measure the growth rate for each k.
    gamma = np.zeros(dNk)
    tInd  = [psu.findNearestIndex(pTimes,sampleTime[0]), psu.findNearestIndex(pTimes,sampleTime[1])]
    for ik in range(1,dNk):  #[ If starting at ik=1, don't compute growth rate for k=0 mode.
      lineFitCoeffs, _ = curve_fit(psu.lineEq, pTimes[tInd[0]:tInd[1]], np.log10(fkr[ik,tInd[0]:tInd[1]]))
      gamma[ik]        = 0.5*np.log(10.0)*lineFitCoeffs[0]

    if outDataFile:
      #[ Save the data to a file.
      if (diagnostic[0] > 0):
        fNamePrefix = iV+'Sq'

      yLabel = r'Growth rate, $\gamma L_n/c_s$'
      if diagnostic[1] == 0:
        psIO.dataToFile(fileAttributes,{'ky':kVar,'gamma':gamma},outDir+'growthRate-'+fNamePrefix+'-kyAv.bp')
        xLabel = r'$k_y\rho_s$'
      elif diagnostic[1] == 1:
        psIO.dataToFile(fileAttributes,{'kx':kVar,'gamma':gamma},outDir+'growthRate-'+fNamePrefix+'-kxAv.bp')
        #psIO.dataToFile(fileAttributes,{'ky':kVar,'gamma':gamma},outDir+'growthRate-'+fNamePrefix+'-kx0p0.bp')
        xLabel = r'$k_x\rho_s$'
      elif diagnostic[1] == 2:
        psIO.dataToFile(fileAttributes,{'k':kVar,'gamma':gamma},outDir+'growthRate-'+fNamePrefix+'-kBandAv.bp')
        xLabel = r'$|k|\rho_s$'

    hpl1a = ax1a.semilogx(kVar, gamma)
    ##[ Load data from Smith's 1997 PoP figure 3 (growth rate).
    #smithGrowthRate = psu.getSmithMaxGrowthRate()
    #hpl1b = ax1a.semilogx(smithGrowthRate[:,0], smithGrowthRate[:,1], linestyle='--')
    #ax1a.legend([r'MuSHrooM',r'SH 1997'], fontsize=legendFontSize, frameon=False)

    #[ Set the x-y limits.
    psu.setXYlim(ax1a, adjustPlotLim, [[None,None],[None,None]], plotLimits, useDefaultX=False, useDefaultY=False)

    ax1a.set_xlabel(xLabel, fontsize=xyLabelFontSize, labelpad=-2)
    ax1a.set_ylabel(yLabel, fontsize=xyLabelFontSize, labelpad=-1)

    if outFigureFile:
      if (diagnostic[0] > 0):
        fNamePrefix = iV+'Sq'

      if diagnostic[1] == 0:
        plt.savefig(outDir+'growthRate-'+fNamePrefix+'-kxAv.'+figureFileFormat, format=figureFileFormat)
#        plt.savefig(outDir+'growthRate-'+fNamePrefix+'-kx0p0.'+figureFileFormat, format=figureFileFormat)
      elif diagnostic[1] == 1:
        plt.savefig(outDir+'growthRate-'+fNamePrefix+'-kyAv.'+figureFileFormat, format=figureFileFormat)
      elif diagnostic[1] == 2:
        plt.savefig(outDir+'growthRate-'+fNamePrefix+'-kBandAv.'+figureFileFormat, format=figureFileFormat)
    else:
      plt.show()
  elif (diagnostic[0] == 3):
    #[ Plot amplitude vs. wavenumber.
    tInd  = [psu.findNearestIndex(pTimes,sampleTime[0]), psu.findNearestIndex(pTimes,sampleTime[1])]
    hpl1a = ax1a.loglog(kVar, np.mean(fkr[:,tInd[0]:tInd[1]],1))

    #[ Set the x-y limits.
    psu.setXYlim(ax1a, adjustPlotLim, [[None,None],[None,None]], plotLimits, useDefaultX=False, useDefaultY=False)

    ##[ Load data from Smith's 1997 PoP figure 3 (Omega(k) spectrum).
    #smithOmegaSpectrum = psu.getSmithOmegaSpectrum()
    #hpl1b = ax1a.loglog(smithOmegaSpectrum[:,0], smithOmegaSpectrum[:,1], linestyle='--')
    #ax1a.legend([r'MuSHrooM',r'SH 1997'], fontsize=legendFontSize, frameon=False)

    if diagnostic[1] == 0:
      xLabel, fNameSuffix = r'$k_y\rho_s$', '-kxAv.'
    elif diagnostic[1] == 1:
      xLabel, fNameSuffix = r'$k_x\rho_s$', '-kyAv.'
    elif diagnostic[1] == 2:
      xLabel, fNameSuffix = r'$|k|\rho_s$', '-kBandAv.'

    ax1a.set_xlabel(xLabel, fontsize=xyLabelFontSize, labelpad=-2)
    ax1a.set_ylabel(yLabel, fontsize=xyLabelFontSize, labelpad=-1)

    if outFigureFile:
      if (diagnostic[0] > 0):
        fNamePrefix = iV+'Sq'

      plt.savefig(outDir+'kSpectrum-'+fNamePrefix+fNameSuffix+figureFileFormat, format=figureFileFormat)
    else:
      plt.show()
  elif (diagnostic[0] == 4):
    #.Plot amplitude of zonal component in x versus time.
    xLabel = r'$c_s t/L_n$'
    yLabel = r'$x/\rho_s$'
    cLabel = r'$|\phi_{k_y=0}|^2$'        #[ zonal potential.

    #[ Plot something to get the colorbar initialized.
    cbar_ax1a = fig1.add_axes(cax1aPos)
    hpl1a = ax1a.pcolormesh(X[0], X[1], pData)
    cbara = plt.colorbar(hpl1a,ax=ax1a,cax=cbar_ax1a)
    cbara.set_label(cLabel, rotation=270, labelpad=14, fontsize=cLabelFontSize)
    ax1a.set_xlabel(xLabel, fontsize=xyLabelFontSize, labelpad=-2)
    ax1a.set_ylabel(yLabel, fontsize=xyLabelFontSize, labelpad=-1)
    #[ If user specifies colorbar limits do not update colorbar.
    if not adjustPlotLim[2]:
      hpl1a.set_clim(plotLimits[2][0],plotLimits[2][1])
      cbara.set_clim(plotLimits[2][0],plotLimits[2][1])
    #[ Set the x-y limits.
    psu.setXYlim(ax1a, adjustPlotLim, [[pTimes[0],pTimes[-1]],[x[0],x[-1]]], plotLimits)

    if outFigureFile:
      plt.savefig(outDir+iV+'ZonalvTime.'+figureFileFormat, format=figureFileFormat)
    else:
      plt.show()

print("    >... post MuSHrooM completed ...<   ")
print("")
