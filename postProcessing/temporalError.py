#!/usr/bin/env python3
#[ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ]#
#[
#[ Compare results from two MuSHrooM runs.  This script is largely based off
#[ of the file pshroom.py, written by Manaure Francisquez.
#[
#[ Dan Reynolds @ SMU
#[ Cody Balos @ LLNL
#[
#[ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ]#

import matplotlib
matplotlib.use('Agg') #[ need this when a display is not available (e.g. ssh w/o X)
from pylab import *
import argparse
import getopt, sys, os.path
import numpy as np
import pshroomUtil as psu
import adios as ad
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from filters import *


class ShroomVariable:

  def __init__(self, dataDir, varName):

    #[ Read variables saved as attributes.
    fileAttributes = psIO.getAttributes(dataDir+varName+fTypeExtension)
    self.Nkx          = fileAttributes['Nkx']            #[ Number of distinct (absolute magnitude) kx modes.
    self.Nky          = fileAttributes['Nky']            #[ Number of distinct (absolute magnitude) kx modes.
    self.kxMin        = fileAttributes['kxMin']          #[ Minimum finite amplitude of kx modes.
    self.kyMin        = fileAttributes['kyMin']          #[ Minimum finite amplitude of ky modes.
    self.omSte        = fileAttributes['omSte']          #[ Parameter in the curvature drift frequency.
    self.omde         = fileAttributes['omde']           #[ Parameter in the diamagnetic drift frequency.
    self.tau          = fileAttributes['tau']            #[ Ratio of ion to electron temperature.
    self.mu           = fileAttributes['muMass']         #[ Square root of ion to electron mass ratio.
    self.deltae       = fileAttributes['deltae']         #[ Electron (T_perp+T_par)/T.
    self.deltaPerpe   = fileAttributes['deltaPerpe']     #[ Electron T_perp/T.
    self.eta_e        = fileAttributes['eta_e']          #[ Electron L_n/L_T.
    self.deltai       = fileAttributes['deltai']         #[ Ion (T_perp+T_par)/T.
    self.deltaPerpi   = fileAttributes['deltaPerpi']     #[ Ion T_perp/T.
    self.eta_i        = fileAttributes['eta_i']          #[ Ion L_n/L_T.
    self.lambdaD      = fileAttributes['lambdaD']        #[ Normalized Debye length.
    self.adiabaticElc = fileAttributes['adiabaticElc']   #[ Adiabatic electrons? =0 no, =else yes.
    self.adiabaticIon = fileAttributes['adiabaticIon']   #[ Adiabatic ions? =0 no, =else yes.
    self.HDmodel      = fileAttributes['HDmodel']        #[ Hyperdiffusion model.

    self.Nk   = [self.Nkx,self.Nky]
    self.kMin = [self.kxMin, self.kyMin]


def PrepareFigure(xyLabelFontSize = 16, cLabelFontSize = 16, xyTickFontSize = 12, legendFontSize = 14,
                  xLabel = r'$x/\rho_s$', yLabel = r'$y/\rho_s$'):
  #[ Prepare figure.
  figProp1a = [5.6,4.8]
  ax1aPos   = [0.12, 0.135, 0.66, 0.77]
  cax1aPos  = [0.8, 0.135, 0.02, 0.77]
  fig1      = plt.figure(figsize=(figProp1a[0],figProp1a[1]))
  ax1a      = fig1.add_axes(ax1aPos)
  ax1a.set_xlabel(xLabel, fontsize=xyLabelFontSize)
  ax1a.set_ylabel(yLabel, fontsize=xyLabelFontSize, labelpad=-1)
  for tick in ax1a.xaxis.get_major_ticks():
    tick.label.set_fontsize(xyTickFontSize)
  for tick in ax1a.yaxis.get_major_ticks():
    tick.label.set_fontsize(xyTickFontSize)
  return fig1, ax1a, cax1aPos


def RealFramePlot(x, Nx, Nek, cField, time, frameNo, outputDir = './',
             figureFileFormat = 'jpeg', showFigureFile = False, titleFontSize = 16,
             xyLabelFontSize = 16, cLabelFontSize = 16, xyTickFontSize = 12, legendFontSize = 14,
             xLabel = r'$x/\rho_s$', yLabel = r'$y/\rho_s$', cLabel = r'$[eL_n/(\rho_sT_e)]\phi(x,y)$'):

  fig1, ax1a, cax1aPos = PrepareFigure(xyLabelFontSize, cLabelFontSize,
                                       xyTickFontSize, legendFontSize,
                                       xLabel, yLabel)

  #[ Array that will hold data to be plotted.
  pData = np.zeros((Nx[0],Nx[1]))

  #[ Mesh grid for plotting.
  X = np.array([np.outer(x[0],np.ones(x[1].shape)), np.outer(np.ones(x[0].shape),x[1])])

  #[ Initialize the plot.
  cbar_ax1a = fig1.add_axes(cax1aPos)
  hpl1a = ax1a.pcolormesh(X[0], X[1], pData)
  cbara = plt.colorbar(hpl1a,ax=ax1a,cax=cbar_ax1a)
  cbara.set_label(cLabel, rotation=270, labelpad=14, fontsize=cLabelFontSize)

  ax1a.set_title(r'Time, $c_s t/L_n = $'+('{:06.4f}'.format(float(time))))
  ax1a.set_xlim(x[0][0],x[0][-1])
  ax1a.set_ylim(x[1][0],x[1][-1])

  #[ Transform to real space and store data.
  pData = np.fft.irfftn(cField, axes=(0,1))

  #[ Plot the actual data
  hpl1a.set_array(pData[:-1,:-1].ravel())
  hpl1a.autoscale()

  if showFigureFile:
    plt.show()
  else:
    plt.savefig(outputDir+('phi_xy-%d' % frameNo)+'.'+figureFileFormat, format=figureFileFormat)
  plt.close()


def LinDecWeights(Nek):
  Nekx, Neky = Nek
  errorWeights = zeros((Nekx, Neky))
  for lckx in range(Nekx):
    for lcky in range(Neky):
      d = np.sqrt(0.5*(lckx-1)**2/(Nekx-1)**2 + 0.5*(lcky-1)**2/(Neky-1)**2)
      errorWeights[lckx,lcky] = (1. - d) + 1e-3*d
  return errorWeights


def WrmsNorm(cField, weightFn=np.ones):
  weights = weightFn(cField.shape)
  flat = np.abs(cField.ravel())
  return np.real( 1./np.sqrt(flat.size) * np.sum(np.multiply(weights.ravel(), flat))**2)


#[ Variables we wish to compare.
varNames = ['phik']

#[ Handle command-line arguments
parser = argparse.ArgumentParser(description='Compare two MuSHrooM runs')

parser.add_argument('dir1', help='data directory for the reference solution')
parser.add_argument('dir2', help='data directory for the solution to compare against the reference')

parser.add_argument('-s', '--showFig', default=False,
                    action='store_true',
                    help='show figures instead of saving them')

parser.add_argument('-o', '--outputDir',
                    default='./',
                    help='output directory')

parser.add_argument('-fmt', '--figFormat',
                    default='png',
                    help='figure output format (png, pdf, ps, eps, svg)')

parser.add_argument('-n', '--numFrames', default=-1,
                    type=int,
                    help='the number of frames to compare (-1 => all)')

parser.add_argument('--plotfourier', default=False,
                    action='store_true',
                    help='output figures showing 2D plane in Fourier space (log of square amplitude of complex)')

parser.add_argument('--plotreal', default=False,
                    action='store_true',
                    help='output figures showing 2D plane in real space')

# parser.add_argument('-f', '--spectralFilter', default=None,
#                     type=str,
#                     choices=['gaussian', 'maeyama'],
#                     help='apply a low-pass SpectralFilter to the k-space data (box or gaussian or maeyama)')
# parser.add_argument('-fdelta', '--filterDelta', default=None,
#                     nargs=2,
#                     type=int,
#                     help='filter width (kx,ky)')

parser.add_argument('-ewt', '--errorWeight', default=None,
                    type=str,
                    choices=['linear', 'gaussian', 'maeyama'],
                    help='error weights to use with WrmsNorm')

args = parser.parse_args()
dataDirs = (args.dir1, args.dir2)
outputDir = args.outputDir
nFrames = args.numFrames
showFigureFile = args.showFig
figureFileFormat = args.figFormat

psu.checkMkdir(args.outputDir)

#[ Determine whether output files are ADIOS or binary (and that they match).
if (os.path.isfile(dataDirs[0]+"/"+varNames[0]+".bp")):
  adiosIO = True
elif (os.path.isfile(dataDirs[0]+"/"+varNames[0]+".bin")):
  adiosIO = False
else:
  sys.exit("\nError: nonexistent output file ",dataDirs[0],"/",varNames[0],"[.bp,.bin]")
if (os.path.isfile(dataDirs[1]+"/"+varNames[0]+".bp")):
  adiosIO2 = True
elif (os.path.isfile(dataDirs[1]+"/"+varNames[0]+".bin")):
  adiosIO2 = False
else:
  sys.exit("\nError: nonexistent output file ",dataDirs[1],"/",varNames[0],"[.bp,.bin]")
if (adiosIO != adiosIO2):
  sys.exit("\nError: non-conformant file types in ",dataDirs[0]," and ",dataDirs[1])

#[ Define IO functions (different for ADIOS and binary).
if adiosIO:
  psIO = psu.psIObp()
  fTypeExtension = '.bp'
else:
  psIO = psu.psIObin()
  fTypeExtension = '.bin'

nVars = len(varNames)    #[ Number of variables to post-process.

phik1 = ShroomVariable(dataDirs[0],varNames[0])
#[ For now alias all the class members to global variables
Nkx          = phik1.Nkx           #[ Number of distinct (absolute magnitude) kx modes.
Nky          = phik1.Nky           #[ Number of distinct (absolute magnitude) ky modes.
kxMin        = phik1.kxMin         #[ Minimum finite amplitude of kx modes.
kyMin        = phik1.kyMin         #[ Minimum finite amplitude of ky modes.
omSte        = phik1.omSte         #[ Parameter in the curvature drift frequency.
omde         = phik1.omde          #[ Parameter in the diamagnetic drift frequency.
tau          = phik1.tau           #[ Ratio of ion to electron temperature.
mu           = phik1.mu            #[ Square root of ion to electron mass ratio.
deltae       = phik1.deltae        #[ Electron (T_perp+T_par)/T.
deltaPerpe   = phik1.deltaPerpe    #[ Electron T_perp/T.
eta_e        = phik1.eta_e         #[ Electron L_n/L_T.
deltai       = phik1.deltai        #[ Ion (T_perp+T_par)/T.
deltaPerpi   = phik1.deltaPerpi    #[ Ion T_perp/T.
eta_i        = phik1.eta_i         #[ Ion L_n/L_T.
lambdaD      = phik1.lambdaD       #[ Normalized Debye length.
adiabaticElc = phik1.adiabaticElc  #[ Adiabatic electrons? =0 no, =else yes.
adiabaticIon = phik1.adiabaticIon  #[ Adiabatic ions? =0 no, =else yes.
HDmodel      = phik1.HDmodel       #[ Hyperdiffusion model.
Nk           = phik1.Nk
kMin         = phik1.kMin

#[ Ensure that second data directory has identical spatial/frequency grid
phik2 = ShroomVariable(dataDirs[1],varNames[0])
if (phik2.Nk[0] != Nk[0]) or (phik2.Nk[1] != Nk[1]):
  sys.exit(sys.argv[0]," error: mismatched Nk (",Nk," vs ",phik2.Nk,"). Terminating...")
if (phik2.kMin[0] != kMin[0]) or (phik2.kMin[1] != kMin[1]):
  sys.exit(sys.argv[0]," error: mismatched kMin (",kMin," vs ",phik2.kMin,"). Terminating...")

#[ Get the k-grid
Nekx, Neky, kx, ky, kxSq, kySq, kSq  = psu.kGrid(Nkx,kxMin,Nky,kyMin)
kMag  = np.sqrt(kSq)                    #[ Magnitude of wavenumber vector k.
Nek   = np.array([Nekx, Neky])          #[ Number of elements in k-space array.
Nx    = np.array([Nekx, (Neky-1)*2])    #[ Number of cells in de-aliased real space.

#[ Define real-space arrays.
if args.plotreal:
  Lx = np.array([2.0*np.pi/kMin[0], 2.0*np.pi/kMin[1]])  #[ Simulation box size.
  dx = np.array([Lx[0]/float(Nx[0]+np.mod(Nx[0],2)), Lx[1]/float(Nx[1]+np.mod(Nx[1],2))])  #[ Cell size.
  xreal  = np.array([np.arange(Nx[0])*dx[0]-Lx[0]/2.0+(1.0-np.mod(Nx[0],2))*0.5*dx[0],
           np.arange(Nx[1])*dx[1]-Lx[1]/2.0+(1.0-np.mod(Nx[1],2))*0.5*dx[1]])

#[ Define k-space arrays.
if args.plotfourier:
  xkspace = np.array([kx, ky])

#[ After establishing the number of elements in each direction
#[ assign value to byte-size variables (used if employing binary IO).
psIO.setByteSizes(Nek)

#[ Determine the number of frames saved in files.
pFrames = psIO.getFrameN(dataDirs[0]+"/"+varNames[0]+fTypeExtension)
pFrames2 = psIO.getFrameN(dataDirs[1]+"/"+varNames[0]+fTypeExtension)
if (pFrames != pFrames2):
  sys.exit("\nError: ",dataDirs[0]," and ",dataDirs[1]," have different numbers of frames (",pFrames," vs ",pFrames2,")")

#[ Determine number of frames to compare
if ((nFrames < 0) or (nFrames > pFrames)):
  nFrames = pFrames-1

print("\n Processing MuSHrooM simulations ",dataDirs[0]," and ",dataDirs[1],":")
print("   Nkx = ",Nkx,"  |  Nky = ",Nky,"  |  nFrames = ",nFrames,"\n")

#[ Extract the array of time stamps for plot
pTimes = psIO.getTimes(dataDirs[0]+"/"+varNames[0]+fTypeExtension,0,nFrames)
pTimes = pTimes[:nFrames]

#[ Initialize temporal 'error' and 'magnitude' arrays
varErrs = np.zeros([len(varNames),len(pTimes)], dtype=float)
varMags0 = np.zeros([len(varNames),len(pTimes)], dtype=float)
varMags1 = np.zeros([len(varNames),len(pTimes)], dtype=float)

#[ Initialize filters.
filt = None
# if args.spectralFilter == 'gaussian':
#   if args.filterDelta:
#     filt = GaussianFilter((Nekx,Neky), (kx,ky), kMin, args.filterDelta)
#   else:
#     filt = GaussianFilter((Nekx,Neky), (kx,ky), kMin)
#   filt.plotResponse(args.outputDir)
# elif args.spectralFilter == 'maeyama':
#   if args.filterDelta:
#     filt = MaeyamaFilter((Nekx,Neky), (kx,ky), kMin, args.filterDelta)
#   else:
#     filt = MaeyamaFilter((Nekx,Neky), (kx,ky), kMin)
#   filt.plotResponse(args.outputDir)

#[ Initialize weights.
weights = None
if args.errorWeight == 'gaussian':
  if filt is None and args.filterDelta:
    filt = GaussianFilter((Nekx,Neky), (kx,ky), kMin, delta=args.filterDelta)
  elif filt is None:
    filt = GaussianFilter((Nekx,Neky), (kx,ky), kMin)
  weights = filt.toEWT()
elif args.errorWeight == 'maeyama':
  if filt is None and args.filterDelta:
    filt = MaeyamaFilter((Nekx,Neky), (kx,ky), kMin, delta=args.filterDelta)
  elif filt is None:
    filt = MaeyamaFilter((Nekx,Neky), (kx,ky), kMin)
  weights = filt.toEWT()

#[ Loop over each variable, accumulating error/magnitude.
totalN = 0   #[ Total data count.
for iV, iVname in enumerate(varNames):

  #[ Open data files.
  fH0 = psIO.fOpen(dataDirs[0]+"/"+iVname+fTypeExtension)
  fH1 = psIO.fOpen(dataDirs[1]+"/"+iVname+fTypeExtension)

  #[ Read attributes (moves the position within the file for binary IO).
  psIO.skipAttributes(fH0)
  psIO.skipAttributes(fH1)

  #[ Squared magnitudes
  sqrs = zeros((nFrames,Nekx,Neky))
  lengths = zeros((nFrames,2))

  #[ Loop over frames.
  for iF in range(0,nFrames):

    #[ Read current frame data and time stamp.
    time0, inData0, _, _, _, _, _, _, _ = psIO.readFrame(fH0,iF)
    time1, inData1, _, _, _, _, _, _, _ = psIO.readFrame(fH1,iF)

    #[ Reorder the complex field in numpy fft format.
    #[ The kx grid is:
    #[ cField[0,] - the zero frequency term
    #[ cField[1:n/2,] - the positive-freq. terms
    #[ cField[n/2+1:,] - the negative-freq. terms in order of decreasingly negative freq.
    #[ The ky grid is:
    #[ cField[,0] - the zero frequency term
    #[ cField[,1:n/2] - the positive-freq. terms
    cField0 = psIO.reorderComplexField(inData0,Nk,Nek)
    cField1 = psIO.reorderComplexField(inData1,Nk,Nek)

    if filt:
      cField0 = filt.apply(cField0)
      cField1 = filt.apply(cField1)

    #[ Plot the frame
    if args.plotfourier:
      FourierFramePlot(spacing, cField0, iF, outputDir=args.outputDir)
    if args.plotreal:
      RealFramePlot(xreal, Nx, Nek, cField0, time0, iF, outputDir=args.outputDir)

    #[ Store data count for an individual frame.
    frameN = np.size(cField0)

    #[ Compute absolute difference, reference solution magnitude, and accumulate counter.
    if args.errorWeight == 'linear':
      varErrs[iV,iF] = WrmsNorm(cField0 - cField1, weightFn=LinDecWeights)
    elif args.errorWeight == 'gaussian' or args.errorWeight == 'maeyama':
      varErrs[iV,iF] = WrmsNorm(cField0 - cField1, weightFn=lambda x: weights)
    else:
      varErrs[iV,iF] = WrmsNorm(cField0 - cField1)
    varMags0[iV,iF] = np.sum(np.abs(cField0)**2)
    varMags1[iV,iF] = np.sum(np.abs(cField1)**2)
    totalN += np.size(cField0)

  #[ close data files for this variable.
  fH0.close()
  fH1.close()

# combine variable-specific errors & magnitudes.
timeErr = np.sum(varErrs, axis=0)
timeMag0 = np.sum(varMags0, axis=0)
timeMag1 = np.sum(varMags1, axis=0)

# convert error & magnitude values to RMS norms
for iV, iVname in enumerate(varNames):
  varErrs[iV] = np.sqrt(varErrs[iV]/frameN)
  varMags0[iV] = np.sqrt(varMags0[iV]/frameN)
  varMags1[iV] = np.sqrt(varMags1[iV]/frameN)
timeErr = np.sqrt(timeErr/frameN/len(varNames))
timeMag0 = np.sqrt(timeMag0/frameN/len(varNames))
timeMag1 = np.sqrt(timeMag1/frameN/len(varNames))
totalErr = np.sqrt(np.sum(timeErr)/totalN)
totalMag0 = np.sqrt(np.sum(timeMag0)/totalN)
totalMag1 = np.sqrt(np.sum(timeMag1)/totalN)

#[ Plot absolute error vs time for each variable
if (len(varNames)>1):
  plt.figure()
  for iV, iVname in enumerate(varNames):
    plt.semilogy(pTimes,varErrs[iV,:],label=iVname+" error")
  plt.xlabel("time"), plt.ylabel("error"), plt.legend()
  plt.title("Absolute difference: "+dataDirs[0]+" vs "+dataDirs[1])
  if showFigureFile:
    plt.show()
  else:
    plt.savefig('absErrs.'+figureFileFormat, format=figureFileFormat)

#[ Plot overall absolute error vs time
plt.figure()
plt.semilogy(pTimes,timeErr)
plt.xlabel("time"), plt.ylabel("error")
plt.title("Overall absolute difference: "+dataDirs[0]+" vs "+dataDirs[1])
if showFigureFile:
  plt.show()
else:
  plt.savefig('absErr.'+figureFileFormat, format=figureFileFormat)

#[ Plot relative error vs time for each variable
if (len(varNames)>1):
  plt.figure()
  for iV, iVname in enumerate(varNames):
    plt.semilogy(pTimes,varErrs[iV,:]/varMags0[iV,:],label=iVname+" error")
  plt.xlabel("time"), plt.ylabel("error"), plt.legend()
  plt.title("Relative difference: "+dataDirs[0]+" vs "+dataDirs[1])
  if showFigureFile:
    plt.show()
  else:
    plt.savefig('relErrs.'+figureFileFormat, format=figureFileFormat)

#[ Plot overall relative error vs time
plt.figure()
plt.semilogy(pTimes,timeErr/timeMag0)
plt.xlabel("time"), plt.ylabel("error")
plt.title("Overall relative difference: "+dataDirs[0]+" vs "+dataDirs[1])
if showFigureFile:
  plt.show()
else:
  plt.savefig('relErr.'+figureFileFormat, format=figureFileFormat)

#[ Plot solution magnitudes vs time for each variable
if (len(varNames)>1):
  plt.figure()
  for iV, iVname in enumerate(varNames):
    plt.semilogy(pTimes,varMags0[iV,:],label=iVname+" |ref|")
    plt.semilogy(pTimes,varMags1[iV,:],label=iVname+" |comp|")
  plt.xlabel("time"), plt.ylabel("|solutions|"), plt.legend()
  plt.title("Magnitudes: ref ("+dataDirs[0]+") vs comp ("+dataDirs[1]+")")
  if showFigureFile:
    plt.show()
  else:
    plt.savefig('mags.'+figureFileFormat, format=figureFileFormat)

#[ Plot solution magnitudes vs time for each variable
plt.figure()
plt.semilogy(pTimes,timeMag0,label="|ref|")
plt.semilogy(pTimes,timeMag1,label="|comp|")
plt.xlabel("time"), plt.ylabel("|solution|"), plt.legend()
plt.title("Overall magnitudes: ref ("+dataDirs[0]+") vs comp ("+dataDirs[1]+")")
if showFigureFile:
  plt.show()
else:
  plt.savefig('mag.'+figureFileFormat, format=figureFileFormat)

# output overall solution errors to screen
print(" Overall absolute error = %.16f" % totalErr)
print(" Overall absolute reference magnitude = %.16f" % totalMag0)
print(" Overall absolute comparison magnitude = %.16f" % totalMag1)
print(" Overall relative error = %.16f" % (totalErr/totalMag0))
