#[ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ]#
#[
#[ Functions (utilities) used by pshroom.py.
#[
#[ Manaure Francisquez.
#[
#[ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ]#

from pylab import *
import numpy as np
import os
import adios as ad
import struct
import matplotlib.pyplot as plt
from shutil import copyfile
from scipy.optimize import curve_fit
import scipy.special as scsp


#[ Function to check existence and/or make directory.
def checkMkdir(dirIn):
  if not os.path.exists(os.path.dirname(dirIn)):
    try:
      os.makedirs(os.path.dirname(dirIn))
    except OSError as exc: # Guard against race condition
      if exc.errno != errno.EEXIST:
        raise
#[ .......... end of checkMkdir function ........... ]#

#[ This function finds the index of the grid point nearest to a given fix value.
def findNearestIndex(array,value):
  return (np.abs(array-value)).argmin()
#[ .......... end of findNearestIndex function ........... ]#

#[ Function computing the line equation (for fitting).
def lineEq(x,a,b):
  return a*x + b
#[ .......... end of lineEq function ........... ]#

#[.. LOADING DATASETS CORRESPONDING TO DATA FROM SMITH'S 1997 PoP FIGURES ..]#

#[ Load the particle flux from Smith's 1997 PoP figure 2.
def getSmithFlux():
  return np.loadtxt(open("./refData/Smith-particleFlux_1997.csv", "rb"), delimiter=",", skiprows=0)

#[ Load the spectrum of the conserved quantity Omega from Smith's 1997 PoP figure 3.
def getSmithOmegaSpectrum():
  return np.loadtxt(open("./refData/Smith-OmegaSpectrum_1997.csv", "rb"), delimiter=",", skiprows=0)

#[ Load the maximum growth rate from Smith's 1997 PoP figure 4.
def getSmithMaxGrowthRate():
  return np.loadtxt(open("./refData/Smith-maxGrowthRate_1997.csv", "rb"), delimiter=",", skiprows=0)

#[ Load the minimum growth rate from Smith's 1997 PoP figure 4.
def getSmithMinGrowthRate():
  return np.loadtxt(open("./refData/Smith-minGrowthRate_1997.csv", "rb"), delimiter=",", skiprows=0)

#[ .......... end of loading data from Smith's figures function ........... ]#


#[ Function that checks the existence of a variable.
#[ If it exists leave it alone, if not assign it a given value.
def checkExistAssign(varIn, valIn):
  try:
    varIn
  except NameError:
    return valIn
  else:
    return varIn
#[ .......... end of checkExistAssign function ........... ]#

#[ Generate k-grids.
def kGrid(Nkx, kxMin, Nky, kyMin):
  Nekx = (Nkx-1)*2+1                     #[ Number of elements along kx in k-space array.
  Neky = Nky                             #[ Number of elements along ky in k-space array.
  dkx  = kxMin
  dky  = kyMin
  kx   = np.zeros(Nekx)
  ky   = np.zeros(Neky)
  for ik in range(Nkx):
    kx[ik] = ik*dkx
  for ik in range(Nkx,Nekx):
    kx[ik] = -(Nkx-1-(ik-Nkx))*dkx
  for jk in range(Neky):
    ky[jk] = jk*dky
  kxSq = kx*kx
  kySq = ky*ky
  kSq  = np.zeros((Neky,Nekx))
  for ik in range(Nekx):
    for jk in range(Neky):
      kSq[jk,ik] = kxSq[ik] + kySq[jk]

  return Nekx, Neky, kx, ky, kxSq, kySq, kSq
#[ .......... End of genKgrid function ........... ]#

#[ Generate k-grid with kx centered at zero.
def kGridCentered(Nkx, kxMin, Nky, kyMin):
  Nekx   = (Nkx-1)*2+1                     #[ Number of elements along kx in k-space array.
  Neky   = Nky                             #[ Number of elements along ky in k-space array.
  dkx = kxMin
  dky = kyMin
  kx  = np.zeros(Nekx)
  ky  = np.zeros(Neky)
  for ik in range(Nekx):
    kx[ik] = (ik-(Nkx-1))*dkx
  for jk in range(Neky):
    ky[jk] = jk*dky
  kxSq = kx*kx
  kySq = ky*ky
  kSq  = np.zeros((Neky,Nekx))
  for ik in range(Nekx):
    for jk in range(Neky):
      kSq[jk,ik] = kxSq[ik] + kySq[jk]

  return Nekx, Neky, kx, ky, kxSq, kySq, kSq
#[ .......... End of genKgridCentered function ........... ]#

#[ Generate x-grids.
def xGridCentered(Nx, dx):
  x, y = np.zeros(Nx[0]), np.zeros(Nx[1])
  Lx   = dx*Nx
  for ik in range(Nx[0]):
    x[ik] = (ik-Nx[0]*0.5)*dx[0]
  for jk in range(Nx[1]):
    y[jk] = (jk-Nx[1]*0.5)*dx[1]

  return [x, y], Lx
#[ .......... End of genKgrid function ........... ]#


#[ .............. FLR functions ....................... ]#
def calcFLR(kxIn,kyIn,spec,tauIn,muIn,**kwargs):
  #[ Calculate the FLR operators.
  #[ Here 'spec' indicates the species (=0 electrons, =1 ions),
  #[ 'tauIn'=Ti0/Te0 and 'muIn'=sqrt(m_i/m_e).
  kxSq    = np.power(kxIn,2)
  kySq    = np.power(kyIn,2)
  kperpSq = np.zeros((np.size(kxIn),np.size(kyIn)))
  for ik in range(np.size(kxIn)):
    for jk in range(np.size(kyIn)):
      kperpSq[ik,jk] = np.add(kxSq[ik],kySq[jk])
  kperp   = np.sqrt(kperpSq)

  if spec == 0:
    krho = kperp/muIn   #[ Electrons.
  elif spec == 1:
    krho = kperp*tauIn  #[ Ions.

  b   = np.power(krho,2)
  bSq = np.power(b,2)

  Gamma0 = scsp.ive(0,b)
  Gamma1 = scsp.ive(1,b)
  avgJ0  = np.sqrt(Gamma0)

  GammaRat  = Gamma1/Gamma0
  hatLap    = b*(GammaRat-1.0)
  hatLapSq  = np.power(hatLap,2)
  hathatLap = b*(0.5*GammaRat-1.0)-0.25*bSq*(3.0+GammaRat)*(GammaRat-1.0)

  Db = 1.0+hathatLap-0.25*hatLapSq
  Nb = 1.0+hathatLap-0.5*hatLapSq

  Sb = (avgJ0*Nb)/Db

  flrDict = {}
  if 'only' in kwargs:
    #[ Output only desired arrays (saves memory).
    for iSt in kwargs['only']:
      if iSt == "Gamma0":
        flrDict["Gamma0"]    = Gamma0
      elif iSt == "avgJ0":
        flrDict["avgJ0"]     = avgJ0
      elif iSt == "hatLap":
        flrDict["hatLap"]    = hatLap
      elif iSt == "hathatLap":
        flrDict["hathatLap"] = hathatLap
      elif iSt == "Db":
        flrDict["Db"]        = Db
      elif iSt == "Nb":
        flrDict["Nb"]        = Nb
      elif iSt == "Sb":
        flrDict["Sb"]        = Sb
  else:
    #[ Output all FLR functions.
    flrDict["Gamma0"]    = Gamma0
    flrDict["avgJ0"]     = avgJ0
    flrDict["hatLap"]    = hatLap
    flrDict["hathatLap"] = hathatLap
    flrDict["Db"]        = Db
    flrDict["Nb"]        = Nb
    flrDict["Sb"]        = Sb
  return flrDict
#[ ........... End of FLR functions ................... ]#

#[ ........... Functions for plot manipulation ................. ]# 
def setXYlim(axisIn, adjustLims, defaultLims, userLims, **opArgs): 
  #[ Set the x-y limits of a plot drawn in axis 'axisIn'. The default limits
  #[ are given in 'defaultLims', while the user-requested limits are in
  #[ 'userLims'. Whether to use default or user lims is dictated by 'adjustLims'.
  #[ Optional arguments indicate whether to enforce 'detaultLims' or not.

  useDefX = opArgs.get('useDefaultX',True)
  useDefY = opArgs.get('useDefaultY',True)
  #[ x limits.
  if useDefX: 
    if adjustLims[0]:
      axisIn.set_xlim(defaultLims[0][0],defaultLims[0][1])
    else:
      axisIn.set_xlim(userLims[0][0],userLims[0][1])
  else:
    if not adjustLims[0]:
      axisIn.set_xlim(userLims[0][0],userLims[0][1])

  #[ y limits.
  if useDefY: 
    if adjustLims[1]:
      axisIn.set_ylim(defaultLims[1][0],defaultLims[1][1])
    else:
      axisIn.set_ylim(userLims[1][0],userLims[1][1])
  else:
    if not adjustLims[1]:
      axisIn.set_ylim(userLims[1][0],userLims[1][1])

#[ ........... End of unctions for plot manipulation ................. ]# 

#[ Define IO functions, depending on whether we are using ADIOS or binary IO.

#[ ............................................................................ ]#
#[ ............................ ADIOS IO CLASS ................................ ]#
#[ ............................................................................ ]#
class psIObp:
  def __init__(self):
    pass

  #[ Open a data file and return a handle to it.
  def fOpen(self,filePathAndName):
    return ad.file(filePathAndName)
  #[ .......... end of fOpen method ........... ]#

  #[ Open a data file and read variables saved as attributes.
  def getAttributes(self,filePathAndName):
    fH = ad.file(filePathAndName)
    kxMaxDyn      = fH.attr['kxMaxDyn'].value        #[ Maximum kx evolved.
    kyMaxDyn      = fH.attr['kyMaxDyn'].value        #[ Maximum ky evolved.
    omSte         = fH.attr['omSte'].value           #[ Parameter in the curvature drift frequency.
    omde          = fH.attr['omde'].value            #[ Parameter in the diamagnetic drift frequency.
    tau           = fH.attr['tau'].value             #[ Ratio of ion to electron temperature.
    muMass        = fH.attr['muMass'].value          #[ Square root of ion to electron mass ratio.
    deltae        = fH.attr['deltae'].value          #[ Electron (T_perp+T_par)/T.
    deltaPerpe    = fH.attr['deltaPerpe'].value      #[ Electron T_perp/T.
    eta_e         = fH.attr['eta_e'].value           #[ Electron L_n/L_T.
    deltai        = fH.attr['deltai'].value          #[ Ion (T_perp+T_par)/T.
    deltaPerpi    = fH.attr['deltaPerpi'].value      #[ Ion T_perp/T.
    eta_i         = fH.attr['eta_i'].value           #[ Ion L_n/L_T.
    lambdaD       = fH.attr['lambdaD'].value         #[ Normalized Debye length.
    kxDiffMin     = fH.attr['kxDiffMin'].value       #[ Minimum kx affected by diffusion.
    kyDiffMin     = fH.attr['kyDiffMin'].value       #[ Minimum ky affected by diffusion.
    tRateOutput   = fH.attr['tRateOutput'].value     #[ Time rate of outputs.
    tRateAdjustHD = fH.attr['tRateAdjustHD'].value   #[ Time rate of hyperdiffusion adjustments.
    tRateAdjustDt = fH.attr['tRateAdjustDt'].value   #[ Time rate of dt adjustments.
    adiabaticElc  = fH.attr['adiabaticElc'].value    #[ Adiabatic electrons? =0 no, else=yes.
    adiabaticIon  = fH.attr['adiabaticIon'].value    #[ Adiabatic ions? =0 no, else=yes.
    useMPI        = fH.attr['useMPI'].value          #[ Was MPI used?
    ioMode        = fH.attr['ioMode'].value          #[ =0 binary, =1 ADIOS.
    modelHD       = fH.attr['HDmodel'].value         #[ Hyperdiffusion model.
    STEPPER       = fH.attr['STEPPER'].value         #[ Time stepping algorithm.
    ADJUSTdt      = fH.attr['ADJUSTdt'].value        #[ Was dt dynamically adjusted?
    if 'Nkx' in fH.attr:
      Nkx         = fH.attr['Nkx'].value             #[ Number of distinct (nonnegative) kx's.
      Nky         = fH.attr['Nky'].value             #[ Number of distinct (nonnegative) ky's.
      kxMin       = fH.attr['kxMin'].value           #[ Minimum finite wavenumber along kx.
      kyMin       = fH.attr['kyMin'].value           #[ Minimum finite wavenumber along ky.
      fH.close()
      return {'Nkx':Nkx, 'Nky':Nky, 'kxMin':kxMin, 'kyMin':kyMin, 'kxMaxDyn':kxMaxDyn, 'kyMaxDyn':kyMaxDyn, \
              'omSte':omSte, 'omde':omde, 'tau':tau, 'muMass':muMass, 'deltae':deltae, \
              'deltaPerpe':deltaPerpe, 'eta_e':eta_e, 'deltai':deltai, 'deltaPerpi':deltaPerpi, 'eta_i':eta_i, 'lambdaD':lambdaD, \
              'kxDiffMin':kxDiffMin, 'kyDiffMin':kyDiffMin, \
              'tRateOutput':tRateOutput, 'tRateAdjustHD':tRateAdjustHD, 'tRateAdjustDt':tRateAdjustDt, 'adiabaticElc':adiabaticElc, \
              'adiabaticIon':adiabaticIon, 'useMPI':useMPI, 'ioMode':ioMode, 'HDmodel':modelHD, 'STEPPER':STEPPER, 'ADJUSTdt':ADJUSTdt}
    elif 'Nx' in fH.attr:
      Nx          = fH.attr['Nx'].value              #[ Number of cells in real space along x.
      Ny          = fH.attr['Ny'].value              #[ Number of cells in real space along y.
      dx          = fH.attr['dx'].value              #[ Real space cell length along x.
      dy          = fH.attr['dy'].value              #[ Real space cell length along y.
      return {'Nx':Nx, 'Ny':Ny, 'dx':dx, 'dy':dy, 'kxMaxDyn':kxMaxDyn, 'kyMaxDyn':kyMaxDyn, \
              'omSte':omSte, 'omde':omde, 'tau':tau, 'muMass':muMass, 'deltae':deltae, \
              'deltaPerpe':deltaPerpe, 'eta_e':eta_e, 'deltai':deltai, 'deltaPerpi':deltaPerpi, 'eta_i':eta_i, 'lambdaD':lambdaD, \
              'kxDiffMin':kxDiffMin, 'kyDiffMin':kyDiffMin, \
              'tRateOutput':tRateOutput, 'tRateAdjustHD':tRateAdjustHD, 'tRateAdjustDt':tRateAdjustDt, 'adiabaticElc':adiabaticElc, \
              'adiabaticIon':adiabaticIon, 'useMPI':useMPI, 'ioMode':ioMode, 'HDmodel':modelHD, 'STEPPER':STEPPER, 'ADJUSTdt':ADJUSTdt}
    elif 'Nxa' in fH.attr:
      Nxa         = fH.attr['Nxa'].value              #[ Number of cells in aliased real space along x.
      Nya         = fH.attr['Nya'].value              #[ Number of cells in aliased real space along y.
      dxa         = fH.attr['dxa'].value              #[ Aliased real space cell length along x.
      dya         = fH.attr['dya'].value              #[ Aliased real space cell length along y.
      return {'Nxa':Nxa, 'Nya':Nya, 'dxa':dxa, 'dya':dya, 'kxMaxDyn':kxMaxDyn, 'kyMaxDyn':kyMaxDyn, \
              'omSte':omSte, 'omde':omde, 'tau':tau, 'muMass':muMass, 'deltae':deltae, \
              'deltaPerpe':deltaPerpe, 'eta_e':eta_e, 'deltai':deltai, 'deltaPerpi':deltaPerpi, 'eta_i':eta_i, 'lambdaD':lambdaD, \
              'kxDiffMin':kxDiffMin, 'kyDiffMin':kyDiffMin, \
              'tRateOutput':tRateOutput, 'tRateAdjustHD':tRateAdjustHD, 'tRateAdjustDt':tRateAdjustDt, 'adiabaticElc':adiabaticElc, \
              'adiabaticIon':adiabaticIon, 'useMPI':useMPI, 'ioMode':ioMode, 'HDmodel':modelHD, 'STEPPER':STEPPER, 'ADJUSTdt':ADJUSTdt}
      fH.close()

  #[ .......... end of getAttributes method ........... ]#

  #[ Given the number of elements in a 2D complex field (and some assumed
  #[ file structure) establish the size of various data in bytes.
  def setByteSizes(self,NeksIn,**kwargs):
    #[ Only needed for binary IO.
    pass
  #[ .......... end of setByteSizes method ........... ]#

  #[ Read 2D complex field.
  #[ Determine how many time frames are saved in file (counting initial condition).
  def getFrameN(self,filePathAndName):
    fH     = self.fOpen(filePathAndName)
    #[ Frames saved (counting initial).
    if 'complexField' in fH.vars:
      frameN = fH['complexField'].nsteps
    elif 'realField' in fH.vars:
      frameN = fH['realField'].nsteps
    elif 'fourierRealField' in fH.vars:
      frameN = fH['fourierRealField'].nsteps
    elif 'realAliasedField' in fH.vars:
      frameN = fH['realAliasedField'].nsteps
    fH.close()
    return frameN
  #[ .......... end of getFrameN method ........... ]#

  #[ Extract simulation time of desired frames (as array).
  def getTimes(self,filePathAndName,iFrame,eFrame):
    fH    = self.fOpen(filePathAndName)
    times = np.array(fH['simulationTime'].read(from_steps=iFrame,nsteps=(eFrame-iFrame+1)))
    fH.close()
    return times
  #[ .......... end of getTimes method ........... ]#

  #[ Skip over the attributes to the position of the first frame.
  def skipAttributes(self,fHandle):
    #[ Only needed for binary IO.
    pass
  #[ .......... end of skipAttributes method ........... ]#

  #[ Read 2D complex field.
  def readComplexField(self,fHandle,tOffset):
    return fHandle['complexField'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readComplexField method ........... ]#

  #[ Read 2D real field.
  def readRealField(self,fHandle,tOffset):
    return fHandle['realField'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readRealField method ........... ]#

  #[ Read 2D Fourier real field.
  def readFourierRealField(self,fHandle,tOffset):
    return fHandle['fourierRealField'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of fourierReadRealField method ........... ]#

  #[ Read 2D real field (defined in aliased real space).
  def readRealAliasedField(self,fHandle,tOffset):
    return fHandle['realAliasedField'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readRealAliasedField method ........... ]#

  #[ Read simulation time of desired frame.
  def readSimulationTime(self,fHandle,tOffset):
    return fHandle['simulationTime'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readSimulationTime method ........... ]#

  #[ Read time step count at the desired frame.
  def readTimeSteps(self,fHandle,tOffset):
    return fHandle['timeSteps'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readTimeSteps method ........... ]#

  #[ Read the frame count at the desired frame.
  def readFramesOut(self,fHandle,tOffset):
    return fHandle['framesOut'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readFramesOut method ........... ]#

  #[ Read the number of hyperdiffusion adjustments at the desired frame.
  def readhdAdjusts(self,fHandle,tOffset):
    return fHandle['hdAdjusts'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readhdAdjusts method ........... ]#

  #[ Read the number of time step adjustments at the desired frame.
  def readdtAdjusts(self,fHandle,tOffset):
    return fHandle['dtAdjusts'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readdtAdjusts method ........... ]#

  #[ Read the number of time step adjustments at the desired frame.
  def readdt(self,fHandle,tOffset):
    return fHandle['dt'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readdt method ........... ]#

  #[ Read the hyperdiffusion order at the desired frame.
  def readhDiffOrder(self,fHandle,tOffset):
    return fHandle['hDiffOrder'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readhDiffOrder method ........... ]#

  #[ Read the electron density hyperdiffusion amplitude at the desired frame.
  def readhDiffne(self,fHandle,tOffset):
    return fHandle['hDiffne'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readhDiffne method ........... ]#

  #[ Read the electron temperature hyperdiffusion amplitude at the desired frame.
  def readhDiffTe(self,fHandle,tOffset):
    return fHandle['hDiffTe'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readhDiffTe method ........... ]#

  #[ Read the ion density hyperdiffusion amplitude at the desired frame.
  def readhDiffni(self,fHandle,tOffset):
    return fHandle['hDiffni'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readhDiffni method ........... ]#

  #[ Read the ion temperature hyperdiffusion amplitude at the desired frame.
  def readhDiffTi(self,fHandle,tOffset):
    return fHandle['hDiffTi'].read(from_steps=tOffset,nsteps=1)
  #[ .......... end of readhDiffTi method ........... ]#

  #[ Read whole frame (2D complex field and simulation time).
  def readFrame(self,fHandle,tOffset,**kwargs):
    fldDatStr  = kwargs.get('fieldType', 'complexField')
    simTime    = fHandle['simulationTime'].read(from_steps=tOffset,nsteps=1)
    fieldData  = fHandle[fldDatStr].read(from_steps=tOffset,nsteps=1)
    steps      = fHandle['timeSteps'].read(from_steps=tOffset,nsteps=1)
    frames     = fHandle['framesOut'].read(from_steps=tOffset,nsteps=1)
    hdChanges  = fHandle['hdAdjusts'].read(from_steps=tOffset,nsteps=1)
    dtChanges  = fHandle['dtAdjusts'].read(from_steps=tOffset,nsteps=1)
    timeStep   = fHandle['dt'].read(from_steps=tOffset,nsteps=1)
    hdOrder    = fHandle['hDiffOrder'].read(from_steps=tOffset,nsteps=1)
    hdAdenElc  = fHandle['hDiffne'].read(from_steps=tOffset,nsteps=1)
    hdAtempElc = fHandle['hDiffTe'].read(from_steps=tOffset,nsteps=1)
    hdAdenIon  = fHandle['hDiffni'].read(from_steps=tOffset,nsteps=1)
    hdAtempIon = fHandle['hDiffTi'].read(from_steps=tOffset,nsteps=1)
    return simTime, fieldData, steps, frames, hdChanges, \
           dtChanges, timeStep, hdOrder, \
           [hdAdenElc, hdAtempElc, hdAdenIon, hdAtempIon]
  #[ .......... end of readFrame method ........... ]#

  #[ Rearrange complex field into the order python's FFT expects.
  def reorderComplexField(self,cFieldIn,Nks,Neks):
    cDataOut = np.empty((Neks[0],Neks[1]), dtype='complex')
    #[ In this case the data is already in the right order.
    cDataOut = cFieldIn
    return cDataOut
  #[ .......... end of reorderComplexField method ........... ]#

  #[ Obtain the square norm of the complex field.
  def getSqNorm(self,cDataIn):
    return cDataIn[:,:].real**2+cDataIn[:,:].imag**2
  #[ .......... end of getSqNorm method ........... ]#

  #[ Obtain the square norm of the complex field element corresponding
  #[ to the iIn-th kx mode and the jIn-th ky mode.
  def getSqNormIJ(self,cDataIn,iIn,jIn):
    return cDataIn[iIn,jIn].real**2+cDataIn[iIn,jIn].imag**2
  #[ .......... end of getSqNormIJ method ........... ]#

  #[ Write post-processed data onto a '.bp' file.
  def dataToFile(self,attrs,fields,filePathAndName):
    #[ 'attrs' must be a dictionary with key-value pairs, each pair corresponding
    #[ to an attribute to be saved in the .bp file (named after their key).
    #[ Likewise, 'fields' must be a key-value dictionary with the variables to be saved.
    print(" Saving data to "+filePathAndName)
    print("")
    if not os.path.isfile(filePathAndName):
      #.ADIOS init.
      ad.init_noxml()
      ad.set_max_buffer_size(1000)
      groupId = ad.declare_group("pshroomOut", "")
      ad.select_method(groupId, "POSIX1", "", "")
      #.Define attributes.
      for attrK, attrV in attrs.items():
        if isinstance(attrV,int):
          ad.define_attribute_byvalue(groupId, attrK, "", attrV)
        elif isinstance(attrV,float):
          ad.define_attribute_byvalue(groupId, attrK, "", attrV)
        elif isinstance(attrV,complex):
          ad.define_attribute_byvalue(groupId, attrK, "", attrV)
      #.Define the ADIOS variables in this file.
      for fieldK, fieldV in fields.items():
        dataShape = np.shape(fieldV)
        if len(dataShape) > 0:
          sSize    = "{:d},".format(int(dataShape[0]))
          for d in range(1,len(dataShape)):
            sSize = sSize+"{:d},".format(int(dataShape[d]))
          sOffsets = "0,"*len(np.shape(fieldV))
        else:
          sSize    = ""
          sOffsets = ""
        if isinstance(fieldV,int) or (isinstance(fieldV,np.ndarray) and (fieldV.dtype=='int')):
          ad.define_var(groupId, fieldK, "", ad.DATATYPE.integer, sSize, sSize, sOffsets)
        elif isinstance(fieldV,float) or (isinstance(fieldV,np.ndarray) and (fieldV.dtype=='double')):
          ad.define_var(groupId, fieldK, "", ad.DATATYPE.double, sSize, sSize, sOffsets)
        elif isinstance(fieldV,complex) or (isinstance(fieldV,np.ndarray) and (fieldV.dtype=='complex')):
          ad.define_var(groupId, fieldK, "", ad.DATATYPE.double_complex, sSize, sSize, sOffsets)
      fh = ad.open("pshroomOut", filePathAndName, 'w')
      #.Write the variables.
      for fieldK, fieldV in fields.items():
        ad.write(fh, fieldK, fieldV)
      ad.close(fh)
      ##.Deal with weird file output where a '.bp.0' file is created.
      ##.IMPT: Actually removing that directory messes up appending.
      #if len(filePathAndName.split('/')) > 1:
      #    nm = filePathAndName.split('/')[-1]
      #else:
      #    nm = filePathAndName
      #shutil.move(filePathAndName + '.dir/' + nm + '.0', filePathAndName)
      #shutil.rmtree(filePathAndName + '.dir')
      ad.finalize()
    else:
      #.Append value to existing file.
      fw = ad.writer(filePathAndName, mode='a')
      fw.declare_group("pshroomOut", method="POSIX1")
      for fieldK, fieldV in fields.items():
        fw.define_var(fieldK,(np.size(fieldV)))
        fw[fieldK] = fieldV
      fw.close()
  #[ .......... end of dataToFile method ........... ]#



#[ ............................................................................ ]#
#[ ............................ BINARY IO CLASS ................................ ]#
#[ ............................................................................ ]#
class psIObin:
  def __init__(self):
    self.Neks                  = []       #[ Number of elements in 2D complex field.
    self.Nes                   = []       #[ Number of elements in 2D real field.
    self.padBytes              = 4        #[ Number of bytes of the side padding in binary files.
    self.intBytes              = 4        #[ Number of bytes of the side padding in binary files.
    self.doubleBytes           = 8        #[ Number of bytes in double precision float.
    self.complexFieldEs        = 0        #[ Total number of elements in kx-ky space.
    self.realFieldEs           = 0        #[ Total number of elements in kx-ky space.
    self.frameBytes            = 0        #[ Size of data written in bytes.
    self.complexFieldBytes     = 0        #[ Size of complex field data written in bytes.
    self.realFieldBytes        = 0        #[ Size of real field data written in bytes.
    self.fourierRealFieldBytes = 0        #[ Size of Fourier real field data written in bytes.
    self.realAliasedFieldBytes = 0        #[ Size of real field data (in aliased real space) written in bytes.
    self.sComplexFieldDs       = ''       #[ Number of doubles for struct.unpack command (complex field).
    self.sRealFieldDs          = ''       #[ Number of doubles for struct.unpack command (real field).
    self.sFourierRealFieldDs   = ''       #[ Number of doubles for struct.unpack command (Fourier real field).
    self.sRealAliasedFieldDs   = ''       #[ Number of doubles for struct.unpack command (aliased real field).
    self.attrBytes             = 9*self.intBytes \
                                +29*self.doubleBytes #[ Bytes taken up by attributes (top of the file).
    #[ Byte size of variables in a frame before the 2D complex field (i.e. simulation time,
    #[ time steps, (absolute) number of frames outputted, (absolute) number of hyperdiffusion
    #[ adjustments, (absolute) number of dt adjustments, dt, hDiffOrder and hDiffs (4 of them)).
    self.pre2DComplexFieldBytes     = 7*self.doubleBytes+4*self.intBytes
    self.pre2DRealFieldBytes        = 7*self.doubleBytes+4*self.intBytes
    self.pre2DFourierRealFieldBytes = 7*self.doubleBytes+4*self.intBytes
    self.pre2DRealAliasedFieldBytes = 7*self.doubleBytes+4*self.intBytes

  #[ Open a data file and return a handle to it.
  def fOpen(self,filePathAndName):
    return open(filePathAndName,'rb')
  #[ .......... end of fOpen method ........... ]#

  #[ Open a data file and read variables saved as attributes.
  def getAttributes(self,filePathAndName):
    fH  = open(filePathAndName,'rb')
    jnk = fH.read(self.padBytes)
    Nks           = np.squeeze(struct.unpack('2i',fH.read(2*self.intBytes)))     #[ Number of distinct (nonnegative) k's.
    kMins         = np.squeeze(struct.unpack('2d',fH.read(2*self.doubleBytes)))  #[ Minimum finite wavenumbers.
    kMaxDyns      = np.squeeze(struct.unpack('2d',fH.read(2*self.doubleBytes)))  #[ Maximum kx,ky evolved.
    omSte         = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Parameter in the curvature drift frequency.
    omde          = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Parameter in the diamagnetic drift frequency.
    tau           = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Ratio of ion to electron temperature.
    muMass        = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Square root of ion to electron mass ratio.
    deltae        = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Electron (T_perp+T_par)/T.
    deltaPerpe    = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Electron T_perp/T.
    eta_e         = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Electron L_n/L_T.
    deltai        = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Ion (T_perp+T_par)/T.
    deltaPerpi    = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Ion T_perp/T.
    eta_i         = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Ion L_n/L_T.
    lambdaD       = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Normalized Debye length.
    kDiffMins     = np.squeeze(struct.unpack('2d',fH.read(2*self.doubleBytes)))  #[ Minimum kx,ky affected by diffusion.
    tRateOutput   = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Time rate of outputs.
    tRateAdjustHD = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Time rate of hyperdiffusion adjustments.
    tRateAdjustDt = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))    #[ Time rate of dt adjustments.
    adiabaticElc  = np.squeeze(struct.unpack('1i',fH.read(self.intBytes)))       #[ Adiabatic electrons? =0 no, =else yes.
    adiabaticIon  = np.squeeze(struct.unpack('1i',fH.read(self.intBytes)))       #[ Adiabatic ions? =0 no, =else yes.
    useMPI        = np.squeeze(struct.unpack('1i',fH.read(self.intBytes)))       #[ Was MPI used?
    ioMode        = np.squeeze(struct.unpack('1i',fH.read(self.intBytes)))       #[ =0 binary, =1 ADIOS.
    modelHD       = np.squeeze(struct.unpack('1i',fH.read(self.intBytes)))       #[ Hyperdiffusion model.
    STEPPER       = np.squeeze(struct.unpack('1i',fH.read(self.intBytes)))       #[ Time stepping algorithm.
    ADJUSTdt      = np.squeeze(struct.unpack('1i',fH.read(self.intBytes)))       #[ Was dt dynamically adjusted?
    jnk = fH.read(self.padBytes)
    fH.close()
    return {'Nkx':Nks[0], 'Nky':Nks[1], 'kxMin':kMins[0], 'kyMin':kMins[1], 'kxMaxDyn':kMaxDyns[0], 'kyMaxDyn':kMaxDyns[1], \
            'omSte':omSte, 'omde':omde, 'tau':tau, 'muMass':muMass, 'deltae':deltae, \
            'deltaPerpe':deltaPerpe, 'eta_e':eta_e, 'deltai':deltai, 'deltaPerpi':deltaPerpi, 'eta_i':eta_i, 'lambdaD':lambdaD, \
            'kxDiffMin':kDiffMins[0], 'kyDiffMin':kDiffMins[1], \
            'tRateOutput':tRateOutput, 'tRateAdjustHD':tRateAdjustHD, 'tRateAdjustDt':tRateAdjustDt, 'adiabaticElc':adiabaticElc, \
            'adiabaticIon':adiabaticIon, 'useMPI':useMPI, 'ioMode':ioMode, 'HDmodel':modelHD, 'STEPPER':STEPPER, 'ADJUSTdt':ADJUSTdt}
  #[ .......... end of getAttributes method ........... ]#

  #[ Given the number of elements in a 2D complex field (and some assumed
  #[ file structure) establish the size of various data in bytes.
  def setByteSizes(self,NeksIn,**kwargs):
    self.Neks                  = NeksIn                                    #[ Number of Fourier-space elements along each direction.
    self.Nes                   = NeksIn                                    #[ Number of elements in real-space along each direction.
    self.complexFieldEs        = self.Neks[0]*self.Neks[1]                 #[ Total number of elements in kx-ky space.
    self.realFieldEs           = self.Nes[0]*self.Nes[1]                   #[ Total number of elements in x-y space.
    self.complexFieldBytes     = self.complexFieldEs*2*self.doubleBytes    #[ Size of complex field in bytes.
    self.realFieldBytes        = self.realFieldEs*self.doubleBytes         #[ Size of real field in bytes.
    self.fourierRealFieldBytes = self.complexFieldEs*self.doubleBytes      #[ Size of Fourier real field in bytes.
    self.realAliasedFieldBytes = self.realFieldEs*self.doubleBytes         #[ Size of real field (in aliased real space) in bytes.
    #[ Byte size of a frame (consisting of variables placed before
    #[ the 2D field, and the 2D field).
    fldDatStr  = kwargs.get('fieldType', 'complexField')
    if fldDatStr=='complexField':
      self.frameBytes = self.pre2DComplexFieldBytes+self.complexFieldBytes
    elif fldDatStr=='realField':
      self.frameBytes = self.pre2DRealFieldBytes+self.realFieldBytes
    elif fldDatStr=='fourierRealField':
      self.frameBytes = self.pre2DFourierRealFieldBytes+self.fourierRealFieldBytes
    elif fldDatStr=='realAliasedField':
      self.frameBytes = self.pre2DRealAliasedFieldBytes+self.realAliasedFieldBytes
    self.sComplexFieldDs     = str(self.complexFieldEs*2)+'d'    #[ Number of doubles for struct.unpack command (complex field).
    self.sRealFieldDs        = str(self.realFieldEs)+'d'         #[ Number of doubles for struct.unpack command (real field).
    self.sFourierRealFieldDs = str(self.complexFieldEs)+'d'      #[ Number of doubles for struct.unpack command (Fourier real field).
    self.sRealAliasedFieldDs = str(self.realFieldEs)+'d'         #[ Number of doubles for struct.unpack command (aliased real field).
  #[ .......... end of setByteSizes method ........... ]#

  #[ Determine how many time frames are saved in file (counting initial condition).
  def getFrameN(self,filePathAndName):
    frameN     = 0
    moreFrames = True

    fH  = open(filePathAndName,'rb')
    jnk = fH.read(2*self.padBytes+self.attrBytes)    #[ Skip the attributes.
    while moreFrames:
      frameData = fH.read(2*self.padBytes+self.frameBytes)    #[ Read this frame's binary data.
      if not frameData:
        moreFrames = False    #[ End of file.
      else:
        frameN = frameN+1
    fH.close()
    return frameN
  #[ .......... end of getFrameN method ........... ]#

  #[ Extract simulation time of desired frames (as array).
  def getTimes(self,filePathAndName,iFrame,eFrame):
    times      = np.array([])
    moreFrames = True

    fH     = open(filePathAndName,'rb')
    jnk    = fH.seek(2*self.padBytes+self.attrBytes                 #[ Skip the attributes.
                     +iFrame*(2*self.padBytes+self.frameBytes),0)  #[ Skip unwanted (tOffset) frames.
    cFrame = 0
    while moreFrames:
      jnk = fH.read(self.padBytes)    #[ Skip this frame's start pad.
      if not jnk:
        moreFrames = False    #[ End of file.
      else:
        cTime = np.squeeze(struct.unpack('1d',fH.read(self.doubleBytes)))
        if (cFrame<=eFrame):
          times = np.append(times,cTime)
          cFrame = cFrame+1
        else:
          moreFrames = False
        jnk   = fH.read(self.pre2DComplexFieldBytes-self.doubleBytes  #[ Skip other scalar variables.
                        +self.complexFieldBytes                       #[ Skip the 2D complex field data.
                        +self.padBytes)                               #[ Skip this frame's end pad.
    fH.close()
    return times
  #[ .......... end of getTimes method ........... ]#

  #[ Skip over the attributes to the position of the first frame.
  def skipAttributes(self,fHandle):
    jnk = fHandle.read(2*self.padBytes+self.attrBytes)    #[ Skip the attributes.
  #[ .......... end of skipAttributes method ........... ]#

  #[ Read 2D complex field.
  def readComplexField(self,fHandle,tOffset):
    jnk       = fHandle.seek(2*self.padBytes+self.attrBytes                #[ Skip the attributes.
                             +tOffset*(2*self.padBytes+self.frameBytes)    #[ Skip unwanted (tOffset) frames.
                             +self.padBytes+self.pre2DComplexFieldBytes,0) #[ Skip variables save in front of 2D complex field.
    fieldData = np.array(np.reshape(struct.unpack(self.sComplexFieldDs,
                fHandle.read(self.complexFieldBytes)),(2,self.Neks[1],self.Neks[0]),order='F'))
    return fieldData
  #[ .......... end of readComplexField method ........... ]#

  #[ Read 2D real field.
  def readRealField(self,fHandle,tOffset):
    jnk       = fHandle.seek(2*self.padBytes+self.attrBytes                #[ Skip the attributes.
                             +tOffset*(2*self.padBytes+self.frameBytes)    #[ Skip unwanted (tOffset) frames.
                             +self.padBytes+self.pre2DRealFieldBytes,0)    #[ Skip variables save in front of 2D real field.
    fieldData = np.array(np.reshape(struct.unpack(self.sRealFieldDs,
                fHandle.read(self.realFieldBytes)),(self.Nes[1],self.Nes[0]),order='F'))
    return fieldData
  #[ .......... end of readRealField method ........... ]#

  #[ Read 2D real field in Fourier space.
  def readFourierRealField(self,fHandle,tOffset):
    jnk       = fHandle.seek(2*self.padBytes+self.attrBytes                #[ Skip the attributes.
                             +tOffset*(2*self.padBytes+self.frameBytes)    #[ Skip unwanted (tOffset) frames.
                             +self.padBytes+self.pre2DFourierRealFieldBytes,0)    #[ Skip variables save in front of 2D Fourier real field.
    fieldData = np.array(np.reshape(struct.unpack(self.sFourierRealFieldDs,
                fHandle.read(self.fourierRealFieldBytes)),(self.Neks[1],self.Neks[0]),order='F'))
    return fieldData
  #[ .......... end of readRealField method ........... ]#

  #[ Read 2D real field (defined in aliased real space).
  def readRealAliasedField(self,fHandle,tOffset):
    jnk       = fHandle.seek(2*self.padBytes+self.attrBytes                #[ Skip the attributes.
                             +tOffset*(2*self.padBytes+self.frameBytes)    #[ Skip unwanted (tOffset) frames.
                             +self.padBytes+self.pre2DRealAliasedFieldBytes,0)    #[ Skip variables save in front of 2D aliased real field.
    fieldData = np.array(np.reshape(struct.unpack(self.sRealAliasedFieldDs,
                fHandle.read(self.realAliasedFieldBytes)),(self.Nes[1],self.Nes[0]),order='F'))
    return fieldData
  #[ .......... end of readRealAliasedField method ........... ]#

  #[ Read simulation time of desired frame.
  def readSimulationTime(self,fHandle,tOffset):
    jnk     = fHandle.seek(2*self.padBytes+self.attrBytes              #[ Skip the attributes.
                           +tOffset*(2*self.padBytes+self.frameBytes)  #[ Skip unwanted (tOffset) frames.
                           +self.padBytes,0)                           #[ Skip pad of this frame.
    simTime = np.squeeze(struct.unpack('1d',fHandle.read(self.doubleBytes)))
    return simeTime
  #[ .......... end of readSimulationTime method ........... ]#

  #[ Read whole frame (2D complex field and simulation time)
  #[ and skip padding at the end of this frame. This method
  #[ assumes attributes have been read and current file position
  #[ is at the beginning of this frame (tOffset not used).
  def readFrame(self,fHandle,tOffset,**kwargs):
    jnk = fHandle.read(self.padBytes)
    simTime    = np.squeeze(struct.unpack('1d',fHandle.read(self.doubleBytes)))
#    simTime    = fHandle.read(self.doubleBytes)
    simTime    = np.squeeze(struct.unpack('1d',simTime))
    steps      = np.squeeze(struct.unpack('1i',fHandle.read(self.intBytes)))
    frames     = np.squeeze(struct.unpack('1i',fHandle.read(self.intBytes)))
    hdChanges  = np.squeeze(struct.unpack('1i',fHandle.read(self.intBytes)))
    dtChanges  = np.squeeze(struct.unpack('1i',fHandle.read(self.intBytes)))
    timeStep   = np.squeeze(struct.unpack('1d',fHandle.read(self.doubleBytes)))
    hdOrder    = np.squeeze(struct.unpack('1d',fHandle.read(self.doubleBytes)))
    hdAdenElc  = np.squeeze(struct.unpack('1d',fHandle.read(self.doubleBytes)))
    hdAtempElc = np.squeeze(struct.unpack('1d',fHandle.read(self.doubleBytes)))
    hdAdenIon  = np.squeeze(struct.unpack('1d',fHandle.read(self.doubleBytes)))
    hdAtempIon = np.squeeze(struct.unpack('1d',fHandle.read(self.doubleBytes)))
    fldDatStr  = kwargs.get('fieldType', 'complexField')
    if fldDatStr=='complexField':
      fieldData  = np.array(np.reshape(struct.unpack(self.sComplexFieldDs,
                    fHandle.read(self.complexFieldBytes)),(2,self.Neks[1],self.Neks[0]),order='F'))
    elif fldDatStr=='realField':
      fieldData  = np.array(np.reshape(struct.unpack(self.sRealFieldDs,
                    fHandle.read(self.realFieldBytes)),(self.Nes[1],self.Nes[0]),order='F'))
    elif fldDatStr=='fourierRealField':
      fieldData  = np.array(np.reshape(struct.unpack(self.sFourierRealFieldDs,
                    fHandle.read(self.complexFieldBytes)),(self.Neks[1],self.Neks[0]),order='F'))
    elif fldDatStr=='realAliasedField':
      fieldData  = np.array(np.reshape(struct.unpack(self.sRealAliasedFieldDs,
                    fHandle.read(self.realFieldBytes)),(self.Nes[1],self.Nes[0]),order='F'))
    jnk = fHandle.read(self.padBytes)
    return simTime, fieldData, steps, frames, hdChanges, \
           dtChanges, timeStep, hdOrder, \
           [hdAdenElc, hdAtempElc, hdAdenIon, hdAtempIon]
  #[ .......... end of readFrame method ........... ]#

  #[ Rearrange complex field into the order python's FFT expects.
  def reorderComplexField(self,cFieldIn,Nks,Neks):
    cDataOut = np.empty((Neks[0],Neks[1]), dtype='complex')
    #[ Place the ky's in the second dimension.
    #[ First store the kx>=0 entries:
    for ik in range(Neks[0]):
      for jk in range(Neks[1]):
        cDataOut[ik,jk] = cFieldIn[0,jk,ik]+1j*cFieldIn[1,jk,ik]
    return cDataOut
  #[ .......... end of reorderComplexField method ........... ]#

  #[ Obtain the square norm of the complex field.
  def getSqNorm(self,cDataIn):
    return np.transpose(cDataIn[0,:,:]**2+cDataIn[1,:,:]**2)
  #[ .......... end of getSqNorm method ........... ]#

  #[ Obtain the square norm of the complex field element corresponding
  #[ to the iIn-th kx mode and the jIn-th ky mode.
  def getSqNormIJ(self,cDataIn,iIn,jIn):
    return cDataIn[0,jIn,iIn]**2+cDataIn[1,jIn,iIn]**2
  #[ .......... end of getSqNormIJ method ........... ]#

#[ ............................ END OF IO CLASSES ................................ ]#
