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
import argparse
import getopt, sys, os.path
import numpy as np
import pshroomUtil as psu
import adios as ad
import matplotlib.colors as colors
import matplotlib.pyplot as plt
# from filters import *


def readFlux(fH, kVar):
  flux = fH['flux'].read(from_steps=0,nsteps=1)
  k = fH[kVar].read(from_steps=0,nsteps=1)
  return k, flux

# pointwise error
def errorFn(reference, other):
  return np.abs(reference-other)


#[ Fluxes to plot
fluxNames = {'heat-xTimeAv': 'heatFluxIon-xTimeAv-f0-1000',
             'heat-yTimeAv': 'heatFluxIon-yTimeAv-f0-1000',
             'n-xTimeAv': 'nFluxIon-xTimeAv-f0-1000',
             'n-yTimeAv': 'nFluxIon-yTimeAv-f0-1000'}

#[ Handle command-line arguments
parser = argparse.ArgumentParser(description='Compare fluxes from MuSHrooM runs')

parser.add_argument('reference_dir', help='data directory for the reference solution')
parser.add_argument('-c', '--compare-dirs', nargs='+', help='data directory for the solution to compare against the reference')
parser.add_argument('-l', '--labels', nargs='+', help='labels for compare dirs')
parser.add_argument('-o', '--outputDir',
                    default='./',
                    help='output directory')
parser.add_argument('-fmt', '--figFormat',
                    default='png',
                    help='figure output format (png, pdf, ps, eps, svg)')


args = parser.parse_args()
dataDirs = [args.reference_dir] + args.compare_dirs
plt_labels = args.labels
outputDir = args.outputDir
figureFileFormat = args.figFormat

psu.checkMkdir(args.outputDir)

heatFlux_xTimeAvs = []
heatFlux_yTimeAvs = []
nFlux_xTimeAvs = []
nFlux_yTimeAvs = []
kx = None
ky = None
maxKX = 0
maxKY = 0
minKX = 1000000
minKY = 1000000
for dataDir in dataDirs:
  #[ Open all the files
  handles = {}
  for key in fluxNames.keys():
    if (os.path.isfile(dataDir+"/"+fluxNames[key]+".bp")):
      adiosIO = True
    elif (os.path.isfile(dataDir+"/"+fluxNames[key]+".bin")):
      adiosIO = False
    else:
      sys.exit("\nError: nonexistent output file %s" % (dataDir+"/"+fluxNames[key]+"[.bp,.bin]"))

    #[ Define IO functions (different for ADIOS and binary).
    if adiosIO:
      psIO = psu.psIObp()
      fTypeExtension = '.bp'
    else:
      psIO = psu.psIObin()
      fTypeExtension = '.bin'

    handles[key] = psIO.fOpen(dataDir+"/"+fluxNames[key]+fTypeExtension)

  #[ Read
  ky, heatFlux_xTimeAv = readFlux(handles['heat-xTimeAv'], 'ky')
  kx, heatFlux_yTimeAv = readFlux(handles['heat-yTimeAv'], 'kx')
  ky, nFlux_xTimeAv = readFlux(handles['n-xTimeAv'], 'ky')
  kx, nFlux_yTimeAv = readFlux(handles['n-yTimeAv'], 'kx')

  maxKX = max(len(kx), maxKX)
  maxKY = max(len(ky), maxKY)
  minKX = min(len(kx), minKX)
  minKY = min(len(ky), minKY)

  heatFlux_xTimeAvs.append(heatFlux_xTimeAv)
  heatFlux_yTimeAvs.append(heatFlux_yTimeAv)
  nFlux_xTimeAvs.append(nFlux_xTimeAv)
  nFlux_yTimeAvs.append(nFlux_yTimeAv)

  #[ Close files
  for fH in handles.values():
    fH.close()

# Adjust array sizes for different grid sizes
for idx, flux in enumerate(list(heatFlux_xTimeAvs)):
  rfac = len(flux)//minKX
  heatFlux_xTimeAvs[idx] = flux[::rfac] 
for idx, flux in enumerate(list(heatFlux_yTimeAvs)):
  rfac = len(flux)//minKX
  heatFlux_yTimeAvs[idx] = flux[::rfac] 
for idx, flux in enumerate(list(nFlux_xTimeAvs)):
  rfac = len(flux)//minKX
  nFlux_xTimeAvs[idx] = flux[::rfac] 
for idx, flux in enumerate(list(nFlux_yTimeAvs)):
  rfac = len(flux)//minKX
  nFlux_yTimeAvs[idx] = flux[::rfac] 

if len(kx) == maxKX:
  rfacX = maxKX//minKX
  kx = kx[::rfacX]
if len(ky) == maxKY:
  rfacY = maxKY//minKY
  ky = ky[::rfacY]

print(kx)
print(ky)

#[ Plot heat x
fig = plt.figure()
for heatFlux_xTimeAv in heatFlux_xTimeAvs:
  plt.loglog(ky, heatFlux_xTimeAv)
plt.legend(['ref']+plt_labels)
plt.xlabel(r'$k_y \rho_s$')
plt.ylabel(r'$\langle\sum_{k_x}\,\hat{\mathcal{Q}}_{i}\rangle_t$')
plt.title('Heat Flux Spectrum')
plt.savefig(fluxNames['heat-xTimeAv']+'.'+figureFileFormat, format=figureFileFormat, bbox_inches="tight")

#[ Plot heat y
fig = plt.figure()
for heatFlux_yTimeAv in heatFlux_yTimeAvs:
  plt.loglog(kx, heatFlux_yTimeAv)
plt.legend(['ref']+plt_labels)
plt.xlabel(r'$k_x \rho_s$')
plt.ylabel(r'$\langle\sum_{k_y}\,\hat{\mathcal{Q}}_{i}\rangle_t$')
plt.title('Heat Flux Spectrum')
plt.savefig(fluxNames['heat-yTimeAv']+'.'+figureFileFormat, format=figureFileFormat, bbox_inches="tight")

#[ Plot particle x
fig = plt.figure()
for nFlux_xTimeAv in nFlux_xTimeAvs:
  plt.loglog(ky, nFlux_xTimeAv)
plt.legend(['ref']+plt_labels)
plt.xlabel(r'$k_y \rho_s$')
plt.ylabel(r'$\langle\sum_{k_x}\,\hat{\mathcal{Q}}_{i}\rangle_t$')
plt.tight_layout()
plt.title('Particle Flux Spectrum')
plt.savefig(fluxNames['n-xTimeAv']+'.'+figureFileFormat, format=figureFileFormat, bbox_inches="tight")

#[ Plot particle y
fig = plt.figure()
for nFlux_yTimeAv in nFlux_yTimeAvs:
  plt.loglog(kx, nFlux_yTimeAv)
plt.legend(['ref']+plt_labels)
plt.xlabel(r'$k_x \rho_s$')
plt.ylabel(r'$\langle\sum_{k_y}\,\hat{\mathcal{Q}}_{i}\rangle_t$')
plt.tight_layout()
plt.title('Particle Flux Spectrum')
plt.savefig(fluxNames['n-yTimeAv']+'.'+figureFileFormat, format=figureFileFormat, bbox_inches="tight")

#[ Plot heat x error
fig = plt.figure()
lg = []
for i, heatFlux_xTimeAv in enumerate(heatFlux_xTimeAvs[1:]):
  plt.loglog(ky, errorFn(heatFlux_xTimeAvs[0], heatFlux_xTimeAv))
  relerr = np.linalg.norm(heatFlux_xTimeAv-heatFlux_xTimeAvs[0])/np.linalg.norm(heatFlux_xTimeAvs[0])
  lg.append('%s, relerr=%.2g' % (plt_labels[i], relerr))
plt.legend(lg)
plt.xlabel(r'$k_y \rho_s$')
plt.ylabel(r'$\langle\sum_{k_x}\,\hat{\mathcal{Q}}_{i}\rangle_t$')
plt.title('Heat Flux Spectrum Error')
plt.savefig(fluxNames['heat-xTimeAv']+'-error'+'.'+figureFileFormat, format=figureFileFormat, bbox_inches="tight")

#[ Plot heat y error
fig = plt.figure()
lg = []
for i, heatFlux_yTimeAv in enumerate(heatFlux_yTimeAvs[1:]):
  plt.loglog(kx, errorFn(heatFlux_yTimeAvs[0], heatFlux_yTimeAv))
  relerr = np.linalg.norm(heatFlux_yTimeAv-heatFlux_yTimeAvs[0])/np.linalg.norm(heatFlux_yTimeAvs[0])
  lg.append('%s, relerr=%.2g' % (plt_labels[i], relerr))
plt.legend(lg)
plt.xlabel(r'$k_x \rho_s$')
plt.ylabel(r'$\langle\sum_{k_y}\,\hat{\mathcal{Q}}_{i}\rangle_t$')
plt.title('Heat Flux Spectrum Error')
plt.savefig(fluxNames['heat-yTimeAv']+'-error'+'.'+figureFileFormat, format=figureFileFormat, bbox_inches="tight")

#[ Plot nFlux x error
fig = plt.figure()
lg = []
for i, nFlux_xTimeAv in enumerate(nFlux_xTimeAvs[1:]):
  plt.loglog(ky, errorFn(nFlux_xTimeAvs[0], nFlux_xTimeAv))
  relerr = np.linalg.norm(nFlux_xTimeAv-nFlux_xTimeAvs[0])/np.linalg.norm(nFlux_xTimeAvs[0])
  lg.append('%s, relerr=%.2g' % (plt_labels[i], relerr))
plt.legend(lg)
plt.xlabel(r'$k_y \rho_s$')
plt.ylabel(r'$\langle\sum_{k_x}\,\hat{\mathcal{Q}}_{i}\rangle_t$')
plt.title('Particle Flux Spectrum Error')
plt.savefig(fluxNames['n-xTimeAv']+'-error'+'.'+figureFileFormat, format=figureFileFormat, bbox_inches="tight")

#[ Plot nFlux y error
fig = plt.figure()
lg = []
for i, nFlux_yTimeAv in enumerate(nFlux_yTimeAvs[1:]):
  plt.loglog(kx, errorFn(nFlux_yTimeAvs[0], nFlux_yTimeAv))
  relerr = np.linalg.norm(nFlux_yTimeAv-nFlux_yTimeAvs[0])/np.linalg.norm(nFlux_yTimeAvs[0])
  lg.append('%s, relerr=%.2g' % (plt_labels[i], relerr))
plt.legend(lg)
plt.xlabel(r'$k_x \rho_s$')
plt.ylabel(r'$\langle\sum_{k_y}\,\hat{\mathcal{Q}}_{i}\rangle_t$')
plt.title('Particle Flux Spectrum Error')
plt.savefig(fluxNames['n-yTimeAv']+'-error'+'.'+figureFileFormat, format=figureFileFormat, bbox_inches="tight")
