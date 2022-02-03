#!/usr/bin/env python
#
# A script which plots the most probable energy loss equation from here: https://lar.bnl.gov/properties/pass.html
#
#   Plots are made for (dE/dx)_MPV as a function of:
#
#   - Kinetic energy
#   - Momentum
#   - BetaGamma
#
#     Note: All of the above plots set the pitch of the tracks to be a single value,
#           to modify what this 'peak' pitch is, change: dxPeak in the parameter input list
#
#   - Pitch
#
#     Note: The pitch distribution sets the kinetic energy of the tracks to be a single value
#           to modify what this 'peak' kinetic energy is, change: EkPeak in the parameter input list
#
#   The 'tag' parameter simply adds a string to the output file names
#
#   There are 2 print statements, the first prints the MP energy loss at a single value of kinetic energy and momentum
#   to change the printed values, modify 'EkPrint' in the parameter input list
#
#   Send any queries to: r.s.jones@sheffield.ac.uk
#

from array import *
from ctypes import *
import ROOT as R
import numpy as np
import math
import os
import sys

lineStyles = array("i",[1,2,5,7])
lineWidths = array("i",[2,3,3,3])
palette    = array("i",[R.kPink+5, R.kAzure+7, R.kTeal+5, R.kOrange+1])

tag = "_010222_MPVs_GeV" # Tag for file naming, include a "_" at the start

# Global plotting settings
R.gStyle.SetLabelFont(132, "X")
R.gStyle.SetLabelFont(132, "Y")
R.gStyle.SetLabelSize(0.04,"X")
R.gStyle.SetLabelSize(0.04,"Y")
R.gStyle.SetTitleFont(132, "X")
R.gStyle.SetTitleFont(132, "Y")
R.gStyle.SetTitleSize(0.05,"X")
R.gStyle.SetTitleSize(0.05,"Y")

# Parameters to input to the MPV energy loss equation
# Equation from here: https://lar.bnl.gov/properties/pass.html
me      = 0.5109989461 # electron mass [MeV/c^2]
mmu     = 105.6583755  # muon mass [MeV/c^2]
k       = 0.307075     # 4*pi*Na*re*me*c^{2} [MeV*cm/mol]
j       = 0.2          # from here: https://pdg.lbl.gov/2016/reviews/rpp2016-rev-passage-particles-matter.pdf (top of page 12)
dxPeak  = 0.55         # 'thickness' = pitch*diffusion [cm], nominal 0.55 (where I see the peak pitch, no diffision)
EkPeak  = 7675         # Approximate peak kinetic energy in the distribution [MeV] from peak true total in my studies
EkAvg   = 292264-105   # Average kinetic energy of the muon sample [MeV] (E-m)
Z       = 18           # atomic number of argon
A       = 39.948       # atomic mass of argon
I       = 188e-6       # mean excitation energy argon [MeV] from here: https://pdg.lbl.gov/2017/AtomicNuclearProperties/HTML/liquid_argon.html
hOmega  = 22.85e-6     # plasma energy, replaces I at betaGamma > 8 (muon momentum > 700) [MeV]
beta    = 1            # relativistic v/c, assume 1
rho     = 1.3973       # argon density in g/cm^3
EkPrint = EkAvg        # The kinetic energy to print

# Now define the pitch and energy ranges to assess
nbins     = 1000                                        # arbitrarily chosen
dxVals    = array('d', np.linspace(0.3,1,nbins))        # Linerly spaced pitch values from the DUNE FD distribution, 0.3 -> 1 cm
EkVals    = array('d', np.logspace(6e-1,6,nbins))       # Log10-spaced kinetic energies, 4 MeV -> 10000 GeV
#EkVals    = array('d', np.logspace(3.602,6.699,nbins))   # Log10-spaced kinetic energies, 4 GeV -> 5000 GeV

# Get the momentum from the kinetic energy
def muonMomentum(Ek):
  eTot = Ek + mmu
  p    = np.sqrt(np.square(eTot)-np.square(mmu))
  return p

# Now calculate beta*gamma for muons
def muonBetaGamma(Ek):
  bg = muonMomentum(Ek)/mmu
  return bg

# Calculate gamma from kinetic energy and mass
def muonGamma(Ek):
  g = (mmu+Ek)/mmu
  return g

# Calculate beta from gamma
def muonBeta(Ek):
  b = np.power(1.0-np.power(muonGamma(Ek),-2.0),0.5)
  return b

# Square beta
def muonBetaSq(Ek):
  bsq = np.power(muonBeta(Ek),2)
  return bsq

# Get eta
def muonXi(Ek,th):
  xi = (k*0.5)*(Z/A)*(th/muonBetaSq(Ek))
  return xi

def muonThickness(pitch):
  t = rho*pitch
  return t

# Calculate the density correction factor, taken from Gray's code: https://github.com/SBNSoftware/sbncode/blob/develop/sbncode/Calibration/notebook/lib/dedx.py
def deltaBetaGamma(Ek):

  # Get betaGamma
  bg = muonBetaGamma(Ek)
  b  = muonBeta(Ek)

  # Medium energy
  delta = 2.0*np.log(10)*np.log10(bg)-5.2146+0.19559*np.power(3.0-np.log10(bg),3.0)

  # low energy
  if np.log10(bg) < 0.2 or b < 1e-6:
    delta = 0.

  # high energy
  if np.log10(bg) > 3.0:
    delta = (2.0*np.log(10)*np.log10(bg)-5.2146)

  return delta

# Function that calculates the energy loss at a given value
def energyLoss(Ek, pitch):

  # Actually calculate the energy loss
  th  = muonThickness(pitch)
  p   = muonMomentum(Ek)
  xi  = muonXi(Ek, th)
  b   = muonBetaSq(Ek)
  bg  = muonBetaGamma(Ek)
  d   = deltaBetaGamma(Ek)

  eLoss = xi*(np.log((2*me*np.square(bg))/I) + np.log((xi)/I) + j - b - d)/pitch
  return eLoss

# Now do some plotting
eLosses = array('d') # Energy loss array
bGammas = array('d') # betaGamma array (for the x axis)
momenta = array('d') # momentum array in GeV (for the x axis)
eGeV    = array('d') # kinetic energy array in GeV (for the x axis)
dPAvg   = array('d') # Energy loss as a function of pitch with average KE (for the x axis)
dPPeak  = array('d') # Energy loss as a function of pitch pitch with peak KE(for the x axis)

for Ek in EkVals:
  eLosses.append(energyLoss(Ek,dxPeak))
  bGammas.append(muonBetaGamma(Ek))
  momenta.append(muonMomentum(Ek)/1000.)
  eGeV.append(Ek/1000.)

for p in dxVals:
  dPPeak.append(energyLoss(EkPeak,p))
  dPAvg.append(energyLoss(EkAvg,p))

# Setup the canvas
c = R.TCanvas("c","",900,900)
c.SetLeftMargin  (0.16)
c.SetRightMargin (0.02)
c.SetBottomMargin(0.12)
c.SetTopMargin   (0.02)

# First, just calculate the energy loss for one of our 'standard' muons: Average energy in truth
pPrint = muonMomentum(EkPrint)
print("Most probable energy loss of a", EkPrint," MeV kinetic energy and ", pPrint, " momentum muon is: ", energyLoss(EkPrint,dxPeak), " [MeV/cm]")

# Now plot the energy-dependent case
c.SetLogx()

g = R.TGraph(nbins,eGeV,eLosses)
g.SetTitle("")
g.SetLineWidth(2)
g.SetLineColor(palette[0])
g.GetXaxis().SetTitle("Kinetic energy, [GeV]")
g.GetYaxis().SetTitle("Most probable energy loss, [MeV/cm]")
g.Draw()

c.SaveAs("kineticEnergy_vs_energyLoss"+tag+".png")
c.SaveAs("kineticEnergy_vs_energyLoss"+tag+".pdf")
c.SaveAs("kineticEnergy_vs_energyLoss"+tag+".root")

c.Clear();

# Now betaGamma
g1 = R.TGraph(nbins,bGammas,eLosses)
g1.SetTitle("")
g1.SetLineWidth(2)
g1.SetLineColor(palette[0])
g1.GetXaxis().SetTitle("#beta#gamma")
g1.GetYaxis().SetTitle("Most probable energy loss, [MeV/cm]")
g1.Draw()

c.SaveAs("betaGamma_vs_energyLoss"+tag+".png")
c.SaveAs("betaGamma_vs_energyLoss"+tag+".pdf")
c.SaveAs("betaGamma_vs_energyLoss"+tag+".root")

c.Clear();

# Now momentum
g2 = R.TGraph(nbins,momenta,eLosses)
g2.SetTitle("")
g2.SetLineWidth(2)
g2.SetLineColor(palette[0])
g2.GetXaxis().SetTitle("Muon momentum [GeV/c]")
g2.GetYaxis().SetTitle("Most probable energy loss, [MeV/cm]")
g2.Draw()

c.SaveAs("momenta_vs_energyLoss"+tag+".png")
c.SaveAs("momenta_vs_energyLoss"+tag+".pdf")
c.SaveAs("momenta_vs_energyLoss"+tag+".root")

c.Clear();

# Now pitch dependence

# Setup the canvas for no logx
c1 = R.TCanvas("c1","",900,900)
c1.SetLeftMargin  (0.16)
c1.SetRightMargin (0.02)
c1.SetBottomMargin(0.12)
c1.SetTopMargin   (0.02)

# Pitch calculated from the peak energy
g3 = R.TGraph(nbins,dxVals,dPPeak)
g3.SetTitle("")
g3.SetLineWidth(2)
g3.SetLineColor(palette[0])
g3.GetXaxis().SetTitle("Pitch [cm]")
g3.GetYaxis().SetTitle("Most probable energy loss, [MeV/cm]")
g3.Draw()

c1.SaveAs("peak_pitch_vs_energyLoss"+tag+".png")
c1.SaveAs("peak_pitch_vs_energyLoss"+tag+".pdf")
c1.SaveAs("peak_pitch_vs_energyLoss"+tag+".root")

c1.Clear();

# Pitch calculated from the average energy
g4 = R.TGraph(nbins,dxVals,dPAvg)
g4.SetTitle("")
g4.SetLineWidth(2)
g4.SetLineColor(palette[0])
g4.GetXaxis().SetTitle("Pitch [cm]")
g4.GetYaxis().SetTitle("Most probable energy loss, [MeV/cm]")
g4.Draw()

c1.SaveAs("avg_pitch_vs_energyLoss"+tag+".png")
c1.SaveAs("avg_pitch_vs_energyLoss"+tag+".pdf")
c1.SaveAs("avg_pitch_vs_energyLoss"+tag+".root")

c1.Clear();

# Calculate the average pitch calculated from the peak energy
totELossPitchPeak = 0
for p in dPPeak:
  totELossPitchPeak += p

avgELossPitchPeak = totELossPitchPeak/len(dPPeak)

# Now get the value of the energy loss at the maximum and minimum pitches
# and calculate how much they vary from the average
minELoss = g3.Eval(dxVals[0])
maxEloss = g3.Eval(dxVals[len(dxVals)-1])
minVar = (avgELossPitchPeak-minELoss)/avgELossPitchPeak
maxVar = (maxEloss-avgELossPitchPeak)/avgELossPitchPeak

print("Average energy loss w.r.t pitches from peak energy: ", avgELossPitchPeak, ", minVar: ", minVar*100, "%, maxVar: ", maxVar*100, "%")

# Calculate the average pitch calculated from the avg energy
totELossPitchAvg = 0
for p in dPAvg:
  totELossPitchAvg += p

avgELossPitchAvg = totELossPitchAvg/len(dPAvg)

# Now get the value of the energy loss at the maximum and minimum pitches
# and calculate how much they vary from the average
minELoss = g3.Eval(dxVals[0])
maxEloss = g3.Eval(dxVals[len(dxVals)-1])
minVar = (avgELossPitchAvg-minELoss)/avgELossPitchAvg
maxVar = (maxEloss-avgELossPitchAvg)/avgELossPitchAvg

print("Average energy loss w.r.t pitches from avg energy: ", avgELossPitchAvg, ", minVar: ", minVar*100, "%, maxVar: ", maxVar*100, "%")

