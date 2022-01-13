#!/usr/bin/env python

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

tag = "_100122_MPVs_GeV" # Tag for file naming, include a "_" at the start

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

me     = 0.5109989461 # electron mass [MeV/c^2]
mmu    = 105.6583755  # muon mass [MeV/c^2]
k      = 0.307075     # 4*pi*Na*re*me*c^{2} [MeV*cm/mol]
j      = 0.2          # from here: https://pdg.lbl.gov/2016/reviews/rpp2016-rev-passage-particles-matter.pdf (top of page 12)
dxPeak = 0.353        # 'thickness' = pitch*diffusion [cm], nominal 0.35 (approximately where I see the peak pitch, no diffision)
EkPeak = 7675         # Approximate peak kinetic energy in the distribution [MeV] from peak true total in my studies
Z      = 18           # atomic number of argon
A      = 39.948       # atomic mass of argon
I      = 188e-6       # mean excitation energy argon [MeV] from here: https://pdg.lbl.gov/2017/AtomicNuclearProperties/HTML/liquid_argon.html
hOmega = 22.85e-6     # plasma energy, replaces I at betaGamma > 8 (muon momentum > 700) [MeV]
beta   = 1            # relativistic v/c, assume 1
rho    = 1.3973       # argon density in g/cm^3

# Now define the energy range to assess
nbins     = 1000 # arbitrarily chosen
EkVals    = array('d', np.logspace(3.602,6.699,nbins))   # Log10-spaced kinetic energies, 4 GeV -> 5000 GeV
dxVals    = array('d', np.linspace(0.3,1,nbins)) # Linerly spaced pitch values from the DUNE FD distribution, 0.3 -> 1 cm

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
eLosses = array('d')
bGammas = array('d')
momenta = array('d')
eGeV    = array('d')
pitches = array('d')

for Ek in EkVals:
  eLosses.append(energyLoss(Ek,dxPeak))
  bGammas.append(muonBetaGamma(Ek))
  momenta.append(muonMomentum(Ek)/1000.)
  eGeV.append(Ek/1000.)

for p in dxVals:
  pitches.append(energyLoss(EkPeak,p))

# Setup the canvas
c = R.TCanvas("c","",900,900)
c.SetLeftMargin  (0.16)
c.SetRightMargin (0.02)
c.SetBottomMargin(0.12)
c.SetTopMargin   (0.02)

# First, just calculate the energy loss for one of our 'standard' muons: Peak energy in truth
eTest = 292264-105
pTest = muonMomentum(eTest)
print("Most probable energy loss of a", eTest," MeV kinetic energy and ", pTest, " momentum muon is: ", energyLoss(eTest,dxPeak), " [MeV/cm]")

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

g3 = R.TGraph(nbins,dxVals,pitches)
g3.SetTitle("")
g3.SetLineWidth(2)
g3.SetLineColor(palette[0])
g3.GetXaxis().SetTitle("Pitch [cm]")
g3.GetYaxis().SetTitle("Most probable energy loss, [MeV/cm]")
g3.Draw()

c1.SaveAs("pitch_vs_energyLoss"+tag+".png")
c1.SaveAs("pitch_vs_energyLoss"+tag+".pdf")
c1.SaveAs("pitch_vs_energyLoss"+tag+".root")

c1.Clear();

# Calculate the average pitch
totELossPitch = 0
for p in pitches:
  totELossPitch += p

avgELossPitch  = totELossPitch/len(pitches)

# Now get the value of the energy loss at the maximum and minimum pitches
# and calculate how much they vary from the average
minELoss = g3.Eval(dxVals[0])
maxEloss = g3.Eval(dxVals[len(dxVals)-1])
minVar = (avgELossPitch-minELoss)/avgELossPitch
maxVar = (maxEloss-avgELossPitch)/avgELossPitch

print("Average energy loss w.r.t pitches: ", avgELossPitch, ", minVar: ", minVar*100, "%, maxVar: ", maxVar*100, "%")

