#!/usr/bin/env python

from array import *
import ROOT as R
import numpy as np
import math
import os
import sys

lineStyles = array("i",[1,2,5,7])
lineWidths = array("i",[2,3,3,3])
palette    = array("i",[R.kPink+5, R.kAzure+7, R.kTeal+5, R.kOrange+1])

tag = "_121221_MPVs" # Tag for file naming, include a "_" at the start

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
delta  = 0            # density correction factor, not applied here - need to find from simulation
x      = 1.0          # 'thickness' = pitch*diffusion [cm], nominal 0.3 (no pitch no diffusion, will make dependent on pitch properly soon)
Z      = 18           # atomic number of argon
A      = 39.948       # atomic mass of argon
I      = 188e-6       # mean excitation energy argon [MeV] from here: https://pdg.lbl.gov/2017/AtomicNuclearProperties/HTML/liquid_argon.html
hOmega = 22.85e-6     # Plasma energy, replaces I at betaGamma > 8 (muon momentum > 700) [MeV]

# Now define the energy range to assess
nbins  = 100 # arbitrarily chosen
EkVals = array('d', np.logspace(2.865,7,nbins)) # Evenly spaced in log10 basis, begin at Ek = 732 MeV, where the plasma energy assumption is valid

# Get the momentum from the kinetic energy
def muonMomentum(Ek):
  eTot = Ek + mmu
  p    = np.sqrt(np.square(eTot)-np.square(mmu))
  return p

# Now calculate beta*gamma for muons
def muonBetaGamma(Ek):
  bg = muonMomentum(Ek)/mmu
  return bg

# Get beta^2 from the kinetic energy
def muonBeta(Ek):
  b = (2*Ek)/mmu
  return b

# Get eta
def eta(Ek):
  e = (k*0.5)*(Z/A)*(1/muonBeta(Ek))
  return e

# Function that calculates the energy loss at a given value
def energyLoss(Ek):

  # Actually calculate the energy loss
  p  = muonMomentum(Ek)
  e  = eta(Ek)
  b  = muonBeta(Ek)
  bg = muonBetaGamma(Ek)
  #print("Ek ", Ek, "p ", p, ", e ", e, ", b ", b, "bg ", bg)

  eLoss = e*(np.log((2*me*np.square(bg))/I) + np.log((e*x)/I) + j + b - delta)
  return eLoss

# Function that calculates the energy loss at a given value
def energyLossSimp(Ek):

  # Actually calculate the energy loss
  e  = eta(Ek)
  c  = 12.325 # constant taken from here: https://meroli.web.cern.ch/lecture_landau_ionizing_particle.html

  eLoss = e*(12.325 + np.log((e*x)/I))
  return eLoss

# Now do some plotting
eLosses     = array('d')
eLossesSimp = array('d')
bGammas     = array('d')

for Ek in EkVals:
  eLossesSimp.append(energyLossSimp(Ek))
  eLosses.append(energyLoss(Ek))
  bGammas.append(muonBetaGamma(Ek))

# Setup the canvas
c = R.TCanvas("c","",900,900)
c.SetLeftMargin  (0.16)
c.SetRightMargin (0.02)
c.SetBottomMargin(0.12)
c.SetTopMargin   (0.02)

# First, just calculate the energy loss for one of our 'standard' muons
eTest = 700
print("Most probable energy loss of a", eTest," MeV kinetic energy muon is: ", energyLoss(eTest), " [MeV/cm]")
print("Simplified MP energy loss of a", eTest," MeV kinetic energy muon is: ", energyLossSimp(eTest), " [MeV/cm]")

# Now plot the energy-dependent case
c.SetLogx()

g = R.TGraph(nbins,EkVals,eLosses)
g.SetTitle("")
g.GetXaxis().SetTitle("Kinetic energy, [MeV]")
g.GetYaxis().SetTitle("Most probable energy loss, [MeV/cm]")
g.Draw()

c.SaveAs("kineticEnergy_vs_energyLoss"+tag+".png")
c.SaveAs("kineticEnergy_vs_energyLoss"+tag+".pdf")
c.SaveAs("kineticEnergy_vs_energyLoss"+tag+".root")

c.Clear();

# Now betaGamma
g1 = R.TGraph(nbins,bGammas,eLosses)
g1.SetTitle("")
g1.GetXaxis().SetTitle("#beta#gamma")
g1.GetYaxis().SetTitle("Most probable energy loss, [MeV/cm]")
g1.Draw()

c.SaveAs("betaGamma_vs_energyLoss"+tag+".png")
c.SaveAs("betaGamma_vs_energyLoss"+tag+".pdf")
c.SaveAs("betaGamma_vs_energyLoss"+tag+".root")

c.Clear();

g2 = R.TGraph(nbins,EkVals,eLossesSimp)
g2.SetTitle("")
g2.GetXaxis().SetTitle("Kinetic energy, [MeV]")
g2.GetYaxis().SetTitle("Most probable energy loss, [MeV/cm]")
g2.Draw()

c.SaveAs("kineticEnergy_vs_energyLoss_Simplified"+tag+".png")
c.SaveAs("kineticEnergy_vs_energyLoss_Simplified"+tag+".pdf")
c.SaveAs("kineticEnergy_vs_energyLoss_Simplified"+tag+".root")

c.Clear();

# Now betaGamma
g3 = R.TGraph(nbins,bGammas,eLossesSimp)
g3.SetTitle("")
g3.GetXaxis().SetTitle("#beta#gamma")
g3.GetYaxis().SetTitle("Most probable energy loss, [MeV/cm]")
g3.Draw()

c.SaveAs("betaGamma_vs_energyLoss_Simplified"+tag+".png")
c.SaveAs("betaGamma_vs_energyLoss_Simplified"+tag+".pdf")
c.SaveAs("betaGamma_vs_energyLoss_Simplified"+tag+".root")

c.Clear();

