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
dx     = 0.35         # 'thickness' = pitch*diffusion [cm], nominal 0.35 (approximately where I see the peak pitch)
x      = 0.35         # mass per unit area, g/cm^2
Z      = 18           # atomic number of argon
A      = 39.948       # atomic mass of argon
I      = 188e-6       # mean excitation energy argon [MeV] from here: https://pdg.lbl.gov/2017/AtomicNuclearProperties/HTML/liquid_argon.html
hOmega = 22.85e-6     # Plasma energy, replaces I at betaGamma > 8 (muon momentum > 700) [MeV]
beta   = 1            # Relativistic v/c, assume 1

# Now define the energy range to assess
nbins  = 100 # arbitrarily chosen
EkVals = array('d', np.logspace(2,7,nbins)) # Evenly spaced in log10 basis, begin at Ek = 732 MeV, where the plasma energy assumption is valid

# Get the momentum from the kinetic energy
def muonMomentum(Ek):
  eTot = Ek + mmu
  p    = np.sqrt(np.square(eTot)-np.square(mmu))
  return p

# Now calculate beta*gamma for muons
def muonBetaGamma(Ek):
  bg = muonMomentum(Ek)/mmu
  return bg

# Get eta
def eta(Ek):
  e = (k*0.5)*(Z/A)*(x/beta)
  return e

# Function that calculates the energy loss at a given value
def energyLoss(Ek):

  # Actually calculate the energy loss
  p  = muonMomentum(Ek)
  e  = eta(Ek)
  b  = beta
  bg = muonBetaGamma(Ek)
  #print("Ek ", Ek, "p ", p, ", e ", e, ", b ", b, "bg ", bg)

  eLoss = e*(np.log((2*me*np.square(bg))/I) + np.log((e)/I) + j - b - delta)
  return eLoss/dx

# Now do some plotting
eLosses     = array('d')
bGammas     = array('d')
momenta     = array('d')

for Ek in EkVals:
  eLosses.append(energyLoss(Ek))
  bGammas.append(muonBetaGamma(Ek))
  momenta.append(muonMomentum(Ek))

# Setup the canvas
c = R.TCanvas("c","",900,900)
c.SetLeftMargin  (0.16)
c.SetRightMargin (0.02)
c.SetBottomMargin(0.12)
c.SetTopMargin   (0.02)

# First, just calculate the energy loss for one of our 'standard' muons
eTest = 9896
pTest = muonMomentum(eTest)
print("Most probable energy loss of a", eTest," MeV kinetic energy and ", pTest, " momentum muon is: ", energyLoss(eTest), " [MeV/cm]")

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

# Now momentum
g2 = R.TGraph(nbins,momenta,eLosses)
g2.SetTitle("")
g2.GetXaxis().SetTitle("Muon momentum [MeV/c]")
g2.GetYaxis().SetTitle("Most probable energy loss, [MeV/cm]")
g2.Draw()

c.SaveAs("momenta_vs_energyLoss"+tag+".png")
c.SaveAs("momenta_vs_energyLoss"+tag+".pdf")
c.SaveAs("momenta_vs_energyLoss"+tag+".root")

c.Clear();

