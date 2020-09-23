from ROOT import TTree, TFile, TH1F, TGraph, TGraphErrors, gStyle, TLine, TF1, TLorentzVector
import glob
from array import array
import os
import numpy as np
import math

def jetdisplay():
	outputfile = "hznb"
	displayfile = TFile(outputfile+".root","RECREATE")

	#for geant4.10.5 FTFPBERT
	inputfiles = ["Results/noBnoX0/2j_0.43/hznb_0.43.root"]
	#end geant4.10.5 FTFPBERT
	
	for counter, inputfile in enumerate(inputfiles):
		inputfile = TFile(inputfile)
		print "Analyzing: "+str(inputfile)+" \n"
		tree = TTree()
		inputfile.GetObject("MyTree", tree)	

		#graphEjet1 = TH1F("energyjet1", "energyjet1", 100, 0., 200.)
		#graphEjet2 = TH1F("energyjet2", "energyjet2", 100, 0., 200.)
		
		#graphEcherjet1 = TH1F("energycher1", "energycherjet", 100, 0., 200.)
		#graph3 = TH1F("energyscinjet", "energyscinjet", 100, 0., 200.)
		
		#graphmass = TH1F("mass_jet", "mass_jet", 100, 0., 200.)
		graphtest = TH1F("test", "test", 80, -40., 40.)
		graphenergy = TH1F("energy", "energy", 100, 60., 160.)
		graphenergytruth = TH1F("energytruth", "energytruth", 100, 60., 160.) 
		graphjs = TH1F("energyjs", "energyjs", 200, 0., 100.) 
		graphjc = TH1F("energyjc", "energyjc", 200, 0., 100.) 
		diff = TH1F("diffmass","diffmass",250,-25.,25.)
		#loop over events
		for Event in range(tree.GetEntries()):		

			tree.GetEntry(Event)	
			#print "Event "+str(Event)
			nmuon = tree.nmuon
			nneu = tree.nneu
			mjjr = tree.mjjr
			mjjt = tree.mjjt
			edep = tree.edep
			muene_che = tree.muene_che
			muene_sci = tree.muene_sci
			emcomp1 = tree.emcomp1
			emcomp2 = tree.emcomp2
			eleak = tree.eleak
			eleakn = tree.eleakn
			   
			j1t_E = tree.j1t_E
			j1t_m = tree.j1t_m
			j1t_theta = tree.j1t_theta
			j1t_pt = tree.j1t_pt
			j1t_eta = tree.j1t_eta
			j1t_phi = tree.j1t_phi
			j2t_E = tree.j2t_E
			j2t_m = tree.j2t_m
			j2t_theta = tree.j2t_theta
			j2t_pt = tree.j2t_pt
			j2t_eta = tree.j2t_eta
			j2t_phi = tree.j2t_phi	

			j1r_E = tree.j1r_E
			j1r_m = tree.j1r_m
			j1r_theta = tree.j1r_theta
			j1r_pt = tree.j1r_pt
			j1r_eta = tree.j1r_eta
			j1r_phi = tree.j1r_phi
			j2r_E = tree.j2r_E
			j2r_m = tree.j2r_m
			j2r_theta = tree.j2r_theta
			j2r_pt = tree.j2r_pt
			j2r_eta = tree.j2r_eta
			j2r_phi = tree.j2r_phi	

			#deltaj1 = 0.04406*j1r_E+0.1158
			#deltaj1 = 0.04135*j1r_E+0.08789
			#deltaj1 = 0.07113*j1r_E+0.5201
			deltaj1 = 0.07211*j1r_E+0.7122
			j1 = TLorentzVector()
			j1.SetPtEtaPhiE(j1r_pt, j1r_eta, j1r_phi, j1r_E)
			#deltaj2 = 0.04406*j2r_E+0.1158
			#deltaj2 = 0.04135*j2r_E+0.08789
			#deltaj2 = 0.07113*j2r_E+0.5201
			deltaj2 = 0.07211*j2r_E+0.7122
			j2 = TLorentzVector()
			j2.SetPtEtaPhiE(j2r_pt, j2r_eta, j2r_phi, j2r_E)
			newmass = (j1+j2).M()

			j1s_E = tree.j1s_E
			j1s_m = tree.j1s_m
			j1s_theta = tree.j1s_theta
			j1s_pt = tree.j1s_pt
			j1s_eta = tree.j1s_eta
			j1s_phi = tree.j1s_phi
			j2s_E = tree.j2s_E
			j2s_m = tree.j2s_m
			j2s_theta = tree.j2s_theta
			j2s_pt = tree.j2s_pt
			j2s_eta = tree.j2s_eta
			j2s_phi = tree.j2s_phi	

			j1c_E = tree.j1c_E
			j1c_m = tree.j1c_m
			j1c_theta = tree.j1c_theta
			j1c_pt = tree.j1c_pt
			j1c_eta = tree.j1c_eta
			j1c_phi = tree.j1c_phi
			j2c_E = tree.j2c_E
			j2c_m = tree.j2c_m
			j2c_theta = tree.j2c_theta
			j2c_pt = tree.j2c_pt
			j2c_eta = tree.j2c_eta
			j2c_phi = tree.j2c_phi	

			cut1 =  nmuon==0 and nneu==2
			#cut1 = True #for tails in H distribution
			cut2 =  abs(j1t_eta)<2.0 and abs(j2t_eta)<2.0
			cut3 = eleak<1.
			#cut3 = True #for tails in H distribution
			cut4 = j1t_E+j2t_E>85.0
			cut4 = True
			cut5 = edep>100
			cut5 = True
		
			if cut1 and cut2 and cut3 and cut4 and cut5:
				graphtest.Fill(j1r_E-j1t_E)
				graphtest.Fill(j2r_E-j2t_E)
				graphenergy.Fill(newmass)
				#graphenergy.Fill(j2r_E+deltaj2)
				#graphenergytruth.Fill(j1t_E)
				#graphenergytruth.Fill(j2t_E+j1t_E)
				graphenergytruth.Fill(mjjt)
				graphjs.Fill(j2s_E)
				graphjs.Fill(j1s_E)
				graphjc.Fill(j2c_E)
				graphjc.Fill(j1c_E)
				diff.Fill(newmass-mjjt)

		displayfile.cd()
		#graphtest.Write()
		scale = 1./graphenergy.Integral()
		graphenergy.Scale(scale)
		graphenergy.Write()
		graphenergytruth.Write()
		#graphjs.Write()
		#graphjc.Write()
		diff.Write()

jetdisplay()