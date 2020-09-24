from ROOT import TTree, TFile, TH2F, TH1F, TGraph, TGraphErrors, gStyle, TLine, TF1, TLorentzVector
import glob
from array import array
import os
import numpy as np
import math

def jetdisplay():
	outputfile = "zjj"
	displayfile = TFile(outputfile+".root","RECREATE")

	energies = [30,50,90, 150, 250]
	cut = [16.72, 31.115, 58.715, 105.335, 180.8]
	
	#for geant4.10.5.p01 FTFPBERT
	inputfiles = ["Results/noBnoX0/jetscan/jetscan_"+str(e)+".root" for e in energies]
	#end geant4.10.5.p01 FTFPBERT

	#for chi scan
	inputfiles = ["Results/noBnoX0/chiscan/chi_0.5/jetscan_"+str(e)+".root" for e in energies]
	#end chi scan

	arrayenergies = array('d', [x/2. for x in energies])
	arraydiffj1 = array('d')
	arraydiffj2 = array('d')
	arrayres = array('d')
	arraysqrtenergies = array('d', [1./((x/2)**0.5) for x in energies])
	
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
		graphtest = TH1F("test"+str(energies[counter]), "test"+str(energies[counter]), 80, -40., 40.)
		graphtest2 = TH1F("test2_"+str(energies[counter]), "test2"+str(energies[counter]), 80, -40., 40.)
		graphenergy = TH1F("energy"+str(energies[counter]), "energy"+str(energies[counter]), 200, 0., 100.)
		graphenergytruth = TH1F("energytruth"+str(energies[counter]), "energytruth"+str(energies[counter]), 200, 0., 100.) 
		graphjs = TH1F("energyjs"+str(energies[counter]), "energyjs"+str(energies[counter]), 200, 0., 100.) 
		graphjc = TH1F("energyjc"+str(energies[counter]), "energyjc"+str(energies[counter]), 200, 0., 100.) 
		histresolution = TH1F("res_"+str(energies[counter]), "res_"+str(energies[counter]), 100, -1., 1.)	
		graph_emcomp02 = TH1F("emcomp02_"+str(energies[counter]), "emcomp02"+str(energies[counter]), 80, -40, 40)
		graph_emcomp04 = TH1F("emcomp04_"+str(energies[counter]), "emcomp04"+str(energies[counter]), 80, -40, 40)
		graph_emcomp06 = TH1F("emcomp06_"+str(energies[counter]), "emcomp06"+str(energies[counter]), 80, -40, 40)
		graph_emcomp08 = TH1F("emcomp08_"+str(energies[counter]), "emcomp08"+str(energies[counter]), 80, -40, 40)
		graph_emcomp1 = TH1F("emcomp1_"+str(energies[counter]), "emcomp1"+str(energies[counter]), 80, -40, 40)
		
		scatterplot = TH2F("diff_"+str(energies[counter]),"diff_"+str(energies[counter]), 70, -20., 50., 70, -50., 20)
		scatterplotedep = TH2F("edep_"+str(energies[counter]), "edep_"+str(energies[counter]), 100, 0.0, 100.0, 100, 0.0, 100.0)
		


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

			cut1 =  nmuon==0 and nneu==0
			cut2 =  abs(j1t_eta)<2.0 and abs(j2t_eta)<2.0
			cut3 = eleak<0.1
			#cut3 = True	
			cut4 = j2s_E+j1s_E>cut[counter]
			#cut4= True
			#cut5 = abs(j1t_E-j2t_E)<5.
			#cut5 = abs(j1t_phi-j2t_phi)>0.1
			cut5 = True
			if cut1 and cut2 and cut3 and cut4 and cut5:
				#deltaj1 = 0.04406*j1r_E+0.1158
				#deltaj2 = 0.04406*j2r_E+0.1158
				#deltaj1 = 0.02825*j1r_E+0.4056
				#deltaj2 = 0.02825*j2r_E+0.4056 
				#deltaj1 = 0.04135*j1r_E+0.08789
				#deltaj2 = 0.04135*j2r_E+0.08789
				#deltaj1 = 0.07113*j1r_E+0.5201
				#deltaj2 = 0.07113*j2r_E+0.5201
				#deltaj1 = 0.07211*j1r_E+0.7122
				#deltaj2 = 0.07211*j2r_E+0.7122
				#deltaj1 = 0.4272+0.07402*j1r_E-(0.00007572*(j1r_E**2.0))
				#deltaj2 = 0.4272+0.07402*j2r_E-(0.00007572*(j2r_E**2.0))
				deltaj1 = 0.0
				deltaj2 = 0.0
				graphtest.Fill(j1r_E+deltaj1-j1t_E)
				graphtest.Fill(j2r_E+deltaj2-j2t_E)
				graphtest2.Fill(j2r_E+deltaj2-j2t_E)
				histresolution.Fill((j1r_E+deltaj1-j1t_E)/j1t_E)
				histresolution.Fill((j2r_E+deltaj2-j2t_E)/j2t_E)
				'''
				if (emcomp1+emcomp2)<0.2*90.:
					graph_emcomp02.Fill(j1r_E+deltaj1-j1t_E)
					graph_emcomp02.Fill(j2r_E+deltaj2-j2t_E)
				if 0.2*90.<emcomp1+emcomp2<0.4*90.:
					graph_emcomp04.Fill(j1r_E+deltaj1-j1t_E)
					graph_emcomp04.Fill(j2r_E+deltaj2-j2t_E)
				if 0.4*90.<emcomp1+emcomp2<0.6*90.:
					graph_emcomp06.Fill(j1r_E+deltaj1-j1t_E)
					graph_emcomp06.Fill(j2r_E+deltaj2-j2t_E)				
				if 0.6*90.<emcomp1+emcomp2<0.8*90.:
					graph_emcomp08.Fill(j1r_E+deltaj1-j1t_E)
					graph_emcomp08.Fill(j2r_E+deltaj2-j2t_E)				
				if 0.8*90.<emcomp1+emcomp2<90.:
					graph_emcomp1.Fill(j1r_E+deltaj1-j1t_E)
					graph_emcomp1.Fill(j2r_E+deltaj2-j2t_E)				
				'''
				'''
				a = np.sin(j1r_theta*180./math.pi)**2.+np.sin(j2r_theta*180./math.pi)**2.
				b = -2.*float(energies[counter])*np.sin(j2r_theta*180./math.pi)**2
				c = -(j1r_m**2.)*np.sin(j1r_theta*180./math.pi)**2.+(float(energies[counter])**2.)*np.sin(j2r_theta*180./math.pi)**2.-(j2r_m**2.)*np.sin(j2r_theta*180./math.pi)**2.
				delta = (b**2.-4.*a*c)
				'''
				graphenergy.Fill(j1r_E+deltaj1)
				graphenergy.Fill(j2r_E+deltaj2)
				graphenergytruth.Fill(j1t_E)
				graphenergytruth.Fill(j2t_E)
				graphjs.Fill(j2s_E)
				graphjs.Fill(j1s_E)
				graphjc.Fill(j2c_E)
				graphjc.Fill(j1c_E)
				scatterplot.Fill(j2r_E+deltaj2-j2t_E, j1r_E+deltaj1-j1t_E)
				scatterplotedep.Fill(edep, j1s_E+j2s_E)

		graphtest.Fit("gaus")
		arraydiffj1.append(graphtest.GetFunction("gaus").GetParameter(1))
		graphtest2.Fit("gaus")
		arraydiffj2.append(graphtest2.GetFunction("gaus").GetParameter(1))
		histresolution.Fit("gaus")
		arrayres.append(histresolution.GetFunction("gaus").GetParameter(2))
		histresolution.GetXaxis().SetTitle("E_{j}^{r} - E_{j}^{t}")
		histresolution.GetYaxis().SetTitle("Events")
		displayfile.cd()
		graphtest.Write()
		graphtest2.Write()
		graphenergy.Write()
		graphenergytruth.Write()
		graphjs.Write()
		graphjc.Write()
		#graph_emcomp02.Write()
		#graph_emcomp04.Write()
		#graph_emcomp06.Write()
		#graph_emcomp08.Write()
		#graph_emcomp1.Write()
		scatterplot.Write()
		scatterplotedep.Write()
		histresolution.Write()
	graphres = TGraph(len(arraysqrtenergies), arraysqrtenergies, arrayres)
	graphres.GetXaxis().SetTitle("1/\sqrt{E_{j}^{norm} (GeV)}}")
	graphres.GetYaxis().SetTitle("\sigma((E_{j}^{r} - E_{j}^{t})/E_{j}^{t})")
	#graphres.GetYaxis().SetRangeUser(0.0,0.12)
	#graphres.GetXaxis().SetRangeUser(0.0,0.28)
	graph = TGraph(len(arraydiffj1), arrayenergies, arraydiffj1)
	graph.Fit("pol2","","",10.,125.)
	graph.GetXaxis().SetTitle("E_{nom}")
	graph.GetYaxis().SetTitle("E_{j1}^{r}-E_{j1}^{t}")
	graph2 = TGraph(len(arraydiffj2), arrayenergies, arraydiffj2)
	graph2.GetXaxis().SetTitle("E_{nom}")
	graph2.GetYaxis().SetTitle("E_{j2}^{r}-E_{j2}^{t}")
	print arrayenergies, arraydiffj1
	graph.Write()
	graph2.Write()
	graphres.Write()

jetdisplay()