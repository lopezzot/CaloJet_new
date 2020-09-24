import os
import glob

energies = [30,50,70,90,150,250]

path = "/home/software/Calo/results/newresults/HepEvents_FTFPBERT/noBnoX0/jetscan/"
'''
for e in energies:
	string = "jetscan_"+str(e)
	f = "jetscan_"+str(e)+".root"
	cmnd1 = "time ./dohist"+" "+str(string)+" "+str(path+string)+" "+str(446.0)
	cmnd2 = "mv bonsai.jetscan_"+str(e)+" chi_0.446/jetscan_"+str(e)+".root"
	os.system(cmnd1)	
	os.system(cmnd2)
'''
types = ["hzjnbn","hznb","wwjl"]

path = "/home/software/Calo/results/newresults/HepEvents_FTFPBERT/noBnoX0/2j/"

for e in types:
	string = "jetscan_"+str(e)
	f = "jetscan_"+str(e)+".root"
	cmnd1 = "time ./dohist"+" "+str(string)+" "+str(path+string)+" "+str(445.0)
	os.system(cmnd1)
