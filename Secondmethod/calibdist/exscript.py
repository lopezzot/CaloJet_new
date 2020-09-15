import os
import glob

energies = [30,50,70,90,150,250]
'''
path = "/home/software/Calo/results/newresults/HepEvents_FTFPBERT/noBnoX0/jetscan/"

for e in energies:
	string = "jetscan_"+str(e)
	f = "jetscan_"+str(e)+".root"
	cmnd1 = "time ./dohist"+" "+str(string)+" "+str(path+string)+" "+str(300.0)
	os.system(cmnd1)	
'''
types = ["hzjnbn","hznb","wwjl"]

path = "/home/software/Calo/results/newresults/HepEvents_FTFPBERT/noBnoX0/2j/"

for e in types:
	string = "jetscan_"+str(e)
	f = "jetscan_"+str(e)+".root"
	cmnd1 = "time ./dohist"+" "+str(string)+" "+str(path+string)+" "+str(300.0)
	os.system(cmnd1)
