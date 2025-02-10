#!/usr/bin/python

#---
# Image Analysis Pipeline: Automates the identification of all .nd2 files for analysis,  
# compiles a structured list, and initiates processing on SciCore’s HPC system.  
# Additionally, it verifies the presence of key result files and generates a streamlined  
# subset .csv file from the main output for further analysis.

#----
import numpy as np
from StringIO import StringIO
from optparse import OptionParser, OptionGroup
import re
import sys
import readline
import subprocess
import os

np.set_printoptions(threshold='nan')

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-t","--timelapse",action="store",dest="nr_timelapse",help="read data from filename",default="")
    parser.add_option("-f","--fields",action="store",dest="fields",help="read fields",default="")

    return parser.parse_args()
	

### Check command line arguments

def check_input_validity(options,args):

	if options.nr_timelapse == "":
		print "Error! You need to specify number of timelapse"
		sys.exit()
	if options.fields == "":
		print "Error! You need to specify fields"
		sys.exit()
		
if __name__ == "__main__":
	
	# Get command line arguments
	
	(options,args) = main()
	check_input_validity(options,args)
	

# IMPORT FILES
numberTL = int(options.nr_timelapse)

pathToFiles2 = subprocess.Popen("pwd", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,).stdout.read().strip().split("/")[-1].split("_directories")[0].split(".")[2]
pathToFiles1 = subprocess.Popen("pwd", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,).stdout.read().strip().split("/")[-1].split("_directories")[0].split(".")[1]
pathToFiles0 = subprocess.Popen("pwd", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,).stdout.read().strip().split("/")[-1].split("_directories")[0].split(".")[0]
pathToFilesLoopA = "../%s.%s.%s_data/%s.%s.%s_data_LoopAt/" %(pathToFiles0,pathToFiles1,pathToFiles2,pathToFiles0,pathToFiles1,pathToFiles2)
print(pathToFilesLoopA)
#print(pathToFilesLoopB)


# CHECK IF NUMBER OF FIELDS IN 48 AND 96 IS ACCURATE 

count_fields_LoopA = subprocess.Popen("ls %s*nd2 | wc | cut -f 3 -d ' ' " %(pathToFilesLoopA,), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,).stdout.read().strip()

# COUNT WELLS AND FIELDS


if options.fields == "all":
	
	file1 = subprocess.Popen("ls %s*.nd2 | grep  '_C'" %(pathToFilesLoopA,), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,).stdout.read().strip()



else: 
	fields2 = subprocess.Popen("cat %s | tr '\t' '\n' | sort | uniq " %(options.fields,), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,).stdout.read().strip()
	grepCode = ""
	for x in fields2.strip().split("\n"):
		grepCode += " -e '_%s_' " %(x,)
	print("Test",grepCode)
	file1 = subprocess.Popen("ls %s*.nd2 | grep %s  " %(pathToFilesLoopA,grepCode), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,).stdout.read().strip()	
	
	
file1A = file1.replace(pathToFilesLoopA, "")
print(file1A)
print("PASSED LINE 102")

file2 = np.genfromtxt(StringIO(file1A), dtype='S', delimiter='_', comments='%%%%%')
print(file2)
print("PASSED LINE 104")

fields = (file2.shape)[0]					#ASSESS NUMBER OF FIELDS
print(fields)

wells1 = file2[:,1]

wells2 = np.unique(wells1).shape[0]			#ASSESS NUMBER OF WELLS

#fiePwell1 = np.sort(file2[:,:3],axis=0)		#ASSESS MIN AND MAX OF FIELDS PER WELL 
fiePwell1 = file2[:,:3]
fiePwell2 = fiePwell1[fiePwell1[:,2].argsort()]
fiePwell3 = fiePwell2[fiePwell2[:,1].argsort()]
#print(fiePwell3)

fiePwell4 = ""
expWell = ""
counter = 0 
for row_fiePwell in fiePwell3: 
	expWellNew = "%s %s" %(row_fiePwell[0],row_fiePwell[1])
	#print(expWellNew)
	if expWellNew  == expWell:
		counter += 1
	else:
		fiePwell4 += "%s\t%d\n" %(expWell,counter)
		expWell = expWellNew
		counter = 1
fiePwell4 += "%s\t%d" %(expWell,counter)
fiePwell5 = np.genfromtxt(StringIO(fiePwell4), dtype=(int,float), delimiter='\t', comments='%%%%%')
fiePwell6 = fiePwell5[1:,1]
minFpW = np.amin(fiePwell6)
maxFpW = np.amax(fiePwell6)

print("\nSTATISTICS OF TIMELAPSE\nNumber of wells: %s\t\tMin/Max number of fields per well: %s/%s\t\tTotal number of fields: %s\t\tNumber of timelapses: %s"  %(wells2,minFpW,maxFpW,fields,numberTL))   
raw_input("Press Enter to continue or control Z to exit")


# CREATE ALL DIRECTORIES

#print(fiePwell1)
numberTL_list = range(numberTL)
#print(numberTL_list)

print("PASSED LINE 149")

print(file2)

#CHECK IF FILES ARE ALREADY PRESENT
count = 0
table1 = ""
for row_Fields in file2:

	fieldName = "%s.%s.%s_%s_%s_%s" %(row_Fields[0].split(".")[0],row_Fields[0].split(".")[1],row_Fields[0].split(".")[2],row_Fields[1],row_Fields[2],row_Fields[3])
	fieldName2 = "%s.%s.%s_%s_%s" %(row_Fields[0].split(".")[0],row_Fields[0].split(".")[1],row_Fields[0].split(".")[2],row_Fields[2],row_Fields[3])
	print(fieldName)
	"""
	try: 
		os.mkdir(fieldName)
		
	except: 
		print(fieldName," already exists")

	"""
	output1A = subprocess.Popen("ls %s/*.o " %(fieldName2,), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,).stdout.read().strip()
	output2A = subprocess.Popen("ls %s/POCbSd-*csv " %(fieldName2,), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,).stdout.read().strip()
	output2B = "%s/POCbSd-%s_table.csv" %(fieldName2,fieldName2)

	
	if output1A != "" and output2A == output2B:
		print("OUTPUT file ALREADY PRESENT:",fieldName)
		
		output1B = subprocess.Popen("ls %s/OCsimple*csv " %(fieldName2,), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,).stdout.read().strip()		
		if output1B == "":
			simp1=subprocess.Popen("cat  %s/POCbSd-*csv | ~/.local/bin/csvcut -c 'timestep','labelimage_oid','Predicted Class','Mean Intensity_1','Size in pixels'"  %(fieldName,) , stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).stdout.read().strip()
			out = open("%s/OCsimple_%s.csv" %(fieldName2,fieldName2)  ,"w")
			out.write(simp1)
			out.close()
			
		
			print("OCsimple calculated")	
	else:
		#GENERATE TABLE
		table1 += "%s %s %s %s %s %s\n" %(row_Fields[0].split(".")[0],row_Fields[0].split(".")[1],row_Fields[0].split(".")[2],row_Fields[1],row_Fields[2],row_Fields[3])	
		
		count += 1
		#print(count)	

out = open("parameters.txt","w")
out.write(table1)
out.close()

print(" Executing commands after line 204")
print(" Executing command A")
commandA = "sbatch --array=1-%d"  %(count,)

print(" Executing command B")
commandB = " /scicore/projects/rinfsci/Jovanovic/SCRIPTS/AJ_scripts/New_Camera/MOC_POC_POC.basicDefault_BF.MOCid_H5/2_launchSLURMjobs_ImageAnaylsis"
command1 = "%s%s" %(commandA,commandB,)

os.system(command1)
print(command1)		
print("CONGRATULATIONS %d Job(s) submitted !!!" %(count,))


