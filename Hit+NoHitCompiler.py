#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 14:54:53 2020

@author: adam
"""
#taxonomy of blast XML
import re


#filenumber = input("what file number would you like to analyse? (2-12):")
for filenumber in range(2,13):
    print("analysing file: " + str(filenumber))
    inputfilename = ("/home/adamh/thrip/BlastResults/E"+str(filenumber)+"BlastResults.txt")
    outputfilename = ("/home/adamh/thrip/taxonomyTXT/namesandlensE"+str(filenumber)+".txt")
    nohitsoutputfilename = ("/home/adamh/thrip/taxonomyTXT/nohitcontigsE"+str(filenumber)+".txt")
    
    #patterns
    querystringpattern = re.compile(r"Query=\s(\w*)\s\w*=\d*\s\w*=\d*.?\d*\slen=(\d*)(?:\n.*){4}\n(?:\*\*\*\*\*\s(No hits found)\s\*\*\*\*\*)?\n(?:\w*\d*.?\d*\s*(?:PREDICTED:\s)?([A-Z]\w*[\s|\W]\w*))?")
    organismsplitter = re.compile(r"([A-Z]\w*)\s(\w*)")
    
    organismlist = []
    
    counter = 0
    nohits = 0
    errors = 0

    #each result is searched using curl https://taxonomy.jgi-psf.org/name/")
    
    writefile = open(outputfilename, "w")
    nohitsfile = open(nohitsoutputfilename, "w")
    
    with open(inputfilename, 'r') as handle:
        contents = handle.read()
        matches = querystringpattern.finditer(contents)
        for match in matches:
            if match.group(3) == "No hits found":
                nohitstring = str(match.group(1) + "\t" + match.group(2) + "\n")
                nohitsfile.write(nohitstring)
                nohits += 1
            elif match.group(3) == None:
                organism = match.group(4)
                temp1 = organism.split()
                if len(temp1) == 1:
                    organism2search = temp1[0]
                else:   
                    organism2search = str(temp1[0] + "_" + temp1[1])
                nameandlen = str(organism2search + "\t" + match.group(2)  + "\t" + match.group(1)  + "\n")
                writefile.write(nameandlen)
                counter += 1
            else:
                errors += 1
    
    
    writefile.close()
    nohitsfile.close()
    
    
    print("filenumber: " + str(filenumber))
    print("Number of organisms identified: " + str(counter))
    print("Number without a hit: " + str(nohits))
    print("Number of errors: " + str(errors))
    print("end of code")




#	curl https://taxonomy.jgi-psf.org/simple/gi/${ginum} >> E${filenumber}taxonomyL.txt
