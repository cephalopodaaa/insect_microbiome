#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 10:49:34 2020

@author: adam
"""
#interpreting taxanomic files

import re
from collections import defaultdict
import csv

counter = 0
errorcounter = 0
#correctpairs = 0

barchartdata = defaultdict(list)
barchartfilename = ("/home/adamh/thrip/taxonomyTXT/species_occurrence/barchartdata3.tsv")
barchartfile = open(barchartfilename, "wt")
barchart_writer = csv.writer(barchartfile, delimiter='\t')
barchart_writer.writerow(["dataset","catagory","length","covXlen","longestcontig","longestcontigcoverage","longestcontigspecies","species"])

LCratio = int(input("select LC ratio: "))

#file handling
#filenumber = input("what file would you like to analyse?")
for filenumber in range(2,13):
    inputfilename = ("/home/adamh/thrip/taxonomyTXT/E"+str(filenumber)+"taxa.txt")
    rawfile = open(inputfilename, "r")
    taxfile = rawfile.read()
    
    
    #regex prep
    pattern = re.compile(r"(\w*)\t(\d*)\t(\w*)\n{\n\s\s\s\"\w*\":\s{\n\s*\"\w*\":\s*\"\w*(?:\s\w*)?\",\n.*\n\s*\"\w*\":\s*\"\w*\",\n\s*\"(\w*)\":\s{\n\s*\"\w*\":\s\"(.*)\",\n\s*\"tax_id\":\s\d*\n\s*},(?:\n\s*\"(\w*)\":\s{\n\s*\"\w*\":\s\"(.*)\",\n\s*\"tax_id\":\s\d*\n\s*},?)?(?:\n\s*\"(\w*)\":\s{\n\s*\"\w*\":\s\"(.*)\",\n\s*\"tax_id\":\s\d*\n\s*},?)?(?:\n\s*\"(\w*)\":\s{\n\s*\"\w*\":\s\"(.*)\",\n\s*\"tax_id\":\s\d*\n\s*},?)?(?:\n\s*\"(\w*)\":\s{\n\s*\"\w*\":\s\"(.*)\",\n\s*\"tax_id\":\s\d*\n\s*},?)?(?:\n\s*\"(\w*)\":\s{\n\s*\"\w*\":\s\"(.*)\",\n\s*\"tax_id\":\s\d*\n\s*},?)?(?:\n\s*\"(\w*)\":\s{\n\s*\"\w*\":\s\"(.*)\",\n\s*\"tax_id\":\s\d*\n\s*},?)?(?:\n\s*\"(\w*)\":\s{\n\s*\"\w*\":\s\"(.*)\",\n\s*\"tax_id\":\s\d*\n\s*},?)?(?:\n\s*\"(\w*)\":\s{\n\s*\"\w*\":\s\"(.*)\",\n\s*\"tax_id\":\s\d*\n\s*},?)?")
    errorpattern = re.compile(r"{\n\s\s\s\"(\d*)\":\s{\n\s*\"error\":\s\"(.*)\",\n.*\s*\"type\":\s(\d*)")
    
    #output file
    outputfilename = ("/home/adamh/thrip/taxonomyTXT/species_occurrence/E"+str(filenumber)+"phyla.tsv")
    outputfile = open(outputfilename, "wt")
  
    
    taxlist = []  
    matches = pattern.findall(taxfile)          #counting matches
    for match in matches:
    #    print("new match:")
     #   print(match)
        taxlist.append(match)
        counter += 1
    
#    print("printing errors:")                #counting errors
    errors = errorpattern.findall(taxfile)
    for error in errors:
#        print("new error:")
#        print(error)
        errorcounter += 1
    
    rawfile.close()
    
    
    specieslist = []
    for i in taxlist:
  #      print(i)
   #     print("")
        temp = []
        temp.append(i[0])
        temp.append(int(i[1]))
        temp.append(i[2])
        specieslist.append(temp)
    



    #in file: species, frequency, length, contigs
    contigdict = defaultdict()
    firstlineskipper = True
    contigfilename = str("/home/adamh/thrip/E"+str(filenumber)+"_assembly/E"+str(filenumber)+"_coverage/E"+str(filenumber)+"_coveragesQ2Q3-COVs.txt")
    contigfile = open(contigfilename, "r")
    for line in contigfile:
        if firstlineskipper == True:
            firstlineskipper = False
            print(line)
        else:
            linesplit = []
            linesplit = line.split()
            contigname = linesplit[0]
            coverage = float(linesplit[1])
            contigdict[contigname] = coverage
    
  #  print(contigdict.keys())

        
    badcontigsfreq = 0
    badcontigslen = 0


    microbiome = defaultdict(list)
    for species,length,contig, a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2, g1, g2, h1, h2, i1, i2 in taxlist:
  #      print(species)
   #     print(length)
    #    print(contig)
        if float(length) > (LCratio*contigdict[contig]):
            contigcoverage = "," + str(contigdict[contig])
            strlength = "," + str(length)
            microbiome[contig] = {}
            microbiome[contig]['contig'] = contig
            microbiome[contig]['length'] = length
            microbiome[contig]['coverage'] = contigdict[contig]
            microbiome[contig]['coverageXlen'] = (float(length)*contigdict[contig])
            microbiome[contig][a1] = a2
            microbiome[contig][b1] = b2
            microbiome[contig][c1] = c2
            microbiome[contig][d1] = d2
            microbiome[contig][e1] = e2
            microbiome[contig][f1] = f2
            microbiome[contig][g1] = g2
            microbiome[contig][h1] = h2
            microbiome[contig][i1] = i2
            if 'species' not in microbiome:
                microbiome[contig]['species'] = species
        else:
            badcontigsfreq += 1
            badcontigslen += int(length)
    
    
    #creating a barchart
    barchartdata = defaultdict(int)
    barchartdata['thrip'] = {}
    barchartdata['thrip']['length'] = 0
    barchartdata['thrip']['covXlen'] = 0
    barchartdata['thrip']['longestcontig'] = 0
    barchartdata['thrip']['longestcontigcov'] = 0
    barchartdata['thrip']['longestcontigspecies'] = ""
    barchartdata['thrip']['species'] = []
    barchartdata['plant'] = {}
    barchartdata['plant']['length'] = 0
    barchartdata['plant']['covXlen'] = 0
    barchartdata['plant']['longestcontig'] = 0
    barchartdata['plant']['longestcontigcov'] = 0
    barchartdata['plant']['longestcontigspecies'] = ""
    barchartdata['plant']['species'] = []
    barchartdata['pseud'] = {}
    barchartdata['pseud']['length'] = 0
    barchartdata['pseud']['covXlen'] = 0
    barchartdata['pseud']['longestcontig'] = 0
    barchartdata['pseud']['longestcontigcov'] = 0
    barchartdata['pseud']['longestcontigspecies'] = ""
    barchartdata['pseud']['species'] = []
    barchartdata['paeni'] = {}
    barchartdata['paeni']['length'] = 0
    barchartdata['paeni']['covXlen'] = 0
    barchartdata['paeni']['longestcontig'] = 0
    barchartdata['paeni']['longestcontigcov'] = 0
    barchartdata['paeni']['longestcontigspecies'] = ""
    barchartdata['paeni']['species'] = []
    barchartdata['rhizo'] = {}
    barchartdata['rhizo']['length'] = 0
    barchartdata['rhizo']['covXlen'] = 0
    barchartdata['rhizo']['longestcontig'] = 0
    barchartdata['rhizo']['longestcontigcov'] = 0
    barchartdata['rhizo']['longestcontigspecies'] = ""
    barchartdata['rhizo']['species'] = []
    barchartdata['Wolbachia'] = {}
    barchartdata['Wolbachia']['length'] = 0
    barchartdata['Wolbachia']['covXlen'] = 0
    barchartdata['Wolbachia']['longestcontig'] = 0
    barchartdata['Wolbachia']['longestcontigcov'] = 0
    barchartdata['Wolbachia']['longestcontigspecies'] = ""
    barchartdata['Wolbachia']['species'] = []

   
    with open(outputfilename, "wt") as outputfile:
        tsv_writer = csv.writer(outputfile, delimiter='\t')
        tsv_writer.writerow(["contig","domain","kingdom","phylum","class","order","family","genus","species"])
        for a, values in microbiome.items():
            temp1 = []
            length = values['length']
            temp1.append(values['contig'])
            temp1.append(values['length'])
            temp1.append(values['coverageXlen'])
            if 'superkingdom' in values:
                temp1.append(values['superkingdom'])
            else:
                temp1.append("none")
            if 'kingdom' in values:
                temp1.append(values['kingdom'])
                if values['kingdom'] == 'Viridiplantae':
                    barchartdata['plant']['length'] += int(length)
                    barchartdata['plant']['covXlen'] += float(length)*contigdict[values['contig']]
                    if 'species' in values:
                        if values['species'] not in barchartdata['plant']['species']:
                            barchartdata['plant']['species'].append(values['species'])
                    if int(values['length']) >= int(barchartdata['plant']['longestcontig']):
                        barchartdata['plant']['longestcontig'] = int(values['length'])
                        barchartdata['plant']['longestcontigcov'] = float(values['coverage'])
                    if 'species' in values:
                        barchartdata['plant']['longestcontigspecies'] = values['species']
                    else:
                        barchartdata['plant']['longestcontigspecies'] = "unknown"
            else:
                temp1.append("none")
            if 'phylum' in values:
                temp1.append(values['phylum'])
            else:
                temp1.append("none")
            if 'class' in values:
                temp1.append(values['class'])
            else:
                temp1.append("none")
            if 'order' in values:
                temp1.append(values['order'])
            else:
                temp1.append("none")
            if 'family' in values:
                temp1.append(values['family'])
                if values['family'] == 'Thripidae':
                    barchartdata['thrip']['length'] += int(length)
                    barchartdata['thrip']['covXlen'] += float(length)*contigdict[values['contig']]
                    if 'species' in values:
                        if values['species'] not in barchartdata['thrip']['species']:
                            barchartdata['thrip']['species'].append(values['species'])
                    if int(values['length']) >= int(barchartdata['thrip']['longestcontig']):
                        barchartdata['thrip']['longestcontig'] = int(values['length'])
                        barchartdata['thrip']['longestcontigcov'] = float(values['coverage'])
                        if 'species' in values:
                            barchartdata['thrip']['longestcontigspecies'] = values['species']
                        else:
                            barchartdata['thrip']['longestcontigspecies'] = "unknown"
                if values['family'] == 'Pseudomonadaceae':
                    barchartdata['pseud']['length'] += int(length)
                    barchartdata['pseud']['covXlen'] += float(length)*contigdict[values['contig']]
                    if 'species' in values:
                        if values['species'] not in barchartdata['pseud']['species']:
                            barchartdata['pseud']['species'].append(values['species'])
                    if int(values['length']) >= int(barchartdata['pseud']['longestcontig']):
                        barchartdata['pseud']['longestcontig'] = int(values['length'])
                        barchartdata['pseud']['longestcontigcov'] = float(values['coverage'])
                        if 'species' in values:
                            barchartdata['pseud']['longestcontigspecies'] = values['species']
                        else:
                            barchartdata['pseud']['longestcontigspecies'] = "unknown"
                if values['family'] == 'Paenibacillaceae':
                    barchartdata['paeni']['length'] += int(length)
                    barchartdata['paeni']['covXlen'] += float(length)*contigdict[values['contig']]
                    if 'species' in values:
                        if values['species'] not in barchartdata['paeni']['species']:
                            barchartdata['paeni']['species'].append(values['species'])
                    if int(values['length']) >= int(barchartdata['paeni']['longestcontig']):
                        barchartdata['paeni']['longestcontig'] = int(values['length'])
                        barchartdata['paeni']['longestcontigcov'] = float(values['coverage'])
                        if 'species' in values:
                            barchartdata['paeni']['longestcontigspecies'] = values['species']
                        else:
                            barchartdata['paeni']['longestcontigspecies'] = "unknown"
                if values['family'] == 'Rhizobiaceae':
                    barchartdata['rhizo']['length'] += int(length)
                    barchartdata['rhizo']['covXlen'] += float(length)*contigdict[values['contig']]
                    if 'species' in values:
                        if values['species'] not in barchartdata['rhizo']['species']:
                            barchartdata['rhizo']['species'].append(values['species'])
                    if int(values['length']) >= int(barchartdata['rhizo']['longestcontig']):
                        barchartdata['rhizo']['longestcontig'] = int(values['length'])
                        barchartdata['rhizo']['longestcontigcov'] = float(values['coverage'])
                        if 'species' in values:
                            barchartdata['rhizo']['longestcontigspecies'] = values['species']
                        else:
                            barchartdata['rhizo']['longestcontigspecies'] = "unknown"
            else:
                temp1.append("none")
            if 'genus' in values:
                temp1.append(values['genus'])
                if values['genus'] == 'Wolbachia':
                    barchartdata['Wolbachia']['length'] += int(length)
                    barchartdata['Wolbachia']['covXlen'] += float(length)*contigdict[values['contig']]
                    if 'species' in values:
                        if values['species'] not in barchartdata['Wolbachia']['species']:
                            barchartdata['Wolbachia']['species'].append(values['species'])
                    if int(values['length']) >= int(barchartdata['Wolbachia']['longestcontig']):
                        barchartdata['Wolbachia']['longestcontig'] = int(values['length'])
                        barchartdata['Wolbachia']['longestcontigcov'] = float(values['coverage'])
                        if 'species' in values:
                            barchartdata['Wolbachia']['longestcontigspecies'] = values['species']
                        else:
                            barchartdata['Wolbachia']['longestcontigspecies'] = "unknown"               
            else:
                temp1.append("none")
            if 'species' in values:
                temp1.append(values['species'])
            else:
                temp1.append("none")
            tsv_writer.writerow(temp1)
           
            
            
            
        
            




    #writing dictionary to file
    
        for catagory, value in barchartdata.items():
            temp2 = []
            temp2.append(filenumber)
            temp2.append(catagory)
            temp2.append(value['length'])
            temp2.append(value['covXlen'])
            temp2.append(value['longestcontig'])
            temp2.append(value['longestcontigcov'])
            temp2.append(value['longestcontigspecies'])
            temp2.append(value['species'])
            barchart_writer.writerow(temp2)
        
    


    
    print("sequences identified: "+str(counter))
    print("number of errors: "+ str(errorcounter))
    print("")
    print("number of contigs with incorrect leng/coverage ratio: " + str(badcontigsfreq))
    print("The cumulative length of these contigs is: " + str(badcontigslen))

    
    

  

outputfile.close()

print("end of code")



