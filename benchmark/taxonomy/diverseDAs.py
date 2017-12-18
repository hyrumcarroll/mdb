#!/usr/bin/env python3

# Description: Determine if a domain architecture (DA) has members (sequences) in multiple taxonomical domains
#              The optional DA subset file has the format of ^DA\s*.*
# Authors: Mileidy Gonzalez and Hyrum Carroll

import sys

# To set __debug__ to True, use python3 -d ...
def DEBUG(printStr):
    if __debug__: print( printStr, file = sys.stderr)

if len(sys.argv) < 3:
    print( "Usage: <DA membership filename> <taxonomy info filename> [<DA subset>]")
    raise SystemExit(1)

daMembershipFilename = sys.argv[1]
taxonomyInfoFilename = sys.argv[2]
daSubsetFilename = None
if len(sys.argv) >= 4:
    daSubsetFilename = sys.argv[3]

# try:
#     daMembershipFile = open( daMembershipFilename)
# except IOError:
#     print('ERROR: Cannot open', daMembershipFilename)
#     raise SystemExit(1)
# 
# daMembershipLines = daMembershipFile.readlines()
# daMembershipFile.close()
# 
# print( "Found", len(daMembershipLines), "lines")

# dictionary for speciesCode to domain
speciesCode2domain = {}
accession2domain = {}

# Using "with ... as", cleans up exceptions after the clause is executed
with open( taxonomyInfoFilename) as taxFO:
    for line in taxFO:
        # Expecting format to be: SPECIES_CODE\tKingdom\tTaxon_Node\tOfficial (scientific) name\tCommon name\tSynonym\tAccession\tsuperkingdom[...]
        fields = line.split('\t')
        speciesCode = fields[0]
        accession = fields[6]
        domain = fields[7]
        speciesCode2domain[ speciesCode] = domain
        accession2domain[ accession] = domain

DEBUG( "Found " + str(len( speciesCode2domain.keys())) + " species codes in " + taxonomyInfoFilename)
DEBUG( "Found " + str(len( accession2domain.keys())) + " accessions in " + taxonomyInfoFilename)

dasOfInterest = {}
if daSubsetFilename != None:
    with open( daSubsetFilename) as subsetFO:
        for line in subsetFO:
            if line.startswith('#'):
                continue
            dasOfInterest[ line.split()[0]] = 1
            
DEBUG( "Found " + str(len( dasOfInterest.keys())) + " DAs in " + daSubsetFilename)

print( "#DA","Breadth", "Archaea","Bacteria","Eukaryota", sep='\t')
with open( daMembershipFilename) as f:
    firstLine = True
    for line in f:
        if firstLine:
            firstLine = False
            continue
            
        # Expecting format to be: DA\tDomains\tCollapsed_Domains\tTaxa
        
        # print( line, end="")
        fields = line.split('\t')
        da = fields[0]
        if daSubsetFilename != None and not dasOfInterest.get(da):
            continue

        taxa = fields[3].split(',')
        # print( da, "has", len(taxa), "taxa");

        diverse = False
        domains = { "Archaea" : 0,
                    "Bacteria" : 0,
                    "Eukaryota" : 0
        }
        for taxon in taxa:
            # Expecting format: (pfam21|up)_ACCESSION_TAXONLABEL_SPECIESCODE
            ids = taxon.split('_')
            # DEBUG( "Found ids: " + str(ids))
            accession = ids[1]
            speciesCode = ids[3].strip()
            if accession in accession2domain:
                domains[ accession2domain[ accession]] = 1
            elif speciesCode in speciesCode2domain:
                domains[ speciesCode2domain[ speciesCode]] = 1
            else:
                DEBUG( "No taxonomy information found for either " + accession + " nor " + speciesCode)
                
            breadth = domains["Archaea"] + domains["Bacteria"] + domains["Eukaryota"]

            # have we already found all 3 domains?
            if breadth == 3:
                break
            
        print( da, breadth, domains["Archaea"], domains["Bacteria"], domains["Eukaryota"], sep='\t')
        #if breadth > 1:
        #    diverse = True
        #elif __debug__:
        #    print( "DEBUGGING:",da, breadth, domains["Archaea"], domains["Bacteria"], domains["Eukaryota"], sep='\t', file = sys.stderr)
            
            

exit(0)            
