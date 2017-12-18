#!/bin/bash

# Based on code from https://www.biostars.org/p/13452/

## Download NCBI's taxonomic data and GI (GenBank ID) taxonomic assignation.

CLEANUP_USELESS_FILES=0


## Variables
NCBI="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
TAXDUMP="taxdump"
#TAXID="gi_taxid_nucl.dmp"
NAMES="names.dmp"
NODES="nodes.dmp"
DMP=$(echo {citations,division,gencode,merged,delnodes}.dmp)
USELESS_FILES="${TAXDUMP}.gz ${DMP} gc.prt readme.txt"

if [ ! -s "$TAXDUMP.tar.gz" ]; then 
    ## Download taxdump
    #if [ $CLEANUP_USELESS_FILES -eq 1 ]; then rm -rf ${USELESS_FILES} "${NODES}" "${NAMES}"; fi
    curl "${NCBI}${TAXDUMP}.tar.gz"
fi
tar -xzf "${TAXDUMP}.tar.gz"
if [ $CLEANUP_USELESS_FILES -eq 1 ]; then rm -rf ${USELESS_FILES}; fi

## Limit search space to scientific names
head -n 3  $NAMES
grep "scientific name" "${NAMES}" > "${NAMES/.dmp/_reduced.dmp}" && \
    rm -f "${NAMES}" && \
    mv "${NAMES/.dmp/_reduced.dmp}" "${NAMES}"
perl -pi -e "s/\t\|\t/,/g" "${NAMES}"
cut -d "," -f 1-2 "${NAMES}" > "${NAMES}.tmp"; mv -f "${NAMES}.tmp" "${NAMES}"
head -n 3  $NAMES

head -n 3 $NODES
perl -pi -e "s/\t\|\t/,/g" "${NODES}"
cut -d "," -f 1-3 "${NODES}" > "${NODES}.tmp"; mv -f "${NODES}.tmp" "${NODES}"
head -n 3 $NODES

# if [ ! -s "$TAXID" ]; then 
#     ## Download gi_taxid_nucl
#     #rm -f "${TAXID}"
#     curl "${NCBI}${TAXID}.gz" && \
# 	gunzip "${TAXID}.gz"
# fi

exit 0
