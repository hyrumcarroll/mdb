#!/bin/bash

#$ -N hardwareSpecs
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -j n
#$ -P unified
#$ -l h_rt=10
#$ -l h_vmem=8.1G
#$ -l mem_free=6G
# use reserve_mem only for large jobs b/c it will block other users from getting access to that node
#$ -l reserve_mem=8G
#$ -m n

hostname
cat /proc/cpuinfo
cat /proc/meminfo
free
