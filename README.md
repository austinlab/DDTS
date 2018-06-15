# DDTS
DNAse-DTS  - Dnase I Direct to Sequencing

usage: DDTS.py -c [undigested.bam files] -e [digested.bam files] 

Options:
  -h                    show this help message and exit
  -c                    undigested background data files in BAM format (WIG format with -w)
  -e                    digested data files in BAM format (WIG format with -w) 
  -z [float]            cut-off value for z-score significance (default:3)
  -o [file]             output file name 
  -d [float]      		  fold change cut-off between undigested and digested (default: 2) 
  -l [int]              fseq feature length (default: 250) 
  -w                    perform analysis with already generated WIG files
  -W [int]              length of scan window (bp) (default: 250) 
  -s [int]              length of step increase (bp) (default: 25) 
