# DDTS
DNAse-DTS  - Dnase I Direct to Sequencing 

usage: DDTS.py -c [undigested.bam files] -e [digested.bam files] 

Options:<br>
  -h                    show this help message and exit <br>
  -c                    undigested background data files in BAM format (WIG format with -w) <br>
  -e                    digested data files in BAM format (WIG format with -w)  <br>
  -z [float]            cut-off value for z-score significance (default:3) <br>
  -o [file]             output file name <br>
  -d [float]      		  fold change cut-off between undigested and digested (default: 2)  <br>
  -l [int]              fseq feature length (default: 250)  <br>
  -w                    perform analysis with already generated WIG files <br>
  -W [int]              length of scan window (bp) (default: 250) <br>
  -s [int]              length of step increase (bp) (default: 25) <br>
