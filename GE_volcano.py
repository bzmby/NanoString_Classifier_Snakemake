#!/usr/bin/env python
#Author: Behzad Moumbeini
# this script has these  outputs:
# 1- CSV files with the list of DEGs 
# 2- volcano plot
# usage: python3 GE_volcano.py




import subprocess
import os

def GE_volcano(r_script):
	subprocess.call (["/usr/bin/Rscript", "--vanilla", r_script])
	return

if __name__ == '__main__':

    r_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'GE_volcano.r') 
    print(GE_volcano(r_script))

