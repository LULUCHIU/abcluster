#!/home/Tina/miniconda3/bin/python

"""
Created On Jan 10, 2022
@author: Tina

"""

import pandas as pd
import numpy as np
import sys
import os
import os.path
from abnumber import Chain

inseq=sys.argv[1]
outpath=sys.argv[2]

seq = inseq
chain = Chain(seq, scheme='kabat')
sys.stdout=open(outpath+"/test.txt","w")
#f=open(outpath+"/test.txt","w")
chain.renumber(cdr_definition='kabat').print_tall(1)
sys.stdout.close()
#f.close()

