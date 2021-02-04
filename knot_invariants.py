#!/bin/python3

# Usage: knot_invariants data_file output_file

import sys
import numpy as np
from pyknotid.spacecurves import Knot
import pyknotid
from sklearn.decomposition import PCA
from sympy import symbols

t=symbols('t')

if(len(sys.argv)==3):
   # Read data from specified file
   points=np.genfromtxt(sys.argv[1],delimiter=',')

   # Use PCA to extract the three most significant components
   pca = PCA(n_components=3)

   # Render as a knot
   k=Knot(pca.fit_transform(points))

   # Write invariants to output file
   with open(sys.argv[2],'wt') as fp:
      fp.write('{}\n'.format(pyknotid.invariants.alexander(k.representation(),variable=t)))
   
   
