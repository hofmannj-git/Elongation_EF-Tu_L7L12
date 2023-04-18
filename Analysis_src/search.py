#import pandas as pd
import sys
import os
sys.path.insert(0,'./')
import h5py as h5
import numpy as np
#import h5py
import freud
from freud import box, locality
#from shape import *
from tqdm import tqdm

# get files
files = np.array([x for x in os.listdir('.') if x.endswith('.h5')])
times = np.array([int(x.split('.')[3][1:]) for x in files])
files = files[np.argsort(times)]
times = times[np.argsort(times)]


# Wrap box, extract first time step from h5 file
b_l = 15.52845

# Create indices for ribosomes vs EF-Tu Complexes
n_rib = 1266
n_tc = 6090
n_crowder = 205610
rib_id = [5*i for i in range(n_rib)]
tc_id = [rib_id[-1]+5+2*i for i in range(n_tc)]

# create binding boolean array (n_tc x n_rib)
max_bond_length = (1+.1433) * (0.5 + 0.4521/2)
box = freud.box.Box.cube(b_l)

# Loop through snapshots and generate binding matrices

for x, f in tqdm(enumerate(files)):
    coords = h5.File(f)['pos']
    bound_idx = np.zeros((len(coords),n_rib,n_tc),dtype=int)
    for i,c in enumerate(coords):
    
        # slice out coordinates of ribosome and TC bodies
        rib_c = c[rib_id]
        tc_c = c[tc_id]
        
        # find periodic distances via freud 
        d = box.compute_all_distances(rib_c,tc_c)
        bound_idx[i,np.where(d <= max_bond_length)[0],np.where(d <= max_bond_length)[1]] = 1
        
    np.save('search/binding%i.npy'%times[x],bound_idx)


