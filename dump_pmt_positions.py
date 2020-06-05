import os, sys, pickle

import matplotlib.pyplot as plt
import collections
#%matplotlib inline
import os

import numpy as np
import ROOT, rat
import rat_misc


# Need to read at least one file to use the rat utility
infile='/home/jpyanez/snoplus/data/snoplus_data/physics_runs/Analysis_r0000109666_s001_p000.root'
ratreader = None
ratreader = rat_misc.openRat(infile, None, ratreader)

ds, run = ratreader.next()

def getBits(arg, loc, n):
    shifted = arg >> loc
    mask = (1 <<n ) -1
    value = shifted & mask
    return value

# Taken from BitManip (following Ben)
pmtid = 49 #pmt.GetID()
crate = getBits(pmtid, 9, 5)
card = getBits(pmtid, 5, 4)
channel = getBits(pmtid, 0, 5)
print crate, card, channel

du = rat.utility()

pmt = du.GetPMTInfo()

pmt.GetPanelNumber(400)


# Obtain and store pmt positions
maxpmts = 10000
pmt_positions = np.zeros([maxpmts,3])
pmt_directions = np.zeros([maxpmts,3])
pmt_type = np.zeros(maxpmts)
pmt_ccc = np.zeros([maxpmts,3], dtype=int)
for i in range(maxpmts):
    try:
        pmt_positions[i,:] = np.array(du.GetPMTInfo().GetPosition(i))
        pmt_directions[i,:] = np.array(du.GetPMTInfo().GetDirection(i))
        pmt_type[i] = du.GetPMTInfo().GetType(i)
        pmt_ccc[i,:] = [getBits(i, 9, 5),
                        getBits(i, 5, 4),
                        getBits(i, 0, 5)]
    except:
        print 'Maximum reached, exiting', i
        exit_i = i
        break
    
pmt_positions = pmt_positions[:exit_i,:]
pmt_directions = pmt_directions[:exit_i,:]
pmt_type = pmt_type[:exit_i]
pmt_ccc = pmt_ccc[:exit_i,:]

pmt_radii = np.sqrt(np.sum(pmt_positions**2,axis=1))
costheta = pmt_positions[:,2]/pmt_radii
phi      = np.arctan2(pmt_positions[:,1],
                      pmt_positions[:,0])
import pickle
pickle.dump({'xyz':pmt_positions, 'dir':pmt_directions, 'phi':phi, 'costheta':costheta, 'type':pmt_type, 'ccc':pmt_ccc},
            open('/home/jpyanez/snoplus/snoplus_python/pmt_positions.pckl','w'))
