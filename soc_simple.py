import os, sys
import numpy as np
import collections
from copy import deepcopy

import jp_analysis as jp

            
def loadTOAhistogram_detailed(infile_name):
    #infile = open(infile_name)
    data = jp.jpickle(infile_name)
    zero_bin = np.argwhere(data['time_edges']>=0)[0][0]

    new_toa = np.zeros_like(data['toa'])
    
    for ipmt in range(data['toa'].shape[0]):
        if data['toa'][ipmt,:].sum() == 0 :
            continue
        # Find the peak
        pmt_max = data['toa'][ipmt,:].argmax()
        # Shift the peak until it reaches the position desired
        peak_diff = zero_bin - pmt_max
        # If positive, move forward.
        if peak_diff == 0:
            new_toa[ipmt, : ] = data['toa'][ipmt]
        elif peak_diff > 0:
            new_toa[ipmt, peak_diff:] = data['toa'][ipmt,:-peak_diff]
        else:
            new_toa[ipmt, :peak_diff] = data['toa'][ipmt,-peak_diff:]
    data['toa'] = new_toa
    #infile.close()
    return data
