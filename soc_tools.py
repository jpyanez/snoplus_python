import rat
import os, sys
import cPickle as pickle

import numpy as np
import collections

from copy import deepcopy

# Only extracting information, always run in airplane mode
# Airplane mode
db = rat.RAT.DB.Get()
db.SetAirplaneModeStatus(True)
db.SetDefaultPlaneLockStatus(False)
print('This is AIRPLANE MODE - be careful!')



def getTOAhistogram_detailed(infile_list = [],
                             outfile = '',
                             bin_width = 1, # Do not move this
                             time_window = 200.,
                             qhs_threshold = 0.,
                             start_time = -50.):

    # Size of the array
    time_edges = np.arange(start_time, start_time+time_window+bin_width, bin_width)
    all_toa = np.zeros([9800, time_edges.size-1])
    t0 = np.zeros(9800)
    
    # Loop over all the files in infile_list
    first_file = True
    for one_file in infile_list:
        print('\nReading ', one_file)
        #try:
        reader = rat.socreader(one_file)
        soc, run = reader.next()
        print('manip position', np.array(soc.calib.GetPos()))
        toa = np.zeros_like(all_toa)
        
        for one_pmt in soc.GetSOCPMTIDs():
            time_array = np.array(soc.GetSOCPMT(one_pmt).GetTimes())
            qhs_array  = np.array(soc.GetSOCPMT(one_pmt).GetQHSs())
            qhs_bool   = qhs_array > qhs_threshold
                
            # Set the offset if the t0 is not set and it is the first file
            if first_file and t0[one_pmt] == 0:
                t0[one_pmt] = np.ceil(time_array[qhs_bool].mean())                       
                
            counts, x = np.histogram(time_array[qhs_bool]-t0[one_pmt], time_edges)
            toa[one_pmt,:] += counts
            if one_pmt % 300 == 0:
                print(one_pmt), time_array.mean(), t0[one_pmt]
                        
        all_toa += toa
        reader.close()
        first_file = False
        #except:
        #    print('Could not open file! skipping it for now...')
        #    try: reader.close()
        #    except: print('No reader closed')
            
        # The pickle dump is done inside the loop in case something goes wrong later on
        myoutfile = open(outfile, 'wb')
        pickle.dump({'time_edges':time_edges,
                     'toa':all_toa,
                     't0':t0},
                    myoutfile, protocol=2)
        myoutfile.close()
            
    print('Done with all the files!')
            
            
def getTOAhistogram(infile_list = [],
                    outfile = '',
                    bin_width = 1., # in ns
                    time_window = 400., # in ns
                    qhs_threshold = 0.,
                    start_time = 300.): # in ns - trigger is at ~340

    # Size of the array
    time_edges = np.arange(start_time, start_time+time_window+bin_width, bin_width)
    all_toa = np.zeros([9800, time_edges.size-1])

    # Loop over all the files in infile_list
    for one_file in infile_list:
        print('\nReading ', one_file)
        try:
            reader = rat.socreader(one_file)
            soc, run = reader.next()
            print('manip position', np.array(soc.calib.GetPos()))
            toa = np.zeros_like(all_toa)

            for one_pmt in soc.GetSOCPMTIDs():
                if one_pmt % 300 == 0:
                    print(one_pmt)
                time_array = np.array(soc.GetSOCPMT(one_pmt).GetTimes())
                qhs_array  = np.array(soc.GetSOCPMT(one_pmt).GetQHSs())
                qhs_bool   = qhs_array > qhs_threshold
                counts, x = np.histogram(time_array[qhs_bool], time_edges)
                toa[one_pmt,:] += counts
                if one_pmt % 300 == 0:
                    print(one_pmt), time_array.mean()
                    
            all_toa += toa
            reader.close()
        except:
            print('Could not open file! skipping it for now...')
            try: reader.close()
            except: print('No reader closed')

        # The pickle dump is done inside the loop in case something goes wrong later on
        pickle.dump({'time_edges':time_edges,
                     'toa':all_toa},
                    open(outfile, 'wb'), protocol=2)
        
    print('Done with all the files!')


def loadTOAhistogram_detailed(infile_name):
    infile = open(infile_name)
    data = pickle.load(infile)
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
    infile.close()
    return data
