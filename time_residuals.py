import ROOT, rat
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
print 'This is AIRPLANE MODE - be careful!'



def getTimeResiduals(infile_list = [],
                     outfile = '',
                     bin_width = 0.5, # in ns
                     time_window = 300., # in ns
                     calibration = True,
                     remove_bad_pmts = False,
                     min_intimehits = 10.,
                     min_nhits = 15,
                     radius_range = [7020, 10000],
                     udotr_range  = [-1., -0.8],
                     emin = 3.,
                     dcmask = -1,
                     start_time = -120.,
                     max_ev = -1): # in ns - trigger is at ~340

    print 'Radius range', radius_range
    print 'UdotR range', udotr_range

    # Size of the array
    time_edges = np.arange(start_time, start_time+time_window+bin_width, bin_width)
    nbins = time_edges.size-1
    all_toa = np.zeros([9800, nbins])

    light_path      = rat.utility().GetLightPathCalculator()
    group_velocity  = rat.utility().GetGroupVelocity()
    pmt_info        = rat.utility().GetPMTInfo()
    ev_counter = 0
    ev_earlyhit10 = 0
    ev_earlyhit20 = 0
    qpdt_flagged = 0
    break_all  = False
    # Loop over all the files in infile_list
    for iFile, one_file in enumerate(infile_list):
        if break_all:
            break
        print iFile,'/', len(infile_list), ' - reading ', one_file
        if True: #try:


            reader = rat.dsreader(one_file)
            this_file_ev = 0
            this_file_qpdt = 0
            for ds, run in reader:
              
                for iEv in range(ds.GetEVCount()):
                    if break_all:
                        break
                
                    event = ds.GetEV(iEv)
                    pmts = event.GetCalPMTs()

                    # Select which events to run on
                    if len(event.GetFitNames()) < 1:
                        continue
                    if event.GetNhitsCleaned() < min_nhits:
                        continue
                    if event.GetInTimeHits100() < min_intimehits:
                        continue

                    event_vtx = event.GetFitResult(event.GetFitNames()[0]).GetVertex(0)
                    if not (event_vtx.ValidDirection()*event_vtx.ValidEnergy()):
                        continue
                    
                    use_event = False
                    # Selecting only FECD tagged
                    if calibration:
                        for iFECD in range(pmts.GetFECDCount()):
                            fecd = pmts.GetFECDPMT(iFECD)
                            if fecd.GetID() != 9188:
                                continue
                            use_event = True
                    else:
                        # For Non calibration runs I'm doing the following cuts
                        # to get the PMT box only
                        # E > 3. MeV
                        # UdotR < -0.8
                        # posR > 7020.
                        if event_vtx.GetEnergy() < emin:
                            continue

                        r = np.array(event_vtx.GetPosition())
                        radius = np.linalg.norm(r)
                        if radius < radius_range[0] or radius > radius_range[1]:
                            continue
                        u = np.array(event_vtx.GetDirection())
                        udotr = np.dot(u,r)/radius
                        if udotr < udotr_range[0] or udotr > udotr_range[1]:
                            continue

                        use_event = True

                        
                    if use_event:

                        this_early10 = 0
                        this_early20 = 0
                        
                        # Data cleaning
                        if dcmask>0:
                            dcapplied = event.GetDataCleaningFlags().GetApplied(0).GetULong64_t(0)
                            dcflagged = event.GetDataCleaningFlags().GetFlags(0).GetULong64_t(0)
                            dcpass = ((dcapplied & dcmask) & dcflagged) == (dcapplied & dcmask)
                            if not dcpass:
                                continue

                        
                        # Fill the PMT info
                        event_pos = event_vtx.GetPosition()
                            
                        for iPMT in range(pmts.GetNormalCount()):
                            one_pmt = pmts.GetNormalPMT(iPMT)
                            pmtid = int(one_pmt.GetID())
                            if remove_bad_pmts:
                                if one_pmt.GetCrate() == 8 and one_pmt.GetCard() == 9:
                                    continue

                            light_path.CalcByPosition(event_pos, pmt_info.GetPosition(pmtid))
                            inner_av_distance = light_path.GetDistInInnerAV()
                            av_distance       = light_path.GetDistInAV()
                            water_distance    = light_path.GetDistInWater()
                            transit_time      = group_velocity.CalcByDistance(inner_av_distance, av_distance, water_distance)

                            time_residual = one_pmt.GetTime() - event_vtx.GetTime() - transit_time
                            tbin = np.digitize(time_residual, time_edges)-1

                            # Getting the rejection probability
                            if time_residual < -10 and time_residual > -100:
                                this_early10 = 1
                            if time_residual < -20 and time_residual > -75:
                                this_early20 = 1

                            # Filling the histogram
                            if tbin >= 0 and tbin < nbins:
                                all_toa[pmtid, tbin] += 1

                        qpdt = event.GetClassifierResult('QPDT:waterFitter').GetClassification('NHits_early')
                        
                        
                        if qpdt > 0:
                            qpdt_flagged += 1
                            this_file_qpdt += 1

                        # Summing up the rejection probability
                        ev_earlyhit10 += this_early10
                        ev_earlyhit20 += this_early20
                        ev_counter += 1
                        this_file_ev += 1
                        if ev_counter%1000 ==0:
                            print 'Counter ', ev_counter
                        if ev_counter == max_ev:
                            break_all = True
                        break
                    if break_all:
                        break
                if break_all:
                    break
            if this_file_ev > 0:
                print '\n\n *** File: ', one_file.split('/')[-1] 
                print 'QPDT this file rej', this_file_qpdt*1./this_file_ev
                print 'QPDT global rejection', qpdt_flagged*1./ev_counter, '\n'
            reader.close()

            # After I'm done with one file, save
            if len(outfile)>0:
                if os.path.isfile(outfile):
                    in_data = pickle.load(open(outfile))
                    in_data['toa'] += all_toa
                    in_data['events'] += ev_counter
                    in_data['ev_early10'] += ev_earlyhit10
                    in_data['ev_early20'] += ev_earlyhit20
                    in_data['qpdt_flagged'] += qpdt_flagged
                    pickle.dump(in_data, open(outfile, 'wb'), protocol=2)
                else:
                    pickle.dump({'time_edges':time_edges,
                                 'toa':all_toa,
                                 'events':ev_counter,
                                 'ev_early10':ev_earlyhit10,
                                 'ev_early20':ev_earlyhit20,
                                 'qpdt_flagged':qpdt_flagged},
                                open(outfile, 'wb'), protocol=2)
            
        else: #except:
            print 'Could not open file! skipping it for now...'
            try: reader.close()
            except: print 'No reader closed'



    return time_edges, all_toa, ev_counter
    print 'Done with all the files!'
