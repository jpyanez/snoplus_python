import ROOT, rat
import numpy as np
import os, sys

hcvalue = 1.239841

def energyToWlen(energy):
    return  1E+3*hcvalue/(np.array(energy)*1E+6)

def getMChits(filename):
    ratreader = rat.dsreader(filename)
    data = []
    counter = 0
    for ds, run in ratreader:
        counter += 1
        for imc in range(0, ds.GetMCEVCount()):
            data.append(ds.GetMCEV(imc).GetMCHits().GetAllCount())
    data = np.array(data)
    ratreader.close()
    return data

# The AV has a radius of ~ 6m
# The belly plates have a dimension (in z) of 1.6 m
# The PSUP is situated at ~ 8.9 m
# The bellyplate projectio onto the PSUP from the detector center is given by
# Angle to the bellyplate edge (alpha): sin(alpha) = (1.6/2)/6
# The z position of the projection:
# tan(alpha) = z/8.9
# z = 8.9 * tan(arcsin(0.8/6))
# z = 1.197 ... rounding to 1.25m (some scattering allowed)
split_middle = 1250. #In mm
# # The neck has a radius of 0.785m. Doing the same exercise
# # Beta is now defined from the vertical to horizontal directions
# # Angle to the neck edge (beta): sin(beta) = (0.785/6)
# # The z position: sin(pi-beta) = z/8.9
# # z = 8.9 * sin(pi/2-arcsin(0.785/6))
# # z = 8.82 - rounding up at 8.6 (for statistics)
# neckz       = 8600.
pmtinfo = rat.utility().Get().GetPMTInfo()
def getHitRegions(filename):
    hits_list = []
    q_list    = []
    ratreader = rat.dsreader(filename)
    for ds, run in ratreader:
        for iEv in range(ds.GetMCEVCount()):
            mcev = ds.GetMCEV(iEv)
            hits = mcev.GetMCHits()
            hits_list.append([0,0,0])
            q_list.append([0,0,0])
            for iHit in range(hits.GetAllCount()):
                pmt = hits.GetAllPMT(iHit)
                pmtz = pmtinfo.GetPosition(pmt.GetID()).z()
                # print pmtinfo.GetType(pmt)
                # Top
                if  split_middle < pmtz:
                    hits_list[-1][0] += 1
                    q_list[-1][0]    += pmt.GetQHS()
                # Middle
                elif -split_middle < pmtz < split_middle:
                    hits_list[-1][1] += 1
                    q_list[-1][1]    += pmt.GetQHS()
                # Bottom
                else:
                    hits_list[-1][2] += 1
                    q_list[-1][2]    += pmt.GetQHS()
    ratreader.close()
    return np.array(hits_list), np.array(q_list)


def getPEcharge(filename):
    '''
    Obtain information on the PEs produced. Returns arrays with:
    - charge per event
    - charge per PMT
    - number of PEs per PMT
    - tracks that contribute to each PE
    '''
    charge = []
    tracks = []
    pexpmt = []
    charge_perev = []
    ratreader = rat.dsreader(filename)
    for ds, run in ratreader:
        mymc = ds.GetMC()
        # Per event
        charge_perev.append(0)
        for pmt_index in range(mymc.GetMCPMTCount()):
            # Per PMT
            mypmt  = mymc.GetMCPMT(pmt_index)
            pexpmt.append(mypmt.GetMCPECount())
            for pe_index in range(pexpmt[-1]):
                this_charge = mypmt.GetMCPE(pe_index).GetCharge()
                charge_perev[-1] += this_charge
                charge.append(this_charge)
                try:
                    tracks.append(np.size(mypmt.GetMCPE(pe_index).GetPhotonTrackID()))
                except:
                    # Photoelectron is noise
                    None

    return np.array(charge_perev), np.array(charge), np.array(pexpmt), np.array(tracks)


def getPEcharge_light(filename):
    '''
    Obtain information on the PEs produced. Returns arrays with:
    - charge per event
    - charge per PMT
    - number of PEs per PMT
    '''
    charge = []
    pexpmt = []
    charge_perev = []
    ratreader = rat.dsreader(filename)
    for ds, run in ratreader:
        mymc = ds.GetMC()
        # Per event
        charge_perev.append(0)
        for pmt_index in range(mymc.GetMCPMTCount()):
            # Per PMT
            mypmt  = mymc.GetMCPMT(pmt_index)
            pexpmt.append(mypmt.GetMCPECount())
            for pe_index in range(pexpmt[-1]):
                this_charge = mypmt.GetMCPE(pe_index).GetCharge()
                charge_perev[-1] += this_charge
                charge.append(this_charge)

    return np.array(charge_perev), np.array(charge), np.array(pexpmt)

# Doesnt work with my simple MC - PMTs not calibrated?
def getHitTimeResiduals_MC(filename):
    data = []
    ratreader = rat.dsreader(filename)
    for ds, run in ratreader:
        light_path      = rat.utility().GetLightPathCalculator()
        group_velocity  = rat.utility().GetGroupVelocity()
        pmt_info        = rat.utility().GetPMTInfo()
        event_position  = ds.GetMC().GetMCParticle(0).GetPosition()
        for iev in range(0, ds.GetEVCount()):
            calibrated_pmts = ds.GetEV(iev).GetCalPMTs()
            #print calibrated_pmts.GetCount()
            for ipmt in range(0, calibrated_pmts.GetCount()):
                pmt_cal = calibrated_pmts.GetPMT(ipmt)
                light_path.CalcByPosition(event_position, pmt_info.GetPosition(pmt_cal.GetID()))
                inner_av_distance = light_path.GetDistInInnerAV()
                av_distance       = light_path.GetDistInAV()
                water_distance    = light_path.GetDistInWater()
                transit_time      = group_velocity.CalcByDistance(inner_av_distance, av_distance, water_distance) 
                # Assumes 400nm photon

                data.append(pmt_cal.GetTime() - transit_time)

    data = np.array(data)
    ratreader.close()
    return data

def getHitTimes_MC(filename):
    data = []
    ratreader = rat.dsreader(filename)
    for ds, run in ratreader:
        for iev in range(0, ds.GetMCEVCount()):
            mc_hits = ds.GetMCEV(iev).GetMCHits()
            for imc in range(0, mc_hits.GetAllCount()):
                data.append(mc_hits.GetAllPMT(imc).GetTime())
    data = np.array(data)
    ratreader.close()
    return data    

def getPETimes_MC(filename):
    data = []
    ratreader = rat.dsreader(filename)
    for ds, run in ratreader:
        mc = ds.GetMC()
        for imcpmt in range(0, mc.GetMCPMTCount()):
            mcpmt = mc.GetMCPMT(imcpmt)
            for imcphotoelectron in range(0, mcpmt.GetMCPECount()):
                data.append(mcpmt.GetMCPE(imcphotoelectron).GetCreationTime())
    data = np.array(data)
    ratreader.close()
    return data

def getPMTWlen(filename):
    elist   = []
    trackid = []
    ratreader = rat.dsreader(filename)
    for ds, run in ratreader:
        mymc = ds.GetMC()
        for pmt_index in range(mymc.GetMCPMTCount()):
            mypmt = mymc.GetMCPMT(pmt_index)
            for photon_index in range(mypmt.GetMCPhotonCount()):
                elist.append(mypmt.GetMCPhoton(photon_index).GetEnergy())
                trackid.append(mypmt.GetMCPhoton(photon_index).GetPhotonTrackID())
    ratreader.close()
    return energyToWlen(elist), np.array(trackid)

def getPEinfo(filename):
    elist   = []
    trackid = []
    ratreader = rat.dsreader(filename)
    for ds, run in ratreader:
        mymc = ds.GetMC()
        for pmt_index in range(mymc.GetMCPMTCount()):
            mypmt = mymc.GetMCPMT(pmt_index)
            for pe_index in range(mypmt.GetMCPECount()):
                elist.append(mypmt.GetMCPE(pe_index).GetCharge())
                trackid.append(mypmt.GetMCPE(pe_index).GetPhotonTrackID())
    ratreader.close()
    return np.array(elist), np.array(trackid)

def getPhotonWlen(filename):
    elist = []
    ratreader = rat.dsreader(filename)
    for ds, run in ratreader:
        mymc = ds.GetMC()
        for itrack in range(mymc.GetMCTrackCount()):
            mytrack = mymc.GetMCTrack(itrack)
            # Take only the first step of each track
            mystep  = mymc.GetMCTrack(itrack).GetMCTrackStep(0)
            elist.append(mystep.GetKineticEnergy())

    ratreader.close()
    return energyToWlen(elist)


def readFile(data, basedir, one_dirname, fname, savename = None):
    if savename == None:
        savename = one_dirname
    if one_dirname in data.keys():
        print 'This directory is already inside:', one_dirname
        return
    else:
        data[savename] = {}
        
    file_name = os.path.join(basedir, one_dirname, fname)

    data[savename]['nhit'], data[savename]['qtot'] = \
        getHitRegions(file_name)
        
    # Make the sum at the end
    data[savename]['nhit'] = \
        np.vstack((data[savename]['nhit'].T, 
                   np.sum(data[savename]['nhit'], axis=1).T)).T
    data[savename]['qtot'] = \
        np.vstack((data[savename]['qtot'].T, 
                   np.sum(data[savename]['qtot'], axis=1).T)).T
    return




