import ROOT, rat
import numpy as np
import os, sys
import matplotlib.pyplot as plt
import jp_mpl as jplot

def getPMToccupancy_data(ratreader, maxsize = np.inf):
    pmt_list = np.zeros(10000)
    counter = 0
    ev_counter = 0
    for ds, run in ratreader:
        ev_counter += 1

        for iEV in range(ds.GetEVCount()):
            event = ds.GetEV(iEV)
            pmts  = event.GetUncalPMTs()

            for iPMT in range(pmts.GetCount()):
                pmt_id = pmts.GetPMT(iPMT).GetID()
                pmt_list[pmt_id] += 1
                counter += 1
                if counter > maxsize:
                    return pmt_list, ev_counter
    return pmt_list, ev_counter


def getPMToccupancy_mc(ratreader, maxsize = np.inf):
    pmt_list = np.zeros(10000)
    counter = 0
    ev_counter = 0
    for ds, run in ratreader:
        ev_counter += 1
        for iEV in range(ds.GetMCEVCount()):
            event = ds.GetMCEV(iEV)
            hits = event.GetMCHits()      
            for iHit in range(hits.GetCount()):
                pmt_id = hits.GetPMT(iHit).GetID()
                pmt_list[pmt_id] += 1
                counter += 1
                if counter > maxsize:
                    return pmt_list
    return pmt_list, ev_counter

def plotOccupancy(pmt_list, pmt_info, hit_lim = [-1,-1]):
    pmt_list = pmt_list[:pmt_info['phi'].size]
    print pmt_list.shape, pmt_info['phi'].size
    pmt_radii = np.sqrt(np.sum(pmt_info['xyz']**2,axis=1))
    pmt_bool  = (pmt_radii < 8500.)*(pmt_list>0.)

    print hit_lim
    if hit_lim[0] < 0:
        aux = pmt_list.mean()-2.5*pmt_list.std()
        if aux < 0:
            hit_lim[0] = 0.
        else:
            hit_lim[0] = aux 
    if hit_lim[1] < 0:
        hit_lim[1] = pmt_list.mean()+1.5*pmt_list.std()
    print hit_lim
    # Print the histogram of the occupancy
    f1 = plt.figure()
    xaxis = np.linspace(hit_lim[0], hit_lim[1], 201)
    b, x = np.histogram(pmt_list[pmt_bool], xaxis)
    jplot.unfilledBar(x,b)
    plt.xlabel('Nhits')
    plt.ylabel('Number of PMTs')
    plt.show()
    
    f2 = plt.figure(figsize=(12,7))
    ax = f2.add_subplot(111)
    plt.scatter(pmt_info['phi'][pmt_bool]/np.pi, pmt_info['costheta'][pmt_bool], 
                c=pmt_list[pmt_bool], 
                cmap='jet',vmin=hit_lim[0], vmax=hit_lim[1],
                marker='o', s=18,lw = 0)
    plt.colorbar()
    plt.xlim([-1,1])
    plt.ylim([-1,1])
    plt.xlabel('Phi / pi')
    plt.ylabel('cos(theta)')
    plt.show()
