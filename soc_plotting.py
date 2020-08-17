import os, sys, pickle

try:
    import matplotlib.pyplot as plt
    import jp_mpl as jplot
except:
    print 'Not using matplotlib'
import numpy as np
import collections

import jp_analysis as jp


try:
    pmt_info = pickle.load(open('/home/users/jpyanez/snoplus/snoplus_python/pmt_positions.pckl'))
except:
    print('Pickle failed in py3, trying py2')
    try:
        pmt_info = pickle.load(open('/home/jpyanez/snoplus/snoplus_python/pmt_positions.pckl'))
    except:
        print('Pickle failed in py2 as well')


def mynorm(c,axis=-1):
    if len(c.shape) > 1 and axis>=0:
        return np.apply_along_axis(np.linalg.norm, axis, c)
    else:
        return np.linalg.norm(c)

#print pmt_info['xyz'], mynorm(pmt_info['xyz'],axis=1)      
pmt_info['xyz_norm'] = pmt_info['xyz']/mynorm(pmt_info['xyz'],axis=1).reshape(pmt_info['xyz'].shape[0],1)
pmt_info['r'] = mynorm(pmt_info['xyz'], axis=1)
pmtbool = pmt_info['type']==1
psup_r = pmt_info['r'][pmtbool].mean()#/10.
#pmt_info.keys()


cs= ['C'+"%i" % x for x in range(0,10)]
for i in range(50):
        cs.append(np.random.rand(3,))


#def rebinHisto(hin, factor):
#    hout = hin.cumsum()[(factor-1)::factor]
        
# Make figure
def plotTOA( data = None,
             position = np.array([0,0,0]),
             theta_bins = 6,
             time_window = 120.,
             time_rebin=1,
             av_reference = False,
             color = None,
             label=None,
             newfig = True):

    time_bin = data['time_edges'][1]-data['time_edges'][0]
    # Time rebinning
    time_nbins = int(np.ceil(time_window/time_bin))
    nbins = 5
    
    theta_edges = np.linspace(0, 180, theta_bins+1)

    if av_reference:
        a = psup_r
        b = mynorm(position)
        c = mynorm(position-pmt_info['xyz'], axis=1)
        mod_ct = (a**2 + b**2 - c**2)/(2*a*b)
        mod_theta = np.rad2deg(np.arccos(mod_ct))
    else:
        mod_pos   = pmt_info['xyz']-position
        mod_ct    = mod_pos[:,2]/pmt_info['r']
        mod_theta = np.rad2deg(np.arccos(mod_ct))
        
    if newfig:
        fig = plt.figure(figsize=(12,7))

    #Rebin?
    time_edges = data['time_edges']
    if time_rebin > 1:
        new_time = time_edges[::time_rebin][:-1]
        time_edges = new_time

    
    for i in range(len(theta_edges)-1):
        mybool = pmtbool*(mod_theta>theta_edges[i])*(mod_theta<theta_edges[i+1])
        
        norm=1.

        
        hits = data['toa'][:len(mybool),:][mybool].sum(axis=0)*1.
        if time_rebin > 1:
            hout = hits.cumsum()[(time_rebin-1)::time_rebin]
            hdiff = hout[1:] - hout[:-1]
            hits = hdiff
        #print hits
        imax = hits.argmax()
        #print imax
        
        hits /= hits[imax-5:imax+5].sum()
        #print hits
        hits /= (2**i)
        
        if color == None:
            this_color= cs[i]
        else:
            this_color = color
        if label == None:
                mylabel = "%i" % theta_edges[i] + ' - ' + "%i" % theta_edges[i+1]
        else:
                mylabel = label
        jplot.unfilledBar(time_edges - time_edges[imax], #+t0[i], ##-data_t0,
                          hits*norm,#*1./real_data['toa'].sum(),
                          color = this_color,
                          label = label)
        
    plt.legend(bbox_to_anchor=(1.15, 1.),
               loc='upper right')
    plt.yscale('log')
    plt.xlim(-10, 100)
    #plt.ylim(10,400000)
    plt.ylabel('Hits')
    plt.xlabel('Corrected time residual (ns)')
    #plt.legend(loc=0)

# Make figure
def plotTOA2D( data = None,
               position = np.array([0,0,0]),
               time_window = 120.,
               theta_bins =30,
               time_rebin = 1,
               av_reference = False,
               renorm = True,
               plot=True,
               plot_time = 15):
    

    time_bin = data['time_edges'][1]-data['time_edges'][0]
    time_nbins = int(np.ceil(time_window/time_bin))
    nbins = 5

    
    theta_edges = np.linspace(0, 180, theta_bins+1)
    h = np.zeros([time_nbins+nbins, theta_bins])
    herr = np.zeros_like(h)

    if av_reference:
        # Theta wrt center of the av
        # Using law of cosines, a and b are psup_r
        #a = psup_r
        #b = mynorm(position)
        #c = mynorm(position-pmt_info['xyz'], axis=1)
        #mod_ct = (a**2 + b**2 - c**2)/(2*a*b)
        #mod_theta = np.rad2deg(np.arccos(mod_ct))

        mod_ct = np.dot(position, pmt_info['xyz'].T)/(mynorm(position)*mynorm(pmt_info['xyz'],axis=1))
        mod_theta =  np.rad2deg(np.arccos(mod_ct))
        
    else:
        # This is the theta angle wrt the source
        #mod_pos   = pmt_info['xyz']-position
        #mod_ct    = mod_pos[:,2]/pmt_info['r']
        #mod_theta = np.rad2deg(np.arccos(mod_ct))

        pmtpos = pmt_info['xyz']-position
        mod_ct = np.dot(position, pmtpos.T)/(mynorm(position)*mynorm(pmtpos,axis=1))
        mod_theta = np.rad2deg(np.arccos(mod_ct))
        
    norm_list = np.zeros(len(theta_edges)-1)
    for i in range(len(theta_edges)-1):
        mybool = pmtbool*(mod_theta>theta_edges[i])*(mod_theta<theta_edges[i+1])
        norm=1.
        hits = data['toa'][:len(mybool),:][mybool].sum(axis=0)*1.
        imax = hits.argmax()
        if imax == 0:
            continue

        if renorm:
            norm = 1./hits[imax+15:imax+10+time_nbins].sum()
        else:
            norm = 1.
        #hits *= norm
        h[:,i]    = hits[(imax-nbins):(imax+time_nbins)]
        h[:,i] *= norm
        norm_list[i] = norm

    # Rebin?
    time_edges = data['time_edges'][:nbins+time_nbins+1]+45
    if time_rebin > 1:
        hout = np.vstack((np.zeros([1,h.shape[1]]),h.cumsum(axis=0)[(time_rebin-1)::time_rebin,:]))
        overflow = h.sum(axis=0)-hout[-1,:]
        hout[-1,:] += overflow

        hdiff = hout[1:,:]-hout[:-1,:]

        new_time = time_edges[::time_rebin]
        if new_time.shape[0] == hdiff.shape[0]:
            new_time = np.concatenate((new_time + [time_edges.max()]))
            
        time_edges = new_time
        h = hdiff
        #print h.shape
    # Get the errors
    #print h.shape, norm_list.shape
    #print h.shape, (h*norm_list).shape
    hcounts = h/norm_list
    herr = np.sqrt(hcounts)/hcounts

    fig =None
    if plot:
        fig = plt.figure(figsize=(12,7))
        tini = np.where(time_edges>plot_time)[0][0]
        refh = np.log10(h[tini:,:])
        plt.pcolor(time_edges,
                   theta_edges, np.log10(h.T),
                   vmax = refh.max(),
                   vmin = refh.max()-2.)
        plt.colorbar()
        
        plt.xlim(plot_time, 120)

    return time_edges, theta_edges, h, herr, fig
    

#def plotSimpleTOA():
#        
