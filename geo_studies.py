import ROOT, rat
import numpy as np
import collections, os, sys
import matplotlib.pyplot as plt
sys.path.append('/home/jp/projects/python_tools')
import jp_mpl as jplot
from copy import deepcopy


color_list = ['r','g','m','y', 'c','k']
generic_labels = ['set_1','set_2','set_3','set_4']
figsize = (14,9)
figsize_small = (9,6)

def doRead(ratreader = None, max_photons = 100000, test=False):
    costheta_dir     = np.zeros(max_photons)
    step_r           = np.zeros(max_photons)
    n_steps          = np.zeros(max_photons)
    end_position     = np.zeros([max_photons, 3])
    end_volume       = ['volume']*max_photons
    end_process      = ['process']*max_photons
    counter = 0

    for ds, run in ratreader:
        mymc = ds.GetMC()
        
        if counter >= max_photons:
            break
        
        # Loop over all the photons
        for itrack in range(1, mymc.GetMCTrackCount()):
            mytrack = mymc.GetMCTrack(itrack)

            
            # End here to study the info in a track
            if test:
                #print mymc.GetMCTrackCount()+1
                return mytrack
            
            max_steps = mymc.GetMCTrack(itrack).GetMCTrackStepCount()
            n_steps[counter] = max_steps
            
            # The step 0 is always the starting point. Start at one.
            starting_point = mymc.GetMCTrack(itrack).GetMCTrackStep(0).GetPosition()
            first_step     = mymc.GetMCTrack(itrack).GetMCTrackStep(1).GetPosition()
            step_r[counter]    = np.linalg.norm(first_step)

            # Get the initial direction of the photon (I'm intersted in theta_dir = [8, 82] degrees)
            costheta_dir[counter] = first_step[2]/step_r[counter]
            
            # Get the last step
            last_step = mymc.GetMCTrack(itrack).GetMCTrackStep(max_steps-1)

            # Get the position, volume and process of the last step
            end_volume[counter] = last_step.GetEndVolume()
            end_position[counter,:] = last_step.GetPosition()
            end_process[counter] = last_step.GetProcess()
            counter += 1
            
            if counter >= max_photons:
                break
                
    ratreader.close()
    
    data = {'step_r':step_r[:counter],
            'costheta_dir':costheta_dir[:counter],
            'n_steps':n_steps[:counter],
            'end_position':end_position[:counter,:],
            'end_volume':np.array(end_volume[:counter]),
            'end_process':np.array(end_process[:counter])}

    data['end_r'] = np.linalg.norm(data['end_position'],axis=1)
    data['end_costheta'] = data['end_position'][:,2]/data['end_r']
        
    return data

def compareTracking(set_list = [], costheta_range = [-1, 1], outdir = '', pltlabels = []):
    reference = set_list[0]

    if len(pltlabels) == 0:
        pltlabels = generic_labels[:len(set_list)]

    set_bool = []
    for one_set in set_list:
        set_bool.append((one_set['costheta_dir'] >= costheta_range[0])*\
                        (one_set['costheta_dir'] <= costheta_range[1]))

    # Do the reference set first
    ref_gen = len(reference['end_volume'][set_bool[0]])
    labels_ref, values_ref = zip(*collections.Counter(reference['end_volume'][set_bool[0]]).items())
    labels_ref = np.array(labels_ref)
    values_ref = np.array(values_ref)
    sort_ref = np.argsort(labels_ref)
    labels_ref = labels_ref[sort_ref]
    values_ref = values_ref[sort_ref]


    values_list = []

    for i in range(1, len(set_list)):
        set_gen = len(set_list[i]['end_volume'][set_bool[i]])
        #print 'Difference in photons produced', ref_gen - set_gen, np.sqrt(ref_gen)

        labels, values = zip(*collections.Counter(set_list[i]['end_volume'][set_bool[i]]).items())
        # Check if there are differences in the fields
        diff_labels = list(set(labels_ref) - set(labels))
        if len(diff_labels)>0:
            #print 'Labels missing ', diff_labels

            labels = list(labels) + diff_labels
            values = list(values) + [0]*len(diff_labels)

        labels = np.array(labels)
        values = np.array(values)
        onesort = np.argsort(labels)
        labels = labels[onesort]
        values = values[onesort]

        values_list.append(values)

    fig1 = plt.figure(figsize=(7,7))
    ax1 = fig1.add_subplot(211)
    plt.title('cos(theta) = [' + "%.2f" % costheta_range[0] + ', '+"%.2f" % costheta_range[1] + ']')
    indices = np.arange(len(labels_ref)+1)
    jplot.unfilledBar(indices, values_ref, label =pltlabels[0])
    plt.xticks(indices + 0.5, labels_ref)
    plt.yscale('log')
    plt.xlim([0, indices[-1]])
    plt.ylabel('Photons')

    for i in range(len(values_list)):
        jplot.unfilledBar(indices, values_list[i], color=color_list[i],label=pltlabels[i+1])
    plt.legend(loc=0)

    ax2 = fig1.add_subplot(212,sharex=ax1)
    #fig2 = plt.figure()
    #plt.title('cos(theta) = [' + "%.2f" % costheta_range[0] + ', '+"%.2f" % costheta_range[1] + ']')
    plt.ylabel('Ratio')
    plt.plot([0, indices[-1]], [1,1], '--k')
    ref_error = np.sqrt(values_ref)/values_ref
    for i in range(len(values_list)):
        ratio = values_list[i]*1./values_ref
        error = np.sqrt(2)*ratio*ref_error
        jplot.unfilledBar(indices, ratio, color=color_list[i])
        jplot.errorMarkVert(indices*1., ratio, yerror=error, color=color_list[i])
        print ratio
    plt.xticks(indices + 0.5, labels_ref)
    plt.ylim([0.9, 1.1])
    #plt.ylim([0.65, 1.35])
    plt.xlim([0, indices[-1]])
    plt.subplots_adjust(hspace=0)

    fig3 = plt.figure(figsize=(7,7))
    ax3 = fig3.add_subplot(211)

    plt.title('cos(theta) = [' + "%.2f" % costheta_range[0] + ', '+"%.2f" % costheta_range[1] + ']')
    plt.ylabel('Photons')
    plt.xlabel('Steps taken')
    myx = np.arange(0.5,17,1)
    #myx[-1] = np.inf

    ref_steps, x = np.histogram(reference['n_steps'][set_bool[0]]-1, myx)
    jplot.unfilledBar(myx, ref_steps, label=pltlabels[0])
    for i in range(1,len(set_list)):
        n_steps, x = np.histogram(set_list[i]['n_steps'][set_bool[i]]-1, myx)
        jplot.unfilledBar(myx, n_steps, color = color_list[i-1],label=pltlabels[i])
    plt.yscale('log')
    plt.legend(loc=0)

    ax4 = fig3.add_subplot(212,sharex=ax3)
    plt.ylabel('Ratio')
    plt.xlabel('Steps taken')

    for i in range(1, len(set_list)):
        n_steps, x = np.histogram(set_list[i]['n_steps'][set_bool[i]]-1, myx)
        ratio = n_steps*1./ref_steps
        error = np.sqrt(2)*ratio*np.sqrt(ref_steps)/ref_steps
        jplot.unfilledBar(myx, ratio,  color=color_list[i-1])
        jplot.errorMarkVert(myx, ratio, yerror=error, color=color_list[i-1])
        print ratio
    plt.ylim([0.9, 1.1])
    #plt.ylim([0.65, 1.35])
    plt.subplots_adjust(hspace=0)

    if len(outdir) > 0:
        theta_str = "%.2f" % costheta_range[0] + '_'+"%.2f" % costheta_range[1]
        theta_str = theta_str.replace('.','')
        fig1.savefig(os.path.join(outdir, 'EndPhotons_'+theta_str+'.pdf'), dpi=200)
        #fig2.savefig(os.path.join(outdir, 'EndPhotonsRatio_'+theta_str+'.pdf'), dpi=200)
        fig3.savefig(os.path.join(outdir, 'Steps_'+theta_str+'.pdf'), dpi=200)
        #fig4.savefig(os.path.join(outdir, 'StepsRatio_'+theta_str+'.pdf'), dpi=200)

    return

def volumeAnalysis(data, volume = [], radius_range = [0, np.inf], 
                   extra_conditions = {},
                   theta_bins = 41, phi_bins = 61,
                   plot_mode = 'histogram',
                   plot_nstd = 2.,
                   radius_plot = True):
    if type(volume) != list:
        volume = [volume]
    radius_range = np.array(radius_range)

    mybool = np.array([False]*len(data['end_volume']))

    for one_volume in volume:
        mybool = mybool | (data['end_volume'] == one_volume)


    # Additional conditions
    for one_key in extra_conditions.keys():
        mybool *= data[one_key] >= extra_conditions[one_key].min()
        mybool *= data[one_key] <= extra_conditions[one_key].max()

    positions = data['end_position'][mybool, :]

    # Converting to spherical coordinates
    r     = np.linalg.norm(positions, axis=1)
    rbool = (r>radius_range.min())*(r<radius_range.max())

    phi   = np.arctan2(positions[rbool,1], positions[rbool,0])
    theta = -1*(np.arccos(positions[rbool,2]/r[rbool]) - np.pi/2.)    

    # Plotting the events as function of r
    raxis = np.linspace(0, 8000., 201)
    raxis_eff = deepcopy(raxis)
    raxis_eff[-1] = np.inf

    if radius_plot:
        b, x = np.histogram(r, raxis_eff)
        
        fig1 = plt.figure(figsize=figsize_small)
        jplot.unfilledBar(raxis, b)
        jplot.unfilledBar(raxis, b.max()*np.cumsum(b)*1./b.sum(), color='green')
        plt.xlabel('Radius (mm)')
        plt.ylabel('Photons absorbed')

    if radius_range.max() < np.inf:
        # Do the vertical lines here
        plt.axvline(x=radius_range.min(), ymax=b.max(), linestyle='--', color='r')
        plt.axvline(x=radius_range.max(), ymax=b.max(), linestyle='--', color='r')

    # Axes for the projection
    paxis= np.linspace(-np.pi, np.pi, phi_bins)
    taxis = np.arccos(np.linspace(-1,1,theta_bins))[::-1] - np.pi/2.

    # Plotting the events in a projection - histogrammed
    fig2 = plt.figure(figsize=figsize)
    plt.title('Photons absorbed')
    ax = fig2.add_subplot(111, projection='mollweide')
    if plot_mode == 'histogram' or plot_mode == 'combo':
        b,x,y = np.histogram2d(theta, phi, [taxis, paxis])
        nonzero = b>0
        bmean = b[nonzero].mean()
        bstd  = b[nonzero].std()

        print 'Bin stats', b[nonzero].mean(), b[nonzero].std()

        plt.pcolor(paxis, taxis, b, 
                   vmin = bmean-plot_nstd*bstd, 
                   vmax = bmean+plot_nstd*bstd)
        plt.colorbar()

        return b, x, y

    if plot_mode == 'scatter' or plot_mode == 'combo':
        plt.plot(phi, theta, '.', color = 'black')


        
