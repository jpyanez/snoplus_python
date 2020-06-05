import os, sys, pickle
import matplotlib.pyplot as plt
import numpy as np
import jp_mpl as jplot
import ROOT, rat

def arrayDiagnostics(values, name, ax_lim = [None, None]):
    print 'Diagnostics', name
    print values.min(),values.max(), values.mean()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if ax_lim[0] == None:
        ax_lim[0] = 0.
    if ax_lim[1] == None:
        ax_lim[1] = np.median((values[values>0]))*10
    xaxis = np.linspace(ax_lim[0], ax_lim[1], 100)
    b, x = np.histogram(values, xaxis)
    
    jplot.unfilledBar(x,b)
    plt.yscale('log')
    plt.xlabel(name)
    return fig, ax

def dumpArray( indir, branch_name, outdir = None,
               show_info = False, ax_lim = [None, None], max_files = -1):
    isdata = False
    print indir
    if 'data' in indir:
        isdata = True
        print '\n\n***Using data!'
        
    file_list = [os.path.join(indir, x) for x in os.listdir(indir) if '.root' in x]
    file_list.sort()
    if max_files >0:
        file_list = file_list[:max_files]
    print 'Going over ', len(file_list), ' files'

    if branch_name == 'energy':
        print 'Using RSP energy correction AND CALIBRATION'
        rc = rat.utility().GetReconCorrector().Get()
        calib = rat.utility().GetReconCalibrator().Get()
        print 'Test RSP correction', rc.CorrectEnergyRSP(5.0)
        print 'Test Energy calibration', calib.CalibrateEnergyRSP(False, 4., 300., 40)
    # Test one file for content
    tfile = ROOT.TFile(file_list[0])
    tree = tfile.Get('output')
    branch = tree.GetBranch(branch_name)
    tot_entries = branch.GetEntries()

    if tot_entries < 100:
        tot_entries = 100
    
    print 'Entries in file ', tot_entries, ' - adding 20% contintengcy'
    print file_list[0]
    
    values = np.zeros(int(tot_entries*len(file_list)*1.2))
    counter = 0
        
    for ifile, one_file in enumerate(file_list):
        if ifile % 400 ==  0:
            print ifile

        try:
            tfile = ROOT.TFile(one_file)
            tree = tfile.Get('output')
            branch = tree.GetBranch(branch_name)
            tot_entries = branch.GetEntries()
        except:
            print 'File cannot be opened ', one_file
            continue

        # Special loop for energy to avoid doing many ifs
        if branch_name == 'energy':
            posx_branch = tree.GetBranch('posx')
            posy_branch = tree.GetBranch('posy')
            posz_branch = tree.GetBranch('posz')
            None
            for i in range(tot_entries):
                branch.GetEntry(i)
                posx_branch.GetEntry(i)
                posy_branch.GetEntry(i)
                posz_branch.GetEntry(i)
                
                corrected =  rc.CorrectEnergyRSP(branch.GetLeaf(branch_name).GetValue(0))
                # Calibrated!
                rho = np.sqrt(posx_branch.GetLeaf('posx').GetValue(0)**2 + posy_branch.GetLeaf('posy').GetValue(0)**2)
                calibrated = calib.CalibrateEnergyRSP(isdata,
                                                      corrected,
                                                      rho,
                                                      posz_branch.GetLeaf('posz').GetValue(0))
                values[counter] = calibrated
                counter += 1
        else:
            for i in range(tot_entries):
                branch.GetEntry(i)
                values[counter] = branch.GetLeaf(branch_name).GetValue(0)
                counter += 1
        tfile.Close()
                    
    values = values[:counter]
                    
    if outdir != None:
        np.save(os.path.join(outdir, branch_name+'.npy'), values)
        print 'Saved!'
                
    if show_info:
        try:
            fig, ax=arrayDiagnostics(values,branch_name, ax_lim)
            fig.savefig(os.path.join(outdir, branch_name+'.png'))
        except:
            print 'Diagnostics failed'
    return values
