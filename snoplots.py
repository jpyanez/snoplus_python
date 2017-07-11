import matplotlib.pyplot as plt
import numpy as np
import os, sys
sys.path.append('/home/jp/projects/python_tools')
import jp_mpl as jplot

mycolors = ['b','r','g','m','c', 'y', '0.6']

def plotComparison(data=None, datasets = [], datakey = '', scale_factor = [], 
                   labels = [], figname = '', xaxis = np.arange(0, 100, 2),
                   xlabel = 'Nhits', ylabel = 'Entries', 
                   outdir = None, verbose=False):

    if len(scale_factor) == 0:
        scale_factor = [1.]*len(datasets)
    if len(labels) == 0:
        labels = datasets


    # Loop over the different levels
    for k, level in enumerate(['top','middle','bottom', 'all']):
        # Reading each of the files one by one
        if verbose:        print '\n****', level , '*****'
        nbins = []
        myfig = plt.figure(figsize=(8,5))    
        
        if level=='all':
            xaxis = np.arange(1.3*xaxis.max(), 2.3*xaxis.max(), xaxis[1]-xaxis[0])
        
        for i, one_set in enumerate(datasets):
            n, x = np.histogram(data[one_set][datakey][:,k], xaxis)
            nbins.append(n)

            if verbose:
                print '\n',labels[i]

                print 'SUM ', data[one_set][datakey][:,k].sum()*scale_factor[i]
                print 'Mean ', data[one_set][datakey][:,k].mean()
                print 'Std  ', data[one_set][datakey][:,k].std()

            if i == 0:
                baseline = data[one_set][datakey][:,k].mean()

            if i > 0:
                print one_set, level, 1-data[one_set][datakey][:,k].mean()/baseline

            jplot.unfilledBar(xaxis, nbins[-1]*scale_factor[i], 
                                  color = mycolors[i])
            jplot.errorMark(xaxis, nbins[-1]*scale_factor[i], 
                            error=np.sqrt(nbins[-1])*scale_factor[i], color=mycolors[i],
                            label =labels[i] + '\n' + \
                                r'$\mu$' + '=' + "%.2f" % data[one_set][datakey][:,k].mean() + '\n'+\
                                r'$\sigma$' + '=' + "%.2f" % data[one_set][datakey][:,k].std()) 

        plt.title(level)



        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        #plt.ylim([0,])
        plt.legend(loc=0,ncol=1)
        if len(figname) > 0:
            myfig.savefig(os.path.join(outdir, figname + '_' + level + '.png'), dpi=300)
            myfig.savefig(os.path.join(outdir, figname + '_' + level + '.pdf'))
