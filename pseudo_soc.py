import numpy as np
from scipy import optimize, ndimage
import jp_analysis as analysis
from detect_peaks import detect_peaks

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import jp_mpl as jplot
from copy import deepcopy




def getGausTime(n = [],
                    ybins = [],
                    prompt_peak_width = 4.,
                    plot = False,
                    debug = False):
    
    # ycenters
    ycenters = (ybins[1:] + ybins[:-1])/2.
    
    # bin width
    bin_width = ybins[1]-ybins[0]
    
    # Find the highest peak
    prompt_peak_index = n.argmax()

    if n[prompt_peak_index] < 20:
        return -1, -1
    
    # First peak range (in bins)
    fit_s = prompt_peak_index - int(prompt_peak_width/bin_width)
    fit_e = prompt_peak_index + int(prompt_peak_width/bin_width)
        
    
    # Fit a gaussian +/- prompt_peak_width ns around the peak
    popt =  [n[prompt_peak_index]/2., ycenters[prompt_peak_index],3.]
    try:
        popt, pcov = optimize.curve_fit(analysis.gaus,
                                        ycenters[fit_s:fit_e+1],
                                        n[fit_s:fit_e+1],
                                        popt,
                                       )
        perr = np.sqrt(np.diag(pcov))
        terr1 = perr[1]
        t1 = popt[1]
    except:
        terr1 = -1
        t1 = -1
    

    if plot:
        fig = plt.figure(figsize=(7,4))
        ax1 = fig.add_subplot(111)
        jplot.unfilledBar(ybins, n, color='0.7')

        newx = np.linspace(ycenters[fit_s]-3., ycenters[fit_e]+3., 101)
        plt.plot(newx, analysis.gaus(newx, *popt), 'r')
        
        plt.yscale('log')
        plt.ylabel('N hits')
        plt.ylim([1., n[prompt_peak_index]*1.3])

        return fig
    
        #raw_input()
    return t1, terr1



class FitLBpos(object):
    def __init__(self,
                 data = None,
                 error = None,
                 pmt_xyz = None,
                 pmtbool = None,
                 psup_radius = 8390.,
                 water_n = 1.38012,
                 print_call = True):
        self.c            = 0.299792458*1000 # mm/ns
        self.data    = deepcopy(data)
        self.pmt_xyz = deepcopy(pmt_xyz)
        self.pmtbool = deepcopy(pmtbool)
        self.water_n = water_n

        
        self.header_done = False


        self.print_call = print_call
        
        self.pmtbool[self.data <= 0] = False
        
        if error == None:
            self.error = np.ones_like(data)
        else:
            self.error = error + 1 # Adding one ns for everything
    
    def print_eval_header(self):
        print 'FCN \t\t t \t u \t v \t w'
        self.header_done = True
        
    def print_eval(self, value, t, u, v, w):
        if not self.header_done:
            self.print_eval_header()
        print value, '\t', t, '\t', u, '\t', v, '\t', w
        
    def __call__(self, t, u, v, w):
        this_pos = np.array([ u, v, w])

        water_c = self.c/self.water_n

        arr_time = t + np.linalg.norm(self.pmt_xyz - np.array([u,v,w]), axis=1)/water_c

        delta = (arr_time - self.data)[self.pmtbool]**2/self.error[self.pmtbool]**2
        delta = np.sum(delta)
        
        if self.print_call:
            self.print_eval(delta, t, u, v, w)
        
        return delta
