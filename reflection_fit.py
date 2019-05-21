import numpy as np
from scipy import optimize, ndimage
import jp_analysis as analysis
from detect_peaks import detect_peaks

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import jp_mpl as jplot
from copy import deepcopy

def getSimplePeakTimes(n = [],
                       ybins = [],
                       expected_tdelay = 76.,
                       refl_peak_width = 5.,
                       plot= False,
                       debug = False):
    # ycenters
    ycenters = (ybins[1:] + ybins[:-1])/2.
    
    # bin width
    bin_width = ybins[1]-ybins[0]
    
    # Find the highest peak
    prompt_peak_index = n.argmax()
    
    t1 = ybins[prompt_peak_index]
    
    # Now find the latest peak 
    refl_min_index = prompt_peak_index + int(expected_tdelay/bin_width) - int(refl_peak_width/bin_width)

    if refl_min_index >= n.size:
        t2 = -1
    else:
        refl_peak_index = refl_min_index + n[refl_min_index:].argmax()

        if n[refl_peak_index] < 1:
            t2 = -1
        else:
            t2 = ybins[refl_peak_index]

    terr1 = terr2 = 1.
    
    return t1, terr1, t2, terr2
    

def getPeakTimes(n = [],
                 ybins = [],
                 expected_tdelay = 76.,
                 refl_peak_width = 5.,
                 plot= False,
                 debug = False):
    # ycenters
    ycenters = (ybins[1:] + ybins[:-1])/2.

    # bin width
    bin_width = ybins[1]-ybins[0]

    # Find the highest peak
    prompt_peak_index = n.argmax()

    t1 = ybins[prompt_peak_index]

    smooth_coeff = 1./bin_width

    nsmooth = ndimage.filters.gaussian_filter1d(n, sigma=smooth_coeff)
    nsmooth[nsmooth<=0] = 1

    # Select a range where the second peak can be
    refl_peak_index = prompt_peak_index + int(expected_tdelay/bin_width)
    refl_si = refl_peak_index - int(refl_peak_width/bin_width)
    refl_ei = refl_peak_index + int((refl_peak_width*2)/bin_width)
    refl_ei = np.min([refl_ei, ybins.size-1])
    


    # Peak finding algorithm
    peaks = detect_peaks(nsmooth[refl_si:refl_ei], 
                         mpd = 10*smooth_coeff,
                         edge='rising',show=False, mph=1.1)
    if peaks.size == 0:
        t2 = -1
    else:
        t2 = ybins[refl_si:refl_ei][peaks].max()

    terr1 = terr2 = 1.

    if plot:
        fig = plt.figure(figsize=(7,4))
        ax1 = fig.add_subplot(111)
        jplot.unfilledBar(ybins, n, color='0.7')
        jplot.unfilledBar(ybins, nsmooth)
        plt.yscale('log')
        plt.plot((ybins[refl_si:refl_ei][peaks]+ybins[refl_si:refl_ei][peaks-1])/2., 
                 nsmooth[refl_si:refl_ei][peaks], 'xg')
        plt.ylabel('N hits')
        plt.axvline(x = t2, ymin=0, ymax = n.max()*1.2, 
                    linestyle='--', color='k')
        plt.axvline(x = t1, ymin=0, ymax = n.max()*1.2, 
                    linestyle='--', color='k')
        print t1, terr1, t2, terr2

    return t1, terr1, t2, terr2


def getDPeakTimes(n = [],
                     ybins = [],
                     expected_tdelay = 76.,
                  refl_peak_width=5.,
                     plot= False,
                     debug = False):

    if n.sum() == 0:
        return -1, -1, -1, -1
    
    # ycenters
    ycenters = (ybins[1:] + ybins[:-1])/2.

    # bin width
    bin_width = ybins[1]-ybins[0]


    smooth_coeff = 1./bin_width
    nsmooth = ndimage.filters.gaussian_filter1d(n, sigma=smooth_coeff)
    nsmooth[nsmooth<=0] = 1

    # Finding the peaks of the derivative
    dx = (ybins[1]-ybins[0])/2.
    ydn = ybins[1:-1]

    # Derivating, and smoothing over it
    dnsmooth = np.diff(np.log10(nsmooth))/dx  
    dnsmooth = ndimage.filters.gaussian_filter1d(dnsmooth, 
                                                 sigma=smooth_coeff)

    # Highest peak is easy
    prompt_minus = int(10./bin_width)
    prompt_region = [n.argmax()-prompt_minus, n.argmax()]

    try:
        prompt_peak_index = prompt_region[0]+dnsmooth[prompt_region[0]:prompt_region[1]].argmax()
    except:
        # Errors appear here when the histogram is too empty
        # print 'Not using this PMT'
        return -2, -2, -2, -2
    
    t1 = ydn[prompt_peak_index]
    refl_peak_index = prompt_peak_index + int(expected_tdelay/bin_width)
    refl_si = refl_peak_index - int(refl_peak_width/bin_width)
    refl_ei = refl_peak_index + int((refl_peak_width*2)/bin_width)
    refl_ei = np.min([refl_ei, ybins.size-1])

    #print refl_si, refl_ei

    # Finding the peaks
    dpeaks = detect_peaks(dnsmooth[refl_si:refl_ei], 
                          mpd = 10*smooth_coeff,
                          edge='rising',show=False, mph=0.001)

    if dpeaks.size == 0:
        t2 = -1
    else:
        t2 = ydn[refl_si:refl_ei][dpeaks].max()

    terr1 = terr2 = 1.


    if plot:
        fig = plt.figure(figsize=(7,4))
        ax1 = fig.add_subplot(111)

        print 'Peaks', dpeaks
        jplot.unfilledBar(ybins, np.log10(n)*dnsmooth.max()/np.log10(n.max()), color='0.7')
        jplot.plot(ydn, dnsmooth)

        jplot.plot([ydn[refl_si],ydn[refl_ei]], [1.,1.], 'or')

        if len(dpeaks) > 0:
            plt.plot((ydn[refl_si:refl_ei][dpeaks]+ydn[refl_si:refl_ei][dpeaks-1])/2.,
                     dnsmooth[refl_si:refl_ei][dpeaks], 'xg')
            plt.ylabel('N hits')
            plt.axvline(x = t2, ymin=0, ymax = dnsmooth.max()*1.2,
                        linestyle='--', color='k')
            plt.axvline(x = t1, ymin=0, ymax = dnsmooth.max()*1.2,
                        linestyle='--', color='k')

        print t1, terr1, t2, terr2


    return t1, terr1, t2, terr2



def getGausTimes(n = [],
                 ybins = [],
                 expected_tdelay = 76.,
                 prompt_peak_width = 4.,
                 refl_peak_width = 4.5,
                 second_gaus = True,
                 plot = False,
                 debug = False):
    
    # ycenters
    ycenters = (ybins[1:] + ybins[:-1])/2.
    
    # bin width
    bin_width = ybins[1]-ybins[0]

    
    # Find the highest peak
    prompt_peak_index = n.argmax()
    
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
    
    # Fit a gaussian around the second peak, expected_tdelay-refl_peak_width+2*refl_peak_width
    refl_peak_index = prompt_peak_index + int(expected_tdelay/bin_width)
    fit_s2 = refl_peak_index - int(refl_peak_width/bin_width)
    fit_e2 = refl_peak_index + int((refl_peak_width)/bin_width)
    fit_e2 = np.min([fit_e2, ybins.size-1])
    
    if refl_peak_index > (ybins.size-1):
        terr2 = -1
        t2 = -1
    else:
        if second_gaus:        
            try:
                popt2 = [n[refl_peak_index]/2.,ycenters[refl_peak_index],3.]
                popt2, pcov2 = optimize.curve_fit(analysis.gaus,
                                                  ycenters[fit_s2:fit_e2+1],
                                                  n[fit_s2:fit_e2+1],
                                                  popt2)
                perr2 = np.sqrt(np.diag(pcov2))
                terr2 = perr2[1]
                t2 = popt2[1]
            except:
                # There are errors here when the ToA histogram is empty
                terr2 = -1
                t2 = -1
        else:
            t2 = ycenters[fit_s2 + n[fit_s2:].argmax()]
            terr2 = -1
    if plot:

        print 'Times', t1, terr1, t2, terr2
        fig = plt.figure(figsize=(7,4))
        ax1 = fig.add_subplot(111)
        jplot.unfilledBar(ybins, n, color='0.7')

        newx = np.linspace(ycenters[fit_s],
                           ycenters[fit_e], 101)
        plt.plot(newx, analysis.gaus(newx, *popt), 'r')
        print t1+expected_tdelay
        plt.axvline(x = t1 + expected_tdelay-refl_peak_width,
                    ymin=0, ymax = n.max()*1.2,
                    linestyle='--', color = 'k')
        plt.axvline(x = t1 + expected_tdelay+refl_peak_width,
                    ymin=0, ymax = n.max()*1.2,
                    linestyle='--', color = 'k')
        if second_gaus:
            try:
                newx = np.linspace(ycenters[fit_s2],
                                   ycenters[fit_e2], 101)
                plt.plot(newx, analysis.gaus(newx, *popt2), 'r')    
                plt.xlim([ycenters[prompt_peak_index]-20., ycenters[refl_peak_index] + 20])
            except:
                print 'Could not plot second gaus (out of bounds)'
        else:
            plt.axvline(x = t2, ymin =0, ymax = n.max()*1.2)
        
        plt.yscale('log')
        plt.ylabel('N hits')
        plt.ylim([1., n[prompt_peak_index]*1.3])

        return fig
    
        #raw_input()
    return t1, terr1, t2, terr2



class FitLBpos(object):
    def __init__(self,
                 data = None,
                 error = None,
                 pmt_xyz = None,
                 pmtbool = None,
                 psup_radius = 8390.,
                 water_n = 1.34389,
                 print_call = True):
        self.c            = 0.299792458*1000 # mm/ns
        self.data    = deepcopy(data)
        self.pmt_xyz = deepcopy(pmt_xyz)
        self.pmtbool = deepcopy(pmtbool)
        self.pmt_r   = np.linalg.norm(pmt_xyz,axis=1)
        self.header_done = False

        self.R = psup_radius
        #self.water_c = c/water_n
        self.print_call = print_call
        
        self.pmtbool[self.data <= 0] = False
        
        if error == None:
            self.error = np.ones_like(data)
        else:
            self.error = error + 1 # Adding one ns for everything
    
    def print_eval_header(self):
        print 'FCN \t\t u \t v \t w \t n'
        self.header_done = True
        
    def print_eval(self, value, u, v, w, n):
        if not self.header_done:
            self.print_eval_header()
        print value, '\t', u, '\t', v, '\t', w, '\t', n
        
    def __call__(self, u, v, w, n):
        this_pos = np.array([ u, v, w])
        water_c = self.c/n
        tdiff = 2*self.R/water_c* \
                (1+np.dot(self.pmt_xyz,this_pos)/(self.pmt_r*self.R))
        
        delta = (tdiff - self.data)[self.pmtbool]**2/self.error[self.pmtbool]**2
        delta = np.sum(delta)
        
        if self.print_call:
            self.print_eval(delta, u, v, w, n)
        
        return delta
