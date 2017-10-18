import numpy as np
from scipy import optimize, ndimage
import jp_analysis as analysis
from detect_peaks import detect_peaks

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import jp_mpl as jplot
from copy import deepcopy

def getGausTime(residual_histo = None,
                residual_axis = None,
                time_guess = None,
                sigma = 4.,
                min_hits = 300):


    # Ask for at least 300 hits
    if residual_histo.sum() < min_hits:
        return -1, -1

    xcenters = (residual_axis[1:] + residual_axis[:-1])/2.
    

    # Mask the bins to use
    binbool = (xcenters>(time_guess-sigma))*(xcenters<(time_guess+sigma))
    if np.sum(residual_histo[binbool]) < 3:
        return -1, -1

    popt    = [residual_histo[binbool].max(), time_guess, sigma/2.]

    try:
        popt, pcov = optimize.curve_fit(analysis.gaus,
                                        xcenters[binbool],
                                        residual_histo[binbool],
                                        popt,
                                        )
        perr = np.sqrt(np.diag(pcov))
        terr = perr[1]
        time = popt[1]
    except:
        terr = -1
        time = -1

    if terr > 20:
        return -1, -1
    
    return time, terr



class FitLBpos(object):
    def __init__(self,
                 data = None, # Should be a dictionary with all the information
                 pmt_xyz = None,
                 pmtbool = None,
                 psup_radius = 8390.,
                 fit_mode    = 'advanced', # can be 'advanced'
                 fix_xy      = False,
                 print_call = True):

        self.data = data

        self.c            = 0.299792458*1000 # mm/ns

        self.pmt_xyz = deepcopy(pmt_xyz)
        self.pmtbool = deepcopy(pmtbool)
        self.pmt_r   = np.linalg.norm(pmt_xyz,axis=1)
        self.fix_xy  = fix_xy
        # In the fit mode "advanced" we need to know which is the PMT that produces the reflection
        if fit_mode == 'advanced':
            self.pmt_mirror = np.zeros_like(self.pmt_r, dtype=int)
            for pmtid in range(self.pmt_r.size):
                residuals = np.linalg.norm(pmt_xyz[pmtid,:] + pmt_xyz, axis=1)
                self.pmt_mirror[pmtid] = residuals.argmin()
                
            self.pmt_r_mirror = np.linalg.norm(pmt_xyz[self.pmt_mirror,:],axis=1)

            self.tdelay_fcn = self.delta_toa_advanced


        elif fit_mode == 'simple':
            self.tdelay_fcn = self.delta_toa_simple


        if fix_xy:
            self.x = None
            self.y = None


        self.header_done = False

        self.R = psup_radius
        self.print_call = print_call
        

        self.first_peak = np.zeros(self.data['time_residuals'].shape[1])
        self.first_peak_error = np.zeros_like(self.first_peak)

        # Calculate the position of the first peak here (all PMTs)
        print 'Calculating first peak position'
        for iPMT in range(self.first_peak_error.size):
            self.first_peak[iPMT],  self.first_peak_error[iPMT] = \
                getGausTime(residual_histo = self.data['time_residuals'][:,iPMT],
                            residual_axis = self.data['residual_axis'],
                            time_guess = 0.,
                            sigma = 4.)

        # PMTs without an error aren't any good - discard them
        self.pmtbool[self.first_peak_error < 0] = False

        # Adding the time delay subtracted when creating histograms
        # self.first_peak += self.data['pmt_delays']
        self.second_peak = []
        self.second_peak_error = []


    # Here lb_pos is a number, not an array
    def delta_toa_xyfixed(self, lb_pos_z, water_c):
        full_lb_pos = np.array([self.x, self.y, lb_pos_z])
        return self.tdelay_aux(full_lb_pos, water_c)

    
    def print_eval_header(self):
        print 'FCN \t\t u \t v \t w \t n \t nbool'
        self.header_done = True
        
    def print_eval(self, value, u, v, w, n, nbool):
        if not self.header_done:
            self.print_eval_header()
        print value, '\t', u, '\t', v, '\t', w, '\t', n, '\t', nbool


    def delta_toa_simple(self,lb_pos, water_c):
        return 2*self.R/water_c* \
            (1+np.dot(self.pmt_xyz, lb_pos)/(self.pmt_r*self.R))

    def delta_toa_advanced(self, lb_pos, water_c):
        return ( np.linalg.norm(self.pmt_xyz - self.pmt_xyz[self.pmt_mirror,:],axis=1) +
                 np.sqrt(self.pmt_r_mirror**2 -
                         2*np.dot(self.pmt_xyz[self.pmt_mirror,:], lb_pos) + np.linalg.norm(lb_pos)**2) -
                 np.sqrt(self.pmt_r**2 -
                         2*np.dot(self.pmt_xyz, lb_pos) + np.linalg.norm(lb_pos)**2) )/water_c

    #def delta_toa_full()
        
    def __call__(self, *args):

        if self.fix_xy:
            u = self.x
            v = self.y
            w = args[0]
            n = args[1]
        else:
            u = args[0]
            v = args[1]
            w = args[2]
            n = args[3]

        this_pos = np.array([u,v,w])

        # Calculated time difference
        water_c = self.c/n
        self.tdiff = self.tdelay_fcn(this_pos, water_c)

        
        # The second peak is first calculated with the initial gues fo the position - it better be good!
        # Run the fit iteratively to remove this issue
        # Only do this once
        if len(self.second_peak_error) == 0:
            print 'Calculating the second peak position'
            self.second_peak = np.zeros_like(self.first_peak)
            self.second_peak_error = np.zeros_like(self.first_peak)
            for iPMT in range(self.second_peak_error.size):
                self.second_peak[iPMT], self.second_peak_error[iPMT] = \
                    getGausTime(residual_histo = self.data['time_residuals'][:,iPMT],
                                residual_axis  = self.data['residual_axis'],
                                time_guess = self.tdiff[iPMT],
                                sigma = 4.)

            # Selecting only PMTs with an error estimate - warning, will change at every step
            self.thisbool = self.pmtbool*(self.second_peak_error > 0)

            # Get the data time diff
            self.data_timediff = self.second_peak - self.first_peak
            self.data_error    = np.sqrt(self.second_peak_error**2 + self.first_peak_error**2)
        
        
        delta = (self.tdiff - self.data_timediff)[self.thisbool]**2/self.data_error[self.thisbool]**2
        delta = np.sum(delta)
        
        if self.print_call:
            self.print_eval(delta, u, v, w, n, np.sum(self.thisbool))
        
        return delta/np.sum(self.thisbool)
