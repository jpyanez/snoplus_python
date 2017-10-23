import numpy as np
import pickle

class PMTsInRing(object):
    def __init__(self,
                 effective_n = 1.4,
                 pmt_positions = '/home/jpyanez/snoplus/snoplus_python/pmt_positions.pckl'):
        if type(pmt_positions) == str:
            self.pmt_positions = pickle.load(open(pmt_positions))['xyz']
        elif type(pmt_positions) == np.ndarray:
            self.pmt_positions = pmt_positions
        else:
            print 'Make sure you pass either the address of a file with PMT positions or the positions themselves'

        self.ch_angle = np.arccos(1./effective_n)

        print 'Using cherenkov angle', np.rad2deg(self.ch_angle)

    def __call__(self,
                 position = None,
                 direction = None,
                 tolerance = 0.015): # In radians):
        normdir = direction/np.linalg.norm(direction)
        pmt_vectors=self.pmt_positions-position
        pmt_vectors=(pmt_vectors/
                     np.reshape(np.linalg.norm(pmt_vectors,axis=1),[pmt_vectors.shape[0],1]))
        angles_to_track = np.arccos(np.dot(pmt_vectors,normdir))
        pmts_cone = (np.abs(angles_to_track-self.ch_angle) < tolerance)
        print 'Selected PMTs:', np.sum(pmts_cone)

        return pmts_cone
        
