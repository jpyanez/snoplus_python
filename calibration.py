import ROOT, rat
import os, sys, pickle


indir = '/home/jp/projects/snoplus/laserball_calibration/all_runs'
outdir = '/home/jp/projects/snoplus/laserball_calibration/all_runs_pckl'

du = rat.utility()
du.LoadDBAndBeginRun()
def wavelengthToEnergy(wavelength): # wl in nm
    hc = 1.2398 #eV * um
    return 1E-3*hc/wavelength # in MeV
light_path = du.GetLightPathCalculator()
shadow_path = du.GetShadowingCalculator()
laserE = wavelengthToEnergy(420.) # This won't matter for direct paths

def isPMTshadowed(source_position,
                 pmt_position,
                 tolerance = 0.):
    shadow_path.SetAllGeometryTolerances(tolerance)
    light_path.CalcByPosition(ROOT.TVector3(source_position),
                              ROOT.TVector3(pmt_position),
                              laserE, 0.)
    return shadow_path.CheckForShadowing(light_path)


def simplifySOCfiles():
    file_list = os.listdir(indir)

    for file_index, file_name in enumerate(file_list):

        print 'Doing file', file_name ,
        outfile_path = os.path.join(outdir, file_name.rstrip('.root')+'.pckl')
        if os.path.isfile(outfile_path):
            print ' ... exists. Skipping it.'
            continue

        qhs       = np.zeros(pmt_radius.size)
        qhl       = np.zeros_like(qhs)
        occupancy = np.zeros_like(qhl)
        not_valid = 0
        reader = rat.socreader(os.path.join(indir, file_name))

        for soc, run in reader:
            # Get the run information
            result  = soc.GetFitResult(soc.GetFitNames()[0])
            vertex  = np.array(result.GetVertex(0).GetPosition())
            calib   = np.array(soc.GetCalib().GetPos())
            npulses = soc.GetNPulsesTriggered()

            # Go over the PMTs of interest
            # Should add some quality criteria here for PMTs
            for one_id in range(pmt_radius.size):
                try:
                    pmt = soc.GetSOCPMT(one_id)
                except:
                    not_valid +=1
                    continue
                qhs[one_id] = np.sum(np.array(pmt.GetQHSs()))
                qhl[one_id] = np.sum(np.array(pmt.GetQHLs()))
                occupancy[one_id] = pmt.GetPromptOccupancy()


        reader.close()
        print '... done! - not valid PMTs ', not_valid

        # Save the file
        ofile = open(outfile_path, 'w')
        pickle.dump({'vertex':vertex, 'calib':calib, 'npulses':npulses,
                     'qhs':qhs, 'qhl':qhl, 'occupancy':occupancy}, ofile)
        ofile.close()
