
# Import tools
from ROOT import *
from math import sqrt
import sys
import re
import array
import os
import string

def get_list_files(directory,pattern = ""):

    dir = []
    
    dir = os.listdir(directory)
    
    lfiles = []
    for f in dir:
        
        if f.find(pattern) != -1 and f.endswith(".root"):
            
            lfiles.append(directory+"/"+f)
            
    return lfiles
                                    

def main():
    
    prefix = "root://cmsxrootd-site.fnal.gov//store/user/iraklis/ggNtuples/"
    prefix2= "root://cmsxrootd-site.fnal.gov//store/user/makouski/"
    
    #files = get_list_files( prefix, "job_summer12_t" )
    #print files
    #files = [#"job_summer12_WWg.root",
             #"job_summer12_Wg.root",
             #"job_summer12_Wgg_FSR.root",
             #"job_summer12_Zg.root",
             #"job_summer12_diphoton_box_10to25.root",
             #"job_summer12_diphoton_box_250toInf.root",
             #"job_summer12_diphoton_box_25to250.root",
             #"job_summer12_t_s.root",
             #"job_summer12_t_t.root",
             #"job_summer12_t_tW.root",
             #"job_summer12_tbar_s.root",
             #"job_summer12_tbar_t.root",
             ##"job_summer12_tbar_tW.root",
             #"job_summer12_ttW.root",
             #"job_summer12_ttZ.root",
             #"job_summer12_ttg.root",
             #"job_summer12_ttjets_1l.root",
             #"job_summer12_DiPhotonBorn_Pt-10To25.root"]

    files_tmp = []
    #for f in files:
    #    files_tmp.append( prefix + f)

    #files_tmp.append (prefix2 + "job_summer12_Wjets.root")
    #files_tmp.append (prefix2 + "job_summer12_DYJetsToLL.root")
    #files_tmp.append (prefix2 + "job_summer12_ttjets_2l.root")
    
    #files = files_tmp
    #files = [ prefix2 + "job_summer12_ttjets_2l.root" ]
    files = ["/eos/uscms/store/user/iraklis/ggNtuples/job_summer12_ttjets_1l.root"]
    #files = ["/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/iraklis/Makouski/CMSSW_5_3_12/src/ggAnalysis/ggNtuplizer/test/W3JetsToLNu_TuneZ2Star_8TeV-madgraph.root"]
    #files = ["/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/iraklis/Makouski/CMSSW_5_3_12/src/ggAnalysis/ggNtuplizer/test/W4JetsToLNu_TuneZ2Star_8TeV-madgraph.root"]
    #files = ["/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/misha/job_summer12_whizard_2to5_ttA.root"]
        #"/eos/uscms/store/user/troy2012/GG_MC_12/TTJets_Hadronic/TTJets_hadronic.root" ]
    
    outname = "MyMCPileupHistogram.root"
    tf_out = TFile( outname, "UPDATE") # update will append new histograms
    
    for f in files:

        print "Open: "+f
        tfile = TFile.Open( f )
        MC_sample = string.replace( f, prefix, "" )
        MC_sample = string.replace( MC_sample, prefix2, "" )
        MC_sample = string.replace( MC_sample, "/", "")
        MC_sample = string.replace( MC_sample, "job_summer12_", "")
        MC_sample = string.replace( MC_sample, ".root", "")
        #MC_sample = "ttjets_0l"
        #MC_sample = "ttgWhizard"
        #MC_sample = "W3jets"
        #MC_sample = "W4jets"
        MC_sample = "ttjets_1l_b"
        
        h = tfile.Get("//ggNtuplizer/hPUTrue")
        hnew_name = "hPUTrue_"+MC_sample
        print " Write: " + hnew_name
        hnew = h.Clone( hnew_name )
    
        tf_out.cd()
        hnew.Write()

        tfile.Close()

    tf_out.Close()

if __name__ =='__main__':
    sys.exit(main())
