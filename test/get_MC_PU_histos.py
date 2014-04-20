
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
    
    prefix = "/eos/uscms/store/user/iraklis/ggNtuples/"
    prefix2= "/eos/uscms/store/user/makouski/"
    
    #files = get_list_files( prefix, "job_summer12_t" )
    #print files
    files = ["job_summer12_WWg.root","job_summer12_Wg.root","job_summer12_Wgg_FSR.root","job_summer12_Zg.root","job_summer12_diphoton_box_10to25.root","job_summer12_diphoton_box_250toInf.root","job_summer12_diphoton_box_25to250.root","job_summer12_t_s.root","job_summer12_t_t.root","job_summer12_t_tW.root","job_summer12_tbar_s.root","job_summer12_tbar_t.root","job_summer12_tbar_tW.root","job_summer12_ttW.root","job_summer12_ttZ.root","job_summer12_ttg.root","job_summer12_ttinclusive.root","job_summer12_DiPhotonBorn_Pt-10To25.root"]

    files_tmp = []
    for f in files:
        files_tmp.append( prefix + f)

    files_tmp.append (prefix2 + "job_summer12_Wjets.root")
    files = files_tmp
    
    outname = "MyMCPileupHistogram.root"
    tf_out = TFile( outname, "RECREATE")
    
    for f in files:

        print "Open: "+f
        tfile = TFile( f )
        MC_sample = string.replace( f, prefix, "" )
        MC_sample = string.replace( f, prefix2, "" )
        MC_sample = string.replace( MC_sample, "/", "")
        MC_sample = string.replace( MC_sample, "job_summer12_", "")
        MC_sample = string.replace( MC_sample, ".root", "")
        
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
