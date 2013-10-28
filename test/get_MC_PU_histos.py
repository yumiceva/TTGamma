
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
    
    prefix = "/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim2/GGNtuMC"

    files = get_list_files( prefix )
    #print files
    #files = ["job_summer12_WW_2l2nu.root"] #,"job_summer12_WWg.root","job_summer12_WZ_2l2q.root"]
    
    outname = "MyMCPileupHistogram.root"
    tf_out = TFile( outname, "RECREATE")
    
    for f in files:
        
        print "Open: "+f
        tfile = TFile( f )
        MC_sample = string.replace( f, prefix, "" )
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
