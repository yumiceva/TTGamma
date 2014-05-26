#! /usr/bin/env python

"""
   Script to create histograms and cut flow table for lepton+jets.
   
   Francisco Yumiceva (yumiceva@gmail.com)
   Florida Institute of Technology, 2013  
"""

class cTerm:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'                              

# External packages
import sys
import os
import math
import re
# Check if pyROOT is available
try:
    #from ROOT import *
    import ROOT
    print cTerm.GREEN+"ROOT module imported"+cTerm.END
except:
    print cTerm.RED+"\nError: Cannot load PYROOT, make sure you have setup ROOT in the path"
    print "and pyroot library is also defined in the variable PYTHONPATH, try:\n"+cTerm.END
    if (os.getenv("PYTHONPATH")):
        print " setenv PYTHONPATH ${PYTHONPATH}:$ROOTSYS/lib\n"
    else:
        print " setenv PYTHONPATH $ROOTSYS/lib\n"
    print "Exit now\n"
    sys.exit()
    
ROOT.PyConfig.IgnoreCommandLineOptions = True
from optparse import OptionParser

def main(argv = None):
    if argv == None:
        argv = sys.argv[1:]
    # OPTIONS
    usage = "usage: %prog [options]\n This script analyzes madgraph trees."
    parser = OptionParser(usage)
    parser.add_option("-b", "--batch",
                      action="store_true",
                      help="run ROOT in batch mode.")
    parser.add_option("-i", "--input",
                      default="input.root",
                      help="filename of input root file.")
    parser.add_option("-q","--quit",
                      action="store_true",
                      help="quit after reading tree otherwise prompt for keyboard to continue.")
    (options, args) = parser.parse_args(sys.argv[1:])
    #print options
    #print args
    return options

if __name__ == '__main__':

    options = main()
    
    if options.batch:
        ROOT.gROOT.SetBatch()
        print cTerm.GREEN+"Run in batch mode."+cTerm.END
    
    # Load ROOT libraries
    ROOT.gSystem.Load('libttgamma.so')

    # Create output root file
    in_name = options.sample
    #inFile = ROOT.TFile(in_name)
    #    print cTerm.RED+"Unknown sample name: "+options.sample+"\nOptions are \"madgraph\" or \"whizard\""+cTerm.RED
    #    sys.exit()

    N_outs = 2

    outFile = []
    outNames = []

    for i in range(0,N_outs):
        outname[i] = in_name[0:in_name.rfind(".root")]
        outname[i] += "_" + str(i+1) + ".root"
        outFile[i] = ROOT.TFile(outname[i],"RECREATE")
    
    # Create chain of root trees
    chain = ROOT.TChain("ggNtuplizer/EventTree")
    maxEntries = -1
    chain.Add(in_name)

    # create two readers, clone, etc.!!!!
    
    # setup ntuple object
    treeReader = ROOT.EventTree(chain)
    # number of entries
    numberOfEntries = chain.GetEntries()
    print "Total number of entries to be processed: " + str(numberOfEntries)

    # clone tree


    maxEntries = int(numberOfEntries/2)
    theNfile = 0
    # Loop over all events
    for entry in xrange(0, numberOfEntries):

        # Load selected branches with data from specified event
        treeReader.GetEntry(entry)

        if entry%1000 == 0:
            print "entry=",entry
            ####
        outFile[theNfile].Fill()
        
        if entry == maxEntries :
            print cTerm.GREEN+"split now"+cTerm.END


            
            
                
    outFile.cd()
    for key in h_nocut.keys():
        
        if h_nocut[key].GetEntries() > 0:
            h_nocut[key].Write()
        if h_cut[key].GetEntries() > 0:
            h_cut[key].Write()
        
    outFile.Close()
    print cTerm.GREEN+"Output file name: "+outname+cTerm.END

    if options.quit == False:
        # Wait
        rep = ''
        while not rep in [ 'q', 'Q', '.q', 'qq' 'p']:
            rep = raw_input( '\nenter: ["q",".q" to quit] ["p" or "print" to print all canvas]: ' )
            if 0<len(rep):
                if rep=='quit': rep = 'q'
            #if rep=='p' or rep=='print':

    #del(treeReader)
    #del(chain)
    #del treeReader
    print cTerm.GREEN+"done."+cTerm.END

            
