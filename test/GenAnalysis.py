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
    parser.add_option("-s", "--sample",
                      default="madgraph",
                      help="input samples. The options are: madgraph or whizard [default: %default]")
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
    outname = "results_gen.root"
    if options.sample == "madgraph":
        outname = "results_gen_madgraph.root"
    elif options.sample == "whizard":
        outname= "results_gen_whizard.root"
    else:
        print cTerm.RED+"Unknown sample name: "+options.sample+"\nOptions are \"madgraph\" or \"whizard\""+cTerm.RED
        sys.exit()
        
    outFile = ROOT.TFile(outname,"RECREATE")
    
    # Create chain of root trees
    chain = ROOT.TChain("ggNtuplizer/EventTree")
    maxEntries = -1

    # MG files
    print "Use dataset: "+options.sample
    
    if options.sample == "madgraph":
        chain.Add("/eos/uscms/store/user/iraklis/ggNtuples/job_summer12_ttinclusive.root")
        #chain.Add("/uscms_data/d2/maravin/TTG_MG5/Two2Seven/ROOT/ttgamma_27_part2.root")
    ##chain.Add("/uscms_data/d2/maravin/TTG_MG5/Two2Seven/ROOT/ttgamma_27_part3.root")
    ##chain.Add("/uscms_data/d2/maravin/TTG_MG5/Two2Seven/ROOT/ttgamma_27_part4.root")

    # Whizard files
    if options.sample == "whizard":
        chain.Add("/uscmst1b_scratch/lpc1/cmsroc/yumiceva/TTGamma/LHE/whizard/TTGamma_Whizard_2to7/ttgamma.root")
        maxEntries = 200000
        
    # setup ntuple object
    treeReader = ROOT.EventTree(chain)
    # number of entries
    numberOfEntries = chain.GetEntries()
    print "Total number of entries to be processed: " + str(numberOfEntries)
    
    # Get pointers to branches used in this analysis
    #Particles = treeReader.UseBranch("Particle")

    # Book histograms
    h_nocut = {}
    h_cut = {}
    h_nocut['pt'] = ROOT.TH1F("pt", "p_{T} [GeV]", 100, 0.0, 100.0)
    h_nocut['PID'] = ROOT.TH1F("PID","Particle ID",50,-25,25)
    h_nocut['ttbar_prod'] = ROOT.TH1F("ttbar_prod","ttbar production",4,0,4)
    h_nocut['ttbar_BR'] = ROOT.TH1F("ttbar_BR","BR",6,0,6)
    h_nocut['photon_mom'] = ROOT.TH1F("photon_mom","photon mother",50,-25,25)
    h_nocut['photon_pt'] = ROOT.TH1F("photon_pt","Photon p_{T} [GeV]",100,0,100)
    h_nocut['photon_eta'] = ROOT.TH1F("photon_eta","Photon #eta",100,-5,5)
    h_nocut['photon_phi'] = ROOT.TH1F("photon_phi","Photon #phi",80,-3.2,3.2)
    h_nocut['photon_deltaRLep'] = ROOT.TH1F('photon_deltaRLep',"#Delta R(#gamma,lepton)",100,0,6)
    h_nocut['N_photons'] = ROOT.TH1F("N_photons","Number of Photons",5,0,5)
    h_nocut['top_pt'] = ROOT.TH1F("top_pt","top p_{T} [GeV]",100,0,200)
    h_nocut['top_eta'] =ROOT.TH1F("top_eta","top #eta",100,-5,5)
    h_nocut['top_phi'] =ROOT.TH1F("top_phi","top #phi",80,-3.2,3.2)
    h_nocut['top_m'] = ROOT.TH1F("top_m","top mass [GeV]",50,150,200)
    h_nocut['W_m'] = ROOT.TH1F("W_m","W mass [GeV]",50,50,100)
    h_nocut['b_m'] = ROOT.TH1F("b_m","b mass [GeV]",50,0,10)
    h_nocut['b_pt'] = ROOT.TH1F("b_pt","b p_{T} [GeV]",100,0,200)
    h_nocut['b_eta'] = ROOT.TH1F("b_eta","b #eta",100,-5,5)
    h_nocut['b_phi'] =ROOT.TH1F("b_phi","b #phi",80,-3.2,3.2)
    h_nocut['q_pt'] = ROOT.TH1F("q_pt","b p_{T} [GeV]",100,0,200)
    h_nocut['q_eta'] = ROOT.TH1F("q_eta","b #eta",100,-5,5)
    h_nocut['q_phi'] =ROOT.TH1F("q_phi","b #phi",80,-3.2,3.2)
    h_nocut['nu_pt'] = ROOT.TH1F("nu_pt","Neutrino p_{T} [GeV]",100,0,200)
    h_nocut['nu_eta'] = ROOT.TH1F("nu_eta","Neutrino #eta",100,-5,5)
    h_nocut['nu_phi'] =ROOT.TH1F("nu_phi","Neutrino #phi",80,-3.2,3.2)
    h_nocut['lep_pt'] = ROOT.TH1F("lep_pt","Lepton p_{T} [GeV]",100,0,200)
    h_nocut['lep_eta'] = ROOT.TH1F("lep_eta","Lepton #eta",100,-5,5)
    h_nocut['lep_phi'] =ROOT.TH1F("lep_phi","Lepton #phi",80,-3.2,3.2)
            
    for key in h_nocut.keys():
        h_cut[key] = h_nocut[key].Clone(h_nocut[key].GetName())
        h_cut[key].Sumw2()
        h_nocut[key].Sumw2()
        h_cut[key].SetTitle( h_cut[key].GetTitle() )
        h_nocut[key].SetTitle( h_nocut[key].GetTitle() )
        
    # Loop over all events
    for entry in xrange(0, numberOfEntries):

        # Load selected branches with data from specified event
        treeReader.GetEntry(entry)

        if entry%2000 == 0:
            print "entry=",entry
            #### 
            if maxEntries!=-1 and maxEntries < entry:
                print cTerm.GREEN+"This sample has a maximum number of entries to process. Stop now."+cTerm.END
                break
        

        # loop over gen particles
        print "number of gen particles "+str(treeReader.nMC)
        for p in range(0,treeReader.nMC):
            h_nocut['PID'].Fill( treeReader.mcPID[p] )
    
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

            
