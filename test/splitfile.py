#! /usr/bin/env python

import os
import math
import sys
import re
import ROOT

# Load ROOT libraries
ROOT.gSystem.Load('libttgamma.so')


newFile1 = ROOT.TFile.Open("root://cmseos.fnal.gov:1094//store/user/yumiceva/ttgamma/skim/v2/muon/skim_ttjets_1l_a.root",'RECREATE')
newFile2 = ROOT.TFile.Open("root://cmseos.fnal.gov:1094//store/user/yumiceva/ttgamma/skim/v2/muon/skim_ttjets_1l_b.root",'RECREATE')


openfile = ROOT.TFile.Open("root://cmseos.fnal.gov:1094//store/user/yumiceva/ttgamma/skim/v2/muon/skim_ttjets_1l.root",'READ')
hnew = openfile.Get("cutflow").Clone()
hnew2 = openfile.Get("h1test").Clone()
newFile1.cd()
hnew2.Write()
hnew.Write()

print " 1st Histogram done"

#openfile = ROOT.TFile.Open("root://cmseos.fnal.gov:1094//store/user/yumiceva/ttgamma/skim/v0/muon/skim_t_t.root",'READ')
hnew2 = openfile.Get("h1test").Clone()
hnew = openfile.Get("cutflow").Clone()
newFile2.cd()
hnew2.Write()
hnew.Write()

print "2nd Histogram done"


print "###########################"

print "Now the trees"

chain = ROOT.TChain("ggNtuplizer/EventTree")
maxEntries = -1
chain.Add("root://cmseos.fnal.gov:1094//store/user/yumiceva/ttgamma/skim/v2/muon/skim_ttjets_1l.root")
treeReader = ROOT.EventTree(chain)
 
# number of entries
numberOfEntries = chain.GetEntries()
print "Total number of entries to be processed: " + str(numberOfEntries)



####Cloning the trees#######
newFile1.cd()
newtree1 = chain.CloneTree(0)
newtree1.CopyAddresses(newtree1)
print "Cloned once"

newFile2.cd()
newtree2 = chain.CloneTree(0)
newtree2.CopyAddresses(newtree2)
print "Cloning done twice"

####looping over Entries and filling trees####


for entry in xrange(0,numberOfEntries):
	#print str(entry)
	treeReader.GetEntry(entry)
	if (entry <  ( (numberOfEntries-1)/2 ) ):
		newtree1.Fill()
	else : 
		newtree2.Fill()

#### Saving new trees in new root files


newFile1.cd()
#newtree1.Print()
print "TREE 1 Entries: "+ str(newtree1.GetEntries())
newFile1.mkdir("ggNtuplizer")
newFile1.cd("ggNtuplizer")
newtree1.Write()

newFile2.cd()
print "TREE 2 Entries: "+ str(newtree2.GetEntries())
newFile2.mkdir("ggNtuplizer")
newFile2.cd("ggNtuplizer")
newtree2.Write()

		
print "Splitting done"

