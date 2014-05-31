#! /usr/bin/env python

from ROOT import *
from decimal import Decimal
import sys,os, math

path = "./"
IsMC = False
Lumi = 19630.0 # in pb-1
XSecFile = "cuyFiles/CrossSections8TeV.txt"

if len(sys.argv) < 2:
    print "usage: mkCutFlowTables.py <path to files> <>"
    sys.exit()

if len(sys.argv) > 1:

    path = sys.argv[1]
    print " path: "+ path
else:
    print "Need to provide path to files. Exit."
    sys.exit()
    
print " Luminosity for MC samples: "+str(Lumi)

Nevents = {}
xsec = {} # xsection in pb
file = {}
label = {}

# Data file
Datafile['data_all'] = path+'/results_data.root'
# MC files
file['ttjets_0l'] = path+'/results_ttjets_0l.root'
file['ttjets_1l'] = path+'/results_ttjets_1l.root'
file['ttjets_2l'] = path+'/results_ttjets_2l.root'
file['Wjets'] = path+'/results_WJets.root'
file['DYJetsToLL'] = path+'/results_DYJetsToLL.root'
file['t_s']   = path+'/results_t_s.root'
file['tbar_s']  = path+'/results_tbar_s.root'
file['t_t']   = path+'/results_t_t.root'
file['tbar_t']  = path+'/results_tbar_t.root'
file['t_tW']   = path+'/results_t_tW.root'
file['tbar_tW']  = path+'/results_tbar_tW.root'
file['ttW']   = path+'/results_ttW.root'
file['ttZ']  = path+'/results_ttZ.root'
file['ttg']  = path+'/results_ttg.root'
file['ttjets_0l_g'] = path+'/results_ttjets_0l_g.root'
file['ttjets_1l_g'] = path+'/results_ttjets_1l_g.root'
file['ttjets_2l_g'] = path+'/results_ttjets_2l_g.root'
#file['WW']    = path+'/results_WW.root'
#file['WZ']    = path+'/results_WZ.root'
    
    
label['ttbar'] = '$t\\overline{t}$'
label['QCD'] = 'QCD'
label['Wjets'] = '$W\\rightarrow l\\nu$'
label['Zjets'] = '$Z\\rightarrow l^{+}l^{-}$'
label['tch'] = 'top t-ch'
label['tWch'] = 'top tW-ch'
label['sch'] = 'top s-ch'
label['tch_bar'] = 'anti-top t-ch'
label['tWch_bar'] = 'anti-top tW-ch'
label['sch_bar'] = 'anti-top s-ch'
label['WW'] = 'WW'
label['WZ'] = 'WZ'
                                    
label['Total'] = 'Total MC'
label['data'] = "Data";

# data files
    #file['Jun14'] = ''
    #file['MB']    = '/uscmst1b_scratch/lpc1/cmsroc/yumiceva/top_prod_Sep22/Jun14/ttmuj_Jun14_MB_PSep23.root'#'/uscms_data/d3/ttmuj/Documents/NtupleMaker/Data/2.88pb-1/ttmuj_MB_Sep3.root'
    #file['Jul16'] = '/uscmst1b_scratch/lpc1/cmsroc/yumiceva/top_prod_Sep22/Jul16/ttmuj_Jul16_PSep23.root'#'/uscms_data/d3/ttmuj/Documents/NtupleMaker/Data/2.88pb-1/ttmuj_Jul16_Sep3.root'
    #file['Prompt']= '/uscms_data/d3/ttmuj/Documents/NtupleMaker/Data/2.88pb-1/ttmuj_Prompt_Sep3.root'

    #label['Jun14'] = 'Jun 14th Mu'
    #label['MB'] = 'Jun 14th MinimumBias'
    #label['Jul16'] = 'Jul 16th'
    #label['Prompt'] = 'PromptReco'

    #if os.path.isfile(path+'/JPT/cutflow_JPT_data.txt'):
    #    file['dataJPT'] = path+'/cutflow_JPT_data.txt'
    #    label['dataJPT'] = 'JPT Ref. Sel.'
    #if os.path.isfile(path+'/calo/cutflow_calo_data.txt'):
    #    file['datacalo'] = path+'/cutflow_calo_data.txt'
    #    label['datacalo'] = 'Calo Ref. Sel.'
    #if os.path.isfile(path+'/PF/cutflow_PF_data.txt'):
    #    file['dataPF'] = path+'/cutflow_PF_data.txt'
    #    label['dataPF'] = 'PF Ref. Sel.'

                     
    
keys = file.keys()
cutflow = {}
cutflowerr = {}
cutlabel = {}
#cutlabel['Processed'] = 'Processed'
#cutlabel['Cleaning'] = 'CleanFilters'
#cutlabel['HLT'] = 'HLT'
#cutlabel['GoodPV'] = 'Skim'
#cutlabel['OneIsoMuon'] = 'One Iso $\mu$'
#cutlabel['LooseMuVeto'] = 'Loose $\mu$ veto'
#cutlabel['LooseEleVeto'] = 'Electron veto'
#cutlabel['1Jets'] = 'jets $> 0$'
#cutlabel['2Jets'] = 'jets $> 1$'
#cutlabel['3Jets'] = 'jets $> 2$'
cutlabel['4Jets'] = 'N jets >= 4'
cutlabel['Onebtag'] = 'btags $> 0$'
cutlabel['PhotonFiducial'] = '
cutlabel['topmass'] = 'top mass'

cutlabelvector = [ 'GoodPV', 'OneIsoMu', 'LooseMuVeto', 'ElectronVeto', 'MET', '1Jet', '2Jet', '3Jet', '4Jet','2Jet1b','2Jet2b','MaxJets','phi','topmass']
SKIPCUTS = ['3Jet','4Jet','MaxJets','phi','topmass']

allmap = {}
allmaperr = {}

weightmap = {}

tablelist = ['ttbar','Wjets','Zjets','QCD','tch','tch_bar','tWch','tWch_bar','sch','sch_bar']


if not IsMC:
    tablelist = ['data']

                            
for sample in tablelist:

    if sample=="Total": continue
    print " processing " + sample

    cutmap = {}
    cutmaperr = {}
    #txtfile = open(file[sample])
    roofile = TFile(file[sample])
    hcutflow = gDirectory.Get("cutflow")
    print "  entries in cutflow histogram: " + str(hcutflow.GetEntries())

    for ibin in xrange(1, hcutflow.GetNbinsX() +1 ):

        skipthiscut = False
        for skipcut in SKIPCUTS:
            if skipcut == cutlabelvector[ibin-1]: skipthiscut = True
        if skipthiscut:
            print "skip counting cut name: "+cutlabelvector[ibin-1]
            continue
        cutname = cutlabelvector[ibin-1]
        acut = hcutflow.GetBinContent(ibin)
        #print cutname
        #print acut
        #if sample=="data":
        if sample.find("data")!=-1:
            cutmap[ cutname ] = int(float(acut))
            cutmaperr[ cutname ] = math.sqrt( cutmap[cutname] )
        else:
            cutmap[ cutname ] = Decimal(str(acut))
            cutmaperr[ cutname ] = cutmap[cutname].sqrt() #math.sqrt( cutmap[cutname] )
            
    roofile.Close()
    
    scale = 1.
    if IsMC and Lumi>0:
        scale = Decimal( str( Lumi * xsec[ sample ] / Nevents[sample] ))
        print "sample weight = "+ str( xsec[ sample ] / Nevents[sample] )
        weightmap[sample] = xsec[ sample ] / Nevents[sample]
    ilabel = 0
    for key in cutmap.keys():

        cutmap[ key ] = scale * cutmap[ key]
        cutmaperr[key]= scale * cutmaperr[key] * scale * cutmaperr[key]
        
        print " cut "+str(key) + " ("+cutlabel[key]+") "+" = "+str( round(cutmap[key],1) )

        if allmap.has_key(key):
            allmap[ key ] += cutmap[ key ]
            allmaperr[ key ] += cutmaperr[key]
        else:
            allmap[ key ] = cutmap[ key ]
            allmaperr[ key ] = cutmap[key]
        ilabel += 1
    cutflow[ sample ] = cutmap
    cutflowerr[ sample ] = cutmaperr
    
print " TOTAL"
ilabel=0
for key in allmap.keys():
    print " cut "+str(key) + " ("+cutlabel[key]+") "+" = "+str( round(allmap[key],1) )
    ilabel += 1
    
cutflow["Total"] = allmap
cutflowerr["Total"] = allmaperr

# write latex
#sortedcutlist = ['CleanFilters','HLT','GoodPV','OneIsoMu','LooseMuVeto','ElectronVeto','MET','1Jet','2Jet','3Jet','4Jet']
sortedcutlist = ['GoodPV','OneIsoMu','LooseMuVeto','ElectronVeto']#,'MET','1Jet','2Jet','2Jet1b']
sortedcutlist2= ['MET','1Jet','2Jet','2Jet1b','2Jet2b']

if IsMC:
    cutlabel['CleanFilters'] = 'Processed'

    
texname = "cutflow_"+JetType+"_data.tex"

if IsMC:
    texname = "cutflow_"+JetType+"_MC.tex"
    if showWprime:
        texname = "cutflow_"+JetType+"_MC_Wp"+WprimeType+".tex"
        
outtex = open(texname,"w")

sss = " & "

# Write Header
outtex.write('''
\\begin{table}[h]
\\begin{centering}
''')
if IsMC:
    aline = '\\begin{tabular}{|l|'+'c|'*(len(sortedcutlist)+1) +'} \hline \n'
    if Lumi<=0:
        aline = '\\begin{tabular}{|l|'+'c|'*len(sortedcutlist) +'} \hline \n'
    outtex.write(aline)
    #outtex.write('Cut & ttbar & Wjets & Zjets & QCD & STtch \hline')
    line = "Sample"
    for icut in sortedcutlist:
        line += sss+cutlabel[icut]
    outtex.write(line+" \\\ \hline \n")
else:
    aline = '\\begin{tabular}{|l|r|r|r|} \hline \n'
    
    outtex.write(aline)
    line = "Sample"
    for icut in sortedcutlist:
        line += sss+cutlabel[icut]
    outtex.write(line+" \\\ \hline \n")

ilabel = 0                   
#for key in allmap.keys():
# Write content
for sample in tablelist:
    
    line = label[sample]
    
    for key in sortedcutlist:
        cutmap = cutflow[sample]
        cutmaperr = cutflowerr[sample]
        acut = cutmap[key]
        acuterr = cutmaperr[key]
        stracut = str(int(acut))
        stracuterr = str(round(math.sqrt(acuterr),1))
        #stracut = stracut.strip('.0')
        if IsMC:
            acut = round(cutmap[key],1)
            stracut = str(acut)
            line += sss + stracut + " $\\pm$ " + stracuterr
        else:
            line += sss + stracut
    if sample=="Total":line = '\\hline \n'+line
    outtex.write(line+" \\\ \n")
    ilabel += 1
    
outtex.write(''' \hline
\end{tabular}
%\caption{MC cutflow}\label{tab:tab}
\end{centering}
\end{table}
''')


if len(sortedcutlist2) > 0:

    
    # Write Header
    outtex.write('''
\\begin{table}[h]
\\begin{centering}
    ''')
    if IsMC:
        aline = '\\begin{tabular}{|l|'+'c|'*(len(sortedcutlist2)+1) +'} \hline \n'
        if Lumi<=0:
            aline = '\\begin{tabular}{|l|'+'c|'*len(sortedcutlist2) +'} \hline \n'
        outtex.write(aline)
        #outtex.write('Cut & ttbar & Wjets & Zjets & QCD & STtch \hline')
        line = "Sample"
        for icut in sortedcutlist2:
            line += sss+cutlabel[icut]
        outtex.write(line+" \\\ \hline \n")
    else:
        aline = '\\begin{tabular}{|l|r|r|r|} \hline \n'
            
        outtex.write(aline)
        line = "Sample"
        for icut in sortedcutlist2:
            line += sss+cutlabel[icut]
        outtex.write(line+" \\\ \hline \n")

    ilabel = 0
    #for key in allmap.keys():
    # Write content
    for sample in tablelist:
        
        line = label[sample]

        for key in sortedcutlist2:
            cutmap = cutflow[sample]
            cutmaperr = cutflowerr[sample]
            acut = cutmap[key]
            acuterr = cutmaperr[key]
            stracut = str(int(acut))
            stracuterr = str(round(math.sqrt(acuterr),1))
            #stracut = stracut.strip('.0')
            if IsMC:
                acut = round(cutmap[key],1)
                stracut = str(acut)
                line += sss + stracut + " $\\pm$ " + stracuterr
            else:
                line += sss + stracut
        if sample=="Total":line = '\\hline \n'+line
        outtex.write(line+" \\\ \n")
        ilabel += 1
            
    outtex.write(''' \hline
\end{tabular}
%\caption{MC cutflow}\label{tab:tab}
\end{centering}
\end{table}
    ''')
    

print "file "+texname+ " has been written."

if IsMC:
    
    print "\n MC Weights"
    tmplistsamples = weightmap.keys()
    tmplistsamples.sort()
    for sample in tmplistsamples:

        print sample+" "+str(weightmap[sample])


    

