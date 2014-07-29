#!/bin/bash

job=$1

echo "Job start time: `date`"

cmsswDir="/uscms_data/d2/yumiceva/CMSSW_5_3_14_patch1/src"
sourceDir="/uscms_data/d2/yumiceva/CMSSW_5_3_14_patch1/src"
#outputDir="root://cmseos.fnal.gov//store/user/yumiceva/ttgamma/skim/v1/muon/"
#outputDir="/eos/uscms/store/user/yumiceva/ttgamma/skim/v1/muon/"
#outputDir="/uscms/home/yumiceva/work/sframe/CMSSW_5_3_3/src/TTGamma/test/histos"
#outputDir= "${_CONDOR_SCRATCH_DIR}/"
#storeDir="/eos/uscms/store/user/yumiceva/ttgamma/skim/v2/muon/"
storeDir="root://cmseos.fnal.gov//store/user/yumiceva/ttgamma/skim/v2/muon/"


echo "1. Current directory is: `pwd`"
echo "condor_scratch_dir = ${_CONDOR_SCRATCH_DIR}"

cd ${cmsswDir}
echo "2. Current directory is: `pwd`"
source /uscmst1/prod/sw/cms/shrc uaf
eval `scram runtime -sh`
cd -

echo "3. Current directory is: `pwd`"

#if [[ ${job} == 0 ]]; then inputDir="/eos/uscms/store/user/makouski/ggNtupleFullSimWHIZARD/";
#else if [[ ${job} == 12 || ${job} == 13 ]]; then inputDir="/eos/uscms/store/user/makouski/ggNtupleTTJets_v11/";
#else inputDir="/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/ggNtuples/";
#fi
#fi

#samples=("GGNtuMC/ttg_WHIZARD_v11.root ${inputDir}ggtree_mc_*" \
#"job_1electron_2012a_Aug6rereco_skim.root ${inputDir}job_1electron_2012a_Aug6rereco_skim.root" \
samples=(\
#"DYJetsToLL" \
#"Wjets" \
#"W3jets" \
#"W4jets"
##"WWg" \
#"Wg" \
##"Wgg_FSR" \
#"Zg" \
##"diphoton_box_10to25" \
##"diphoton_box_25to250" \
##"diphoton_box_250toInf" \
#"t_s" \
#"t_t" \
#"t_tW" \
#"tbar_s" \
#"tbar_t" \
#"tbar_tW" \
#"ttW" \
#"ttZ" \
##"ttg" \
#"ttgWhizard"
#"ttjets_1l" \
#"ttjets_2l" \
#"ttjets_0l"
"data_mu_a" \
#"data_mu_b" \
#"data_mu_c" \
#"data_mu_d"
)

# change to scratch area
#cd ${_CONDOR_SCRATCH_DIR}

# change to source directory to access helper files
#cd ${sourceDir}
#echo "Current directory is: `pwd`"

options="skim muon outdir=${_CONDOR_SCRATCH_DIR}"

echo "${sourceDir}/ttgamma ${samples[job]} ${options}"
${sourceDir}/ttgamma ${samples[job]} ${options}
echo "ttgamma executable finshed at: `date`"

echo "Now transfer file out of node to ${storeDir}"
echo ${storeDir}
#echo "cp skim_${samples[job]}.root ${storeDir}/."
#cp skim_${samples[job]}.root ${storeDir}/.
echo "xrdcp skim_${samples[job]}.root ${storeDir}/skim_${samples[job]}.root -f"
xrdcp skim_${samples[job]}.root ${storeDir}/skim_${samples[job]}.root -f -d 1
echo "transfer finish at: `date`"
echo "Copy done."

echo "Job done."
