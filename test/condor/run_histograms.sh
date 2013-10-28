#!/bin/bash

job=$1
#job=37

cmsswDir="/uscms/home/yumiceva/work/sframe/CMSSW_5_3_3/src"
sourceDir="/uscms/home/yumiceva/work/sframe/CMSSW_5_3_3/src/TTGamma/test"
outputDir="/eos/uscms/store/user/yumiceva/ttgamma/histos/v0/"

cd ${cmsswDir}
echo "Current directory is: `pwd`"
source /uscmst1/prod/sw/cms/shrc uaf
eval `scram runtime -sh`


#if [[ ${job} == 0 ]]; then inputDir="/eos/uscms/store/user/makouski/ggNtupleFullSimWHIZARD/";
#else if [[ ${job} == 12 || ${job} == 13 ]]; then inputDir="/eos/uscms/store/user/makouski/ggNtupleTTJets_v11/";
#else inputDir="/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/ggNtuples/";
#fi
#fi

#samples=("GGNtuMC/ttg_WHIZARD_v11.root ${inputDir}ggtree_mc_*" \
#"job_1electron_2012a_Aug6rereco_skim.root ${inputDir}job_1electron_2012a_Aug6rereco_skim.root" \
samples=("WW_2l2nu" \
"WWg" \
"WZ_2l2q" \
"WZ_3lnu" \
"Wg" \
"ZZ_2e2mu" \
"ZZ_2e2tau" \
"ZZ_2mu2tau" \
"ZZ_4e" \
"ZZ_4mu" \
"ZZ_4tau" \
"Zg" \
"diphoton_box_10to25" \
"diphoton_box_25to250" \
"t_s" \
"t_t" \
"t_tW" \
"tbar_s" \
"tbar_t" \
"tbar_tW" \
"ttW" \
"ttZ" \
"ttg" \
"ttgWz" \
"data")

# change to scratch area
cd ${_CONDOR_SCRATCH_DIR}

options="outdir=${outputDir}"

echo "${sourceDir}/ttgamma ${samples[job]} ${options}"
${sourceDir}/ttgamma ${samples[job]} ${options}
