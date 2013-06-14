void runttgamma(TString sample="all", bool getLogs=false)  
{

  TString desdir = "/uscms/home/yumiceva/lpc1/4tops2013_templates/";

  TProof *p = TProof::Open("lite://");

  //p->AddDynamicPath("");
  //p->Exec("gSystem->Load(\"/uscms/home/yumiceva/work/CMSSW_5_3_3/lib/slc5_amd64_gcc462/libYumicevaTop7TeV.so\")");
  //p->AddIncludePath("/uscms/home/yumiceva/work/CMSSW_5_3_3/src/");
  //p->AddIncludePath("-I/uscmst1/prod/sw/cms/slc5_amd64_gcc462/external/boost/1.47.0/include/");

  p->Archive(" ",desdir);

  //p->AddInput(new TNamed("PROOF_OUTPUTFILE_LOCATION", "LOCAL"));

  if (sample=="MC"||sample=="2w2l"||sample=="all")
    {
      TDSet *mc_2w2l = new TDSet("EventTree","*","/ggNtuplizer");
      mc_2w2l->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_WW_2l2nu.root"); 
      mc_2w2l->Process("ttgamma.C+","sample=2w2l");
    }
    
  if (sample=="MC"||sample=="wwg"||sample=="all")
    {
      TDSet *mc_wwg = new TDSet("EventTree","*","/ggNtuplizer");
      mc_wwg->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_WWg.root"); 
      mc_wwg->Process("ttgamma.C+","sample=wwg");
    }    
    
  if (sample=="MC"||sample=="wz2l2q"||sample=="all")
    {
      TDSet *mc_wz2l2q = new TDSet("EventTree","*","/ggNtuplizer");
      mc_wz2l2q->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_WZ_2l2q.root"); 
      mc_wz2l2q->Process("ttgamma.C+","sample=wz2l2q");
    }
    
  if (sample=="MC"||sample=="wz3l"||sample=="all")
    {
      TDSet *mc_wz3l = new TDSet("EventTree","*","/ggNtuplizer");
      mc_wz3l->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_WZ_3lnu.root"); 
      mc_wz3l->Process("ttgamma.C+","sample=wz3l");
    }      
    
  if (sample=="MC"||sample=="wg"||sample=="all")
    {
      TDSet *mc_wg = new TDSet("EventTree","*","/ggNtuplizer");
      mc_wg->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_Wg.root"); 
      mc_wg->Process("ttgamma.C+","sample=wg");
    }       
    
  if (sample=="MC"||sample=="zz2e2mu"||sample=="all")
    {
      TDSet *mc_zz2e2mu = new TDSet("EventTree","*","/ggNtuplizer");
      mc_zz2e2mu->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_ZZ_2e2mu.root"); 
      mc_zz2e2mu->Process("ttgamma.C+","sample=zz2e2mu");
    }     
    
  if (sample=="MC"||sample=="zz2e2tau"||sample=="all")
    {
      TDSet *mc_zz2e2tau = new TDSet("EventTree","*","/ggNtuplizer");
      mc_zz2e2tau->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_ZZ_2e2tau.root"); 
      mc_zz2e2tau->Process("ttgamma.C+","sample=zz2e2tau");
    }       
    
  if (sample=="MC"||sample=="zz2mu2tau"||sample=="all")
    {
      TDSet *mc_zz2mu2tau = new TDSet("EventTree","*","/ggNtuplizer");
      mc_zz2mu2tau->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_ZZ_2mu2tau.root"); 
      mc_zz2mu2tau->Process("ttgamma.C+","sample=zz2mu2tau");
    }       
    
  if (sample=="MC"||sample=="zz4e"||sample=="all")
    {
      TDSet *mc_zz4e = new TDSet("EventTree","*","/ggNtuplizer");
      mc_zz4e->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_ZZ_4e.root"); 
      mc_zz4e->Process("ttgamma.C+","sample=zz4e");
    }   
    
  //  if (sample=="MC"||sample=="zz4mu"||sample=="all")
  //    {
  //      TDSet *mc_zz4mu = new TDSet("EventTree","*","/ggNtuplizer");
  //      mc_zz4mu->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_ZZ_4mu.root"); 
  //      mc_zz4mu->Process("ttgamma.C+","sample=zz4mu");
  //    }   
    
  if (sample=="MC"||sample=="zz4tau"||sample=="all")
    {
      TDSet *mc_zz4tau = new TDSet("EventTree","*","/ggNtuplizer");
      mc_zz4tau->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_ZZ_4tau.root"); 
      mc_zz4tau->Process("ttgamma.C+","sample=zz4tau");
    }   
    
  if (sample=="MC"||sample=="zg"||sample=="all")
    {
      TDSet *mc_zg = new TDSet("EventTree","*","/ggNtuplizer");
      mc_zg->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_Zg.root"); 
      mc_zg->Process("ttgamma.C+","sample=zg");
    }     
    
  if (sample=="MC"||sample=="diphoton1"||sample=="all")
    {
      TDSet *mc_diphoton1 = new TDSet("EventTree","*","/ggNtuplizer");
      mc_diphoton1->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_diphoton_box_10to25.root"); 
      mc_diphoton1->Process("ttgamma.C+","sample=diphoton1");
    }        
    
 if (sample=="MC"||sample=="diphoton2"||sample=="all")
    {
      TDSet *mc_diphoton2 = new TDSet("EventTree","*","/ggNtuplizer");
      mc_diphoton2->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_diphoton_box_25to250.root"); 
      mc_diphoton2->Process("ttgamma.C+","sample=diphoton2");
    }     
    
 if (sample=="MC"||sample=="ts"||sample=="all")
    {
      TDSet *mc_ts = new TDSet("EventTree","*","/ggNtuplizer");
      mc_ts->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_t_s.root"); 
      mc_ts->Process("ttgamma.C+","sample=ts");
    }    
    
 if (sample=="MC"||sample=="tt"||sample=="all")
    {
      TDSet *mc_tt = new TDSet("EventTree","*","/ggNtuplizer");
      mc_tt->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_t_t.root"); 
      mc_tt->Process("ttgamma.C+","sample=tt");
    } 
    
 if (sample=="MC"||sample=="2tw"||sample=="all")
    {
      TDSet *mc_2tw = new TDSet("EventTree","*","/ggNtuplizer");
      mc_2tw->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_t_tW.root"); 
      mc_2tw->Process("ttgamma.C+","sample=2tw");
    }    
    
 if (sample=="MC"||sample=="tbars"||sample=="all")
    {
      TDSet *mc_tbars = new TDSet("EventTree","*","/ggNtuplizer");
      mc_tbars->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_tbar_s.root"); 
      mc_tbars->Process("ttgamma.C+","sample=tbars");
    }     
    
 if (sample=="MC"||sample=="tbart"||sample=="all")
    {
      TDSet *mc_tbart = new TDSet("EventTree","*","/ggNtuplizer");
      mc_tbart->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_tbar_t.root"); 
      mc_tbart->Process("ttgamma.C+","sample=tbart");
    }        
    
 if (sample=="MC"||sample=="tbartw"||sample=="all")
    {
      TDSet *mc_tbartw = new TDSet("EventTree","*","/ggNtuplizer");
      mc_tbartw->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_tbar_tW.root"); 
      mc_tbartw->Process("ttgamma.C+","sample=tbartw");
    }   
    
 if (sample=="MC"||sample=="ttw"||sample=="all")
    {
      TDSet *mc_ttw = new TDSet("EventTree","*","/ggNtuplizer");
      mc_ttw->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_ttW.root"); 
      mc_ttw->Process("ttgamma.C+","sample=ttw");
    }    
    
 if (sample=="MC"||sample=="ttz"||sample=="all")
    {
      TDSet *mc_ttz = new TDSet("EventTree","*","/ggNtuplizer");
      mc_ttz->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_ttZ.root"); 
      mc_ttz->Process("ttgamma.C+","sample=ttz");
    }     
    
 if (sample=="MC"||sample=="ttg"||sample=="all")
    {
      TDSet *mc_ttg = new TDSet("EventTree","*","/ggNtuplizer");
      mc_ttg->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_ttg.root"); 
      mc_ttg->Process("ttgamma.C+","sample=ttg");
    }  
    
 if (sample=="MC"||sample=="ttjets1"||sample=="all")
    {
      TDSet *mc_ttjets1 = new TDSet("EventTree","*","/ggNtuplizer");
      mc_ttjets1->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_ttjets_1l.root"); 
      mc_ttjets1->Process("ttgamma.C+","sample=ttjets1");
    }   
    
 if (sample=="MC"||sample=="ttjets2"||sample=="all")
    {
      TDSet *mc_ttjets2 = new TDSet("EventTree","*","/ggNtuplizer");
      mc_ttjets2->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/GGNtuMC/job_summer12_ttjets_2l.root"); 
      mc_ttjets2->Process("ttgamma.C+","sample=ttjets2");
    }    
    
  
  if (sample=="data"||sample=="all")
    {
      TDSet *data = new TDSet("EventTree","*","/ggNtuplizer");
      data->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012a_Jul13rereco_skim.root");
      data->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012b_Jul13rereco_skim.root");
      data->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012a_Aug6rereco_skim.root");
      data->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012c_Aug24rereco_skim.root");
      data->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012c_Dec11rereco_skim.root");
      data->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012c_PRv2_skim.root");
      data->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012d_PRv1_part1_skim.root");
      data->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012d_PRv1_part2_skim.root");
      data->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012d_PRv1_part3_skim.root");
      data->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012d_PRv1_part4_skim.root"); 
      data->Add("/eos/uscms/store/user/makouski/ggNtupleElePhoJetSkim/job_1electron_2012d_PRv1_part5_skim.root"); 
      
      data->Process("ttgamma.C+","verbose sample=data");
      // get log files
      if (getLogs)
	{
	  logList = p->GetManager()->GetSessionLogs()->GetListOfLogs();
	  for (int i=1; i< logList->GetSize(); ++i)
	    {
	      logElem = ( TProofLogElem* ) logList->At( i );
	      macro = logElem->GetMacro();
	      macro->SaveSource("data_ttgamma_"+TString(Form("%i",i))+".stdout");
	    }
	}
    }
}
