
#include <stdlib.h>
#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCint.h>
#include <TMath.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <fstream>
#include "stdio.h"
#include "string"
#include "TStyle.h"
#include "Riostream.h"
#include "TAxis.h"
#include <TRandom3.h>
#include <cctype>
#include <cmath>
#include <vector>
#include <TArrow.h>
#include <TH1F.h>
#include <TH1.h>
#include <TEntryList.h>
#include <TLegend.h>
#include<THStack.h>
#include<TColor.h>
#include<Rtypes.h>
#include<TPaveStats.h>
#include<TVirtualHistPainter.h>
using std::cout;
using std::endl;

void stack() 
{
  //Data histograms
  //TFile *f1 = new TFile("results_data.root","READ");
  //TH1F *h1 = (TH1F*)f1->Get("h14test");
  
  //MC Histogram
  //TFile *f2 = new TFile("results_2tw.root","READ");
  //TH1F *h2 = (TH1F*)f2->Get("h14test");  
  
  //THStack stack("stack","final");
  //stack.Add(&h1);
  //stack.Add(&h2);
  
  //h1->SetfillColor(kBlue);
  //h2->SetfillColor(kRed);
  
  //stack->Draw();

  TCanvas *c0 = new TCanvas("c0","",0,0,500,600);
  c0->Divide(1,2);
  TCanvas *c1 = new TCanvas("c1","",0,0,500,600);
  THStack *hs = new THStack("hs","Number of Jets");

  TFile *f1 = new TFile("results_data.root");
  TFile *f2 = new TFile("results_ttjets2.root");
    
  //TFile *f1 = new TFile("results_data.root");
  f1->cd();
  TH1F *h1 = (TH1F*)f1->Get("h14test");
  //h1->SetDirectory(0);
  //h1->Print();
  //h1->SetFillColor(kBlue);
  c0->cd(1);
  h1->SetMarkerStyle(2);
  h1->SetMarkerSize(1.2);
  h1->Draw("P");
  hs->Add(h1);
  //h1->Draw();          
  
  //TFile *f2 = new TFile("results_tt.root");
  f2->cd();
  TH1F *h2 = (TH1F*)f2->Get("h14test");
  //h2->SetDirectory(0);
  //h2->Print();
  h2->SetFillColor(kRed);
  h2->SetLineColor(kRed);
  //h2->Scale(500);
  c0->cd(2);
  h2->Draw();
  hs->Add(h2);
  //h2->Draw("SAME");  
    
  //h2.SetFillColor(kRed);
  c1->cd();
  hs->Print();
  hs->Draw();
  
  
  
  
}
  
  
