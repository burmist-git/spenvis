//root
#include <TROOT.h>
#include <TStyle.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TString.h>
#include <TFile.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLine.h>
#include <TMath.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMultiGraph.h>

//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>

using namespace std;

void load_data_file(TString fileInName, TGraph *gr_tot, TGraph *gr_diff, Double_t &hMin, Double_t &hMax);
void norm_graph(TGraph *gr, TGraph *gr_norm);
void load_data_file_solar_prot(TString fileInName, TGraph *gr_tot, TGraph *gr_diff);

Int_t plots_prot(){
  //
  Double_t hMin;
  Double_t hMax;
  //
  //
  //535
  TGraph *gr_tot_535km_p = new TGraph();
  TGraph *gr_dif_535km_p = new TGraph();
  TGraph *gr_dif_535km_p_norm = new TGraph();
  load_data_file("./data/535_km_p.dat", gr_tot_535km_p, gr_dif_535km_p, hMin, hMax);
  //norm_graph(gr_dif_535km_e_norm, gr_dif_550km_e);  
  //
  //
  TGraph *gr_tot_solar_prot = new TGraph();
  gr_tot_solar_prot->SetNameTitle("gr_tot_solar_prot","gr_tot_solar_prot");
  TGraph *gr_diff_solar_prot = new TGraph();
  gr_diff_solar_prot->SetNameTitle("gr_diff_solar_prot","gr_diff_solar_prot");
  load_data_file_solar_prot("./data/535km_solar_prot_peakflux.dat", gr_tot_solar_prot, gr_diff_solar_prot);
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  c1->SetRightMargin(0.01);
  c1->SetLeftMargin(0.12);
  c1->SetTopMargin(0.01);
  c1->SetBottomMargin(0.08);
  //
  c1->SetLogy();
  c1->SetLogx();
  c1->SetGridx();
  c1->SetGridy();
  //
  gr_dif_535km_p->SetLineColor(kBlack);
  gr_dif_535km_p->SetLineWidth(3.0);
  //
  gr_diff_solar_prot->SetLineColor(kBlue+2);
  gr_diff_solar_prot->SetLineWidth(3.0);
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr_dif_535km_p);
  //mg->Add(gr_flux_p_PDG);
  mg->Add(gr_diff_solar_prot);
  //mg->SetMinimum(1.0e-2);
  //mg->SetMaximum(10000.0);
  //
  //
  mg->Draw("APL");
  mg->GetXaxis()->SetTitle("Energy, GeV");
  mg->GetYaxis()->SetTitle("Differential Flux, 1/m/m/s/GeV");
  //mg->GetYaxis()->SetTitle("Relative differential flux");
  mg->GetXaxis()->SetRangeUser(1.0e-3,400);
  //gr_dif_535km_p->GetXaxis()->SetTitle("Energy, GeV");
  //gr_dif_535km_p->GetYaxis()->SetTitle("Differential Flux, 1/m/m/s/sr/GeV");
  //gr_dif_535km_p->Draw("APL");
  //
  //
  TLegend *leg = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg->AddEntry(gr_dif_535km_p, "Trapped protons (SPENVIS 535 km)","l");
  leg->AddEntry(gr_diff_solar_prot, "Solar protons (SPENVIS 535 km)","l");
  leg->Draw();
  return 0;
}

void load_data_file_solar_prot(TString fileInName, TGraph *gr_tot, TGraph *gr_diff){
  ifstream fileIn(fileInName.Data());
  string mot;
  Double_t ekin;
  Double_t ftot;
  Double_t fdiff;
  if(fileIn.is_open()){
    fileIn>>mot;
    fileIn>>mot;
    fileIn>>mot;
    while(fileIn>>ekin>>ftot>>fdiff){
      gr_tot->SetPoint(gr_tot->GetN(),ekin,ftot/10000.0);
      if(fdiff>0.0)
	gr_diff->SetPoint(gr_diff->GetN(),ekin,fdiff/10000.0);
    }
    fileIn.close();
  }
  else{
    cout<<"Unable to open file \n";
  }
}

void load_data_file(TString fileInName, TGraph *gr_tot, TGraph *gr_diff, Double_t &hMin, Double_t &hMax){
  ifstream fileIn(fileInName.Data());
  string mot;
  Double_t ekin;
  Double_t ftot;
  Double_t fdiff;
  if(fileIn.is_open()){
    fileIn>>mot;
    fileIn>>hMax;
    fileIn>>hMin;
    fileIn>>mot>>mot>>mot;
    while(fileIn>>ekin>>ftot>>fdiff){
      gr_tot->SetPoint(gr_tot->GetN(),ekin,ftot/TMath::Pi()/2.0);
      if(fdiff>0.0)
	gr_diff->SetPoint(gr_diff->GetN(),ekin,fdiff/TMath::Pi()/2.0);
    }
    fileIn.close();
  }
  else{
    cout<<"Unable to open file \n";
  }
}

void norm_graph(TGraph *gr, TGraph *gr_norm){
  Double_t x,y;
  Double_t norm;
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i,x,y);
    gr_norm->GetPoint(i,x,norm);
    if(norm>0.0)
      gr->SetPoint(i,x,y/norm);
  }
}
