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

void load_data_file(TString fileInName, TGraphErrors *gr, Int_t &np, Double_t &eMin, Double_t &eMax, Double_t &integral);
void load_data_file(TString fileInName, TGraph *gr, Int_t &np, Double_t &eMin, Double_t &eMax, Double_t &integral);
void get_eff(TGraph *gr_sim, TGraph *gr_trg, TGraph *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV);
void get_eff(TGraphErrors *gr_sim, TGraphErrors *gr_trg, TGraphErrors *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV);
void get_eff_area(TGraph *gr_sim, TGraph *gr_trg, TGraph *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV, Double_t area);
void get_eff_area(TGraphErrors *gr_sim, TGraphErrors *gr_trg, TGraphErrors *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV, Double_t area);

Int_t plots_LST_Crab_paper_eff(){
  //
  Double_t eMin_TeV = 0.02;
  Double_t eMax_TeV = 20.0;  
  //
  Int_t np_sim;
  Double_t eMin_sim;
  Double_t eMax_sim;
  Double_t integral_sim;
  //
  Int_t np_trg;
  Double_t eMin_trg;
  Double_t eMax_trg;
  Double_t integral_trg;
  //
  Int_t np_sel;
  Double_t eMin_sel;
  Double_t eMax_sel;
  Double_t integral_sel;
  //
  TGraph *gr_sim = new TGraph();
  TGraph *gr_trg = new TGraph();
  TGraph *gr_sel = new TGraph();
  TGraph *gr_trg_eff = new TGraph();
  TGraph *gr_sel_eff = new TGraph();
  TGraph *gr_trg_effArea = new TGraph();
  //
  Int_t np_sim_LST_AdvCam, np_trg_LST_AdvCam;
  Double_t eMin_sim_LST_AdvCam, eMax_sim_LST_AdvCam, integral_sim_LST_AdvCam;
  Double_t eMin_trg_LST_AdvCam, eMax_trg_LST_AdvCam, integral_trg_LST_AdvCam;
  TGraphErrors *gr_sim_LST_AdvCam = new TGraphErrors();
  TGraphErrors *gr_trg_LST_AdvCam = new TGraphErrors();
  TGraphErrors *gr_trg_eff_AdvCam = new TGraphErrors();
  TGraphErrors *gr_trg_effArea_AdvCam = new TGraphErrors();
  TGraphErrors *gr_trg_effArea_ratio = new TGraphErrors();
  //
  load_data_file("./data/LST_Crab_paper_Zenith_23.63deg_sim.csv", gr_sim, np_sim, eMin_sim, eMax_sim, integral_sim);
  load_data_file("./data/LST_Crab_paper_Zenith_23.63deg_trg.csv", gr_trg, np_trg, eMin_trg, eMax_trg, integral_trg);
  load_data_file("./data/LST_Crab_paper_Zenith_23.63deg_sel.csv", gr_sel, np_sel, eMin_sel, eMax_sel, integral_sel);
  //
  get_eff(gr_sim, gr_trg, gr_trg_eff, 16, 0.02, 20.0);
  get_eff(gr_sim, gr_sel, gr_sel_eff, 16, 0.02, 20.0);
  get_eff_area(gr_sim, gr_trg, gr_trg_effArea, 16, 0.02, 20.0, TMath::Pi()*764.06453*764.06453);
  //
  cout<<"np_sim       "<<np_sim<<endl
      <<"eMin_sim     "<<eMin_sim<<endl
      <<"eMax_sim     "<<eMax_sim<<endl
      <<"integral_sim "<<integral_sim<<endl;
  //
  //
  load_data_file("./data/LST_AdvCam_Zenith_20.00deg_sim.csv", gr_sim_LST_AdvCam, np_sim_LST_AdvCam, eMin_sim_LST_AdvCam, eMax_sim_LST_AdvCam, integral_sim_LST_AdvCam);
  cout<<"np_sim_LST_AdvCam       "<<np_sim_LST_AdvCam<<endl
      <<"eMin_sim_LST_AdvCam     "<<eMin_sim_LST_AdvCam<<endl
      <<"eMax_sim_LST_AdvCam     "<<eMax_sim_LST_AdvCam<<endl
      <<"integral_sim_LST_AdvCam "<<integral_sim_LST_AdvCam<<endl;
  load_data_file("./data/LST_AdvCam_Zenith_20.00deg_trg.csv", gr_trg_LST_AdvCam, np_trg_LST_AdvCam, eMin_trg_LST_AdvCam, eMax_trg_LST_AdvCam, integral_trg_LST_AdvCam);
  cout<<"np_trg_LST_AdvCam       "<<np_trg_LST_AdvCam<<endl
      <<"eMin_trg_LST_AdvCam     "<<eMin_trg_LST_AdvCam<<endl
      <<"eMax_trg_LST_AdvCam     "<<eMax_trg_LST_AdvCam<<endl
      <<"integral_trg_LST_AdvCam "<<integral_trg_LST_AdvCam<<endl;  
  //
  get_eff(gr_sim_LST_AdvCam, gr_trg_LST_AdvCam, gr_trg_eff_AdvCam, 50, 0.005, 50.0);  
  get_eff_area(gr_sim_LST_AdvCam, gr_trg_LST_AdvCam, gr_trg_effArea_AdvCam, 50, 0.005, 50.0, TMath::Pi()*800.0*800.0);    
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1200,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  c1->SetRightMargin(0.01);
  c1->SetLeftMargin(0.12);
  c1->SetTopMargin(0.01);
  c1->SetBottomMargin(0.08);
  c1->Divide(2,1);
  //
  c1->cd(1);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  gr_sim->SetLineColor(kBlack);
  gr_sim->SetLineWidth(3.0);
  gr_sim->SetMarkerStyle(20);
  //
  gr_trg->SetLineColor(kBlue+2);
  gr_trg->SetLineWidth(3.0);
  gr_trg->SetMarkerColor(kBlue+2);
  //
  gr_sel->SetLineColor(kRed+2);
  gr_sel->SetLineWidth(3.0);
  gr_sel->SetMarkerColor(kRed+2);
  //
  TMultiGraph *mg01 = new TMultiGraph();
  mg01->Add(gr_sim);
  mg01->Add(gr_trg);
  mg01->Add(gr_sel);
  mg01->SetMinimum(1.0e3);
  mg01->SetMaximum(1.0e8);
  mg01->Draw("APL");
  mg01->GetXaxis()->SetTitle("Energy, TeV");
  mg01->GetYaxis()->SetTitle("Number of events");
  //
  TLegend *leg01 = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg01->AddEntry(gr_sim, "Simulated (Zenith 23.63 deg)","l");
  leg01->AddEntry(gr_trg, "Triggered (Zenith 23.63 deg)","l");
  leg01->AddEntry(gr_sel, "Selected  (Zenith 23.63 deg)","l");
  leg01->Draw();
  //
  //
  //
  c1->cd(2);
  //
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  //
  gr_trg_eff->SetLineColor(kBlue+2);
  gr_trg_eff->SetLineWidth(3.0);
  gr_trg_eff->SetMarkerColor(kBlue+2);
  //
  gr_sel_eff->SetLineColor(kRed+2);
  gr_sel_eff->SetLineWidth(3.0);
  gr_sel_eff->SetMarkerColor(kRed+2);
  //
  gr_trg_eff_AdvCam->SetLineColor(kBlack);
  gr_trg_eff_AdvCam->SetLineWidth(3.0);
  gr_trg_eff_AdvCam->SetMarkerColor(kBlack);
  //
  TMultiGraph *mg02 = new TMultiGraph();
  mg02->Add(gr_trg_eff);
  mg02->Add(gr_sel_eff);
  //mg02->Add(gr_trg_eff_AdvCam);
  //mg02->SetMinimum(1.0e-2);
  //mg02->SetMaximum(1.0e+0);
  mg02->Draw("APL");
  mg02->GetXaxis()->SetTitle("Energy, TeV");
  mg02->GetYaxis()->SetTitle("Efficiency");
  //
  //
  gr_trg_eff->SetMinimum(1.0e-2);
  gr_trg_eff->SetMaximum(1.0e0);
  gr_trg_eff->Draw("APL");
  mg02->Draw("APL"); 
  //
  gr_trg_eff->GetXaxis()->SetTitle("Energy, TeV"); 
  gr_trg_eff->GetYaxis()->SetTitle("Efficiency");
  //
  //mg->GetYaxis()->SetTitle("Relative differential flux");
  //
  TLegend *leg02 = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg02->AddEntry(gr_trg_eff, "Triggered (Zenith 23.63 deg)","l");
  leg02->AddEntry(gr_sel_eff, "Selected  (Zenith 23.63 deg)","l");
  leg02->Draw();
  //
  //
  //
  //
  //
  //

  //
  TCanvas *c2 = new TCanvas("c2","c2",10,10,1200,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  c2->SetRightMargin(0.01);
  c2->SetLeftMargin(0.12);
  c2->SetTopMargin(0.01);
  c2->SetBottomMargin(0.08);
  c2->Divide(2,1);
  //
  c2->cd(1);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  gr_sim_LST_AdvCam->SetLineColor(kBlack);
  gr_sim_LST_AdvCam->SetLineWidth(3.0);
  gr_sim_LST_AdvCam->SetMarkerStyle(20);
  //
  gr_trg_LST_AdvCam->SetLineColor(kBlue+2);
  gr_trg_LST_AdvCam->SetLineWidth(3.0);
  gr_trg_LST_AdvCam->SetMarkerColor(kBlue+2);
  //
  TMultiGraph *mg03 = new TMultiGraph();
  mg03->Add(gr_sim_LST_AdvCam);
  mg03->Add(gr_trg_LST_AdvCam);
  //mg03->SetMinimum(1.0e3);
  //mg03->SetMaximum(1.0e8);
  mg03->Draw("APL");
  mg03->GetXaxis()->SetTitle("Energy, TeV");
  mg03->GetYaxis()->SetTitle("Number of events");
  //
  TLegend *leg03 = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg03->AddEntry(gr_sim_LST_AdvCam, "Simulated LST AdvCam. (Zenith 20.00 deg)","l");
  leg03->AddEntry(gr_trg_LST_AdvCam, "Triggered LST AdvCam. (Zenith 20.00 deg)","l");
  leg03->Draw();
  //
  //
  //
  c2->cd(2);
  //
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  //
  gr_trg_eff_AdvCam->SetLineColor(kBlack);
  gr_trg_eff_AdvCam->SetLineWidth(3.0);
  gr_trg_eff_AdvCam->SetMarkerColor(kBlack);
  //
  TMultiGraph *mg04 = new TMultiGraph();
  mg04->Add(gr_trg_eff_AdvCam);
  mg04->Add(gr_trg_eff);
  //mg04->SetMinimum(1.0e-2);
  //mg04->SetMaximum(1.0e+0);
  mg04->Draw("APL");
  mg04->GetXaxis()->SetTitle("Energy, TeV");
  mg04->GetYaxis()->SetTitle("Efficiency");
  //
  //
  //gr_trg_eff->SetMinimum(1.0e-2);
  //gr_trg_eff->SetMaximum(1.0e0);
  //gr_trg_eff->Draw("APL");
  mg04->Draw("APL"); 
  //
  mg04->GetXaxis()->SetTitle("Energy, TeV"); 
  mg04->GetYaxis()->SetTitle("Efficiency");
  //
  //mg->GetYaxis()->SetTitle("Relative differential flux");
  //
  TLegend *leg04 = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg04->AddEntry(gr_trg_eff,        "Triggered LST (Zenith 23.63 deg)","l");
  leg04->AddEntry(gr_trg_eff_AdvCam, "Triggered LST AdvCam (Zenith 20.00 deg)","l");
  leg04->Draw();
  //
  //
  //
  TCanvas *c3 = new TCanvas("c3","c3",10,10,1200,600);
  c3->Divide(2,1);
  c3->cd(1);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  c3->SetRightMargin(0.01);
  c3->SetLeftMargin(0.12);
  c3->SetTopMargin(0.01);
  c3->SetBottomMargin(0.08);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  gr_trg_effArea_AdvCam->SetLineColor(kBlack);
  gr_trg_effArea_AdvCam->SetLineWidth(3.0);
  gr_trg_effArea_AdvCam->SetMarkerStyle(20);
  //
  gr_trg_effArea->SetLineColor(kBlue+2);
  gr_trg_effArea->SetLineWidth(3.0);
  gr_trg_effArea->SetMarkerColor(kBlue+2);
  //
  TMultiGraph *mg05 = new TMultiGraph();
  mg05->Add(gr_trg_effArea_AdvCam);
  mg05->Add(gr_trg_effArea);
  mg05->Draw("APL");
  mg05->GetXaxis()->SetTitle("Energy, TeV");
  mg05->GetYaxis()->SetTitle("Trigger effective collection area/m/m");
  //
  //
  c3->cd(2);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  c3->SetRightMargin(0.01);
  c3->SetLeftMargin(0.12);
  c3->SetTopMargin(0.01);
  c3->SetBottomMargin(0.08);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  gPad->SetLogx();
  //
  //
  //
  Double_t eTev_effArea;
  Double_t effArea;
  for(Int_t i = 0; i<gr_trg_effArea->GetN(); i++){
    gr_trg_effArea->GetPoint( i, eTev_effArea, effArea);
    gr_trg_effArea_ratio->SetPoint(i,eTev_effArea, gr_trg_effArea_AdvCam->Eval(eTev_effArea)/effArea);
    gr_trg_effArea_ratio->SetPointError(i,eTev_effArea*0.05, 0.1);
  }
  gr_trg_effArea_ratio->SetLineColor(kBlack);
  gr_trg_effArea_ratio->SetLineWidth(3.0);
  gr_trg_effArea_ratio->SetMarkerStyle(20);
  gr_trg_effArea_ratio->SetMarkerColor(kBlue);
  //
  gr_trg_effArea_ratio->Draw("APL");
  gr_trg_effArea_ratio->GetXaxis()->SetTitle("Energy, TeV");
  gr_trg_effArea_ratio->GetYaxis()->SetTitle("Trigger effective collection area ratio");
  //
  //
  //
  return 0;
}

void load_data_file(TString fileInName, TGraphErrors *gr, Int_t &np, Double_t &eMin, Double_t &eMax, Double_t &integral){
  ifstream fileIn(fileInName.Data());
  string mot;
  Double_t ekin;
  Double_t ftot;
  Int_t npoints_tot = 0;
  integral = 0.0;
  if(fileIn.is_open()){
    while(fileIn>>ekin>>ftot){
      npoints_tot = gr->GetN();
      gr->SetPoint(npoints_tot,ekin,ftot);
      if(gr->GetN()==1)
	eMin = ekin;
      gr->SetPointError(npoints_tot,ekin*0.05,TMath::Sqrt(ftot));
      integral += ftot;
      eMax = ekin;
    }
    fileIn.close();
  }
  else{
    cout<<"Unable to open file \n";
  }
  np=gr->GetN();
}

void load_data_file(TString fileInName, TGraph *gr, Int_t &np, Double_t &eMin, Double_t &eMax, Double_t &integral){
  ifstream fileIn(fileInName.Data());
  string mot;
  Double_t ekin;
  Double_t ftot;
  integral = 0.0;
  if(fileIn.is_open()){
    while(fileIn>>ekin>>ftot){
      gr->SetPoint(gr->GetN(),ekin,ftot);
      integral += ftot;
      if(gr->GetN()==1)
	eMin = ekin;
      eMax = ekin;
    }
    fileIn.close();
  }
  else{
    cout<<"Unable to open file \n";
  }
  np=gr->GetN();
}

void get_eff(TGraph *gr_sim, TGraph *gr_trg, TGraph *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV){
  //
  Double_t log10_eTeV_min = TMath::Log10(eMin_TeV);
  Double_t log10_eTeV_max = TMath::Log10(eMax_TeV);  
  Double_t log10_eTeV_R;
  Double_t log10_eTeV_L;
  Double_t eTeV;
  Double_t F;
  //
  for(Int_t i = 1; i< np;i++){
    log10_eTeV_R = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(np-1)*i;
    log10_eTeV_L = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(np-1)*(i-1);
    eTeV = TMath::Power(10.0,(log10_eTeV_R + log10_eTeV_L)/2.0);
    gr_trg_eff->SetPoint(gr_trg_eff->GetN(),eTeV,gr_trg->Eval(eTeV)/gr_sim->Eval(eTeV));
  }
}

void get_eff(TGraphErrors *gr_sim, TGraphErrors *gr_trg, TGraphErrors *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV){
  //
  Double_t log10_eTeV_min = TMath::Log10(eMin_TeV);
  Double_t log10_eTeV_max = TMath::Log10(eMax_TeV);  
  Double_t log10_eTeV_R;
  Double_t log10_eTeV_L;
  Double_t eTeV;
  Double_t F;
  Double_t eff_val;
  Double_t eff_val_err;
  Int_t npoints_tot = 0;
  //
  for(Int_t i = 1; i< np;i++){
    log10_eTeV_R = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(np-1)*i;
    log10_eTeV_L = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(np-1)*(i-1);
    eTeV = TMath::Power(10.0,(log10_eTeV_R + log10_eTeV_L)/2.0);
    npoints_tot = gr_trg_eff->GetN();
    eff_val=gr_trg->Eval(eTeV)/gr_sim->Eval(eTeV);
    //eff_val_err=TMath::Sqrt(eff_val*(1.0 - eff_val)*gr_sim->Eval(eTeV));
    eff_val_err=2.0/gr_sim->Eval(eTeV)*TMath::Sqrt(gr_trg->Eval(eTeV)*(1.0 - eff_val));
    gr_trg_eff->SetPoint(npoints_tot,eTeV,eff_val);
    gr_trg_eff->SetPointError(npoints_tot,eTeV*0.05,eff_val_err);
  }
}

void get_eff_area(TGraph *gr_sim, TGraph *gr_trg, TGraph *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV, Double_t area){
  //
  Double_t log10_eTeV_min = TMath::Log10(eMin_TeV);
  Double_t log10_eTeV_max = TMath::Log10(eMax_TeV);  
  Double_t log10_eTeV_R;
  Double_t log10_eTeV_L;
  Double_t eTeV;
  Double_t F;
  //
  for(Int_t i = 1; i< np;i++){
    log10_eTeV_R = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(np-1)*i;
    log10_eTeV_L = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(np-1)*(i-1);
    eTeV = TMath::Power(10.0,(log10_eTeV_R + log10_eTeV_L)/2.0);
    gr_trg_eff->SetPoint(gr_trg_eff->GetN(),eTeV,gr_trg->Eval(eTeV)/gr_sim->Eval(eTeV)*area);
  }
}

void get_eff_area(TGraphErrors *gr_sim, TGraphErrors *gr_trg, TGraphErrors *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV, Double_t area){
  //
  Double_t log10_eTeV_min = TMath::Log10(eMin_TeV);
  Double_t log10_eTeV_max = TMath::Log10(eMax_TeV);  
  Double_t log10_eTeV_R;
  Double_t log10_eTeV_L;
  Double_t eTeV;
  Double_t F;
  Double_t eff_val;
  Double_t eff_val_err;
  Int_t npoints_tot = 0;
  //
  for(Int_t i = 1; i< np;i++){
    log10_eTeV_R = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(np-1)*i;
    log10_eTeV_L = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(np-1)*(i-1);
    eTeV = TMath::Power(10.0,(log10_eTeV_R + log10_eTeV_L)/2.0);
    npoints_tot = gr_trg_eff->GetN();
    eff_val=gr_trg->Eval(eTeV)/gr_sim->Eval(eTeV);
    //eff_val_err=TMath::Sqrt(eff_val*(1.0 - eff_val)*gr_sim->Eval(eTeV));
    eff_val_err=2.0/gr_sim->Eval(eTeV)*TMath::Sqrt(gr_trg->Eval(eTeV)*(1.0 - eff_val));
    gr_trg_eff->SetPoint(npoints_tot,eTeV,eff_val*area);
    gr_trg_eff->SetPointError(npoints_tot,eTeV*0.05,eff_val_err*area);
  }
}
