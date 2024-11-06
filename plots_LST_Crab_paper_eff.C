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
#include <vector>
#include <assert.h>

using namespace std;

void load_data_file(TString fileInName, TGraphErrors *gr, Int_t &np, Double_t &eMin, Double_t &eMax, Double_t &integral);
void load_data_file(TString fileInName, TGraph *gr, Int_t &np, Double_t &eMin, Double_t &eMax, Double_t &integral);
void get_eff(TGraph *gr_sim, TGraph *gr_trg, TGraph *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV);
void get_eff(TGraphErrors *gr_sim, TGraphErrors *gr_trg, TGraphErrors *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV);
void get_eff_area(TGraph *gr_sim, TGraph *gr_trg, TGraph *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV, Double_t area);
void get_eff_area(TGraphErrors *gr_sim, TGraphErrors *gr_trg, TGraphErrors *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV, Double_t area);
Double_t modiﬁed_log_parabola(Double_t eGeV);
Double_t get_crab_rate(TGraphErrors *gr_trg_effArea_msq_vs_TeV, Double_t eMin_GeV, Double_t eMax_GeV);
Double_t get_crab_rate(TGraph *gr_trg_effArea_msq_vs_TeV, Double_t eMin_GeV, Double_t eMax_GeV);
void get_crab_gamma_rate(TGraphErrors *gr_gamma_rate_LST, TGraph *gr_trg_effArea_msq_vs_TeV, TGraphErrors *gr_gamma_rate_LST_AdvCam);
//
void fill_log_test(TGraph *gr, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV, Double_t AA, Double_t BB, Double_t CC, Double_t DD);
void logalyze(TGraphErrors *gr, TGraphErrors *gr_log);
void logalyze(TGraph *gr, TGraphErrors *gr_log);
//
Double_t eff_area_f(Double_t *x, Double_t *par);

template <class T>
void save_gr_to_csv( T *gr_wf, TString csv_file_out, bool csv_format = true);

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
  TGraphErrors *gr_trg_eff_log = new TGraphErrors();
  TGraph *gr_sel_eff = new TGraph();
  TGraph *gr_trg_effArea = new TGraph();
  TGraphErrors *gr_trg_effArea_log = new TGraphErrors();
  //
  Int_t np_sim_LST_AdvCam, np_trg_LST_AdvCam;
  Double_t eMin_sim_LST_AdvCam, eMax_sim_LST_AdvCam, integral_sim_LST_AdvCam;
  Double_t eMin_trg_LST_AdvCam, eMax_trg_LST_AdvCam, integral_trg_LST_AdvCam;
  TGraphErrors *gr_sim_LST_AdvCam = new TGraphErrors();
  TGraphErrors *gr_trg_LST_AdvCam = new TGraphErrors();
  TGraphErrors *gr_trg_eff_AdvCam = new TGraphErrors();
  TGraphErrors *gr_trg_effArea_AdvCam = new TGraphErrors();
  TGraphErrors *gr_trg_effArea_AdvCam_log = new TGraphErrors();
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
  //
  leg04->Draw();
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
  //
  //
  TCanvas *c4 = new TCanvas("c4","c4",10,10,1200,600);
  c4->Divide(2,1);
  c4->cd(1);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  c4->SetRightMargin(0.01);
  c4->SetLeftMargin(0.12);
  c4->SetTopMargin(0.01);
  c4->SetBottomMargin(0.08);
  //
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  //
  //
  TGraph *gr_trg_effArea_fit = new TGraph();
  fill_log_test(gr_trg_effArea_fit, 100, 0.02, 17.0, 5.99535e+00, 9.22654e+01, 4.68694e+00, -3.70726e+00);
  //  
  gr_trg_effArea_fit->SetLineColor(kRed+2);
  gr_trg_effArea_fit->SetLineWidth(3.0);
  gr_trg_effArea_fit->SetMarkerColor(kRed+2);
  //

  //
  //logalyze( gr_trg_effArea, gr_trg_effArea_log);
  //
  //Fit
  //const Int_t npar = 4;
  //Double_t inParameters[npar];
  //inParameters[0] = 5.99535e+00;
  //inParameters[1] = 9.22654e+01;
  //inParameters[2] = 4.68694e+00;
  //inParameters[3] = -3.70726e+00;
  //
  //TF1 *f_eff_area_f = new TF1("f_eff_area_f", eff_area_f, -1.6, 1.2, npar);
  //f_eff_area_f->SetParameters(inParameters);
  //f_eff_area_f->FixParameter(0,inParameters[0]);
  //f_eff_area_f->FixParameter(1,inParameters[1]);
  //f_eff_area_f->FixParameter(2,inParameters[2]);
  //f_eff_area_f->FixParameter(3,inParameters[3]);  
  //
  //gr_trg_effArea_log->Fit("f_eff_area_f","","",-1.6, 1.2);
  //gr_trg_effArea_log->Draw("APL");
  //
  TMultiGraph *mg006 = new TMultiGraph();
  mg006->Add(gr_trg_effArea);
  mg006->Add(gr_trg_effArea_fit);
  //
  mg006->Draw("APL");
  mg006->GetXaxis()->SetTitle("Energy, TeV");
  mg006->GetYaxis()->SetTitle("Trigger effective collection area/m/m");
  //
  TLegend *leg005 = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg005->AddEntry(gr_trg_effArea,     "Triggered LST PMT (Zenith 23.63 deg) fit.","l");
  leg005->AddEntry(gr_trg_effArea_fit, "Triggered LST PMT (Zenith 23.63 deg) sim.","l");
  leg005->Draw();

  
  //
  //
  //
  c4->cd(2);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  c4->SetRightMargin(0.01);
  c4->SetLeftMargin(0.12);
  c4->SetTopMargin(0.01);
  c4->SetBottomMargin(0.08);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  TGraph *gr_log = new TGraph();
  fill_log_test(gr_log, 100, 0.005, 50.0, 6.01082e+00, 8.43731e+00, 3.88802e+00,-2.54485e+00);
  //
  //logalyze(gr_trg_effArea_AdvCam, gr_trg_effArea_AdvCam_log);
  //
  //Fit
  //const Int_t npar = 4;
  //Double_t inParameters[npar];
  //inParameters[0] = 5.99534e+00;
  //inParameters[1] = 1.65013e+01;
  //inParameters[2] = 4.19283e+00;
  //inParameters[3] = -2.92505e+00;
  //
  //TF1 *f_eff_area_f = new TF1("f_eff_area_f", eff_area_f, -2.2, 2.0, npar);
  //f_eff_area_f->SetParameters(inParameters);
  //f_eff_area_f->FixParameter(0,inParameters[0]);
  //f_eff_area_f->FixParameter(1,inParameters[1]);
  //f_eff_area_f->FixParameter(2,inParameters[2]);
  //f_eff_area_f->FixParameter(3,inParameters[3]);  
  //gr_trg_effArea_AdvCam_log->Fit("f_eff_area_f","","",-2.3, 2.0);
  //
  TMultiGraph *mg06 = new TMultiGraph();
  //mg06->Add(gr_trg_effArea_AdvCam_log);
  mg06->Add(gr_trg_effArea_AdvCam);
  gr_log->SetLineColor(kRed+2);
  gr_log->SetLineWidth(3.0);
  gr_log->SetMarkerColor(kRed+2);
  //
  mg06->Add(gr_log);
  mg06->Draw("APL");
  mg06->GetXaxis()->SetTitle("Energy, TeV");
  mg06->GetYaxis()->SetTitle("Trigger effective collection area/m/m");
  //
  TLegend *leg05 = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg05->AddEntry(gr_log,            "Triggered LST AdvCam (Zenith 20.00 deg) fit.","l");
  leg05->AddEntry(gr_trg_eff_AdvCam, "Triggered LST AdvCam (Zenith 20.00 deg) sim.","l");
  leg05->Draw();
  //
  //
  //
  cout<<get_crab_rate(gr_trg_effArea_AdvCam, 5.0, 50000.0)<<endl
      <<get_crab_rate(gr_log, 5.0, 50000.0)<<endl
      <<get_crab_rate(gr_trg_effArea, 20.0, 1000.0)<<endl;
  cout<<get_crab_rate(gr_trg_effArea_AdvCam, 5.0, 300.0)<<endl
      <<get_crab_rate(gr_log, 5.0, 300.0)<<endl;
  //
  //
  TCanvas *c5 = new TCanvas("c5","c5",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  c5->SetRightMargin(0.01);
  c5->SetLeftMargin(0.12);
  c5->SetTopMargin(0.01);
  c5->SetBottomMargin(0.08);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  //gPad->SetLogx();  
  //
  //
  //
  TGraphErrors *gr_gamma_rates_trg_LST = new TGraphErrors();
  TGraphErrors *gr_gamma_rates_trg_LST_AdvCam = new TGraphErrors();
  Int_t np_gr_trg_LST;
  Double_t eMin_gr_trg_LST, eMax_gr_trg_LST, integral_gr_trg_LST;
  //
  load_data_file("./data/LST_Crab_paper_gamma_rate_trg.csv", gr_gamma_rates_trg_LST, np_gr_trg_LST, eMin_gr_trg_LST, eMax_gr_trg_LST, integral_gr_trg_LST);
  cout<<"np_gr_trg_LST       "<<np_gr_trg_LST<<endl
      <<"eMin_gr_trg_LST     "<<eMin_gr_trg_LST<<endl
      <<"eMax_gr_trg_LST     "<<eMax_gr_trg_LST<<endl
      <<"integral_gr_trg_LST "<<integral_gr_trg_LST<<endl;
  //
  //
  //
  get_crab_gamma_rate( gr_gamma_rates_trg_LST, gr_log, gr_gamma_rates_trg_LST_AdvCam);
  TMultiGraph *mg07 = new TMultiGraph();
  mg07->Add(gr_gamma_rates_trg_LST);
  mg07->Add(gr_gamma_rates_trg_LST_AdvCam);
  mg07->Draw("APL");
  mg07->GetXaxis()->SetTitle("Energy, TeV");
  mg07->GetYaxis()->SetTitle("Gamma rate, 1/s 1/TeV");
  //
  mg07->SetMinimum(0.0);
  mg07->SetMaximum(160);
  //
  gr_gamma_rates_trg_LST->SetLineColor(kBlack);
  gr_gamma_rates_trg_LST->SetLineWidth(3.0);
  gr_gamma_rates_trg_LST->SetMarkerColor(kBlack);
  //
  gr_gamma_rates_trg_LST_AdvCam->SetLineColor(kRed + 2);
  gr_gamma_rates_trg_LST_AdvCam->SetLineWidth(3.0);
  gr_gamma_rates_trg_LST_AdvCam->SetMarkerColor(kRed + 2.0);  
  //
  TLegend *leg06 = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg06->AddEntry(gr_gamma_rates_trg_LST,        "Triggered LST PMT cam. (Zenith 10.00 deg) 4.35 1/s","l");
  leg06->AddEntry(gr_gamma_rates_trg_LST_AdvCam, "Triggered LST Adv cam. (Zenith 20.00 deg) 7.10 1/s","l");
  leg06->Draw();
  //
  //  
  //
  TCanvas *c6 = new TCanvas("c6","c6",10,10,1200,600);
  c6->Divide(2,1);
  c6->cd(1);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  c6->SetRightMargin(0.01);
  c6->SetLeftMargin(0.12);
  c6->SetTopMargin(0.01);
  c6->SetBottomMargin(0.08);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  TGraph *gr_trg_effArea_AdvCam_fit = new TGraph();
  TGraph *gr_trg_effArea_AdvCam_vs_PMT_fit = new TGraph();
  fill_log_test(gr_trg_effArea_AdvCam_fit, 100, 0.005, 50.0, 6.01082e+00, 8.43731e+00, 3.88802e+00,-2.54485e+00);
  //
  gr_trg_effArea_AdvCam_fit->SetLineColor(kBlack);
  gr_trg_effArea_AdvCam_fit->SetLineWidth(3.0);
  gr_trg_effArea_AdvCam_fit->SetMarkerColor(kBlack);
  //
  TMultiGraph *mg08 = new TMultiGraph();
  mg08->Add(gr_trg_effArea_fit);
  mg08->Add(gr_trg_effArea_AdvCam_fit);
  mg08->Draw("APL");
  mg08->GetXaxis()->SetTitle("Energy, TeV");
  mg08->GetYaxis()->SetTitle("Trigger effective collection area, m^{2}");
  //
  //
  TLegend *leg07 = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg07->AddEntry(gr_trg_effArea_fit,        "Triggered LST PMT cam. (Zenith 23.63 deg) fit","l");
  leg07->AddEntry(gr_trg_effArea_AdvCam_fit, "Triggered LST Adv cam. (Zenith 20.00 deg) fit","l");
  leg07->Draw();
  //
  //
  c6->cd(2);
  //
  //
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  gPad->SetLogx();
  //
  for(Int_t i = 0; i<gr_trg_effArea->GetN(); i++){
    gr_trg_effArea->GetPoint( i, eTev_effArea, effArea);    
    gr_trg_effArea_AdvCam_vs_PMT_fit->SetPoint(i,eTev_effArea, gr_trg_effArea_AdvCam_fit->Eval(eTev_effArea)/gr_trg_effArea_fit->Eval(eTev_effArea));
  }
  gr_trg_effArea_AdvCam_vs_PMT_fit->SetLineColor(kBlack);
  gr_trg_effArea_AdvCam_vs_PMT_fit->SetLineWidth(3.0);
  gr_trg_effArea_AdvCam_vs_PMT_fit->SetMarkerStyle(20);
  gr_trg_effArea_AdvCam_vs_PMT_fit->SetMarkerColor(kBlue);
  //
  gr_trg_effArea_AdvCam_vs_PMT_fit->Draw("APL");
  gr_trg_effArea_AdvCam_vs_PMT_fit->GetXaxis()->SetTitle("Energy, TeV");
  gr_trg_effArea_AdvCam_vs_PMT_fit->GetYaxis()->SetTitle("Trigger effective collection area ratio"); 
  //
  //
  //
  //
  //
  save_gr_to_csv<TGraph>(gr_trg_effArea_fit, "gr_trg_effArea_fit.csv");
  save_gr_to_csv<TGraph>(gr_trg_effArea_AdvCam_fit, "gr_trg_effArea_AdvCam_fit.csv");
  save_gr_to_csv<TGraph>(gr_trg_effArea_AdvCam_vs_PMT_fit, "gr_trg_effArea_AdvCam_vs_PMT_fit.csv");
  save_gr_to_csv<TGraphErrors>(gr_gamma_rates_trg_LST,"gr_gamma_rates_trg_LST.csv");
  save_gr_to_csv<TGraphErrors>(gr_gamma_rates_trg_LST_AdvCam,"gr_gamma_rates_trg_LST_AdvCam.csv");
  //  
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

Double_t modiﬁed_log_parabola(Double_t eGeV){
  Double_t A    = 2.97375e-11;
  Double_t B    = 1.00000e+03;
  Double_t C    = 2.73558e+00;
  Double_t D    = 1.82440e-01;
  Double_t Emax = 2.20294e+01;
  Double_t aa   = 9.94035e-01;
  Double_t logE = TMath::Abs(TMath::Log10(TMath::Abs(eGeV/B/Emax)));
  Double_t val  = A*(TMath::Power(eGeV/B, -C + D*TMath::Power(logE,aa)))/1000.0;
  return val;
}

Double_t modiﬁed_log_parabola_TeV(Double_t eTeV){
  Double_t A    = 2.97375e-11;
  Double_t B    = 1.00000e+00;
  Double_t C    = 2.73558e+00;
  Double_t D    = 1.82440e-01;
  Double_t Emax = 2.20294e+01;
  Double_t aa   = 9.94035e-01;
  Double_t logE = TMath::Abs(TMath::Log10(TMath::Abs(eTeV/B/Emax)));
  Double_t val  = A*(TMath::Power(eTeV/B, -C + D*TMath::Power(logE,aa)));
  return val;
}

Double_t get_crab_rate(TGraphErrors *gr_trg_effArea_msq_vs_TeV, Double_t eMin_GeV, Double_t eMax_GeV){
  Int_t nn = 10000000;
  Double_t ee;
  Double_t dx = (eMax_GeV - eMin_GeV)/nn;
  Double_t dx_half = dx/2.0;
  Double_t integral = 0.0;
  for(Int_t i = 0; i<nn; i++){
    ee = eMin_GeV + dx_half + dx*i;
    integral+=modiﬁed_log_parabola(ee)*gr_trg_effArea_msq_vs_TeV->Eval(ee/1000.0)*10000.0*dx;
  }
  return integral;
}

Double_t get_crab_rate(TGraph *gr_trg_effArea_msq_vs_TeV, Double_t eMin_GeV, Double_t eMax_GeV){
  Int_t nn = 10000000;
  Double_t ee;
  Double_t dx = (eMax_GeV - eMin_GeV)/nn;
  Double_t dx_half = dx/2.0;
  Double_t integral = 0.0;
  for(Int_t i = 0; i<nn; i++){
    ee = eMin_GeV + dx_half + dx*i;
    integral+=modiﬁed_log_parabola(ee)*gr_trg_effArea_msq_vs_TeV->Eval(ee/1000.0)*10000.0*dx;
  }
  return integral;
}

void fill_log_test(TGraph *gr, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV, Double_t AA, Double_t BB, Double_t CC, Double_t DD){
  //
  Double_t log10_eTeV_min = TMath::Log10(eMin_TeV);
  Double_t log10_eTeV_max = TMath::Log10(eMax_TeV);  
  Double_t log10_eTeV_R;
  Double_t log10_eTeV_L;
  //Double_t eTeV, eTeV_L, eTeV_R;
  Double_t eTeV;
  Double_t val;
  Double_t val_err;
  Int_t npoints_tot = 0;
  //
  for(Int_t i = 1; i< np;i++){
    //
    //
    //eTeV_R = eMin_TeV + (eMax_TeV-eMin_TeV)/(np-1)*i;
    //eTeV_L = eMin_TeV + (eMax_TeV-eMin_TeV)/(np-1)*(i-1);
    //eTeV = (eTeV_R + eTeV_L)/2.0;
    //val=(6.01082e+00 - 8.43731e+00*(TMath::Power((eTeV + 3.88802e+00),-2.54485e+00)));
    //
    //
    log10_eTeV_R = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(np-1)*i;
    log10_eTeV_L = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(np-1)*(i-1);
    eTeV = TMath::Power(10.0,(log10_eTeV_R + log10_eTeV_L)/2.0);
    //val=(6.01082e+00 - 8.43731e+00*(TMath::Power((TMath::Log10(eTeV) + 3.88802e+00),-2.54485e+00)));
    //val=(-1.97117e-01 - 5.95699e+00*(TMath::Power((TMath::Log10(eTeV) + 3.46346e+00),-2.22730e+00)));    
    val=(AA-BB*(TMath::Power((TMath::Log10(eTeV) + CC),DD)));
    val=TMath::Power(10,val);
    //
    //val=1.0e+4*TMath::Power(TMath::Log10(TMath::Power(eTeV*1000.0*0.2,1.0)),4.0);
    //val=1.0e+5*TMath::Power(eTeV*100.0,0.7)*TMath::Log10(eTeV*1000.0*0.2)*TMath::Log10(eTeV*1000.0*0.2);
    //val=1.0e+0*TMath::Power(10,1.0e+0*TMath::Log10(eTeV*1000.0) - 1.0e+0*TMath::Power(TMath::Log10(1.0 + eTeV),-2.0));
    //val=3.66 + 0.536*TMath::Log10(eTeV*1000) - TMath::Log10(1.0 + TMath::Power(eTeV/0.128,-2.07));
    //val=3.66 + 0.536*TMath::Log10(eTeV*1000) - 0.8e+0*TMath::Log10(1.0 + TMath::Power(eTeV/0.0128,-3.0));
    //val=1.0e+0/(51.0 - TMath::Power(TMath::Log10(eTeV*1000), 1.0));
    //val=3.5*TMath::Power((eTeV + 2.5),0.25) + 1.0;
    //eff_val_err=TMath::Sqrt(eff_val*(1.0 - eff_val)*gr_sim->Eval(eTeV));
    //val_err=2.0/gr_sim->Eval(eTeV)*TMath::Sqrt(gr_trg->Eval(eTeV)*(1.0 - eff_val));
    //
    npoints_tot = gr->GetN();
    gr->SetPoint(npoints_tot,eTeV,val);
  }
}

void logalyze(TGraphErrors *gr, TGraphErrors *gr_log){
  Double_t xx, yy;
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i,xx,yy);
    gr_log->SetPoint(i,TMath::Log10(xx),TMath::Log10(yy));
    gr_log->SetPointError(i,0.01,0.05);
  }
}

void logalyze(TGraph *gr, TGraphErrors *gr_log){
  Double_t xx, yy;
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i,xx,yy);
    gr_log->SetPoint(i,TMath::Log10(xx),TMath::Log10(yy));
    gr_log->SetPointError( i, 0.01, 0.03);
    //if(TMath::Log10(xx)<0.0)
    //gr_log->SetPointError(i, 0.01,0.03);
    //else
    //gr_log->SetPointError(i, 0.01,0.01);
  }
}

Double_t eff_area_f(Double_t *x, Double_t *par){
  Double_t A  = par[0];
  Double_t B  = par[1];
  Double_t C  = par[2];
  Double_t D  = par[3];
  Double_t val = A - B*TMath::Power((x[0] + C),D);
  return val;
}

void get_crab_gamma_rate(TGraphErrors *gr_gamma_rate_LST, TGraph *gr_trg_effArea_msq_vs_TeV, TGraphErrors *gr_gamma_rate_LST_AdvCam){
  //
  Double_t xx, yy;  
  for( Int_t i = 0; i < gr_gamma_rate_LST->GetN(); i++){
    //
    gr_gamma_rate_LST->GetPoint(i,xx,yy);
    gr_gamma_rate_LST->SetPointError(i,0.0,0.0);
    gr_gamma_rate_LST_AdvCam->SetPoint(i,xx,gr_trg_effArea_msq_vs_TeV->Eval(xx)*10000.0*modiﬁed_log_parabola_TeV(xx));    
  }
}

template <class T>
void save_gr_to_csv( T *gr_wf, TString csv_file_out, bool csv_format){
  //if (typeid(T) == typeid(TGraph))
  //cout<<"T = TGraph"<<endl;
  //if (typeid(T) == typeid(TGraphErrors))
  //cout<<"T = TGraphErrors"<<endl;
  Double_t xx;
  Double_t yy;
  ofstream csvfile;
  csvfile.open (csv_file_out.Data());
  if(csv_format){
    csvfile<<"xx,yy"<<std::endl;
    for(Int_t i = 0;i<gr_wf->GetN();i++){
      gr_wf->GetPoint(i,xx,yy);
      csvfile<<xx<<","<<yy<<std::endl;
    }
  }
  else{
    for(Int_t i = 0;i<gr_wf->GetN();i++){
      gr_wf->GetPoint(i,xx,yy);
      csvfile<<xx<<" "<<yy<<std::endl;
    }
  }
  csvfile.close();
}
