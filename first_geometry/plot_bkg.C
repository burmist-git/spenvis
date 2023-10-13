//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>

#include <time.h>

using namespace std;

void plot_bkg_p(TGraph *gr_nPhot, TGraph *gr_nPhot_more, TGraph *gr_npe, TGraph *gr_dedx);
void plot_bkg_e(TGraph *gr_nPhot, TGraph *gr_nPhot_more, TGraph *gr_npe, TGraph *gr_dedx);
void load_data_file(TString fileInName, TGraph *gr_tot, TGraph *gr_diff, Double_t &hMin, Double_t &hMax);
Double_t getProbOfmoreThenThreshold( TH1D *h1, Double_t thresh);

Int_t plot_bkg(){

  //
  TGraph *gr_nPhot_p = new TGraph();
  TGraph *gr_nPhot_more_p = new TGraph();
  TGraph *gr_npe_p = new TGraph();
  TGraph *gr_dedx_p = new TGraph();
  //
  TGraph *gr_nPhot_e = new TGraph();
  TGraph *gr_nPhot_more_e = new TGraph();
  TGraph *gr_npe_e = new TGraph();
  TGraph *gr_dedx_e = new TGraph();

  //
  plot_bkg_p(gr_nPhot_p, gr_nPhot_more_p, gr_npe_p, gr_dedx_p);
  plot_bkg_e(gr_nPhot_e, gr_nPhot_more_e, gr_npe_e, gr_dedx_e);

  //  
  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  //gStyle->SetOptStat(kFALSE);
  //
  c1->SetRightMargin(0.01);
  c1->SetLeftMargin(0.12);
  c1->SetTopMargin(0.01);
  c1->SetBottomMargin(0.1);
  //
  gPad->SetGridx();
  gPad->SetGridy();

  //
  gr_npe_e->SetMarkerStyle(20);
  gr_npe_e->SetLineColor(kBlack);
  gr_npe_e->SetLineWidth(3.0);
  gr_npe_e->SetMarkerColor(kBlack);
  gr_npe_e->GetXaxis()->SetTitle("Energy, MeV");
  gr_npe_e->GetYaxis()->SetTitle("# p.e. per particle");
  //
  gr_nPhot_e->SetMarkerStyle(20);
  gr_nPhot_e->SetLineColor(kBlack);
  gr_nPhot_e->SetLineWidth(3.0);
  gr_nPhot_e->SetMarkerColor(kBlack);
  gr_nPhot_e->GetXaxis()->SetTitle("Energy, MeV");
  gr_nPhot_e->GetYaxis()->SetTitle("# of photons per particle");
  //
  gr_dedx_e->SetMarkerStyle(20);
  gr_dedx_e->SetLineColor(kBlack);
  gr_dedx_e->SetLineWidth(3.0);
  gr_dedx_e->SetMarkerColor(kBlack);
  gr_dedx_e->GetXaxis()->SetTitle("Energy, MeV");
  gr_dedx_e->GetYaxis()->SetTitle("dE/dx per particle, MeV");

  //
  gr_npe_p->SetMarkerStyle(20);
  gr_npe_p->SetLineColor(kBlack);
  gr_npe_p->SetLineWidth(3.0);
  gr_npe_p->SetMarkerColor(kBlack);
  gr_npe_p->GetXaxis()->SetTitle("Energy, MeV");
  gr_npe_p->GetYaxis()->SetTitle("# p.e. per particle");
  //
  gr_nPhot_p->SetMarkerStyle(20);
  gr_nPhot_p->SetLineColor(kBlack);
  gr_nPhot_p->SetLineWidth(3.0);
  gr_nPhot_p->SetMarkerColor(kBlack);
  gr_nPhot_p->GetXaxis()->SetTitle("Energy, MeV");
  gr_nPhot_p->GetYaxis()->SetTitle("# of photons per particle");
  //
  gr_dedx_p->SetMarkerStyle(20);
  gr_dedx_p->SetLineColor(kBlack);
  gr_dedx_p->SetLineWidth(3.0);
  gr_dedx_p->SetMarkerColor(kBlack);
  gr_dedx_p->GetXaxis()->SetTitle("Energy, MeV");
  gr_dedx_p->GetYaxis()->SetTitle("dE/dx per particle, MeV");
  
  //gr_nPhot_e->Draw("APL");
  //gr_npe_e->Draw("APL");
  //gr_dedx_e->Draw("APL");
  
  //gr_nPhot_p->Draw("APL");
  //gr_npe_p->Draw("APL");
  //gr_dedx_p->SetMinimum(1.0e-10);
  //gr_dedx_p->Draw("APL");

  ////////////////////////////////////
  TFile *f1_e = new TFile("./hist_root_bkg/hist_e-_Ekin_7.00MeV.dat.root");
  TH1D *h1_dedx_sipm_zoom_zoom_e = (TH1D*)f1_e->Get("h1_dedx_sipm_zoom_zoom");
  TH1D *h1_nPhot_e = (TH1D*)f1_e->Get("h1_nPhot");
  //
  TFile *f1_p = new TFile("./hist_root_bkg/hist_proton_Ekin_300.0MeV.dat.root");
  TH1D *h1_dedx_sipm_zoom_zoom_p = (TH1D*)f1_p->Get("h1_dedx_sipm_zoom_zoom");
  TH1D *h1_nPhot_p = (TH1D*)f1_p->Get("h1_nPhot");
  //
  h1_dedx_sipm_zoom_zoom_e->SetLineColor(kBlack);
  h1_dedx_sipm_zoom_zoom_e->SetLineWidth(3.0);
  //
  h1_nPhot_e->SetLineColor(kBlack);
  h1_nPhot_e->SetLineWidth(3.0);
  //
  h1_dedx_sipm_zoom_zoom_p->SetLineColor(kBlack);
  h1_dedx_sipm_zoom_zoom_p->SetLineWidth(3.0);
  //
  h1_nPhot_p->SetLineColor(kBlack);
  h1_nPhot_p->SetLineWidth(3.0);
  //
  h1_dedx_sipm_zoom_zoom_e->Rebin(30);
  h1_dedx_sipm_zoom_zoom_p->Rebin(30);
  //
  h1_dedx_sipm_zoom_zoom_e->SetTitle("");
  h1_dedx_sipm_zoom_zoom_e->GetXaxis()->SetTitle("dE/dx per electron, MeV");
  //
  h1_nPhot_e->SetTitle("");
  h1_nPhot_e->GetXaxis()->SetTitle("# of photons per electron");
  //
  h1_dedx_sipm_zoom_zoom_p->SetTitle("");
  h1_dedx_sipm_zoom_zoom_p->GetXaxis()->SetTitle("dE/dx per proton, MeV");
  //
  h1_nPhot_p->SetTitle("");
  h1_nPhot_p->GetXaxis()->SetTitle("# of photons per proton");
  //
  //
  //
  //h1_dedx_sipm_zoom_zoom_e->Draw();
  //h1_nPhot_e->Draw();
  //h1_dedx_sipm_zoom_zoom_p->Draw();
  //h1_nPhot_p->Draw();
  //
  ////////////////////////////////////
  //
  TGraph *gr_tot_e = new TGraph();
  TGraph *gr_diff_e = new TGraph();
  TGraph *gr_tot_p = new TGraph();
  TGraph *gr_diff_p = new TGraph();
  Double_t hMin;
  Double_t hMax;
  load_data_file("../data/550_km_e.dat", gr_tot_e, gr_diff_e, hMin, hMax);
  load_data_file("../data/550_km_p.dat", gr_tot_p, gr_diff_p, hMin, hMax);
  //
  //gr_diff_e->Draw("APL");
  //gr_diff_p->Draw("APL");
  //
  //
  TGraph *gr_dedx_m_diff_e  = new TGraph();
  TGraph *gr_dedx_m_diff_p  = new TGraph();
  TGraph *gr_nPhot_m_diff_e = new TGraph();
  TGraph *gr_nPhot_m_diff_p = new TGraph();
  //
  //TGraph *gr_nPhot_more_m_diff_e = new TGraph();
  //TGraph *gr_nPhot_more_m_diff_p = new TGraph();  
  //
  Int_t i = 0;
  Int_t n = 1000;
  //
  //gr_diff_e
  //gr_dedx_e
  Double_t eMin_e = 0.2;
  Double_t eMax_e = 7.0;
  Double_t e = 7.0;
  Double_t R_sim = 50.0;//cm
  Double_t sTot = 4*TMath::Pi()*R_sim*R_sim;
  Double_t tTot = 24.0*60.0*60.0*365.0*3;
  Double_t volumeSiPM = 10*0.05*2.5*2.5; //cm^3
  Double_t rho = 2.7/1000.0;             //kg/cm^3
  Double_t tot_mass = volumeSiPM*rho;    //kg
  // 1 MeV = 1.60218e-13 J
  Double_t Mev_to_J = 1.60218e-13;       //J/MeV
  Double_t gy_per_MeV;
  Double_t photon_per_MeV;
  Double_t photon_more_per_MeV;
  Double_t dose_e = 0.0;
  Double_t dose_p = 0.0;
  Double_t photon_rate_e = 0.0;
  Double_t photon_rate_p = 0.0;
  Double_t photon_more_rate_e = 0.0;
  Double_t photon_more_rate_p = 0.0;
  for(i = 0;i<n;i++){
    e = (eMax_e - eMin_e)/(n-1)*i + eMin_e;
    gy_per_MeV = gr_diff_e->Eval(e)*gr_dedx_e->Eval(e)*sTot*tTot*Mev_to_J/tot_mass;
    photon_per_MeV = gr_diff_e->Eval(e)*gr_nPhot_e->Eval(e)*sTot;
    photon_more_per_MeV = gr_diff_e->Eval(e)*gr_nPhot_more_e->Eval(e)*sTot;
    gr_dedx_m_diff_e->SetPoint(gr_dedx_m_diff_e->GetN(),e,gy_per_MeV);
    dose_e += gy_per_MeV*(eMax_e - eMin_e)/(n-1);
    photon_rate_e += photon_per_MeV*(eMax_e - eMin_e)/(n-1);
    photon_more_rate_e += photon_more_per_MeV*(eMax_e - eMin_e)/(n-1);
  }
  //
  Double_t eMin_p = 0.3;
  Double_t eMax_p = 400.0;
  for(i = 0;i<n;i++){
    e = (eMax_p - eMin_p)/(n-1)*i + eMin_e;
    gy_per_MeV = gr_diff_p->Eval(e)*gr_dedx_p->Eval(e)*sTot*tTot*Mev_to_J/tot_mass;
    photon_per_MeV = gr_diff_p->Eval(e)*gr_nPhot_p->Eval(e)*sTot;
    photon_more_per_MeV = gr_diff_p->Eval(e)*gr_nPhot_more_p->Eval(e)*sTot;
    gr_dedx_m_diff_p->SetPoint(gr_dedx_m_diff_p->GetN(),e,gy_per_MeV);
    dose_p += gy_per_MeV*(eMax_p - eMin_p)/(n-1);
    photon_rate_p += photon_per_MeV*(eMax_p - eMin_p)/(n-1);
    photon_more_rate_p += photon_more_per_MeV*(eMax_p - eMin_p)/(n-1);
  }
  //
  //gr_dedx_m_diff_p->Draw("APL");
  gr_dedx_m_diff_e->SetMarkerStyle(6);
  gr_dedx_m_diff_e->SetLineColor(kBlack);
  gr_dedx_m_diff_e->SetLineWidth(3.0);
  gr_dedx_m_diff_e->SetMarkerColor(kBlack);
  gr_dedx_m_diff_e->GetXaxis()->SetTitle("Energy, MeV");
  gr_dedx_m_diff_e->GetYaxis()->SetTitle("Gy/MeV");

  gr_dedx_m_diff_p->SetMarkerStyle(6);
  gr_dedx_m_diff_p->SetLineColor(kBlack);
  gr_dedx_m_diff_p->SetLineWidth(3.0);
  gr_dedx_m_diff_p->SetMarkerColor(kBlack);
  gr_dedx_m_diff_p->GetXaxis()->SetTitle("Energy, MeV");
  gr_dedx_m_diff_p->GetYaxis()->SetTitle("Gy/MeV");

  gr_dedx_m_diff_e->Draw("APL");
  //gr_dedx_m_diff_p->Draw("APL");

  cout<<"dose_e        = "<<dose_e<<endl
      <<"dose_p        = "<<dose_p<<endl
      <<"photon_rate_e = "<<photon_rate_e<<endl
      <<"photon_rate_p = "<<photon_rate_p<<endl
      <<"photon_more_rate_e = "<<photon_more_rate_e<<endl
      <<"photon_more_rate_p = "<<photon_more_rate_p<<endl;

  //  
  TCanvas *c2 = new TCanvas("c2","c2",10,10,800,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  //gStyle->SetOptStat(kFALSE);
  //
  c2->SetRightMargin(0.01);
  c2->SetLeftMargin(0.12);
  c2->SetTopMargin(0.01);
  c2->SetBottomMargin(0.1);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  
  gr_nPhot_more_e->Draw();
  //gr_nPhot_more_p->Draw();
  
  return 0;

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
      gr_tot->SetPoint(gr_tot->GetN(),ekin,ftot);
      gr_diff->SetPoint(gr_diff->GetN(),ekin,fdiff);
    }
    fileIn.close();
  }
  else{
    cout<<"Unable to open file \n";
  }
}

void plot_bkg_p(TGraph *gr_nPhot, TGraph *gr_nPhot_more, TGraph *gr_npe, TGraph *gr_dedx){
  const Int_t n = 27;
  TString fileNarr[n];
  fileNarr[0]="./hist_root_bkg/hist_proton_Ekin_0.200MeV.dat.root";
  fileNarr[1]="./hist_root_bkg/hist_proton_Ekin_0.300MeV.dat.root";
  fileNarr[2]="./hist_root_bkg/hist_proton_Ekin_0.400MeV.dat.root";
  fileNarr[3]="./hist_root_bkg/hist_proton_Ekin_0.500MeV.dat.root";
  fileNarr[4]="./hist_root_bkg/hist_proton_Ekin_0.600MeV.dat.root";
  fileNarr[5]="./hist_root_bkg/hist_proton_Ekin_0.700MeV.dat.root";
  fileNarr[6]="./hist_root_bkg/hist_proton_Ekin_1.000MeV.dat.root";
  fileNarr[7]="./hist_root_bkg/hist_proton_Ekin_1.500MeV.dat.root";
  fileNarr[8]="./hist_root_bkg/hist_proton_Ekin_2.000MeV.dat.root";
  fileNarr[9]="./hist_root_bkg/hist_proton_Ekin_3.000MeV.dat.root";
  fileNarr[10]="./hist_root_bkg/hist_proton_Ekin_4.000MeV.dat.root";
  fileNarr[11]="./hist_root_bkg/hist_proton_Ekin_5.000MeV.dat.root";
  fileNarr[12]="./hist_root_bkg/hist_proton_Ekin_6.000MeV.dat.root";
  fileNarr[13]="./hist_root_bkg/hist_proton_Ekin_7.000MeV.dat.root";
  fileNarr[14]="./hist_root_bkg/hist_proton_Ekin_10.00MeV.dat.root";
  fileNarr[15]="./hist_root_bkg/hist_proton_Ekin_15.00MeV.dat.root";
  fileNarr[16]="./hist_root_bkg/hist_proton_Ekin_20.00MeV.dat.root";
  fileNarr[17]="./hist_root_bkg/hist_proton_Ekin_30.00MeV.dat.root";
  fileNarr[18]="./hist_root_bkg/hist_proton_Ekin_40.00MeV.dat.root";
  fileNarr[19]="./hist_root_bkg/hist_proton_Ekin_50.00MeV.dat.root";
  fileNarr[20]="./hist_root_bkg/hist_proton_Ekin_60.00MeV.dat.root";
  fileNarr[21]="./hist_root_bkg/hist_proton_Ekin_70.00MeV.dat.root";
  fileNarr[22]="./hist_root_bkg/hist_proton_Ekin_100.0MeV.dat.root";
  fileNarr[23]="./hist_root_bkg/hist_proton_Ekin_150.0MeV.dat.root";
  fileNarr[24]="./hist_root_bkg/hist_proton_Ekin_200.0MeV.dat.root";
  fileNarr[25]="./hist_root_bkg/hist_proton_Ekin_300.0MeV.dat.root";
  fileNarr[26]="./hist_root_bkg/hist_proton_Ekin_400.0MeV.dat.root";
  //
  Double_t ekin[n] = { 0.200, 0.300, 0.400, 0.500, 0.600,
		       0.700, 1.000, 1.500, 2.000, 3.000,
		       4.000, 5.000, 6.000, 7.000, 10.00,
		       15.00, 20.00, 30.00, 40.00, 50.00,
		       60.00, 70.00, 100.0, 150.0, 200.0,
		       300.0, 400.0};
  //
  //
  TFile *f1[n];
  TH1D *h1_nPhot[n];
  TH1D *h1_npe[n];
  TH1D *h1_dedx_sipm[n];
  Double_t nPhot_mean[n];
  Double_t npe_mean[n];
  Double_t dedx_mean[n];
  Double_t nPhot_more;
  for(Int_t i = 0;i<n;i++){
   f1[i] = new TFile(fileNarr[i].Data());
   h1_nPhot[i] = (TH1D*)f1[i]->Get("h1_nPhot");
   h1_npe[i] = (TH1D*)f1[i]->Get("h1_npe");
   h1_dedx_sipm[i] = (TH1D*)f1[i]->Get("h1_dedx_sipm");
   //
   nPhot_mean[i] = h1_nPhot[i]->GetMean();
   npe_mean[i] = h1_npe[i]->GetMean()*h1_npe[i]->GetEntries()/h1_nPhot[i]->GetEntries();
   dedx_mean[i] = h1_dedx_sipm[i]->GetMean();
   //
   gr_nPhot->SetPoint(gr_nPhot->GetN(),ekin[i],nPhot_mean[i]);
   nPhot_more = getProbOfmoreThenThreshold( h1_nPhot[i], 20.0);
   gr_nPhot_more->SetPoint(gr_nPhot_more->GetN(),ekin[i],nPhot_more);
   gr_npe->SetPoint(gr_npe->GetN(),ekin[i],npe_mean[i]);
   gr_dedx->SetPoint(gr_dedx->GetN(),ekin[i],dedx_mean[i]);
  }

  //
  //h1[0]->Draw();
  //h1[0]->SetMaximum(0.5);
  //h1[0]->SetMinimum(1.0e-8);
  //
  /*
  for(Int_t i = 0;i<n;i++){
    gr1[i]->SetLineWidth(3.0);
    gr1[i]->SetLineColor(kBlack);
    gr1[i]->SetMarkerStyle(20);
    gr1[i]->SetMarkerColor(kBlack);
  }

  //
  gr1[0]->SetMarkerStyle(6);
  gr1[0]->SetLineColor(kBlack);
  gr1[0]->SetMarkerColor(kBlack);
  //
  gr1[1]->SetMarkerStyle(6);
  gr1[1]->SetLineColor(kRed);
  gr1[1]->SetMarkerColor(kRed);
  //
  gr1[2]->SetMarkerStyle(6);
  gr1[2]->SetLineColor(kBlue);
  gr1[2]->SetMarkerColor(kBlue);
  //
  gr1[3]->SetMarkerStyle(6);
  gr1[3]->SetLineColor(kGreen+2);
  gr1[3]->SetMarkerColor(kGreen+2);
  //
  gr1[4]->SetMarkerStyle(20);
  gr1[4]->SetLineColor(kBlack);
  gr1[4]->SetMarkerColor(kBlack);
  //
  gr1[5]->SetMarkerStyle(20);
  gr1[5]->SetLineColor(kRed);
  gr1[5]->SetMarkerColor(kRed);
  //
  gr1[6]->SetMarkerStyle(20);
  gr1[6]->SetLineColor(kBlue);
  gr1[6]->SetMarkerColor(kBlue);
  //
  gr1[7]->SetMarkerStyle(20);
  gr1[7]->SetLineColor(kGreen+2);
  gr1[7]->SetMarkerColor(kGreen+2);

  TMultiGraph *mg = new TMultiGraph();
  for(Int_t i = 0;i<n;i++)
    mg->Add(gr1[i]);

  //mg->GetXaxis()->SetLimits(xMin,xMax);
  mg->Draw("APL");
  mg->SetMinimum(1.0e-8);
  mg->SetMaximum(1.0);
  //mg->SetMaximum(6.0);
  //mg->SetMaximum(5.0);
  mg->GetXaxis()->SetTitle("# p.e.");
  mg->GetYaxis()->SetTitle("efficiency");

  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr1[0], "10   PeV", "pl");
  leg->AddEntry(gr1[1], "35   PeV", "pl");
  leg->AddEntry(gr1[2], "70   PeV", "pl");
  leg->AddEntry(gr1[3], "100  PeV", "pl");
  leg->AddEntry(gr1[4], "350  PeV", "pl");
  leg->AddEntry(gr1[5], "700  PeV", "pl");
  leg->AddEntry(gr1[6], "1000 PeV", "pl");
  leg->AddEntry(gr1[7], "3000 PeV", "pl");
  leg->Draw();

  //
  TCanvas *c2 = new TCanvas("c2","c2",10,10,1200,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  c2->SetRightMargin(0.01);
  c2->SetLeftMargin(0.08);
  c2->SetTopMargin(0.01);
  c2->SetBottomMargin(0.08);
  //
  //
  TGraphErrors *gr_eff_vs_e = new TGraphErrors();
  gr_eff_vs_e->SetName("gr_eff_vs_e");
  gr_eff_vs_e->SetTitle("gr_eff_vs_e");
  int pointID;
  //
  gr_eff_vs_e->SetMarkerStyle(20);
  gr_eff_vs_e->SetLineColor(kBlack);
  gr_eff_vs_e->SetLineWidth(3.0);
  gr_eff_vs_e->SetMarkerColor(kBlack);
  //
  pointID = gr_eff_vs_e->GetN();
  gr_eff_vs_e->SetPoint(pointID,10,1.0e-9);
  //gr_eff_vs_e->SetPointError(pointID,0.0,1.0e-9/10.0);
  gr_eff_vs_e->SetPointError(pointID,0.0,0.0);
  //
  pointID = gr_eff_vs_e->GetN();
  gr_eff_vs_e->SetPoint(pointID,35,h1[1]->GetBinContent(8));
  //gr_eff_vs_e->SetPointError(pointID,0.0,h1[1]->GetBinContent(8)/3.0);
  gr_eff_vs_e->SetPointError(pointID,0.0,0.0);
  //
  pointID = gr_eff_vs_e->GetN();
  gr_eff_vs_e->SetPoint(pointID,70,h1[2]->GetBinContent(8));
  //gr_eff_vs_e->SetPointError(pointID,0.0,h1[2]->GetBinContent(8)/3.0);
  gr_eff_vs_e->SetPointError(pointID,0.0,0.0);
  //
  pointID = gr_eff_vs_e->GetN();
  gr_eff_vs_e->SetPoint(pointID,100,h1[3]->GetBinContent(8));
  //gr_eff_vs_e->SetPointError(pointID,0.0,h1[3]->GetBinContent(8)/3.0);
  gr_eff_vs_e->SetPointError(pointID,0.0,0.0);
  //
  pointID = gr_eff_vs_e->GetN();
  gr_eff_vs_e->SetPoint(pointID,350,h1[4]->GetBinContent(8));
  //gr_eff_vs_e->SetPointError(pointID,0.0,h1[4]->GetBinContent(8)/100.0);
  gr_eff_vs_e->SetPointError(pointID,0.0,0.0);
  //
  pointID = gr_eff_vs_e->GetN();
  gr_eff_vs_e->SetPoint(pointID,700,h1[5]->GetBinContent(8));
  //gr_eff_vs_e->SetPointError(pointID,0.0,h1[5]->GetBinContent(8)/100.0);
  gr_eff_vs_e->SetPointError(pointID,0.0,0.0);
  //
  pointID = gr_eff_vs_e->GetN();
  gr_eff_vs_e->SetPoint(pointID,1000,h1[6]->GetBinContent(8));
  //gr_eff_vs_e->SetPointError(pointID,0.0,h1[6]->GetBinContent(8)/100.0);
  gr_eff_vs_e->SetPointError(pointID,0.0,0.0);
  //
  pointID = gr_eff_vs_e->GetN();
  gr_eff_vs_e->SetPoint(pointID,3000,h1[7]->GetBinContent(8));
  //gr_eff_vs_e->SetPointError(pointID,0.0,h1[7]->GetBinContent(8)/100.0);
  gr_eff_vs_e->SetPointError(pointID,0.0,0.0);
  //
  gr_eff_vs_e->SetTitle("");
  gr_eff_vs_e->Draw("APL");

  gr_eff_vs_e->GetXaxis()->SetTitle("Energy, PeV");
  gr_eff_vs_e->GetYaxis()->SetTitle("efficiency");

  //gPad->SetGridx();
  //gPad->SetGridy();
  */
  //////////

  /*
  Double_t e_min01 = 100;
  Double_t e_max01 = 355;
  Int_t npar_eff_Fit_func = 2;
  TF1 *f_eff_vs_e01 = new TF1("f_eff_vs_e01", eff_Fit_func, e_min01, e_max01, npar_eff_Fit_func);
  f_eff_vs_e01->SetParameter(0,0.0);
  f_eff_vs_e01->SetParameter(1,1.0);  
  
  Double_t e_min02 = 65;
  Double_t e_max02 = 105;
  npar_eff_Fit_func = 2;
  TF1 *f_eff_vs_e02 = new TF1("f_eff_vs_e02", eff_Fit_func, e_min01, e_max01, npar_eff_Fit_func);
  f_eff_vs_e02->SetParameter(0,0.0);
  f_eff_vs_e02->SetParameter(1,1.0);  

  gr_eff_vs_e->Fit("f_eff_vs_e01","","",e_min01,e_max01);
  gr_eff_vs_e->Fit("f_eff_vs_e02","","",e_min02,e_max02);

  f_eff_vs_e01->Draw("same");
  f_eff_vs_e02->Draw("same");
  */
  /*  
  ///////////////////
  Double_t e_min = 10;
  Double_t e_max = 3000;
  Double_t loge_min = TMath::Log10(e_min);
  Double_t loge_max = TMath::Log10(e_max);
  Double_t loge;
  Double_t el;
  ///////////////////


  ///////////////////
  Int_t nn = 1000;
  TGraph *gr_eff_vs_e_extrapolated = new TGraph();
  gr_eff_vs_e_extrapolated->SetName("gr_eff_vs_e_extrapolated");
  gr_eff_vs_e_extrapolated->SetTitle("gr_eff_vs_e_extrapolated");
  TGraph *gr_eff_vs_e_extrapolated_log = new TGraph();
  gr_eff_vs_e_extrapolated_log->SetName("gr_eff_vs_e_extrapolated_log");
  gr_eff_vs_e_extrapolated_log->SetTitle("gr_eff_vs_e_extrapolated_log");
  Double_t logLog_spline;
  for(Int_t i = 0;i<nn;i++){
    loge = loge_min + (loge_max-loge_min)/(nn-1)*i;
    el = TMath::Power(10,loge);
    gr_eff_vs_e_extrapolated->SetPoint(gr_eff_vs_e_extrapolated->GetN(),el,gr_eff_vs_e->Eval(el));
  }
  //
  e_min = 10;
  e_max = 2999;
  nn = 30;
  loge_min = TMath::Log10(e_min);
  loge_max = TMath::Log10(e_max);
  for(Int_t i = 0;i<nn;i++){
    loge = loge_min + (loge_max-loge_min)/(nn-1)*i;
    el = TMath::Power(10,loge);
    logLog_spline = getLogLog_spline(gr_eff_vs_e,el);
    gr_eff_vs_e_extrapolated_log->SetPoint(gr_eff_vs_e_extrapolated_log->GetN(),el,logLog_spline);
  }

  //
  gr_eff_vs_e_extrapolated->SetMarkerStyle(6);
  gr_eff_vs_e_extrapolated->SetLineColor(kRed);
  gr_eff_vs_e_extrapolated->SetLineWidth(3.0);
  gr_eff_vs_e_extrapolated->SetMarkerColor(kRed);
  //
  gr_eff_vs_e_extrapolated_log->SetMarkerStyle(20);
  gr_eff_vs_e_extrapolated_log->SetLineColor(kBlack);
  gr_eff_vs_e_extrapolated_log->SetLineWidth(3.0);
  gr_eff_vs_e_extrapolated_log->SetMarkerColor(kBlue);
  //
  
  TCanvas *c3 = new TCanvas("c3","c3",10,10,1200,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  c3->SetRightMargin(0.01);
  c3->SetLeftMargin(0.08);
  c3->SetTopMargin(0.01);
  c3->SetBottomMargin(0.08);
  //
  
  TMultiGraph *mg2 = new TMultiGraph();
  mg2->Add(gr_eff_vs_e);
  mg2->Add(gr_eff_vs_e_extrapolated);
  mg2->Add(gr_eff_vs_e_extrapolated_log);
  mg2->Draw("APL");
  mg2->SetMinimum(1.0e-9);
  mg2->SetMaximum(1.0);
  
  TLegend *leg03 = new TLegend(0.6,0.3,0.9,0.6,"","brNDC");
  leg03->AddEntry(gr_eff_vs_e, "MC.", "pl");
  leg03->AddEntry(gr_eff_vs_e_extrapolated, "Linear interpolation", "pl");
  leg03->AddEntry(gr_eff_vs_e_extrapolated_log, "Log log interpolation", "pl");
  leg03->Draw();

  mg2->GetXaxis()->SetTitle("Energy, PeV");
  mg2->GetYaxis()->SetTitle("Efficiency");
  //----------------------------------
  //
  //
  Double_t eff;
  Double_t e_In_GeV;
  Double_t e_In_PeV;
  for(Int_t i = 1;i<=h1_particle_per_year->GetNbinsX();i++){
    e_In_GeV = h1_particle_per_year->GetBinCenter(i);
    e_In_PeV = e_In_GeV/1.0e+6;
    //
    //cout<<"e_In_GeV "<<e_In_GeV<<endl
    //	<<"e_In_PeV "<<e_In_PeV<<endl;
    //
    if(e_In_PeV<3000){
      eff = getLogLog_spline(gr_eff_vs_e,e_In_PeV);
      h1_particle_per_year_eff->SetBinContent(i,h1_particle_per_year->GetBinContent(i)*eff);
      eff = gr_eff_vs_e->Eval(e_In_PeV);
      h1_particle_per_year_eff_lin->SetBinContent(i,h1_particle_per_year->GetBinContent(i)*eff);
    }
    else{
      eff = getLogLog_spline(gr_eff_vs_e,2999);
      h1_particle_per_year_eff->SetBinContent(i,h1_particle_per_year->GetBinContent(i)*eff);
      eff = gr_eff_vs_e->Eval(2999);
      h1_particle_per_year_eff_lin->SetBinContent(i,h1_particle_per_year->GetBinContent(i)*eff);
    }
  }
  //
  //
  //----------------------------------
  //
  TCanvas *c4 = new TCanvas("c4","c4",10,10,1200,800);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  //gStyle->SetOptStat(kFALSE);
  //
  c4->SetRightMargin(0.01);
  c4->SetLeftMargin(0.08);
  c4->SetTopMargin(0.01);
  c4->SetBottomMargin(0.08);
  //h1_particle_per_year_eff->Draw();
  ////////////////////////////////////

  h1_particle_per_year_tot->SetTitle("");

  h1_particle_per_year_tot->SetName("Total");
  h1_particle_per_year_eff->SetName("Log-log ");
  h1_particle_per_year_eff_lin->SetName("Linear");

  h1_particle_per_year_tot->SetLineColor(kBlack);
  h1_particle_per_year_eff->SetLineColor(kRed);
  h1_particle_per_year_eff_lin->SetLineColor(kBlue);

  h1_particle_per_year_tot->SetLineWidth(3.0);
  h1_particle_per_year_eff->SetLineWidth(3.0);
  h1_particle_per_year_eff_lin->SetLineWidth(3.0);

  h1_particle_per_year_tot->Draw();
  h1_particle_per_year_eff->Draw("sames");
  h1_particle_per_year_eff_lin->Draw("sames");

  h1_particle_per_year_tot->GetXaxis()->SetTitle("Energy, GeV");
  h1_particle_per_year_tot->GetYaxis()->SetTitle("Number of events during 3 years");
  
  TLegend *leg04 = new TLegend(0.6,0.3,0.9,0.6,"","brNDC");
  leg04->AddEntry(h1_particle_per_year_tot, "Total", "l");
  leg04->AddEntry(h1_particle_per_year_eff, "Log log interpolation","l");
  leg04->AddEntry(h1_particle_per_year_eff_lin, "Linear interpolation", "l");
  leg04->Draw();
  */  

}

void plot_bkg_e(TGraph *gr_nPhot, TGraph *gr_nPhot_more, TGraph *gr_npe, TGraph *gr_dedx){
  const Int_t n = 28;
  TString fileNarr[n];
  //
  fileNarr[0]="./hist_root_bkg/hist_e-_Ekin_0.20MeV.dat.root";
  fileNarr[1]="./hist_root_bkg/hist_e-_Ekin_0.30MeV.dat.root";
  fileNarr[2]="./hist_root_bkg/hist_e-_Ekin_0.40MeV.dat.root";
  fileNarr[3]="./hist_root_bkg/hist_e-_Ekin_0.50MeV.dat.root";
  fileNarr[4]="./hist_root_bkg/hist_e-_Ekin_0.60MeV.dat.root";
  //
  fileNarr[5]="./hist_root_bkg/hist_e-_Ekin_0.70MeV.dat.root";
  fileNarr[6]="./hist_root_bkg/hist_e-_Ekin_0.80MeV.dat.root";
  fileNarr[7]="./hist_root_bkg/hist_e-_Ekin_1.00MeV.dat.root";
  fileNarr[8]="./hist_root_bkg/hist_e-_Ekin_1.25MeV.dat.root";
  fileNarr[9]="./hist_root_bkg/hist_e-_Ekin_1.50MeV.dat.root";
  //
  fileNarr[10]="./hist_root_bkg/hist_e-_Ekin_1.75MeV.dat.root";
  fileNarr[11]="./hist_root_bkg/hist_e-_Ekin_2.00MeV.dat.root";
  fileNarr[12]="./hist_root_bkg/hist_e-_Ekin_2.25MeV.dat.root";
  fileNarr[13]="./hist_root_bkg/hist_e-_Ekin_2.50MeV.dat.root";
  fileNarr[14]="./hist_root_bkg/hist_e-_Ekin_2.75MeV.dat.root";
  //
  fileNarr[15]="./hist_root_bkg/hist_e-_Ekin_3.00MeV.dat.root";
  fileNarr[16]="./hist_root_bkg/hist_e-_Ekin_3.25MeV.dat.root";
  fileNarr[17]="./hist_root_bkg/hist_e-_Ekin_3.50MeV.dat.root";
  fileNarr[18]="./hist_root_bkg/hist_e-_Ekin_3.75MeV.dat.root";
  fileNarr[19]="./hist_root_bkg/hist_e-_Ekin_4.00MeV.dat.root";
  //
  fileNarr[20]="./hist_root_bkg/hist_e-_Ekin_4.25MeV.dat.root";
  fileNarr[21]="./hist_root_bkg/hist_e-_Ekin_4.50MeV.dat.root";
  fileNarr[22]="./hist_root_bkg/hist_e-_Ekin_4.75MeV.dat.root";
  fileNarr[23]="./hist_root_bkg/hist_e-_Ekin_5.00MeV.dat.root";
  fileNarr[24]="./hist_root_bkg/hist_e-_Ekin_5.50MeV.dat.root";
  //
  fileNarr[25]="./hist_root_bkg/hist_e-_Ekin_6.00MeV.dat.root";
  fileNarr[26]="./hist_root_bkg/hist_e-_Ekin_6.50MeV.dat.root";
  fileNarr[27]="./hist_root_bkg/hist_e-_Ekin_7.00MeV.dat.root";
  //
  Double_t ekin[n] = { 0.20, 0.30, 0.40, 0.50, 0.60,
		       0.70, 0.80, 1.00, 1.25, 1.50,
		       1.75, 2.00, 2.25, 2.50, 2.75,
		       3.00, 3.25, 3.50, 3.75, 4.00,
		       4.25, 4.50, 4.75, 5.00, 5.50,
		       6.00, 6.50, 7.00};
  //
  //
  TFile *f1[n];
  TH1D *h1_nPhot[n];
  TH1D *h1_npe[n];
  TH1D *h1_dedx_sipm[n];
  Double_t nPhot_mean[n];
  Double_t npe_mean[n];
  Double_t dedx_mean[n];
  Double_t nPhot_more;
  for(Int_t i = 0;i<n;i++){
   f1[i] = new TFile(fileNarr[i].Data());
   h1_nPhot[i] = (TH1D*)f1[i]->Get("h1_nPhot");
   h1_npe[i] = (TH1D*)f1[i]->Get("h1_npe");
   h1_dedx_sipm[i] = (TH1D*)f1[i]->Get("h1_dedx_sipm");
   //
   nPhot_mean[i] = h1_nPhot[i]->GetMean();
   npe_mean[i] = h1_npe[i]->GetMean()*h1_npe[i]->GetEntries()/h1_nPhot[i]->GetEntries();
   dedx_mean[i] = h1_dedx_sipm[i]->GetMean();
   //
   gr_nPhot->SetPoint(gr_nPhot->GetN(),ekin[i],nPhot_mean[i]);
   nPhot_more = getProbOfmoreThenThreshold( h1_nPhot[i], 20.0);
   gr_nPhot_more->SetPoint(gr_nPhot_more->GetN(),ekin[i],nPhot_more);
   gr_npe->SetPoint(gr_npe->GetN(),ekin[i],npe_mean[i]);
   gr_dedx->SetPoint(gr_dedx->GetN(),ekin[i],dedx_mean[i]);
  }
}

Double_t getProbOfmoreThenThreshold( TH1D *h1, Double_t thresh){
  Double_t nTot = h1->Integral(1,h1->GetNbinsX());
  Double_t nMore = h1->Integral(h1->FindBin(thresh),h1->GetNbinsX());
  if(nTot>0.0)
    return nMore/nTot;
  return 0.0;
}
