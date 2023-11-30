//https://arxiv.org/pdf/1406.6892.pdf
//Measurement of the Crab Nebula spectrum over three decades in energy with the MAGIC telescopes

//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

void get_log_parabola_fit_MAGIC_tel_crab(TGraph *gr_tev, Int_t nn, Double_t e_min, Double_t e_max);
void load_log_parabola_fit_MAGIC_tel_crab(TGraphErrors *gr_tev, TString name_dat_file);
Double_t function_log_parabola_fit_MAGIC_tel_crab(Double_t e);
Double_t integral_log_parabola(Int_t nn, Double_t e_min, Double_t e_max);
Double_t test_f(Double_t e);
Double_t test_int(Double_t e_min, Double_t e_max);
//
Double_t function_log_parabola_fit_MAGIC_tel_crab_GeV(Double_t e);
Double_t integral_log_parabola_GeV(Int_t nn, Double_t e_min, Double_t e_max);
//
Double_t log_parabola(Double_t *x, Double_t *par){
  Double_t A  = par[0];
  Double_t B  = par[1];
  Double_t C  = par[2];
  Double_t D  = par[3];
  Double_t val = A*(TMath::Power(x[0]/B,-C + D*TMath::Log10(x[0]/B)));
  return val;
}


Int_t log_parabola_fit_MAGIC_tel_crab(){
  //
  Int_t nn = 1000;
  Double_t e_min = 0.05; //TeV 50 GeV
  Double_t e_max = 30.0; //TeV 30000 GeV
  //
  TGraphErrors *gr_tev = new TGraphErrors();
  gr_tev->SetNameTitle("gr_tev","gr_tev");
  TGraph *gr_tev_f = new TGraph();
  gr_tev_f->SetNameTitle("gr_tev_f","gr_tev_f");
  //
  get_log_parabola_fit_MAGIC_tel_crab(gr_tev_f, nn, e_min, e_max);
  load_log_parabola_fit_MAGIC_tel_crab(gr_tev, "./data/data_crab_digi_TeV.csv");
  //


  //Fit
  const Int_t npar = 4;
  Double_t inParameters[npar];
  inParameters[0] = 3.23*1.0e-11;
  inParameters[1] = 1.0;
  inParameters[2] = 2.47;
  inParameters[3] = -0.24;
  TF1 *f_log_parabola = new TF1( "f_log_parabola", log_parabola, e_min, e_max, npar);
  f_log_parabola->SetParameters(inParameters);
  f_log_parabola->FixParameter(0,inParameters[0]);
  f_log_parabola->FixParameter(1,inParameters[1]);
  f_log_parabola->FixParameter(2,inParameters[2]);
  f_log_parabola->FixParameter(3,inParameters[3]);
  gr_tev->Fit("f_log_parabola","","",e_min, e_max);
  //
  TCanvas *c1 = new TCanvas("c1","https://arxiv.org/pdf/1406.6892.pdf",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  gPad->SetLogy();
  //
  //gr_proton->Draw("APL");
  //gr_DAMPE->Draw("APL");
  //
  TMultiGraph *mg = new TMultiGraph();
  //
  gr_tev->SetMarkerStyle(20);
  gr_tev->SetMarkerColor(kBlack);
  gr_tev->SetLineWidth(2);
  gr_tev->SetLineColor(kBlack);
  //
  mg->Add(gr_tev);
  mg->Add(gr_tev_f);
  //
  mg->Draw("APL");
  //
  /*
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr_DAMPE, "DAMPE Collaboration,Sci. Adv.2019;5:eaax3793 (Proton)", "apl");
  leg->AddEntry(gr_DAMPE_He, "DAMPE Collaboration (Helium)", "apl");
  leg->AddEntry(gr_flux_p_PDG, "PDG data", "apl");
  //
  leg->AddEntry(gr_TAUP2023_p, "TAUP2023 p", "apl");
  leg->AddEntry(gr_TAUP2023_he, "TAUP2023 he", "apl");
  leg->AddEntry(gr_TAUP2023_O, "TAUP2023 O", "apl");
  leg->AddEntry(gr_TAUP2023_Fe, "TAUP2023 Fe", "apl");
  //
  leg->AddEntry(gr_TAUP2023_gamma_bkg_galactic, "Gamma bkg. (diffuse) galactic", "apl");
  leg->AddEntry(gr_TAUP2023_gamma_bkg_extragalactic, "Gamma bkg. (diffuse) extragalactic", "apl");
  //
  leg->AddEntry(gr_TAUP2023_epem, "e+e-", "apl");
  leg->AddEntry(gr_TAUP2023_ep, "e+", "apl");
  //
  leg->Draw();
  //
  */
  mg->GetXaxis()->SetTitle("E, TeV");
  mg->GetYaxis()->SetTitle("dN/dE, cm^{-2} s^{-1} TeV^{-1}");

  //
  cout<<integral_log_parabola(10000,  e_min,  e_max)<<endl;
  cout<<integral_log_parabola(100000,  e_min,  e_max)<<endl;
  cout<<integral_log_parabola(1000000,  e_min,  e_max)<<endl;
  cout<<integral_log_parabola(100000000,  e_min,  e_max)<<endl;
  //
  cout<<test_int(e_min, e_max)<<endl;
  //
  cout<<integral_log_parabola_GeV(10000000,  e_min*1000.0,  e_max*1000.0)<<endl;
  
  return 0;
}

void get_log_parabola_fit_MAGIC_tel_crab(TGraph *gr_tev, Int_t nn, Double_t e_min, Double_t e_max){
  Double_t f;
  Double_t e;
  for(Int_t i = 0;i<nn;i++){
    e = (e_max - e_min)/(nn-1)*i + e_min;
    f=function_log_parabola_fit_MAGIC_tel_crab(e);
    gr_tev->SetPoint(i,e,f);
  }  
}

Double_t function_log_parabola_fit_MAGIC_tel_crab(Double_t e){
  return ((3.31327)*1.0e-11)*TMath::Power((e/1.00729),(-2.49570 - 0.130707*TMath::Log(e/1.00729)));
}

Double_t function_log_parabola_fit_MAGIC_tel_crab_GeV(Double_t e){
  return ((3.31327)*1.0e-14)*TMath::Power((e/1.00729/1000.0),(-2.49570 - 0.130707*TMath::Log(e/1.00729/1000.0)));
}

Double_t integral_log_parabola_GeV(Int_t nn, Double_t e_min, Double_t e_max){
  Double_t de = (e_max - e_min)/(nn-1);
  Double_t f;
  Double_t e;
  Double_t integral = 0.0;
  for(Int_t i = 0;i<nn;i++){
    e = de*i + e_min;
    f=function_log_parabola_fit_MAGIC_tel_crab_GeV(e);
    //f=test_f(e);
    integral += f*de;
  }
  return integral;
}

Double_t integral_log_parabola(Int_t nn, Double_t e_min, Double_t e_max){
  Double_t de = (e_max - e_min)/(nn-1);
  Double_t f;
  Double_t e;
  Double_t integral = 0.0;
  for(Int_t i = 0;i<nn;i++){
    e = de*i + e_min;
    f=function_log_parabola_fit_MAGIC_tel_crab(e);
    //f=test_f(e);
    integral += f*de;
  }
  return integral;
}

Double_t test_f(Double_t e){
  return e*e;
}

Double_t test_int(Double_t e_min, Double_t e_max){
  return 1.0/3.0*(e_max*e_max*e_max - e_min*e_min*e_min);
}
  
void load_log_parabola_fit_MAGIC_tel_crab(TGraphErrors *gr, TString name_dat_file){
  Double_t Ekin;
  Double_t F;
  string mot;
  //  
  ifstream indata;
  indata.open(name_dat_file.Data());
  assert(indata.is_open());  
  Int_t np;
  while(indata>>Ekin>>F){
    cout<<Ekin<<endl;
    np=gr->GetN();
    gr->SetPoint(np,Ekin,F);
    gr->SetPointError(np,0.01,F/10.0);
  }
  indata.close();
}
