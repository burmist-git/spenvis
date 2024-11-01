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

void get_MAGIC_crab_data(TString name_dat_file, TGraphErrors *gr, Int_t &nn, Double_t &e_min, Double_t &e_max);
void get_MAGIC_Fermi_LAT_crab_data(TString name_dat_file, TGraphErrors *gr, Int_t &nn, Double_t &e_min, Double_t &e_max);
Double_t log_parabola_Crab_Nebula_MAGIC (Double_t eGeV);
void fill_gr_flux_crab_from_fit(TGraph *gr, Int_t nn, Double_t eGeV_min, Double_t eGeV_max);
//Double_t modiﬁed_log_parabola_function_Crab_Nebula_MAGIC_Fermi_LAT(Double_t eGeV);
void fill_gr_modiﬁed_log_parabola_function_Crab_Nebula_MAGIC_Fermi_LAT_fit(TGraph *gr, Int_t nn, Double_t eGeV_min, Double_t eGeV_max);
Double_t modiﬁed_log_parabola(Double_t *x, Double_t *par);
Double_t modiﬁed_log_parabola(Double_t eGeV);
void fill_gr_modiﬁed_log_parabola_fit(TGraph *gr, Int_t nn, Double_t eGeV_min, Double_t eGeV_max);
Double_t modiﬁed_log_parabola_TeV(Double_t eTeV);
Double_t modiﬁed_log_parabola_TeV_integral(Double_t eTeV_min, Double_t eTeV_max);
Double_t modiﬁed_log_parabola_GeV_integral(Double_t eGeV_min, Double_t eGeV_max);

Int_t flux_crab_old(){
  //
  TGraphErrors *gr_flux_crab = new TGraphErrors();
  TGraphErrors *gr_flux_crab_MAGIC_Fermi_LAT = new TGraphErrors();
  TGraph *gr_flux_crab_from_fit = new TGraph();
  TGraph *gr_flux_crab_MAGIC_Fermi_LAT_from_fit = new TGraph();
  TGraph *gr_flux_crab_MAGIC_Fermi_LAT_fit = new TGraph();
  Int_t nn;
  Double_t e_min;
  Double_t e_max;
  //
  Int_t nn_all;
  Double_t e_min_all;
  Double_t e_max_all;
  //
  //
  //
  get_MAGIC_crab_data("data/MAGICdata.txt", gr_flux_crab, nn, e_min, e_max);  
  get_MAGIC_Fermi_LAT_crab_data("data/Crab_Nebul_MAGIC_stereo_Fermi_LAT_data.dat", gr_flux_crab_MAGIC_Fermi_LAT, nn_all, e_min_all, e_max_all);
  fill_gr_flux_crab_from_fit(gr_flux_crab_from_fit, 300,  1.0, 30000.0);
  fill_gr_modiﬁed_log_parabola_fit(gr_flux_crab_MAGIC_Fermi_LAT_fit, 300,  1.0, 30000.0);
  //
  //fill_gr_modiﬁed_log_parabola_function_Crab_Nebula_MAGIC_Fermi_LAT_fit(gr_flux_crab_MAGIC_Fermi_LAT_from_fit, 1000,  1.0, 30000.0);
  //
  cout<<"nn        "<<nn<<endl
      <<"e_min     "<<e_min<<endl
      <<"e_max     "<<e_max<<endl;
  //
  cout<<"nn_all    "<<nn_all<<endl
      <<"e_min_all "<<e_min_all<<endl
      <<"e_max_all "<<e_max_all<<endl;
  //

    //Fit
  const Int_t npar = 6;
  Double_t inParameters[npar];
  //inParameters[0] = 3.23*1.0e-11;
  //inParameters[1] = 1000.0;
  //inParameters[2] = 2.47;
  //inParameters[3] = 0.24;
  //inParameters[4] = 48.0/48.0;
  //inParameters[5] = 1.0;  
  //
  inParameters[0] = 2.76203e-11;
  inParameters[1] = 1.00000e+03;
  inParameters[2] = 2.72388e+00;
  inParameters[3] = 1.73776e-01;
  inParameters[4] = 2.62158e+01;
  inParameters[5] = 1.00000e+00;
  //
  TF1 *f_modiﬁed_log_parabola = new TF1("f_modiﬁed_log_parabola", modiﬁed_log_parabola, e_min, e_max, npar);
  f_modiﬁed_log_parabola->SetParameters(inParameters);
  //f_modiﬁed_log_parabola->FixParameter(0,inParameters[0]);
  f_modiﬁed_log_parabola->FixParameter(1,inParameters[1]);
  //f_modiﬁed_log_parabola->FixParameter(2,inParameters[2]);
  //f_modiﬁed_log_parabola->FixParameter(3,inParameters[3]);  
  //f_modiﬁed_log_parabola->FixParameter(4,inParameters[4]);
  //f_modiﬁed_log_parabola->FixParameter(5,inParameters[5]);  
  //gr_flux_crab->Fit("f_modiﬁed_log_parabola","","",e_min, e_max);
  gr_flux_crab_MAGIC_Fermi_LAT->Fit("f_modiﬁed_log_parabola","","",e_min_all, e_max_all);
  //

  cout<<modiﬁed_log_parabola_TeV_integral(0.001, 20.0)*TMath::Pi()*80000.0*80000.0<<endl;
  cout<<modiﬁed_log_parabola_GeV_integral(1, 20000.0)*TMath::Pi()*80000.0*80000.0<<endl;
  
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
  //
  gr_flux_crab->SetMarkerStyle(20);
  gr_flux_crab->SetMarkerColor(kBlack);
  gr_flux_crab->SetLineWidth(2);
  gr_flux_crab->SetLineColor(kBlack);
  //
  gr_flux_crab_from_fit->SetMarkerStyle(7);
  gr_flux_crab_from_fit->SetMarkerColor(kRed);
  gr_flux_crab_from_fit->SetLineWidth(2);
  gr_flux_crab_from_fit->SetLineColor(kRed);
  //  
  //
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr_flux_crab);
  mg->Add(gr_flux_crab_MAGIC_Fermi_LAT);
  mg->Add(gr_flux_crab_MAGIC_Fermi_LAT_fit);  
  //mg->Add(gr_flux_crab_from_fit);
  //mg->Add(gr_flux_crab_MAGIC_Fermi_LAT_from_fit);
  mg->Draw("AP");
  //
  mg->GetXaxis()->SetTitle("E, GeV");
  mg->GetYaxis()->SetTitle("dN/dE/dA/dt, TeV^{-1} cm^{-2} s^{-1}");
  //
  //

  /*
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
*/
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
  /*

  //
  cout<<integral_log_parabola(10000,  e_min,  e_max)<<endl;
  cout<<integral_log_parabola(100000,  e_min,  e_max)<<endl;
  cout<<integral_log_parabola(1000000,  e_min,  e_max)<<endl;
  cout<<integral_log_parabola(100000000,  e_min,  e_max)<<endl;
  //
  cout<<test_int(e_min, e_max)<<endl;
  //
  cout<<integral_log_parabola_GeV(10000000,  e_min*1000.0,  e_max*1000.0)<<endl;
  */
  return 0;
}

void get_MAGIC_crab_data(TString name_dat_file, TGraphErrors *gr, Int_t &nn, Double_t &e_min, Double_t &e_max){
  //
  Double_t Ekin;
  Double_t EkinErr;
  Double_t F;
  string mot;
  //
  nn = 0;
  e_min = 100;
  e_max = 100;
  //  
  ifstream indata;
  indata.open(name_dat_file.Data());
  assert(indata.is_open());  
  Int_t np = 0;
  indata>>mot>>mot>>mot>>mot;
  //cout<<mot<<endl;
  while(indata>>F>>Ekin>>EkinErr){
    //cout<<Ekin<<endl;
    //F=F*(Ekin/1000.0)*(Ekin/1000.0);
    F=F/1000.0;
    np=gr->GetN();
    if(np==0)
      e_min = Ekin;
    gr->SetPoint(np,Ekin,F);
    gr->SetPointError(np,Ekin*EkinErr/100.0,F*0.01 + 1.0*1.0e-15);
  }
  nn = gr->GetN();
  e_max = Ekin;
  indata.close();
}

void get_MAGIC_Fermi_LAT_crab_data(TString name_dat_file, TGraphErrors *gr, Int_t &nn, Double_t &e_min, Double_t &e_max){
  //
  Double_t Ekin;
  Double_t F;
  string mot;
  //
  nn = 0;
  e_min = 100;
  e_max = 100;
  //  
  ifstream indata;
  indata.open(name_dat_file.Data());
  assert(indata.is_open());  
  Int_t np = 0;
  while(indata>>Ekin>>F){
    F=F/(Ekin/1000.0)/(Ekin/1000.0)/1000.0;
    np=gr->GetN();
    if(np==0)
      e_min = Ekin;
    gr->SetPoint(np,Ekin,F);
    gr->SetPointError(np,Ekin*0.10,F*0.10);
  }
  nn = gr->GetN();
  e_max = Ekin;
  indata.close();
}

//Differential energy spectrum of the Crab Nebula obtained with data recorded by the MAGIC stereoscopic system
//1/TeV x 1/cm^2 1/s
//From : 65.200 GeV
//To   : 21.500 TeV
Double_t log_parabola_Crab_Nebula_MAGIC(Double_t eGeV){
  Double_t A  = 3.23*1.0e-11; // (3.23 +/- 0.03)*1.0e-11
  Double_t B  = 1000.0;       // (GeV->TeV)
  Double_t C  = 2.47;         // (2.47 +/- 0.01)
  Double_t D  = 0.24;         // (0.24 +/- 0.01)
  Double_t val = A*(TMath::Power(eGeV/B, -C - D*TMath::Log10(eGeV/B)));
  return val;
}

Double_t modiﬁed_log_parabola(Double_t *x, Double_t *par){
  Double_t eGeV = x[0];
  Double_t A  = par[0];
  Double_t B  = par[1];
  Double_t C  = par[2];
  Double_t D  = par[3];
  Double_t Emax  = par[4];
  Double_t aa  = par[5];
  //Double_t val = A*(TMath::Power(eGeV/B, -C - D*TMath::Log10(TMath::Abs(eGeV/B))));
  //Double_t logE = TMath::Abs(TMath::Log10(TMath::Abs(eGeV/B/Emax)));
  Double_t logE = TMath::Abs(TMath::Log10(TMath::Abs(eGeV/B/Emax)));
  Double_t val = A*(TMath::Power(eGeV/B, -C + D*TMath::Power(logE,aa)));
  return val;
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
  return modiﬁed_log_parabola(eTeV*1000.0);
}
  
//Double_t modiﬁed_log_parabola_function_Crab_Nebula_MAGIC_Fermi_LAT(Double_t eGeV){
//}
//Double_t modiﬁed_log_parabola_function_Crab_Nebula_MAGIC_Fermi_LAT(Double_t eGeV){
//Double_t logf0 = -0.120;  // (  0.120 +/- 0.008)
//Double_t C     = -10.248; // (-10.248 +/- 0.006)
//Double_t EIC   = 48.0;    // ( 48.000 +/- 2.000)
//Double_t a     = 2.5;     // (  2.500 +/- 0.100)
//Double_t eTeV  = eGeV/1000.0;
//Double_t log_a = TMath::Power(TMath::Abs(TMath::Log10(eTeV/(EIC/1000.0))),a);
//cout<<TMath::Log10(eGeV/EIC)<<endl;
//Double_t val   = TMath::Power(10.0, logf0 + C/10*log_a);
//Double_t val   = logf0 + C*log_a;
//return val;
//}

void fill_gr_flux_crab_from_fit(TGraph *gr, Int_t nn, Double_t eGeV_min, Double_t eGeV_max){
  //
  Double_t log10_eGeV_min = TMath::Log10(eGeV_min);
  Double_t log10_eGeV_max = TMath::Log10(eGeV_max);  
  Double_t log10_eGeV;
  Double_t eGeV;
  Double_t F;
  //
  for(Int_t i = 0; i< nn;i++){
    log10_eGeV = log10_eGeV_min + (log10_eGeV_max-log10_eGeV_min)/(nn-1)*i;
    eGeV = TMath::Power(10.0,log10_eGeV);
    F = log_parabola_Crab_Nebula_MAGIC(eGeV);
    gr->SetPoint(gr->GetN(),eGeV,F);
  }  
}

void fill_gr_modiﬁed_log_parabola_function_Crab_Nebula_MAGIC_Fermi_LAT_fit(TGraph *gr, Int_t nn, Double_t eGeV_min, Double_t eGeV_max){
  //
  Double_t log10_eGeV_min = TMath::Log10(eGeV_min);
  Double_t log10_eGeV_max = TMath::Log10(eGeV_max);  
  Double_t log10_eGeV;
  Double_t eGeV;
  Double_t F;
  //
  for(Int_t i = 0; i< nn;i++){
    log10_eGeV = log10_eGeV_min + (log10_eGeV_max-log10_eGeV_min)/(nn-1)*i;
    eGeV = TMath::Power(10.0,log10_eGeV);
    //F = modiﬁed_log_parabola_function_Crab_Nebula_MAGIC_Fermi_LAT(eGeV);
    gr->SetPoint(gr->GetN(),eGeV,F);
    cout<<F<<endl;
  }  
}

void fill_gr_modiﬁed_log_parabola_fit(TGraph *gr, Int_t nn, Double_t eGeV_min, Double_t eGeV_max){
  //
  Double_t log10_eGeV_min = TMath::Log10(eGeV_min);
  Double_t log10_eGeV_max = TMath::Log10(eGeV_max);  
  Double_t log10_eGeV;
  Double_t eGeV;
  Double_t F;
  //
  for(Int_t i = 0; i< nn;i++){
    log10_eGeV = log10_eGeV_min + (log10_eGeV_max-log10_eGeV_min)/(nn-1)*i;
    eGeV = TMath::Power(10.0,log10_eGeV);
    F = modiﬁed_log_parabola(eGeV);
    gr->SetPoint(gr->GetN(),eGeV,F);
  }  
}

Double_t modiﬁed_log_parabola_TeV_integral(Double_t eTeV_min, Double_t eTeV_max){
  Int_t nn = 1000000;
  Double_t ee;
  Double_t dx = (eTeV_max - eTeV_min)/nn;
  Double_t dx_half = dx/2.0;
  Double_t integral = 0.0;
  for(Int_t i = 0; i<nn; i++){
    ee = eTeV_min + dx_half + dx*i;
    integral+=modiﬁed_log_parabola_TeV(ee)*1000.0*dx;
    //cout<<ee<<endl;
  }
  return integral;
}

Double_t modiﬁed_log_parabola_GeV_integral(Double_t eGeV_min, Double_t eGeV_max){
  Int_t nn = 1000000;
  Double_t ee;
  Double_t dx = (eGeV_max - eGeV_min)/nn;
  Double_t dx_half = dx/2.0;
  Double_t integral = 0.0;
  for(Int_t i = 0; i<nn; i++){
    ee = eGeV_min + dx_half + dx*i;
    integral+=modiﬁed_log_parabola(ee)*dx;
  }
  return integral;
}
