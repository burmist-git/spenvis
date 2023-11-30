//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

void get_DAMPE(TGraphErrors *gr, TString name_dat_file);
void get_DAMPE_fit(TGraph *gr);
void get_flux_p_PDG(TGraph *gr, TString name_dat_file);
void get_TAUP2023_data(TGraph *gr, TString name_dat_file);
void get_MAGIC_crab_data(TGraph *gr, TString name_dat_file);
void get_MAGIC_crab_data_digi(TGraph *gr, TString name_dat_file);
Double_t function_log_parabola_fit_MAGIC_tel_crab_GeV(Double_t e);
Double_t function_log_parabola_fit_MAGIC_tel_crab_GeV(Double_t e);
void get_crab_log_parabola(TGraph *gr);
Double_t get_proton_diff_flux_DAMPE(Double_t e);
void get_ratio_crab_DAMPE(TGraph *gr);
void get_ratio_crab_epem(TGraph *gr, TGraph *gr_TAUP2023_epem);

Double_t acceptance_frac = 2*TMath::Pi()*(1.0-TMath::Cos(0.5/180.0*TMath::Pi()));
const Int_t nn_DAMPE_crab_fit = 300;

Int_t TAUP2023_ZhenCao_2deg(){
  //
  TString fileN;
  //
  const Int_t n = 10000;
  Double_t I[n];
  Double_t etot[n];
  //
  Double_t etot_min = 1;      //in GeV (5 GeV)
  Double_t etot_max = 1000000; //in GeV (100 TeV)
  //
  Double_t N = 1;
  //
  TGraphErrors *gr_DAMPE = new TGraphErrors();
  gr_DAMPE->SetNameTitle("gr_DAMPE","gr_DAMPE");
  //
  TGraphErrors *gr_r_crab_DAMPE_fit = new TGraphErrors();
  gr_r_crab_DAMPE_fit->SetNameTitle("gr_r_crab_DAMPE_fit","gr_r_crab_DAMPE_fit");
  //
  TGraphErrors *gr_r_crab_epem = new TGraphErrors();
  gr_r_crab_epem->SetNameTitle("gr_r_crab_epem","gr_r_crab_epem");
  //
  TGraphErrors *gr_DAMPE_fit = new TGraphErrors();
  gr_DAMPE_fit->SetNameTitle("gr_DAMPE_fit","gr_DAMPE_fit");  
  //
  TGraphErrors *gr_DAMPE_He = new TGraphErrors();
  gr_DAMPE_He->SetNameTitle("gr_DAMPE_He","gr_DAMPE_He");
  //
  TGraph *gr_flux_p_PDG = new TGraph();
  gr_flux_p_PDG->SetNameTitle("gr_flux_p_PDG","gr_flux_p_PDG");  
  //
  TGraphErrors *gr_TAUP2023_p = new TGraphErrors();
  gr_TAUP2023_p->SetNameTitle("gr_TAUP2023_p","gr_TAUP2023_p");
  //
  TGraphErrors *gr_TAUP2023_he = new TGraphErrors();
  gr_TAUP2023_he->SetNameTitle("gr_TAUP2023_he","gr_TAUP2023_he");
  //
  TGraphErrors *gr_TAUP2023_O = new TGraphErrors();
  gr_TAUP2023_O->SetNameTitle("gr_TAUP2023_O","gr_TAUP2023_O");
  //
  TGraphErrors *gr_TAUP2023_Fe = new TGraphErrors();
  gr_TAUP2023_Fe->SetNameTitle("gr_TAUP2023_Fe","gr_TAUP2023_Fe");
  //
  TGraphErrors *gr_TAUP2023_gamma_bkg_galactic = new TGraphErrors();
  gr_TAUP2023_gamma_bkg_galactic->SetNameTitle("gr_TAUP2023_gamma_bkg_galactic","gr_TAUP2023_gamma_bkg_galactic");
  //
  TGraphErrors *gr_TAUP2023_gamma_bkg_extragalactic = new TGraphErrors();
  gr_TAUP2023_gamma_bkg_extragalactic->SetNameTitle("gr_TAUP2023_gamma_bkg_extragalactic","gr_TAUP2023_gamma_bkg_extragalactic");
  //
  TGraphErrors *gr_TAUP2023_epem = new TGraphErrors();
  gr_TAUP2023_epem->SetNameTitle("gr_TAUP2023_epem","gr_TAUP2023_epem");
  //
  TGraphErrors *gr_TAUP2023_ep = new TGraphErrors();
  gr_TAUP2023_ep->SetNameTitle("gr_TAUP2023_ep","gr_TAUP2023_ep");
  //
  TGraphErrors *gr_gamma_crab_digi = new TGraphErrors();
  gr_gamma_crab_digi->SetNameTitle("gr_gamma_crab_digi","gr_gamma_crab_digi");
  //
  TGraphErrors *gr_gamma_crab = new TGraphErrors();
  gr_gamma_crab->SetNameTitle("gr_gamma_crab","gr_gamma_crab");
  //
  TGraphErrors *gr_gamma_crab_par = new TGraphErrors();
  gr_gamma_crab_par->SetNameTitle("gr_gamma_crab_par","gr_gamma_crab_par");
  //
  get_DAMPE(gr_DAMPE,"./data/DAMPE.dat");
  get_DAMPE_fit(gr_DAMPE_fit);
  get_DAMPE(gr_DAMPE_He,"./data/DAMPE_He.dat");
  get_flux_p_PDG(gr_flux_p_PDG, "./data/flux_P_PDG.csv");  
  //
  get_TAUP2023_data(gr_TAUP2023_p, "data/TAUP2023_ZhenCao_p.csv");
  get_TAUP2023_data(gr_TAUP2023_he, "data/TAUP2023_ZhenCao_he.csv");
  get_TAUP2023_data(gr_TAUP2023_O, "data/TAUP2023_ZhenCao_O.csv");
  get_TAUP2023_data(gr_TAUP2023_Fe, "data/TAUP2023_ZhenCao_Fe.csv");
  get_TAUP2023_data(gr_TAUP2023_gamma_bkg_galactic, "data/TAUP2023_ZhenCao_gamma_bkg_galactic.csv");
  get_TAUP2023_data(gr_TAUP2023_gamma_bkg_extragalactic, "data/TAUP2023_ZhenCao_gamma_bkg_extragalactic.csv");
  get_TAUP2023_data(gr_TAUP2023_epem, "data/TAUP2023_ZhenCao_e+e-.csv");
  get_TAUP2023_data(gr_TAUP2023_ep, "data/TAUP2023_ZhenCao_e+.csv");
  //
  get_MAGIC_crab_data(gr_gamma_crab, "data/MAGICdata.txt");
  get_MAGIC_crab_data_digi(gr_gamma_crab_digi, "data/MAGICdata_crab_digi.txt");
  //
  get_crab_log_parabola(gr_gamma_crab_par);
  //
  //
  get_ratio_crab_DAMPE(gr_r_crab_DAMPE_fit);
  //
  //
  get_ratio_crab_epem(gr_r_crab_epem, gr_TAUP2023_epem);
  //
  //
  TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,600,600);
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
  gr_DAMPE->SetMarkerStyle(20);
  gr_DAMPE->SetMarkerColor(kBlack+2);
  gr_DAMPE->SetLineWidth(2);
  gr_DAMPE->SetLineColor(kBlack+2);
  //
  gr_DAMPE_fit->SetMarkerStyle(20);
  gr_DAMPE_fit->SetMarkerColor(kBlack+2);
  gr_DAMPE_fit->SetLineWidth(2);
  gr_DAMPE_fit->SetLineColor(kBlack+2);
  //
  gr_DAMPE_He->SetMarkerStyle(20);
  gr_DAMPE_He->SetMarkerColor(kRed+2);
  gr_DAMPE_He->SetLineWidth(2);
  gr_DAMPE_He->SetLineColor(kRed+2);
  //
  gr_flux_p_PDG->SetLineWidth(2);
  gr_flux_p_PDG->SetMarkerColor(kBlue+2);
  gr_flux_p_PDG->SetLineColor(kBlue+2);
  //
  gr_TAUP2023_p->SetLineWidth(2);
  gr_TAUP2023_p->SetMarkerColor(kBlue+3);
  gr_TAUP2023_p->SetLineColor(kBlue+3);
  //
  gr_TAUP2023_he->SetLineWidth(2);
  gr_TAUP2023_he->SetMarkerColor(kRed+3);
  gr_TAUP2023_he->SetLineColor(kRed+3);
  //  
  gr_TAUP2023_O->SetLineWidth(2);
  gr_TAUP2023_O->SetMarkerColor(kRed);
  gr_TAUP2023_O->SetLineColor(kRed);
  //
  gr_TAUP2023_Fe->SetLineWidth(2);
  gr_TAUP2023_Fe->SetMarkerColor(kBlue);
  gr_TAUP2023_Fe->SetLineColor(kBlue);
  //
  gr_TAUP2023_gamma_bkg_galactic->SetMarkerStyle(20);
  gr_TAUP2023_gamma_bkg_galactic->SetLineWidth(2);
  gr_TAUP2023_gamma_bkg_galactic->SetMarkerColor(kMagenta+2);
  gr_TAUP2023_gamma_bkg_galactic->SetLineColor(kMagenta+2);
  //
  gr_TAUP2023_gamma_bkg_extragalactic->SetMarkerStyle(21);
  gr_TAUP2023_gamma_bkg_extragalactic->SetLineWidth(2);
  gr_TAUP2023_gamma_bkg_extragalactic->SetMarkerColor(kMagenta+3);
  gr_TAUP2023_gamma_bkg_extragalactic->SetLineColor(kMagenta+3);
  //
  gr_TAUP2023_epem->SetMarkerStyle(21);
  gr_TAUP2023_epem->SetLineWidth(2);
  gr_TAUP2023_epem->SetMarkerColor(kRed+1);
  gr_TAUP2023_epem->SetLineColor(kRed+1);
  //
  gr_TAUP2023_ep->SetMarkerStyle(20);
  gr_TAUP2023_ep->SetLineWidth(2);
  gr_TAUP2023_ep->SetMarkerColor(kBlue+1);
  gr_TAUP2023_ep->SetLineColor(kBlue+1);
  //
  gr_gamma_crab_par->SetMarkerStyle(22);
  gr_gamma_crab_par->SetLineWidth(3);
  gr_gamma_crab_par->SetMarkerColor(kMagenta);
  gr_gamma_crab_par->SetLineColor(kMagenta);
  //
  mg->Add(gr_DAMPE);
  //mg->Add(gr_DAMPE_fit);
  mg->Add(gr_DAMPE_He);
  mg->Add(gr_flux_p_PDG);
  mg->Add(gr_TAUP2023_p);
  mg->Add(gr_TAUP2023_he);
  mg->Add(gr_TAUP2023_O);
  mg->Add(gr_TAUP2023_Fe);
  mg->Add(gr_TAUP2023_gamma_bkg_galactic);
  mg->Add(gr_TAUP2023_gamma_bkg_extragalactic);
  mg->Add(gr_TAUP2023_epem);
  mg->Add(gr_TAUP2023_ep);
  //
  //mg->Add(gr_gamma_crab_digi);
  mg->Add(gr_gamma_crab_par);
  //mg->Add(gr_gamma_crab);
  //
  mg->Draw("APL");
  //
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
  leg->AddEntry(gr_gamma_crab_par, "gamma from crab MAGIC tel.", "apl");
  //
  leg->Draw();
  //
  mg->GetXaxis()->SetTitle("Ekin, GeV");
  mg->GetYaxis()->SetTitle("dN/dE, m^{-2} s^{-1} GeV^{-1}");
  //
  //
  //
  //
  //  
  TCanvas *c2 = new TCanvas("c2",fileN.Data(),10,10,1200,600);
  c2->Divide(2,1);
  c2->cd(1);
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
  TMultiGraph *mg02 = new TMultiGraph();
  //
  gr_DAMPE_fit->SetMarkerStyle(20);
  gr_DAMPE_fit->SetMarkerColor(kBlack);
  gr_DAMPE_fit->SetLineWidth(2);
  gr_DAMPE_fit->SetLineColor(kBlack);
  //
  gr_gamma_crab_par->SetMarkerStyle(22);
  gr_gamma_crab_par->SetLineWidth(3);
  gr_gamma_crab_par->SetMarkerColor(kMagenta+2);
  gr_gamma_crab_par->SetLineColor(kMagenta+2);
  //
  mg02->Add(gr_DAMPE_fit);
  mg02->Add(gr_gamma_crab_par);
  //
  mg02->Draw("APL");
  //
  mg02->GetXaxis()->SetTitle("Ekin, GeV");
  mg02->GetYaxis()->SetTitle("dN/dE, m^{-2} s^{-1} GeV^{-1}");
  //
  c2->cd(2);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  gPad->SetLogy();
  //  
  gr_r_crab_DAMPE_fit->SetTitle("");
  gr_r_crab_DAMPE_fit->Draw("APL");
  gr_r_crab_DAMPE_fit->GetXaxis()->SetTitle("Ekin, GeV");
  gr_r_crab_DAMPE_fit->GetYaxis()->SetTitle("dN/dE, m^{-2} s^{-1} GeV^{-1}");
  //
  //
  //
  //
  TCanvas *c3 = new TCanvas("c3",fileN.Data(),10,10,1200,600);
  c3->Divide(2,1);
  c3->cd(1);
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
  TMultiGraph *mg03 = new TMultiGraph();
  //
  gr_DAMPE_fit->SetMarkerStyle(20);
  gr_DAMPE_fit->SetMarkerColor(kBlack);
  gr_DAMPE_fit->SetLineWidth(2);
  gr_DAMPE_fit->SetLineColor(kBlack);
  //
  gr_gamma_crab_par->SetMarkerStyle(22);
  gr_gamma_crab_par->SetLineWidth(3);
  gr_gamma_crab_par->SetMarkerColor(kMagenta+2);
  gr_gamma_crab_par->SetLineColor(kMagenta+2);
  //
  mg03->Add(gr_TAUP2023_epem);
  mg03->Add(gr_gamma_crab_par);
  //
  mg03->Draw("APL");
  //
  mg03->GetXaxis()->SetTitle("Ekin, GeV");
  mg03->GetYaxis()->SetTitle("dN/dE, m^{-2} s^{-1} GeV^{-1}");
  //
  //
  c3->cd(2);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  gPad->SetLogy();
  //  
  gr_r_crab_epem->SetTitle("");
  gr_r_crab_epem->Draw("APL");
  gr_r_crab_epem->GetXaxis()->SetTitle("Ekin, GeV");
  gr_r_crab_epem->GetYaxis()->SetTitle("dN/dE, m^{-2} s^{-1} GeV^{-1}");
  //
  //
  //
  return 0;
}

//E     Emin  Emax     F    s_stat s_ana s_had (GeV-1_m-2_s-1_sr-1)
//49.8  39.8  63.1     2.97 0.00 0.14 0.20  1.0e-1
void get_DAMPE(TGraphErrors *gr, TString name_dat_file){
  //
  Double_t E;
  Double_t Emin;
  Double_t Emax;
  Double_t F;
  Double_t s_stat;
  Double_t s_ana;
  Double_t s_had;
  Double_t frac;
  Double_t Flux;
  Double_t Flux_err;
  //
  Int_t np = 0;
  //
  ifstream indata;
  string mot;
  indata.open(name_dat_file.Data());
  assert(indata.is_open());  
  for(Int_t i = 0;i<8;i++)
    indata  >> mot;
  while(indata>>E>>Emin>>Emax>>F>>s_stat>>s_ana>>s_had>>frac){
    //cout<<E<<endl;
    Flux = F*frac*acceptance_frac;
    Flux_err = TMath::Sqrt(s_stat*s_stat*frac*frac + s_ana*s_ana*frac*frac + s_had*s_had*frac*frac);    
    gr->SetPoint(np,E,Flux);
    gr->SetPointError(np,(Emax - Emin)/2.0,Flux_err*acceptance_frac);
    np++;
    //cout<<mot<<endl;
  }
  indata.close();
}

void get_flux_p_PDG(TGraph *gr, TString name_dat_file){
  //
  Double_t Ekin;
  Double_t F;
  //  
  ifstream indata;
  indata.open(name_dat_file.Data());
  assert(indata.is_open());  
  while(indata>>Ekin>>F){
    //cout<<E<<endl;
    gr->SetPoint(gr->GetN(),Ekin,F*acceptance_frac);
  }
  indata.close();
}

void get_TAUP2023_data(TGraph *gr, TString name_dat_file){
  //
  Double_t Ekin;
  Double_t F;
  //  
  ifstream indata;
  indata.open(name_dat_file.Data());
  assert(indata.is_open());
  while(indata>>Ekin>>F){
    //cout<<E<<endl;
    gr->SetPoint(gr->GetN(),Ekin,F/Ekin/Ekin*acceptance_frac);
  }
  indata.close();
}

//dN/dEdAdt(TeV-1cm-2s-1) Energy(GeV) Energy errors
//1.2e-8     6.52e1  10
void get_MAGIC_crab_data(TGraph *gr, TString name_dat_file){
  //
  Double_t Ekin;
  Double_t F;
  string mot;
  //  
  ifstream indata;
  indata.open(name_dat_file.Data());
  assert(indata.is_open());  
  indata>>mot>>mot>>mot>>mot;
  while(indata>>F>>Ekin>>mot){
    cout<<Ekin<<endl;
    gr->SetPoint(gr->GetN(),Ekin,F*10000.0/(Ekin/1000.0)/(Ekin/1000.0));
  }
  indata.close();
}

void get_MAGIC_crab_data_digi(TGraph *gr, TString name_dat_file){
  //
  Double_t Ekin;
  Double_t F;
  string mot;
  //  
  ifstream indata;
  indata.open(name_dat_file.Data());
  assert(indata.is_open());  
  while(indata>>Ekin>>F){
    //cout<<E<<endl;
    gr->SetPoint(gr->GetN(),Ekin*1000.0,F*10000.0/Ekin/Ekin);
  }
  indata.close();
}

Double_t function_log_parabola_fit_MAGIC_tel_crab_GeV(Double_t e){
  //return 10000.0*((3.31327)*1.0e-14)*TMath::Power((e/1.00729/1000.0),(-2.49570 - 0.130707*TMath::Log(e/1.00729/1000.0)));
  return 10000.0*((3.23)*1.0e-14)*TMath::Power((e/1.00000/1000.0),(-2.47 - 0.24*TMath::Log(e/1.00000/1000.0)));
}

Double_t cpf(Double_t *x, Double_t *par){
  Double_t A  = par[0];
  Double_t B  = par[1];
  Double_t C  = 0.0;
  Double_t val = A/(TMath::Power(x[0],B)) + C;
  return val;
}

void get_DAMPE_fit(TGraph *gr){
  Double_t e_min = 50.0;
  Double_t e_max = 3000.0;
  Double_t e;
  Int_t i = 0;
  Int_t nn = nn_DAMPE_crab_fit;
  for(Int_t i = 0;i<nn;i++){
    e = (e_max - e_min)/(nn-1)*i + e_min;
    gr->SetPoint(gr->GetN(),e,get_proton_diff_flux_DAMPE(e)*acceptance_frac);
  }
}

Double_t get_proton_diff_flux_DAMPE(Double_t e){
  return 9.67871e+03*TMath::Power(e,-2.674);
}

void get_crab_log_parabola(TGraph *gr){
  Double_t e_min = 50.0;
  Double_t e_max = 3000.0;
  Double_t e;
  Int_t i = 0;
  Int_t nn = nn_DAMPE_crab_fit;
  for(Int_t i = 0;i<nn;i++){
    e = (e_max - e_min)/(nn-1)*i + e_min;
    gr->SetPoint(gr->GetN(),e,function_log_parabola_fit_MAGIC_tel_crab_GeV(e));
  }
}

void get_ratio_crab_DAMPE(TGraph *gr){
  Double_t e_min = 50.0;
  Double_t e_max = 3000.0;
  Double_t e;
  Int_t i = 0;
  Int_t nn = nn_DAMPE_crab_fit;
  for(Int_t i = 0;i<nn;i++){
    e = (e_max - e_min)/(nn-1)*i + e_min;
    gr->SetPoint(gr->GetN(),e,function_log_parabola_fit_MAGIC_tel_crab_GeV(e)/(get_proton_diff_flux_DAMPE(e)*acceptance_frac));
  }
}

void get_ratio_crab_epem(TGraph *gr, TGraph *gr_TAUP2023_epem){
  Double_t e_min = 50.0;
  Double_t e_max = 3000.0;
  Double_t e;
  Double_t f;
  Int_t i = 0;
  for(Int_t i = 0;i<gr_TAUP2023_epem->GetN();i++){
    gr_TAUP2023_epem->GetPoint(i,e,f);
    if( e >= e_min && e <= e_max)
      gr->SetPoint(gr->GetN(),e,function_log_parabola_fit_MAGIC_tel_crab_GeV(e)/f);
  }
}
