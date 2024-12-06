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

template <class T>
void get_spectrum( Double_t(*spectrum)(Double_t), T *gr, Int_t nn, Double_t ekin_TeV_min, Double_t ekin_TeV_max);

Double_t particles_spectrum_PDG_rate(Double_t ekin_TeV);
Double_t DAMPE_hard(Double_t ekin_TeV);
Double_t DAMPE_soft(Double_t ekin_TeV);
Double_t DAMPE_all(Double_t ekin_TeV);

Int_t smoothly_broken_power_law_SBPL_proton_spectrum(){
  //
  TGraphErrors *gr_DAMPE = new TGraphErrors();
  gr_DAMPE->SetNameTitle("gr_DAMPE","gr_DAMPE");
  TGraph *gr_spectrum_PDG = new TGraph();
  gr_spectrum_PDG->SetNameTitle("gr_spectrum_PDG","gr_spectrum_PDG");
  TGraph *gr_DAMPE_hard = new TGraph();
  gr_DAMPE_hard->SetNameTitle("gr_DAMPE_hard","gr_DAMPE_hard");
  TGraph *gr_DAMPE_soft = new TGraph();
  gr_DAMPE_soft->SetNameTitle("gr_DAMPE_soft","gr_DAMPE_soft");
  TGraph *gr_DAMPE_all = new TGraph();
  gr_DAMPE_all->SetNameTitle("gr_DAMPE_all","gr_DAMPE_all");
  //
  get_spectrum<TGraph>( particles_spectrum_PDG_rate, gr_spectrum_PDG,  100, 0.01, 100.0);
  get_spectrum<TGraph>( DAMPE_hard, gr_DAMPE_hard,  100, 0.01, 100.0);
  get_spectrum<TGraph>( DAMPE_soft, gr_DAMPE_soft,  100, 0.01, 100.0);
  get_spectrum<TGraph>( DAMPE_all, gr_DAMPE_all,  100, 0.01, 100.0);
  //
  get_DAMPE(gr_DAMPE,"./data/DAMPE.dat");
  //  
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1200,600);
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
  TMultiGraph *mg = new TMultiGraph();
  gr_DAMPE->SetMarkerStyle(20);
  gr_DAMPE->SetMarkerColor(kRed+2);
  gr_DAMPE->SetLineWidth(2);
  gr_DAMPE->SetLineColor(kRed+2);
  //
  mg->Add(gr_DAMPE);
  //mg->Add(gr_spectrum_PDG);
  //mg->Add(gr_DAMPE_hard);
  //mg->Add(gr_DAMPE_soft);
  mg->Add(gr_DAMPE_all);
  //
  mg->Draw("APL");
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(gr_DAMPE, "DAMPE Collaboration,Sci. Adv.2019;5:eaax3793", "apl");
  leg->AddEntry(gr_DAMPE_all, "DAMPE fit (all)", "apl");
  leg->Draw();  
  //
  mg->GetXaxis()->SetTitle("E, GeV");
  mg->GetYaxis()->SetTitle("dN/dE, m^{-2} s^{-1} sr^{-1} GeV^{-1}");
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
    Flux = F*frac;
    Flux_err = TMath::Sqrt(s_stat*s_stat*frac*frac + s_ana*s_ana*frac*frac + s_had*s_had*frac*frac);    
    gr->SetPoint(np,E,Flux);
    gr->SetPointError(np,(Emax - Emin)/2.0,Flux_err);
    np++;
    //cout<<mot<<endl;
  }
  indata.close();
}

//All particles spectrum from PDG
Double_t particles_spectrum_PDG_rate(Double_t ekin_TeV){
  Double_t result = 0.096*TMath::Power(ekin_TeV, -2.7);
  if(result < 0.0)
    result = 0.0;
  return result;
}

//Hard DAMPE spectrum (0.1 - 6.3 TeV) for protons
Double_t DAMPE_hard(Double_t ekin_TeV){
  Double_t result = 0.0758*TMath::Power(ekin_TeV, -2.772) * TMath::Power((1 + TMath::Power((ekin_TeV/0.48),5.0)), 0.0346);
  if(result < 0.0)
    result = 0.0;
  return result;
}

//Soft DAMPE spectrum (6.3 - 100 TeV) for protons
Double_t DAMPE_soft(Double_t ekin_TeV){
  Double_t result = 0.0868*TMath::Power(ekin_TeV, -2.6) * TMath::Power((1 + TMath::Power((ekin_TeV/13.6),5.0)), -0.05);
  if(result < 0.0)
    result = 0.0;
  return result;
}

Double_t DAMPE_all(Double_t ekin_TeV){
  if(ekin_TeV<6.3)
    return DAMPE_hard(ekin_TeV); 
  return DAMPE_soft(ekin_TeV);  
}

template <class T>
void get_spectrum( Double_t(*spectrum)(Double_t), T *gr, Int_t nn, Double_t ekin_TeV_min, Double_t ekin_TeV_max){
  //
  Double_t log10_ekin_TeV_min = TMath::Log10(ekin_TeV_min);
  Double_t log10_ekin_TeV_max = TMath::Log10(ekin_TeV_max);  
  Double_t log10_ekin_TeV;
  Double_t ekin_TeV;
  Double_t ekin_GeV;
  for( Int_t i = 0; i < nn; i++){
    log10_ekin_TeV = log10_ekin_TeV_min + (log10_ekin_TeV_max-log10_ekin_TeV_min)/(nn-1)*i;
    ekin_TeV = TMath::Power(10.0,log10_ekin_TeV);
    ekin_GeV = ekin_TeV*1000.0;
    gr->SetPoint(i,ekin_GeV,spectrum(ekin_TeV)/1000.0);
  }
}
