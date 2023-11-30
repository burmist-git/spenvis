//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

//
void get_TAUP2023_data(TGraphErrors *gr, TString name_dat_file);
//
Double_t log_parabola(Double_t *x, Double_t *par){
  Double_t A  = par[0];
  Double_t B  = par[1];
  Double_t C  = par[2];
  Double_t D  = par[3];
  Double_t val = A*(TMath::Power(x[0]/B,-C + D*TMath::Log10(x[0]/B)));
  return val;
}

Double_t power_law(Double_t *x, Double_t *par){
  Double_t A1  = par[0];
  Double_t B1  = par[1];
  Double_t C1  = par[2];
  Double_t A2 = par[3];
  Double_t B2 = par[4];
  Double_t C2 = par[5];
  Double_t val = A1*(TMath::Power(x[0]/B1,-C1)) + A2*(TMath::Power(x[0]/B2,-C2));
  return val;
}

Int_t log_parabola_fit_gamma_bkg_galactic(){
  //
  Int_t nn = 1000;
  Double_t e_min = 0.1;       // GeV
  Double_t e_max = 1000000.0; // GeV
  //
  TGraphErrors *gr_gamma_bkg_galactic = new TGraphErrors();
  gr_gamma_bkg_galactic->SetNameTitle("gr_gamma_bkg_galactic","gr_gamma_bkg_galactic");
  get_TAUP2023_data(gr_gamma_bkg_galactic,"./data/TAUP2023_ZhenCao_gamma_bkg_galactic.csv");
  //

  //Fit
  const Int_t npar = 6;
  //const Int_t npar = 3;
  Double_t inParameters[npar];
  //
  //1  p0           7.79234e-06   1.09762e-06   5.65735e-09  -1.71915e+02
  //2  p1           6.56007e+01   2.87790e+00   1.49294e-02  -3.76915e-05
  //3  p2           2.56356e+00   4.38212e-02  -1.90181e-05  -1.77122e-02
  //4  p3          -6.23527e-02   1.83214e-02  -8.33510e-05   4.66514e-02
  //
  // 1  p0          -2.19865e-09   6.60968e-12   3.57510e-16   8.21451e+07
  // 2  p1           3.87098e+03   3.49209e+00   2.50024e-04  -1.12088e-04
  // 3  p2           2.51752e+00   5.56789e-04   1.20045e-07  -5.60599e-01
  // 4  p3           1.93562e-08   5.91072e-11   3.09916e-15  -3.33282e+06
  // 5  p4           1.67809e+03   2.25870e+00   1.06288e-04  -1.10139e-04
  // 6  p5           2.52789e+00   5.27030e-04   1.20539e-07   5.07717e-01
  //
  // -2.19865e-09/(x/3.87098e+03)^2.51752e+00 + 1.93562e-08/(x/1.67809e+03)^2.52789e+00
  //
  //
  inParameters[0] = -2.20431e-09;
  inParameters[1] =  3.87496e+03;
  inParameters[2] =  2.49803e+00;
  inParameters[3] =  1.92880e-08;
  inParameters[4] =  1.67573e+03;
  inParameters[5] =  2.50983e+00;
  //
  //
  //TF1 *f_log_parabola = new TF1( "f_log_parabola", log_parabola, e_min, e_max, npar);
  TF1 *f_log_parabola = new TF1( "f_log_parabola", power_law, e_min, e_max, npar);
  f_log_parabola->SetParameters(inParameters);
  //f_log_parabola->FixParameter(0,inParameters[0]);
  //f_log_parabola->FixParameter(1,inParameters[1]);
  //f_log_parabola->FixParameter(2,inParameters[2]);
  //f_log_parabola->FixParameter(3,inParameters[3]);
  //f_log_parabola->FixParameter(6,inParameters[6]);
  gr_gamma_bkg_galactic->Fit("f_log_parabola","","",e_min, e_max);
  //
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
  gr_gamma_bkg_galactic->SetMarkerStyle(20);
  gr_gamma_bkg_galactic->SetMarkerColor(kBlack);
  gr_gamma_bkg_galactic->SetLineWidth(2);
  gr_gamma_bkg_galactic->SetLineColor(kBlack);
  //
  mg->Add(gr_gamma_bkg_galactic);
  //mg->Add(gr_tev_f);
  //
  mg->SetTitle("p0*1/(e/p1)^(-p2) + p3*1/(e/p4)^(-p5) + p6");
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
  mg->GetXaxis()->SetTitle("E, GeV");
  mg->GetYaxis()->SetTitle("dN/dE, m^{-2} s^{-1} sr^{-1} GeV^{-1}");
  
  return 0;
}

void get_TAUP2023_data(TGraphErrors *gr, TString name_dat_file){
  //
  Double_t Ekin;
  Double_t F;
  //
  Int_t ii = 0;
  //
  ifstream indata;
  indata.open(name_dat_file.Data());
  assert(indata.is_open());  
  while(indata>>Ekin>>F){
    //cout<<E<<endl;
    ii = gr->GetN();
    gr->SetPoint(ii,Ekin,F/Ekin/Ekin);
    gr->SetPointError(ii,Ekin/50.0,F/Ekin/Ekin/300.0);
  }
  indata.close();
}
