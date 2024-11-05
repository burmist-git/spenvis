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

void load_data_file(TString fileInName, TGraph *gr, Int_t &np, Double_t &eMin, Double_t &eMax, Double_t &integral);
void fill_hist(TH1D *h1, Int_t nn, Double_t eMin, Double_t eMax);
void extend_gr(TGraph *gr, TGraph *gr_ext, Int_t np, Double_t eMin, Double_t eMax);
void get_trg_eff(TH1D *h1_sim, TGraph *gr_trg, TGraph *gr_trg_eff);

Double_t modiﬁed_log_parabola(Double_t eGeV);
Double_t modiﬁed_log_parabola_TeV(Double_t eTeV);

void get_crab_gamma_rate(TGraph *gr_gamma_rates, TGraph *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV, Double_t gen_area_cmsq);

Double_t integrate_gr(TGraph *gr, Double_t eMin, Double_t eMax);

Int_t plots_LST_Crab_paper(){
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
  TGraph *gr_sim_ext = new TGraph();
  TGraph *gr_trg_ext = new TGraph();
  TGraph *gr_sel_ext = new TGraph();
  //
  Int_t np_gr_trg;
  Double_t eMin_gr_trg, eMax_gr_trg, integral_gr_trg;  
  //
  TGraph *gr_gamma_rates_trg = new TGraph();
  TGraph *gr_gamma_rates_trg_reco = new TGraph();
  TGraph *gr_trg_eff = new TGraph();
  //
  TH1D *h1_sim = new TH1D();
  TH1D *h1_sim_m2 = new TH1D();
  //
  load_data_file("./data/LST_Crab_paper_Zenith_23.63deg_sim.csv", gr_sim, np_sim, eMin_sim, eMax_sim, integral_sim);
  load_data_file("./data/LST_Crab_paper_Zenith_23.63deg_trg.csv", gr_trg, np_trg, eMin_trg, eMax_trg, integral_trg);
  load_data_file("./data/LST_Crab_paper_Zenith_23.63deg_sel.csv", gr_sel, np_sel, eMin_sel, eMax_sel, integral_sel);
  load_data_file("./data/LST_Crab_paper_gamma_rate_trg.csv", gr_gamma_rates_trg, np_gr_trg, eMin_gr_trg, eMax_gr_trg, integral_gr_trg);
  //
  extend_gr( gr_sim, gr_sim_ext, 100, eMin_sim, eMax_sim);
  extend_gr( gr_trg, gr_trg_ext, 100, eMin_trg, eMax_trg);
  extend_gr( gr_sel, gr_sel_ext, 100, eMin_sel, eMax_sel);
  //
  cout<<"np_sim       "<<np_sim<<endl
      <<"eMin_sim     "<<eMin_sim<<endl
      <<"eMax_sim     "<<eMax_sim<<endl
      <<"integral_sim "<<integral_sim<<endl;
  //
  cout<<"np_gr_trg       "<<np_gr_trg<<endl
      <<"eMin_gr_trg     "<<eMin_gr_trg<<endl
      <<"eMax_gr_trg     "<<eMax_gr_trg<<endl
      <<"integral_gr_trg "<<integral_gr_trg<<endl;
  //
  cout<<"integral "<< integrate_gr(gr_gamma_rates_trg, eMin_gr_trg, eMax_gr_trg)<<endl;
  cout<<"integral "<< integrate_gr(gr_gamma_rates_trg, 0.005, 0.02)<<endl;
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
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  //
  //gr_dif_550km_e->SetLineColor(kRed);
  //gr_dif_550km_e->SetLineWidth(3.0);
  //gr_dif_550km_e_norm->SetLineColor(kRed);
  //gr_dif_550km_e_norm->SetLineWidth(3.0);
  //
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr_sim);
  mg->Add(gr_trg);
  mg->Add(gr_sel);
  mg->Add(gr_sim_ext);
  mg->Add(gr_trg_ext);
  mg->Add(gr_sel_ext);
  //mg->SetMinimum(1.0e-2);
  //mg->SetMinimum(0.0);
  //mg->SetMaximum(1.1);
  //mg->SetMaximum(100.0);
  mg->Draw("APL");
  //mg->GetXaxis()->SetTitle("Energy, MeV");
  //mg->GetYaxis()->SetTitle("Differential Flux, 1/cm/cm/s/MeV");
  //mg->GetYaxis()->SetTitle("Relative differential flux");
  //
  //
  TCanvas *c2 = new TCanvas("c2","c2",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  c2->SetRightMargin(0.01);
  c2->SetLeftMargin(0.12);
  c2->SetTopMargin(0.01);
  c2->SetBottomMargin(0.08);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  Int_t nx=16;  
  Double_t *xBins = new Double_t[nx]; 
  Double_t log10_eTeV_min = TMath::Log10(eMin_TeV);
  Double_t log10_eTeV_max = TMath::Log10(eMax_TeV);  
  Double_t log10_eTeV;
  Double_t eTeV;
  for(Int_t i = 0; i < nx; i++){
    log10_eTeV = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(nx-1)*i;
    eTeV = TMath::Power(10.0,log10_eTeV);
    xBins[i]=eTeV;
  }
  h1_sim->SetBins(nx-1,xBins);
  h1_sim_m2->SetBins(nx-1,xBins);
  for(Int_t i = 0; i < (nx-1); i++){
    Double_t xx, yy;
    gr_sim->GetPoint(i,xx,yy);
    h1_sim->SetBinContent(i+1,yy);
  }
  h1_sim->Draw();  
  //
  //
  //
  fill_hist(h1_sim_m2, integral_sim, eMin_TeV, eMax_TeV);

  get_trg_eff( h1_sim_m2, gr_trg, gr_trg_eff);

  //
  TMultiGraph *mg2 = new TMultiGraph();
  mg2->Add(gr_sim);
  mg2->Draw("APL");
  h1_sim_m2->Draw("sames");


  TCanvas *c3 = new TCanvas("c3","c3",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  c3->SetRightMargin(0.01);
  c3->SetLeftMargin(0.12);
  c3->SetTopMargin(0.01);
  c3->SetBottomMargin(0.08);
  //
  //
  //gr_gamma_rates_trg->Draw("APL");
  get_crab_gamma_rate(gr_gamma_rates_trg_reco, gr_trg_eff, 70, 0.02, 0.3, TMath::Pi()*70000.0*70000.0);

  TMultiGraph *mg3 = new TMultiGraph();
  mg3->Add(gr_gamma_rates_trg);
  mg3->Add(gr_gamma_rates_trg_reco);
  mg3->Draw("APL");

  

  TCanvas *c4 = new TCanvas("c4","c4",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  c4->SetRightMargin(0.01);
  c4->SetLeftMargin(0.12);
  c4->SetTopMargin(0.01);
  c4->SetBottomMargin(0.08);
  //
  gr_trg_eff->Draw("APL");

  //
  //
  //
  //
  //
  //TLegend *leg = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  //leg->AddEntry(gr_dif_550km_e, "550 km","l");
  //leg->AddEntry(gr_dif_545km_e, "545 km","l");
  //leg->AddEntry(gr_dif_540km_e, "540 km","l");
  //leg->AddEntry(gr_dif_535km_e, "535 km","l");
  //leg->AddEntry(gr_dif_525km_e, "525 km","l");
  //leg->Draw();
  //
  /*
  //
  TCanvas *c2 = new TCanvas("c2","c2",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  //
  c2->SetRightMargin(0.01);
  c2->SetLeftMargin(0.12);
  c2->SetTopMargin(0.01);
  c2->SetBottomMargin(0.08);
  //
  gr_dif_550km_e_copy->SetLineColor(kBlack);
  gr_dif_550km_e_copy->SetLineWidth(3.0);
  gr_dif_550km_e_copy->Draw("APL");
  //
  gr_dif_550km_e_copy->GetXaxis()->SetTitle("Energy, MeV");
  gr_dif_550km_e_copy->GetYaxis()->SetTitle("Differential Flux, 1/cm/cm/s/MeV");

  */
  //
  //
  /*
  c1->Divide(2,1);
  //
  c1->cd(1);  
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  gr_al_mirror_refl_Unpolarized_nm->SetLineColor(kBlack);
  gr_al_mirror_refl_Unpolarized_nm->SetLineWidth(2.0);
  gr_al_mirror_refl_Spolarized_nm->SetLineColor(kRed);
  gr_al_mirror_refl_Spolarized_nm->SetLineWidth(2.0);
  gr_al_mirror_refl_Ppolarized_nm->SetLineColor(kBlue);
  gr_al_mirror_refl_Ppolarized_nm->SetLineWidth(2.0);
  //
  TMultiGraph *mg_nm = new TMultiGraph();
  mg_nm->Add(gr_al_mirror_refl_Unpolarized_nm);
  mg_nm->Add(gr_al_mirror_refl_Ppolarized_nm);
  mg_nm->Add(gr_al_mirror_refl_Spolarized_nm);
  mg_nm->SetMinimum(50.0);
  mg_nm->SetMaximum(100.0);
  mg_nm->Draw("AL");
  mg_nm->GetXaxis()->SetTitle("#lambda, nm");
  mg_nm->GetYaxis()->SetTitle("Reflectance, %");
  //
  TLegend *leg = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg->AddEntry(gr_al_mirror_refl_Unpolarized_nm, "Unpolarized","l");
  leg->AddEntry(gr_al_mirror_refl_Spolarized_nm, "S-Polarized","l");
  leg->AddEntry(gr_al_mirror_refl_Ppolarized_nm, "P-Polarized","l");
  leg->Draw();
  //
  c1->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  //
  gr_al_mirror_refl_Unpolarized_ev->SetLineColor(kBlack);
  gr_al_mirror_refl_Unpolarized_ev->SetLineWidth(2.0);
  gr_al_mirror_refl_Spolarized_ev->SetLineColor(kRed);
  gr_al_mirror_refl_Spolarized_ev->SetLineWidth(2.0);
  gr_al_mirror_refl_Ppolarized_ev->SetLineColor(kBlue);
  gr_al_mirror_refl_Ppolarized_ev->SetLineWidth(2.0);
  //
  TMultiGraph *mg_ev = new TMultiGraph();
  mg_ev->Add(gr_al_mirror_refl_Unpolarized_ev);
  mg_ev->Add(gr_al_mirror_refl_Ppolarized_ev);
  mg_ev->Add(gr_al_mirror_refl_Spolarized_ev);
  mg_ev->SetMinimum(50.0);
  mg_ev->SetMaximum(100.0);
  mg_ev->Draw("AL");
  mg_ev->GetXaxis()->SetTitle("Photon energy, eV");
  mg_ev->GetYaxis()->SetTitle("Reflectance, %");
  //
  leg->Draw();
  //
  save_data_file("al_mirror_Reflectance_Unpolarized_ev.dat", gr_al_mirror_refl_Unpolarized_ev);
  save_data_file("al_mirror_Reflectance_P-polarized_ev.dat", gr_al_mirror_refl_Ppolarized_ev);
  save_data_file("al_mirror_Reflectance_S-polarized_ev.dat", gr_al_mirror_refl_Spolarized_ev);
  //FILE *fp01;
  //fp01 = fopen(outData_ABSORB_vs_eV.Data(), "w+");
  //fprintf(fp01, " ev        m \n");
  //for(int i=0; i<71; i++) {
  //Double_t x_nm,x,y;
  //x_nm = 1000 - 10*i;
  //x = f_nm_to_eV->Eval(x_nm);
  //y = gr_abs_ev->Eval(x);
  //fprintf(fp01, "%10.7f %10.7f \n",x,y);
  //}
  //fclose(fp01);
  //
  return gr_al_mirror_refl_Unpolarized_ev;
  */
  return 0;
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

void fill_hist(TH1D *h1, Int_t nn, Double_t eMin, Double_t eMax){
  TRandom3 *rnd = new TRandom3(123123);
  for(Int_t i = 0;i<nn;i++)
    h1->Fill(1.0/(rnd->Uniform(1.0/eMax,1.0/eMin)));
}

void extend_gr(TGraph *gr, TGraph *gr_ext, Int_t np, Double_t eMin, Double_t eMax){
  Double_t log10_e_min = TMath::Log10(eMin);
  Double_t log10_e_max = TMath::Log10(eMax);  
  Double_t log10_e;
  Double_t e;
  for(Int_t i = 0; i < np; i++){
    log10_e = log10_e_min + (log10_e_max-log10_e_min)/(np-1)*i;
    e = TMath::Power(10.0,log10_e);
    gr_ext->SetPoint(i,e,gr->Eval(e));
  }
}

Double_t integrate_gr(TGraph *gr, Double_t eMin, Double_t eMax){
  Int_t np = 1000000;
  Double_t e;
  Double_t de = (eMax - eMin)/np;
  Double_t de_half = de/2.0;
  Double_t integral = 0.0;
  for(Int_t i = 0; i < np; i++){
    e = de_half + de*i;
    integral += gr->Eval(e)*de;
  }
  return integral;
}

void get_trg_eff(TH1D *h1_sim, TGraph *gr_trg, TGraph *gr_trg_eff){
  //
  Double_t eTeV;
  for(Int_t i = 1; i <= h1_sim->GetNbinsX(); i++){
    cout<<h1_sim->GetNbinsX()<<endl;
    eTeV=h1_sim->GetBinCenter(i);
    gr_trg_eff->SetPoint(gr_trg_eff->GetN(),eTeV,gr_trg->Eval(eTeV)/h1_sim->GetBinContent(i));
    //log10_eTeV = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(nx-1)*i;
    //eTeV = TMath::Power(10.0,log10_eTeV);
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
  Double_t B    = 1.00000;
  Double_t C    = 2.73558e+00;
  Double_t D    = 1.82440e-01;
  Double_t Emax = 2.20294e+01;
  Double_t aa   = 9.94035e-01;
  Double_t logE = TMath::Abs(TMath::Log10(TMath::Abs(eTeV/B/Emax)));
  Double_t val  = A*(TMath::Power(eTeV/B, -C + D*TMath::Power(logE,aa)));
  return val;
}

void get_crab_gamma_rate(TGraph *gr_gamma_rates, TGraph *gr_trg_eff, Int_t np, Double_t eMin_TeV, Double_t eMax_TeV, Double_t gen_area_cmsq){
  Double_t log10_eTeV_min = TMath::Log10(eMin_TeV);
  Double_t log10_eTeV_max = TMath::Log10(eMax_TeV);  
  Double_t log10_eTeV;
  Double_t eTeV;
  Double_t F;
  //
  for(Int_t i = 0; i< np;i++){
    log10_eTeV = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/np*i;
    eTeV = TMath::Power(10.0,log10_eTeV);
    //F = modiﬁed_log_parabola_TeV(eTeV)*(log10_eTeV_max-log10_eTeV_min)/np*gen_area_cmsq*gr_trg_eff->Eval(eTeV)*100;
    F = modiﬁed_log_parabola(eTeV*1000.0)*(log10_eTeV_max-log10_eTeV_min)/np*1000.0*gen_area_cmsq*gr_trg_eff->Eval(eTeV);
    gr_gamma_rates->SetPoint(gr_gamma_rates->GetN(),eTeV - (log10_eTeV_max-log10_eTeV_min)/np/2.0,F);
  }
}
