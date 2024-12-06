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

template <class T>
void load_data_file(TString fileInName, T *gr, Int_t &np, Double_t &eMin, Double_t &eMax, Double_t &integral);
template <class T>
void logalyze(T *gr, T *gr_log);

template <class T1, class T2>
void operation_graph( Double_t (*operation)(Double_t, Double_t), T1 *gr1, T1 *gr2, T2 *grs, Int_t n, Double_t xmin, Double_t xmax, bool if_linearScale = false, Double_t area = 1.0);
Double_t add(Double_t a, Double_t b) {return a + b;}
Double_t subtract(Double_t a, Double_t b) {return a - b;}
Double_t multiply(Double_t a, Double_t b) {return a * b;}
Double_t divide(Double_t a, Double_t b) { if (b != 0.0) return a / b; return 0.0;}
Double_t eff_area_f(Double_t *x, Double_t *par);

Double_t particles_spectrum_PDG_rate(Double_t ekin_TeV);
//Hard DAMPE spectrum (0.1 - 6.3 TeV) for protons
Double_t DAMPE_hard(Double_t ekin_TeV);
//Soft DAMPE spectrum (6.3 - 100 TeV) for protons
Double_t DAMPE_soft(Double_t ekin_TeV);
Double_t DAMPE_all(Double_t ekin_TeV);
template <class T>
void get_spectrum( Double_t(*spectrum)(Double_t), T *gr, Int_t nn, Double_t ekin_TeV_min, Double_t ekin_TeV_max, Double_t solid_angle);
//
template <class T>
void get_eff_fit( TF1 *f1log, T *gr, Int_t nn, Double_t ekin_TeV_min, Double_t ekin_TeV_max);
//
Double_t get_solid_angle(Double_t thetaDeg);
//
//
template <class T>
Double_t get_integral(T *gr, Double_t eMin_TeV, Double_t eMax_TeV);
//
//
//
Int_t plots_LST_proton_effArea_rate(){
  //
  Int_t np;
  Double_t eMin;
  Double_t eMax;
  Double_t integral;
  //
  TGraph *gr_sim = new TGraph();
  TGraph *gr_trg = new TGraph();
  TGraphErrors *gr_eff = new TGraphErrors();
  TGraph *gr_eff_fit = new TGraph();
  TGraph *gr_DAMPE_all = new TGraph();
  TGraphErrors *gr_proton_rate = new TGraphErrors();
  //
  TGraphErrors *gr_eff_log = new TGraphErrors();
  //
  load_data_file<TGraph>("../CTA/DBscan_on_simtel_data/LST_AdvCam_Zenith_20.00deg_proton_sim.csv", gr_sim, np, eMin, eMax, integral);
  cout<<"np       "<<np<<endl
      <<"eMin     "<<eMin<<endl
      <<"eMax     "<<eMax<<endl
      <<"integral "<<integral<<endl;
  load_data_file<TGraph>("../CTA/DBscan_on_simtel_data/LST_AdvCam_Zenith_20.00deg_proton_trg.csv", gr_trg, np, eMin, eMax, integral);
  cout<<"np       "<<np<<endl
      <<"eMin     "<<eMin<<endl
      <<"eMax     "<<eMax<<endl
      <<"integral "<<integral<<endl;
  //
  Double_t area = TMath::Pi()*1500.0*1500.0;
  operation_graph< TGraph, TGraphErrors>( divide, gr_trg, gr_sim, gr_eff, np, eMin, eMax, false, area);
  logalyze<TGraphErrors>(gr_eff, gr_eff_log);
  //
  //Fit
  const Int_t npar = 5;
  Double_t inParameters[npar];
  //
  inParameters[0] = -2.09;
  inParameters[1] = -1.49;
  inParameters[2] = 14.12;
  inParameters[3] =  3.03;  
  inParameters[4] = -7.05e-03;
  //
  TF1 *f_eff_area_f = new TF1("f_eff_area_f", eff_area_f, -1.97, 1.97, npar);
  f_eff_area_f->SetParameters(inParameters);
  f_eff_area_f->FixParameter(0,inParameters[0]);
  f_eff_area_f->FixParameter(1,inParameters[1]);
  f_eff_area_f->FixParameter(2,inParameters[2]);
  f_eff_area_f->FixParameter(3,inParameters[3]);
  f_eff_area_f->FixParameter(4,inParameters[4]);
  //
  gr_eff_log->Fit("f_eff_area_f","","",-1.97, 1.97);
  get_eff_fit( f_eff_area_f, gr_eff_fit, 100, 0.01, 100.0);
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
  TMultiGraph *mg01 = new TMultiGraph();
  mg01->Add(gr_sim);
  mg01->Add(gr_trg);
  mg01->SetMinimum(1.0e1);
  mg01->SetMaximum(1.0e8);
  mg01->Draw("APL");
  mg01->GetXaxis()->SetTitle("Energy, TeV");
  mg01->GetYaxis()->SetTitle("Number of events");
  //
  TLegend *leg01 = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg01->AddEntry(gr_sim, "Simulated (Zenith 20 deg)","l");
  leg01->AddEntry(gr_trg, "Triggered (Zenith 20 deg)","l");
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
  gr_eff->SetLineColor(kBlack);
  gr_eff->SetLineWidth(3.0);
  gr_eff->SetMarkerColor(kBlack);
  //
  TMultiGraph *mg02 = new TMultiGraph();
  mg02->Add(gr_eff);
  //mg02->Add(gr_trg_eff_AdvCam);
  //mg02->SetMinimum(1.0e-2);
  //mg02->SetMaximum(1.0e+0);
  mg02->Draw("APL");
  mg02->GetXaxis()->SetTitle("Energy, TeV");
  mg02->GetYaxis()->SetTitle("Trigger effective collection area, m^{2}");
  mg02->Draw("APL"); 
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
  //
  gr_eff_log->SetLineColor(kBlack);
  gr_eff_log->SetLineWidth(3.0);
  gr_eff_log->SetMarkerStyle(20);
  //
  gr_eff_log->Draw("APL");
  //
  c2->cd(2);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  //
  gr_eff_fit->SetLineColor(kBlue+2);
  gr_eff_fit->SetLineWidth(3.0);
  gr_eff_fit->SetMarkerColor(kBlue+2);
  //
  //
  TMultiGraph *mg03 = new TMultiGraph();
  mg03->Add(gr_eff);
  mg03->Add(gr_eff_fit);
  mg03->SetMinimum(1.0e1);
  mg03->SetMaximum(1.0e6);
  mg03->Draw("APL");
  mg03->GetXaxis()->SetTitle("Energy, TeV");
  mg03->GetYaxis()->SetTitle("Trigger effective collection area, m^{2}");
  //
  //
  //
  TCanvas *c3 = new TCanvas("c3","c3",10,10,700,700);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetOptStat(kFALSE);
  c3->SetRightMargin(0.01);
  c3->SetLeftMargin(0.12);
  c3->SetTopMargin(0.01);
  c3->SetBottomMargin(0.08);
  //
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetLogx();
  //
  //
  Double_t solid_angle = get_solid_angle(10.0);
  get_spectrum<TGraph>( DAMPE_all, gr_DAMPE_all, 100, 0.01, 100, solid_angle);
  operation_graph<TGraph,TGraphErrors>( multiply, gr_DAMPE_all, gr_eff, gr_proton_rate, 1000, 0.01, 100, false, 1.0);
  //
  gr_proton_rate->Draw("APL");
  //
  //
  //
  //
  cout<<get_integral(gr_proton_rate, 0.01, 100)<<endl;
  //
  return 0;
}

template <class T>
void load_data_file(TString fileInName, T *gr, Int_t &np, Double_t &eMin, Double_t &eMax, Double_t &integral){
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
      //gr->SetPointError(npoints_tot,ekin*0.05,TMath::Sqrt(ftot));
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

template <class T1, class T2>
void operation_graph( Double_t (*operation)(Double_t, Double_t), T1 *gr1, T1 *gr2, T2 *grs, Int_t n, Double_t xmin, Double_t xmax, bool if_linearScale, Double_t area){
  //
  Double_t log10_eTeV_min = TMath::Log10(xmin);
  Double_t log10_eTeV_max = TMath::Log10(xmax);  
  Double_t log10_eTeV_R;
  Double_t log10_eTeV_L;
  Double_t eTeV;
  Int_t npoints_tot = 0;
  Int_t i_zero = 1;
  Int_t i_last = n;
  Int_t ip = 0;
  if(if_linearScale){
    i_zero = 0;
    i_last = n;
  }
  //
  Double_t x, y1, y2, val;
  for( Int_t i = i_zero; i < i_last; i++){
    if(if_linearScale){
      x = xmin + (xmax-xmin)/(n-1)*i;
    }
    else{
      log10_eTeV_R = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(n-1)*i;
      log10_eTeV_L = log10_eTeV_min + (log10_eTeV_max-log10_eTeV_min)/(n-1)*(i-1);
      eTeV = TMath::Power(10.0,(log10_eTeV_R + log10_eTeV_L)/2.0);
      x = eTeV;
    }
    //
    //cout<<x<<endl;
    //  
    y1 = gr1->Eval(x);
    y2 = gr2->Eval(x);
    if(y2 != 0.0)
      val = operation(y1,y2);
    else
      val = 0.0;
    val = val*area;
    //cout<<val<<endl;
    ip=grs->GetN();
    grs->SetPoint(ip,x,val);
    TString typeid_t2 = typeid(TGraph).name();
    TString typeid_gre = typeid(TGraphErrors).name();
    if (typeid_t2 == typeid_gre){
      cout<<typeid(TGraphErrors).name()<<endl;
      cout<<typeid(TGraph).name()<<endl;
      grs->SetPointError(ip,0.0,val*0.1);
    }
  }
}

template <class T>
void logalyze(T *gr, T *gr_log){
  Double_t xx, yy;
  Double_t errxx = 0.00001;
  Double_t erryy = 0.00001;
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i,xx,yy);
    gr_log->SetPoint(i,TMath::Log10(xx),TMath::Log10(yy));
    if (typeid(T) == typeid(TGraphErrors))
      gr_log->SetPointError(i,errxx,erryy);
  }
}

// EXT PARAMETER
// NO.   NAME     VALUE      ERROR
// 1     p0      -2.09     9.55e-06
// 2     p1      -1.49     7.07e-06
// 3     p2      14.12     1.52e-05
// 4     p3       3.03     4.00e-06
// 5     p4      -7.05e-03 9.11e-06
Double_t eff_area_f(Double_t *x, Double_t *par){
  //
  Double_t A = par[0];
  Double_t B = par[1];
  Double_t C = par[2];
  Double_t D = par[3];
  Double_t E = par[4];
  Double_t x_truncation = E;
  Double_t K;
  Double_t BB;
  //
  Double_t val;
  Double_t val_L;
  Double_t val_R;
  Double_t xx_L;
  Double_t xx_R;
  Double_t dxx = 0.01;
  //
  if(x[0]<=x_truncation){
    val = A + B*x[0] + C*TMath::Log10(x[0]+D);
  }
  else{
    xx_R = x_truncation;
    xx_L = xx_R - dxx;
    //
    val_R = A + B*xx_R + C*TMath::Log10(xx_R+D);
    val_L = A + B*xx_L + C*TMath::Log10(xx_L+D);
    //
    K = (val_R - val_L)/(xx_R - xx_L);
    BB = val_R - K*xx_R;
    //
    val = BB + K*x[0];
  }
  //
  return val;
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

//Hard and soft DAMPE spectrum (0.1 - 100 TeV) for protons
Double_t DAMPE_all(Double_t ekin_TeV){
  if(ekin_TeV<6.3)
    return DAMPE_hard(ekin_TeV); 
  return DAMPE_soft(ekin_TeV);  
}

template <class T>
void get_spectrum( Double_t(*spectrum)(Double_t), T *gr, Int_t nn, Double_t ekin_TeV_min, Double_t ekin_TeV_max, Double_t solid_angle){
  //
  Double_t log10_ekin_TeV_min = TMath::Log10(ekin_TeV_min);
  Double_t log10_ekin_TeV_max = TMath::Log10(ekin_TeV_max);  
  Double_t log10_ekin_TeV;
  Double_t ekin_TeV;
  //Double_t ekin_GeV;
  for( Int_t i = 0; i < nn; i++){
    log10_ekin_TeV = log10_ekin_TeV_min + (log10_ekin_TeV_max-log10_ekin_TeV_min)/(nn-1)*i;
    ekin_TeV = TMath::Power(10.0,log10_ekin_TeV);
    //ekin_GeV = ekin_TeV*1000.0;
    //gr->SetPoint(i,ekin_GeV,spectrum(ekin_TeV)/1000.0);
    gr->SetPoint(i,ekin_TeV,spectrum(ekin_TeV)*solid_angle);
  }
}

template <class T>
void get_eff_fit( TF1 *f1log, T *gr, Int_t nn, Double_t ekin_TeV_min, Double_t ekin_TeV_max){
  //
  Double_t log10_ekin_TeV_min = TMath::Log10(ekin_TeV_min);
  Double_t log10_ekin_TeV_max = TMath::Log10(ekin_TeV_max);  
  Double_t log10_ekin_TeV;
  Double_t ekin_TeV;
  //Double_t ekin_GeV;
  for( Int_t i = 0; i < nn; i++){
    log10_ekin_TeV = log10_ekin_TeV_min + (log10_ekin_TeV_max-log10_ekin_TeV_min)/(nn-1)*i;
    ekin_TeV = TMath::Power(10.0,log10_ekin_TeV);
    gr->SetPoint(i,ekin_TeV,TMath::Power(10,f1log->Eval(log10_ekin_TeV)));
  }
}

Double_t get_solid_angle(Double_t thetaDeg){
  return 2*TMath::Pi()*(1.0-TMath::Cos(thetaDeg/180.0*TMath::Pi()));
}

template <class T>
Double_t get_integral(T *gr, Double_t eMin_TeV, Double_t eMax_TeV){
  Int_t nn = 1000000;
  Double_t ee;
  Double_t dx = (eMax_TeV - eMin_TeV)/nn;
  Double_t dx_half = dx/2.0;
  Double_t integral = 0.0;
  for(Int_t i = 0; i<nn; i++){
    ee = eMin_TeV + dx_half + dx*i;
    integral+=gr->Eval(ee)*dx;
  }
  return integral;
}
