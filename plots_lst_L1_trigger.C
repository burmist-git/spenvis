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
void load_data_file( TString fileInName, T *gr, Int_t &np, Double_t &thMin, Double_t &thMax, Double_t &integral);

Int_t plots_lst_L1_trigger(){
  //
  Int_t np;
  Double_t thMin;
  Double_t thMax;
  Double_t integral;
  //
  TGraph *gr_lst_L1_trigger_rate = new TGraph();
  TGraph *gr_P_DAMPE_lst_L1_trigger_rate = new TGraph();
  TGraph *gr_He_DAMPE_lst_L1_trigger_rate = new TGraph();
  //
  load_data_file<TGraph>( "./data/lst_L1_trigger_simulation_P_He_DAMPE.dat", gr_lst_L1_trigger_rate, np, thMin, thMax, integral);
  load_data_file<TGraph>("./data/P_DAMPE_lst_L1_trigger.dat", gr_P_DAMPE_lst_L1_trigger_rate, np, thMin, thMax, integral);
  load_data_file<TGraph>("./data/He_DAMPE_lst_L1_trigger.dat", gr_He_DAMPE_lst_L1_trigger_rate, np, thMin, thMax, integral);
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
  //c1->Divide(2,1);
  //
  //c1->cd(1);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  //gPad->SetLogx();
  //
  /*
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
  */
  //
  TMultiGraph *mg01 = new TMultiGraph();
  mg01->Add(gr_lst_L1_trigger_rate);
  mg01->Add(gr_P_DAMPE_lst_L1_trigger_rate);
  mg01->Add(gr_He_DAMPE_lst_L1_trigger_rate);
  //mg01->SetMinimum(1.0e3);
  //mg01->SetMaximum(1.0e8);
  mg01->Draw("APL");
  //mg01->GetXaxis()->SetTitle("Energy, TeV");
  //mg01->GetYaxis()->SetTitle("Number of events");
  //
  //
  cout<<"gr_lst_L1_trigger_rate->Eval(230)          "<<gr_lst_L1_trigger_rate->Eval(230)<<endl
      <<"gr_P_DAMPE_lst_L1_trigger_rate->Eval(230)  "<<gr_P_DAMPE_lst_L1_trigger_rate->Eval(230)<<endl
      <<"gr_He_DAMPE_lst_L1_trigger_rate->Eval(230) "<<gr_He_DAMPE_lst_L1_trigger_rate->Eval(230)<<endl;
  cout<<"gr_lst_L1_trigger_rate->Eval(270)          "<<gr_lst_L1_trigger_rate->Eval(270)<<endl
      <<"gr_P_DAMPE_lst_L1_trigger_rate->Eval(270)  "<<gr_P_DAMPE_lst_L1_trigger_rate->Eval(270)<<endl
      <<"gr_He_DAMPE_lst_L1_trigger_rate->Eval(270) "<<gr_He_DAMPE_lst_L1_trigger_rate->Eval(270)<<endl;
  //
  //
  /*
  TLegend *leg01 = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg01->AddEntry(gr_sim, "Simulated (Zenith 23.63 deg)","l");
  leg01->AddEntry(gr_trg, "Triggered (Zenith 23.63 deg)","l");
  leg01->AddEntry(gr_sel, "Selected  (Zenith 23.63 deg)","l");
  leg01->Draw();
  //
  //
  //
  */
  /*
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
  */
  //
  return 0;
}

//  x  y
//  9.769094138543517  344545.07758890773
// 19.538188277087034  190506.45979290613
// 29.307282415630553  112114.40064156917
// 39.964476021314390   73973.21462878345
template <class T>
void load_data_file(TString fileInName, T *gr, Int_t &np, Double_t &thMin, Double_t &thMax, Double_t &integral){
  ifstream fileIn(fileInName.Data());
  string mot;
  Double_t threshold;
  Double_t rate;
  Int_t npoints_tot = 0;
  integral = 0.0;
  if(fileIn.is_open()){
    fileIn>>mot>>mot;
    while(fileIn>>threshold>>rate){
      npoints_tot = gr->GetN();
      gr->SetPoint(npoints_tot,threshold,rate);
      if(gr->GetN()==1)
	thMin = threshold;
      //gr->SetPointError(npoints_tot,ekin*0.05,TMath::Sqrt(ftot));
      integral += rate;
      thMax = threshold;
    }
    fileIn.close();
  }
  else{
    cout<<"Unable to open file \n";
  }
  np=gr->GetN();
}
