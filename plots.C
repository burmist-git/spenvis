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

void load_data_file(TString fileInName, TGraph *gr_tot, TGraph *gr_diff, Double_t &hMin, Double_t &hMax);
void norm_graph(TGraph *gr, TGraph *gr_norm);

Int_t plots(){
  //
  Double_t hMin;
  Double_t hMax;
  //
  TGraph *gr_tmp = new TGraph();
  //
  //max
  TGraph *gr_tot_max_e = new TGraph();
  TGraph *gr_dif_max_e = new TGraph();
  TGraph *gr_dif_max_e_norm = new TGraph();
  load_data_file("./data/max_e.dat", gr_tot_max_e, gr_dif_max_e,hMin,hMax);
  load_data_file("./data/max_e.dat", gr_tmp, gr_dif_max_e_norm,hMin,hMax);
  //
  //550
  TGraph *gr_tot_550km_e = new TGraph();
  TGraph *gr_dif_550km_e = new TGraph();
  TGraph *gr_dif_550km_e_norm = new TGraph();
  TGraph *gr_dif_550km_e_copy = new TGraph();
  load_data_file("./data/550_km_e.dat", gr_tot_550km_e, gr_dif_550km_e,hMin,hMax);
  load_data_file("./data/550_km_e.dat", gr_tmp, gr_dif_550km_e_norm,hMin,hMax);
  load_data_file("./data/550_km_e.dat", gr_tmp, gr_dif_550km_e_copy,hMin,hMax);
  norm_graph(gr_dif_550km_e_norm, gr_dif_550km_e);
  norm_graph(gr_dif_max_e_norm, gr_dif_550km_e);
  //
  //545
  TGraph *gr_tot_545km_e = new TGraph();
  TGraph *gr_dif_545km_e = new TGraph();
  TGraph *gr_dif_545km_e_norm = new TGraph();
  load_data_file("./data/545_km_e.dat", gr_tot_545km_e, gr_dif_545km_e,hMin,hMax);
  load_data_file("./data/545_km_e.dat", gr_tmp, gr_dif_545km_e_norm,hMin,hMax);
  norm_graph(gr_dif_545km_e_norm, gr_dif_550km_e);
  //
  //540
  TGraph *gr_tot_540km_e = new TGraph();
  TGraph *gr_dif_540km_e = new TGraph();
  TGraph *gr_dif_540km_e_norm = new TGraph();
  load_data_file("./data/540_km_e.dat", gr_tot_540km_e, gr_dif_540km_e,hMin,hMax);
  load_data_file("./data/540_km_e.dat", gr_tmp, gr_dif_540km_e_norm,hMin,hMax);
  norm_graph(gr_dif_540km_e_norm, gr_dif_550km_e);
  //
  //535
  TGraph *gr_tot_535km_e = new TGraph();
  TGraph *gr_dif_535km_e = new TGraph();
  TGraph *gr_dif_535km_e_norm = new TGraph();
  load_data_file("./data/535_km_e.dat", gr_tot_535km_e, gr_dif_535km_e,hMin,hMax);
  load_data_file("./data/535_km_e.dat", gr_tmp, gr_dif_535km_e_norm,hMin,hMax);
  norm_graph(gr_dif_535km_e_norm, gr_dif_550km_e);
  //
  //525
  TGraph *gr_tot_525km_e = new TGraph();
  TGraph *gr_dif_525km_e = new TGraph();
  TGraph *gr_dif_525km_e_norm = new TGraph();
  load_data_file("./data/525_km_e.dat", gr_tot_525km_e, gr_dif_525km_e,hMin,hMax);
  load_data_file("./data/525_km_e.dat", gr_tmp, gr_dif_525km_e_norm,hMin,hMax);
  norm_graph(gr_dif_525km_e_norm, gr_dif_550km_e);
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
  TMultiGraph *mg = new TMultiGraph();
  //mg->Add(gr_dif_max_e);
  mg->Add(gr_dif_550km_e_norm);
  mg->Add(gr_dif_545km_e_norm);
  mg->Add(gr_dif_540km_e_norm);
  mg->Add(gr_dif_535km_e_norm);
  mg->Add(gr_dif_525km_e_norm);
  //mg->Add(gr_dif_max_e_norm);
  //mg->SetMinimum(1.0e-2);
  mg->SetMinimum(0.0);
  mg->SetMaximum(1.1);
  //mg->SetMaximum(100.0);
  mg->Draw("AL");
  mg->GetXaxis()->SetTitle("Energy, MeV");
  //mg->GetYaxis()->SetTitle("Differential Flux, 1/cm/cm/s/MeV");
  mg->GetYaxis()->SetTitle("Relative differential flux");

  gr_dif_max_e->SetLineColor(kBlack);
  gr_dif_max_e->SetLineWidth(3.0);

  gr_dif_550km_e->SetLineColor(kRed);
  gr_dif_550km_e->SetLineWidth(3.0);
  gr_dif_550km_e_norm->SetLineColor(kRed);
  gr_dif_550km_e_norm->SetLineWidth(3.0);

  gr_dif_545km_e->SetLineColor(kGreen+2);
  gr_dif_545km_e->SetLineWidth(2.0);
  gr_dif_545km_e_norm->SetLineColor(kGreen+2);
  gr_dif_545km_e_norm->SetLineWidth(2.0);

  gr_dif_540km_e->SetLineColor(kMagenta+2);
  gr_dif_540km_e->SetLineWidth(2.0);
  gr_dif_540km_e_norm->SetLineColor(kMagenta+2);
  gr_dif_540km_e_norm->SetLineWidth(2.0);

  gr_dif_535km_e->SetLineColor(kBlue+2);
  gr_dif_535km_e->SetLineWidth(2.0);
  gr_dif_535km_e_norm->SetLineColor(kBlue+2);
  gr_dif_535km_e_norm->SetLineWidth(2.0);

  gr_dif_525km_e->SetLineColor(kBlue);
  gr_dif_525km_e->SetLineWidth(3.0);
  gr_dif_525km_e_norm->SetLineColor(kBlue);
  gr_dif_525km_e_norm->SetLineWidth(3.0);

  TLegend *leg = new TLegend(0.7,0.2,0.9,0.4,"","brNDC");
  leg->AddEntry(gr_dif_550km_e, "550 km","l");
  leg->AddEntry(gr_dif_545km_e, "545 km","l");
  leg->AddEntry(gr_dif_540km_e, "540 km","l");
  leg->AddEntry(gr_dif_535km_e, "535 km","l");
  leg->AddEntry(gr_dif_525km_e, "525 km","l");
  leg->Draw();


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

void norm_graph(TGraph *gr, TGraph *gr_norm){
  Double_t x,y;
  Double_t norm;
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i,x,y);
    gr_norm->GetPoint(i,x,norm);
    if(norm>0.0)
      gr->SetPoint(i,x,y/norm);
  }
}
