#include <iostream>
#include <cmath>

#include "TCanvas.h"
#include "TTree.h"
#include "TCut.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TGraphErrors.h"
#include <sstream>
#include <string>

using namespace std;

void anchi(){
#include "PlotStyle.C+"
SetPlotStyle();
  TFile* file1 = TFile::Open("bonsai.hzjnbn30k.root","OLD");
//  
  TH2D* jetchi = new TH2D("jetchi", "jetchi", 160., -30.,30.,6,0.6,0.9);
  TH2D* jetchi1 = new TH2D("jetchi1", "jetchi1", 160., -30.,30.,6,0.6,0.9);
  TH2D* jetchi2 = new TH2D("jetchi2", "jetchi2", 160., -30.,30.,6,0.6,0.9);
  TH2D* jetchi3 = new TH2D("jetchi3", "jetchi3", 160., -30.,30.,6,0.6,0.9);
//
  TTree* tree1 = (TTree*) file1->Get("MyTree");
  string cut0s;
  cut0s="nmuon==0 && nneu==0 && abs(j1t_eta)<2 && abs(j2t_eta)<2  && abs(j1t_phi-j2t_phi)>0.1 && j1t_E+j2t_E>89.95 &&  j1r_E+j2r_E>85";
  cout << " string: " << endl;
  cout << cut0s << endl;
  TCut cut0 = cut0s.c_str();
//
  tree1->Project("jetchi", "j1c_E/j1s_E:(j1s_E-0.4*j1c_E)/(1-0.4)-j1t_E",cut0);
  tree1->Project("jetchi1", "j1c_E/j1s_E:(j1s_E-0.34*j1c_E)/(1-0.34)-j1t_E",cut0);
  tree1->Project("jetchi2", "j1c_E/j1s_E:(j1s_E-0.3*j1c_E)/(1-0.3)-j1t_E",cut0);
  tree1->Project("jetchi3", "j1c_E/j1s_E:(j1s_E-0.25*j1c_E)/(1-0.25)-j1t_E",cut0);
  tree1->Project("jetchi", "j2c_E/j2s_E:(j2s_E-0.4*j2c_E)/(1-0.4)-j2t_E",cut0);
  tree1->Project("jetchi1", "j2c_E/j2s_E:(j2s_E-0.34*j2c_E)/(1-0.34)-j2t_E",cut0);
  tree1->Project("jetchi2", "j2c_E/j2s_E:(j2s_E-0.3*j2c_E)/(1-0.3)-j2t_E",cut0);
  tree1->Project("jetchi3", "j2c_E/j2s_E:(j2s_E-0.25*j2c_E)/(1-0.25)-j2t_E",cut0);
//  tree1->Project("jetchi", "j2c_E/j2s_E:j2r_E-j2t_E",cut0);
  jetchi->FitSlicesX(0,0,-1,0);
  jetchi1->FitSlicesX(0,0,-1,0);
  jetchi2->FitSlicesX(0,0,-1,0);
  jetchi3->FitSlicesX(0,0,-1,0);
//
  TH1D *jetchi_1 = (TH1D*)gDirectory->Get("jetchi_1");
  TH1D *jetchi1_1 = (TH1D*)gDirectory->Get("jetchi1_1");
  TH1D *jetchi2_1 = (TH1D*)gDirectory->Get("jetchi2_1");
  TH1D *jetchi3_1 = (TH1D*)gDirectory->Get("jetchi3_1");
  jetchi_1->GetYaxis()->SetTitle("E(jrec)-E(jtru)");
  jetchi_1->GetXaxis()->SetTitle("Cerenkov/Scintillation");

  jetchi_1->SetMarkerStyle(20);
  jetchi_1->SetMarkerSize(1.2);
  jetchi1_1->SetMarkerStyle(21);
  jetchi1_1->SetMarkerSize(1.2);
  jetchi2_1->SetMarkerStyle(22);
  jetchi2_1->SetMarkerSize(1.2);
  jetchi3_1->SetMarkerStyle(23);
  jetchi3_1->SetMarkerSize(1.2);
//
//  TCanvas *c1 = new TCanvas("","",600.,650.);
  TLegend *legend = new TLegend(0.75,0.65,0.95,0.90);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->SetFillColor(0);
  legend->SetLineColor(0);
  legend->AddEntry(jetchi_1,"c=0.4","p");
  legend->AddEntry(jetchi1_1,"c=0.34","p");
  legend->AddEntry(jetchi2_1,"c=0.3","p");
  legend->AddEntry(jetchi3_1,"c=0.25","p");

  jetchi_1->SetMaximum(3.5);
  jetchi_1->SetMinimum(-4.5);
  jetchi_1->Draw();
  jetchi1_1->Draw("same");
  jetchi2_1->Draw("same");
  jetchi3_1->Draw("same");
  legend->Draw();
//  c1->Print("copt.jpg");
//  m_jetchi->Draw("colz");
}
