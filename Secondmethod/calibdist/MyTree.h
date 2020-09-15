//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 21 18:40:04 2013 by ROOT version 5.34/04
// from TTree MyTree/Ntuple
// found on file: StopSignal.root
//////////////////////////////////////////////////////////

#ifndef MyTree_h
#define MyTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MyTree {
public :
   TTree           *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // ---- Declaration of leaf types ----
   Int_t nmuon ;
   Int_t nneu ;
   Int_t nele ;
   Double_t mjjr;
   Double_t mjjt;
   Double_t edep;
   Double_t muene_sci;
   Double_t muene_che;
   Double_t emcomp1;
   Double_t emcomp2;
   Double_t chcomp1;
   Double_t chcomp2;
   Double_t chcomp1sm;
   Double_t chcomp2sm;
   Double_t elcomp1;
   Double_t elcomp2;
   Double_t neucomp1;
   Double_t neucomp2;
   Double_t etotjt1;
   Double_t etotjt2;
   Double_t eleak;
   Double_t eleakn;
   Double_t drmmu;
   Double_t enumu;
   Double_t mnumu;
   Double_t emu;
   Double_t mbos_noc;
   Double_t ej1b;
   Double_t ej2b;

   Double_t j1t_E;
   Double_t j1t_pt;
   Double_t j1t_eta;
   Double_t j1t_phi;
   Double_t j1t_m;
   Double_t j1t_theta;
   Double_t j2t_E;
   Double_t j2t_pt;
   Double_t j2t_eta;
   Double_t j2t_phi;
   Double_t j2t_m;
   Double_t j2t_theta;
   Double_t j1r_E;
   Double_t j1r_pt;
   Double_t j1r_eta;
   Double_t j1r_phi;
   Double_t j1r_m;
   Double_t j1r_theta;
   Double_t j2r_E;
   Double_t j2r_pt;
   Double_t j2r_eta;
   Double_t j2r_phi;
   Double_t j2r_m;
   Double_t j2r_theta;
   Double_t j1s_E;
   Double_t j1s_pt;
   Double_t j1s_eta;
   Double_t j1s_phi;
   Double_t j1s_m;
   Double_t j1s_theta;
   Double_t j2s_E;
   Double_t j2s_pt;
   Double_t j2s_eta;
   Double_t j2s_phi;
   Double_t j2s_m;
   Double_t j2s_theta;
   Double_t j1c_E;
   Double_t j1c_pt;
   Double_t j1c_eta;
   Double_t j1c_phi;
   Double_t j1c_m;
   Double_t j1c_theta;
   Double_t j2c_E;
   Double_t j2c_pt;
   Double_t j2c_eta;
   Double_t j2c_phi;
   Double_t j2c_m;
   Double_t j2c_theta;
//

   TBranch *b_nmuon ;
   TBranch *b_nneu ;
   TBranch *b_nele ;
   TBranch *b_mjjr;
   TBranch *b_mjjt;
   TBranch *b_edep;
   TBranch *b_muene_sci;
   TBranch *b_muene_che;
   TBranch *b_emcomp1;
   TBranch *b_emcomp2;
   TBranch *b_chcomp1;
   TBranch *b_chcomp2;
   TBranch *b_chcomp1sm;
   TBranch *b_chcomp2sm;
   TBranch *b_elcomp1;
   TBranch *b_elcomp2;
   TBranch *b_neucomp1;
   TBranch *b_neucomp2;
   TBranch *b_etotjt1;
   TBranch *b_etotjt2;
   TBranch *b_eleak;
   TBranch *b_eleakn;
   TBranch *b_drmmu;
   TBranch *b_enumu;
   TBranch *b_mnumu;
   TBranch *b_emu;
   TBranch *b_mbos_noc;
   TBranch *b_ej1b;
   TBranch *b_ej2b;

   TBranch *b_j1t_E;
   TBranch *b_j1t_pt;
   TBranch *b_j1t_eta;
   TBranch *b_j1t_phi;
   TBranch *b_j1t_m;
   TBranch *b_j1t_theta;
   TBranch *b_j2t_E;
   TBranch *b_j2t_pt;
   TBranch *b_j2t_eta;
   TBranch *b_j2t_phi;
   TBranch *b_j2t_m;
   TBranch *b_j2t_theta;
   TBranch *b_j1r_E;
   TBranch *b_j1r_pt;
   TBranch *b_j1r_eta;
   TBranch *b_j1r_phi;
   TBranch *b_j1r_m;
   TBranch *b_j1r_theta;
   TBranch *b_j2r_E;
   TBranch *b_j2r_pt;
   TBranch *b_j2r_eta;
   TBranch *b_j2r_phi;
   TBranch *b_j2r_m;
   TBranch *b_j2r_theta;
   TBranch *b_j1s_E;
   TBranch *b_j1s_pt;
   TBranch *b_j1s_eta;
   TBranch *b_j1s_phi;
   TBranch *b_j1s_m;
   TBranch *b_j1s_theta;
   TBranch *b_j2s_E;
   TBranch *b_j2s_pt;
   TBranch *b_j2s_eta;
   TBranch *b_j2s_phi;
   TBranch *b_j2s_m;
   TBranch *b_j2s_theta;
   TBranch *b_j1c_E;
   TBranch *b_j1c_pt;
   TBranch *b_j1c_eta;
   TBranch *b_j1c_phi;
   TBranch *b_j1c_m;
   TBranch *b_j1c_theta;
   TBranch *b_j2c_E;
   TBranch *b_j2c_pt;
   TBranch *b_j2c_eta;
   TBranch *b_j2c_phi;
   TBranch *b_j2c_m;
   TBranch *b_j2c_theta;

   MyTree();
   virtual ~MyTree();
   void     Init(); 
   virtual void     Write();
   void     Fill();
   void     Reset();
   
};

#endif

#ifdef MyTree_cxx
MyTree::MyTree() 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   //if (tree == 0) {
   /*
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("StopSignal.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("StopSignal.root");
      }
      f->GetObject("MyTree",tree);
   */
   //}
   Init();
}

MyTree::~MyTree()
{
}

#endif // #ifdef MyTree_cxx
