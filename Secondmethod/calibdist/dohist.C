/***********************************************************************
	This is the main program to process the CBNT Root Ntuple from
Athena with SUSYtup. See SUSYtup.h for more information.
	For version 7.0.0++ using name susy for Ntuple.
***********************************************************************/
#include <TROOT.h> 
#include "TTree.h" 
#include "TBranch.h" 
#include <TFile.h>
#include "MyTree.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <tuple>
#include <TRandom3.h>
using namespace std;
//
//  setup for particle charge
//
TDatabasePDG *fPDG;
TParticlePDG *pdgParticle;

TFile* ftree;
MyTree bonsaiTree;

vector<fastjet::PseudoJet> inputparticles_tru;
vector<fastjet::PseudoJet> inputparticles_scin;
vector<fastjet::PseudoJet> inputparticles_cher;
vector<fastjet::PseudoJet> inputparticles_all;
vector<fastjet::PseudoJet> jetexc;
vector<fastjet::PseudoJet> jet_all;
vector<fastjet::PseudoJet> jet_scin;
vector<fastjet::PseudoJet> jet_cher_t;
vector<fastjet::PseudoJet> jet_cher;
vector<fastjet::PseudoJet> jet_rec;
vector<fastjet::PseudoJet> jet_tru;
vector<fastjet::PseudoJet>  sci_comp;
vector<fastjet::PseudoJet> cher_comp;
vector<fastjet::PseudoJet> rec_comp;
vector<TLorentzVector> jet_comp;
vector<TLorentzVector> muvec;
vector<TLorentzVector> nuvec;
vector<TLorentzVector> numuvec;

vector<fastjet::PseudoJet> trackghost;
vector<fastjet::PseudoJet> truthghost;
vector<fastjet::PseudoJet> tracks;

vector <double> emcomp;
vector <double> elcomp;
vector <double> chcomp;
vector <double> chcompsm;
vector <double> neucomp;
vector <double> etotjt;

vector<double> Calib_VectorScinR;
vector<double> Calib_VectorScinL;
vector<double> Calib_VectorCherR;
vector<double> Calib_VectorCherL;	

TRandom3   m_random3;

double GeV=1000.;
double pi=3.14159265;
// threshold for including particle/cell in pseudojet
double threshold=0.01;
double threshold_c=0.01;
// eta limit of calorimeter
double etalim=5; //5.0
//double chi=0.46;
//////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) { 
  vector<double> calibscin(std::vector<double> vectorscin);
  vector<double> calibcher(std::vector<double> vectorcher);
  tuple<double, double, double> maptower(int index, string side);
  fastjet::PseudoJet matchjet(fastjet::PseudoJet jet_in, vector<fastjet::PseudoJet> testvec);
  fastjet::PseudoJet mergejet(fastjet::PseudoJet jet_scin, fastjet::PseudoJet jet_cher, double); 
  double ptsmeatr(double, double);
  pair<double,double>ejconst(double k, double th1, double th2, double mass1, double mass2); 
  fPDG = TDatabasePDG::Instance();

  cout << "-------------------------------------------" << std::endl;
  cout << " N arguments " << argc << std::endl;
  if(argc<2){
    cout << "Please give output.root input.root" << endl;
    return 0;
  }
  // Output tree file
  string histName=argv[1];
  std::cout << " Output file name: " << std::endl;
  std::cout << "      " << histName << std::endl;

  // Open input file
  std::string fn = argv[2];
  std::cout << " Read file name: " << std::endl;
  std::cout << "      " << fn << std::endl;

  std::string chis = argv[3];
  std::cout << " chi value * 1000" << std::endl;
  std::cout << "      " << chis << std::endl;
  double chi=atof(chis.c_str())/1000.;
  cout << " chi " << chi << endl;

  std::string filetru = fn+"_truth.root";
  std::string filesim = fn+".root";
  std::cout<< "Reading file " <<filetru <<std::endl;
  TFile* f = new TFile( filetru.c_str() );
  std::cout<< "Reading file " <<filesim <<std::endl;
  TFile* f1 = new TFile( filesim.c_str() );

  int pos_st=histName.rfind("/");
  string stma=histName.substr(pos_st+1);
  string newfile="bonsai."+stma;
  cout << " bonsai file " << newfile << endl;
  ftree = new TFile(newfile.c_str(), "RECREATE");
  bonsaiTree.Init();
#include "truthdec.h"
#include "B4dec.h"
//
  TTree* tree1 = (TTree*)f->Get("truth");
  TTree* tree2 = (TTree*)f1->Get("B4");
  tree1->AddFriend(tree2);
#include "truthset.h"
#include "B4set.h"
  tree1->GetEntry(0);
//
  if (tree1== 0) return 1;
  Int_t nentries = Int_t(tree1->GetEntries());
  Int_t nbytes= 0, nb = 0;
//
  cout << " Number of events " << nentries << endl;
  for (Int_t jentry=0; jentry<nentries;jentry++) {
//    cout << " jentry " << jentry << endl;
    nb = tree1->GetEntry(jentry);   nbytes += nb;
    bonsaiTree.Reset();
    inputparticles_tru.clear();
    muvec.clear();
    nuvec.clear();
    tracks.clear();
    trackghost.clear();
    truthghost.clear();
    int nmuon=0;
    int nele=0;
    int nneu=0;
    int nmun=0;
    double muene_sci=0.;
    double muene_che=0.;
    double etott=0;
    double ecalt=0;
    int ncha=0;
//  loop on truth particles
    for(uint itru=0;itru<mcs_n;itru++){
      int partid = mcs_pdgId->at(itru);
      pdgParticle = fPDG->GetParticle(partid);
      float charge=pdgParticle ? int(pdgParticle->Charge()/3.0) : -999;
      double parteta = mcs_eta->at(itru);
      etott+=mcs_E->at(itru);
//    only select particles fully absorbed in calo
//    reject muons, neutrinos and neutralinos
      if(abs(partid) ==11)nele++;
      if(abs(partid) != 13 &&  
         abs(partid) !=12  && abs(partid) != 14 && abs(partid) != 16 &&
         abs(partid) != 1000022){
//    require that they are within calo rapidity	 
        if(abs(parteta)<etalim){
          TLorentzVector trup;
          trup.SetPtEtaPhiM(mcs_pt->at(itru), mcs_eta->at(itru), 
			    mcs_phi->at(itru), mcs_m->at(itru));
//    zero the masses of particles for consistency with cells
	  double escal=trup.E()/trup.P();
//    build vector for stnadlaone truth jet reco
          fastjet::PseudoJet fj(trup.Px()*escal, trup.Py()*escal, 
			        trup.Pz()*escal, trup.E());
          ecalt+=mcs_E->at(itru);
          fj.set_user_index(itru);
          inputparticles_tru.push_back(fj);
//   define pseudojet of ghost particle for 
//   association of truth particles to jets
          fastjet::PseudoJet fjr(trup.Px()*1e-18, 
	                         trup.Py()*1e-18, 
	                 	 trup.Pz()*1e-18, trup.E()*1e-18);
          fjr.set_user_index(itru+200000);
//   store for ghosts only tru particles with energy 
//   above the threshold used for calo cells
          if(mcs_E->at(itru)>threshold)truthghost.push_back(fjr);
//   create vectors of charged tracks, full and
//   scaled to ghost for possible use
          if(abs(charge)>0.0001) {
            tracks.push_back(fj);
            if(mcs_E->at(itru)>threshold)trackghost.push_back(fjr);
            ncha++;
          }
        }
      }
//    end if on absorbed particles
//    count neutrinos and save vector
      if(abs(partid) ==12  ||  abs(partid) == 14 || abs(partid) == 16){
        TLorentzVector nuall;
        nuall.SetPtEtaPhiM(mcs_pt->at(itru), mcs_eta->at(itru), 
			  mcs_phi->at(itru), mcs_m->at(itru));
        nuvec.push_back(nuall);
	nneu++;
      }
//    build vector of muons
      if(abs(partid) == 13){
        TLorentzVector muon;
        muon.SetPtEtaPhiM(mcs_pt->at(itru), mcs_eta->at(itru), 
			  mcs_phi->at(itru), mcs_m->at(itru));
        muvec.push_back(muon);
        nmuon++;
      }
//    build vector of muon neutrinos 
      if(abs(partid) == 14){
        TLorentzVector numu;
        numu.SetPtEtaPhiM(mcs_pt->at(itru), mcs_eta->at(itru), 
			  mcs_phi->at(itru), mcs_m->at(itru));
        numuvec.push_back(numu);
        nmun++;
      }
    } // loop on truth particles    
//  create standalone truth jets    
//    jetexc.clear();
//    fastjet::JetDefinition jet_def(fastjet::ee_genkt_algorithm, 2.*pi, 1.);
//    fastjet::ClusterSequence clust_seq(inputparticles_tru, jet_def); 
//    jetexc = clust_seq.exclusive_jets(int(2));
//
//  now the rec part
//
    Calib_VectorScinR.clear();
    Calib_VectorScinL.clear();
    Calib_VectorCherR.clear();
    Calib_VectorCherL.clear();

//  fill vectors of cells calibrated at em scale for 
//  both scintillator and cerenkov component  
    Calib_VectorScinR = calibscin(*VectorSignalsR);
    Calib_VectorScinL = calibscin(*VectorSignalsL);
    Calib_VectorCherR = calibcher(*VectorSignalsCherR);
    Calib_VectorCherL = calibcher(*VectorSignalsCherL);

    double energy=0;
    double energyc=0;
//
//  calculate total energy 
//
    for(uint i=0; i<Calib_VectorScinR.size(); i++) {
      energy+=Calib_VectorScinR.at(i)+Calib_VectorScinL.at(i);
      energyc+=Calib_VectorCherR.at(i)+Calib_VectorCherL.at(i);
    } 
//  create  pseudojet vectors from calo cells
    if(energy>0){
      inputparticles_scin.clear();
      inputparticles_cher.clear();
      inputparticles_all.clear();
      TLorentzVector scin_bos(0.,0.,0.,0.);
      TLorentzVector cher_bos(0.,0.,0.,0.);
// right side
      for(int towerindex=1; towerindex<=75*36; towerindex++) {
        auto thphieta=maptower(towerindex, "right");
        double theta=get<0>(thphieta);
        double phi=get<1>(thphieta);
        double eta=get<2>(thphieta);
        double energy_scin = Calib_VectorScinR[towerindex];
        double pt_scin = energy_scin*sin(theta*pi/180.);

        double energy_cher = Calib_VectorCherR[towerindex];
        double pt_cher = energy_cher*sin(theta*pi/180.);

        TLorentzVector towerscin;
        towerscin.SetPtEtaPhiM(pt_scin, eta, phi*pi/180., 0.);
        TLorentzVector towercher;
        towercher.SetPtEtaPhiM(pt_cher, eta, phi*pi/180., 0.);
//  for each cell find the minimum distance with a muon
        double deltamumin=999999.;
        for(uint i=0;i<muvec.size();i++) {
           double deltaR=abs(towerscin.DeltaR(muvec[i]));
           if(deltaR<deltamumin)deltamumin=deltaR;
        }
// build pseudojet vector for scintillator  
        if(energy_scin > threshold) {
          if(deltamumin<0.1){
            muene_sci = muene_sci+towerscin.E();
          }
          if(deltamumin>0.1){
            scin_bos=scin_bos+towerscin;	
            fastjet::PseudoJet cellsci(towerscin.Px(), towerscin.Py(), 
			               towerscin.Pz(), towerscin.E());
            cellsci.set_user_index(towerindex);
            inputparticles_scin.push_back(cellsci);
            inputparticles_all.push_back(cellsci);
          }
        }  
// build pseudojet vector for cerenkov
        if(energy_cher > threshold_c) {
          if(deltamumin<0.1){
            muene_che = muene_che+towercher.E();
          }
          if(deltamumin>0.1){
            cher_bos=cher_bos+towercher;
            fastjet::PseudoJet cellcher(towercher.Px(), towercher.Py(), 
			                towercher.Pz(), towercher.E());
            cellcher.set_user_index(towerindex+100000);
            inputparticles_cher.push_back(cellcher);
            inputparticles_all.push_back(cellcher);
          }
        }  
      }
// left side
      for(int towerindex=1; towerindex<=75*36; towerindex++) {
        auto thphieta=maptower(towerindex, "left");
        double theta=get<0>(thphieta);
        double phi=get<1>(thphieta);
        double eta=get<2>(thphieta);
        double energy_scin = Calib_VectorScinL[towerindex];
        double pt_scin = energy_scin*sin(theta*pi/180.);

        double energy_cher = Calib_VectorCherL[towerindex];
        double pt_cher = energy_cher*sin(theta*pi/180.);

        TLorentzVector towerscin;
        towerscin.SetPtEtaPhiM(pt_scin, eta, phi*pi/180., 0.);
        TLorentzVector towercher;
        towercher.SetPtEtaPhiM(pt_cher, eta, phi*pi/180., 0.);
        double deltamumin=999999.;
        for(uint i=0;i<muvec.size();i++) {
           double deltaR=abs(towerscin.DeltaR(muvec[i]));
           if(deltaR<deltamumin)deltamumin=deltaR;
        }
        if(energy_scin > threshold) {
          if(deltamumin<0.1){
            muene_sci = muene_sci+towerscin.E();
          }
          if(deltamumin>0.1){
            scin_bos=scin_bos+towerscin;	
            fastjet::PseudoJet cellsci(towerscin.Px(), towerscin.Py(), 
			               towerscin.Pz(), towerscin.E());
            cellsci.set_user_index(towerindex);
            inputparticles_scin.push_back(cellsci);
            inputparticles_all.push_back(cellsci);
          }
        }  
        if(energy_cher > threshold_c) {
          if(deltamumin<0.1){
            muene_che = muene_che+towercher.E();
          }
          if(deltamumin>0.1){
            cher_bos=cher_bos+towercher;
            fastjet::PseudoJet cellcher(towercher.Px(), towercher.Py(), 
			                towercher.Pz(), towercher.E());
            cellcher.set_user_index(towerindex+100000);
            inputparticles_cher.push_back(cellcher);
            inputparticles_all.push_back(cellcher);
          }
        }  
      }
      fastjet::PseudoJet merge_bos =mergejet(
            fastjet::PseudoJet(scin_bos.Px(), scin_bos.Py(), 
		               scin_bos.Pz(), scin_bos.E()),
            fastjet::PseudoJet(cher_bos.Px(), cher_bos.Py(), 
		               cher_bos.Pz(), cher_bos.E()), chi);

      double mbos_noc=merge_bos.m();
//
//  add ghost truth
//
      for(uint ig=0;ig<truthghost.size();ig++) {
         inputparticles_scin.push_back(truthghost.at(ig));
         inputparticles_all.push_back(truthghost.at(ig));
      }

      fastjet::JetDefinition jet_defs(fastjet::ee_genkt_algorithm, 2.*pi, 1.);
//  create jet using the full set of cells
      fastjet::ClusterSequence clust_seq_all(inputparticles_all, jet_defs);
//  clear vectors of jets
      jet_all.clear();
      jet_rec.clear();
      jet_tru.clear();
//   create vector of jet_all 
      jet_all.push_back(clust_seq_all.exclusive_jets(int(2))[0]);
      jet_all.push_back(clust_seq_all.exclusive_jets(int(2))[1]);
//
//   muiso for W
//
     double drminmu=9999.;
/*
     if(muvec.size()>0) {	      
       for(uint jt=0; jt<jet_tru.size(); jt++) {
          TLorentzVector jj;
          jj.SetPtEtaPhiM(jet_tru[jt].pt(), jet_tru[jt].eta(), jet_tru[jt].phi(),
                        jet_tru[jt].m());
          double deltaR=abs(jj.DeltaR(muvec[0]));
	  if(deltaR<drminmu)drminmu=deltaR;
       }
     }
*/
//
//  tru jets from ghost components
//  and calculation of charged, em, electron, and hadneutral 
//  from truth partilcles
//  for charged components also perform smearing based 
//  on official delphes card
// 
     emcomp.clear();
     elcomp.clear();
     chcomp.clear();
     chcompsm.clear();
     neucomp.clear();
     etotjt.clear();
     jet_comp.clear();
     sci_comp.clear();
     cher_comp.clear();
     rec_comp.clear();
//
//   loop on jets built with with sci and cher cells
//   and rebuild sci and cher jets 
//
     for(uint jn=0; jn<jet_all.size();jn++) {
        vector<fastjet::PseudoJet> constituents = jet_all[jn].constituents();
	TLorentzVector jcomp;
	fastjet::PseudoJet scicomp(0.,0.,0.,0.);
	fastjet::PseudoJet chercomp(0.,0.,0.,0.);
//
        double eem=0;
        double eel=0;
        double ech=0;
        double echsm=0;
        double eneu=0;
        double etotj=0;
//  loop on components
        for (unsigned j = 0; j < constituents.size(); j++) {
          int ui=constituents[j].user_index();
// rebuild scintillator jets
          if(ui<100000) {
            fastjet::PseudoJet scij(constituents[j].px(), constituents[j].py(),
			            constituents[j].pz(), constituents[j].e());
            scicomp+=scij; 
          } 
// rebuild cerenkov jets
          if(ui>=100000 && ui<200000) {
            fastjet::PseudoJet cherj(constituents[j].px(), constituents[j].py(),                                     constituents[j].pz(), constituents[j].e());
            chercomp+=cherj; 
          }
// build truth jets and store different components of jet
          if(ui>=200000) {
            int itru=ui-200000;
            TLorentzVector trup;
            trup.SetPtEtaPhiM(mcs_pt->at(itru),  mcs_eta->at(itru), 
			      mcs_phi->at(itru), mcs_m->at(itru));
            jcomp+=trup; 
            int partid = mcs_pdgId->at(itru);
            etotj+=mcs_E->at(itru);
            if(abs(partid)==22 || abs(partid)==11) eem+=mcs_E->at(itru);
            if(abs(partid)==11) eel+=mcs_E->at(itru);
            pdgParticle = fPDG->GetParticle(partid);
            float charge=pdgParticle ? int(pdgParticle->Charge()/3.0) : -999;
            if(abs(charge)>0.01) {
              ech+=mcs_E->at(itru);
// save charge component smeared as for official Delphes datacard
              double pptr=ptsmeatr(mcs_pt->at(itru), mcs_eta->at(itru));
	      echsm+=pptr*mcs_E->at(itru)/mcs_pt->at(itru);
	    }
            if(abs(charge)<0.01 && abs(partid) !=22) eneu+=mcs_E->at(itru);
          }
        }
	jet_comp.push_back(jcomp);
	sci_comp.push_back(scicomp);
	cher_comp.push_back(chercomp);
        rec_comp.push_back(mergejet(scicomp,chercomp,chi));
        emcomp.push_back(eem);
        elcomp.push_back(eel);
        chcomp.push_back(ech);
        chcompsm.push_back(echsm);
        neucomp.push_back(eneu);
        etotjt.push_back(etotj);
     }	
      double emu=0.;
      double enumu=0; 
      double mnumu=0;
      if(nmuon==1)emu=muvec[0].E();
      if(nmuon==1 && nmun==1) {
        enumu=muvec[0].E()+nuvec[0].E();
        mnumu=(muvec[0]+nuvec[0]).M();
      }
//
//    save in ntuple for two-jets case
//
      if(jet_comp.size()==2 && rec_comp.size()==2) {
        fastjet::PseudoJet jetrec=rec_comp[0]+rec_comp[1];
        TLorentzVector jettruth=jet_comp[0]+jet_comp[1];
	double kk=round(jet_comp[0].E()+jet_comp[1].E());
	double th1=rec_comp[0].theta();
	double th2=rec_comp[1].theta();
	double mass1=rec_comp[0].m();
	double mass2=rec_comp[1].m();
	auto ejbeam = ejconst(kk, th1, th2, mass1, mass2);
	double ej1beam=ejbeam.first;
        double ej2beam=ejbeam.second;
//	cout << " kk " << kk << " j1b " << ej1beam << " " << jet_comp[0].E() <<
//	        " j2b " << ej2beam << " " << jet_comp[1].E() << endl;
//
	bonsaiTree.nmuon = nmuon;
        bonsaiTree.nneu = nneu;
        bonsaiTree.nele = nele;
        bonsaiTree.mjjr= jetrec.m();
        bonsaiTree.mjjt= jettruth.M();
        bonsaiTree.edep= EnergyTot/1000.;
        bonsaiTree.muene_sci=muene_sci;
        bonsaiTree.muene_che=muene_che;
	bonsaiTree.emcomp1=emcomp[0];
	bonsaiTree.emcomp2=emcomp[1];
	bonsaiTree.chcomp1=chcomp[0];
	bonsaiTree.chcomp2=chcomp[1];
	bonsaiTree.chcomp1sm=chcompsm[0];
	bonsaiTree.chcomp2sm=chcompsm[1];
	bonsaiTree.neucomp1=neucomp[0];
	bonsaiTree.neucomp2=neucomp[1];
	bonsaiTree.elcomp1=elcomp[0];
	bonsaiTree.elcomp2=elcomp[1];
	bonsaiTree.etotjt1=etotjt[0];
	bonsaiTree.etotjt2=etotjt[1];
	bonsaiTree.eleak=leakage/1000.;
	bonsaiTree.eleakn=neutrinoleakage/10000.;
	bonsaiTree.mbos_noc=mbos_noc;
	bonsaiTree.drmmu=drminmu;
	bonsaiTree.enumu=enumu;
	bonsaiTree.mnumu=mnumu;
	bonsaiTree.emu=emu;
	bonsaiTree.ej1b=ej1beam;
	bonsaiTree.ej2b=ej2beam;
//	
//      four-momenta of jets for further processing 
//      on ntuples
//
//      truth jets from ghost components
//
        bonsaiTree.j1t_E=jet_comp[0].E();
        bonsaiTree.j1t_pt=jet_comp[0].Pt();
        bonsaiTree.j1t_eta=jet_comp[0].Eta();
        bonsaiTree.j1t_phi=jet_comp[0].Phi();
        bonsaiTree.j1t_m=jet_comp[0].M();
        bonsaiTree.j1t_theta=jet_comp[0].Theta();
        bonsaiTree.j2t_E=jet_comp[1].E();
        bonsaiTree.j2t_pt=jet_comp[1].Pt();
        bonsaiTree.j2t_eta=jet_comp[1].Eta();
        bonsaiTree.j2t_phi=jet_comp[1].Phi();
        bonsaiTree.j2t_m=jet_comp[1].M();
        bonsaiTree.j2t_theta=jet_comp[1].Theta();
//
//      reco jets based on DR formula with fixed chi
//
//
        bonsaiTree.j1r_E=rec_comp[0].E();
        bonsaiTree.j1r_pt=rec_comp[0].pt();
        bonsaiTree.j1r_eta=rec_comp[0].eta();
        bonsaiTree.j1r_phi=rec_comp[0].phi();
        bonsaiTree.j1r_m=rec_comp[0].m();
        bonsaiTree.j1r_theta=rec_comp[0].theta();
        bonsaiTree.j2r_E=rec_comp[1].E();
        bonsaiTree.j2r_pt=rec_comp[1].pt();
        bonsaiTree.j2r_eta=rec_comp[1].eta();
        bonsaiTree.j2r_phi=rec_comp[1].phi();
        bonsaiTree.j2r_m=rec_comp[1].m();
        bonsaiTree.j2r_theta=rec_comp[1].theta();

//
        bonsaiTree.j1s_E=sci_comp[0].E();
        bonsaiTree.j1s_pt=sci_comp[0].pt();
        bonsaiTree.j1s_eta=sci_comp[0].eta();
        bonsaiTree.j1s_phi=sci_comp[0].phi();
        bonsaiTree.j1s_m=sci_comp[0].m();
        bonsaiTree.j1s_theta=sci_comp[0].theta();
        bonsaiTree.j2s_E=sci_comp[1].E();
        bonsaiTree.j2s_pt=sci_comp[1].pt();
        bonsaiTree.j2s_eta=sci_comp[1].eta();
        bonsaiTree.j2s_phi=sci_comp[1].phi();
        bonsaiTree.j2s_m=sci_comp[1].m();
        bonsaiTree.j2s_theta=sci_comp[1].theta();
        bonsaiTree.j1c_E=cher_comp[0].E();
        bonsaiTree.j1c_pt=cher_comp[0].pt();
        bonsaiTree.j1c_eta=cher_comp[0].eta();
        bonsaiTree.j1c_phi=cher_comp[0].phi();
        bonsaiTree.j1c_m=cher_comp[0].m();
        bonsaiTree.j1c_theta=cher_comp[0].theta();
        bonsaiTree.j2c_E=cher_comp[1].E();
        bonsaiTree.j2c_pt=cher_comp[1].pt();
        bonsaiTree.j2c_eta=cher_comp[1].eta();
        bonsaiTree.j2c_phi=cher_comp[1].phi();
        bonsaiTree.j2c_m=cher_comp[1].m();
        bonsaiTree.j2c_theta=cher_comp[1].theta();
  	}
    }// energy>0          
    
// fill output tree
    bonsaiTree.Fill();
  } // loop on events
  
  ftree->cd();
  bonsaiTree.Write();
  delete f;
  delete f1;

}
std::vector<double> calibscin(std::vector<double> vectorscin){
	std::vector<double> s_cont;
	s_cont = {408.21638950554075, 408.3954472740771, 407.1870232421094, 406.63875945884087, 404.8060585388971, 403.97304819147996, 403.3691105878475, 403.49367909804056, 404.55647780600043, 405.58591491094637, 403.9575182245898, 404.4757730162475, 404.72249522199195, 405.272159576985, 404.74332809708255, 404.83205898107536, 405.23195412471205, 404.9766105533868, 404.9085068798063, 404.9314555180952, 404.67532710488985, 404.58364980855805, 405.012793566413, 405.0007315500301, 404.30902206187204, 405.6974274788762, 405.2261341502687, 405.63975175649347, 404.90683641527, 404.37034541526305, 405.67260217215875, 405.5109490861691, 404.2898135363692, 405.07073526391474, 405.58981257625425, 405.3751447994642, 405.36549518339785, 405.3332161707569, 404.88956759976287, 405.37027184803094, 404.8980725551248, 405.34774082392767, 405.2984093045488, 405.14372480308344, 405.19187487160525, 405.03757034167137, 405.16280927227615, 404.7829216539207, 405.03107640207867, 404.7292557576276, 404.8025372723253, 403.9177916263665, 404.7460239584375, 403.96821450150077, 404.1905949169899, 404.1704924951662, 403.16496315846314, 402.2360298379118, 403.3863719919289, 402.9762332238292, 403.15699339382735, 403.4020052256797, 402.3032561236677, 402.8453577277423, 401.11356268338346, 401.3504783424065, 400.94087925309395, 400.29569405733, 400.0328154316862, 399.5130445431503, 398.66148407548866, 399.83880015591535, 398.96289406538807, 398.42261837089694, 391.76612693948175};
	
	int loop1 = s_cont.size();
	int loop2 = 36;

	for(int b=1; b<loop2; b++){
		for(int i=0; i<loop1; i++){
			s_cont.push_back(s_cont[i]);
		}
	}

	double c = 0.1;
	s_cont.insert(s_cont.begin(),c);

	if(s_cont.size() != vectorscin.size()){cout<<"ERROR in calibration!"<<endl;}

	std::vector<double> Calib_vectorScin;

	for(uint i=0; i<vectorscin.size(); i++){
		Calib_vectorScin.push_back(vectorscin[i]*(1.0/s_cont[i]));
	}

	return Calib_vectorScin;
}

std::vector<double> calibcher(std::vector<double> vectorcher){
	std::vector<double> c_cont;
        c_cont = {103.08779161895677, 102.91302749597065, 102.69865952763615, 102.61869191270468, 102.54928716539662, 102.48068194031679, 102.49984890080964, 102.35556540203991, 102.47969263317724, 102.6281510005559, 102.43322742473204, 102.47810836409134, 102.55371034296142, 102.67118096060427, 102.67297232291142, 102.48284061965019, 102.5649981010228, 102.56155933915096, 102.67809243921879, 102.56067521092992, 102.60224889784466, 102.63726587197354, 102.63191774143888, 102.76496337880408, 102.6929637252195, 102.60491403169074, 102.85913301772406, 102.741217657914, 102.69546934772463, 102.67035622618218, 102.69304228926421, 102.75886941001674, 102.75976221892324, 102.731492956408, 102.7188845221274, 102.77429845330465, 102.78649420797491, 102.75140309520445, 102.70051794706535, 102.68996042906552, 102.78365100098196, 102.8153738834064, 102.71292597825087, 102.73146416207084, 102.6450394621172, 102.61404003462839, 102.66675609739092, 102.60991640602225, 102.750246685674, 102.62575682868824, 102.42720794074478, 102.51305416968992, 102.52098979376447, 102.59751750679058, 102.45780037787654, 102.53083482963227, 102.47068539942974, 102.5721049950492, 102.56599170316093, 102.46469174495641, 102.19238017547394, 102.28148980648412, 102.19817435184497, 102.1330715125064, 102.09230341456059, 102.05765775486448, 101.9644426420847, 101.96014956820567, 101.85273676485993, 101.93311307596035, 101.96637882465569, 101.68716060542853, 101.55050000833062, 101.67603040894112, 99.77195006099979};

	int loop1 = c_cont.size();
	int loop2 = 36;

	for(int b=1; b<loop2; b++){
		for(int i=0; i<loop1; i++){
			c_cont.push_back(c_cont[i]);
		}
	}

	double c = 0.1;
	c_cont.insert(c_cont.begin(),c);

	if(c_cont.size() != vectorcher.size()){cout<<"ERROR in calibration!"<<endl;}

	std::vector<double> Calib_vectorCher;

	for(uint i=0; i<vectorcher.size(); i++){
		Calib_vectorCher.push_back(vectorcher[i]*(1.0/c_cont[i]));
	}

	return Calib_vectorCher;
}
//
std::tuple<double, double, double> maptower(int index, string side){
//Function to return tower angles (theta and phi) given index
  int NbOfBarrel=40;
  int NbOfEndcap=35;
  int NZrot=36;
  int TotTower=NbOfBarrel+NbOfEndcap;
  index = index-1;
  int sliceindex = index/TotTower;
  int towerindex = index-(sliceindex*TotTower);
  double deltatheta = 45./(NbOfBarrel);
//get theta
  double theta = towerindex*deltatheta+deltatheta/2.;
//  cout << " thetap " << theta << endl;
//get phi
  double phi_unit = 360./NZrot;
  double phi = (sliceindex)*phi_unit;
  
  if (side == "right"){
//     cout << " thetai " << theta+90. << " phii " << phi << " etai " <<  -log(tan(((90.-theta)*pi/180./2.))) << endl;
    return std::make_tuple(theta+90., phi, -log(tan(((90.-theta)*pi/180./2.))));
  }
  if (side == "left"){
    return std::make_tuple(90.-theta, phi, log(tan(((90.-theta)*pi/180./2.))));
  }
  return std::make_tuple(0.,0.,0.);
}
fastjet::PseudoJet mergejet(fastjet::PseudoJet jet_scin, fastjet::PseudoJet jet_cher, double c) {


  double jetPx = (jet_scin.px()-c*jet_cher.px())/(1-c);
  double jetPy = (jet_scin.py()-c*jet_cher.py())/(1-c);
  double jetPz = (jet_scin.pz()-c*jet_cher.pz())/(1-c);
  double jetE = (jet_scin.e()-c*jet_cher.e())/(1.-c);
  return fastjet::PseudoJet(jetPx, jetPy, jetPz, jetE);
}
fastjet::PseudoJet matchjet(fastjet::PseudoJet jet_in, vector<fastjet::PseudoJet> testvec) {

  int imin=-1;
  double deltarmin=99999.;
  for(uint i=0; i<testvec.size(); i++){
    double deltar=jet_in.delta_R(testvec.at(i));
    if(deltar<deltarmin) {
      deltarmin=deltar;
      imin=i;
    }
  }
  if(imin != -1) return testvec.at(imin);
  else
  return fastjet::PseudoJet(0., 0., 0., 0.);
}
double ptsmeatr(double pttr, double etatr) {
  double resotr(double, double);
  if(pttr>0.) {
    double sigma=resotr(pttr,etatr);
    double sigma1p=sigma/pttr/pttr;
    double rpttr=1./pttr;
    while(1) {
      double rptsmear=rpttr+m_random3.Gaus(0,sigma1p);
      if(rptsmear>0.0000001) {
        double ptsmear=1./rptsmear;
        return ptsmear;
//	cout << " pt " << pttr << " ptsmear " << ptsmear << endl;
      }
    }
  }
  else return 0.;
}
// DELPHES card
//sqrt(0.0001145^2 + 0.0002024^2*pt + (pt*2.093e-005)^2)
double resotr(double pttr, double etatr) {
  double sigma=0.;
  double etalow=0.;
  double etahigh=3.;
  double a1=0.0001145;
  double a2=0.0002024;
  double a3=2.093e-005;
  if(fabs(etatr)>=etalow && fabs(etatr)<etahigh) {
      double sigmaid=pttr*sqrt(a1*a1+a2*a2*pttr+a3*a3*pttr*pttr);
      sigma=sigmaid;
  }
  return sigma;
}
//
//  calculate jet energies in 2-jet events from beam energy constraint
//  and reconstructed jet masses and angles
//
std::pair<double,double>ejconst(double k, double th1, double th2, double mass1, double mass2) {
   double e1=0.;
   double e2=0.;
   double s1=sin(th1)*sin(th1);
   double s2=sin(th2)*sin(th2);
   double m1=mass1*mass1;
   double m2=mass2*mass2;
   double a=s1-s2;
   double b=2*k*s2;
   double c=m2*s2-m1*s1-k*k*s2;
   double disc=b*b-4*a*c;
   if(disc>=0){
     double e1a=(-b+sqrt(disc))/2/a;
     double e1b=(-b-sqrt(disc))/2/a;
     if(abs(e1a)<abs(e1b) && e1a>0 && e1a<k){
       e1=e1a;
       e2=k-e1;
     }
   }
   return pair<double,double>(e1,e2);
}


