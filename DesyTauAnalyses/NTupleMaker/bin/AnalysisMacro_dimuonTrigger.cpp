#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"

#include "TLorentzVector.h"

#include "TRandom.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

#include "DesyTauAnalyses/NTupleMaker/interface/RunLumiReader.h"


const float MuMass = 0.105658367;

int binNumber(float x, int nbins, float * bins) {

  int binN = 0;

  for (int iB=0; iB<nbins; ++iB) {
    if (x>=bins[iB]&&x<bins[iB+1]) {
      binN = iB;
      break;
    }
  }

  return binN;

}

float effBin(float x, int nbins, float * bins, float * eff) {

  int bin = binNumber(x, nbins, bins);

  return eff[bin];

}

double cosRestFrame(TLorentzVector boost, TLorentzVector vect) {

  double bx = -boost.Px()/boost.E();
  double by = -boost.Py()/boost.E();
  double bz = -boost.Pz()/boost.E();

  vect.Boost(bx,by,bz);
  double prod = -vect.Px()*bx-vect.Py()*by-vect.Pz()*bz;
  double modBeta = TMath::Sqrt(bx*bx+by*by+bz*bz); 
  double modVect = TMath::Sqrt(vect.Px()*vect.Px()+vect.Py()*vect.Py()+vect.Pz()*vect.Pz());
  
  double cosinus = prod/(modBeta*modVect);

  return cosinus;

}

double QToEta(double Q) {
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;
}

double EtaToQ(double Eta) {
  double Q = 2.0*TMath::ATan(TMath::Exp(-Eta));
  if (Q<0.0) Q += TMath::Pi();
  return Q;
}

double PtoEta(double Px, double Py, double Pz) {

  double P = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
  double cosQ = Pz/P;
  double Q = TMath::ACos(cosQ);
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;

}

double PtoPhi(double Px, double Py) {
  return TMath::ATan2(Py,Px);
}

double PtoPt(double Px, double Py) {
  return TMath::Sqrt(Px*Px+Py*Py);
}

double dPhiFrom2P(double Px1, double Py1,
		  double Px2, double Py2) {


  double prod = Px1*Px2 + Py1*Py2;
  double mod1 = TMath::Sqrt(Px1*Px1+Py1*Py1);
  double mod2 = TMath::Sqrt(Px2*Px2+Py2*Py2);
  
  double cosDPhi = prod/(mod1*mod2);
  
  return TMath::ACos(cosDPhi);

}

double deltaEta(double Px1, double Py1, double Pz1,
		double Px2, double Py2, double Pz2) {

  double eta1 = PtoEta(Px1,Py1,Pz1);
  double eta2 = PtoEta(Px2,Py2,Pz2);

  double dEta = eta1 - eta2;

  return dEta;

}

double deltaR(double Eta1, double Phi1,
	      double Eta2, double Phi2) {

  double Px1 = TMath::Cos(Phi1);
  double Py1 = TMath::Sin(Phi1);

  double Px2 = TMath::Cos(Phi2);
  double Py2 = TMath::Sin(Phi2);

  double dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2);
  double dEta = Eta1 - Eta2;

  double dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);

  return dR;

}

double PtEtaToP(double Pt, double Eta) {

  //  double Q = EtaToQ(Eta);

  //double P = Pt/TMath::Sin(Q);
  double P = Pt*TMath::CosH(Eta);

  return P;
}
double Px(double Pt, double Phi){

  double Px=Pt*TMath::Cos(Phi);
  return Px;
}
double Py(double Pt, double Phi){

  double Py=Pt*TMath::Sin(Phi);
  return Py;
}
double Pz(double Pt, double Eta){

  double Pz=Pt*TMath::SinH(Eta);
  return Pz;
}
double InvariantMass(double energy,double Px,double Py, double Pz){

  double M_2=energy*energy-Px*Px-Py*Py-Pz*Pz;
  double M=TMath::Sqrt(M_2);
  return M;


}
double EFromPandM0(double M0,double Pt,double Eta){

  double E_2=M0*M0+PtEtaToP(Pt,Eta)*PtEtaToP(Pt,Eta);
  double E =TMath::Sqrt(E_2);
  return E;

}

bool electronMvaIdTight(float eta, float mva) {

  float absEta = fabs(eta);

  bool passed = false;
  if (absEta<0.8) {
    if (mva>0.73) passed = true;
  }
  else if (absEta<1.479) {
    if (mva>0.57) passed = true;
  }
  else {
    if (mva>0.05) passed = true;
  }

  return passed;

}

bool electronMvaIdLoose(float eta, float mva) {

  float absEta = fabs(eta);

  bool passed = false;
  if (absEta<0.8) {
    if (mva>0.35) passed = true;
  }
  else if (absEta<1.479) {
    if (mva>0.20) passed = true;
  }
  else {
    if (mva>-0.52) passed = true;
  }

  return passed;

}

bool electronMvaIdWP80(float pt, float eta, float mva) {

  float absEta = fabs(eta);
  bool passed = false;
  if (absEta<0.8) {
    if (pt<10) 
      passed = mva > -0.253;
    else 
      passed = mva > 0.965;
  }
  else if (absEta<1.479) {
    if (pt<10)
      passed = mva > 0.081;
    else
      passed = mva > 0.917;
  }
  else {
    if (pt<10)
      passed = mva > -0.081;
    else
      passed = mva > 0.683;
  }

  return passed;

}

bool electronMvaIdWP90(float pt, float eta, float mva) {

  float absEta = fabs(eta);
  bool passed = false;
  if (absEta<0.8) {
    if (pt<10) 
      passed = mva > -0.483;
    else 
      passed = mva > 0.933;
  }
  else if (absEta<1.479) {
    if (pt<10)
      passed = mva > -0.267;
    else
      passed = mva > 0.825;
  }
  else {
    if (pt<10)
      passed = mva > -0.323;
    else
      passed = mva > 0.337;
  }

  return passed;

}

struct myclass {
  bool operator() (int i,int j) { return (i<j);}
} myobject;

const float electronMass = 0;
const float muonMass = 0.10565837;
const float pionMass = 0.1396;

int main(int argc, char * argv[]) {

  // first argument - config file 
  // second argument - filelist

  using namespace std;

  // **** configuration
  Config cfg(argv[1]);

  // Data
  const bool isData = cfg.get<bool>("IsData");
  const bool applyGoodRunSelection = cfg.get<bool>("ApplyGoodRunSelection"); 
  const string jsonFile = cfg.get<string>("jsonFile");

  // kinematic cuts on muons
  const float ptMuonLowCut   = cfg.get<float>("ptMuonLowCut");
  const float ptMuonHighCut  = cfg.get<float>("ptMuonHighCut");
  const float etaMuonLowCut  = cfg.get<float>("etaMuonLowCut");
  const float etaMuonHighCut = cfg.get<float>("etaMuonHighCut");
  const float dxyMuonCut     = cfg.get<float>("dxyMuonCut");
  const float dzMuonCut      = cfg.get<float>("dzMuonCut");
  const float isoMuonCut     = cfg.get<float>("isoMuonCut");
  const bool applyMuonId     = cfg.get<bool>("ApplyMuonId");
  const bool applyMuonIso    = cfg.get<bool>("ApplyMuonIso");

  // topological cuts
  const float dRleptonsCut   = cfg.get<float>("dRleptonsCut");
  const float dZleptonsCut   = cfg.get<float>("dZleptonsCut");
  const float DRTrigMatch    = cfg.get<float>("DRTrigMatch");
  const bool oppositeSign    = cfg.get<bool>("OppositeSign");

  // triggers
  const bool applyTrigger = cfg.get<bool>("ApplyTrigger");
  const string hltMu17Mu8           = cfg.get<string>("HLTMu17Mu8");
  const string hltMu17Mu8DZ         = cfg.get<string>("HLTMu17Mu8DZ");
  const string hltMu17Mu8SameSignDZ = cfg.get<string>("HLTMu17Mu8SameSignDZ");

  // HLT filters
  const string hltMu17Leg     = cfg.get<string>("HLTMu17Leg");
  const string hltMu8Leg      = cfg.get<string>("HLTMu8Leg");
  const string dzFilter       = cfg.get<string>("dzFilter"); 
  const string sameSignFilter = cfg.get<string>("sameSignFilter");

  TString HLTMu17Mu8(hltMu17Mu8);
  TString HLTMu17Mu8DZ(hltMu17Mu8DZ);
  TString HLTMu17Mu8SameSignDZ(hltMu17Mu8SameSignDZ);
  TString HLTMu17Leg(hltMu17Leg);
  TString HLTMu8Leg(hltMu8Leg);
  TString DZFilter(dzFilter);
  TString SameSignFilter(sameSignFilter);
  
  // vertex cuts
  const float ndofVertexCut  = cfg.get<float>("NdofVertexCut");   
  const float zVertexCut     = cfg.get<float>("ZVertexCut");
  const float dVertexCut     = cfg.get<float>("DVertexCut");


  // **** end of configuration

  string cmsswBase = (getenv ("CMSSW_BASE"));

  // Run-lumi selector
  std::vector<std::string> jsonFiles;
  jsonFiles.push_back(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/test/"+jsonFile);

  RunLumiSelector runLumiSelector;
  runLumiSelector = RunLumiSelector(jsonFiles);

  // file name and tree name
  std::string rootFileName(argv[2]);
  std::ifstream fileList(argv[2]);
  std::ifstream fileList0(argv[2]);
  std::string ntupleName("makeroottree/AC1B");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  

  // output fileName with histograms
  TFile * file = new TFile(TStrName+TString(".root"),"recreate");
  file->cd("");
  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);
  TH1F * weightsH = new TH1F("weightsH","",1,-0.5,0.5);


  // J/Psi ->
  TH1F * JPsiMassDZFilterPassH =  new TH1F("JPsiMassDZFilterPassH","",200,2,4);
  TH1F * JPsiMassDZFilterFailH =  new TH1F("JPsiMassDZFilterFailH","",200,2,4);

  TH1F * JPsiMassDZFilterDz0to1PassH =  new TH1F("JPsiMassDZFilterDz0to1PassH","",200,2,4);
  TH1F * JPsiMassDZFilterDz0to1FailH =  new TH1F("JPsiMassDZFilterDz0to1FailH","",200,2,4);

  TH1F * JPsiMassDZFilterDz1to2PassH =  new TH1F("JPsiMassDZFilterDz1to2PassH","",200,2,4);
  TH1F * JPsiMassDZFilterDz1to2FailH =  new TH1F("JPsiMassDZFilterDz1to2FailH","",200,2,4);

  TH1F * JPsiMassDZFilterDzGt2PassH =  new TH1F("JPsiMassDZFilterDzGt2PassH","",200,2,4);
  TH1F * JPsiMassDZFilterDzGt2FailH =  new TH1F("JPsiMassDZFilterDzGt2FailH","",200,2,4);

  TH1F * JPsiMassSameSignFilterPassH =  new TH1F("JPsiMassSameSignFilterPassH","",200,2,4);
  TH1F * JPsiMassSameSignFilterFailH =  new TH1F("JPsiMassSameSignFilterFailH","",200,2,4);

  // Z ->
  TH1F * ZMassDZFilterPassH =  new TH1F("ZMassDZFilterPassH","",60,60,120);
  TH1F * ZMassDZFilterFailH =  new TH1F("ZMassDZFilterFailH","",60,60,120);

  TH1F * ZMassDZFilterDz0to1PassH =  new TH1F("ZMassDZFilterDz0to1PassH","",60,60,120);
  TH1F * ZMassDZFilterDz0to1FailH =  new TH1F("ZMassDZFilterDz0to1FailH","",60,60,120);

  TH1F * ZMassDZFilterDz1to2PassH =  new TH1F("ZMassDZFilterDz1to2PassH","",60,60,120);
  TH1F * ZMassDZFilterDz1to2FailH =  new TH1F("ZMassDZFilterDz1to2FailH","",60,60,120);

  TH1F * ZMassDZFilterDzGt2PassH =  new TH1F("ZMassDZFilterDzGt2PassH","",60,60,120);
  TH1F * ZMassDZFilterDzGt2FailH =  new TH1F("ZMassDZFilterDzGt2FailH","",60,60,120);

  TH1F * ZMassSameSignFilterPassH =  new TH1F("ZMassSameSignFilterPassH","",60,60,120);
  TH1F * ZMassSameSignFilterFailH =  new TH1F("ZMassSameSignFilterFailH","",60,60,120);


  unsigned int iRun;
  unsigned int iEvent;
  TTree * eventTree = new TTree("eventTree","eventTree");
  eventTree->Branch("Run",&iRun,"Run/i");
  eventTree->Branch("Event",&iEvent,"Event/i");


  int nFiles = 0;
  int nEvents = 0;

  int selEvents = 0;
  int selEventsHLTMu17Mu8 = 0;
  int selEventsHLTMu17Mu8DZ = 0;
  int selEventsHLTMu17Mu8SameSignDZ = 0;

  int selPairs = 0;
  int selPairsHLTMu17Mu8 = 0;
  int selPairsHLTMu17Mu8DZ = 0;
  int selPairsHLTMu17Mu8SameSignDZ = 0;

  unsigned int minRun = 99999999;
  unsigned int maxRun = 0;

  std::vector<unsigned int> allRuns; allRuns.clear();

  int nTotalFiles = 0;
  std::string dummy;
  // count number of files --->
  while (fileList0 >> dummy) nTotalFiles++;

  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    std::cout << "file " << iF+1 << " out of " << nTotalFiles << " filename : " << filen << std::endl;
    TFile * file_ = TFile::Open(TString(filen));
    
    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));
  
    if (_tree==NULL) continue;
    
    TH1D * histoInputEvents = NULL;
   
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    
    if (histoInputEvents==NULL) continue;
    
    int NE = int(histoInputEvents->GetEntries());
    
    std::cout << "      number of input events    = " << NE << std::endl;
    
    for (int iE=0;iE<NE;++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree);
    
    Long64_t numberOfEntries = analysisTree.GetEntries();
    
    std::cout << "      number of entries in Tree = " << numberOfEntries << std::endl;
    
    for (Long64_t iEntry=0; iEntry<numberOfEntries; iEntry++) { 
    
      analysisTree.GetEntry(iEntry);
      nEvents++;
      
      if (nEvents%10000==0) 
	cout << "      processed " << nEvents << " events" << endl; 

      float weight = 1;

      if (analysisTree.event_run>maxRun)
	maxRun = analysisTree.event_run;

      if (analysisTree.event_run<minRun)
	minRun = analysisTree.event_run;


      bool isNewRun = true;
      if (allRuns.size()>0) {
	for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
	  if (analysisTree.event_run==allRuns.at(iR)) {
	    isNewRun = false;
	    break;
	  }
	}
      }

      if (isNewRun) 
	allRuns.push_back(analysisTree.event_run);

      if (isData) {
	if (applyGoodRunSelection && !runLumiSelector.accept(analysisTree.event_run, analysisTree.event_luminosityblock))
	  continue;
      }

      weightsH->Fill(0.0,weight);

      // triggers
      bool isHLTMu17Mu8 = false;
      bool isHLTMu17Mu8DZ = false;
      bool isHLTMu17Mu8SameSignDZ = false;
      for (std::map<string,int>::iterator it=analysisTree.hltriggerresults->begin(); it!=analysisTree.hltriggerresults->end(); ++it) {
	TString trigName(it->first);
	if (trigName.Contains(HLTMu17Mu8)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  if (it->second==1)
	    isHLTMu17Mu8 = true;
	}
	if (trigName.Contains(HLTMu17Mu8DZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  if (it->second==1)
            isHLTMu17Mu8DZ = true;
	}
	if (trigName.Contains(HLTMu17Mu8SameSignDZ)) {
	  //	  std::cout << trigName << " : " << it->second << std::endl;
	  if (it->second==1)
            isHLTMu17Mu8SameSignDZ = true;
	}
      }

      unsigned int nMu17Leg = 0;
      bool isMu17Leg = false;
      unsigned int nMu8Leg = 0;
      bool isMu8Leg = false;
      unsigned int nDZFilter = 0;
      bool isDZFilter = false;
      unsigned int nSameSignFilter = 0;
      bool isSameSignFilter = false;
      unsigned int nfilters = analysisTree.run_hltfilters->size();
      for (unsigned int i=0; i<nfilters; ++i) {
        TString HLTFilter(analysisTree.run_hltfilters->at(i));
	if (HLTFilter==HLTMu17Leg) {
	  nMu17Leg = i;
	  isMu17Leg = true;
	}
	if (HLTFilter==HLTMu8Leg) {
	  nMu8Leg = i;
	  isMu8Leg = true;
	}
	if (HLTFilter==DZFilter) {
	  nDZFilter = i;
	  isDZFilter = true;
	}
	if (HLTFilter==SameSignFilter) {
	  nSameSignFilter = i;
	  isSameSignFilter = true;
	}
      }
      if (!isMu17Leg) {
	std::cout << "HLT filter " << HLTMu17Leg << " not found" << std::endl;
	exit(-1);
      }
      if (!isMu8Leg) {
	std::cout << "HLT filter " << HLTMu8Leg << " not found" << std::endl;
	exit(-1);
      }
      if (!isDZFilter) {
	std::cout << "HLT filter " << DZFilter << " not found" << std::endl;
	exit(-1);
      }
      if (!isSameSignFilter) {
	std::cout << "HLT filter " << SameSignFilter << " not found" << std::endl;
	exit(-1);
      }

      // vertex cuts
      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      float dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
		       analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;

      
      // muon selection
      vector<int> muons; muons.clear();
      for (unsigned int im = 0; im<analysisTree.muon_count; ++im) {
	if (analysisTree.muon_pt[im]<ptMuonLowCut) continue;
	if (fabs(analysisTree.muon_eta[im])>etaMuonLowCut) continue;
	if (fabs(analysisTree.muon_dxy[im])>dxyMuonCut) continue;
	if (fabs(analysisTree.muon_dz[im])>dzMuonCut) continue;
	if (applyMuonId && !analysisTree.muon_isMedium[im]) continue;
	float absIso = analysisTree.muon_chargedHadIso[im];
	float relIso = absIso/analysisTree.muon_pt[im];
	if (relIso>isoMuonCut&&applyMuonIso) continue;
	muons.push_back(im);
      }

      if (muons.size()<2) continue;

      bool isPairSelected = false;
      bool isPairSelectedHLTMu17Mu8 = false;
      bool isPairSelectedHLTMu17Mu8DZ = false;
      bool isPairSelectedHLTMu17Mu8SameSignDZ = false;

      // selecting muon pair
      for (unsigned int im1=0; im1<muons.size()-1; ++im1) {
	//	  std::cout << "Muon " << im << std::endl;
	int  mu1Index = muons[im1];
	bool mu1MatchMu17 = false;
	bool mu1MatchMu8  = false;
	bool mu1MatchDz   = false;
	bool mu1MatchSS   = false;
	for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  float dRtrig = deltaR(analysisTree.muon_eta[mu1Index],analysisTree.muon_phi[mu1Index],
				analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	  if (dRtrig>DRTrigMatch) continue;
	  if (analysisTree.trigobject_filters[iT][nMu17Leg]) // Muon17 Leg
	      mu1MatchMu17 = true;
	  if (analysisTree.trigobject_filters[iT][nMu8Leg]) // Muon8 Leg
	      mu1MatchMu8 = true;
	  if (analysisTree.trigobject_filters[iT][nDZFilter]) // DZ filter
	      mu1MatchDz = true;
	  if (analysisTree.trigobject_filters[iT][nSameSignFilter]) // same-sign filter
	      mu1MatchSS = true;
	   
	}
	bool mu1Mu17 = mu1MatchMu17 && analysisTree.muon_pt[mu1Index]>ptMuonHighCut && fabs(analysisTree.muon_eta[mu1Index])<etaMuonHighCut;
	bool mu1Mu8  = mu1MatchMu8 && analysisTree.muon_pt[mu1Index]>ptMuonLowCut && fabs(analysisTree.muon_eta[mu1Index])<etaMuonLowCut;

	float q1 = analysisTree.muon_charge[mu1Index];

	for (unsigned int im2=im1+1; im2<muons.size(); ++im2) {

	  int  mu2Index = muons[im2];

	  float q2 = analysisTree.muon_charge[mu2Index];
	  if (oppositeSign && (q1*q2>0)) continue;

	  float dRmumu = deltaR(analysisTree.muon_eta[mu1Index],analysisTree.muon_phi[mu1Index],
				analysisTree.muon_eta[mu2Index],analysisTree.muon_phi[mu2Index]);
	  if (dRmumu>dRleptonsCut) continue;

	  bool mu2MatchMu17 = false;
	  bool mu2MatchMu8  = false;
	  bool mu2MatchDz   = false;
	  bool mu2MatchSS   = false;
	  for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	    float dRtrig = deltaR(analysisTree.muon_eta[mu2Index],analysisTree.muon_phi[mu2Index],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig>DRTrigMatch) continue;
	    if (analysisTree.trigobject_filters[iT][nMu17Leg]) // Muon17 Leg
              mu2MatchMu17 = true;
	    if (analysisTree.trigobject_filters[iT][nMu8Leg]) // Muon8 Leg
              mu2MatchMu8 = true;
	    if (analysisTree.trigobject_filters[iT][nDZFilter]) // DZ filter
              mu2MatchDz = true;
	    if (analysisTree.trigobject_filters[iT][nSameSignFilter]) // same-sign filter
              mu2MatchSS = true;

	  }
	  bool mu2Mu17 = mu2MatchMu17 && analysisTree.muon_pt[mu2Index]>ptMuonHighCut && fabs(analysisTree.muon_eta[mu2Index])<etaMuonHighCut;
	  bool mu2Mu8  = mu2MatchMu8 && analysisTree.muon_pt[mu2Index]>ptMuonLowCut && fabs(analysisTree.muon_eta[mu2Index])<etaMuonLowCut;

	  bool triggerMatch = (mu1Mu17&&mu2Mu8) || (mu1Mu8&&mu2Mu17);
	  bool triggerMatchDz = triggerMatch && mu1MatchDz && mu2MatchDz;
	  bool triggerMatchSS = triggerMatchDz && mu1MatchSS && mu2MatchSS;

	  float dZ = fabs(analysisTree.muon_dz[mu1Index]-analysisTree.muon_dz[mu2Index]);

	  TLorentzVector mu1lv; mu1lv.SetXYZM(analysisTree.muon_px[mu1Index],
					      analysisTree.muon_py[mu1Index],
					      analysisTree.muon_pz[mu1Index],
					      muonMass);
	  TLorentzVector mu2lv; mu2lv.SetXYZM(analysisTree.muon_px[mu2Index],
                                              analysisTree.muon_py[mu2Index],
                                              analysisTree.muon_pz[mu2Index],
                                              muonMass);

	  float mass = (mu1lv+mu2lv).M();

	  isPairSelected = true;
	  selPairs++;
	  
	  if (isHLTMu17Mu8 && triggerMatch) { // pass HLT_Mu17_Mu8
	    isPairSelectedHLTMu17Mu8 = true;
	    selPairsHLTMu17Mu8++;
	    if (triggerMatchDz) { // pass HLT_Mu17_Mu8_DZ
	      if (dZ<dZleptonsCut) {
		JPsiMassDZFilterPassH->Fill(mass,weight);
		ZMassDZFilterPassH->Fill(mass,weight);
	      }
	      if (dZ<0.1) {
		JPsiMassDZFilterDz0to1PassH->Fill(mass,weight);
		ZMassDZFilterDz0to1PassH->Fill(mass,weight);
	      }
	      else if (dZ<0.2) {
		JPsiMassDZFilterDz1to2PassH->Fill(mass,weight); 
                ZMassDZFilterDz1to2PassH->Fill(mass,weight);
	      }
	      else {
		JPsiMassDZFilterDzGt2PassH->Fill(mass,weight);
                ZMassDZFilterDzGt2PassH->Fill(mass,weight);
	      }
	    }
	    else { // fail HLT_Mu17_Mu8_DZ
	      if (dZ<dZleptonsCut) {
		JPsiMassDZFilterFailH->Fill(mass,weight);
                ZMassDZFilterFailH->Fill(mass,weight);
              }
	      if (dZ<0.1) {
                JPsiMassDZFilterDz0to1FailH->Fill(mass,weight);
                ZMassDZFilterDz0to1FailH->Fill(mass,weight);
              }
              else if (dZ<0.2) {
                JPsiMassDZFilterDz1to2FailH->Fill(mass,weight);
                ZMassDZFilterDz1to2FailH->Fill(mass,weight);
              }
              else {
                JPsiMassDZFilterDzGt2FailH->Fill(mass,weight);
                ZMassDZFilterDzGt2FailH->Fill(mass,weight);
              }
	    }
	  }

	  if (isHLTMu17Mu8DZ && triggerMatchDz) { // pass HLT_Mu17_Mu8_DZ
	    isPairSelectedHLTMu17Mu8DZ = true;
            selPairsHLTMu17Mu8DZ++;
	    if (triggerMatchSS) { // pass HLT_Mu17_Mu8_SameSign_DZ
              if (dZ<dZleptonsCut) {
                JPsiMassDZFilterPassH->Fill(mass,weight);
                ZMassDZFilterPassH->Fill(mass,weight);
              }
            }
            else { // fail HLT_Mu17_Mu8_SameSign_DZ
              if (dZ<dZleptonsCut) {
                JPsiMassDZFilterFailH->Fill(mass,weight);
                ZMassDZFilterFailH->Fill(mass,weight);
              }
            }
	  }

	  if (isHLTMu17Mu8SameSignDZ && triggerMatchSS) {
	    isPairSelectedHLTMu17Mu8SameSignDZ = true;
            selPairsHLTMu17Mu8SameSignDZ++;
	  }

	}
      }
    
      if (isPairSelected)
	selEvents++;
      if (isPairSelectedHLTMu17Mu8)
	selEventsHLTMu17Mu8++;
      if (isPairSelectedHLTMu17Mu8DZ)
	selEventsHLTMu17Mu8DZ++;
      if (isPairSelectedHLTMu17Mu8SameSignDZ)
	selEventsHLTMu17Mu8SameSignDZ++;

      
    } // end of file processing (loop over events in one file)
    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }
  std::cout << std::endl;
  int allEvents = int(inputEventsH->GetEntries());
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << "Total number of selected pairs  = " << selPairs << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events (HLT_Mu17_Mu8) = " << selEventsHLTMu17Mu8 << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Mu17_Mu8) = " << selPairsHLTMu17Mu8  << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events (HLT_Mu17_Mu8_DZ) = " << selEventsHLTMu17Mu8DZ << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Mu17_Mu8_DZ) = " << selPairsHLTMu17Mu8DZ << std::endl;
  std::cout << std::endl;
  std::cout << "Total number of selected events (HLT_Mu17_Mu8_SameSign_DZ) = " << selEventsHLTMu17Mu8SameSignDZ << std::endl;
  std::cout << "Total number of selected pairs  (HLT_Mu17_Mu8_SameSign_DZ) = " << selPairsHLTMu17Mu8SameSignDZ << std::endl;
  std::cout << std::endl;
  std::cout << "Run range " << minRun << ":" << maxRun << std::endl;
  std::cout << std::endl;
  // using object as comp
  std::sort (allRuns.begin(), allRuns.end(), myobject);
  std::cout << "Runs : ";
  for (unsigned int iR=0; iR<allRuns.size(); ++iR)
    std::cout << " " << allRuns.at(iR);
  std::cout << std::endl;
  file->cd("");
  file->Write();
  file->Close();
  delete file;
  
}



