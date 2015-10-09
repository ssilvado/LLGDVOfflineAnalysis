#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <stdlib.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TColor.h>
#include <math.h>
#include <THnSparse.h>
#include <TChain.h>
#include <map>
#include <string>
#include <vector>
#include <TRandom3.h>
#include "THnSparse.h"
#include "TF1.h"
#include "TSystem.h"

using namespace std;

class LLGAnalysis {
    public:
        static LLGAnalysis* GetInstance( char* configFileName );
        ~LLGAnalysis() {}
        
        vector<double> CalculateVertex( vector<double> x, vector<double> y, vector<double> z, vector<double> weight, vector<int> charge, vector<double> distance, unsigned int &nConsidered, double &weightednConsidered, vector<double> &error ); 
        
        bool Init();
        void RunEventLoop( int nEventsMax = -1);
        void FinishRun();
    
    private:
        static LLGAnalysis* _instance;

        void RunObjectID();

        void makeHist( string nametitle, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, string xtitle, string ytitle, string ztitle, string drawOption = "", double xAxisOffset = 1., double yAxisOffset = 1.2, double zAxisOffset = 1. ); 
        void makeHist( string nametitle, int nbins, double xmin, double xmax, string xtitle, string ytitle, string drawOption = "", double xAxisOffset = 1., double yAxisOffset = 1.2 );
        void setStyle(double ytoff = 1.0, bool marker = true, double left_margin = 0.15); 
        void MakeEfficiencyPlot( TH1D hpass, TH1D htotal, TCanvas *c, string triggerName = "");
        void FillEfficiencyHistograms();
        
        void SetupSignalRegion();
        void SignalRegionSelection();
        void SetupTTJetsCR();
        void TTJetsCRSelection();
        // INSERT YOUR SELECTION HERE


    private:
        LLGAnalysis() {}
        LLGAnalysis( char* configFileName );
        
        map<string, int>        _cutFlow;
        map<string, TH1D>       _histograms1D;
        map<string, TH2D>       _histograms2D;
        map<string, string>     _histograms1DDrawOptions;
        map<string, string>     _histograms2DDrawOptions;
        vector<string>          _inputFileNames;
        string                  _inputTreeName;
        TChain                  *_inputTree;
        string                  _outputFileName;
        string 					outputHistos;
        TTree                   *_outputTree;
        TFile                   *_outputFile;
        bool                    _writeOutputTree;

        vector<double> *recoJet_pt; 
        vector<double> *recoJet_phi; 
        vector<double> *recoJet_eta; 
        vector<double> *recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags;
        vector<double> *recoJet_btag_jetBProbabilityBJetTags;
        vector<double> *recoJet_btag_jetProbabilityBJetTags;
        vector<double> *recoJet_btag_trackCountingHighPurBJetTags;
        vector<double> *recoJet_btag_trackCountingHighEffBJetTags;
        vector<double> *recoJet_vertex_x;
        vector<double> *recoJet_vertex_y;
        vector<double> *recoJet_vertex_z;
        vector<double> *muon_px; 
        vector<double> *muon_py; 
        vector<double> *muon_pz; 
        vector<double> *muon_phi; 
        vector<double> *muon_eta; 
        vector<double> *muon_iso;
        vector<bool>   *muon_isLooseMuon;
        vector<bool>   *muon_isTightMuon;
        vector<double> *electron_px; 
        vector<double> *electron_py; 
        vector<double> *electron_pz; 
        vector<double> *electron_phi; 
        vector<double> *electron_eta; 
        vector<double> *electron_iso;
        vector<bool> *electron_isVeto;
        vector<bool> *electron_isLoose;
        vector<bool> *electron_isMedium;
        vector<bool> *electron_isTight;
        vector<bool> *electron_isHEEP;
        vector<int> *triggerBits; 
        vector<string> *triggerNames; 
        /*
        vector<vector<double> > *recoJet_constVertex_x; 
        vector<vector<double> > *recoJet_constVertex_y; 
        vector<vector<double> > *recoJet_constVertex_z; 
        vector<vector<double> > *recoJet_const_pt; 
        vector<vector<double> > *recoJet_const_closestVertex_dxy; 
        vector<vector<double> > *recoJet_const_closestVertex_dz; 
        vector<vector<double> > *recoJet_const_closestVertex_d;
        vector<vector<int> > *recoJet_const_charge; 
        */
        vector<bool> *recoJet_isLeptonLike;

        vector<double> *vertex_x;
        vector<double> *vertex_y;
        vector<double> *vertex_z; 
        vector<double> *vertex_nTracks; 
        vector<double> *vertex_pt; 
        vector<double> *vertex_ndof;
        vector<double> *vertex_dx;
        vector<double> *vertex_dy;
        vector<double> *vertex_dz;

        vector<double> *secVertex_x;
        vector<double> *secVertex_y; 
        vector<double> *secVertex_z;
        vector<double> *secVertex_ndof;
        vector<double> *secVertex_pt;
        vector<double> *secVertex_dx;
        vector<double> *secVertex_dy; 
        vector<double> *secVertex_dz; 

        
        vector<int> vetoElectrons;
        vector<int> looseElectrons;
        vector<int> mediumElectrons;
        vector<int> tightElectrons;
        vector<int> heepElectrons;
        vector<int> vetoMuons;
        vector<int> selectedJets;

        double met;
        double met_x;
        double met_y;


        double evtWeight;
        double JET_PT_CUT_SV;
        double JET_PT_CUT_PV;
        double JET_ETA_CUT;
        double MUON_PT_CUT;
        double ELECTRON_PT_CUT;
        double MET_CUT;
        double LEADING_SV_JET_CUT;
        double MJJ_CUT;
        std::string SELECTION;
        std::string metadataFileName;
        std::string datasetName;
        std::vector<std::vector<double> > _yields2DOptimisation;

        // the total number of events per sample (take it from DAS!)
        double PROC_NTOT;
        // the total cross section for the process in pb
        double PROC_XSEC;
        // the target luminosity in fb-1
        double TARGET_LUMI;

        bool applyEventWeights;
        
    
        vector<string> _plotFormats;

};

