#include "LLGAnalysis.h"
#include "TGraphAsymmErrors.h"
#include <fstream>
#include <iomanip>

LLGAnalysis* LLGAnalysis::GetInstance( char *configFileName ) {
    if( !_instance ) {
        _instance = new LLGAnalysis( configFileName );
    }
    return _instance;
}

LLGAnalysis::LLGAnalysis( char *configFileName ) {
    
    
    // Setup the default values for the cuts:
    JET_PT_CUT_SV = 30;
    JET_PT_CUT_PV = 75;
    JET_ETA_CUT = 5.0;
    LEADING_SV_JET_CUT = 10000.;
    MJJ_CUT = 10000;
    MUON_PT_CUT = 15.;
    ELECTRON_PT_CUT = 15.;
    MET_CUT = 210.;
    PROC_XSEC = 1.;
    PROC_NTOT = 1.;
    applyEventWeights = false;
    TARGET_LUMI = 1000.;
    //SELECTION = "SignalRegion";
    SELECTION = "TTJetsCR";
    metadataFileName = "Configuration/DatasetMetadata.txt";
    datasetName = "Signal_500_60";
    _writeOutputTree = true;

    ifstream configFile( configFileName, ios::in );
    while( configFile.good() ) {
        string key, value;
        configFile >> key >> ws >> value;
        if( configFile.eof() ) break;
        if( key == "InputFile"         ) _inputFileNames.push_back( value ); 
        if( key == "InputTree"         ) _inputTreeName = value; 
        if( key == "JET_PT_CUT_SV"     ) JET_PT_CUT_SV = atof(value.c_str());
        if( key == "JET_PT_CUT_PV"     ) JET_PT_CUT_PV = atof(value.c_str());
        if( key == "JET_ETA_CUT"       ) JET_ETA_CUT = atof(value.c_str());
        if( key == "MUON_PT_CUT"       ) MUON_PT_CUT = atof(value.c_str()); 
        if( key == "ELECTRON_PT_CUT"   ) ELECTRON_PT_CUT = atof(value.c_str()); 
        if( key == "MET_CUT"           ) MET_CUT = atof(value.c_str()); 
        //if( key == "PROC_XSEC"         ) PROC_XSEC = atof(value.c_str());
        //if( key == "PROC_NTOT"         ) PROC_NTOT = atof(value.c_str());
        if( key == "TARGET_LUMI"       ) TARGET_LUMI = atof(value.c_str());
        if( key == "ApplyEventWeights" ) applyEventWeights = ( atoi(value.c_str()) == 1 );
        if( key == "Selection"         ) SELECTION = value;
        if( key == "MetadataFileName"  ) metadataFileName = value;
        if( key == "DatasetName"       ) datasetName = value;
        if( key == "WriteOutputTree"   ) _writeOutputTree = (bool)(atoi(value.c_str()));
        if( key == "LEADING_SV_JET_CUT" ) LEADING_SV_JET_CUT = atof(value.c_str());
        if( key == "MJJ_CUT"           ) MJJ_CUT = atof(value.c_str());
    }
    outputHistos = datasetName + "_histos.root";
    //_outputHistos = datasetName + "_histos.root";
    _outputFileName = datasetName + "_tree.root";
    _outputFile = new TFile( _outputFileName.c_str(), "RECREATE" );
    _outputTree = new TTree( _inputTreeName.c_str(), _inputTreeName.c_str() );

    bool foundMetadata = false;
    ifstream metadataFile( metadataFileName.c_str(), ios::in );
    while( metadataFile.good() ) {
      string name, xs, ntot;
      metadataFile >> name >> ws >> xs >> ws >> ntot;
      if( metadataFile.eof() ) break;
      if( name == datasetName ) {
        PROC_XSEC = atof( xs.c_str() );
        PROC_NTOT = atof( ntot.c_str() );
        foundMetadata = true;
        break;
      }
    }
    if( !foundMetadata ) {
      cout << "Did not find dataset " << datasetName << " in " << metadataFileName << ". Using standard values for xsec and ntot: " << PROC_XSEC << " " << PROC_NTOT << endl;
    }
    evtWeight = applyEventWeights ? PROC_XSEC/PROC_NTOT * TARGET_LUMI : 1.;
    std::cout << "evtWeight is " << evtWeight << std::endl;
}

void LLGAnalysis::MakeEfficiencyPlot( TH1D hpass, TH1D htotal, TCanvas *c, string triggerName ) {
    
    TGraphAsymmErrors geff;
    geff.BayesDivide( &hpass, &htotal );
    geff.GetXaxis()->SetTitle( hpass.GetXaxis()->GetTitle() );
    string ytitle = "#varepsilon (" + triggerName + ")";
    geff.GetYaxis()->SetTitle( ytitle.c_str() );
    string efftitle = "efficiency_" + triggerName;
    geff.SetNameTitle(efftitle.c_str(), efftitle.c_str());
    geff.SetMarkerColor(kBlue);
    geff.Draw("APZ"); 
    for( vector<string>::iterator itr_f = _plotFormats.begin(); itr_f != _plotFormats.end(); ++itr_f ) {
        string thisPlotName = efftitle + (*itr_f);
        c->Print( thisPlotName.c_str() );
    }

}

vector<double> LLGAnalysis::CalculateVertex( vector<double> x, vector<double> y, vector<double> z, vector<double> weight, vector<int> charge, vector<double> distance, unsigned int &nConsidered, double &weightednConsidered, vector<double> &error ) {
   
   nConsidered = 0;
   vector<double> diff_x;
   vector<double> diff_y;
   vector<double> diff_z;
   vector<double> score;
   
   for( unsigned int i = 0; i < x.size(); ++i ) {
      if( charge.at(i) == 0 ) continue;
      nConsidered += 1;
      bool knownPoint = false;
      int iKnown = -1;
      for( unsigned int i2 = 0; i2 < diff_x.size(); ++i2 ) {
        if( fabs( diff_x.at(i2) - x.at(i) ) < 1.e-10 && fabs( diff_y.at(i2) - y.at(i) ) < 1.e-10 && fabs( diff_z.at(i2) - z.at(i) ) < 1.e-10 ) {
            knownPoint = true;
            iKnown = i2;
        }
      }
    
      if( knownPoint ) {
        score.at(iKnown) += weight.at(i)/distance.at(i);
      }
      else {
        diff_x.push_back( x.at(i) );
        diff_y.push_back( y.at(i) );
        diff_z.push_back( z.at(i) );
        score.push_back( weight.at(i)/distance.at(i) );
      }
   }
    
   double scoreMax = 0.;
   vector<double> mean(3, -10000.);
   for( unsigned int i = 0; i < diff_x.size(); ++i ) {
        if ( score.at(i) > scoreMax ) {
            scoreMax = score.at(i);
            mean.at(0) = diff_x.at(i);
            mean.at(1) = diff_y.at(i);
            mean.at(2) = diff_z.at(i);
        }
    }
    return mean;
}

void LLGAnalysis::makeHist( string nametitle, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, string xtitle, string ytitle, string ztitle, string drawOption, double xAxisOffset, double yAxisOffset, double zAxisOffset ) {

    TH2D hist(nametitle.c_str(), nametitle.c_str(), nbinsx, xmin, xmax, nbinsy, ymin, ymax );
    hist.GetXaxis()->SetTitle( xtitle.c_str() );
    hist.GetYaxis()->SetTitle( ytitle.c_str() );
    hist.GetZaxis()->SetTitle( ztitle.c_str() );
    hist.GetXaxis()->SetTitleOffset( xAxisOffset );
    hist.GetYaxis()->SetTitleOffset( yAxisOffset );
    hist.GetZaxis()->SetTitleOffset( zAxisOffset );
    _histograms2D.insert( pair<string, TH2D>( nametitle, hist ) );
    _histograms2DDrawOptions.insert( pair<string,string>( nametitle, drawOption ) );

}

void LLGAnalysis::makeHist( string nametitle, int nbins, double xmin, double xmax, string xtitle, string ytitle, string drawOption, double xAxisOffset, double yAxisOffset ) {

    TH1D hist(nametitle.c_str(), nametitle.c_str(), nbins, xmin, xmax );
    hist.GetXaxis()->SetTitle( xtitle.c_str() );
    hist.GetYaxis()->SetTitle( ytitle.c_str() );
    hist.GetYaxis()->SetTitleOffset( yAxisOffset );
    hist.GetXaxis()->SetTitleOffset( xAxisOffset );
    _histograms1D.insert( pair<string, TH1D>( nametitle, hist ) );
    _histograms1DDrawOptions.insert( pair<string,string>( nametitle, drawOption ) );
}

bool LLGAnalysis::Init() {
    
    gROOT->ProcessLine(".L Loader.C+");
    gROOT->ProcessLine("#include <vector>");
    gSystem->Load("Loader_C.so");
    setStyle(1.0,true,0.15);
   
    _inputTree = new TChain(_inputTreeName.c_str());
    for( vector<string>::iterator itr = _inputFileNames.begin(); itr != _inputFileNames.end(); ++itr ) {
        _inputTree -> Add( (*itr).c_str() );  
    }

    // create the histograms and add them to the list
    makeHist( "MET_allEvents", 50, 0., 2000., "MET [GeV]", "Number of Events" );
    makeHist( "MET_HLT_PFMET170_NoiseCleaned_v1", 50, 0., 2000., "MET [GeV]", "Number of Events" );
    makeHist( "jet1_pt_allEvents", 50, 0., 1000., "Leading Jet p_{T} [GeV]", "Number of Events");
    makeHist( "jet1_pt_HLT_PFJet260_v1", 50, 0., 1000., "Leading Jet p_{T} [GeV]", "Number of Events");

    // allocate memory for all the variables
    recoJet_pt = new vector<double>;
    recoJet_phi = new vector<double>;
    recoJet_eta = new vector<double>;
    recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags = new vector<double>;
    recoJet_vertex_x = new vector<double>;
    recoJet_vertex_y = new vector<double>;
    recoJet_vertex_z = new vector<double>;

    muon_px = new vector<double>;
    muon_py = new vector<double>;
    muon_pz = new vector<double>;
    muon_phi = new vector<double>;
    muon_eta = new vector<double>;
    muon_iso = new vector<double>;
    muon_isTightMuon = new vector<bool>;
    muon_isLooseMuon = new vector<bool>;
    
    electron_px = new vector<double>;
    electron_py = new vector<double>;
    electron_pz = new vector<double>;
    electron_phi = new vector<double>;
    electron_eta = new vector<double>;
    electron_iso = new vector<double>;
    electron_isVeto = new vector<bool>;
    electron_isLoose = new vector<bool>;
    electron_isMedium = new vector<bool>;
    electron_isTight = new vector<bool>;
    electron_isHEEP = new vector<bool>;

    triggerBits = new vector<int>;
    triggerNames = new vector<string>;
    /* 
    recoJet_constVertex_x = new vector<vector<double> >;
    recoJet_constVertex_y = new vector<vector<double> >;
    recoJet_constVertex_z = new vector<vector<double> >;
    recoJet_const_pt = new vector<vector<double> >;
    recoJet_const_closestVertex_dxy = new vector<vector<double> >;
    recoJet_const_closestVertex_dz = new vector<vector<double> >;
    recoJet_const_closestVertex_d = new vector<vector<double> >;
    recoJet_const_charge = new vector<vector<int> >;
    */
    recoJet_isLeptonLike = new vector<bool>;
    vertex_x = new vector<double>;
    vertex_y = new vector<double>;
    vertex_z = new vector<double>;
    vertex_dx = new vector<double>;
    vertex_dy = new vector<double>;
    vertex_dz = new vector<double>;
    vertex_nTracks = new vector<double>;
    vertex_pt = new vector<double>;
    vertex_ndof = new vector<double>;
    secVertex_x = new vector<double>;
    secVertex_y = new vector<double>;
    secVertex_z = new vector<double>;
    secVertex_dx = new vector<double>;
    secVertex_dy = new vector<double>;
    secVertex_dz = new vector<double>;
    secVertex_pt = new vector<double>;
    secVertex_ndof = new vector<double>;
    
    // and set the branch addresses
    _inputTree->SetBranchAddress("RecoMuon_px", &muon_px );
    _inputTree->SetBranchAddress("RecoMuon_py", &muon_py );
    _inputTree->SetBranchAddress("RecoMuon_pz", &muon_pz );
    _inputTree->SetBranchAddress("RecoMuon_eta", &muon_eta );
    _inputTree->SetBranchAddress("RecoMuon_phi", &muon_phi );
    _inputTree->SetBranchAddress("RecoMuon_iso", &muon_iso );
    _inputTree->SetBranchAddress("RecoMuon_isLooseMuon", &muon_isLooseMuon );
    _inputTree->SetBranchAddress("RecoMuon_isTightMuon", &muon_isTightMuon );
    _inputTree->SetBranchAddress("RecoElectron_px", &electron_px );
    _inputTree->SetBranchAddress("RecoElectron_py", &electron_py );
    _inputTree->SetBranchAddress("RecoElectron_pz", &electron_pz );
    _inputTree->SetBranchAddress("RecoElectron_eta", &electron_eta );
    _inputTree->SetBranchAddress("RecoElectron_phi", &electron_phi );
    _inputTree->SetBranchAddress("RecoElectron_iso", &electron_iso );
    _inputTree->SetBranchAddress("RecoElectron_isVeto", &electron_isVeto );
    _inputTree->SetBranchAddress("RecoElectron_isLoose", &electron_isLoose );
    _inputTree->SetBranchAddress("RecoElectron_isMedium", &electron_isMedium );
    _inputTree->SetBranchAddress("RecoElectron_isTight", &electron_isTight );
    _inputTree->SetBranchAddress("RecoElectron_isHEEP", &electron_isHEEP );
    _inputTree->SetBranchAddress("TriggerNames", &triggerNames );
    _inputTree->SetBranchAddress("TriggerBits", &triggerBits );
    _inputTree->SetBranchAddress("RecoJet_pt", &recoJet_pt );
    _inputTree->SetBranchAddress("RecoJet_eta", &recoJet_eta );
    _inputTree->SetBranchAddress("RecoJet_phi", &recoJet_phi );
    _inputTree->SetBranchAddress("RecoJet_btag_pfCombinedInclusiveSecondaryVertexV2BJetTags", &recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags );
    _inputTree->SetBranchAddress("RecoJet_btag_pfJetBProbabilityBJetTags", &recoJet_btag_jetBProbabilityBJetTags );
    _inputTree->SetBranchAddress("RecoJet_btag_pfJetProbabilityBJetTags", &recoJet_btag_jetProbabilityBJetTags );
    _inputTree->SetBranchAddress("RecoJet_btag_pfTrackCountingHighPurBJetTags", &recoJet_btag_trackCountingHighPurBJetTags );
    _inputTree->SetBranchAddress("RecoJet_btag_pfTrackCountingHighEffBJetTags", &recoJet_btag_trackCountingHighEffBJetTags );
    _inputTree->SetBranchAddress("RecoJet_Vertex_x", &recoJet_vertex_x );
    _inputTree->SetBranchAddress("RecoJet_Vertex_y", &recoJet_vertex_y );
    _inputTree->SetBranchAddress("RecoJet_Vertex_z", &recoJet_vertex_z );
    /*
    _inputTree->SetBranchAddress("RecoJet_constVertex_x", &recoJet_constVertex_x );
    _inputTree->SetBranchAddress("RecoJet_constVertex_y", &recoJet_constVertex_y );
    _inputTree->SetBranchAddress("RecoJet_constVertex_z", &recoJet_constVertex_z );
    _inputTree->SetBranchAddress("RecoJet_const_pt", &recoJet_const_pt );
    _inputTree->SetBranchAddress("RecoJet_const_charge", &recoJet_const_charge );
    _inputTree->SetBranchAddress("RecoJet_const_closestVertex_dxy", &recoJet_const_closestVertex_dxy );
    _inputTree->SetBranchAddress("RecoJet_const_closestVertex_dz", &recoJet_const_closestVertex_dz );
    _inputTree->SetBranchAddress("RecoJet_const_closestVertex_d", &recoJet_const_closestVertex_d );
    */
    _inputTree->SetBranchAddress("RecoVertex_x", &vertex_x );
    _inputTree->SetBranchAddress("RecoVertex_y", &vertex_y );
    _inputTree->SetBranchAddress("RecoVertex_z", &vertex_z );
    _inputTree->SetBranchAddress("RecoVertex_xError", &vertex_dx );
    _inputTree->SetBranchAddress("RecoVertex_yError", &vertex_dy );
    _inputTree->SetBranchAddress("RecoVertex_zError", &vertex_dz );
    _inputTree->SetBranchAddress("RecoVertex_nTracks", &vertex_nTracks );
    _inputTree->SetBranchAddress("RecoVertex_pt", &vertex_pt );
    _inputTree->SetBranchAddress("RecoVertex_ndof", &vertex_ndof );
    _inputTree->SetBranchAddress("RecoSecVertex_x", &secVertex_x );
    _inputTree->SetBranchAddress("RecoSecVertex_y", &secVertex_y );
    _inputTree->SetBranchAddress("RecoSecVertex_z", &secVertex_z );
    _inputTree->SetBranchAddress("RecoSecVertex_xError", &secVertex_dx );
    _inputTree->SetBranchAddress("RecoSecVertex_yError", &secVertex_dy );
    _inputTree->SetBranchAddress("RecoSecVertex_zError", &secVertex_dz );
    _inputTree->SetBranchAddress("RecoSecVertex_pt", &secVertex_pt );
    _inputTree->SetBranchAddress("RecoSecVertex_ndof", &secVertex_ndof );
    _inputTree->SetBranchAddress("MET", &met );
    _inputTree->SetBranchAddress("MET_x", &met_x );
    _inputTree->SetBranchAddress("MET_y", &met_y );
    
    _outputTree->Branch("RecoMuon_px", &muon_px );
    _outputTree->Branch("RecoMuon_py", &muon_py );
    _outputTree->Branch("RecoMuon_pz", &muon_pz );
    _outputTree->Branch("RecoMuon_eta", &muon_eta );
    _outputTree->Branch("RecoMuon_phi", &muon_phi );
    _outputTree->Branch("RecoMuon_iso", &muon_iso );
    _outputTree->Branch("RecoMuon_isLooseMuon", &muon_isLooseMuon );
    _outputTree->Branch("RecoMuon_isTightMuon", &muon_isTightMuon );
    _outputTree->Branch("RecoElectron_px", &electron_px );
    _outputTree->Branch("RecoElectron_py", &electron_py );
    _outputTree->Branch("RecoElectron_pz", &electron_pz );
    _outputTree->Branch("RecoElectron_eta", &electron_eta );
    _outputTree->Branch("RecoElectron_phi", &electron_phi );
    _outputTree->Branch("RecoElectron_iso", &electron_iso );
    _outputTree->Branch("RecoElectron_isVeto", &electron_isVeto );
    _outputTree->Branch("RecoElectron_isLoose", &electron_isLoose );
    _outputTree->Branch("RecoElectron_isMedium", &electron_isMedium );
    _outputTree->Branch("RecoElectron_isTight", &electron_isTight );
    _outputTree->Branch("RecoElectron_isHEEP", &electron_isHEEP );
    _outputTree->Branch("TriggerNames", &triggerNames );
    _outputTree->Branch("TriggerBits", &triggerBits );
    _outputTree->Branch("RecoJet_pt", &recoJet_pt );
    _outputTree->Branch("RecoJet_eta", &recoJet_eta );
    _outputTree->Branch("RecoJet_phi", &recoJet_phi );
    _outputTree->Branch("RecoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags", &recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags );
    _outputTree->Branch("RecoJet_btag_jetBProbabilityBJetTags", &recoJet_btag_jetBProbabilityBJetTags );
    _outputTree->Branch("RecoJet_btag_jetProbabilityBJetTags", &recoJet_btag_jetProbabilityBJetTags );
    _outputTree->Branch("RecoJet_btag_trackCountingHighPurBJetTags", &recoJet_btag_trackCountingHighPurBJetTags );
    _outputTree->Branch("RecoJet_btag_trackCountingHighEffBJetTags", &recoJet_btag_trackCountingHighEffBJetTags );
    /*
    _outputTree->Branch("RecoJet_constVertex_x", &recoJet_constVertex_x );
    _outputTree->Branch("RecoJet_constVertex_y", &recoJet_constVertex_y );
    _outputTree->Branch("RecoJet_constVertex_z", &recoJet_constVertex_z );
    _outputTree->Branch("RecoJet_const_pt", &recoJet_const_pt );
    _outputTree->Branch("RecoJet_const_charge", &recoJet_const_charge );
    _outputTree->Branch("RecoJet_const_closestVertex_dxy", &recoJet_const_closestVertex_dxy );
    _outputTree->Branch("RecoJet_const_closestVertex_dz", &recoJet_const_closestVertex_dz );
    _outputTree->Branch("RecoJet_const_closestVertex_d", &recoJet_const_closestVertex_d );
    */
    _outputTree->Branch("RecoVertex_x", &vertex_x );
    _outputTree->Branch("RecoVertex_y", &vertex_y );
    _outputTree->Branch("RecoVertex_z", &vertex_z );
    _outputTree->Branch("RecoVertex_dx", &vertex_dx );
    _outputTree->Branch("RecoVertex_dy", &vertex_dy );
    _outputTree->Branch("RecoVertex_dz", &vertex_dz );
    _outputTree->Branch("RecoVertex_nTracks", &vertex_nTracks );
    _outputTree->Branch("RecoVertex_pt", &vertex_pt );
    _outputTree->Branch("RecoVertex_ndof", &vertex_ndof );
    _outputTree->Branch("RecoSecVertex_x", &secVertex_x );
    _outputTree->Branch("RecoSecVertex_y", &secVertex_y );
    _outputTree->Branch("RecoSecVertex_z", &secVertex_z );
    _outputTree->Branch("RecoSecVertex_pt", &secVertex_pt );
    _outputTree->Branch("RecoSecVertex_ndof", &secVertex_ndof );
    _outputTree->Branch("RecoSecVertex_dx", &secVertex_dx );
    _outputTree->Branch("RecoSecVertex_dy", &secVertex_dy );
    _outputTree->Branch("RecoSecVertex_dz", &secVertex_dz );
    _outputTree->Branch("MET", &met );
    _outputTree->Branch("MET_x", &met_x );
    _outputTree->Branch("MET_y", &met_y );
  

    // crate eps, png and pdf in the end
    //_plotFormats.push_back(".eps");
    //_plotFormats.push_back(".png");
    //_plotFormats.push_back(".pdf");

    // finally set the style
    setStyle();

    return true;
}

void LLGAnalysis::RunEventLoop( int nEntriesMax ) {

    if( nEntriesMax < 0 ) nEntriesMax = _inputTree -> GetEntries();


    if( SELECTION == "SignalRegion" ) SetupSignalRegion();
    else if( SELECTION == "TTJetsCR" ) SetupTTJetsCR();
    // SETUP YOUR SELECTION HERE

    else {
      std::cout << "Unknown selection requested. Exiting. " << std::endl;
      return;
    }

    for( int i = 0; i < nEntriesMax; ++i ) {
        
        //cout << "NOW RUNNING EVENT " << i << "\r"; fflush(stdout);
        //cout << "====================" << endl;
    
        _inputTree->GetEntry(i);
        RunObjectID(); 
        FillEfficiencyHistograms();

        if( SELECTION == "SignalRegion" ) SignalRegionSelection();
        else if( SELECTION == "TTJetsCR" ) TTJetsCRSelection();
        // CALL YOUR SELECTION HERE

    }
    cout << endl;
    return;

}   

void LLGAnalysis::FillEfficiencyHistograms() {
        
    // don' fill these histograms with weights as they are used for the efficiency plots!
    _histograms1D.at("MET_allEvents").Fill( met );
    int leadingJet = -1;
    double leadingJetPt = 0.;
    for( unsigned int iJet = 0; iJet < recoJet_pt->size(); ++iJet ) {
        if( recoJet_pt->at(iJet) > leadingJetPt ) {
            leadingJetPt = recoJet_pt->at(iJet);
            leadingJet = iJet;
        }
    }
    // don' fill these histograms with weights as they are used for the efficiency plots!
    if( leadingJet >= 0 ) _histograms1D.at("jet1_pt_allEvents").Fill( recoJet_pt->at(leadingJet) );


    for( unsigned int iTrig = 0; iTrig < triggerNames->size(); ++iTrig ) {
            
        // don' fill these histograms with weights as they are used for the efficiency plots!
        if( triggerNames->at(iTrig) == "HLT_PFMET170_NoiseCleaned_v1" && triggerBits->at(iTrig) == 1 ) {
            _histograms1D.at("MET_HLT_PFMET170_NoiseCleaned_v1").Fill( met );
        }
        if( triggerNames->at(iTrig) == "HLT_PFJet260_v1" && triggerBits->at(iTrig) == 1 ) {
            _histograms1D.at("jet1_pt_HLT_PFJet260_v1" ).Fill( recoJet_pt->at(leadingJet) );
        }
    }
}

void LLGAnalysis::FinishRun() {

    
    TCanvas c("","");
    for( map<string,TH1D>::iterator itr_h = _histograms1D.begin(); itr_h != _histograms1D.end(); ++itr_h ) {
        (*itr_h).second.Draw( (_histograms1DDrawOptions.at((*itr_h).first)).c_str() );
        for( vector<string>::iterator itr_f = _plotFormats.begin(); itr_f != _plotFormats.end(); ++itr_f ) {
            string thisPlotName = (*itr_h).first + (*itr_f);
            c.Print( thisPlotName.c_str() );
        }
    }
    
    for( map<string,TH2D>::iterator itr_h = _histograms2D.begin(); itr_h != _histograms2D.end(); ++itr_h ) {
        (*itr_h).second.Draw( (_histograms2DDrawOptions.at((*itr_h).first)).c_str()  );
        for( vector<string>::iterator itr_f = _plotFormats.begin(); itr_f != _plotFormats.end(); ++itr_f ) {
            string thisPlotName = (*itr_h).first + (*itr_f);
            c.Print( thisPlotName.c_str() );
        }
    }
    
    MakeEfficiencyPlot( _histograms1D.at("MET_HLT_PFMET170_NoiseCleaned_v1"), _histograms1D.at("MET_allEvents"), &c, "HLT_PFMET170_NoiseCleaned_v1" ); 
    MakeEfficiencyPlot( _histograms1D.at("jet1_pt_HLT_PFJet260_v1"), _histograms1D.at("jet1_pt_allEvents"), &c, "HLT_PFJet260_v1" ); 

    cout << endl << "RECO CUT FLOW " << endl;
    cout << "-----------------------------" << endl;
    for( map<string,int>::iterator itr = _cutFlow.begin(); itr != _cutFlow.end(); ++itr ) {
        cout << (*itr).first << " " << (*itr).second * evtWeight << endl;
    }

    TFile *fHistOutput = new TFile( outputHistos.c_str() , "RECREATE" );
    for( map<string,TH1D>::iterator itr_h = _histograms1D.begin(); itr_h != _histograms1D.end(); ++itr_h ) {
      (*itr_h).second.Write();
    }
    for( map<string,TH2D>::iterator itr_h = _histograms2D.begin(); itr_h != _histograms2D.end(); ++itr_h ) {
      (*itr_h).second.Write();
    }
    fHistOutput->Close();

    gDirectory = _outputFile;
    _outputTree->Write();
    _outputFile->Close();

    /* 
    cout << "PRINTING 2D OPTIMISATION MATRIX : " << endl;
    for( unsigned int i = 0; i < _yields2DOptimisation.size(); ++i ) {
      for( unsigned int j = 0; j < _yields2DOptimisation.at(i).size(); ++j ) {
        cout << std::setw(8) << _yields2DOptimisation.at(i).at(j) << " ";
      }
      cout << endl;
    }
    */

    delete _inputTree;
    delete recoJet_pt;
    delete recoJet_phi;
    delete recoJet_eta;
    delete recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags;
    delete recoJet_btag_jetBProbabilityBJetTags;
    delete recoJet_btag_jetProbabilityBJetTags;
    delete recoJet_btag_trackCountingHighPurBJetTags;
    delete recoJet_btag_trackCountingHighEffBJetTags;
    delete recoJet_isLeptonLike;
    delete recoJet_vertex_x;
    delete recoJet_vertex_y;
    delete recoJet_vertex_z;
    delete muon_px;
    delete muon_py;
    delete muon_pz;
    delete muon_phi;
    delete muon_eta;
    delete muon_iso;
    delete muon_isTightMuon;
    delete muon_isLooseMuon;
    delete electron_px;
    delete electron_py;
    delete electron_pz;
    delete electron_phi;
    delete electron_eta;
    delete electron_iso;
    delete electron_isVeto;
    delete electron_isLoose;
    delete electron_isMedium;
    delete electron_isTight;
    delete electron_isHEEP;
    delete triggerBits;
    delete triggerNames;
    /*
    delete recoJet_constVertex_x;
    delete recoJet_constVertex_y;
    delete recoJet_constVertex_z;
    delete recoJet_const_pt;
    delete recoJet_const_closestVertex_dxy;
    delete recoJet_const_closestVertex_dz;
    delete recoJet_const_closestVertex_d;
    delete recoJet_const_charge;
    */
    delete vertex_x;
    delete vertex_y;
    delete vertex_z;
    delete vertex_dx;
    delete vertex_dy;
    delete vertex_dz;
    delete vertex_nTracks;
    delete vertex_pt;
    delete vertex_ndof;
    delete secVertex_x;
    delete secVertex_y;
    delete secVertex_z;
    delete secVertex_dx;
    delete secVertex_dy;
    delete secVertex_dz;
    delete secVertex_ndof;
    delete secVertex_pt;

} 


void LLGAnalysis::setStyle( double ytoff, bool marker, double left_margin ) {
// use plain black on white colors
Int_t icol=0;
gStyle->SetFrameBorderMode(icol);
gStyle->SetCanvasBorderMode(icol);
gStyle->SetPadBorderMode(icol);
gStyle->SetPadColor(icol);
gStyle->SetCanvasColor(icol);
gStyle->SetStatColor(icol);
gStyle->SetTitleFillColor(icol);
// set the paper & margin sizes
gStyle->SetPaperSize(20,26);
gStyle->SetPadTopMargin(0.10);
gStyle->SetPadRightMargin(0.15);
gStyle->SetPadBottomMargin(0.16);
gStyle->SetPadLeftMargin(0.15);
// use large fonts
Int_t font=62;
Double_t tsize=0.04;
gStyle->SetTextFont(font);
gStyle->SetTextSize(tsize);
gStyle->SetLabelFont(font,"x");
gStyle->SetTitleFont(font,"x");
gStyle->SetLabelFont(font,"y");
gStyle->SetTitleFont(font,"y");
gStyle->SetLabelFont(font,"z");
gStyle->SetTitleFont(font,"z");
gStyle->SetLabelSize(tsize,"x");
gStyle->SetTitleSize(tsize,"x");
gStyle->SetLabelSize(tsize,"y");
gStyle->SetTitleSize(tsize,"y");
gStyle->SetLabelSize(tsize,"z");
gStyle->SetTitleSize(tsize,"z");
gStyle->SetTitleBorderSize(0);
//use bold lines and markers
if ( marker ) {
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
}
gStyle->SetHistLineWidth(Width_t(3.));
// postscript dashes
gStyle->SetLineStyleString(2,"[12 12]");
gStyle->SetOptStat(0);
gStyle->SetOptFit(1111);
// put tick marks on top and RHS of plots
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
// DLA overrides
gStyle->SetPadLeftMargin(left_margin);
gStyle->SetPadBottomMargin(0.13);
gStyle->SetTitleYOffset(ytoff);
gStyle->SetTitleXOffset(1.0);
gStyle->SetOptTitle(0);
//gStyle->SetStatStyle(0);
//gStyle->SetStatFontSize();
gStyle->SetStatW(0.17);

}

LLGAnalysis *LLGAnalysis::_instance = 0;
