#include "LLGAnalysis.h"
#include "TLorentzVector.h"

void LLGAnalysis::SetupTTJetsCR() {

    // setup the cutflow
    _cutFlow.insert(pair<string,int>("0_NoCut", 0) );
    _cutFlow.insert(pair<string,int>("1_Trigger", 0) );
    _cutFlow.insert(pair<string,int>("2_AtLeastOneLepton", 0) );
    _cutFlow.insert(pair<string,int>("3_MET", 0) );
    _cutFlow.insert(pair<string,int>("4_NJets", 0) );
    _cutFlow.insert(pair<string,int>("5_NBJets_CSVv2L", 0) );
    _cutFlow.insert(pair<string,int>("5_NBJets_CSVv2M", 0) );
    _cutFlow.insert(pair<string,int>("5_NBJets_CSVv2T", 0) );
    _cutFlow.insert(pair<string,int>("5_NBJets_JPL", 0) );
    _cutFlow.insert(pair<string,int>("5_NBJets_JPM", 0) );
    _cutFlow.insert(pair<string,int>("5_NBJets_JPT", 0) );

    makeHist( "met", 600, 0., 600., "MET [GeV]", "Number of Events" );
    makeHist( "nJet", 40, 0, 40, "# Jets", "Number of Events" ); 
    makeHist( "nBJet_CSVv2L", 10, 0, 10, "# B-Jets", "Number of Events" ); 
    makeHist( "nBJet_CSVv2M", 10, 0, 10, "# B-Jets", "Number of Events" );
    makeHist( "nBJet_CSVv2T", 10, 0, 10, "# B-Jets", "Number of Events" );
    makeHist( "nBJet_JPL", 10, 0, 10, "# B-Jets", "Number of Events" ); 
    makeHist( "nBJet_JPM", 10, 0, 10, "# B-Jets", "Number of Events" );
    makeHist( "nBJet_JPT", 10, 0, 10, "# B-Jets", "Number of Events" );

    makeHist( "AtLeastOneLepton_met", 600, 0., 600, "MET [GeV]", "Number of Events" );
    makeHist( "AtLeastOneLepton_nJet", 40, 0, 40, "# Jets", "Number of Events" );

    makeHist( "NJets_met", 600, 0., 600, "MET [GeV]", "Number of Events" );
    makeHist( "NJets_nJet", 40, 0, 40, "# Jets", "Number of Events" );

    makeHist( "MET_met", 600, 0., 600, "MET [GeV]", "Number of Events" );
    makeHist( "MET_nJet", 40, 0, 40, "# Jets", "Number of Events" );
    
    makeHist( "NBJets_CSVv2L_met", 600, 0., 600, "MET [GeV]", "Number of Events" );     
    makeHist( "NBJets_CSVv2L_nJet", 40, 0, 40, "# Jets", "Number of Events"); 
    
    makeHist( "NBJets_CSVv2M_met", 600, 0., 600, "MET [GeV]", "Number of Events" );     
    makeHist( "NBJets_CSVv2M_nJet", 40, 0, 40, "# Jets", "Number of Events");  
    
    makeHist( "NBJets_CSVv2T_met", 600, 0., 600, "MET [GeV]", "Number of Events" );     
    makeHist( "NBJets_CSVv2T_nJet", 40, 0, 40, "# Jets", "Number of Events"); 
    
    makeHist( "NBJets_JPL_met", 600, 0., 600, "MET [GeV]", "Number of Events" );     
    makeHist( "NBJets_JPL_nJet", 40, 0, 40, "# Jets", "Number of Events"); 
    
    makeHist( "NBJets_JPM_met", 600, 0., 600, "MET [GeV]", "Number of Events" );     
    makeHist( "NBJets_JPM_nJet", 40, 0, 40, "# Jets", "Number of Events"); 
    
    makeHist( "NBJets_JPT_met", 600, 0., 600, "MET [GeV]", "Number of Events" );     
    makeHist( "NBJets_JPT_nJet", 40, 0, 40, "# Jets", "Number of Events"); 
    
    return;

}

void LLGAnalysis::TTJetsCRSelection() {
	
		//vector<bool> recoJet_isLeptonLike;
        
    	_cutFlow.at("0_NoCut") += 1;

    	bool passTrigger = false;
    	for( unsigned int iTrig = 0; iTrig < triggerNames->size(); ++iTrig ) {
        	if( triggerNames->at(iTrig) == "HLT_PFMET170_NoiseCleaned_v1" && triggerBits->at(iTrig) == 1 ) passTrigger = true;
        	//if( (triggerNames->at(iTrig) == "HLT_PFJet260_v1" || triggerNames->at(iTrig) == "HLT_PFMET170_NoiseCleaned_v1") && triggerBits->at(iTrig) == 1 ) passTrigger = true;
        }

    	if( !passTrigger ) return; 
    	
    	_cutFlow.at("1_Trigger") += 1;
    	
    	int nJets = 0;
    	for( unsigned int iselJet = 0; iselJet < selectedJets.size(); ++iselJet ) {
        	nJets += 1;
        }
    	        
        _histograms1D.at("met").Fill(met, evtWeight );
        _histograms1D.at("nJet").Fill(nJets, evtWeight );
        
        //_outputTree->Fill();
    
        // lepton:
        if( (vetoMuons.size() + vetoElectrons.size()) < 1 ) return;
        
        _cutFlow.at("2_AtLeastOneLepton") += 1;

        _histograms1D.at("AtLeastOneLepton_met").Fill(met, evtWeight );
        _histograms1D.at("AtLeastOneLepton_nJet").Fill(nJets, evtWeight );
        
        //MET cut
        if( met > MET_CUT ) {
                
        	_cutFlow.at("3_MET") += 1;
        			
        	_histograms1D.at("MET_nJet").Fill(nJets, evtWeight );
        	_histograms1D.at("MET_met").Fill(met, evtWeight ); 
        
        	if ( nJets >= 2 ) { 
        	
        		_cutFlow.at("4_NJets") += 1;
        
        		_histograms1D.at("NJets_nJet").Fill(nJets, evtWeight );
        		_histograms1D.at("NJets_met").Fill(met, evtWeight );

        		int nBjets_CSVv2L = 0;
        		int nBjets_CSVv2M = 0;
        		int nBjets_CSVv2T = 0;
        		int nBjets_JPL = 0;
        		int nBjets_JPM = 0;
        		int nBjets_JPT = 0;
        	
        		//https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation74X
        		//JPL=0.275; JPM=0.545; JPT=0.790; CVSVv2L=0.605; CSVv2M=0.89; CSVv2T=0.97
        	
        		//for( unsigned int iselJet = 0; iselJet < selectedJets.size(); ++iselJet ) {
        		//int iJet = selectedJets.at(iselJet);
        		for( unsigned int iiJet = 0; iiJet < recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->size(); ++iiJet ){ 
        			if ( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iiJet) > 0.605 ) nBjets_CSVv2L += 1;
        			if ( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iiJet) > 0.89 ) nBjets_CSVv2M += 1;
        			if ( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iiJet) > 0.97 ) nBjets_CSVv2T += 1;
        		}
        		for( unsigned int iiJet = 0; iiJet < recoJet_btag_jetProbabilityBJetTags ->size(); ++iiJet ){ 
        			if ( recoJet_btag_jetProbabilityBJetTags ->at(iiJet) > 0.275 ) nBjets_JPL += 1;
        			if ( recoJet_btag_jetProbabilityBJetTags ->at(iiJet) > 0.545 ) nBjets_JPM += 1;
        			if ( recoJet_btag_jetProbabilityBJetTags ->at(iiJet) > 0.790 ) nBjets_JPT += 1;
        		}
        		
        	
        		_histograms1D.at("nBJet_CSVv2L").Fill(nBjets_CSVv2L, evtWeight );
        		_histograms1D.at("nBJet_CSVv2M").Fill(nBjets_CSVv2M, evtWeight );
        		_histograms1D.at("nBJet_CSVv2T").Fill(nBjets_CSVv2T, evtWeight );
        		_histograms1D.at("nBJet_JPL").Fill(nBjets_JPL, evtWeight );
        		_histograms1D.at("nBJet_JPM").Fill(nBjets_JPM, evtWeight );
        		_histograms1D.at("nBJet_JPT").Fill(nBjets_JPT, evtWeight );
        	
        		if( nBjets_CSVv2L >= 2 ){ 
        			_cutFlow.at("5_NBJets_CSVv2L") += 1; 
        			_histograms1D.at("NBJets_CSVv2L_met").Fill(met, evtWeight );
        			_histograms1D.at("NBJets_CSVv2L_nJet").Fill(nJets, evtWeight );
        		}
        		
        		if( nBjets_CSVv2M >= 2 ){ 
        			_cutFlow.at("5_NBJets_CSVv2M") += 1; 
        			_histograms1D.at("NBJets_CSVv2M_met").Fill(met, evtWeight );
        			_histograms1D.at("NBJets_CSVv2M_nJet").Fill(nJets, evtWeight );
        		}
        		if( nBjets_CSVv2T >= 2 ){ 
        			_cutFlow.at("5_NBJets_CSVv2T") += 1; 
        			_histograms1D.at("NBJets_CSVv2T_met").Fill(met, evtWeight );
                	_histograms1D.at("NBJets_CSVv2T_nJet").Fill(nJets, evtWeight );
                }
                if( nBjets_JPL >= 2 ){ 
                	_cutFlow.at("5_NBJets_JPL") += 1; 
                	_histograms1D.at("NBJets_JPL_met").Fill(met, evtWeight );
                	_histograms1D.at("NBJets_JPL_nJet").Fill(nJets, evtWeight );
                } 
                if( nBjets_JPM >= 2 ){ 
                	_cutFlow.at("5_NBJets_JPM") += 1; 
                	_histograms1D.at("NBJets_JPM_met").Fill(met, evtWeight );
                	_histograms1D.at("NBJets_JPM_nJet").Fill(nJets, evtWeight );
                }
                if( nBjets_JPT >= 2 ){ 
                	_cutFlow.at("5_NBJets_JPT") += 1; 
                	_histograms1D.at("NBJets_JPT_met").Fill(met, evtWeight );
                	_histograms1D.at("NBJets_JPT_nJet").Fill(nJets, evtWeight );
                }
            }
        }

        return;
}