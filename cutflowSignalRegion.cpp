#include "LLGAnalysis.h"
#include "TLorentzVector.h"

void LLGAnalysis::SetupSignalRegion() {

    // setup the cutflow
    _cutFlow.insert(pair<string,int>("0_NoCut", 0) );
    _cutFlow.insert(pair<string,int>("1_Trigger", 0) );
    _cutFlow.insert(pair<string,int>("2_MuonVeto", 0) );
    _cutFlow.insert(pair<string,int>("3_ElectronVeto", 0) );
    _cutFlow.insert(pair<string,int>("4_HasPVJet", 0) );
    _cutFlow.insert(pair<string,int>("5_HasSV20", 0) );
    _cutFlow.insert(pair<string,int>("6_MET", 0) );
    _cutFlow.insert(pair<string,int>("6a_SVJetPT", 0 ) );
    _cutFlow.insert(pair<string,int>("7_DiJetMass", 0 ) );
    _cutFlow.insert(pair<string,int>("8_BVeto", 0) );
    _cutFlow.insert(pair<string,int>("9_SVPVDistance", 0) );
    
    // and the histograms 
    makeHist( "nBjetAtSV", 5, -0.5, 4.5, "Number of b-jets associated to SV", "Number of SV" );
    makeHist( "mJJSV", 100, 0., 500., "DiJet mass at SV", "Number of Jet Pairs" );
    makeHist( "nJetsSV", 7, -0.5, 6.5, "Number of Jets associated to SV", "Number of SV" );
    makeHist( "SVJet1Pt", 50, 0., 500., "SV Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "SVJet2Pt", 50, 0., 500., "SV 2^{nd} Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "SVJet3Pt", 50, 0., 500., "SV 3^{rd} Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "SVJet4Pt", 50, 0., 500., "SV 4^{th} Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "PVJet1Pt", 50, 0., 1000., "PV Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "PVJet2Pt", 50, 0., 500., "PV 2^{nd} Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "PVJet3Pt", 50, 0., 500., "PV 3^{rd} Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "PVJet4Pt", 50, 0., 500., "PV 4^{th} Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "nJetsTotal", 26, -0.5, 25.5, "Number of Jets with p_{T} > 30 GeV", "Number of Events" );
    makeHist( "BJet1Pt", 50, 0., 500., "Leading B-tagged Jet p_{T} [GeV]", "Number of Events" );
    makeHist( "BJet2Pt", 50, 0., 500., "Subleading B-tagged Jet p_{T} [GeV]", "Number of Events" );
    makeHist( "BJet3Pt", 50, 0., 500., "3^{rd} leading B-tagged Jet p_{T} [GeV]", "Number of Events" );
    makeHist( "BJet4Pt", 50, 0., 500., "4^{th} leading B-tagged Jet p_{T} [GeV]", "Number of Events" );
    makeHist( "JetLeptonDr", 50, 0., 6., "#DeltaR(jet, electron)", "# Jet-Electron Pairs" ); 
    makeHist( "nPVWithJet75", 10, -0.5, 9.5, "Number of PV with >= 1 Jet > 75 GeV", "# events" );
    makeHist( "nSVWith2Jets30", 10, -0.5, 9.5, "Number of SV with >= 2 Jets > 30 GeV", "# events" );
    makeHist( "mjjvsleadingjetpt", 100, 0., 500., 100, 0., 500., "DiJet Mass at SV [GeV]", "Leading Jet p_{T} at SV [GeV]", "Number of Events", "COLZ" );
    makeHist( "distancePVSV", 40, 0., 40., "Distance between leading PV and SV [mm]", "Number of PV-SV pairs" );

    for( double jptpv = 30.; jptpv <= 401.; jptpv += 5. ) {
      std::vector<double> yield;
      //for( double jptsv = 40.; jptsv <= 401.; jptsv += 10. ) {
      //for( double jptsv = 15.; jptsv <= 201.; jptsv += 5. ) {
      for( double jptsv = 0.; jptsv <= 50.; jptsv += 2. ) {
        yield.push_back( 0. );
      }
      _yields2DOptimisation.push_back( yield );
    }
    return;
}

void LLGAnalysis::SignalRegionSelection() {


    int IDXFIRSTCUT = -1;
    int IDXSECONDCUT = -1;


    _cutFlow.at("0_NoCut") += 1;

    
    
    bool passTrigger = false;
    for( unsigned int iTrig = 0; iTrig < triggerNames->size(); ++iTrig ) {
        if( triggerNames->at(iTrig) == "HLT_PFMET170_NoiseCleaned_v1" && triggerBits->at(iTrig) == 1 ) passTrigger = true;
    }

    if( !passTrigger ) return; 
    _cutFlow.at("1_Trigger") += 1;

    // lepton veto:
    if( vetoMuons.size() > 0 ) return; 
    _cutFlow.at("2_MuonVeto") += 1;
        
    
    if( vetoElectrons.size() > 0 ) return;
    _cutFlow.at("3_ElectronVeto") += 1;


    // now assign jets to the vertices:
    vector<int> nJetsToPV( vertex_x->size(), 0 );
    vector<int> nJetsToSV( secVertex_x->size(), 0 );
    vector<vector<int> > idJetsToPV;
    vector<vector<int> > idJetsToSV;
    for( unsigned int iVtx = 0; iVtx < vertex_x->size(); ++iVtx ) {
        vector<int> idx;
        idJetsToPV.push_back( idx );
    }
    for( unsigned int iVtx = 0; iVtx < secVertex_x->size(); ++iVtx ) {
        vector<int> idx;
        idJetsToSV.push_back( idx );
    }
    

    for( unsigned int iselJet = 0; iselJet < selectedJets.size(); ++iselJet ) {
        int iJet = selectedJets.at(iselJet);
        //calculate jet vertex position:
        //vector<double> error(3,0.);
        vector<double> position(3,0.);
        position.at(0) = recoJet_vertex_x->at(iJet);
        position.at(1) = recoJet_vertex_y->at(iJet);
        position.at(2) = recoJet_vertex_z->at(iJet);
        //vector<double> position = CalculateVertex( recoJet_constVertex_x->at(iJet), recoJet_constVertex_y->at(iJet), recoJet_constVertex_z->at(iJet), recoJet_const_pt->at(iJet), recoJet_const_charge->at(iJet), recoJet_const_closestVertex_d->at(iJet), nCons, weightednCons, error );
        int nMatch = 0;
        
        //std::cout << "this is jet # " << iJet << " at position " << position.at(0) << " " << position.at(1) << " " << position.at(2) << std::endl; 
        for( unsigned int iVtx = 0; iVtx < vertex_x->size(); ++iVtx ) {
            if( fabs(position.at(0) - vertex_x->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(1) - vertex_y->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(2) - vertex_z->at(iVtx) ) < 1.e-10 ) {
                nJetsToPV.at(iVtx) += 1;
                idJetsToPV.at(iVtx).push_back( iJet );
                nMatch += 1;
            }
        }
        for( unsigned int iVtx = 0; iVtx < secVertex_x->size(); ++iVtx ) {
            //std::cout << "checking SV with " << secVertex_x->at(iVtx) << " " << secVertex_y->at(iVtx) << " " << secVertex_z->at(iVtx) << std::endl;
            if( fabs(position.at(0) - secVertex_x->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(1) - secVertex_y->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(2) - secVertex_z->at(iVtx) ) < 1.e-10 ) {
                
                nJetsToSV.at(iVtx) += 1;
                if( recoJet_pt->at(iJet) > JET_PT_CUT_SV ) { 
                    idJetsToSV.at(iVtx).push_back( iJet );
                }
                nMatch += 1;

            }
        }
        if( nMatch > 1 ) {
            cout << "WARNING! ASSOCIATED JET TO MORE THAN 1 VERTEX ?!" << endl;
        }
    }

    // now count the number of vertices with jets:
    vector<int> PVWithJet;
    vector<int> SVWithJets;
    vector<int> SVWith2Jets;

 
    double allPVLeadingJetPt = -1.;
    double allSVLeadingJetPt = -1.;
    double allSVLeadingmJJ = -1.;
    for( unsigned int iPV = 0; iPV < vertex_x -> size(); ++iPV ) {
      bool hasJetPV = false;
      double leadingJetPt = 0.;
      for( unsigned int iiJet = 0; iiJet < idJetsToPV.at(iPV).size(); ++iiJet ) {
          int iJet = idJetsToPV.at(iPV).at(iiJet);
          if( recoJet_pt->at(iJet) > JET_PT_CUT_PV ) hasJetPV = true;
          if( recoJet_pt->at(iJet) > leadingJetPt ) leadingJetPt = recoJet_pt->at(iJet);
          if( leadingJetPt > allPVLeadingJetPt ) allPVLeadingJetPt = leadingJetPt;
      }
      if( hasJetPV ) {
        PVWithJet.push_back( iPV );
      }
    }
    
    if( PVWithJet.size() > 0 ) {
    int counter = -1;
    for( double jptpv = 30.; jptpv <= 401.; jptpv += 5. ) {
      counter++;
      if( allPVLeadingJetPt > jptpv && allPVLeadingJetPt <= jptpv + 5. ) break;
    }
    IDXFIRSTCUT = counter;
    
    for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
        if( idJetsToSV.at(iSV).size() > 0 ) SVWithJets.push_back( iSV );
        if( idJetsToSV.at(iSV).size() >= 2 ) SVWith2Jets.push_back( iSV );
    }
    }

    if( SVWith2Jets.size() > 0 ) _outputTree->Fill();


    // do the n-1 plots here:
    // 1st: Leading Lepton pT from PV:
    if( SVWith2Jets.size() > 0 && met > MET_CUT ) {
      for( unsigned int iPV = 0; iPV < vertex_x -> size(); ++iPV ) {
        int idxLeadingJetPV = -1;
        double ptLeadingJetPV = -1.;
        int idxSubLeadingJetPV = -1;
        double ptSubLeadingJetPV = -1;
        int idxThirdLeadingJetPV = -1;
        double ptThirdLeadingJetPV = -1;
        int idxFourthLeadingJetPV = -1;
        double ptFourthLeadingJetPV = -1;
        bool hasJetPV = false;
        for( unsigned int iiJet = 0; iiJet < idJetsToPV.at(iPV).size(); ++iiJet ) {
            int iJet = idJetsToPV.at(iPV).at(iiJet);
            if( recoJet_pt->at(iJet) > JET_PT_CUT_PV ) hasJetPV = true;
            
            if( recoJet_pt->at(iJet) > ptLeadingJetPV ) {
                idxFourthLeadingJetPV = idxThirdLeadingJetPV;
                ptFourthLeadingJetPV = ptThirdLeadingJetPV;
                idxThirdLeadingJetPV = idxSubLeadingJetPV;
                ptThirdLeadingJetPV = ptSubLeadingJetPV;
                idxSubLeadingJetPV = idxLeadingJetPV;
                ptSubLeadingJetPV = ptLeadingJetPV;
                idxLeadingJetPV = iJet;
                ptLeadingJetPV = recoJet_pt->at(iJet);
            }
            else if( recoJet_pt->at(iJet) > ptSubLeadingJetPV ) {
                idxFourthLeadingJetPV = idxThirdLeadingJetPV;
                ptFourthLeadingJetPV = ptThirdLeadingJetPV;
                idxThirdLeadingJetPV = idxSubLeadingJetPV;
                ptThirdLeadingJetPV = ptSubLeadingJetPV;
                idxSubLeadingJetPV = iJet;
                ptSubLeadingJetPV = recoJet_pt->at(iJet);
            }
            else if( recoJet_pt->at(iJet) > ptThirdLeadingJetPV ) {
                idxFourthLeadingJetPV = idxThirdLeadingJetPV;
                ptFourthLeadingJetPV = ptThirdLeadingJetPV;
                idxThirdLeadingJetPV = iJet;
                ptThirdLeadingJetPV = recoJet_pt->at(iJet);
            }
            else if( recoJet_pt->at(iJet) > ptThirdLeadingJetPV ) {
                idxFourthLeadingJetPV = iJet;
                ptFourthLeadingJetPV = recoJet_pt->at(iJet);
            }
        }
        _histograms1D.at("PVJet1Pt").Fill( ptLeadingJetPV, evtWeight );
        _histograms1D.at("PVJet2Pt").Fill( ptSubLeadingJetPV, evtWeight );
        _histograms1D.at("PVJet3Pt").Fill( ptThirdLeadingJetPV, evtWeight );
        _histograms1D.at("PVJet4Pt").Fill( ptFourthLeadingJetPV, evtWeight );
      }
      _histograms1D.at("nPVWithJet75").Fill( PVWithJet.size(), evtWeight ); 
    }
    
    if( PVWithJet.size() >= 1 && met > MET_CUT ) {
      
      for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
        _histograms1D.at("nJetsSV").Fill( idJetsToSV.at(iSV).size(), evtWeight );
      }
      _histograms1D.at("nSVWith2Jets30").Fill( SVWith2Jets.size(), evtWeight );
    }
    

    // and run the selection:
    if( PVWithJet.size() >= 1 ) {
        _cutFlow.at("4_HasPVJet") += 1;
        
        if( SVWith2Jets.size() > 0 ) {
            _cutFlow.at("5_HasSV20") += 1;
                
            vector<double> allDistances;
            double maxDist = 0.;
            for( unsigned int iPV = 0; iPV < PVWithJet.size(); ++iPV ) {
                double thispv_x = vertex_x->at(PVWithJet.at(iPV));
                double thispv_y = vertex_y->at(PVWithJet.at(iPV));
                double thispv_z = vertex_z->at(PVWithJet.at(iPV));
                for( unsigned int iSV = 0; iSV < SVWithJets.size(); ++iSV ) {
                    double thissv_x = secVertex_x->at(SVWithJets.at(iSV));
                    double thissv_y = secVertex_y->at(SVWithJets.at(iSV));
                    double thissv_z = secVertex_z->at(SVWithJets.at(iSV));
                    double dx = thissv_x - thispv_x;
                    double dy = thissv_y - thispv_y;
                    double dz = thissv_z - thispv_z;
                    double dist = 10.*sqrt(dx*dx + dy*dy + dz*dz );
                    allDistances.push_back( dist );
                    if( dist > maxDist ) maxDist = dist;
                }
            }
            if( met > MET_CUT ) {
                _cutFlow.at("6_MET") += 1;
                
                for( unsigned int iDist = 0; iDist < allDistances.size(); ++iDist ) {
                  _histograms1D.at("distancePVSV").Fill(allDistances.at(iDist), evtWeight );
                }

                bool hasDiJetPair100 = false;
                bool hasSVJetCut = false;

                for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
                  if( idJetsToSV.at(iSV).size() <= 1 ) continue;
                  _histograms1D.at("nJetsSV").Fill( idJetsToSV.at(iSV).size(), evtWeight );
                  
                  int idxLeadingJet = -1;
                  double ptLeadingJet = -1.;
                  int idxSubLeadingJet = -1;
                  double ptSubLeadingJet = -1;
                  int idxThirdLeadingJet = -1;
                  double ptThirdLeadingJet = -1;
                  int idxFourthLeadingJet = -1;
                  double ptFourthLeadingJet = -1;

                  for( unsigned int iJToSV = 0; iJToSV < idJetsToSV.at(iSV).size(); ++iJToSV ) {
                    int jIdx = idJetsToSV.at(iSV).at(iJToSV);
                    if( recoJet_pt->at(jIdx) > ptLeadingJet ) {
                      idxFourthLeadingJet = idxThirdLeadingJet;
                      ptFourthLeadingJet = ptThirdLeadingJet;
                      idxThirdLeadingJet = idxSubLeadingJet;
                      ptThirdLeadingJet = ptSubLeadingJet;
                      idxSubLeadingJet = idxLeadingJet;
                      ptSubLeadingJet = ptLeadingJet;
                      idxLeadingJet = jIdx;
                      ptLeadingJet = recoJet_pt->at(jIdx);
                    }
                    else if ( recoJet_pt->at(jIdx) > ptSubLeadingJet ) {
                      idxFourthLeadingJet = idxThirdLeadingJet;
                      ptFourthLeadingJet = ptThirdLeadingJet;
                      idxThirdLeadingJet = idxSubLeadingJet;
                      ptThirdLeadingJet = ptSubLeadingJet;
                      idxSubLeadingJet = jIdx;
                      ptSubLeadingJet = recoJet_pt->at(jIdx);
                    }
                    else if ( recoJet_pt->at(jIdx) > ptThirdLeadingJet ) {
                      idxFourthLeadingJet = idxThirdLeadingJet;
                      ptFourthLeadingJet = ptThirdLeadingJet;
                      idxThirdLeadingJet = jIdx;
                      ptThirdLeadingJet = recoJet_pt->at(jIdx);
                    }
                    else if ( recoJet_pt->at(jIdx) > ptFourthLeadingJet ) {
                      idxFourthLeadingJet = jIdx;
                      ptFourthLeadingJet = recoJet_pt->at(jIdx);
                    }
                  }

                  _histograms1D.at("SVJet1Pt").Fill( ptLeadingJet, evtWeight );
                  _histograms1D.at("SVJet2Pt").Fill( ptSubLeadingJet, evtWeight );
                  _histograms1D.at("SVJet3Pt").Fill( ptThirdLeadingJet, evtWeight );
                  _histograms1D.at("SVJet4Pt").Fill( ptFourthLeadingJet, evtWeight );
                  TLorentzVector p4Jet1, p4Jet2;
                  p4Jet1.SetPtEtaPhiM( recoJet_pt->at(idxLeadingJet), recoJet_eta->at(idxLeadingJet), recoJet_phi->at(idxLeadingJet), 0. );
                  p4Jet2.SetPtEtaPhiM( recoJet_pt->at(idxSubLeadingJet), recoJet_eta->at(idxSubLeadingJet), recoJet_phi->at(idxSubLeadingJet), 0. );
                  TLorentzVector p4DiJet = p4Jet1 + p4Jet2;
                  _histograms1D.at("mJJSV").Fill( p4DiJet.M(), evtWeight );
                  
                  if( p4DiJet.M() > MJJ_CUT ) hasDiJetPair100 = true;
                  if( ptLeadingJet > LEADING_SV_JET_CUT ) hasSVJetCut = true;
                  _histograms2D.at("mjjvsleadingjetpt").Fill( p4DiJet.M(), ptLeadingJet, evtWeight ); 
                  if( ptLeadingJet > allSVLeadingJetPt ) allSVLeadingJetPt = ptLeadingJet;
                  if( p4DiJet.M() > allSVLeadingmJJ ) allSVLeadingmJJ = p4DiJet.M();
                }
                int counter = -1;
                /*
                for( double jptsv = 40.; jptsv <= 401.; jptsv += 10. ) {
                  counter++;
                  if( (allSVLeadingJetPt < jptsv && allSVLeadingJetPt >= jptsv - 10.) || allSVLeadingJetPt < 40. ) break;
                }*/
                /*for( double jptsv = 15.; jptsv <= 201.; jptsv += 5. ) {
                  counter ++;
                  if( (allSVLeadingmJJ < jptsv && allSVLeadingmJJ >= jptsv - 5. ) || allSVLeadingmJJ < 15. ) break;
                }*/
                for( double jptsv = 0.; jptsv <= 50.; jptsv += 2. ) {
                  counter++;
                  if( maxDist >= jptsv && maxDist < jptsv+2. ) break;
                }
                IDXSECONDCUT = counter;
                for( unsigned int firstCut = IDXFIRSTCUT; firstCut >= 0; --firstCut ) {
                  //for( unsigned int secondCut = IDXSECONDCUT; secondCut < _yields2DOptimisation.at(firstCut).size(); ++secondCut ) {
                  for( unsigned int secondCut = IDXSECONDCUT; secondCut >=0; --secondCut ) {
                    //std::cout << "attempting to attac at " << firstCut << " and " << secondCut << endl;
                    _yields2DOptimisation.at(firstCut).at(secondCut) += evtWeight;
                    if( secondCut == 0 ) break;
                  }
                  if( firstCut == 0  ) break;
                }
                if( hasSVJetCut )  return;
                _cutFlow.at("6a_SVJetPT") += 1;
                if( hasDiJetPair100 ) return;
                _cutFlow.at("7_DiJetMass") += 1;

                bool hasBjetFromSV = false;
                for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
                  if( idJetsToSV.at(iSV).size() <= 1 ) continue;
                  int nBjets = 0;
                  for( unsigned int iJToSV = 0; iJToSV < idJetsToSV.at(iSV).size(); ++iJToSV ) {
                    //int jIdx = idJetsToSV.at(iSV).at(iJToSV);
                    //if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(jIdx) > 0.814 ) {
                    /*
                    if( recoJet_btag_jetProbabilityBJetTags -> at(jIdx) > 0.790 ) { 
                      nBjets += 1;
                      hasBjetFromSV = true;
                    */
                    //}
                    
                  }
                  _histograms1D.at("nBjetAtSV").Fill(nBjets, evtWeight );
                }
                if( !hasBjetFromSV ) {
                  _cutFlow.at("8_BVeto") += 1;
                }
            }
        }
    }   
    return;
}
