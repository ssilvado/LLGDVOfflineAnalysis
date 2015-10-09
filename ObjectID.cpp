#include "LLGAnalysis.h"

void LLGAnalysis::RunObjectID() {


    vetoElectrons.clear();
    looseElectrons.clear();
    mediumElectrons.clear();
    tightElectrons.clear();
    heepElectrons.clear();
    vetoMuons.clear();
    selectedJets.clear();
    recoJet_isLeptonLike->clear();
    // first the muons 
    for( unsigned int im = 0; im < muon_px->size(); ++im ) {
        double pt = sqrt(muon_px->at(im)*muon_px->at(im) + muon_py->at(im)*muon_py->at(im));
        if( muon_iso->size() > 0 ) {
          if( muon_iso->at(im) / pt  > 0.2 ) continue;
        }
        if( pt < MUON_PT_CUT ) continue;
        
        vetoMuons.push_back(im);
    }
    
    // now the electrons
    for( unsigned int ie = 0; ie < electron_px->size(); ++ie ) {
        double pt = sqrt(electron_px->at(ie)*electron_px->at(ie) + electron_py->at(ie)*electron_py->at(ie));
        if( pt < ELECTRON_PT_CUT ) continue;
        if( electron_isVeto->at(ie) )   vetoElectrons.push_back(ie);
        if( electron_isLoose->at(ie) )  looseElectrons.push_back(ie);
        if( electron_isMedium->at(ie) ) mediumElectrons.push_back(ie);
        if( electron_isTight->at(ie) )  tightElectrons.push_back(ie);
        if( electron_isHEEP->at(ie) )   heepElectrons.push_back(ie);
    }
    
    // and the jets
    // first, identify lepton-like jets:
    for( unsigned int iJet = 0; iJet < recoJet_pt->size(); ++iJet ) {
      double jeta = recoJet_eta->at(iJet);
      double jphi = recoJet_phi->at(iJet);
      double drMin = 10000;
      
      // remove jets overlapping with electrons
      for( unsigned int ivetoEle = 0; ivetoEle < vetoElectrons.size(); ++ivetoEle ) {
        int iEle = vetoElectrons.at(ivetoEle);
        double eeta = electron_eta->at(iEle);
        double ephi = electron_phi->at(iEle);
        double deta = fabs(eeta - jeta);
        double dphi = fabs( ephi - jphi );
        if( dphi > M_PI ) dphi = 2*M_PI - dphi;
        double dr = sqrt( deta*deta + dphi*dphi );
        if( dr < drMin ) drMin = dr;
      }
      
      //remove jets overlapping with muons:
      for( unsigned int ivetoMuon = 0; ivetoMuon < vetoMuons.size(); ++ivetoMuon ) {
        int iMuon = vetoMuons.at(ivetoMuon);
        double eeta = muon_eta->at(iMuon);
        double ephi = muon_phi->at(iMuon);
        double deta = fabs(eeta - jeta);
        double dphi = fabs( ephi - jphi );
        if( dphi > M_PI ) dphi = 2*M_PI - dphi;
        double dr = sqrt( deta*deta + dphi*dphi );
        if( dr < drMin ) drMin = dr;
      }

      // and fill the recoJet variable:
      recoJet_isLeptonLike->push_back( (drMin < 0.4 ) ? true : false );
    }

    // now fill the signal jets:
    for( unsigned int iJet = 0; iJet < recoJet_pt->size(); ++iJet ) {
        if( recoJet_isLeptonLike->at(iJet) ) continue;
        if( fabs(recoJet_eta->at(iJet)) > JET_ETA_CUT ) continue;
        selectedJets.push_back(iJet);
    }
}
