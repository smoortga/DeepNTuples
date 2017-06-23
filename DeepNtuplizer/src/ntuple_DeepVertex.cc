/*
 * ntuple_DeepVertex.cc
 *
 *  Created on: 23 June 2017
 *      Author: Seth Moortgat
 */


#include "../interface/ntuple_DeepVertex.h"
// for ivf
// #include "DataFormats/VertexReco/interface/VertexFwd.h"
// #include "DataFormats/VertexReco/interface/Vertex.h"
// #include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
// #include "DataFormats/TrackReco/interface/Track.h"
// #include "TLorentzVector.h"
// #include "DataFormats/Math/interface/deltaR.h"
// #include "TrackingTools/IPTools/interface/IPTools.h"
// #include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"


// system include files
//#include <memory>
#include "DataFormats/GeometrySurface/interface/Line.h"

// user include files
//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/one/EDAnalyzer.h"

//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"

//#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "TTree.h"


#include "DataFormats/JetReco/interface/PFJet.h"
#include <DataFormats/TrackReco/interface/Track.h>
//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

//#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"

// #include "RecoVertex/VertexPrimitives/interface/BasicVertexState.h"
// #include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
// #include "RecoVertex/VertexPrimitives/interface/VertexState.h"
// #include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
// #include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
// #include "RecoVertex/VertexTools/interface/VertexDistance.h"
// #include "RecoVertex/VertexTools/interface/SharedTracks.h"

#include <DataFormats/Candidate/interface/Candidate.h>
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/CandidatePtrTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "../interface/trackVars2.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
//#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

//#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"






ntuple_DeepVertex::ntuple_DeepVertex(double jetR):ntuple_content(jetR){
}
ntuple_DeepVertex::~ntuple_DeepVertex(){}


void ntuple_DeepVertex::getInput(const edm::ParameterSet& iConfig){

}

void ntuple_DeepVertex::initBranches(TTree* tree){
    addBranch(tree,"n_seed",&n_seed, "n_seed/I");

    addBranch(tree,"seed_pt",&seed_pt, "seed_pt[10]/D");
    addBranch(tree,"seed_eta",&seed_eta, "seed_eta[10]/D");
    addBranch(tree,"seed_phi",&seed_phi, "seed_phi[10]/D");
    addBranch(tree,"seed_mass",&seed_mass, "seed_mass[10]/D");
    
    addBranch(tree,"seed_dz", &seed_dz, "seed_dz[10]/D");
    addBranch(tree,"seed_dxy", &seed_dxy, "seed_dxy[10]/D");
    addBranch(tree,"seed_3D_ip", &seed_3D_ip, "seed_3D_ip[10]/D");
    addBranch(tree,"seed_3D_sip", &seed_3D_sip, "seed_3D_sip[10]/D");
    addBranch(tree,"seed_2D_ip", &seed_2D_ip, "seed_2D_ip[10]/D");
    addBranch(tree,"seed_2D_sip", &seed_2D_sip, "seed_2D_sip[10]/D");
    
    addBranch(tree,"seed_3D_signedIp", &seed_3D_signedIp, "seed_3D_signedIp[10]/D");
    addBranch(tree,"seed_3D_signedSip", &seed_3D_signedSip, "seed_3D_signedSip[10]/D");
    addBranch(tree,"seed_2D_signedIp", &seed_2D_signedIp, "seed_2D_signedIp[10]/D");
    addBranch(tree,"seed_2D_signedSip", &seed_2D_signedSip, "seed_2D_signedSip[10]/D");
    
    //addBranch(tree,"seed_JetMatch", &seed_JetMatch, "seed_JetMatch[10]/D");
    addBranch(tree,"seed_chi2reduced",&seed_chi2reduced, "seed_chi2reduced[10]/D");
    addBranch(tree,"seed_nPixelHits",&seed_nPixelHits, "seed_nPixelHits[10]/D");
    addBranch(tree,"seed_nHits",&seed_nHits, "seed_nHits[10]/D");
    addBranch(tree,"seed_jetAxisDistance",&seed_jetAxisDistance, "seed_jetAxisDistance[10]/D");
    addBranch(tree,"seed_jetAxisDlength",&seed_jetAxisDlength, "seed_jetAxisDlength[10]/D");
    
    
    
    // near Tracks
    addBranch(tree,"nearTracks_Nvtx", &nearTracks_Nvtx);
    addBranch(tree,"nearTracks_nTracks", &nearTracks_nTracks);
    addBranch(tree,"nearTracks_pt", &nearTracks_pt, "nearTracks_pt[200]/D");
    addBranch(tree,"nearTracks_eta", &nearTracks_eta, "nearTracks_eta[200]/D");
    addBranch(tree,"nearTracks_phi", &nearTracks_phi, "nearTracks_phi[200]/D");
    addBranch(tree,"nearTracks_mass", &nearTracks_mass, "nearTracks_mass[200]/D");
    addBranch(tree,"nearTracks_dz", &nearTracks_dz, "nearTracks_dz[200]/D");
    addBranch(tree,"nearTracks_dxy", &nearTracks_dxy, "nearTracks_dxy[200]/D");
    addBranch(tree,"nearTracks_3D_ip", &nearTracks_3D_ip, "nearTracks_3D_ip[200]/D");
    addBranch(tree,"nearTracks_3D_sip", &nearTracks_3D_sip, "nearTracks_3D_sip[200]/D");
    addBranch(tree,"nearTracks_2D_ip", &nearTracks_2D_ip, "nearTracks_2D_ip[200]/D");
    addBranch(tree,"nearTracks_2D_sip", &nearTracks_2D_sip, "nearTracks_2D_sip[200]/D");

    addBranch(tree,"nearTracks_PCAdist", &nearTracks_PCAdist, "nearTracks_PCAdist[200]/D");
    addBranch(tree,"nearTracks_PCAdsig", &nearTracks_PCAdsig, "nearTracks_PCAdsig[200]/D");
    
    addBranch(tree,"nearTracks_PCAonSeed_x", &nearTracks_PCAonSeed_x, "nearTracks_PCAonSeed_x[200]/D");
    addBranch(tree,"nearTracks_PCAonSeed_y", &nearTracks_PCAonSeed_y, "nearTracks_PCAonSeed_y[200]/D");
    addBranch(tree,"nearTracks_PCAonSeed_z", &nearTracks_PCAonSeed_z, "nearTracks_PCAonSeed_z[200]/D");

    addBranch(tree,"nearTracks_PCAonSeed_xerr", &nearTracks_PCAonSeed_xerr, "nearTracks_PCAonSeed_xerr[200]/D");
    addBranch(tree,"nearTracks_PCAonSeed_yerr", &nearTracks_PCAonSeed_yerr, "nearTracks_PCAonSeed_yerr[200]/D");
    addBranch(tree,"nearTracks_PCAonSeed_zerr", &nearTracks_PCAonSeed_zerr, "nearTracks_PCAonSeed_zerr[200]/D");

    addBranch(tree,"nearTracks_PCAonTrack_x", &nearTracks_PCAonTrack_x, "nearTracks_PCAonTrack_x[200]/D");
    addBranch(tree,"nearTracks_PCAonTrack_y", &nearTracks_PCAonTrack_y, "nearTracks_PCAonTrack_y[200]/D");
    addBranch(tree,"nearTracks_PCAonTrack_z", &nearTracks_PCAonTrack_z, "nearTracks_PCAonTrack_z[200]/D");

    addBranch(tree,"nearTracks_PCAonTrack_xerr", &nearTracks_PCAonTrack_xerr, "nearTracks_PCAonTrack_xerr[200]/D");
    addBranch(tree,"nearTracks_PCAonTrack_yerr", &nearTracks_PCAonTrack_yerr, "nearTracks_PCAonTrack_yerr[200]/D");
    addBranch(tree,"nearTracks_PCAonTrack_zerr", &nearTracks_PCAonTrack_zerr, "nearTracks_PCAonTrack_zerr[200]/D"); 

    addBranch(tree,"nearTracks_dotprodTrack", &nearTracks_dotprodTrack, "nearTracks_dotprodTrack[200]/D");
    addBranch(tree,"nearTracks_dotprodSeed", &nearTracks_dotprodSeed, "nearTracks_dotprodSeed[200]/D");
    addBranch(tree,"nearTracks_dotprodTrackSeed2D", &nearTracks_dotprodTrackSeed2D, "nearTracks_dotprodTrackSeed2D[200]/D");
    addBranch(tree,"nearTracks_dotprodTrackSeed3D", &nearTracks_dotprodTrackSeed3D, "nearTracks_dotprodTrackSeed3D[200]/D");
    addBranch(tree,"nearTracks_dotprodTrackSeedVectors2D", &nearTracks_dotprodTrackSeedVectors2D, "nearTracks_dotprodTrackSeedVectors2D[200]/D");
    addBranch(tree,"nearTracks_dotprodTrackSeedVectors3D", &nearTracks_dotprodTrackSeedVectors3D, "nearTracks_dotprodTrackSeedVectors3D[200]/D");
    
    addBranch(tree,"nearTracks_PCAonSeed_pvd", &nearTracks_PCAonSeed_pvd, "nearTracks_PCAonSeed_pvd[200]/D");
    addBranch(tree,"nearTracks_PCAonTrack_pvd", &nearTracks_PCAonTrack_pvd, "nearTracks_PCAonTrack_pvd[200]/D");
    addBranch(tree,"nearTracks_PCAjetAxis_dist",&nearTracks_PCAjetAxis_dist,"nearTracks_PCAjetAxis_dist[200]/D");
    addBranch(tree,"nearTracks_PCAjetMomenta_dotprod",&nearTracks_PCAjetMomenta_dotprod,"nearTracks_PCAjetMomenta_dotprod[200]/D");
    addBranch(tree,"nearTracks_PCAjetDirs_DEta",&nearTracks_PCAjetDirs_DEta,"nearTracks_PCAjetDirs_DEta[200]/D");
    addBranch(tree,"nearTracks_PCAjetDirs_DPhi",&nearTracks_PCAjetDirs_DPhi,"nearTracks_PCAjetDirs_DPhi[200]/D");


}


void ntuple_DeepVertex::readEvent(const edm::Event& iEvent){

    iEvent.getByToken(CandidateToken, tracks);
    iEvent.getByToken(BeamSpotToken, beamSpot);

}


void ntuple_DeepVertex::readSetup(const edm::EventSetup& iSetup){

    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

}




bool ntuple_DeepVertex::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){

    // pv info
    const reco::Vertex &pv = vertices()->at(0);
    GlobalPoint pvp(pv.x(),pv.y(),pv.z());

    
    std::vector<reco::TransientTrack> selectedTracks;
    std::vector<float> masses;
    
    
    //for(typename edm::View<pat::PackedCandidate>::const_iterator track = tracks->begin(); track != tracks->end(); ++track) {
    for(size_t k = 0; k<tracks->size(); ++k) {
        if((*tracks)[k].bestTrack() != 0 &&  (*tracks)[k].pt()>0.5 && std::fabs(pvp.z()-builder->build(tracks->ptrAt(k)).track().vz())<0.5) {
            selectedTracks.push_back(builder->build(tracks->ptrAt(k)));
        }
    }
    
    std::cout << "TEST1" << std::endl;

    double jet_radius = jetR();
    GlobalVector direction(jet.px(), jet.py(), jet.pz());
    
    unsigned int n_seeds = 0;
    for(std::vector<reco::TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){
        std::cout << "TEST2 in loop" << std::endl;
        
        //is the track in the jet cone?
        float angular_distance=std::sqrt(std::pow(jet.eta()-it->track().eta(),2) + std::pow(jet.phi()-it->track().phi(),2) );
        if (angular_distance>jet_radius) { continue; }
        
        // is it a seed track?
        
        std::pair<bool,Measurement1D> ip = IPTools::absoluteImpactParameter3D(*it, pv);        
        std::pair<bool,Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(*it, pv);
		std::pair<double, Measurement1D> jet_dist =IPTools::jetTrackDistance(*it, direction, pv);                   
        TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(it->impactPointState(),pv, direction,it->field());
		float length=999;
        if (closest.isValid()) length=(closest.globalPosition() - pvp).mag();
        
        bool is_seed_candidate = (ip.first && ip.second.value() >= 0.0 && ip.second.significance() >= 1.0 &&
            ip.second.value() <= max3DIPValue && ip.second.significance() <= max3DIPSignificance &&
            it->track().normalizedChi2()<5. && std::fabs(it->track().dxy(pv.position())) < 2 &&
            std::fabs(it->track().dz(pv.position())) < 17  && jet_dist.second.value()<0.07 && length<5. );
        
        if (!is_seed_candidate){continue;}
        
        std::pair<bool,Measurement1D> ipSigned = IPTools::signedImpactParameter3D(*it,direction, pv); 
        n_seeds++;
        
        nearTracks.clear();
        //now that we found a seed, loop over all other tracks and look for neighbours
        for(std::vector<reco::TransientTrack>::const_iterator tt = selectedTracks.begin();tt!=selectedTracks.end(); ++tt ) {
            std::cout << "TEST3 in loop in loop" << std::endl;
            VertexDistance3D distanceComputer;
            TwoTrackMinimumDistance dist;
            if(*tt==*it) continue;
            if(std::fabs(pvp.z()-tt->track().vz())>0.1) continue;
            if(dist.calculate(tt->impactPointState(),it->impactPointState())) {
                std::cout << "TEST3 in loop in loop condition forfilled" << std::endl;
                GlobalPoint ttPoint          = dist.points().first;
                GlobalError ttPointErr       = tt->impactPointState().cartesianError().position();
                GlobalPoint seedPosition     = dist.points().second;
                GlobalError seedPositionErr  = it->impactPointState().cartesianError().position();
                Measurement1D m = distanceComputer.distance(VertexState(seedPosition,seedPositionErr), VertexState(ttPoint, ttPointErr));
                GlobalPoint cp(dist.crossingPoint()); 
                
                std::cout << "TEST31" << std::endl;
                
                GlobalVector PairMomentum(it->track().px()+tt->track().px(), it->track().py()+tt->track().py(), it->track().pz()+tt->track().pz());
                GlobalVector  PCA_pv(cp-pvp);

                 std::cout << "TEST32" << std::endl;

                float PCAseedFromPV =  (dist.points().second-pvp).mag();
                float PCAtrackFromPV =  (dist.points().first-pvp).mag();               
                float distance = dist.distance();

                 std::cout << "TEST33" << std::endl;
                
                GlobalVector trackDir2D(tt->impactPointState().globalDirection().x(),tt->impactPointState().globalDirection().y(),0.); 
                GlobalVector seedDir2D(it->impactPointState().globalDirection().x(),it->impactPointState().globalDirection().y(),0.); 
                GlobalVector trackPCADir2D(dist.points().first.x()-pvp.x(),dist.points().first.y()-pvp.y(),0.); 
                GlobalVector seedPCADir2D(dist.points().second.x()-pvp.x(),dist.points().second.y()-pvp.y(),0.); 
                
                 std::cout << "TEST34" << std::endl;

                float dotprodTrack = (dist.points().first-pvp).unit().dot(tt->impactPointState().globalDirection().unit());
                float dotprodSeed = (dist.points().second-pvp).unit().dot(it->impactPointState().globalDirection().unit());                    
                float dotprodTrackSeed2D = trackDir2D.unit().dot(seedDir2D.unit());
                float dotprodTrackSeed3D = it->impactPointState().globalDirection().unit().dot(tt->impactPointState().globalDirection().unit());
                float dotprodTrackSeed2DV = trackPCADir2D.unit().dot(seedPCADir2D.unit());
                float dotprodTrackSeed3DV = (dist.points().second-pvp).unit().dot((dist.points().first-pvp).unit());
                
                 std::cout << "TEST35" << std::endl;

                std::pair<bool,Measurement1D> t_ip = IPTools::absoluteImpactParameter3D(*tt,pv);        
                std::pair<bool,Measurement1D> t_ip2d = IPTools::absoluteTransverseImpactParameter(*tt,pv);
                
                 std::cout << "TEST36" << std::endl;

                myTrack.set_values(tt->track().pt(), tt->track().eta(), tt->track().phi(),  tt->track().dz(pv.position()), tt->track().dxy(pv.position()), distance,  m.significance(), seedPosition.x(), seedPosition.y(), seedPosition.z(), seedPositionErr.cxx(), seedPositionErr.cyy(), seedPositionErr.czz(),  ttPoint.x(),  ttPoint.y(),  ttPoint.z(),  ttPointErr.cxx(),  ttPointErr.cyy(),  ttPointErr.czz(), dotprodTrack, dotprodSeed );
                std::cout << "TEST361" << std::endl;
                myTrack.set_index(-1);
                std::cout << "TEST362" << std::endl;
                myTrack.set_distances(PCAseedFromPV, PCAtrackFromPV);
                std::cout << "TEST363" << std::endl;
                myTrack.set_vars(masses[tt-selectedTracks.begin()],t_ip2d.second.value() , t_ip2d.second.significance(), t_ip.second.value() , t_ip.second.significance(), dotprodTrackSeed2D, dotprodTrackSeed3D, dotprodTrackSeed2DV, dotprodTrackSeed3DV ); 
                
                 std::cout << "TEST37" << std::endl;
                
                Line::PositionType pos(pvp);
                Line::DirectionType dir(direction);
                Line::DirectionType pairMomentumDir(PairMomentum);
                Line jetLine(pos,dir);   
                Line PCAMomentumLine(cp,pairMomentumDir);
                float PCA_JetAxis_dist=jetLine.distance(cp).mag();
                float dotprodMomenta=PairMomentum.unit().dot(direction.unit());
                float dEta=std::fabs(PCA_pv.eta()-jet.eta());
                float dPhi=std::fabs(PCA_pv.phi()-jet.phi());
                
                 std::cout << "TEST38" << std::endl;
            
                myTrack.setSeedMass(masses[it-selectedTracks.begin()]);                    
                myTrack.set_JetAxisVars(PCA_JetAxis_dist,dotprodMomenta,dEta,dPhi);
                nearTracks.push_back(myTrack);
                
                 std::cout << "TEST39" << std::endl;
          
            
            
            }
            std::cout << "TEST3 in loop in loop at the end" << std::endl;
        }            
         
        std::sort (nearTracks.begin(), nearTracks.end(), sortfunction2());
        nearTracks.resize(20);
        SortedSeedsMap.insert(std::make_pair(-ipSigned.second.significance(), std::make_pair(&(*it), nearTracks)));
            
    }
    
    std::cout << "TEST4 before multimap loop" << std::endl;
       
    unsigned int seeds_max_counter=0;
    for(std::multimap<double,std::pair<const reco::TransientTrack*,const std::vector<trackVars2> > >::const_iterator im = SortedSeedsMap.begin(); im != SortedSeedsMap.end(); im++){
        
        if(seeds_max_counter>=10) break;
        
        std::cout << "TEST5 in multimap loop" << std::endl;
        
        std::cout<<-(im->first)<< "  sorted sips" <<std::endl; 			
        std::cout << "check the seed again before filling" << im->second.first->track().phi() << " " << im->second.first->track().pt() << " " << im->second.first->track().eta() << std::endl;
        std::pair<bool,Measurement1D> ipSigned = IPTools::signedImpactParameter3D(*im->second.first,direction, pv);        
        std::pair<bool,Measurement1D> ip2dSigned = IPTools::signedTransverseImpactParameter(*im->second.first,direction, pv);  
        std::pair<bool,Measurement1D> ip = IPTools::absoluteImpactParameter3D(*im->second.first, pv);        
        std::pair<bool,Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(*im->second.first, pv);	
        
        seed_pt[seeds_max_counter]=im->second.first->track().pt();
        seed_eta[seeds_max_counter]=im->second.first->track().eta();
        seed_phi[seeds_max_counter]=im->second.first->track().phi();
        seed_mass[seeds_max_counter]=im->second.second.at(0).seedMass;
        seed_dz[seeds_max_counter]=im->second.first->track().dz(pv.position());
        seed_dxy[seeds_max_counter]=im->second.first->track().dxy(pv.position());
        seed_3D_ip[seeds_max_counter]=ip.second.value();
        seed_3D_sip[seeds_max_counter]=ip.second.significance();
        seed_2D_ip[seeds_max_counter]=ip2d.second.value();
        seed_2D_sip[seeds_max_counter]=ip2d.second.significance();
        seed_3D_signedIp[seeds_max_counter]=ipSigned.second.value();
        seed_3D_signedSip[seeds_max_counter]=ipSigned.second.significance();
        seed_2D_signedIp[seeds_max_counter]=ip2dSigned.second.value();
        seed_2D_signedSip[seeds_max_counter]=ip2dSigned.second.significance();
        //seed_JetMatch[seeds_max_counter]=seed_jet_match_i;			
        seed_chi2reduced[seeds_max_counter]=im->second.first->track().normalizedChi2();
        seed_nPixelHits[seeds_max_counter]=im->second.first->track().hitPattern().numberOfValidPixelHits();
        seed_nHits[seeds_max_counter]=im->second.first->track().hitPattern().numberOfValidHits();

        std::pair<double, Measurement1D> jet_distance =IPTools::jetTrackDistance(*im->second.first, direction, pv);
        seed_jetAxisDistance[seeds_max_counter]=std::fabs(jet_distance.second.value());

        TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(im->second.first->impactPointState(),pv, direction,im->second.first->field());
        if (closest.isValid()) seed_jetAxisDlength[seeds_max_counter]=(closest.globalPosition() - pvp).mag(); 
        else seed_jetAxisDlength[seeds_max_counter]= -99;
        
        // TODO: FILL NEAREAST VARIABLES
        
        seeds_max_counter++; 
    }
    
    std::cout << "TEST6" << std::endl;
    SortedSeedsMap.clear();
    nearTracks.clear();

    return true;
}









///helpers seldomly touched



// Measurement1D ntuple_SV::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
//     VertexDistanceXY dist;
//     reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
//     reco::Vertex svtx(svcand.vertex(), csv);
//     return dist.distance(svtx, pv);
// }
// 
// Measurement1D ntuple_SV::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
//     VertexDistance3D dist;
//     reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
//     reco::Vertex svtx(svcand.vertex(), csv);
//     return dist.distance(svtx, pv);
// }
// 
// float ntuple_SV::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv)  {
//     reco::Candidate::Vector p = sv.momentum();
//     reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
//     return p.Unit().Dot(d.Unit());
// }

