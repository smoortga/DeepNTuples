// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//ROOT includes
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include <map>
#include "TSystem.h"
#include <iterator>

//CMSSW includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
//#include "RecoBTau/JetTagComputer/interface/TagInfoMVACategorySelector.h"
//#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"

using namespace std;
class DeepNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit DeepNtuplizer(const edm::ParameterSet&);
  ~DeepNtuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data --------------------------- 
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::JetCollection>     jetToken_;
  edm::EDGetTokenT<double>     rhoToken_;
  std::string   CSVTagComputer_;
  std::string   ipTagInfoName_;
  std::string   svTagInfoName_;

  TFile *file_ = new TFile("output.root","recreate");
  TTree *tree_ = new TTree("tree","tree");

  unsigned int npv_;

  // jet variables
  float jet_pt_, jet_eta_, jet_phi_, rho_;
	map<string, float> discriminators_;
	// ********************************
	//            TAG INFOS
	// ********************************
	//jet general
	float trackJetPt_;              // track-based jet transverse momentum
	float jetNTracks_;              // tracks associated to jet
	float jetNTracksEtaRel_;        // tracks associated to jet for which trackEtaRel is calculated
	float jetNSecondaryVertices_;   // number of secondary vertices associated to the jet
	float trackSumJetEtRatio_;      // ratio of track sum transverse energy over jet energy
	float trackSumJetDeltaR_;       // pseudoangular distance between jet axis and track fourvector sum
	float trackSip2dValAboveCharm_; // track 2D signed impact parameter of first track lifting mass above charm
	float trackSip2dSigAboveCharm_; // track 2D signed impact parameter significance of first track lifting mass above charm
	float trackSip3dValAboveCharm_; // track 3D signed impact parameter of first track lifting mass above charm
	float trackSip3dSigAboveCharm_; // track 3D signed impact parameter significance of first track lifting mass above charm
	float vertexCategory_;          // category of secondary vertex (Reco, Pseudo, No)
	//track info
	vector<float> trackMomentum_;    // track momentum
	vector<float> trackEta_;         // track pseudorapidity
	vector<float> trackPhi_;         // track polar angle
	vector<float> trackPtRel_;       // track transverse momentum, relative to the jet axis
	vector<float> trackPPar_;        // track parallel momentum, along the jet axis
	vector<float> trackDeltaR_;      // track pseudoangular distance from the jet axis
	vector<float> trackPtRatio_;     // track transverse momentum, relative to the jet axis, normalized to its energy
	vector<float> trackPParRatio_;   // track parallel momentum, along the jet axis, normalized to its energy
	vector<float> trackSip2dVal_;    // track 2D signed impact parameter
	vector<float> trackSip2dSig_;    // track 2D signed impact parameter significance
	vector<float> trackSip3dVal_;    // track 3D signed impact parameter
	vector<float> trackSip3dSig_;    // track 3D signed impact parameter significance
	vector<float> trackDecayLenVal_; // track decay length
	vector<float> trackDecayLenSig_; // track decay length significance
	vector<float> trackJetDistVal_;  // minimum track approach distance to jet axis
	vector<float> trackJetDistSig_;  // minimum track approach distance to jet axis significance
	vector<float> trackEtaRel_;      // track pseudorapidity, relative to the jet axis
	//SV info
	vector<float> vertexMass_;          // mass of track sum at secondary vertex
	vector<float> vertexNTracks_;       // number of tracks at secondary vertex
	vector<float> vertexEnergyRatio_;   // ratio of energy at secondary vertex over total energy
	vector<float> vertexJetDeltaR_;     // pseudoangular distance between jet axis and secondary vertex direction
	vector<float> flightDistance2dVal_; // transverse distance between primary and secondary vertex
	vector<float> flightDistance2dSig_; // transverse distance significance between primary and secondary vertex
	vector<float> flightDistance3dVal_; // distance between primary and secondary vertex
	vector<float> flightDistance3dSig_; // distance significance between primary and secondary vertex
};


DeepNtuplizer::DeepNtuplizer(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  CSVTagComputer_(iConfig.getParameter<std::string>("CSVcomputer")),
  ipTagInfoName_(iConfig.getParameter<std::string>("ipTagInfos")),
  svTagInfoName_(iConfig.getParameter<std::string>("svTagInfos"))
{
  //now do what ever initialization is needed
  usesResource("TFileService");

  tree_->Branch("npv"    ,&npv_    ,"npv/i"    );

  // jet variables
  tree_->Branch("jet_pt", &jet_pt_);
	tree_->Branch("jet_eta", &jet_eta_);
	tree_->Branch("jet_phi", &jet_phi_);
	tree_->Branch("rho", &rho_);
	vector<string> disc_names = iConfig.getParameter<vector<string> >("bDiscriminators");
	for(auto& name : disc_names) {
		discriminators_[name] = 0.;
		tree_->Branch(name.c_str(), &discriminators_[name], (name+"/F").c_str());
	}

	// ********************************
	//            TAG INFOS
	// ********************************
	//jet general
	tree_->Branch("trackJetPt"             , &trackJetPt_             , "trackJetPt/F"             );
	tree_->Branch("jetNTracks"             , &jetNTracks_             , "jetNTracks/F"             );
	tree_->Branch("jetNTracksEtaRel"       , &jetNTracksEtaRel_       , "jetNTracksEtaRel/F"       );
	tree_->Branch("jetNSecondaryVertices"  , &jetNSecondaryVertices_  , "jetNSecondaryVertices/F"  );
	tree_->Branch("trackSumJetEtRatio"     , &trackSumJetEtRatio_     , "trackSumJetEtRatio/F"     );
	tree_->Branch("trackSumJetDeltaR"      , &trackSumJetDeltaR_      , "trackSumJetDeltaR/F"      );
	tree_->Branch("vertexCategory"         , &vertexCategory_         , "vertexCategory/F"         );
	tree_->Branch("trackSip2dValAboveCharm", &trackSip2dValAboveCharm_, "trackSip2dValAboveCharm/F");
	tree_->Branch("trackSip2dSigAboveCharm", &trackSip2dSigAboveCharm_, "trackSip2dSigAboveCharm/F");
	tree_->Branch("trackSip3dValAboveCharm", &trackSip3dValAboveCharm_, "trackSip3dValAboveCharm/F");
	tree_->Branch("trackSip3dSigAboveCharm", &trackSip3dSigAboveCharm_, "trackSip3dSigAboveCharm/F");
	//track info
	tree_->Branch("trackMomentum"   , &trackMomentum_   );
	tree_->Branch("trackEta"        , &trackEta_        );
	tree_->Branch("trackPhi"        , &trackPhi_        );
	tree_->Branch("trackPtRel"      , &trackPtRel_      );
	tree_->Branch("trackPPar"       , &trackPPar_       );
	tree_->Branch("trackDeltaR"     , &trackDeltaR_     );
	tree_->Branch("trackPtRatio"    , &trackPtRatio_    );
	tree_->Branch("trackPParRatio"  , &trackPParRatio_  );
	tree_->Branch("trackSip2dVal"   , &trackSip2dVal_   );
	tree_->Branch("trackSip2dSig"   , &trackSip2dSig_   );
	tree_->Branch("trackSip3dVal"   , &trackSip3dVal_   );
	tree_->Branch("trackSip3dSig"   , &trackSip3dSig_   );
	tree_->Branch("trackDecayLenVal", &trackDecayLenVal_);
	tree_->Branch("trackDecayLenSig", &trackDecayLenSig_);
	tree_->Branch("trackJetDistVal" , &trackJetDistVal_ );
	tree_->Branch("trackJetDistSig" , &trackJetDistSig_ );
	tree_->Branch("trackEtaRel"     , &trackEtaRel_     );
	//SV info
    tree_->Branch("vertexMass"         , &vertexMass_         );
	tree_->Branch("vertexNTracks"      , &vertexNTracks_      );
	tree_->Branch("vertexEnergyRatio"  , &vertexEnergyRatio_  );
	tree_->Branch("vertexJetDeltaR"    , &vertexJetDeltaR_    );
	tree_->Branch("flightDistance2dVal", &flightDistance2dVal_);
	tree_->Branch("flightDistance2dSig", &flightDistance2dSig_);
	tree_->Branch("flightDistance3dVal", &flightDistance3dVal_);
	tree_->Branch("flightDistance3dSig", &flightDistance3dSig_);
}


DeepNtuplizer::~DeepNtuplizer()
{
  file_->Close();
}


// ------------ method called for each event  ------------
void
DeepNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();

    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);

    edm::Handle<double> rho;
    iEvent.getByToken(rhoToken_, rho);
    rho_ = *rho;

    // clear vectors
    npv_ = vertices->size();


    edm::ESHandle<JetTagComputer> computerHandle;
    iSetup.get<JetTagComputerRecord>().get(CSVTagComputer_, computerHandle);
    const GenericMVAJetTagComputer *computer = dynamic_cast<const GenericMVAJetTagComputer*>(computerHandle.product());
        

  // loop over the jets
    for (const pat::Jet &jet : *jets) {
        jet_pt_  = jet.correctedJet("Uncorrected").pt();
        jet_eta_ = jet.eta();
        jet_phi_ = jet.phi();
        for(auto& entry : discriminators_) {
            entry.second = jet.bDiscriminator(entry.first);
        }


        std::vector<const reco::BaseTagInfo* > tagInfoPtrs;//(nTagInfos)
        JetTagComputer::TagInfoHelper helper(tagInfoPtrs);
        tagInfoPtrs.push_back(jet.tagInfoCandIP(ipTagInfoName_));
        tagInfoPtrs.push_back(jet.tagInfoCandSecondaryVertex(svTagInfoName_));

        reco::TaggingVariableList vars = computer->taggingVariables(helper);
        
        //*******************
        //
        //  jet general
        //
        //*******************
        trackJetPt_                 = ( vars.checkTag(reco::btau::trackJetPt) ? vars.get(reco::btau::trackJetPt) : -999 );
        jetNSecondaryVertices_      = ( vars.checkTag(reco::btau::jetNSecondaryVertices) ? vars.get(reco::btau::jetNSecondaryVertices) : -1 );
        trackSumJetEtRatio_         = ( vars.checkTag(reco::btau::trackSumJetEtRatio) ? vars.get(reco::btau::trackSumJetEtRatio) : -999 );
        trackSumJetDeltaR_          = ( vars.checkTag(reco::btau::trackSumJetDeltaR) ? vars.get(reco::btau::trackSumJetDeltaR) : -999 );
        vertexCategory_             = ( vars.checkTag(reco::btau::vertexCategory) ? vars.get(reco::btau::vertexCategory) : -999 );
        trackSip2dValAboveCharm_    = ( vars.checkTag(reco::btau::trackSip2dValAboveCharm) ? vars.get(reco::btau::trackSip2dValAboveCharm) : -999 );
        trackSip2dSigAboveCharm_    = ( vars.checkTag(reco::btau::trackSip2dSigAboveCharm) ? vars.get(reco::btau::trackSip2dSigAboveCharm) : -999 );
        trackSip3dValAboveCharm_    = ( vars.checkTag(reco::btau::trackSip3dValAboveCharm) ? vars.get(reco::btau::trackSip3dValAboveCharm) : -999 );
        trackSip3dSigAboveCharm_    = ( vars.checkTag(reco::btau::trackSip3dSigAboveCharm) ? vars.get(reco::btau::trackSip3dSigAboveCharm) : -999 );

        
        
        std::vector<float> tagValList;
        //*******************
        //
        //  track info  //FIXME: check that the order is the same of the charged components!
        //  FIXME: right now there are no default values in the vectors... Is this something we want? or do we want fixed size vectors with zero padding in the end?
        //
        //*******************
        // jetNTracks --> FIXME: what do we want here? number of selected tracks for the TagVars? of total number of tracks associated to the jet?
        tagValList = vars.getList(reco::btau::trackSip2dSig,false);
        jetNTracks_ = tagValList.size();

        // jetNTracksEtaRel_
        tagValList = vars.getList(reco::btau::trackEtaRel,false);
        jetNTracksEtaRel_ = tagValList.size();
        

        tagValList = vars.getList(reco::btau::trackMomentum,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackMomentum_) );
        tagValList = vars.getList(reco::btau::trackEta,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackEta_) );
        tagValList = vars.getList(reco::btau::trackPhi,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackPhi_) );
        tagValList = vars.getList(reco::btau::trackPtRel,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackPtRel_) );
        tagValList = vars.getList(reco::btau::trackPPar,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackPPar_) );
        tagValList = vars.getList(reco::btau::trackDeltaR,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackDeltaR_) );
        tagValList = vars.getList(reco::btau::trackPtRatio,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackPtRatio_) );
        tagValList = vars.getList(reco::btau::trackPParRatio,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackPParRatio_) );
        tagValList = vars.getList(reco::btau::trackSip2dVal,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackSip2dVal_) );
        tagValList = vars.getList(reco::btau::trackSip2dSig,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackSip2dSig_) );
        tagValList = vars.getList(reco::btau::trackSip3dVal,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackSip3dVal_) );
        tagValList = vars.getList(reco::btau::trackSip3dSig,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackSip3dSig_) );
        tagValList = vars.getList(reco::btau::trackDecayLenVal,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackDecayLenVal_) );
        tagValList = vars.getList(reco::btau::trackDecayLenSig,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackDecayLenSig_) );
        tagValList = vars.getList(reco::btau::trackJetDistVal,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackJetDistVal_) );
        tagValList = vars.getList(reco::btau::trackJetDistSig,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackJetDistSig_) );
        tagValList = vars.getList(reco::btau::trackEtaRel,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(trackEtaRel_) );


        //*******************
        //
        //  vertex info
        //
        //*******************
        tagValList = vars.getList(reco::btau::vertexMass,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(vertexMass_) );
        tagValList = vars.getList(reco::btau::vertexNTracks,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(vertexNTracks_) );
        tagValList = vars.getList(reco::btau::vertexEnergyRatio,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(vertexEnergyRatio_) );
        tagValList = vars.getList(reco::btau::vertexJetDeltaR,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(vertexJetDeltaR_) );
        tagValList = vars.getList(reco::btau::flightDistance2dVal,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(flightDistance2dVal_) );
        tagValList = vars.getList(reco::btau::flightDistance2dSig,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(flightDistance2dSig_) );
        tagValList = vars.getList(reco::btau::flightDistance3dVal,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(flightDistance3dVal_) );
        tagValList = vars.getList(reco::btau::flightDistance3dSig,false);
        if(tagValList.size()>0) std::copy( tagValList.begin(), tagValList.end(), back_inserter(flightDistance3dSig_) );

        
        
        //*******************
        //
        //  fill the tree
        //
        //*******************
        tree_->Fill();
        
        
        //*******************
        //
        //  clear all vectors!!!
        //
        //*******************
        
        trackMomentum_.clear();    // track momentum
        trackEta_.clear();         // track pseudorapidity
        trackPhi_.clear();         // track polar angle
        trackPtRel_.clear();       // track transverse momentum, relative to the jet axis
        trackPPar_.clear();        // track parallel momentum, along the jet axis
        trackDeltaR_.clear();      // track pseudoangular distance from the jet axis
        trackPtRatio_.clear();     // track transverse momentum, relative to the jet axis, normalized to its energy
        trackPParRatio_.clear();   // track parallel momentum, along the jet axis, normalized to its energy
        trackSip2dVal_.clear();    // track 2D signed impact parameter
        trackSip2dSig_.clear();    // track 2D signed impact parameter significance
        trackSip3dVal_.clear();    // track 3D signed impact parameter
        trackSip3dSig_.clear();    // track 3D signed impact parameter significance
        trackDecayLenVal_.clear(); // track decay length
        trackDecayLenSig_.clear(); // track decay length significance
        trackJetDistVal_.clear();  // minimum track approach distance to jet axis
        trackJetDistSig_.clear();  // minimum track approach distance to jet axis significance
        trackEtaRel_.clear();      // track pseudorapidity, relative to the jet axis

        vertexMass_.clear();          // mass of track sum at secondary vertex
        vertexNTracks_.clear();       // number of tracks at secondary vertex
        vertexEnergyRatio_.clear();   // ratio of energy at secondary vertex over total energy
        vertexJetDeltaR_.clear();     // pseudoangular distance between jet axis and secondary vertex direction
        flightDistance2dVal_.clear(); // transverse distance between primary and secondary vertex
        flightDistance2dSig_.clear(); // transverse distance significance between primary and secondary vertex
        flightDistance3dVal_.clear(); // distance between primary and secondary vertex
        flightDistance3dSig_.clear(); // distance significance between primary and secondary vertex
        
        
        
  }
}


// ------------ method called once each job just before starting event loop  ------------
void
DeepNtuplizer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DeepNtuplizer::endJob()
{
  file_->cd();
  tree_->Write();
  file_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DeepNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepNtuplizer);
