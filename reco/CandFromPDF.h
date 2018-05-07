//File: CandFromPDF.h
//Brief: A Reconstructor that reads in a TTree from EDepNeutrons and creates NeutronCands from the MCClusters it contains.  
//       Tries to stitch together MCClusters that could have come from the same FS neutron based on a PDF of beta between clusters 
//       versus energy deposited.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//EdepNeutrons includes
#include "reco/Reconstructor.h"
#include "persistency/NeutronCand.h"
#include "persistency/MCCluster.h"

//ROOT includes
#include "TTreeReaderArray.h"

//c++ includes
#include <memory>

#ifndef RECO_CANDFROMPDF_H
#define RECO_CANDFROMPDF_H

class TH2D;

namespace reco
{
  class CandFromPDF: public plgn::Reconstructor
  {
    public:
      CandFromPDF(const plgn::Reconstructor::Config& config);
      virtual ~CandFromPDF() = default;

    protected:
      virtual bool DoReconstruct() override; //Look at what is already in the tree and do your own reconstruction.

      //Location of MCClusters that will be written to tree
      std::vector<pers::NeutronCand> fCands;

      //Location from which MCClusters will be read
      TTreeReaderArray<pers::MCCluster> fClusters;

      std::string fClusterAlgName; //Name of the cluster algorithm to be stitched

      //Configuration data
      double fTimeRes; //Time resolution for 3DST in ns
      double fPosRes; //Position resolution for 3DST in mm
      std::unique_ptr<TH2D> fBetaVsEDep; //PDF of Beta between clusters versus energy deposited in a cluster
  };
}

#endif //RECO_CANDFROMPDF_H
