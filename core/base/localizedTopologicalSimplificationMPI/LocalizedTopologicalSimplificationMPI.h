/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::LocalizedTopologicalSimplificationMPI
/// \author Sylvain Gerbaud <sylvain.gerbaud@lip6.fr>
/// \date October 2024
///
/// This module defines the %LocalizedTopologicalSimplificationMPI 
///
/// \b Related \b publication: \n
/// 'LocalizedTopologicalSimplificationMPI'
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <LocalizedTopologicalSimplification.h>

namespace ttk {
  namespace lts {
    /**
     * The LocalizedTopologicalSimplificationMPI class provides methods to
     * compute for each vertex of a triangulation the average scalar value of
     * itself and its direct neighbors.
     */
    class LocalizedTopologicalSimplificationMPI
      : virtual public Debug,
        public LocalizedTopologicalSimplification {

    public:
      LocalizedTopologicalSimplificationMPI(){
        this->setDebugMsgPrefix("LTS_MPI version");
      }

      ~LocalizedTopologicalSimplificationMPI() override = default;

      int preconditionTriangulation(
        ttk::AbstractTriangulation *triangulation) const {
        this->printMsg(
          "Precondition", 0, 0, this->threadNumber_, debug::LineMode::REPLACE);
        if(triangulation) {
          triangulation->preconditionVertexNeighbors();
//#ifdef TTK_ENABLE_MPI
//          triangulation->preconditionExchangeGhostVertices();
//#endif // TTK_ENABLE_MPI
        }
        return 0;
      }


    }; // LocalizedTopologicalSimplificationMPI class
  } // namespace lts
} // namespace ttk
