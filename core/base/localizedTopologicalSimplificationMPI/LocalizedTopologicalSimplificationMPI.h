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
#ifdef TTK_ENABLE_MPI
        hasMPISupport_ = true;
#endif
        this->setDebugMsgPrefix("LTS_MPI version");
      }

      ~LocalizedTopologicalSimplificationMPI() override = default;

      int preconditionTriangulation(
        ttk::AbstractTriangulation *triangulation) const {
        this->printMsg(
          "Precondition", ttk::debug::Separator::L2);
        if(triangulation) {
          triangulation->preconditionVertexNeighbors();
#ifdef TTK_ENABLE_MPI
          triangulation->preconditionExchangeGhostVertices();
#endif // TTK_ENABLE_MPI
        }
        return 0;
      }


      /// Modification of existing initializePropagations in the non-MPI version
      /// 1. Global maximum can not be in all blocks (triangulation), so no-need to check for it now
      /// 2. 
      template <typename IT, class TT>
      int initializePropagations(
        std::vector<Propagation<IT>> &propagations,
        IT *authorizationMask, // assumes preservation mask is initialized as -1
        IT *maximaBuffer,

        const IT *authorizedExtremaIndices,
        const IT &nAuthorizedExtremaIndices,
        const IT *order,
        const TT *triangulation) const {

        ttk::Timer timer;
        this->printMsg("Initializing Propagations", 0, 0, this->threadNumber_,
                       debug::LineMode::REPLACE);

        const IT nVertices = triangulation->getNumberOfVertices();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads( \
  this->threadNumber_) if(nAuthorizedExtremaIndices > 1000)
#endif // TTK_ENABLE_OPENMP
        for(IT i = 0; i < nAuthorizedExtremaIndices; i++)
          authorizationMask[authorizedExtremaIndices[i]] = -2;

        IT writeIndex = 0;
// find discarded maxima
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(IT v = 0; v < nVertices; v++) {

          // if v needs to be preserved then skip
          if(authorizationMask[v] == -2)
            continue;

          // check if v has larger neighbors
          bool hasLargerNeighbor = false;
          const IT &vOrder = order[v];
          const IT nNeighbors = triangulation->getVertexNeighborNumber(v);
          for(IT n = 0; n < nNeighbors; n++) {
            IT u{-1};
            triangulation->getVertexNeighbor(v, n, u);
            if(vOrder < order[u]) {
              hasLargerNeighbor = true;
              break;
            }
          }

          // if v has larger neighbors then v can not be maximum
          if(hasLargerNeighbor)
            continue;

          // get local write index for this thread
          IT localWriteIndex = 0;
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic capture
#endif // TTK_ENABLE_OPENMP
          localWriteIndex = writeIndex++;

// atomic write to prevent cache line conflicts
#ifdef TTK_ENABLE_OPENMP
#pragma omp atomic write
#endif // TTK_ENABLE_OPENMP
          maximaBuffer[localWriteIndex] = v;
        }

        // sort maxima in ascending order
        std::sort(maximaBuffer, maximaBuffer + writeIndex,
                  [=](const IT &a, const IT &b) -> bool {
                    return order[a] < order[b];
                  });

        // if there are no authorized maxima then always skip most persistent
        // prop
        if(nAuthorizedExtremaIndices < 1)
          writeIndex--;

        // init propagations
        propagations.resize(writeIndex);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(IT p = 0; p < writeIndex; p++) {
          propagations[p].criticalPoints.push_back(maximaBuffer[p]);
        }

        this->printMsg("Initializing Propagations ("
                         + std::to_string(writeIndex) + "|"
                         + std::to_string(nVertices) + ")",
                       1, timer.getElapsedTime(), this->threadNumber_);

        return 0;
      }


      template <typename DT, typename IT, class TT>
      int removeUnauthorizedExtrema(DT *scalars,
                                    IT *order,

                                    const TT *triangulation,
                                    const IT *authorizedExtremaIndices,
                                    const IT &nAuthorizedExtremaIndices,
                                    const bool &computePerturbation) const {
        ttk::Timer globalTimer;
#ifdef TTK_ENABLE_MPI
        if(!ttk::hasInitializedMPI())
          return -1;
#endif // TTK_ENABLE_MPI        
        IT nVertices = triangulation->getNumberOfVertices();
        
#ifdef TTK_ENABLE_MPI
        this->printMsg(std::to_string(nVertices) + " number of vertices on this block",
                         ttk::debug::Separator::L2);

        this->printMsg(std::to_string(MPIrank_) + " id block on " + std::to_string(MPIsize_),
                         ttk::debug::Separator::L2);
#endif // TTK_ENABLE_MPI
        

        return 0;

        bool hasMinima{false};
        bool hasMaxima{false};
        for(IT i = 0; i < nAuthorizedExtremaIndices; i++) {
          const auto &extremum{authorizedExtremaIndices[i]};
          const auto nNeigh{triangulation->getVertexNeighborNumber(extremum)};
          if(nNeigh > 0) {
            // look at the first neighbor to determine if minimum or maximum
            SimplexId neigh{};
            triangulation->getVertexNeighbor(extremum, 0, neigh);
            if(order[extremum] > order[neigh]) {
              hasMaxima = true;
            } else {
              hasMinima = true;
            }
          }
          // don't look further if we have both minima and maxima
          if(hasMaxima && hasMinima) {
            this->printMsg("----------- [BLOC VIDE OF EXTREMA]",
                         ttk::debug::Separator::L2);
          }
        }

        this->printMsg(debug::Separator::L1);

        return 0;
      }

    }; // LocalizedTopologicalSimplificationMPI class
  } // namespace lts
} // namespace ttk
