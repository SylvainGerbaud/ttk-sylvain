/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::BettiNumbers
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %BettiNumbers class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'BettiNumbers'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <UnionFind.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The BettiNumbers class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class BettiNumbers : virtual public Debug {

  public:
    BettiNumbers();

    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      triangulation->preconditionEdges();
      triangulation->preconditionBoundaryEdges();
      triangulation->preconditionBoundaryVertices();
      return triangulation->preconditionVertexNeighbors();
    }

    /**
     * TODO 3: Implementation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */
    template <class triangulationType = ttk::AbstractTriangulation>
    int execute(const triangulationType *triangulation) const {
      // ---------------------------------------------------------------------
      // EXERCICE TTK
      // ---------------------------------------------------------------------

      std::cout << std::endl << "Hello World!" << std::endl << std::endl;

      // ---------------------------------------------------------------------
      // B0
      // ---------------------------------------------------------------------
      size_t const nVertices = triangulation->getNumberOfVertices();
      std::vector<UnionFind> uf(nVertices);
      std::vector<UnionFind*> ufSeed(nVertices);

      for(SimplexId i = 0; i < (SimplexId)nVertices; i++) {
        ufSeed[i] = &(uf[i]);
      }
       
      // compute each edge of M
      size_t const nEdges = triangulation->getNumberOfEdges();
      for(SimplexId i = 0; i < (SimplexId)nEdges; i++) {
        SimplexId a, b;
        triangulation->getEdgeVertex(i, 0, a);
        triangulation->getEdgeVertex(i, 1, b);
        UnionFind::makeUnion(ufSeed[a], ufSeed[b]);
      } 
  
      std::vector<UnionFind *> cc;
      for(auto _uf : ufSeed)
      {
        auto parent = _uf->find();
        bool isfound = false;
        for(auto _ufb : cc)
        {
          if(parent == _ufb)
            isfound = true;
        }
        if(!isfound)
          cc.push_back(parent);
      }

     
      // Compute B2

      // Extract boundary edges from the surface
      std::vector<SimplexId> boundaryV(nVertices, -1);
      int v_boundary = 0;

      for (int i = 0 ; i < nVertices ; i++){
        if(triangulation->isVertexOnBoundary(i)){           
          boundaryV[i] = v_boundary++;
        }
      }

      std::vector<SimplexId> boundaryEdges;
      for (int i = 0 ; i < nEdges ; i++){
        if(triangulation->isEdgeOnBoundary(i)){           
          boundaryEdges.push_back(i);
        }
      }

      //Counting connected components of bounding surface

      std::vector<UnionFind> boundaryVTmp(v_boundary);
      std::vector<UnionFind *> boundaryVSeeds(v_boundary);
      for (int i = 0 ; i < v_boundary ; i++){
        boundaryVSeeds[i] = &(boundaryVTmp[i]);
      }

      for (int i = 0; i < boundaryEdges.size(); i++){
        SimplexId a,b;
        triangulation->getEdgeVertex(boundaryEdges[i] , 0 , a);
        triangulation->getEdgeVertex(boundaryEdges[i] , 1, b);
        SimplexId surfaceId1=boundaryV[a];
        SimplexId surfaceId2=boundaryV[b];
        UnionFind::makeUnion(boundaryVSeeds[surfaceId1], boundaryVSeeds[surfaceId2]);
      }

      std::vector<UnionFind *> boundaryCC;
      for(auto seed : boundaryVSeeds) {
        UnionFind* parent = seed->find();
        bool found = false;
        for (auto d : boundaryCC){
          if(d==parent){
            found=true;
            break;
          }
        }
        if(!found){
          boundaryCC.push_back(parent);
        }
      }

      int B2 = boundaryCC.size()-cc.size();


      // Compute Euler characteristic
      int x = nVertices - nEdges + triangulation->getNumberOfTriangles() - triangulation->getNumberOfCells();
      std::cout << "Euler characteristic x "<< x << std::endl; 
      std::cout << "Number of B0 "<< cc.size() << std::endl; 
      std::cout << "Number of B1 "<< cc.size() + B2 - x << std::endl; 
      std::cout << "Number of B2 "<< B2 << std::endl; 


      return 1; // return success
    }

  }; // BettiNumbers class

} // namespace ttk
