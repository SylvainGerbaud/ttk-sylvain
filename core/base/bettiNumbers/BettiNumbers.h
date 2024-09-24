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
      std::vector<UnionFind*> uf;
      //std::vector<UnionFind *> ufList;
      
      for(SimplexId i = 0; i < (SimplexId)nVertices; i++) {
        //ufList[i] = &(uf[i]);
        uf.push_back(new UnionFind());
      }
       
      // compute each edge of M
      size_t const nEdges = triangulation->getNumberOfEdges();
      for(SimplexId i = 0; i < (SimplexId)nEdges; i++) {
        SimplexId a, b;
        triangulation->getEdgeVertex(i, 0, a);
        triangulation->getEdgeVertex(i, 1, b);
        UnionFind::makeUnion(
          uf[a], 
          uf[b]);
      } 
  
      std::vector<UnionFind *> res_uf;
      for(auto _uf : uf)
      {
        auto parent = _uf->find();
        bool isfound = false;
        for(auto _ufb : res_uf)
        {
          if(parent == _ufb)
          {
            isfound = true;
          }
        }
        if(!isfound)
        {
          res_uf.push_back(parent);
        }
      }

      std::cout << std::endl << "Number of B0 "<< res_uf.size() << std::endl << std::endl; 

      // Compute B2

      // Extract edges from the surface
      std::vector<SimplexId> boundaryVertices(nVertices, -1);
      int nvertices_boundary = 0;

      for (int i = 0 ; i < nVertices ; i++){
        if(triangulation->isVertexOnBoundary(i)){           
          boundaryVertices[i]=nvertices_boundary;
          nvertices_boundary++;
        }
      }

      std::vector<SimplexId> boundaryEdges;

      for (int i = 0 ; i < nEdges ; i++){
        if(triangulation->isEdgeOnBoundary(i)){           
          boundaryEdges.push_back(i);
        }
      }

      //Counting connected components of bounding surface

      std::vector<UnionFind *> boundaryVerticeSeeds(nvertices_boundary);
      for (int i = 0 ; i < nvertices_boundary ; i++){
        boundaryVerticeSeeds[i] = new UnionFind();
      }

      for (int i = 0 ; i < boundaryEdges.size() ; i++){
        SimplexId a,b;
        triangulation->getEdgeVertex(boundaryEdges[i] , 0 , a);
        triangulation->getEdgeVertex(boundaryEdges[i] , 1, b);
        SimplexId surfaceId1=boundaryVertices[a];
        SimplexId surfaceId2=boundaryVertices[b];
        UnionFind::makeUnion(boundaryVerticeSeeds[surfaceId1], boundaryVerticeSeeds[surfaceId2]);
      }

      std::vector<UnionFind *> boundaryCC;
      for(auto seed : boundaryVerticeSeeds) {
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

      int B_2 = boundaryCC.size()-res_uf.size();

      std::cout << std::endl << "Number of B2 "<< B_2 << std::endl << std::endl; 

      // Calculer caraacteristique euler
      int x = nVertices - nEdges + triangulation->getNumberOfTriangles() - triangulation->getNumberOfCells();
      std::cout << std::endl << "Euler characteristics x "<< x << std::endl << std::endl; 

      std::cout << std::endl << "Number of B1 "<< res_uf.size() + B_2 - x << std::endl << std::endl; 

      return 1; // return success
    }

  }; // BettiNumbers class

} // namespace ttk
