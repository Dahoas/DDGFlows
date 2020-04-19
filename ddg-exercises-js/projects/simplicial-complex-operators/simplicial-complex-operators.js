"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
              let i = 0;
              for(let v of mesh.vertices){
                v.index = i;
                i++;
              }
              i=0;
              for(let e of mesh.edges){
                e.index = i;
                i++;
              }
              i=0;
              for(let f of mesh.faces){
                f.index = i;
                i++;
              }
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
              var m = mesh.edges.length;
              var n = mesh.vertices.length;
              let T = new Triplet(m,n);
              for(let v of mesh.vertices){
                for(let e of v.adjacentEdges()){
                  T.addEntry(1,e.index,v.index);
                }
              }
              let a1 = SparseMatrix.fromTriplet(T);
              return a1;
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                var m = mesh.faces.length;
                var n = mesh.edges.length;
                let T = new Triplet(m,n);
                for(let f of mesh.faces){
                  for(let e of f.adjacentEdges()){
                    T.addEntry(1,f.index,e.index);
                  }
                }
                let a2 = SparseMatrix.fromTriplet(T);
                return a2;
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.SparseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
               let vertices = subset.vertices;
               var m = this.mesh.vertices.length;
               var n = 1;
               let M = DenseMatrix.zeros(m,n);
               for(let v of vertices){
                 M.set(1,v,0);
               }
               return M;
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.SparseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
          let edges = subset.edges;
          var m = this.mesh.edges.length;
          var n = 1;
          let M = DenseMatrix.zeros(m,n);
          for(let e of edges){
           M.set(1,e,0);
          }
          return M;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.SparseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
          let faces = subset.faces;
          var m = this.mesh.faces.length;
          var n = 1;
          let M = DenseMatrix.zeros(m,n);
          //for(let f of faces){
           //M.set(1,f,0);
          //}
          for(let i = 0; i < subset.faces.length; i++){
            M.set(1,subset.faces[i],0);
          }
          return M;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
              //let strs = Core.MeshSubset.deepCopy(subset);
              let vs = this.buildVertexVector(subset);
              let es = this.buildEdgeVector(subset);

              let inEdges = this.A0.timesDense(vs);
              let inFacesOne = this.A1.timesDense(inEdges);
              let inFacesTwo = this.A1.timesDense(es);

              for(let i = 0; i < inEdges.nRows();i++){
                if(inEdges.get(i,0) > 0){
                  subset.addEdge(i);
                }
              }

              for(let i = 0; i < inFacesOne.nRows();i++){
                if(inFacesOne.get(i,0) > 0){
                  subset.addFace(i);
                }
              }

              for(let i = 0; i < inFacesTwo.nRows();i++){
                if(inFacesTwo.get(i,0) > 0){
                  subset.addFace(i);
                }
              }
              return subset;
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {

          let fs = this.buildFaceVector(subset);
          let es = this.buildEdgeVector(subset);

          let inEdges = (this.A1.transpose()).timesDense(fs);
          let inVerticesOne = (this.A0.transpose()).timesDense(inEdges);
          let inVerticesTwo = (this.A0.transpose()).timesDense(es);

          for(let i = 0; i < inEdges.nRows();i++){
            if(inEdges.get(i,0) > 0){
              subset.addEdge(i);
            }
          }

          for(let i = 0; i < inVerticesOne.nRows();i++){
            if(inVerticesOne.get(i,0) > 0){
              subset.addVertex(i);
            }
          }

          for(let i = 0; i < inVerticesTwo.nRows();i++){
            if(inVerticesTwo.get(i,0) > 0){
              subset.addVertex(i);
            }
          }

                return subset; // placeholder
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                let temp = MeshSubset.deepCopy(subset);
                let cloOfStar = this.closure(this.star(subset));
                let starOfClo = this.star(this.closure(temp));
                cloOfStar.deleteSubset(starOfClo);
                return cloOfStar; // placeholder
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
          let temp = MeshSubset.deepCopy(subset);
          let closed = this.closure(temp);
          closed.deleteSubset(subset);
          if(closed.vertices.size == 0 && closed.edges.size == 0 && closed.faces.size == 0){
            return true;
          }
          else{
            return false;
          }
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
             if(this.isComplex(subset)){
                if(subset.faces.size> 0){
                  let temp = MeshSubset.deepCopy(subset);
                  temp.deleteVertices(subset.vertices);
                  temp.deleteEdges(subset.edges);
                  let closed = this.closure(temp);
                  subset.deleteSubset(closed);
                  if(subset.vertices.size == 0 && subset.edges.size == 0 && subset.faces.size == 0){
                    return 2;
                  }
                  else{
                    return -1;
                  }
                }
                else if(subset.edges.size > 0){
                  let temp = MeshSubset.deepCopy(subset);
                  temp.deleteVertices(subset.vertices);
                  let closed = this.closure(temp);
                  subset.deleteSubset(closed);
                  if(subset.vertices.size == 0 && subset.edges.size == 0 && subset.faces.size == 0){
                    return 1;
                  }
                  else{
                    return -1;
                  }
                }
                else {
                  return 0;
                }
              }
              else {
                return -1;
              }
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
          if(subset.faces.size> 0){
            let fs = this.buildFaceVector(subset);
            let inEdges = (this.A1.transpose()).timesDense(fs);
            for(let i = 0; i < inEdges.nRows();i++){
              if(inEdges.get(i,0) > 1){
                subset.deleteEdge(i);
              }
            }
            subset.deleteFaces(subset.faces);
          }
          else if(subset.edges.size > 0){
            let es = this.buildEdgeVector(subset);
            let inVertices = (this.A0.transpose()).timesDense(es);
            for(let i = 0; i < inVertices.nRows();i++){
              if(inVertices.get(i,0) > 1){
                subset.deleteVertex(i);
              }
            }
            subset.deleteEdges(subset.edges);
          }
        else {
          subset.deleteVertices(subset.vertices);
        }
                this.closure(subset);
                return subset; // placeholder
        }
}
