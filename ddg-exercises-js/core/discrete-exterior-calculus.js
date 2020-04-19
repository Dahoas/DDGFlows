"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		// TODO
		var m = geometry.mesh.vertices.length;
		var n = geometry.mesh.vertices.length;
		let T = new Triplet(m,n);
		for(let v of geometry.mesh.vertices){
			let i = vertexIndex[v];
			T.addEntry(geometry.barycentricDualArea(v),i,i);
		}
		let hs0 = SparseMatrix.fromTriplet(T);
		return hs0; // placeholder
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		// TODO
		var m = geometry.mesh.edges.length;
		var n = geometry.mesh.edges.length;
		let T = new Triplet(m,n);
		for(let e of geometry.mesh.edges){
			let i = edgeIndex[e];
			let a1 = geometry.cotan(e.halfedge);
			let a2 = geometry.cotan(e.halfedge.twin);
			let r = (a1+a2)/2;
			T.addEntry(r,i,i);
		}
		let hs1 = SparseMatrix.fromTriplet(T);
		return hs1;
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		// TODO
		var m = geometry.mesh.faces.length;
		var n = geometry.mesh.faces.length;
		let T = new Triplet(m,n);
		for(let f of geometry.mesh.faces){
			let i = faceIndex[f];
			let s = 0;
			let lens = [];
			for(let e of f.adjacentEdges()){
				s = s + geometry.length(e);
				lens.push(geometry.length(e));
			}
			s = s/2;
			let r = Math.sqrt(s*(s-lens[0])*(s-lens[1])*(s-lens[2]));
			r = 1/r;
			T.addEntry(r,i,i);
		}
		let hs2 = SparseMatrix.fromTriplet(T);
		return hs2;
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		var m = geometry.mesh.edges.length;
		var n = geometry.mesh.vertices.length;
		let T = new Triplet(m,n);
		for(let v of geometry.mesh.vertices){
			let i = vertexIndex[v];
			for(let e of v.adjacentEdges()){
				let h = e.halfedge;
				let v1 = vertexIndex[h.vertex];
				let e1 = edgeIndex[e];
				if(v1 == i){
					T.addEntry(-1,e1,i);
				}
				else{
					T.addEntry(1,e1,i);
				}
			}
		}
		let d0 = SparseMatrix.fromTriplet(T);
		return d0;
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
		var m = geometry.mesh.faces.length;
		var n = geometry.mesh.edges.length;
		let T = new Triplet(m,n);
		for(let f of geometry.mesh.faces){
			let j = faceIndex[f];
			for(let e of f.adjacentEdges()){
				let i = edgeIndex[e];
				let h = e.halfedge;
				let j1 = faceIndex[h.face];
				if(j1 == j){
					T.addEntry(1,j,i);
				}
				else{
					T.addEntry(-1,j,i);
				}
			}
		}
		let d0 = SparseMatrix.fromTriplet(T);
		return d0;
	}
}
