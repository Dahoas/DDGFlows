"use strict";

class MeanCurvatureFlow {
	/**
	 * This class performs {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf mean curvature flow} on a surface mesh.
	 * @constructor module:Projects.MeanCurvatureFlow
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);
	}

	/**
	 * Builds the mean curvature flow operator.
	 * @private
	 * @method module:Projects.MeanCurvatureFlow#buildFlowOperator
	 * @param {module:LinearAlgebra.SparseMatrix} M The mass matrix of the input mesh.
	 * @param {number} h The timestep.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	buildFlowOperator(M, h) {
		// TODO
		let n = M.nCols();
  	let L = this.geometry.laplaceMatrix(this.vertexIndex);
		L.scaleBy(h);
		//Add not subtract b/c positive semi-definite
		let final = M.plus(L);
		return final; // placeholder
	}

	/**
	 * Performs mean curvature flow on the input mesh with timestep h.
	 * @method module:Projects.MeanCurvatureFlow#integrate
	 * @param {number} h The timestep.
	 */
	integrate(h) {
		// TODO
		let vertices = this.geometry.mesh.vertices;
		let n = vertices.length;
		let M = this.geometry.massMatrix(this.vertexIndex);
		let op = this.buildFlowOperator(M,h);
		let f = DenseMatrix.ones(n,3);
		//console.log(n);
		for(let v of vertices){
			let i = this.vertexIndex[v];
			let pos = this.geometry.positions[v];
			f.set(pos.x,i,0);
			f.set(pos.y,i,1);
			f.set(pos.z,i,2);
			//console.log("try1");
		}
		let ch = op.chol();
		let Mp = M.timesDense(f);
		//Mp.scaleBy(-1);
		let res = ch.solvePositiveDefinite(Mp);
		// center mesh positions around origin
		//let verts = this.geometry.positions;
		for(let v of vertices){
			let i = this.vertexIndex[v];
			this.geometry.positions[v].x = res.get(i,0);
			this.geometry.positions[v].y = res.get(i,1);
			this.geometry.positions[v].z = res.get(i,2);
		}
		normalize(this.geometry.positions, vertices, false);
	}
}
