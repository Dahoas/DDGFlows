"use strict";

class HarmonicBases {
	/**
	 * This class computes the {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf harmonic bases} of a surface mesh.
	 * @constructor module:Projects.HarmonicBases
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 */
	constructor(geometry) {
		this.geometry = geometry;
	}

	/**
	 * Builds a closed, but not exact, primal 1-form ω.
	 * @private
	 * @method module:Projects.HarmonicBases#buildClosedPrimalOneForm
	 * @param {module:Core.Halfedge[]} generator An array of halfedges representing a
	 * {@link https://en.wikipedia.org/wiki/Homology_(mathematics)#Surfaces homology generator}
	 * of the input mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of the input mesh
	 * to a unique index.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	buildClosedPrimalOneForm(generator, edgeIndex) {
		// TODO
		let res = DenseMatrix.zeros(this.geometry.mesh.edges.length,1);
		let alt = true;
		for(let h of generator){
			let e = h.edge;
			let i = edgeIndex[e];
			if(alt){
				res.set(1,i,0);
				alt = false;
			}
			else{
				res.set(-1,i,0);
				alt = true;
			}
		}
		return res;
	}

	/**
	 * Computes the harmonic bases [γ1, γ2 ... γn] of the input mesh.
	 * @method module:Projects.HarmonicBases#compute
	 * @param {module:Projects.HodgeDecomposition} hodgeDecomposition A hodge decomposition object that
	 * can be used to compute the exact component of the closed, but not exact, primal
	 * 1-form ω.
	 * @returns {module:LinearAlgebra.DenseMatrix[]}
	 */
	compute(hodgeDecomposition) {
		// TODO
		let tct = new TreeCotree(this.geometry.mesh);
		tct.buildGenerators();
		let edgesIndex = indexElements(this.geometry.mesh.edges);
		let basis = [];

		for(let loop of tct.mesh.generators){
			let omega = this.buildClosedPrimalOneForm(loop,edgesIndex);
			let dalpha = hodgeDecomposition.computeExactComponent(omega);
			let dbeta = hodgeDecomposition.computeCoExactComponent(omega);
			let res = hodgeDecomposition.computeHarmonicComponent(omega, dalpha, dbeta)
			basis.push(res);
		}
		return basis; // placeholder
	}
}
