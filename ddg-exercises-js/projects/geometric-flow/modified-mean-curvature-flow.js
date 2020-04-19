"use strict";

class ModifiedMeanCurvatureFlow extends MeanCurvatureFlow {
	/**
	 * This class performs a {@link http://cs.jhu.edu/~misha/MyPapers/SGP12.pdf modified version} of {@link https://www.cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf mean curvature flow} on a surface mesh.
	 * @constructor module:Projects.ModifiedMeanCurvatureFlow
	 * @augments module:Projects.MeanCurvatureFlow
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:LinearAlgebra.SparseMatrix} A The laplace matrix of the input mesh.
	 */
	constructor(geometry) {
		super(geometry);
		this.vertexIndex = indexElements(geometry.mesh.vertices);
		// TODO: build the laplace matrix
		this.A = geometry.laplaceMatrix(this.vertexIndex); // placeholder
	}

	/**
	 * @inheritdoc
	 */
	buildFlowOperator(M, h) {
		// TODO

		let n = M.nCols();
  	let L = this.A;
		let Lp = L.timesReal(h);
		//Add not subtract b/c positive semi-definite
		let final = M.plus(Lp);
		return final; // placeholder
	}
}
