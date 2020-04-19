"use strict";

class HeatMethod {
	/**
	 * This class implements the {@link http://cs.cmu.edu/~kmcrane/Projects/HeatMethod/ heat method} to compute geodesic distance
	 * on a surface mesh.
	 * @constructor module:Projects.HeatMethod
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 * @property {module:LinearAlgebra.SparseMatrix} A The laplace matrix of the input mesh.
	 * @property {module:LinearAlgebra.SparseMatrix} F The mean curvature flow operator built on the input mesh.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);

		// TODO: build laplace and flow matrices
		this.A = this.geometry.laplaceMatrix(this.vertexIndex).timesReal(-1);
		let h = this.geometry.meanEdgeLength();
		let M = this.geometry.massMatrix(this.vertexIndex);
		this.F =  M.minus(this.A.timesReal(h*h));// placeholder
	}

	/**
	 * Computes the vector field X = -∇u / |∇u|.
	 * @private
	 * @method module:Projects.HeatMethod#computeVectorField
	 * @param {module:LinearAlgebra.DenseMatrix} u A dense vector (i.e., u.nCols() == 1) representing the
	 * heat that is allowed to diffuse on the input mesh for a brief period of time.
	 * @returns {Object} A dictionary mapping each face of the input mesh to a {@link module:LinearAlgebra.Vector Vector}.
	 */
	computeVectorField(u) {
		// TODO
		let dict = {};
		for(let f of this.geometry.mesh.faces){
			let N = this.geometry.faceNormal(f);
			N = N.unit();
			let area = this.geometry.area(f);
			let vect = new Vector(0,0,0);
			//How do we know orienation of face?
			for(let h of f.adjacentHalfedges()){
			//	let h = e.halfedge;
				let vert = h.corner.vertex;
				let crs = N.cross(this.geometry.vector(h));
				crs.scaleBy(u.get(this.vertexIndex[vert],0));
				vect.incrementBy(crs);
			}
			vect.divideBy(2*area);
			vect = vect.unit();
			vect = vect.negated();
			dict[f] = vect;
		}
		return dict;
	}

	/**
	 * Computes the integrated divergence ∇.X.
	 * @private
	 * @method module:Projects.HeatMethod#computeDivergence
	 * @param {Object} X The vector field -∇u / |∇u| represented by a dictionary
	 * mapping each face of the input mesh to a {@link module:LinearAlgebra.Vector Vector}.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	computeDivergence(X) {
		// TODO
		let n = this.geometry.mesh.vertices.length;
		let ret = DenseMatrix.zeros(n,1);
		for(let v of this.geometry.mesh.vertices){
			let i = this.vertexIndex[v];
			let sum = 0;

			for(let h of v.adjacentHalfedges()){

				if(!h.onBoundary){
				let f = h.face;
				let coa = this.geometry.cotan(h);
				let coaT = this.geometry.cotan(h.twin);
				if(this.vertexIndex[h.vertex] != i){
					h=h.twin;
				}
				let vect = this.geometry.vector(h);
				let Xi = X[f];
				let cdot = vect.dot(Xi);
				let res = cdot*coa;
				sum = sum+res;

				h = h.twin;
				f = h.face;
				Xi = X[f];
				cdot = vect.dot(Xi);
				res = cdot*coaT;

				sum = sum + res;
			}


			}
		/*	for(let f of v.adjacentFaces()){
				let vect = X[f];

				for(let h of f.adjacentHalfedges()){
					if(this.vertexIndex[h.vertex] == i){
						let a = this.geometry.cotan(h);
						let dot = vect.dot(this.geometry.vector(h));
						sum = sum + a*dot;
					}
					else if(this.vertexIndex[h.twin.vertex] == i){
						//Do I need to adjust this vector???
						h = h.twin;
						let a = this.geometry.cotan(h);
						let dot = vect.dot(this.geometry.vector(h));
					}
				}

			}*/
			sum = sum/2;
			ret.set(sum,i,0);
		}
		return ret;
	}

	/**
	 * Shifts φ such that its minimum value is zero.
	 * @private
	 * @method module:Projects.HeatMethod#subtractMinimumDistance
	 * @param {module:LinearAlgebra.DenseMatrix} phi The (minimum 0) solution to the poisson equation Δφ = ∇.X.
	 */
	subtractMinimumDistance(phi) {
		let min = Infinity;
		for (let i = 0; i < phi.nRows(); i++) {
			min = Math.min(phi.get(i, 0), min);
		}

		for (let i = 0; i < phi.nRows(); i++) {
			phi.set(phi.get(i, 0) - min, i, 0);
		}
	}

	/**
	 * Computes the geodesic distances φ using the heat method.
	 * @method module:Projects.HeatMethod#compute
	 * @param {module:LinearAlgebra.DenseMatrix} delta A dense vector (i.e., delta.nCols() == 1) containing
	 * heat sources, i.e., u0 = δ(x).
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	compute(delta) {
		// TODO
		let Lc = this.A;
		let flow = this.F;
		let solve1 = flow.chol();
		let sol1 = solve1.solvePositiveDefinite(delta);

		let phi1 = this.computeDivergence(this.computeVectorField(sol1));

		let solve2 = Lc.chol();
		let phi = solve2.solvePositiveDefinite(phi1);
		// since φ is unique up to an additive constant, it should
		// be shifted such that the smallest distance is zero
		this.subtractMinimumDistance(phi);

		return phi;
	}
}
