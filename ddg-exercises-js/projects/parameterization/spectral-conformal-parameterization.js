"use strict";

class SpectralConformalParameterization {
	/**
	 * This class implements the {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf spectral conformal parameterization} algorithm to flatten
	 * surface meshes with boundaries conformally.
	 * @constructor module:Projects.SpectralConformalParameterization
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);
	}

	/**
	 * Builds the complex conformal energy matrix EC = ED - A.
	 * @private
	 * @method module:Projects.SpectralConformalParameterization#buildConformalEnergy
	 * @returns {module:LinearAlgebra.ComplexSparseMatrix}
	 */
	buildConformalEnergy() {

		let n = this.geometry.mesh.vertices.length;
		let lp = this.geometry.complexLaplaceMatrix(this.vertexIndex);
		lp = lp.timesComplex(new Complex(1/2,0));
		let T = new ComplexTriplet(n,n);

		for(let f of this.geometry.mesh.boundaries){
			let h1 = f.halfedge;
			let v1 = h1.vertex;
			let n0 = this.vertexIndex[v1];
			do{
				let h2 = h1.twin;
				let v2 = h2.vertex;
				let i = this.vertexIndex[v1];
				let j = this.vertexIndex[v2];
				let ret = (1/4)*(1);
				let num = new Complex(0,ret);
				let num2 = new Complex(0,ret*(-1));
				T.addEntry(num,i,j);
				T.addEntry(num2,j,i);
				h1 = h1.next;
				v1 = h1.vertex;
			}while(n0 != this.vertexIndex[v1]);
		}

		let A = ComplexSparseMatrix.fromTriplet(T);
		let res = lp.minus(A);
		return res; // placeholder
	}

	/**
	 * Flattens the input surface mesh with 1 or more boundaries conformally.
	 * @method module:Projects.SpectralConformalParameterization#flatten
	 * @returns {Object} A dictionary mapping each vertex to a vector of planar coordinates.
	 */
	flatten() {
		// TODO
		let vertices = this.geometry.mesh.vertices;
		//let flattening = this.geometry.positions; // placeholder
		let A = this.buildConformalEnergy();
		let vect = Solvers.solveInversePowerMethod(A);
		let n = vect.nRows();
		let dict = [];
		for(let i = 0; i < n; i++){
			let c = vect.get(i,0);
			let x = c.re;
			let y = c.im;
			let v = new Vector(x,y,0);
			dict[vertices[i]] = v;
		}
		// normalize flattening
		normalize(dict, vertices);
		return dict;
	}
}
