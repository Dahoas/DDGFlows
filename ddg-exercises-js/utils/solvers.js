"use strict";

/**
 * This class implements frequently used numerical algorithms such as the inverse power method.
 * @memberof module:Utils
 */
class Solvers {
	/**
	 * Computes the residual of Ax - λx, where x has unit norm and λ = x.Ax.
	 * @param {module:LinearAlgebra.ComplexSparseMatrix} A The complex sparse matrix whose eigen decomposition
	 * is being computed.
	 * @param {module:LinearAlgebra.ComplexDenseMatrix} x The current guess for the smallest eigenvector
	 * (corresponding to the smallest eigenvalue λ) of A.
	 * @returns {number}
	 */
	static residual(A, x) {
		// TODO
		let eVect = A.timesDense(x);
		let n = eVect.nRows();
		let sum = new Complex(0,0);
		for(let i = 0; i < n; i++){
			let a = x.get(i);
			let b = eVect.get(i);
			let bCong = b.conjugate();
			let res = a.timesComplex(bCong);
			sum = sum.plus(res);
		}
		sum = sum.overReal(x.norm(2)*x.norm(2));
		let eVal = x.timesComplex(sum);
		let ret = eVect.minus(eVal);
		return ret.norm(2)/(x.norm(2)); // placeholder
	}

	/**
	 * Solves Ax = λx, where λ is the smallest nonzero eigenvalue of A and x is the
	 * corresponding eigenvector. x should be initialized to a random complex dense
	 * vector (i.e., x.nCols() == 1).
	 * @param {module:LinearAlgebra.ComplexSparseMatrix} A The complex positive definite sparse matrix
	 * whose eigen decomposition needs to be computed.
	 * @returns {module:LinearAlgebra.ComplexDenseMatrix} The smallest eigenvector (corresponding to the
	 * smallest eigenvalue λ) of A.
	 */
	static solveInversePowerMethod(A) {

		let n = A.nRows();
		let y0 = ComplexDenseMatrix.random(n,1);
		while(this.residual(A,y0) > 10e-10){
			let sys = A.chol();
			y0 = sys.solvePositiveDefinite(y0);
			y0.scaleBy(new Complex(1/y0.norm(2),0));
		
			let tot = y0.sum();
			console.log(tot.re,tot.im);
			tot = tot.overReal(n);
			let con = ComplexDenseMatrix.constant(tot,n,1);
			y0.decrementBy(con);
			console.log(y0.get(0,0).re,y0.get(0,0).im);
			console.log(y0.get(1,0).re,y0.get(1,0).im);
			console.log(y0.get(2,0).re,y0.get(2,0).im);
		}
		return y0;
	}

	/**
	 * Inverts a 2x2 matrix.
	 * @param {module:LinearAlgebra.DenseMatrix} m The matrix to be inverted.
	 * @returns {module:LinearAlgebra.DenseMatrix}
	 */
	static invert2x2(m) {
		let m00 = m.get(0, 0);
		let m01 = m.get(0, 1);
		let m10 = m.get(1, 0);
		let m11 = m.get(1, 1);

		let det = m00 * m11 - m01 * m10;
		m.set(m11, 0, 0);
		m.set(m00, 1, 1);
		m.set(-m01, 0, 1);
		m.set(-m10, 1, 0);
		m.scaleBy(1.0 / det);

		return m;
	}
}
