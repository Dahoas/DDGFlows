"use strict";

class TreeCotree {
	/**
	 * This class computes the {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf tree cotree} decomposition of a surface mesh
	 * to build its {@link https://en.wikipedia.org/wiki/Homology_(mathematics)#Surfaces homology generators}.
	 * @constructor module:Projects.TreeCotree
	 * @param {module:Core.Mesh} mesh The input mesh this class acts on.
	 * @property {module:Core.Mesh} mesh The input mesh this class acts on.
	 * @property {vertexParent} vertexParent A dictionary mapping each vertex of the input mesh to
	 * its parent in the primal spanning tree.
	 * @property {faceParent} faceParent A dictionary mapping each face of the input mesh to
	 * its parent in the dual spanning tree.
	 */
	constructor(mesh) {
		this.mesh = mesh;
		this.vertexParent = {};
		this.faceParent = {};
	}

	/**
	 * Builds a primal spanning tree on a boundaryless mesh.
	 * @private
	 * @method module:Projects.TreeCotree#buildPrimalSpanningTree
	 */
	buildPrimalSpanningTree() {
		let root = this.mesh.vertices[0];
		//let marks = Array(this.mesh.vertices.length).fill(false);
		this.vertexParent[root] = root;
		let stack = [];
		for(let v of root.adjacentVertices()){
			stack.push([v,root]);
		}
		while(stack.length != 0){
			let vertex = stack.pop();
			this.vertexParent[vertex[0]] = vertex[1];
			for(let v of vertex[0].adjacentVertices()){

				if(this.vertexParent[v] == undefined){
					stack.push([v,vertex[0]]);
				}
			}
		}
	}

	/**
	 * Checks whether a halfedge is in the primal spanning tree.
	 * @private
	 * @method module:Projects.TreeCotree#inPrimalSpanningTree
	 * @param {module:Core.Halfedge} h A halfedge on the input mesh.
	 * @returns {boolean}
	 */
	inPrimalSpanningTree(h) {
		let v1 = h.vertex;
		let top =h.twin.vertex;
		return top == this.vertexParent[v1] || v1 == this.vertexParent[top];
	}

	/**
	 * Builds a dual spanning tree on a boundaryless mesh.
	 * @private
	 * @method module:Projects.TreeCotree#buildDualSpanningCotree
	 */
	buildDualSpanningCotree() {
		let root = this.mesh.faces[0];
		//let marks = Array(this.mesh.vertices.length).fill(false);
		this.faceParent[root] = root;
		let stack = [];
		for(let f of root.adjacentFaces()){
			let h = this.sharedHalfedge(f,root);
			if(!this.inPrimalSpanningTree(h)){
				stack.push([f,root]);
			}
		}
		while(stack.length != 0){
			let face = stack.pop();
			this.faceParent[face[0]] = face[1];
			for(let f of face[0].adjacentFaces()){
				let h = this.sharedHalfedge(f,face[0]);

				if(this.faceParent[f] == undefined && !this.inPrimalSpanningTree(h)){
					stack.push([f,face[0]]);
				}
			}
		}
	}

	/**
	 * Checks whether a halfedge is in the dual spanning tree.
	 * @private
	 * @method module:Projects.TreeCotree#inDualSpanningTree
	 * @param {module:Core.Halfedge} h A halfedge on the input mesh.
	 * @returns {boolean}
	 */
	inDualSpanningTree(h) {
		let f1 = h.face;
		let f2 = h.twin.face;

		return f1 === this.faceParent[f2] || f2 === this.faceParent[f1];
	}

	/**
	 * Returns a halfedge lying on the shared edge between face f and g.
	 * @private
	 * @method module:Projects.TreeCotree#sharedHalfedge
	 * @param {module:Core.Face} f A face on the input mesh.
	 * @param {module:Core.Face} g A neighboring face to f on the input mesh.
	 * @returns {module:Core.Halfedge}
	 */
	sharedHalfedge(f, g) {
		for (let h of f.adjacentHalfedges()) {
			if (h.twin.face === g) {
				return h;
			}
		}

		alert("Code should not reach here!");
		return new Halfedge();
	}

	/**
	 * Computes the {@link https://en.wikipedia.org/wiki/Homology_(mathematics)#Surfaces homology generators} of the input mesh and stores them
	 * in the {@link module:Core.Mesh Mesh}'s generators property.
	 * @method module:Projects.TreeCotree#buildGenerators
	 */
	buildGenerators() {
		// build spanning trees
		this.buildPrimalSpanningTree();
		this.buildDualSpanningCotree();

		let list = [];

		let marks = {};

		for(let f1 of this.mesh.faces){
			for(let f2 of f1.adjacentFaces()){
				let loop  = [];
				let h1 = this.sharedHalfedge(f1,f2);
				if(!this.inDualSpanningTree(h1) && !this.inPrimalSpanningTree(h1) && marks[[f1,f2]] == undefined){
					marks[[f2,f1]] = true;
					let bot = f1;
					let top = this.faceParent[f1];
					while(top != bot){
						let h = this.sharedHalfedge(bot,top);
						loop.push(h);
						bot = top;
						top = this.faceParent[top];
					}

					let temp = [];
					bot = f2;
					top = this.faceParent[f2];
					while(top != bot){
						let h = this.sharedHalfedge(bot,top);
						temp.unshift(h);
						bot = top;
						top = this.faceParent[top];
					}

					loop = loop.concat(temp);
					let h = this.sharedHalfedge(f1,f2);
					loop.unshift(h);
					list.push(loop);
				}
			}
		}

		this.mesh.generators = list;
		// TODO
	}
}
