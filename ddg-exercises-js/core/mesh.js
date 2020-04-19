"use strict";

class Mesh {
	/**
	 * This class represents a Mesh.
	 * @constructor module:Core.Mesh
	 * @property {module:Core.Vertex[]} vertices The vertices contained in this mesh.
	 * @property {module:Core.Edge[]} edges The edges contained in this mesh.
	 * @property {module:Core.Face[]} faces The faces contained in this mesh.
	 * @property {module:Core.Corner[]} corners The corners contained in this mesh.
	 * @property {module:Core.Halfedge[]} halfedges The halfedges contained in this mesh.
	 * @property {module:Core.Face[]} boundaries The boundary loops contained in this mesh.
	 * @property {Array.<module:Core.Halfedge[]>} generators An array of halfedge arrays, i.e.,
	 * [[h11, h21, ..., hn1], [h12, h22, ..., hm2], ...] representing this mesh's
	 * {@link https://en.wikipedia.org/wiki/Homology_(mathematics)#Surfaces homology generators}.
	 */
	constructor() {
		this.vertices = [];
		this.edges = [];
		this.faces = [];
		this.corners = [];
		this.halfedges = [];
		this.boundaries = [];
		this.generators = [];
	}

	/**
	 * Computes the euler characteristic of this mesh.
	 * @method module:Core.Mesh#eulerCharacteristic
	 * @returns {number}
	 */
	eulerCharacteristic() {
		return this.vertices.length - this.edges.length + this.faces.length;
	}

	/**
	 * Constructs this mesh.
	 * @method module:Core.Mesh#build
	 * @param {Object} polygonSoup A polygon soup mesh containing vertex positions and indices.
	 * @param {module:LinearAlgebra.Vector[]} polygonSoup.v The vertex positions of the polygon soup mesh.
	 * @param {number[]} polygonSoup.f The indices of the polygon soup mesh.
	 * @returns {boolean} True if this mesh is constructed successfully and false if not
	 * (when this mesh contains any one or a combination of the following - non-manifold vertices,
	 *  non-manifold edges, isolated vertices, isolated faces).
	 */
	build(polygonSoup) {
		// preallocate elements
		//Modified for curves
		let positions = polygonSoup["v"];
		let indices = polygonSoup["f"];
		this.preallocateElements(positions, indices);

		// create and insert vertices
		//Do I want an extra halfedge between start/end vertices?
		let indexToVertex = new Map();
		for (let i = 0; i < positions.length; i++) {
			let v = new Vertex();
			this.vertices[i] = v;
			indexToVertex.set(i, v);
			let e = new Edge();
			this.edges[i] = e;
			let h = new Halfedge();
			h.vertex = v;
			v.halfedge = h;
			this.halfedges[i] = h;
			e.halfedge = h;
			h.edge = e;
			let hp = new Halfedge();
			this.halfedges[positions.length+i] = hp;
			h.twin = hp;
			hp.twin = h;
			hp.edge = e;
		}

		//Let's add dummy halfedges at endpoints for ease of computation
		//But the algorithm must operate as if they do not exist otherwise

		let h = this.halfedges[0];
		h.next = this.halfedges[1];
		h.prev = this.halfedges[positions.length-1];
		h = this.halfedges[positions.length];
		h.prev = this.halfedges[positions.length+1];
		h.next = this.halfedges[2*positions.length-1];
		for(let i = 1; i < positions.length-1;i++){
			h = this.halfedges[i];
			h.next = this.halfedges[i+1];
			h.prev = this.halfedges[i-1];
			h = this.halfedges[positions.length+i];
			h.next = this.halfedges[positions.length+i-1];
			h.prev = this.halfedges[positions.length+i+1];
		}
		h = this.halfedges[positions.length -1];
		h.prev = this.halfedges[positions.length - 2];
		h.next = this.halfedges[0];
		h = this.halfedges[2*positions.length-1];
		h.next = this.halfedges[2*positions.length-2];
		h.prev = this.halfedges[positions.length];


		//Need to modify halfedge construction process to rely on vertices

		// create and insert halfedges, edges and non boundary loop faces
		
		this.indexElements();

		return true;
	}

	/**
	 * Preallocates mesh elements.
	 * @private
	 * @method module:Core.Mesh#preallocateElements
	 * @param {module:LinearAlgebra.Vector[]} positions The vertex positions of a polygon soup mesh.
	 * @param {number[]} indices The indices of a polygon soup mesh.
	 */
	preallocateElements(positions, indices) {
		let nBoundaryHalfedges = 0;
		let sortedEdges = new Map();
		for (let I = 0; I < indices.length; I += 3) {
			for (let J = 0; J < 3; J++) {
				let K = (J + 1) % 3;
				let i = indices[I + J];
				let j = indices[I + K];

				// swap if i > j
				if (i > j) j = [i, i = j][0];

				let value = [i, j]
				let key = value.toString();
				if (sortedEdges.has(key)) {
					nBoundaryHalfedges--;

				} else {
					sortedEdges.set(key, value);
					nBoundaryHalfedges++;
				}
			}
		}

		let nVertices = positions.length;
		let nEdges = positions.length;
		let nFaces = 0;
		let nHalfedges = 2 * nEdges;
		let nInteriorHalfedges = 0;

		// clear arrays
		this.vertices.length = 0;
		this.edges.length = 0;
		this.faces.length = 0;
		this.halfedges.length = 0;
		this.corners.length = 0;
		this.boundaries.length = 0;
		this.generators.length = 0;

		// allocate space
		this.vertices = new Array(nVertices);
		this.edges = new Array(nEdges);
		this.faces = new Array(nFaces);
		this.halfedges = new Array(nHalfedges);
		this.corners = new Array(nInteriorHalfedges);
	}

	/**
	 * Checks whether this mesh has isolated vertices.
	 * @private
	 * @method module:Core.Mesh#hasIsolatedVertices
	 * @returns {boolean}
	 */
	hasIsolatedVertices() {
		for (let v of this.vertices) {
			if (v.isIsolated()) {
				alert("Mesh has isolated vertices!");
				return true;
			}
		}

		return false;
	}

	/**
	 * Checks whether this mesh has isolated faces.
	 * @private
	 * @method module:Core.Mesh#hasIsolatedFaces
	 * @returns {boolean}
	 */
	hasIsolatedFaces() {
		for (let f of this.faces) {
			let boundaryEdges = 0;
			for (let h of f.adjacentHalfedges()) {
				if (h.twin.onBoundary) boundaryEdges++;
			}

			if (boundaryEdges === 3) {
				alert("Mesh has isolated faces!");
				return true;
			}
		}

		return false;
	}

	/**
	 * Checks whether this mesh has non-manifold vertices.
	 * @private
	 * @method module:Core.Mesh#hasNonManifoldVertices
	 * @returns {boolean}
	 */
	hasNonManifoldVertices() {
		let adjacentFaces = new Map();
		for (let v of this.vertices) {
			adjacentFaces.set(v, 0);
		}

		for (let f of this.faces) {
			for (let v of f.adjacentVertices()) {
				adjacentFaces.set(v, adjacentFaces.get(v) + 1);
			}
		}

		for (let b of this.boundaries) {
			for (let v of b.adjacentVertices()) {
				adjacentFaces.set(v, adjacentFaces.get(v) + 1);
			}
		}

		for (let v of this.vertices) {
			if (adjacentFaces.get(v) !== v.degree()) {
				return true;
			}
		}

		return false;
	}

	/**
	 * Assigns indices to this mesh's elements.
	 * @private
	 * @method module:Core.Mesh#indexElements
	 */
	indexElements() {
		let index = 0;
		for (let v of this.vertices) {
			v.index = index++;
		}

		index = 0;
		for (let e of this.edges) {
			e.index = index++;
		}

		index = 0;
		for (let f of this.faces) {
			f.index = index++;
		}

		index = 0;
		for (let h of this.halfedges) {
			h.index = index++;
		}

		index = 0;
		for (let c of this.corners) {
			c.index = index++;
		}

		index = 0;
		for (let b of this.boundaries) {
			b.index = index++;
		}
	}
}


/**
 * Assigns an index to each element in elementList. Indices can be accessed by using
 * elements as keys in the returned dictionary.
 * @global
 * @function module:Core.indexElements
 * @param {Object[]} elementList An array of any one of the following mesh elements -
 * vertices, edges, faces, corners, halfedges, boundaries.
 * @returns {Object} A dictionary mapping each element in elementList to a unique index
 * between 0 and |elementList|-1.
 * @example
 * let vertexIndex = indexElements(mesh.vertices);
 * let v = mesh.vertices[0];
 * let i = vertexIndex[v];
 * console.log(i); // prints 0
 */
function indexElements(elementList) {
	let i = 0;
	let index = {};
	for (let element of elementList) {
		index[element] = i++;
	}

	return index;
}
