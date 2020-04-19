"use strict";

describe("NonisometricCurvatureFlow", function() {
	let steps, h;

	describe("integrate", function() {
		it("performs NonisometricCurvatureFlow with a fixed timestep", function() {

			let polygonSoup = MeshIO.readOBJ("cycle.obj");
			let mesh = new Mesh();
			mesh.build(polygonSoup);
			let geometry = new Geometry(mesh, polygonSoup["v"], false);

			let nonisometricCurvatureFlow = new NonisometricCurvatureFlow(geometry);

			

			let success = true;
			for (let i = 0; i < geometry.positions.length; i++) {
				if (!geometry.positions[i].isValid()) {
					success = false;
					break;
				}
			}

			chai.assert.strictEqual(success, true);
			memoryManager.deleteExcept([]);
		});
	});
});
