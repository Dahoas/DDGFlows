"use strict";

class NonisometricCurvatureFlow {
	//Not sure if this is the proper way to create mesh
	
	constructor(geometry){
		this.geometry = geometry;
		//console.log(JSON.stringify(this.geometry.positions));
		this.n = geometry.mesh.vertices.length;
	}

	//Postitive Curvature is counterclockwise rotation
	angle(vec1, vec2){
		let dot =  vec1.dot(vec2);
		let res = dot/(vec1.norm()*vec2.norm());
		//Temp fix?
		if(res > 1){
			res = 1;
		}
		else if(res < -1){
			res = -1;
		}
		let red = Math.acos(res);
		return red;
	}

	//For arrays
	norm(vec){
		let sum = 0;
		for(let i = 0; i < vec.length; i++){
			sum = sum + vec[i]*vec[i];
		}
		sum = Math.sqrt(sum);
		return sum;
	}

	//For arrays
	dot(vec1, vec2){
		let sum = 0; 
		for(let i = 0; i < vec1.length; i++){
			sum = sum+vec1[i]*vec2[i];
		}
		return sum;
	}

	//For arrays
	minus(vec1,vec2){
		for(let i = 0; i < vec1.length; i++){
			vec1[i] = vec1[i] - vec2[i];
		}
		return vec1;
	}

	build_coordinates(){
		//logging objects seem to cache without actually printing
		console.log(JSON.stringify(this.geometry.positions));
		//console.log(this.geometry.mesh.boundaries[0]);
		let bound = this.geometry.mesh.boundaries[0].halfedge;
		this.edges = new Array(this.n);
		this.curvatures = new Array(this.n);
		this.lengths = new Array(this.n);
		let v0 = bound.vertex.index;
		this.vert_ordering = [v0];
		bound = bound.next;
		while(v0 != bound.vertex){
			this.vert_ordering.push(bound.vertex.index);
			bound = bound.next;
		}

		console.log(v0);

		for(let i = 0; i < this.n; i++){
			//Iterate around boundary
			//First bound is halfedge coming out of first vertex
			this.edges[i] = bound;
			//Length is length of halfedge extending from
			//vertex i
			this.lengths[i] = this.geometry.length(bound.edge);
			bound = bound.next;
		}
		//Test parity of halfedge: Want to be going in counterclockwise direction

		//Set up so that each length/curvature corresponds to
		//halfedge of same index(halfedge points to corner with curvature ki)
		for(let i =0;i<this.n;i++){
			//this.lengths[i] = this.geometry.length(mesh.halfedges[i].edge);
			//Note that following is consistent with how we specified halfedges of vertices
			//console.log("Computing " + i);
			//Problem: Curvatures are not signed
			let vec1 = this.geometry.vector(this.edges[i]);
			let vec2 = this.geometry.vector(this.edges[i].next);
			let dot =  vec1.dot(vec2);
			let res = dot/(vec1.norm()*vec2.norm());
			//Temp fix?
			if(res > 1){
				res = 1;
			}
			else if(res < -1){
				res = -1;
			}
			let red = Math.acos(res);
			/*console.log(dot);
			console.log(vec1.norm());
			console.log(vec2.norm());
			console.log(res);
			console.log(red);*/
			//To account for sign, check the counterclockwise rotation
			let vec = new Vector(vec1.x,vec1.y,0);
			vec.normalize();
			let targ = new Vector(vec2.x,vec2.y,0);
			targ.normalize();
			let tempx = vec.x;
			let tempy = vec.y;
			//Switching rotation direction fixes parity flipping
			vec.x = Math.cos(red)*tempx - Math.sin(red)*tempy;
			vec.y = Math.sin(red)*tempx + Math.cos(red)*tempy;
			//Floating point issues...
			if(Math.abs(targ.x - vec.x) > 10e-3 || Math.abs(targ.y - vec.y) > 10e-3){
				red = -red;
				//console.log(targ);
				//console.log(vec);
			}
			this.curvatures[i] = red;
		}
		//console.log(JSON.stringify(this.geometry.positions));
		console.log(JSON.stringify(this.lengths));
		console.log(JSON.stringify(this.curvatures));
	}

	print(i){
		let x = this.geometry.positions[i].x;
		let y = this.geometry.positions[i].y;
		let z = this.geometry.positions[i].z;
		console.log("Set vertex " + i + " at " + x + " " + y + " " + z);
	}

	reconstruct(){
		//Reconstruction:
		//I think figure gets rotated because building in opposite direction
		console.log(JSON.stringify(this.vert_ordering));
		this.geometry.positions[this.vert_ordering[0]].x = this.geometry.positions[this.vert_ordering[0]].x;
		this.geometry.positions[this.vert_ordering[0]].y = this.geometry.positions[this.vert_ordering[0]].y;
		this.geometry.positions[this.vert_ordering[0]].z = this.geometry.positions[this.vert_ordering[0]].z;
		//this.print(0);
		//Compute initial curvature
		let init = new Vector(1,0,0);
		let first = this.geometry.vector(this.edges[0]);
		first.normalize();
		let res = init.dot(first);
		if(res > 1){
			res = 1;
		}
		else if(res < -1){
			res = -1;
		}
		let crv = Math.acos(res);
		init.x = Math.cos(crv)*1-Math.sin(crv)*0;
		init.y = Math.sin(crv)*1+Math.cos(crv)*0;
		if(Math.abs(init.x - first.x) > 10e-3 || Math.abs(init.y - first.y) > 10e-3){
			crv = -crv;
			//console.log(targ);
			//console.log(vec);
		}
		for(let i = 0; i < this.n-1; i++){
			//Determine location of i+1 index from i
			//Center vertex 0 at 0
			//Make initial angle theta 0, I think we can always do this
			//What is our orientation?
			if(i == 0){
				//Don't include next curvature
				let vec = new Vector(0,0,0);
				vec.x = Math.cos(crv);
				vec.y = Math.sin(crv);
				vec.scaleBy(this.lengths[i]);
				this.geometry.positions[this.vert_ordering[i+1]].x = this.geometry.positions[this.vert_ordering[i]].x + vec.x;
				this.geometry.positions[this.vert_ordering[i+1]].y = this.geometry.positions[this.vert_ordering[i]].y + vec.y;
				this.geometry.positions[this.vert_ordering[i+1]].z = 0;
			}
			else{
				//This seems wrong?
				//Last curvature is irrelevant since
				//forcing second vertex on x-axis
				//Second to last curvature irrelevant since fixing
				//First vertex at origin
				//location of vert2 determined by 0th curvature
				//Adjusting i+1st position
				//Made this negative to try to eliminate flipping
				let crvi = this.curvatures[i-1];
				crv = (crvi + crv);
				//Make sure rotating right way
				//let half = this.edges[i-1];
				//let vec = this.geometry.vector(half);
				//vec.normalize();
				//vec.scaleBy(this.lengths[i]);
				//Right direction?
				//console.log(crvi);
				//console.log(crv);
				let vec = new Vector(0,0,0);
				vec.x = Math.cos(crv);
				vec.y = Math.sin(crv);
				vec.scaleBy(this.lengths[i]);
				this.geometry.positions[this.vert_ordering[i+1]].x = this.geometry.positions[this.vert_ordering[i]].x + vec.x;
				this.geometry.positions[this.vert_ordering[i+1]].y = this.geometry.positions[this.vert_ordering[i]].y + vec.y;
				this.geometry.positions[this.vert_ordering[i+1]].z = 0;
			}
			this.print(this.vert_ordering[i+1]);
		}
		console.log(JSON.stringify(this.geometry.positions));
		//console.log(JSON.stringify(this.lengths));
		//console.log(JSON.stringify(this.curvatures));
	}

	compute(){

		this.l_grad_energies = new Array(this.n);
		this.k_grad_energies = new Array(this.n);
		for(let i = 0; i < this.n; i++){
			//Our formula will be reindexed with l_i + l_(i+1)
			//because of how mesh initialized
			this.k_grad_energies[i] = 4*this.curvatures[i]/(this.lengths[i]+this.lengths[(i+1)%this.n]);
			//Formula for edges will contain i-1 and i curvatures
			let l_grad1 = -2*this.curvatures[(i-1+this.n)%this.n]/((this.lengths[(i-1+this.n)%this.n]+this.lengths[i])**2);
			let l_grad2 = -2*this.curvatures[i]/((this.lengths[i]+this.lengths[(i+1)%this.n])**2);
			this.l_grad_energies[i] = l_grad1+l_grad2;
		}
		console.log(JSON.stringify(this.l_grad_energies));
		console.log(JSON.stringify(this.k_grad_energies));

		//Computing constraints
		//Need to make sure curve is centered at origin and 
		//first tangent is along x-axis
		//Or compute theta
		//Also need right orientation(this is determined by mesh input)
		let zero = new Vector(1,0,0);
		let first = this.geometry.vector(this.edges[0]);
		let theta = this.angle(zero,first);
		let theta_const = theta;
		let sum = 0;
		let c1 = [];
		let c2 = [];
		let c3 = [];
		//First length constraints
		for(let i = 0; i < this.n; i++){
			let vec = this.geometry.vector(this.edges[i]);
			c1.push(vec.x);
		}
		//First curvature constraints
		for(let i = 0; i < this.n; i++){
			//i==0 determined by theta?
			if(i == 0){
				c1.push(0);
			}
			else{
				let res = -1.0*Math.sin(theta);
				res = res*this.lengths[i];
				c1.push(res);
				theta = theta + this.curvatures[i];
			}
		}

		//Second length constraints
		theta = theta_const;
		for(let i = 0; i < this.n; i++){
			let vec = this.geometry.vector(this.edges[i]);
			c2.push(vec.y);
		}
		//Second curvature constraints
		for(let i = 0; i < this.n; i++){
			//i==0 determined by theta?
			if(i == 0){
				c2.push(0);
			}
			else{
				let res = Math.cos(theta);
				res = res*this.lengths[i];
				c2.push(res);
				theta = theta + this.curvatures[i];
			}
		}

		//Third length constraints
		for(let i = 0; i < this.n; i++){
			c3.push(0);
		}
		//Third curvature constraints
		for(let i = 0; i < this.n; i++){
			c3.push(1);
		}

		let n1 = this.norm(c1);
		let n2 = this.norm(c2);
		let n3 = this.norm(c3);
		c1 = c1.map(x => x/n1);
		c2 = c2.map(x => x/n2);
		c3 = c3.map(x => x/n3);

		let state = [];
		for(let i = 0; i < this.n; i++){
			state.push(this.l_grad_energies[i]);
		}
		for(let i = 0; i < this.n; i++){
			state.push(this.k_grad_energies[i]);
		}

		console.log(state);

		let p1 = this.dot(state,c1);
		let p2 = this.dot(state,c2);
		let p3 = this.dot(state,c3);
		let proj1 = c1.map(x => x*p1);
		let proj2 = c2.map(x => x*p2);
		let proj3 = c3.map(x => x*p3);

		state = this.minus(state,proj1);
		state = this.minus(state,proj2);
		state = this.minus(state,proj3);

		for(let i = 0; i < this.n; i++){
			this.l_grad_energies[i] = state[i];
			this.k_grad_energies[i] = state[i+this.n];
		}
	}

	update(h){
		h = 0.02;
		//Computing updates
		for(let i = 0; i < this.n; i++){
			//Recall that we want to go in the opposite direction of 
			//our gradient
			//Make sure sign is correct
			this.lengths[i] = this.lengths[i] - h*this.l_grad_energies[i];
			this.curvatures[i] = this.curvatures[i] - h*this.k_grad_energies[i];
		}
	}

	integrate(h){
		//Why on earth is console.log nonlinear

		//Reconstruct and build_coords do not work correctly
		
		this.build_coordinates();

		this.compute();
		
		this.update(h);

		this.reconstruct();
		//May need to do some normalization

		//Then need normalize curve(via least squares), address disc. error
		//Pick vertex positions and then minimize diff
		//in edge lengths between actual and result

	}
}
