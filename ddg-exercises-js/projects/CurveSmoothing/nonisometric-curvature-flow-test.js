"use strict";

class NonisometricCurvatureFlow {
	//Not sure if this is the proper way to create mesh

	/*Think about using line slide test instead of rotatitins to test direction of curvature*/
	
	constructor(geometry){
		this.geometry = geometry;
		//console.log(JSON.stringify(this.geometry.positions));
		//Subtract 1 to eliminate imaginary last edge
		this.n = geometry.mesh.vertices.length-1;
		this.step_count = 0;
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

	normalize(vec){
		let res = this.norm(vec);
		vec = vec.map(x => x/res);
		return vec;
	}

	//Projects vec1 onto vec2
	//Assumes vec2 is unit vector
	proj(vec1,vec2){
		let res = this.dot(vec1,vec2);
		let proj = vec2.map(x => x*res);
		return proj;
	}

	dist(i,j){
		return Math.sqrt((this.geometry.positions[i].x-this.geometry.positions[j].x)**2 + (this.geometry.positions[i].y-this.geometry.positions[j].y)**2 + (this.geometry.positions[i].z-this.geometry.positions[j].z)**2);
	}

	compute_energy(type){
		let energy = 0;
		for(let i = 0; i < this.n; i++){
			let term = 0;
			//Maybe should change so that curvature of/at a vertex i is not the ith-1 curvature
			if(type == "wilmore"){
				term = 2*(this.curvatures[(i-1+this.n)%this.n]**2)*(1/(this.lengths[(i-1+this.n)%this.n]+this.lengths[i]));
			}
			else if(type == "squared_lengths"){
				term = 2*(this.curvatures[(i-1+this.n)%this.n]**2)*(1/(this.lengths[(i-1+this.n)%this.n]+this.lengths[i])) + this.lengths[i]**2;
			}
			else if(type == "squared_curvatures"){
				term = (this.curvatures[(i-1+this.n)%this.n]**2);
			}
			energy = energy + term;
		}
		return energy;
	}

	build_coordinates(){

		let bound = this.geometry.mesh.vertices[0].halfedge;
		this.edges = new Array(this.n);
		this.curvatures = new Array(this.n);
		this.lengths = new Array(this.n);
	

		for(let i = 0; i < this.n; i++){
			//Iterate around boundary
			//First bound is halfedge coming out of first vertex
			this.edges[i] = bound;
			//Length is length of halfedge extending from
			//vertex i
			this.lengths[i] = this.geometry.length(bound.edge);
			bound = bound.next;
		}

		for(let i = 0;i<this.n;i++){
		
			let vec1 = this.geometry.vector(this.edges[i]);
			let vec2 = undefined;
			if(i == this.n-1){
				vec2 = this.geometry.vector(this.edges[i].next.next);
			}
			else{
				vec2 = this.geometry.vector(this.edges[i].next);
			}
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
		//console.log(JSON.stringify(this.lengths));
		//console.log(JSON.stringify(this.curvatures));
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
		//console.log(JSON.stringify(this.vert_ordering));
		this.geometry.positions[0].x = this.geometry.positions[0].x;
		this.geometry.positions[0].y = this.geometry.positions[0].y;
		this.geometry.positions[0].z = this.geometry.positions[0].z;
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

		this.tangents = [];
		for(let i = 0; i < this.n; i++){
		
			if(i == 0){
				//Don't include next curvature
				let vec = new Vector(0,0,0);
				vec.x = Math.cos(crv);
				vec.y = Math.sin(crv);
				vec.scaleBy(this.lengths[i]);
				this.tangents.push(vec);

				this.geometry.positions[i+1].x = this.geometry.positions[i].x + vec.x;
				this.geometry.positions[i+1].y = this.geometry.positions[i].y + vec.y;
				this.geometry.positions[i+1].z = 0;
			}
			//Last curvature only arises because of closure condition
			else{

				let crvi = this.curvatures[i-1];
				crv = (crvi + crv);
			
				let vec = new Vector(0,0,0);
				vec.x = Math.cos(crv);
				vec.y = Math.sin(crv);
				vec.scaleBy(this.lengths[i]);
				this.tangents.push(vec);

				this.geometry.positions[i+1].x = this.geometry.positions[i].x + vec.x;
				this.geometry.positions[i+1].y = this.geometry.positions[i].y + vec.y;
				this.geometry.positions[i+1].z = 0;
		
			}
			//this.print(this.vert_ordering[i+1]);
		}

	

		//console.log(JSON.stringify(this.geometry.positions));
		//console.log(JSON.stringify(this.lengths));
		//console.log(JSON.stringify(this.curvatures));
	}

	compute_wilmore(){

		this.l_grad_energies = new Array(this.n);
		this.k_grad_energies = new Array(this.n);
		for(let i = 0; i < this.n; i++){
			//Our formula will be reindexed with l_i + l_(i+1)
			//because of how mesh initialized
			this.k_grad_energies[i] = 4*this.curvatures[i]/(this.lengths[i]+this.lengths[(i+1)%this.n]);
			/*console.log(this.k_grad_energies[i]);
			console.log(this.curvatures[i]);
			console.log(this.lengths[i]);
			console.log(this.lengths[(i+1)%this.n]);
			*///Formula for edges will contain i-1 and i curvatures
			let l_grad1 = -2*(this.curvatures[(i-1+this.n)%this.n]**2)/((this.lengths[(i-1+this.n)%this.n]+this.lengths[i])**2);
			let l_grad2 = -2*(this.curvatures[i]**2)/((this.lengths[i]+this.lengths[(i+1)%this.n])**2);
			this.l_grad_energies[i] = l_grad1+l_grad2;
			//console.log(this.l_grad_energies[i]);
		}
		//console.log(JSON.stringify(this.l_grad_energies));
		//console.log(JSON.stringify(this.k_grad_energies));	
	}

	compute_squared_lengths(){

		this.l_grad_energies = new Array(this.n);
		this.k_grad_energies = new Array(this.n);
		for(let i = 0; i < this.n; i++){
			//Our formula will be reindexed with l_i + l_(i+1)
			//because of how mesh initialized
			this.k_grad_energies[i] = 4*this.curvatures[i]/(this.lengths[i]+this.lengths[(i+1)%this.n]);
			//Formula for edges will contain i-1 and i curvatures
			let l_grad1 = -2*(this.curvatures[(i-1+this.n)%this.n]**2)/((this.lengths[(i-1+this.n)%this.n]+this.lengths[i])**2);
			let l_grad2 = -2*(this.curvatures[i]**2)/((this.lengths[i]+this.lengths[(i+1)%this.n])**2);
			//Adding sum of squared lengths
			let l_grad3 = 2*this.lengths[i];
			this.l_grad_energies[i] = l_grad1+l_grad2+l_grad3;
		}		
	}

	compute_squared_curvatures(){
		this.l_grad_energies = new Array(this.n);
		this.k_grad_energies = new Array(this.n);
		for(let i = 0; i < this.n; i++){
			//Our formula will be reindexed with l_i + l_(i+1)
			//because of how mesh initialized
			this.k_grad_energies[i] = 2*this.curvatures[i];
			this.l_grad_energies[i] = 0;
		}		
	}

	compute(type){
		if(type == "wilmore"){
			this.compute_wilmore();
		}
		else if(type == "squared_lengths"){
			this.compute_squared_lengths();
		}
		else if(type == "squared_curvatures"){
			this.compute_squared_curvatures();
		}
	}

	compute_constraints(){

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

		let state = [];
		for(let i = 0; i < this.n; i++){
			state.push(this.l_grad_energies[i]);
		}
		for(let i = 0; i < this.n; i++){
			state.push(this.k_grad_energies[i]);
		}

		c1 = this.normalize(c1.slice());
		let c2Onc1 = this.proj(c2.slice(),c1.slice());
		c2 = this.minus(c2.slice(),c2Onc1.slice());
		c2 = this.normalize(c2);
		let c3onc1 = this.proj(c3.slice(),c1.slice());
		c3 = this.minus(c3.slice(),c3onc1.slice());
		let c3onc2 = this.proj(c3.slice(),c2.slice());
		c3 = this.minus(c3.slice(),c3onc2.slice());
		c3 = this.normalize(c3.slice());

		let proj1 = this.proj(state.slice(),c1);
		let proj2 = this.proj(state.slice(),c2);
		let proj3 = this.proj(state.slice(),c3);

		state = this.minus(state.slice(),proj1.slice());
		state = this.minus(state.slice(),proj2.slice());
		state = this.minus(state.slice(),proj3.slice());

		for(let i = 0; i < this.n; i++){
			this.l_grad_energies[i] = state[i];
			this.k_grad_energies[i] = state[i+this.n];
		}
		//console.log(this.l_grad_energies);
		//console.log(this.k_grad_energies);
		//console.log(JSON.stringify(this.l_grad_energies));
		//console.log(JSON.stringify(this.k_grad_energies));	
	}

	update(h,g){
		//Note: Attains convergence for fixed curvature
		/*Fucking up somehow because one edge getting too large(seems to be the last one?)
			Perhaps I have implemented constraints correctly: Last edge is not fixed length so
			we seek to balance flattening of curvature(a straight line) against its blowup
		*/
		/*Need to add two endpoints for start/end to ensure use of last curvature condition
		*/
		/*How much of cycle breakage is discretization error vs
		incorrect constraints?
			Constraints do seem to help with breakage*/
			/*Edge lengths are made large as a result of trying to
			maximize energy*/
		//Computing updates
		//h=0.02;
		for(let i = 0; i < this.n; i++){
			//Recall that we want to go in the opposite direction of 
			//our gradient
			//Make sure sign is correct
			//Note: not changing edge lengths should correspond to wilmore flow
			this.lengths[i] = this.lengths[i] - h*g*this.l_grad_energies[i];
			this.curvatures[i] = this.curvatures[i] - h*this.k_grad_energies[i];
			/*if(i == this.n-1){
				console.log(this.curvatures[i]);
				console.log(this.k_grad_energies[i]);
			}
			else if(i == 0){
				console.log(this.curvatures[i]);
				console.log(this.k_grad_energies[i]);
			}*/
		}
	}

	correct(){

		let offset = 10e-8;
		let T = new Triplet(this.n,this.n);
		for(let i = 0; i < this.n; i++){
			for(let j = 0; j < this.n; j++){
				let res = 0;
				if(i == j){
					if(i == 0){
						res = 1/this.lengths[this.n-1] + 1/this.lengths[0]+offset;
					}
					else{
						res = 1/this.lengths[i-1] + 1/this.lengths[i]+offset;
					}
				}
				else if(i+1 == j){
					res = -1/this.lengths[i];
				}
				else if(i == this.n-1 && j == 0){
					res = -1/this.lengths[i];
				}
				else if(i == j+1){
					res = -1/this.lengths[j];
				}
				else if(j == this.n-1 && i == 0){
					res = -1/this.lengths[j];
				}
				T.addEntry(res,i,j);
			}
		}

		let L = SparseMatrix.fromTriplet(T);

		let bx = DenseMatrix.zeros(this.n,1);
		for(let i = 0; i < this.n; i++){
			let res = 0;
			if(i == 0){
				res = this.tangents[this.n-1].x/this.lengths[this.n-1] - this.tangents[0].x/this.lengths[0];
			}
			else{
				res = this.tangents[i-1].x/this.lengths[i-1]-this.tangents[i].x/this.lengths[i];
			}
			bx.set(res,i,0);
		}

		let by = DenseMatrix.zeros(this.n,1);
		for(let i = 0; i < this.n; i++){
			let res = 0;
			if(i == 0){
				res = this.tangents[this.n-1].y/this.lengths[this.n-1] - this.tangents[0].y/this.lengths[0];
			}
			else{
				res = this.tangents[i-1].y/this.lengths[i-1]-this.tangents[i].y/this.lengths[i];
			}
			by.set(res,i,0);
		}

		let lu = L.lu();
		//Error caused by next line
		//debugger;
		let xs = lu.solveSquare(bx);
		let ys = lu.solveSquare(by);

		for(let i = 0; i< this.n+1;i++){
			if(i == this.n){
				let x = xs.get(0,0);
				let y = ys.get(0,0);
				this.geometry.positions[i].x = x;
				this.geometry.positions[i].y = y;
			}
			else{
				let x = xs.get(i,0);
				let y = ys.get(i,0);
				this.geometry.positions[i].x = x;
				this.geometry.positions[i].y = y;
			}
		}
		//console.log(JSON.stringify(this.geometry.positions));
	}

	invariant_check(){
		let total_curv = 0;
		for(let i = 0; i < this.n; i++){
			total_curv = total_curv+this.curvatures[i];
		}
		//Only considering 2pi curvature curves
		if(Math.abs(Math.abs(total_curv)-(2*Math.PI)) > 10e-2){
			console.log(Math.abs(Math.abs(total_curv)-(2*Math.PI)));
			chai.assert.fail("Curvature not 2 pi");
		}

		for(let i = 0; i < this.n; i++){
			chai.assert(this.lengths[i] >= 0);
		}
	}

	measure(type){
		let energy = this.compute_energy(type);
		

		let edge_length_sum = 0;
		for(let i =0; i< this.n; i++){
			edge_length_sum = edge_length_sum + this.lengths[i];
		}
		let mean = edge_length_sum/this.n;

		let variance = 0;
		for(let i =0;i < this.n;i++){
			variance = variance + (this.lengths[i]-mean)**2;
		}
		variance = variance/(this.n-1);

		let irregular_edges = [];
		for(let i = 0; i < this.n; i++){
			let temp = [];
			if(Math.abs(this.lengths[i] - mean) >= Math.sqrt(variance)*3){
				temp.push(i);
				temp.push(this.lengths[i]);
				temp.push(this.l_grad_energies[i]);
				irregular_edges.push(temp);
			}
		}

		let irregular_curvatures = [];
		for(let i = 0; i < this.n; i++){
			let temp = [];
			if(Math.abs(this.curvatures[i]) >= Math.PI/8){
				temp.push(i);
				temp.push(this.curvatures[i]);
				temp.push(this.k_grad_energies[i]);
			}
		}

		console.log(type);
		console.log("Energy");
		console.log(energy);
		console.log("Length Sum");
		console.log(edge_length_sum);
		console.log("Edge Mean");
		console.log(mean);
		console.log("Edge Variance");
		console.log(variance);
		console.log("Irregular Edges");
		console.log(JSON.stringify(irregular_edges));
		console.log("Irregular Curvatures");
		console.log(JSON.stringify(irregular_curvatures));
		console.log("Steps");
		console.log(this.step_count);
	}

	integrate(h,g,type){

		//Decide on Metrics:
		//Speed of convergence
		//Step size taken
		//Stability(dependent on step size)
		//edge length change
		//Aesthetic
		
		this.build_coordinates();

		this.invariant_check();

		this.compute(type);

		this.compute_constraints();

		this.measure(type);
		
		this.update(h,g);

		this.reconstruct();
		//May need to do some normalization
		//this.correct();

		//Center curve around origin
		normalize(this.geometry.positions,this.geometry.mesh.vertices,false);

		this.step_count = this.step_count +1;

	}
}
