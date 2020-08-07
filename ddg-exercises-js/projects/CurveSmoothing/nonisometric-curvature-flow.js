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
		this.g = 1;
		console.log("Size");
		console.log(this.n);
		console.log(this.g);

		this.build_coordinates();
		this.compute_cycles();
	}

	compute_cycles(){
		let curvature = 0;
		for(let i = 0; i < this.n; i++){
			curvature = this.curvatures[i]+curvature;
		}
		this.cycles = curvature;
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

	//I think a counterclockwise rotation
	rotate(vec,curv){
		let x = Math.cos(curv)*vec.x - Math.sin(curv)*vec.y;
		let y = Math.sin(curv)*vec.x + Math.cos(curv)*vec.y;
		let ret = new Vector(x,y,0);
		return ret;
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

	compute_wilmore(){

		this.l_grad_energies = new Array(this.n);
		this.k_grad_energies = new Array(this.n);
		for(let i = 0; i < this.n; i++){
			//Our formula will be reindexed with l_i + l_(i+1)
			//because of how mesh initialized
			this.k_grad_energies[i] = 4*this.curvatures[i]/(this.lengths[i]+this.lengths[(i+1)%this.n]);
			let l_grad1 = -2*(this.curvatures[(i-1+this.n)%this.n]**2)/((this.lengths[(i-1+this.n)%this.n]+this.lengths[i])**2);
			let l_grad2 = -2*(this.curvatures[i]**2)/((this.lengths[i]+this.lengths[(i+1)%this.n])**2);
			this.l_grad_energies[i] = l_grad1+l_grad2;
			//console.log(this.l_grad_energies[i]);
		}
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

	build_coordinates(){

		let bound = this.geometry.mesh.vertices[0].halfedge;
		this.edges = new Array(this.n);
		this.curvatures = new Array(this.n);
		this.lengths = new Array(this.n);

		//Check to make sure not on degenerate halfedge(only there for graphical reasons)
		if(this.geometry.length(bound.edge) == 0){
			bound = bound.twin;
			bound = bound.next;
		}
	

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
				vec2 = this.geometry.vector(this.edges[0]);
			}
			else{
				vec2 = this.geometry.vector(this.edges[i+1]);
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
		
			//Checking sign of curvature, why do I have to do this?
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

		this.old_lengths = this.lengths.slice(0);
		this.old_curvatures = this.curvatures.slice(0);
	}

	compute_constraints(){

		let zero = new Vector(1,0,0);
		let first = this.geometry.vector(this.edges[0]);
		let theta = this.angle(zero,first);
		let theta_const = theta;
		//???Really shouldn't be using edges(unmutated by gradient)
		let sum = 0;
		let c1 = [];
		let c2 = [];
		let c3 = [];

		//Letting last curvature control first edge
		theta = this.curvatures[this.n-1];

		//Constraints don't seem to care about last curvature; maybe last curvature should control first component of sum
		//Computing cumulative angles: kth angle determines kth tangent
		let angles = [];
		angles.push(theta);
		for(let k = 1; k < this.n; k++){
			angles.push(angles[k-1]+this.curvatures[k-1]);
		}

		//First length constraints
		for(let i = 0; i < this.n; i++){
			let vec = this.geometry.vector(this.edges[i]);
			c1.push(vec.x);
		}
		//First componenent curvature constraints
		sum =0;
		let c1tmp = [];
		//Last curvature(does not determine constraints?)
		c1tmp.unshift(0);
		for(let k = this.n-1; k >= 1; k--){
			//k is kth curvature
			sum = sum-this.lengths[k]*Math.sin(angles[k]);
			c1tmp.unshift(sum);
		}
		c1 = c1.concat(c1tmp);

		//Second length constraints
		for(let i = 0; i < this.n; i++){
			let vec = this.geometry.vector(this.edges[i]);
			c2.push(vec.y);
		}
		//Second curvature constraints
		sum =0;
		let c2tmp = [];
		//Last curvature(does not determine constraints?)
		c2tmp.unshift(0);
		for(let k = this.n-1; k >= 1; k--){
			//k is kth curvature
			sum = sum+this.lengths[k]*Math.cos(angles[k]);
			c2tmp.unshift(sum);
		}
		c2 = c2.concat(c2tmp);
		
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
	}

	update(h){
		for(let i = 0; i < this.n; i++){
			this.lengths[i] = this.lengths[i] - h*this.g*this.l_grad_energies[i];
			this.curvatures[i] = this.curvatures[i] - h*this.k_grad_energies[i];
		}
	}

	reset(){
		this.lengths = this.old_lengths.slice(0);
		this.curvatures = this.old_curvatures.slice(0);
	}

	reconstruct(){
		this.geometry.positions[0].x = this.geometry.positions[0].x;
		this.geometry.positions[0].y = this.geometry.positions[0].y;
		this.geometry.positions[0].z = this.geometry.positions[0].z;

		this.tangents = [];
		//crv = this.curvatures[this.n-1]/2;
		//Im not convinced reconstruction is exactly right
		for(let i = 0; i < this.n; i++){
		
			if(i == 0){
				//Don't include next curvature
				let vec = new Vector(0,0,0);
				//vec.x = Math.cos(crv);
				//vec.y = Math.sin(crv);
				vec.scaleBy(this.lengths[i]);

				//New
				let turn = this.curvatures[this.n-1];
				let dir = this.geometry.vector(this.edges[this.n-1]);
				dir.normalize();
				vec.x = Math.cos(turn)*dir.x - Math.sin(turn)*dir.y;
				vec.y = Math.sin(turn)*dir.x + Math.cos(turn)*dir.y;
				vec.scaleBy(this.lengths[i]);
				//this.edges[i] = vec;
				this.tangents.push(vec);	
			}
			//Last curvature only arises because of closure condition
			else{
				let crvi = this.curvatures[i-1];
				//crv = (crvi + crv);

				//Not building on self
				//To build on self use tangents
				let vec = this.rotate(this.tangents[i-1],crvi);
				/*let vec = new Vector(0,0,0);
				//New
				vec.x = Math.cos(crvi)*this.edges[i-1].x - Math.sin(crvi)*this.edges[i-1].y;
				vec.y = Math.sin(crvi)*this.edges[i-1].x + Math.cos(crvi)*this.edges[i-1].y;*/
				vec.normalize();
				vec.scaleBy(this.lengths[i]);

				this.tangents.push(vec);
		
			}
			//this.print(this.vert_ordering[i+1]);
			this.geometry.positions[i+1].x = this.geometry.positions[i].x+this.tangents[i].x;
			this.geometry.positions[i+1].y = this.geometry.positions[i].y+this.tangents[i].y;
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
		if(Math.abs(total_curv-(this.cycles)) > 10e-2){
			console.log(Math.abs(total_curv-this.cycles));
			chai.assert.fail("Curvature not constant");
		}

		for(let i = 0; i < this.n; i++){
			chai.assert(this.lengths[i] >= 0);
		}

		//For closedness
		let start = this.geometry.positions[0];
		let end = this.geometry.positions[this.n];
		if(start.x != end.x && start.y != end.y){
			chai.assert.fail("Curve not closed");
		}

	}

	measure(type){
		let energy = this.compute_energy(type);
		

		let edge_length_sum = 0;
		let curvature_sum = 0;
		for(let i =0; i< this.n; i++){
			curvature_sum = curvature_sum + this.curvatures[i];
			edge_length_sum = edge_length_sum + this.lengths[i];
		}
		let mean = edge_length_sum/this.n;
		let mean_curv = curvature_sum/this.n;

		let variance = 0;
		let curv_var = 0;
		for(let i =0;i < this.n;i++){
			variance = variance + (this.lengths[i]-mean)**2;
			curv_var = curv_var + (this.curvatures[i]-mean_curv)**2;
		}
		variance = variance/(this.n-1);
		curv_var = curv_var/(this.n-1);

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
			if(Math.abs(this.curvatures[i] - mean_curv) >= Math.sqrt(curv_var)*3){
				temp.push(i);
				temp.push(this.curvatures[i]);
				temp.push(this.k_grad_energies[i]);
				irregular_curvatures.push(temp);
			}
		}

		console.log(type);
		console.log("Energy");
		console.log(energy);
		console.log("Length Sum");
		console.log(edge_length_sum);
		//console.log("Edge Mean");
		//console.log(mean);
		//console.log("Edge Variance");
		//console.log(variance);
		console.log("Irregular Edges");
		console.log(JSON.stringify(irregular_edges));
		console.log("Irregular Curvatures");
		console.log(JSON.stringify(irregular_curvatures));
	}

	line_search(type){
		//Energy might not be concave
		let energy = this.compute_energy(type);
		let h=.04;
		let alpha = .1;
		let beta = .5;
		let new_energy = 0;
		let n1 = this.l_grad_energies.map((x)=>x);
		let n2 = this.k_grad_energies.map((x)=>x);
		let dir = n1.concat(n2);
		let dot = this.dot(dir,dir);
		let lin = -alpha*dot;
		let lin_appx = 0;
		do{
			h=h*beta;
			this.update(h);
			new_energy = this.compute_energy(type);
			lin_appx = lin*h+energy;
			this.reset();
			/*console.log("New energy");
			console.log(new_energy);
			console.log("Lin appx");
			console.log(lin_appx);*/
		}while(new_energy > lin_appx);
		console.log(h);
		return h;
	}

	integrate(h,g,type){
		
		this.build_coordinates();

		this.invariant_check();

		this.compute(type);

		this.compute_constraints();

		this.measure(type);
		
		h = this.line_search(type);

		this.update(h);

		this.reconstruct();
		//May need to do some normalization
		this.correct();

		//Center curve around origin
		normalize(this.geometry.positions,this.geometry.mesh.vertices,false);

		this.step_count = this.step_count + 1;

	}

	run(h,g,type){
		this.build_coordinates();
		this.energy = this.compute_energy(type);
		let energy = this.energy;
		do{
			this.energy = energy;
			this.integrate(h,g,type);
			energy = this.compute_energy(type);
		}while(energy <= this.energy && (this.energy-energy) > 10e-2);
	}

	run25(h,g,type){
		let cnt = 0;
		this.build_coordinates();
		this.energy = this.compute_energy(type);
		let energy = this.energy;
		do{
			this.energy = energy;
			this.integrate(h,g,type);
			energy = this.compute_energy(type);
			cnt++;
		}while(cnt < 25);
	}

}
