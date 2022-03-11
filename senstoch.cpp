#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


/**************************************************
*
*    	statdist: compute column vector with top left eigenvector of input matrix
*
**************************************************/

// [[Rcpp::export]]
vec statdist(mat P){
    cx_vec eigval;
    cx_mat eigvecL;
    eig_gen(eigval,eigvecL,trans(P)); // Note: Left e-vectors are still given as columns
    vec abseigval=abs(eigval);
    uword maxind=abseigval.index_max();  // maxind is now the index of the
    // top eigenvalue
    vec topvec=conv_to<vec>::from(eigvecL.col(maxind));
    return topvec/sum(topvec);
}


/**************************************************
*
*		reverseMC: compute transitions for backward MC
*
**************************************************/



// [[Rcpp::export]]
mat reverseMC(mat P){
	vec nu=statdist(P);
	uword K=P.n_rows;
	mat Prev(K,K);
	for (uword i=0; i<K;i++){
		for (uword j=0; j<K;j++){
			Prev(i,j)=P(j,i)*nu[j]/nu[i];
		}
	}
	return Prev;	
}

/*
MATMAX=function to input matrix and a number, output same matrix with
 	each element maxed with the number
 	
*/

// [[Rcpp::export]]
mat matmax(mat M, double x){
	mat MM(M);
	double xx(x);
	uword K=MM.n_rows;
	uword J=MM.n_cols;
	for (uword k=0; k<K;k++){
		for (uword j=0; j<J;j++){
			if(MM(k,j)< xx){ MM(k,j)=xx;}
		}
	}
	return(MM);
}

/*
MATROUND=function to input matrix and a number n, output same matrix with
 	each element rounded to n places
 	
*/

// [[Rcpp::export]]
mat matround(mat M, double n){
	double x=pow(10,n);
	mat rm=round(M*x)/x;
	return(rm);
}

/*
ROWSUM1=function to input matrix, output same with rows normalised to sum to 1
 	
*/

// [[Rcpp::export]]
mat rowsum1(mat M){
	mat MM(M);
	double s;
	uword K=MM.n_rows;
	for (uword k=0; k<K;k++){
		s=sum(MM.row(k));
		if (s==0){s=1;}
		MM.row(k)=MM.row(k)/s;
		}
	return MM;
}

/*
Rowsums
    Input matrices as vector of all entries mats
        dxd matrices
        We're picking out the *whichmat* matrix
*/

// [[Rcpp::export]]
vec rowsums(vec mats, uword whichmat, uword d){
    vec init=zeros<vec>(d);
    uword index0=whichmat*d*d;
    for (uword i=0; i<d; i++){
        for (uword j=0; j<d; j++){
            init[i]+=mats[index0+d*j+i]; // sum along the i-th row
        }
    }
    return init;
}

/*
matcum=function to input matrix, output same with rows as cumulative sum
 	
*/

// [[Rcpp::export]]
mat matcum(mat M){
	uword J=M.n_rows;
	uword K=M.n_cols;
	mat MM(M.begin(), J,K,true);
	for (uword j=0; j<J; j++){
		for (uword k=1; k<K;k++){
			MM(j,k) = MM(j,k) + MM(j,k-1);
		}
	}
	return MM;
}

// [[Rcpp::export]]
double hilbert(vec vek1,vec vek2){
	double H;
	if (min(vek1)<=0||min(vek2)<=0){H=INFINITY;}
	else {vec lvek1=log(vek1);
		vec lvek2=log(vek2);
		H=max(lvek1-lvek2)+max(lvek2-lvek1);}
	return H;
	}
	
// [[Rcpp::export]]
uword firsttrue(uvec LV){
	uword i=0;
	uword L=LV.n_elem;
	while (i<L && LV[i]==0){i++;}
	return i;}	

 // [[Rcpp::export]]
uword furthestpoint(mat matpoints,vec v1,uvec whichpoints){
	mat M(matpoints);
	vec V(v1);
	uvec wp(whichpoints);
	int L=M.n_cols;
	vec distances(L);
	for (int i=0;i<L; i++){
		if(wp[i]==1){distances[i]=hilbert(V,M.col(i));}
		else{distances[i]=0;}
	}
	uword index=distances.index_max();
	return index;
	}
	

// Input matrix whose columns are vectors, a point, and a radius;
//    Output unsigned integer vector of which points are outside the ball

// [[Rcpp::export]]
uvec outsideball(mat veks, vec point, double radius, uvec whichpoints){
	// Initialize variables 
	mat vekvec=veks;
	uvec wp(whichpoints);
	uword L=wp.n_elem;
	//
	for(int i=0; i<L ; i++){
		if (wp[i]==1){
			if (hilbert(point,vekvec.col(i))<= radius){wp[i]=0;}
			}
		}
	return wp;
}

/********************************************************************************
 * 
 *  doublenormal: Auxiliary function to diam to find a double normal in a set of points;
 *          that is, a pair of points such that each is the farthest from the other
 *          as measured by the distance measure hilbert
 * 
 *  Input: matpoints= matrix whose column vectors are the set of points
 *          whichpoints= binary vector of length= number of points
 *              indicating which points are still under consideration
 *              1 means point is still being considered
 *          startpoint= index of point that we start with.
 *********************************************************************************/


// whichpoints is a binary vector, 
// [[Rcpp::export]]
uvec doublenormal(mat matpoints,uvec whichpoints,uword startpoint){
	//declare variables
	double maxdist(0),testdist(0);
	uvec wp(whichpoints);
	uvec testpoint(3);
	uvec wpreturn(wp.begin(),wp.n_elem); //clone so it won't change the input
	testpoint[0]=0;
	testpoint[1]=startpoint;
	wpreturn[testpoint[1]]=false;  // Remove the point just chosen.
	vec v2=matpoints.col(testpoint[1]);
	vec v1=v2;
	uword stop(0);
	testpoint[2]=furthestpoint(matpoints,v2,wpreturn);// Find a furthest point
		vec v3=matpoints.col(testpoint[2]);
		testdist=hilbert(v3,v2);
	wpreturn[testpoint[2]]=0;   // Remove the point just chosen.
	// At this moment, v3 is the farthest from v2
	do {
		v1=v2;
		testpoint[0]=testpoint[1];
		v2=v3;
		testpoint[1]=testpoint[2];// Now v2 is the farthest from v1.
		// We need to test if v1 is the farthest from v2.
		maxdist=testdist;
		testpoint[2]=furthestpoint(matpoints,v2,wpreturn);
		v3=matpoints.col(testpoint[2]);
		testdist=hilbert(v2,v3);
		if (testdist>maxdist){wpreturn[testpoint[2]]=0;}
			else {stop=1;}
		}
		while (stop==0); // When the new distance (to 3) isn't longer, points 1 and 2 are a DN.
	testpoint.resize(2);
	return join_cols(testpoint,wpreturn);
}



// [[Rcpp::export]]
double shortdiam(mat vek){
	// Initialize variables 
	uword npoints=vek.n_cols; //npoints= number of points
	double dmax=0, delta=1;
	for (uword i=0; i<npoints-1;i++){
		for (uword j=i+1; j<npoints;j++){
			delta=hilbert(vek.col(i),vek.col(j));
			if (delta>dmax){dmax=delta;}
		}
	}
	return dmax;
}


/* Input matrix whose columns are vectors, output
		diameter of set of vectors and indices of vectors at max distance */

// [[Rcpp::export]]
double diamlong(mat vek){
	// Initialize variables 
	int k=vek.n_rows;   //k= dimension of space
	int npoints=vek.n_cols; //npoints= number of points
	double dmax=0, delta=1, delta2=0;
	vec maxdistance(npoints);
	vec vekvek=vec(vek.begin(), k*npoints,false);
	vec midpoint(k);
	uvec dnsreturn(npoints+2);  // used to receive output from doublenormal routine
	uvec stillin=ones<uvec>(npoints);    // vector describes whether points are under consideration or not
			// start with all points in
	uvec Q=stillin;
	uword stop=0,sQ=0;
	umat storedns(npoints+2,0);
	mat vekmod=mat(vek.begin(),k,npoints);  //the matrix of points for internal consumption
//	IntegerVector Q=seq(0,npoints-1);
	uvec candidate(2);
		// place 0 and 1 say which points
		// as we choose more double normals, we add columns to storedns
	if (min(vekvek)<=0){
		uvec outlist(2);  // vector =(distance, first point, second point)
		candidate[0]=0;
		candidate[1]=0;
		dmax=1000000;
		storedns=join_cols(outlist,stillin);}
 else{ 
		while (sum(Q)>0 &&stop==0)
		{
			dnsreturn=doublenormal(vek,stillin,firsttrue(Q));
			candidate[0]=dnsreturn[0];
			candidate[1]=dnsreturn[1];
			// Calculate midpoint of double normal
				for (int i=0;i<k ; i++){
					midpoint[i]=sqrt(vekmod(i,candidate[0])*vekmod(i,candidate[1]));}
			delta=hilbert(vekmod.col(candidate[1]),vekmod.col(candidate[0]));
			stillin=dnsreturn.subvec(2,npoints+1);
			Q=outsideball(vekmod,midpoint,delta/2,stillin);
			storedns=join_rows(join_cols(candidate,Q),storedns);
			if (dmax >= delta){stop=1;}
				else {dmax=delta;}
		}
	if (stop==1) {
		Q=storedns.col(1).subvec(2,npoints+1);
		uvec fQ=find(Q==1);
		uvec fP=find(stillin==1);
		for (uword i=0; i<sQ; i++){
			for(uword j=0; j<sum(stillin); j++){
				delta2=hilbert(vekmod.col(fQ[i]),vekmod.col(fP[j]));
				if (delta2>delta){delta=delta2;}
				}
			}
		if (delta> dmax){dmax=delta;}
		}
	}
	return dmax;
}

// [[Rcpp::export]]
double diam(mat vek){
	uword npoints=vek.n_cols;
	double diameter;
	if (npoints < 15){diameter=shortdiam(vek);}
		else {diameter=diamlong(vek);}
	return diameter;
}

//     lmatprod
//
// Input matrices as the vector mx, and a sequence. Produce the product of successive
//     left multiplications.
// Matrices are assumed all dimmats x dimmats.
// sequence is a sequence of states, each corresponding to 

// [[Rcpp::export]]
cube lmatprod(vec mx, mat startmatrix,umat sequence,uword dimmats){
uword d(dimmats);
umat ss(sequence);
int slices = mx.n_elem/(d*d);
// Turn mx into cube 
	cube matrixcube(mx.begin(),d,d,slices,false);
//
int N=ss.n_rows;  // number of products to make
int n=ss.n_cols;  // length of each product
int outcols=startmatrix.n_cols;
cube uitmat=cube(d,outcols,N);
mat currentprod=mat(startmatrix);
for (int whichain=0; whichain<N;whichain++){
	currentprod=startmatrix;
	// now multiply 
	for (int time=0; time<n;time++){
		currentprod=matrixcube.slice(ss(whichain,time))*currentprod;
	}
	// put the result into the output cube 
/*	for(int i =0; i<K;i++){
		for(int j =0; j<outcols;j++){
			uitmat(i,j,whichain)=currentprod(i,j);
		}
	}*/
	uitmat.slice(whichain)=currentprod;
}
return uitmat;
}

// Input matrices as the vector mx, and a sequence. Produce the product of successive
//     RIGHT multiplications.
// Matrices are assumed all dimmats x dimmats.
// sequence is a sequence of states, each corresponding to 

// [[Rcpp::export]]
cube rmatprod(vec mx, mat startmatrix,umat sequence,uword dimmats){
// uword whichside(0); // =0 for left multiplication;
uword d(dimmats);
umat ss(sequence);
uword slices = mx.n_elem/(d*d);
// Turn mx into cube 
	cube matrixcube(mx.begin(),d,d,slices,false);
//
int N=ss.n_rows;  // number of products to make
int n=ss.n_cols;  // length of each product
int outrows=startmatrix.n_rows;
cube uitmat=cube(outrows,d,N);
mat currentprod=mat(startmatrix);
for (int whichain=0; whichain<N;whichain++){
	currentprod=startmatrix;
	// now multiply 
	for (int time=0; time<n;time++){
		currentprod=currentprod*matrixcube.slice(ss(whichain,time));
	}
	// put the result into the output cube 
/*	for(int i =0; i<d;i++){
		for(int j =0; j<outcols;j++){
			uitmat(i,j,whichain)=currentprod(i,j);
		}
	}*/
	uitmat.slice(whichain)=currentprod;
}
return uitmat;
}


//     lvecprod
//
// Input matrices as the vector mx, and a sequence. Produce the product of successive
//     left multiplications acting on the column vector initial.
// Matrices are assumed all dimmats x dimmats.
// sequence is a sequence of states, each corresponding to 

// [[Rcpp::export]]
mat lvecprod(vec mx, vec initial,umat sequence,uword dimmats){
uword d(dimmats);
umat ss(sequence);
int slices = mx.n_elem/(d*d);
// Turn mx into cube 
	cube matrixcube(mx.begin(),d,d,slices,false);
//
int N=ss.n_rows;  // number of products to make
int n=ss.n_cols;  // length of each product
mat uitmat(d,N);
vec currentprod=initial;
for (int whichain=0; whichain<N;whichain++){
	currentprod=initial;
	// now multiply 
	for (int time=0; time<n;time++){
		currentprod=matrixcube.slice(ss(whichain,time))*currentprod;
	}
	// put the result into the output cube 
/*	for(int i =0; i<K;i++){
		for(int j =0; j<outcols;j++){
			uitmat(i,j,whichain)=currentprod(i,j);
		}
	}*/
	uitmat.col(whichain)=currentprod;
}
return uitmat;
}

/*
		sampledist: Draw one sample from a discrete distribution
		
	Input variables:
		cumdist= cumulative distribution of states
		R = uniform random variable
*/
// [[Rcpp::export]]
uword sampledist(vec cumdist,double R){
	uword teststate=0;
	while(cumdist[teststate]<R) {
			teststate++;}
	return teststate;
}

/*
 sample_iid: Draw multiple samples from a discrete distribution

Input variables:
cumdist= 1 x K cumulative distribution of states
R = uniform random variable
n = number of samples
*/

// [[Rcpp::export]]
uvec sample_iid(vec cumdist, uword n){
    uvec test_samples= zeros<uvec>(n);
    vec Random = randu<vec>(n);
    for (uword i=0; i<n; i++){
        test_samples[i] = sampledist(cumdist,Random[i]);}
    return test_samples;
}

/*  	MC: Simulate realisations of Markov chain

	Input variables:
		P is the transition matrix
		K  is the number of states
		q is the length K vector giving the initial distribution
		numchains is the number of times the chain is run
		T is length of each chain
		rnd is a vector of N*(T+1) iid uniform random variables
		
	Output numstates x T+1 matrix of states (numbered 0..K-1)
*/


// [[Rcpp::export]]
umat MC(uword numchains,uword T,mat P,rowvec q, vec rnd){
//Initialise local variables
//	uword K=P.n_rows;
	uword N(numchains); 
    umat mc(N,T+1);  // storage for the output chains
	vec Random(rnd);
// turn qcum and P into cumulative probabilities
	rowvec qq(q);
	mat Pcum=matcum(P);
	rowvec qcum=matcum(qq);
	uword currentstate;  // Will be used in looping each chain to store the current state
	uword basecount;
	for (uword whichain=0; whichain<N; whichain++){
		basecount=whichain*T;
		currentstate=sampledist(trans(qcum),Random[basecount]);  //trans because sampledist takes column vectors
		mc(whichain,0)=currentstate;
	for (uword time=1; time<T+1; time++){
		currentstate=sampledist(trans(Pcum.row(currentstate)),Random[basecount+time]);
		mc(whichain,time)=currentstate;
}
    }
return mc;
}

/**********************************************************************
 * 
 *  inhomseq1: Input a sequence of transition matrices;
 *          Output simulated sequences from the inhomegeneous Markov chain they define
 * 
 ***********************************************************************/


// simulate n steps from an inhomogeneous MC with transitions P[[1]],P[[2]],...P[[T-1]]
// N sequences
// Output is NxT array

// [[Rcpp::export]]
umat inhomseq1(uword Nin,vec Pinput,rowvec qin, vec rnd,uvec timesin){
/*  Input variables:
		Pvec is a vector with length K^2*(T-1), giving the T different transition matrices
		numstates  is the number of states
		q is the length K vector giving the initial distribution
		numchains is the number of times the chain is run
		rnd is a vector of N*(T+1) iid uniform random variables
        times (timesin) = 2-place vector of starting and ending times
*/
//Initialise local variables
	vec Random(rnd), Pvec(Pinput);
	uvec times(timesin);
	uword N(Nin);
	rowvec q(qin);
    uword K=q.n_elem;
	uword T=times[1]-times[0]+1;  //length of each chain
	umat mc=zeros<umat>(N,T);  // storage for the output chains
	cube Pcube(Pvec.begin()+K*K*times[0], K,K,T,true); //slices of this cube are transition matrices
	rowvec qcum=matcum(q);
	uword currentstate;  // Will be used in looping each chain to store the current state
	uword basecount;
	double R;
// turn Pcube into cumulative probabilities
	for(uword t=0;t<T-1; t++){
		Pcube.slice(t)=matcum(Pcube.slice(t));
		Pcube.slice(t).col(K-1).fill(1);  // In case we got a row of all 0's,
								// put a 1 on the end; it shouldn't matter,
								// because any such row should correspond to q=0.
		}
		for (uword whichain=0; whichain<N; whichain++){
		basecount=whichain*T;
		uword teststate=0;
		while(qcum[teststate]<Random[basecount]) {
			teststate++;}
		mc(whichain,0)=teststate;
		currentstate=teststate;
	for (uword time=1; time<T; time++){
		teststate=0;
		R=Random[basecount+time];
		while(Pcube(currentstate,teststate,time-1)<R && teststate<K) {
			teststate++;}
		mc(whichain,time)=teststate;
		currentstate=teststate;
		}
	}
return mc;
}

/*
				SimulateU: Output samples from the stationary distribution
		Input variables:
			mx=matrix elements (as vector, in column normal form)
			P=transition matrix
			d= dimensions of the matrices
			N= number of samples
			T= number steps in the first round
			start= starting state (numbered from 0)
			maxerror= maximum allowable error
		
		Output:
			(K+1) x N matrix. Top K x N is a collection of N samples of U;
				bottom row is the diameters
*/

// [[Rcpp::export]]
mat SimulateU(vec matsin, mat Pin,uword din,uword Nin, uword Tin, uword startin,double maxerrin){
	vec mx(matsin);
	mat P(Pin);
	uword d(din),N(Nin),start(startin),T(Tin);
	double maxerror(maxerrin);
	uword slices = mx.n_elem/(d*d);
	uword K=P.n_rows;
	uword newstart;
	urowvec testsequence(T); // This will store one MC sequence
	mat storeU(d,N); // Store all the output
	mat init(d,d);
	mat revP=reverseMC(P);
	init.eye();  // init=Identity matrix for initialising rmatprod
	colvec testvec(d);
	rowvec storediam(N);
	cube matrixcube(mx.begin(),d,d,slices,false);
	rowvec q(K),qnew(K);
	q.zeros();
	q[start]=1;
	vec rnd=randu<vec>(N*(T+1)), morernd(T+1);
	umat sequences=MC(N,T,revP,q,rnd);
		sequences.shed_col(0);
	cube testproduct=rmatprod(mx,init,sequences,d);
	for (uword whichain=0;whichain<N; whichain++){
		storediam[whichain]=diam(testproduct.slice(whichain));
		while (storediam[whichain]> maxerror){
			testsequence=sequences.row(whichain);
			morernd=randu<vec>(T+1);
			newstart=testsequence[T-1];
			qnew.zeros();
			qnew[newstart]=1;
			testsequence=MC(1,T,revP,qnew,morernd);
			testsequence=testsequence.subvec(1,T);
			testproduct.slice(whichain)=rmatprod(mx,testproduct.slice(whichain),testsequence,d);
			storediam[whichain]=diam(testproduct.slice(whichain));
		}
		testvec=sum(testproduct.slice(whichain),1);
		storeU.col(whichain) = testvec/sum(testvec);
	} 
	return join_cols(storeU, storediam);
}


/*
				SimulateV: Output samples from the stationary distribution
		Input variables:
			mx=matrix elements (as vector, in column normal form)
			P=reverse transition matrix
			d= dimensions of the matrices
			N= number of samples
			T= number steps in the first round
			start= starting state (numbered from 0)
			maxerror= maximum allowable error
		
		Output:
			(K+1) x N matrix. Top K x N is a collection of N samples of U;
				bottom row is the diameters
*/

// [[Rcpp::export]]
mat SimulateV(vec matsin, mat Pin,uword din,uword Nin, uword Tin, uword startin,double maxerrin){
	vec mx(matsin);
	mat P(Pin);
	uword d(din),N(Nin),start(startin),T(Tin);
	double maxerror(maxerrin);
	uword slices = mx.n_elem/(d*d);
	uword K=P.n_elem;
	uword newstart;
	urowvec testsequence(T); // This will store one MC sequence
			// Note for later: This doesn't work if it's a row vector,
			// since matrix subsetting needs a matrix output
	mat storeV(N,d); // Store all the output
	mat init(d,d);
	init.eye();  // init=Identity matrix for initialising rmatprod
	rowvec testvec(d);
	vec storediam(N);
	cube matrixcube(mx.begin(),d,d,slices,false);
	rowvec q(K),qnew(K);
	q.zeros();
	q[start]=1;
	vec rnd=randu<vec>(N*(T+1)), morernd(T+1);
	umat sequences=MC(N,T,P,q,rnd);
		sequences.shed_col(0); // Get rid of starting state (which will be always the same)
	cube testproduct=lmatprod(mx,init,sequences,d);
	 for (uword whichain=0;whichain<N; whichain++){
		storediam[whichain]=diam(testproduct.slice(whichain).t());
		while (storediam[whichain]> maxerror){
			testsequence=sequences.row(whichain);
			morernd=randu<vec>(T+1);
			newstart=testsequence[T-1];
			qnew.zeros();
			qnew[newstart]=1;
			testsequence=MC(1,T,P,qnew,morernd);
			testsequence=testsequence.subvec(1,T);
			testproduct.slice(whichain)=lmatprod(mx,testproduct.slice(whichain),testsequence,d);
			storediam[whichain]=diam(testproduct.slice(whichain));
		}
		testvec=sum(testproduct.slice(whichain),0);
		storeV.row(whichain) = testvec/sum(testvec);
	}
	return join_rows(storeV, storediam);
}

/*
 SimulateV_iid: Output samples from the stationary distribution
    when sampling iid from a given distribution
Input variables:
mx=matrix elements (as vector, in column normal form)
nuin=distribution sampling from
d= dimensions of the matrices
N= number of samples
T= number steps in the first round
maxerror= maximum allowable error

Output:
(K+1) x N matrix. Top K x N is a collection of N samples of U;
bottom row is the diameters
*/

// [[Rcpp::export]]
mat SimulateV_iid(vec matsin, vec nuin,uword din,uword Nin, uword Tin,double maxerrin){
    vec mx(matsin);
    vec nu(nuin);
    uword d(din),N(Nin),T(Tin);
    double maxerror(maxerrin);
    uword slices = mx.n_elem/(d*d);
    urowvec testsequence(T); // This will store one iid sequence
    mat storeV(N,d); // Store all the output
    mat init(d,d);
    init.eye();  // init=Identity matrix for initialising rmatprod
    rowvec testvec(d);
    vec storediam(N);
    cube matrixcube(mx.begin(),d,d,slices,false);
    uvec sequences=sample_iid(cumsum(nu), N*T);
    umat sequencemat;
    sequencemat.insert_cols(0,sequences);
    sequencemat.reshape(N,T);
    cube testproduct=lmatprod(mx,init,sequencemat,d);
    for (uword whichain=0;whichain<N; whichain++){
        storediam[whichain]=diam(testproduct.slice(whichain).t());
        while (storediam[whichain]> maxerror){
            testsequence=sample_iid(cumsum(nu), T);
            testproduct.slice(whichain)=lmatprod(mx,testproduct.slice(whichain),testsequence,d);
            storediam[whichain]=diam(testproduct.slice(whichain));
        }
        testvec=sum(testproduct.slice(whichain),0);
        storeV.row(whichain) = testvec/sum(testvec);
    }
    return join_rows(storeV, storediam);
}

// [[Rcpp::export]]
mat SimulateU_iid(vec matsin, vec nuin,uword din,uword Nin, uword Tin,double maxerrin){
    vec mx(matsin);
    vec nu(nuin);
    uword d(din),N(Nin),T(Tin);
    double maxerror(maxerrin);
    uword slices = mx.n_elem/(d*d);
    uvec testsequence(T); // This will store one iid sequence
    mat storeU(d,N); // Store all the output
    mat init(d,d);
    init.eye();  // init=Identity matrix for initialising rmatprod
    vec testvec(d);
    rowvec storediam(N);
    cube matrixcube(mx.begin(),d,d,slices,false);
    uvec sequences=sample_iid(cumsum(nu), N*T);
    umat sequencemat;
    sequencemat.insert_cols(0,sequences);
    sequencemat.reshape(N,T);
    cube testproduct=lmatprod(mx,init,sequencemat,d);
    for (uword whichain=0;whichain<N; whichain++){
        storediam[whichain]=diam(testproduct.slice(whichain).t());
        while (storediam[whichain]> maxerror){
            testsequence=sample_iid(cumsum(nu), T);
            testproduct.slice(whichain)=rmatprod(mx,testproduct.slice(whichain),testsequence,d);
            storediam[whichain]=diam(testproduct.slice(whichain));
        }
        testvec=sum(testproduct.slice(whichain),1);
        storeU.col(whichain) = testvec/sum(testvec);
    }
    return join_cols(storeU, storediam);
}



/**********************************************************************

				SimulateY

************************************************************************/
// We need to give both K (number of states) and d (dimensions of the X matrices)
// so that the function can parse mx (the matrices) and transitions.
// Input variables:
//		matsin=vector of all matrix entries (for the products)
//		transin= transition matrix elements
//		couplein=state at which coupling occurs; starts at 0
//		Nin= number of chains
//		din= dimension of matrices
//		timesin= # time steps. So length of trans is K^2*(T-1)

// Product includes the coupling state, but not the starting state


// [[Rcpp::export]]
cube SimulateYmat(vec matsin, vec transin,uword Kin,uword din,uword Nin, uword couplein,uvec timesin){
	vec matrices(matsin),trans(transin);
	uword K(Kin),d(din),N(Nin),couple(couplein);
	uvec times(timesin);
		rowvec q=zeros<rowvec>(K);
		q[couple]=1;
		uword T=trans.n_elem/(K*K);
		vec rnd=randu<vec>(N*(T+1));
		umat chainruns=inhomseq1(N,trans,q,rnd,times);
		mat init(d,d);
		init.eye();
		cube Y=rmatprod(matrices,init,chainruns,d);
		return Y;
}

/*
			SimulateYsums
	
	Just like SimulateY, but we output only the Y*init vectors
    init= X_e * vector of all 1s

*/

// [[Rcpp::export]]
mat SimulateYsums(vec matsin, vec transin,uword Kin,uword din,uword Nin, uword couplein,uvec timesin,vec init){
	vec matrices(matsin),trans(transin);
	uword K(Kin),d(din),N(Nin),couple(couplein);
	uvec times(timesin);
		rowvec q=zeros<rowvec>(K);
		q[couple]=1;
		uword T=trans.n_elem/(K*K);
		vec rnd=randu<vec>(N*(T+1));
		umat chainruns=inhomseq1(N,trans,q,rnd,times);
		umat chainrunsflip=fliplr(chainruns); // We want to multiply on the left,
				// so reverse the order of the sequences
		mat Y=lvecprod(matrices,init,chainrunsflip,d);
		return Y;
}

// [[Rcpp::export]]
umat SimulateYsumstest(vec matrices, vec trans,uword K,uword d,uword N, uword couple,uvec times){
		rowvec q(K);
		q[couple]=1;
		uword T=trans.n_elem/(K*K);
		vec rnd=randu<vec>(N*(T+1));
		umat chainruns=inhomseq1(N,trans,q,rnd,times);
		vec init=ones<vec>(d);
		umat chainrunsflip=fliplr(chainruns); // We want to multiply on the left,
				// so reverse the order of the sequences
		mat Y=lvecprod(matrices,init,fliplr(chainruns),d);
		return chainruns;
}


/*  

	Compute transition matrices for backward coupled process
	Input 
		MCtrans= transition matrix for the Markov chain
		T= max coupling time
		statesin1,statesin2 = starting states (numbered 0..K-1)
	Variables
		Q = K x T matrix of probabilities for coupling
	Output
		vector of length K*K*(2T-2)+K*Tt=
			First two pieces are (T-1) transition matrices
			Last piece is Q.
		
*/


// [[Rcpp::export]]
vec coupletransitions(mat Pin,uword Tin,uvec statesin){
	mat P(Pin);
	uword T(Tin);
	uvec states(statesin);
	uword K=P.n_rows;
	cube transitions(K,K,2*T-2);
// time 0 is the time when the paths first split, corresponding to state1 and state2
// We have alpha for time 0 through time T, allowing us to compute Q for
// time 1 through T. This gives us just T-1 transition matrices, corresponding 
// time  T -> T-1, ... , 2 -> 1
//  The first T will correspond to alphaplus; the next T to alphaminus
	mat tP=P.t();
	mat alpha(K,T+1),Q(K,T);  // Note: Q(i,t) corresponds to time t+1
	alpha.zeros();
	Q.zeros();
	alpha(states[0],0)=1;
	alpha(states[1],0)=-1;
	for (uword t=0; t<T;t++){
		alpha.col(t+1)=tP*alpha.col(t);}
	mat alphaplus=matmax(alpha,0);//
	mat alphaminus=matmax(-alpha,0);
	for (uword t=0; t<T;t++){
		Q.col(t)=matmax(tP*alphaplus.col(t)-alphaplus.col(t+1),0);
		}
	// Now compute transition matrices for backward coupled chain
	for(uword t=0;t<T-1;t++){
		for(uword j=0;j<K;j++){
			for(uword i=0;i<K;i++){
				transitions(i,j,t)=P(j,i)*alphaplus(j,T-t-1);
				transitions(i,j,t+T-1)=P(j,i)*alphaminus(j,T-t-1);
				}
			}
			transitions.slice(t)=rowsum1(transitions.slice(t));
			transitions.slice(t+T-1)=rowsum1(transitions.slice(t+T-1));
		}
	vec transvec(transitions.begin(),K*K*(2*T-2));
	vec Qvec(Q.begin(),K*T);
	vec outputvector=join_cols(transvec,Qvec);
	return outputvector;
}



		
//    maxV
//
//			auxiliary function to compute the maximum value of ||XV||
//
//		Input matrices=vector of matrix values
//				d=matrix dimensions
//				V=n x d matrix of V samples

// [[Rcpp::export]]
double maxV(vec matrices, mat V){
	uword d=V.n_cols;
	uword n=V.n_rows;
	uword K=matrices.n_elem/(d*d);
	vec maxes(2*K);
	cube X(matrices.begin(),d,d,K,false);
	vec VX(n);
	for (uword i=0;i<K;i++){
		VX=sum(V*X.slice(i),1);
		maxes[i]=log(max(VX));
		maxes[K+i]=-log(min(VX));
	}
	return max(maxes);
}


// [[Rcpp::export]]
double maxU(vec matrices, mat U){
	uword d=U.n_rows;
	uword n=U.n_cols;
	uword K=matrices.n_elem/(d*d);
	vec maxes(2*K);
	cube X(matrices.begin(),d,d,K,false);
	rowvec XU(n);
	for (uword i=0;i<K;i++){
		XU=sum(X.slice(i)*U,0);
		maxes[i]=log(max(XU));
		maxes[K+i]=-log(min(XU));
	}
	return max(maxes);
}


//				LinReg			//
//
/*		routine to compute trends

			Input vector y
			
			output alpha, beta best fit to y[i]=alpha=beta*i
*/

// [[Rcpp::export]]
vec LinReg(vec y){
	vec beta(2);
	uword k=y.n_elem;
	double ybar=mean(y);
	for (uword i=0; i<k; i++){
		y[i]=y[i]*(i+1);}
	double yibar=mean(y);
	beta[0]=((2*k+1)*ybar-3*yibar)*2/(k-1);
	beta[1]=(-3*ybar+6*yibar/(k+1))*2/(k-1);
	return beta;
}


//  maxdiff= auxiliary function for QTbound that takes a vector as input, outputs
// 			maximum difference between coefficients

// [[Rcpp::export]]
double maxdiff(cx_vec vin){
	cx_vec v(vin);
	uword K=v.n_elem;
	double absdiff;
	double vijmax; // We're going to run through components to find max (v_i-v_j)
	for (uword i=0; i<K-1; i++){
		for (uword j=i+1; j<K; j++){
			absdiff=abs(v[i]-v[j]);
			if (absdiff>vijmax){vijmax=absdiff;}
		}
	}
	return vijmax;
}

/*******************************************************************************
 *  QTbound: Given a transition matrix P and depth T,
 *          compute a bound on the error in estimating zeta that can arise from truncating
 *          the sum at T.
 * 
 *  output = 2-place vector. First place is bound on sum of q_i,t * C_i,t for |C|<1
 *                          Second place is bound for |C_i,t| < T-t.
 ******************************************************************************/

// [[Rcpp::export]]
vec QTbound(mat P,uword T){
    vec boundoutput(2);  // First element will be the bound on sum(q(t);t>T),
    //  Second element is bound on sum(q(t)(T-t))
    cx_mat Pc=conv_to<cx_mat>::from(P);
    uword K=P.n_rows;
    cx_vec eigval;
    vec abseigval;
    cx_mat eigvecR;
    eig_gen(eigval,eigvecR,Pc); // right eigenvectors  
    cx_mat eigvecL=eigvecR.i();  // left eigenvectors
    abseigval=abs(eigval);
    vec vijmax(K),wsum(K),powers(K),recip(K);
    for (uword i=0; i<K; i++){
        vijmax[i]=maxdiff(eigvecR.col(i)); // Largest possible entry of (q -pi)vi
        wsum[i]=sum(abs(eigvecL.col(i)));
    }
    uword maxind=abseigval.index_max();  // maxind is now the index of the
    // top eigenvalue
    abseigval[maxind]=0;
    powers=pow(abseigval,T+1);
    for (uword i=0; i<K;i++){
        recip[i]=1/(1-abseigval[i]);}
    boundoutput[0]=sum(vijmax%wsum%powers); // first missing term
    boundoutput[1]=sum(vijmax%wsum%powers%recip);  // sum of all missing terms
    return boundoutput;
}
//******************************************************************
// 			QTchooseT
//
//		Input: transition matrix P
//				desired maximum error Qbound
//
//		Output: T=depth to which we need to go in computing zeta
//******************************************************************

// [[Rcpp::export]]
uword QTchooseT(mat PP,double Qboundin){
    mat P(PP);
    double Qbound(Qboundin);
    vec boundoutput(2);  // First element will be the bound on sum(q(t);t>T),
    //  Second element is bound on sum(q(t)(T-t))
    cx_mat Pc=conv_to<cx_mat>::from(P);
    uword K=P.n_rows;
    cx_vec eigval;
    vec abseigval;
    cx_mat eigvecR;
    eig_gen(eigval,eigvecR,Pc);
    cx_mat eigvecL=eigvecR.i();
    abseigval=abs(eigval);
    vec vijmax(K),wsum(K),powers(K),recip(K),error_bounds(K);
    for (uword i=0; i<K; i++){
        vijmax[i]=maxdiff(eigvecR.col(i));
        wsum[i]=sum(abs(eigvecL.col(i)));
    }
    uword maxind=abseigval.index_max(); 
    abseigval[maxind]=0;  // top eigenvalue/vector doesn't get used in the bound
    uword T=1;
    for (uword i=0; i<K;i++){
        recip[i]=1/(1-abseigval[i]);}
    error_bounds=vijmax % wsum % abseigval % recip;
    while (sum(error_bounds)>Qbound){
        T++;
        error_bounds=error_bounds % abseigval;}
    // Loop by multiplying by eigenvalues until error is small enough
    return T;
}

/*************************************************************************
 * 
 * arraycovariance: Input a matrix of values W_ij=F(V_i,U_j)
 * 
 * Output 4-place vector with estimates of variance(W_ij),
 *          cov(W_ij,W_i'j)     (same U)
 *          cov(W_ij,W_ij')     (same V)
 *          variance(mean(W_ij))
 * 
 *************************************************************************/

// [[Rcpp::export]]
vec arraycovariance(mat WW){
	mat W(WW);
	uword m=W.n_rows; // Number of V's
	uword n=W.n_cols; // Number of U's
	double meanW=mean(mean(W));
	double C1(0),C2(0),C3(0);
	vec hatsigma=zeros<vec>(4);
	rowvec meanU=mean(W,0);
	vec meanV=mean(W,1);
	C1=sum(sum(square(W-meanW)));
for (uword i=0; i<m ;i++){
		for (uword j=0; j< n ;j++){
			C3+= pow(W(i,j)-meanU[j],2);
			C2+= pow(W(i,j)-meanV[i],2);
		}
	}
	hatsigma[0]=(C1-C2/m-C3/n)/(m-1)/(n-1);
	hatsigma[1]=(C1-C2/m-C3)/(m-1)/(n-1);  //hatgammaU=covariance between F's with same U
	hatsigma[2]=(C1-C2-C3/n)/(m-1)/(n-1) ; //hatgammaV=covariance between F's with same V
	hatsigma[3]=hatsigma[0]/m/n+hatsigma[1]*(m-1)/m/n+hatsigma[2]*(n-1)/n/m;
	return hatsigma;}

/******************************************************************************************
 zetaintern

		Inputs the transition matrix, the projection matrices, and the two starting states,
          and outputs zeta for starting from those states for all possible coupling states
          
          transin (subtrans1)= vector of all transitions for the first reverse MC = K^2*(TY-2) vector
          transinprime (subtrans2)= vector of all transitions for the second reverse MC
          PP (P)=input transition matrix
          dd (d)= dimensions of projection matrices
          Nin (N)= Kx(TY) with number of samples (sequences of length 1 automatically zero)
          ToV (TV)= length of sequence in each try of sampling V
          ostarts (starts)= pair of starting states
          maxostep (maxstep)= estimated maximum log change due to one matrix product
      Output is a list of KxTY matrices: mean estimate of zeta, variance of the estimate, and bound on the systematic error.
            place (i,t) is the appropriate quantity for coupling in state i at time t

    The function zeta will compute the appropriate TY (the length of the longest coupling sequence)
           and a vector of N's of length TY-1
    Note: Every coupling sequence starts with the same matrices. So the sequence of length 1 has no random component,
        and the longest sequence has random component of length TY-1, which is why there are TY-2 transition matrices

******************************************************************************************/

// [[Rcpp::export]]
List zetaintern(vec mats,vec transin,vec transinprime, mat PP,uword dd,umat Nin,uword ToV,uvec ostarts,double maxoerror,double maxobias,double maxostep){
    mat P(PP);
    uword K=P.n_rows;
    uword TY=Nin.n_cols;     // We simulate up to coupling at time TY after the split.
						// Note that coupling at time 1 contributes 0 to zeta.
    uword d(dd),  TV(ToV);
	vec matrices(mats);
    cube matrixcube(matrices.begin(),d,d,K,false);  // Need access to matrices as matrices
    mat currentmatrix(d,d);
  umat N(Nin);
	uvec starts(ostarts);
// We need initial vectors = X_e * 1 to feed to the Y simulator
    vec init0=sum(matrixcube.slice(starts[0]),1);
    vec init1=sum(matrixcube.slice(starts[1]),1);
	double maxerror(maxoerror), maxstep(maxostep);
  vec QTb=QTbound(P,TY);
  double zetabias2=QTb[1]*maxstep; // This is one piece of the bias; this piece
                    // arises from truncating at T.
					// There's still the bias component due to differences in zeta, and
          // a matrix of biases corresponding to the errors in sampling V.
	mat zetamean=zeros<mat>(K,TY); // Store estimated mean values of components of zeta
	mat zetaVar=zeros<mat>(K,TY); // Store estimated variance of components of zeta
	mat zetabias=zeros<mat>(K,TY);
    // Note: We may be able to speed this up by defining the size of these matrices
    uword maxN=max(max(N));
	mat Yvecs(d,maxN);
	mat Yprimevecs(d,maxN);
	mat Vsamples(d,maxN);
    vec Vsums(maxN);
	vec subtrans;
	mat F;
	uvec times(2);
	times[1]=TY-1;  // We're always running up to the last transition.
	vec subtrans1(transin),subtrans2(transinprime);
  vec cutoff(K),cutoffprime(K),meanlogvy;
  double meanvy;
	vec acov(4); //Storing return of variance/covariance estimates for F
	for (uword j=0; j<K; j++){ // Simulate Y corresponding to coupling at state j
        currentmatrix=matrixcube.slice(j);
	    // N(j,t) is the number of simulations to make for coupling
	    //      at state j at time t.
        if (N(j,0)>1){
            Vsamples=SimulateV(matrices,P,d,N(j,0),TV,j,maxerror); //Take #Vsamples to be N(j,0), though that's actually for 2-step coupling
            zetabias(j,0)=mean(Vsamples.col(d));
            Vsamples.shed_col(d);
            Vsums=log(Vsamples*currentmatrix*init0)-log(Vsamples*currentmatrix*init1);
            // Note: no need to simulate in this case, because we've gone directly from the split to the coupling state,
            //   where we're already conditioning on it.
            zetamean(j,0)=sum(Vsums)/N(j,0);
            zetaVar(j,0)=var(Vsums)/N(j,0);}
    for (uword t=1; t<TY; t++){ // Simulate Y's t+1 steps long, corresponding to coupling after t+1 steps. (Final step takes us 
                                // back to the splitting state.)
                                // Note that the last transition matrix only takes us back to 1 step after the split
		  times[0]=TY-t-1; 
			if(N(j,t)>0){
            Yvecs=SimulateYsums(matrices,subtrans1,K,d,N(j,t),j,times,init0);
			  Yprimevecs=SimulateYsums(matrices,subtrans2,K,d,N(j,t),j,times,init1);
			  Vsamples=SimulateV(matrices,P,d,N(j,t),TV,j,maxerror);
			  zetabias(j,t)=mean(Vsamples.col(d));  //Last column gives the diameters
			    // These are bounds on the finite-sample error.
			  Vsamples.shed_col(d);
			  F=log(Vsamples*Yvecs)-log(Vsamples*Yprimevecs);
			  zetamean(j,t)=sum(sum(F))/N(j,t)/N(j,t);
			  acov=arraycovariance(F);
			  zetaVar(j,t)=acov[3];
			}
		}
    meanlogvy=mean(log(Vsamples*Yvecs),1);
    meanvy=mean(meanlogvy);
    cutoff[j]=max(abs(meanlogvy-meanvy));
    meanlogvy=mean(log(Vsamples*Yprimevecs),1);
    meanvy=mean(meanlogvy);
    cutoffprime[j]=max(abs(meanlogvy-meanvy));
	}
  cx_vec lastcol=conv_to<cx_vec>::from(zetamean.col(TY-2));
  zetabias2 +=QTb[0]*maxdiff(lastcol);
  vec zetabias3(2);
  zetabias3[0]=max(cutoff);
  zetabias3[1]=max(cutoffprime);
  List ret; 
  ret["Mean"] = zetamean; 
  ret["Variance"] = zetaVar;
  ret["Bias"] = zetabias;  // errors in stationary samples
  ret["Cutoff"] = zetabias2; // errors from stopping at T
  ret["Cutoff2"] = zetabias3; // errors from differences in zeta
  return(ret);
}

// [[Rcpp::export]]
mat zetatest(vec mats, mat PP,uword dd,uword Nin,uword ToV,uvec ostarts){
  mat P(PP);
    vec matrices(mats);
    uword K=P.n_rows;
  uword N(Nin);
  mat outsizes(2,N); // This is going to store the output norms
	uword d(dd),  TV(ToV);
    vec rnd=randu<vec>(N*(TV+1));
    rowvec q=zeros<rowvec>(K);
    q[ostarts[0]]=1;
    umat sequences=MC(N,TV,P,q,rnd);
    mat init=ones<mat>(d,1);
	cube testproduct=lmatprod(matrices,init,sequences,d);
	 for (uword whichain=0;whichain<N; whichain++){
	    outsizes(0,whichain)=log(sum(sum(testproduct.slice(whichain))));
	 }
     q[ostarts[0]]=0;
     q[ostarts[1]]=1;
     sequences=MC(N,TV,P,q,rnd);
     testproduct=lmatprod(matrices,init,sequences,d);
     for (uword whichain=0;whichain<N; whichain++){
        outsizes(1,whichain)=log(sum(sum(testproduct.slice(whichain))));
	 }
     mat outp(3,2);
     outp(0,0)=mean(outsizes.row(0));
     outp(1,0)=mean(outsizes.row(1));
     outp(0,1)=stddev(outsizes.row(0));
     outp(1,1)=stddev(outsizes.row(1));
     outp(2,0)=outp(0,0)-outp(1,0);
     outp(2,1)=stddev(outsizes.row(0)-outsizes.row(1))/sqrt(N);
     return outp;
}
    

/****************************************************************
						xi
	Compute the factors xi(e,e'), the second component in sensitivity
		to change in driving process
	
	xi(e',e)= E[log VX_e X_e' U/ ||VX_e'||]
	
	Input variables:
		matsin (mats) = vector of all X matrix entries
		PP (P) = transition matrix
		din (d) = dimension of X
		NUin(NU) and NVin (NV) = number of U and V samples
		startsin (starts) = (e',e)
		maxerrorin (maxerror) = diameter of matrices producing U and V *3
				(because the actual error is 3 times as big)

allxi takes all inputs except startsin

    outputs K x 3K matrix. 3 blocks of KxK giving mean, variance, bias for xi(i,j)
****************************************************************/

// [[Rcpp::export]]
vec xi(vec matsin, mat PP,uword din,uword NUin,uword NVin,uword TVin,uvec startsin,double maxerrorin){
	vec mats(matsin);
    vec xioutput=zeros<vec>(3); // (mean estimates; Vars ; bias)
	mat P(PP);
	uword d(din),NU(NUin),NV(NVin),TV(TVin);
	rowvec vrow(d);
	uvec starts(startsin);
	double maxerror(maxerrorin);
	uword K=P.n_rows;
	cube X(mats.begin(),d,d,K,false);
	mat oneV(NV,d+1),oneU(d+1,NU);
	mat uvmat(NV,NU);
	vec hsigma(4);
			oneV=SimulateV(mats,P,d,NV,TV,starts[1],maxerror);//Simulate vectors from stationary distribution of left vectors from second starting pt
			oneU=SimulateU(mats,P,d,NU,TV,starts[0],maxerror);// Simulate vectors from stationary distribution of right vectors from first starting pt
			xioutput[2]=mean(oneV.col(d))+mean(oneU.row(d)); // These are the diameters
			oneV.shed_col(d);
            oneU.shed_row(d); // don't need the diameters anymore; cut them out; we're left with just the simulated vectors
			oneV=oneV*(X.slice(starts[1]));  // V's get multiplied by the second X
            oneV=rowsum1(oneV); //normalise all the rows; we now have V X/||VX||
			uvmat=log(oneV*X.slice(starts[0])*oneU); 
			xioutput[0]=sum(sum(uvmat))/NU/NV;
			hsigma=arraycovariance(uvmat); // Variances and covariances
			xioutput[1]=hsigma[3]; 
	return xioutput;
}

// [[Rcpp::export]]
mat allxi(vec matsin, mat PP,uword din,uword NUin,uword NVin,uword TVin,double maxerrorin){
	vec mats(matsin);
	mat P(PP);
	uword d(din),NU(NUin),NV(NVin),TV(TVin);
	rowvec vrow(d);
	uvec starts(2);
	double maxerror(maxerrorin);
	uword K=P.n_rows;
	cube X(mats.begin(),d,d,K,false);
	mat oneV(NV,d+1),oneU(d+1,NU);
	vec diamV(NV),diamU(NU);
	mat xiout=zeros<mat>(K,3*K); // (mean estimates; SEs ; bias)
	mat uvmat(NV,NU);
	vec retxi(3);
	for (uword eprime=0;eprime<K; eprime++){
		for (uword e=0;e<K; e++){
			starts[0]=e;
			starts[1]=eprime;
			retxi=xi(mats,P,d,NU,NV,TV,starts,maxerror);
			xiout(e,eprime)=retxi[0];
			xiout(e,eprime+K)=retxi[1];
			xiout(e,eprime+K*2)=retxi[2];
		}
	}
	return xiout;
}

/****************************************************************
 allxi_iid
Compute the factors xi(e), the sensitivity
to change in driving process in the iid case

xi(e)= E[log || VX_e ||]

Input variables:
matsin (mats) = vector of all X matrix entries
nu = probabilities
din (d) = dimension of X
NUin(NU) and NVin (NV) = number of U and V samples
startsin (starts) = (e',e)
maxerrorin (maxerror) = diameter of matrices producing U and V *3
(because the actual error is 3 times as big)

allxi_iid takes all inputs except startsin

outputs K x 3K matrix. 3 blocks of KxK giving mean, variance, bias for xi(i,j)
****************************************************************/

// [[Rcpp::export]]
mat allxi_iid(vec matsin, vec nuin,uword din,uword NUin,uword NVin,uword TVin,double maxerrorin){
    vec mats(matsin);
    vec nu(nuin);
    uword K=nu.n_elem;
    vec xiout=zeros<vec>(2*K+1); // (mean estimates; Vars ; bias)
    uword d(din),NU(NUin),NV(NVin),TV(TVin);
    double maxerror(maxerrorin);
    cube X(mats.begin(),d,d,K,false);
    mat oneV(NV,d+1),oneU(d+1,NU);
    mat uvmat(NV,NU);
    vec hsigma(4);
    oneV=SimulateV_iid(mats,nu,d,NV,TV,maxerror);//Simulate vectors from stationary distribution of left vectors from second starting pt
    oneU=SimulateU_iid(mats,nu,d,NU,TV,maxerror);// Simulate vectors from stationary distribution of right vectors from first starting pt
    xiout(2*K) = mean(oneV.col(d)) +mean(oneU.row(d)); // These are the diameters
    oneV.shed_col(d);
    oneU.shed_row(d); // don't need the diameters anymore; cut them out; we're left with just the simulated vectors
    for (uword e=0;e<K; e++){
//        oneV=oneV*(X.slice(0));  // V's get multiplied by the second X
        //oneV=rowsum1(oneV); //normalise all the rows; we now have V X/||VX||
        uvmat=log(oneV*X.slice(e)*oneU);
        xiout(e)=sum(sum(uvmat))/NU/NV;
        hsigma=arraycovariance(uvmat); // Variances and covariances
        xiout(e+K)=hsigma[3];
    }
    return xiout;
}


/***************************************************************
            sensitivityij: Compute sensitivity wrt changing the projection matrices
        mats(matsin)= Projection matrices as vector
        P(PP) = transition matrix
        d(din)= dimensions of projection matrices
        NU(NUin)= number of U samples to take
        NV(NVin)= number of V samples to take
        TV(TVin)= depth for initial try for U and V samples
		ij(coords)= (i,j,e) where we are computing the derivative wrt
				changing the (i,j) component of X_e
        maxerror(maxerrorin)= Maximum allowed systematic error in U and V samples
    
    Output: 3-place vector (mean estimate, variance, bias)

***************************************************************/
			
// [[Rcpp::export]]
vec sensitivityij(vec matsin, mat PP,uword din,uword NUin,uword NVin,uword TVin,uvec ij,double maxerrorin){
	vec outputs(3);
    vec mats(matsin);
	uvec coords(ij);
	mat P(PP);
	uword d(din),NU(NUin),NV(NVin),TV(TVin);
	mat vXu(NV,NU);
	rowvec vrow(d);
	double maxerror(maxerrorin);
	uword K=P.n_rows;
	vec q0=zeros<vec>(K);
	vec q1=zeros<vec>(K);
	q0[coords[0]]=1;
	q1[coords[1]]=1;
	mat X(mats.begin()+d*d*coords[2u],d,d,true);
	mat oneV(NV,d+1),oneU(d+1,NU);
	vec diamV(NV),diamU(NU);
	oneV=SimulateV(mats,P,d,NV,TV,coords[2],maxerror);
	oneU=SimulateU(mats,P,d,NU,TV,coords[2],maxerror);
	diamV=oneV.col(d);
	diamU=oneU.row(d).t();// turn into column vector
	oneV.shed_col(d);
	oneU.shed_row(d);
	vXu=(oneV.col(coords[0])*oneU.row(coords[1]))/(oneV*X*oneU);
	vec CUvec=sum(vXu,1);
	rowvec CVvec=sum(vXu,0);
	double CU=1+max(CUvec)/NU;  // upperbound for mean over U_i
	double CV=1+max(CVvec)/NV;  // upperbound for mean over V_i
	outputs[0]=sum(sum(vXu))/NV/NU;
	vec hsigma=arraycovariance(vXu);
	outputs[1]=hsigma[3];
	outputs[2]=CU*mean(diamU)+CV*mean(diamV);
	return outputs;
}


/*********************************************************************
 *
 *          allsens: Estimate all sensitivities wrt any change in any projection matrix
 *  
 *  Input variables as in sensitivities1
 * 
 * Output cube d x 3d x K
 *      Each slice (choice of last coordinate) defines a choice of matrix being changed
 *      Within a slice, position (i,j) is the estimate for sensitivity wrt change at (i,j)
 *                      position (i,j+d) is the variance for sensitivity wrt change at (i,j)
 *                      position (i,j+2d) is the bias for sensitivity wrt change at (i,j)
 * 
 **********************************************************************/

// [[Rcpp::export]]				
cube allsens1(vec matsin, mat PP,uword din,uword NUin,uword NVin,uword TVin,double maxerrorin){
	vec mats(matsin);
	mat P(PP);
	uword d(din),NU(NUin),NV(NVin),TV(TVin);
	mat vXu(NV,NU);
	rowvec vrow(d);
	double maxerror(maxerrorin);
	uword K=P.n_rows;
	uvec coords(3);
	vec sensij(3);
	cube outcube(d,3*d,K);
	for (uword e=0; e<K; e++){ 
	for (uword i=0; i<d; i++){
		for (uword j=0; j<d; j++){
			coords[0]=i;
			coords[1]=j;
			coords[2]=e;
			sensij=sensitivityij(mats,P,d,NU,NV,TV,coords,maxerror);
			outcube(i,j,e)=sensij[0];
			outcube(i,j+d,e)=sensij[1];
			outcube(i,j+2*d,e)=sensij[2];
		}
	}
	}
	return outcube;
}

/*****************************************************************************
 * 
 * stochgrowthU/stochgrowthV: Compute stochastic growth rates, based on 
 *                              left/right multiplication of column/row vectors
 * 
 *  Input:         
 *      mats(matsin)= Projection matrices as vector
 *      P(PP) = transition matrix
 *      d(din)= dimensions of projection matrices
 *      N(Nin)= number of U or V samples to take
 *      T(Tin)= depth for initial try for U and V samples
 *      start = starting state
 *      maxerror(maxerrorin)= Maximum allowed systematic error in U and V samples
 * 
 *  Output: 3-place vector giving mean estimate, variance, bias
 * 
 * Note: averaged wrt stationary distribution stochgrowthU and stochgrowthV must be the
 *          same, up to random error and bias
 *****************************************************************************/


// [[Rcpp::export]]				
vec stochgrowthU(vec matsin, mat PP,uword din,uword Nin,uword Tin,uword start,double maxerrorin){
	vec outs(3);
	vec mats(matsin);
	mat P(PP);
	uword d(din),N(Nin),T(Tin);
	mat X(mats.begin()+d*d*start,d,d,true);
	double maxerror(maxerrorin);
	mat U=SimulateU(mats,P,d,N,T,start,maxerror);
	rowvec Udiam=U.row(d);
	U.shed_row(d);
	mat XU=X*U;
	rowvec XUlog=log(sum(XU,0));
	outs[0]=sum(XUlog)/N;
	outs[1]=var(XUlog)/N;
	outs[2]=mean(Udiam);
	return outs;
}

// [[Rcpp::export]]				
vec stochgrowthV(vec matsin, mat PP,uword din,uword Nin,uword Tin,uword start,double maxerrorin){
	vec outs(3);
	vec mats(matsin);
	mat P(PP);
	uword d(din),N(Nin),T(Tin);
	mat X(mats.begin()+d*d*start,d,d,true);
	double maxerror(maxerrorin);
	mat V=SimulateV(mats,P,d,N,T,start,maxerror);
	vec Vdiam=V.col(d);
	V.shed_col(d);
	mat XV=V*X;
	vec XVlog=log(sum(XV,1));
	outs[0]=sum(XVlog)/N;
	outs[1]=var(XVlog)/N;
	outs[2]=mean(Vdiam);
	return outs;
}

/*****************************************************************************
 * 
 * stochgrowthallU/V: Compute stochastic growth rates for all starting states, based on 
 *                              left/right multiplication of column/row vectors
 * 
 *  Input:         
 *      mats(matsin)= Projection matrices as vector
 *      P(PP) = transition matrix
 *      d(din)= dimensions of projection matrices
 *      N(Nin)= number of U or V samples to take
 *      T(Tin)= depth for initial try for U and V samples
 *      maxerror(maxerrorin)= Maximum allowed systematic error in U and V samples
 *      whichside = 0 for call to stochgrowthV, otherwise stochgrowthU
 * 
 *  Output: 3 x K matrix giving mean estimate, variance, bias for all K possible initial states
 * 
 * Note: averaged wrt stationary distribution stochgrowthU and stochgrowthV must be the
 *          same, up to random error and bias
 *****************************************************************************/
// [[Rcpp::export]]				
mat stochgrowthall(vec matsin, mat PP,uword din,uvec Nin,uword Tin,double maxerrorin,uword whichside){
	vec mats(matsin);
	mat P(PP);
	double maxerror(maxerrorin);
	uword d(din),T(Tin);
  uvec N(Nin);
	uword K=P.n_rows;
	mat outputs(3,K);
    if (whichside==0){
	    for (uword i=0; i<K ; i++){
		    outputs.col(i)=stochgrowthV(mats,P,d,N[i],T,i,maxerror);
	    }
    }
    else{
        for (uword i=0; i<K ; i++){
		    outputs.col(i)=stochgrowthU(mats,P,d,N[i],T,i,maxerror);
	    }
    }
	return outputs;
}