#pragma once

#include "linpack/linpack_d.hpp"

#define FMT_HEADER_ONLY 
#include "../fmt/include/fmt/format.h"
#include "../fmt/include/fmt/ranges.h"
#include "../fmt/include/fmt/color.h"

#include <iostream>
#include <vector>
#include <cassert>
#include <iomanip>
#include <fstream>

#define assertm(exp, msg) assert(((void)msg, exp))


//**********************************************************************
//  this package solves boundary value problems for
//  ordinary differential equations with constraints,
//  as described below. for more see [1].
//
//  COLDAE is a modification of the package COLNEW by bader and
//  ascher [4], which in turn is a modification of COLSYS by ascher,
//  christiansen and russell [3]. it does what colnew does plus
//  optionally solving semi-explicit differential-algebraic equations
//  with index at most 2.
//**********************************************************************
//----------------------------------------------------------------------
//                            p a r t  1
//        main storage allocation and program control subroutines
//----------------------------------------------------------------------
//
//**********************************************************************
//
//     written by
//                  uri ascher and ray spiteri,
//                            department of computer science,
//                            university of british columbia,
//                            vancouver, b. c., canada   v6t 1z2
//
//**********************************************************************
//
//     purpose
//
//     this package solves a multi-point boundary value
//     problem for a mixed order system of ode-s with constraints
//     given by
//
//          (m(i))
//         u       =  f  ( x; z(u(x)), y(x) )   i = 1, ... ,ncomp
//          i          i
//
//             0   =  f  ( x; z(u(x)), y(x) )   i = ncomp+1,...,ncomp+ny
//                     i
//
//                                          aleft .lt. x .lt. aright,
//
//
//         g  ( zeta(j); z(u(zeta(j))) ) = 0   j = 1, ... ,mstar
//          j
//                                    mstar = m(1)+m(2)+...+m(ncomp),
//
//
//         where                          t                       t
//               u = (u , u , ... ,u     ) , y = (y ,y , ... ,y  )
//                     1   2        ncomp          1  2        ny
//
//         is the exact solution vector: u are the differential solution
//         components and y are the algebraic solution components.
//
//                (mi)
//               u     is the mi=m(i) th  derivative of u
//                i                                      i
//
//                                  (1)        (m1-1)       (mncomp-1)
//               z(u(x)) = ( u (x),u  (x),...,u    (x),...,u      (x) )
//                            1     1          1            ncomp
//
//                f (x,z(u),y)   is a (generally) nonlinear function of
//                 i
//                             z(u)=z(u(x)) and y=y(x).
//
//                g (zeta(j);z(u))  is a (generally) nonlinear function
//                 j
//                               used to represent a boundary condition.
//
//         the boundary points satisfy
//               aleft .le. zeta(1) .le. .. .le. zeta(mstar) .le. aright
//
//         the orders mi of the differential equations satisfy
//                            1 .le. m(i) .le. 4.
//
//         regarding the dae, note that:
//         i)  with ny=0, the code is essentially identical to the ode
//             code  colnew
//         ii) no explicit checking of the index of the problem
//             is provided. if the index is > 2 then the code will
//             not work well.
//         iii) the constraints are treated like de-s of order 0 and
//             correspondingly the approximation to y is sought in
//             a piecewise discontinuous polynomial space.
//         iv) the number of boundary conditions required is independent
//             of the index. it is the user's responsibility to ensure
//             that these conditions are consistent with the constraints.
//             the conditions at the left end point aleft  must include
//             a subset equivalent to specifying the index-2
//             constraints there.
//         v)  for an index-2 problem in hessenberg form, the projected
//             collocation method of ascher and petzold [2] is used.
//         vi) if the constraints are of a mixed type (and possibly
//             a mixed index) then coldae can transform and project
//             appropriately -- see description of ipar(12).
//
//**********************************************************************
//
//     method
//
//        the method used to approximate the solution u is
//     collocation at gaussian points, requiring m(i)-1 continuous
//     derivatives in the i-th component, i = 1, ..., ncomp.
//     here, k is the number of collocation points (stages) per
//     subinterval and is chosen such that k .ge. max m(i).
//     a runge-kutta-monomial solution representation is utilized.
//     for hessenberg index-2 daes, a projection on the constraint
//     manifold at each interval's end is used [2].
//
//     references
//
//     [1] u. ascher and r. spiteri,
//         collocation software for boundary value differential-
//         algebraic equations,
//         siam j. scient. stat. comput., to appear.
//
//     [2] u. ascher and l. petzold,
//         projected implicit runge-kutta methods for differential-
//         algebraic equations,
//         siam j. num. anal. 28 (1991), 1097-1120.
//
//     [3] u. ascher, j. christiansen and r.d. russell,
//         collocation software for boundary-value odes,
//         acm trans. math software 7 (1981), 209-222.
//
//     [4] g. bader and u. ascher,
//         a new basis implementation for a mixed order
//         boundary value ode solver,
//         siam j. scient. stat. comput. 8 (1987), 483-500.
//
//     [5] u. ascher, j. christiansen and r.d. russell,
//         a collocation solver for mixed order
//         systems of boundary value problems,
//         math. comp. 33 (1979), 659-679.
//
//     [6] u. ascher, j. christiansen and r.d. russell,
//         colsys - a collocation code for boundary
//         value problems,
//         lecture notes comp.sc. 76, springer verlag,
//         b. childs et. al. (eds.) (1979), 164-185.
//
//     [7] c. deboor and r. weiss,
//         solveblok: a package for solving almost block diagonal
//         linear systems,
//         acm trans. math. software 6 (1980), 80-87.
//


//**********************************************************************
//
//     *************    user supplied subroutines   *************
//
//
//     the following subroutines must be declared external in the
//     main program which calls coldae.
//
//
//     fsub  - name of subroutine for evaluating f(x,z(u(x)),y(x)) =
//                               t
//             (f ,...,f        )  at a point x in (aleft,aright).  it
//               1      ncomp+ny
//             should have the heading
//
//                       subroutine fsub (x , z , y, f)
//
//             where f is the vector containing the value of fi(x,z(u),y)
//             in the i-th component and y(x) = (y(1),...y(ny)) and
//                                       z(u(x))=(z(1),...,z(mstar))
//             are defined as above under  purpose .
//
//
//     dfsub - name of subroutine for evaluating the jacobian of
//             f(x,z(u),y) at a point x.  it should have the heading
//
//                       subroutine dfsub (x , z , y, df)
//
//             where z(u(x)) and y(x) are defined as for fsub and the
//             (ncomp+ny) by (mstar+ny)
//             array  df  should be filled by the partial deriv-
//             atives of f, viz, for a particular call one calculates
//                          df(i,j) = dfi / dzj, i=1,...,ncomp+ny
//                                               j=1,...,mstar,
//                          df(i,mstar+j) = dfi / dyj, i=1,...,ncomp+ny
//                                               j=1,...,ny.
//
//
//     gsub  - name of subroutine for evaluating the i-th component of
//             g(x,z(u(x))) = g (zeta(i),z(u(zeta(i)))) at a point x =
//                             i
//             zeta(i) where 1.le.i.le.mstar. it should have the heading
//
//                       subroutine gsub (i , z , g)
//
//             where z(u) is as for fsub, and i and g=g  are as above.
//                                                     i
//             note that in contrast to f in  fsub , here
//             only one value per call is returned in g.
//             also, g is independent of y and, in case of a higher
//             index dae (index = 2), should include the constraints
//             sampled at aleft  plus independent additional constraints,
//             or an equivalent set.
//
//
//     dgsub - name of subroutine for evaluating the i-th row of
//             the jacobian of g(x,z(u(x))).  it should have the heading
//
//                       subroutine dgsub (i , z , dg)
//
//             where z(u) is as for fsub, i as for gsub and the mstar-
//             vector dg should be filled with the partial derivatives
//             of g, viz, for a particular call one calculates
//                   dg(i,j) = dgi / dzj      j=1,...,mstar.
//
//
//     guess - name of subroutine to evaluate the initial
//             approximation for  z(u(x)), y(x) and for dmval(u(x))= vector
//             of the mj-th derivatives of u(x). it should have the
//             heading
//
//                       subroutine guess (x , z , y, dmval)
//
//             note that this subroutine is needed only if using
//             ipar(9) = 1, and then all  mstar  components of z,
//             ny components of  y
//             and  ncomp  components of  dmval  should be specified
//             for any x,  aleft .le. x .le. aright .
//


//**********************************************************************
//
//     ************   use of output from coldae   ************
//
//                 ***   solution evaluation   ***
//
//     on return from coldae, the arrays fspace and ispace
//     contain information specifying the approximate solution.
//     the user can produce the solution vector  z( u(x) ), y(x)  at
//     any point x, aleft .le. x .le. aright, by the statement,
//
//           call appsln (x, z, y, fspace, ispace)
//
//     when saving the coefficients for later reference, only
//     ispace(1),...,ispace(8+ncomp)    and
//     fspace(1),...,fspace(ispace(8))    need to be saved as
//     these are the quantities used by appsln.
//
//
//                 ***   simple continuation   ***
//
//
//     a formerly obtained solution can easily be used as the
//     first approximation for the nonlinear iteration for a
//     new problem by setting   (iguess =) ipar(9) = 2, 3 or 4.
//
//     if the former solution has just been obtained then the
//     values needed to define the first approximation are
//     already in ispace and fspace.
//     alternatively, if the former solution was obtained in a
//     previous run and its coefficients were saved then those
//     coefficients must be put back into
//     ispace(1),..., ispace(8+ncomp)    and
//     fspace(1),..., fspace(ispace(8)).
//
//     for ipar(9) = 2 or 3 set ipar(3) = ispace(1) ( = the
//     size of the previous mesh ).
//
//     for ipar(9) = 4 the user specifies a new mesh of n subintervals
//     as follows.
//     the values in  fspace(1),...,fspace(ispace(8))  have to be
//     shifted by n+1 locations to  fspace(n+2),..,fspace(ispace(8)+n+1)
//     and the new mesh is then specified in fspace(1),..., fspace(n+1).
//     also set ipar(3) = n.
//


//**********************************************************************
//
//     ***************      package subroutines      ***************
//
//     the following description gives a brief overview of how the
//     procedure is broken down into the subroutines which make up
//     the package called  coldae . for further details the
//     user should refer to documentation in the various subroutines
//     and to the references cited above.
//
//     the subroutines fall into four groups:
//
// part 1 - the main storage allocation and program control subr
//
//     coldae - tests input values, does initialization and breaks up
//              the work areas, fspace and ispace, into the arrays
//              used by the program.
//
//     contrl - is the actual driver of the package. this routine
//              contains the strategy for nonlinear equation solving.
//
//     skale  - provides scaling for the control
//              of convergence in the nonlinear iteration.
//
//
// part 2 - mesh selection and error estimation subroutines
//
//     consts - is called once by  coldae  to initialize constants
//              which are used for error estimation and mesh selection.
//
//     newmsh - generates meshes. it contains the test to decide
//              whether or not to redistribute a mesh.
//
//     errchk - produces error estimates and checks against the
//              tolerances at each subinterval
//
//
// part 3 - collocation system set-up subroutines
//
//     lsyslv - controls the set-up and solution of the linear
//              algebraic systems of collocation equations which
//              arise at each newton iteration.
//
//     gderiv - is used by lsyslv to set up the equation associated
//              with a side condition point.
//
//     vwblok - is used by lsyslv to set up the equation(s) associated
//              with a collocation point.
//
//     gblock - is used by lsyslv to construct a block of the global
//              collocation matrix or the corresponding right hand
//              side.
//
//
// part 4 - service subroutines
//
//     appsln - sets up a standard call to  approx .
//
//     approx - evaluates a piecewise polynomial solution.
//
//     rkbas  - evaluates the mesh independent runge-kutta basis
//
//     vmonde - solves a vandermonde system for given right hand
//              side
//
//     horder - evaluates the highest order derivatives of the
//              current collocation solution used for mesh refinement.
//
//
// part 5 - linear algebra  subroutines
//
//     to solve the global linear systems of collocation equations
//     constructed in part 3,  coldae  uses a column oriented version
//     of the package  solveblok originally due to de boor and weiss.
//
//     to solve the linear systems for static parameter condensation
//     in each block of the collocation equations, the linpack
//     routines  dgefa and  dgesl  are included. but these
//     may be replaced when solving problems on vector processors
//     or when solving large scale sparse jacobian problems.
//----------------------------------------------------------------------



//------------------------------------------------------------------------------------------------------


template<typename T>
class arrRef1 {
	template<typename T> friend class arrRef2;
	template<typename T> friend class arrData1;
protected:
	T* data;
	int size;
	arrRef1(){}
public:
	void assertDim(int dim0) {
		if (dim0 == 1) return; // we assume that DIMENSION(1) does not give the actual dimension and only makes sure it is an array

		if (size > dim0) return; // allow downsizing

		assertm(size == dim0, "Dimension 1 does not match!");
	}

	inline const T& operator()(int idx) const {
		assertm(idx >= 1 && idx <= size, "Dimension 1 index out of range!");
		return data[idx - 1];
	}
	inline T& operator()(int idx) {
		assertm(idx>=1 && idx <= size,"Dimension 1 index out of range!");
		return data[idx - 1];
	}
	arrRef1 sub(int idx) {
		arrRef1 s;
		s.data = data + (idx - 1);
		s.size = size - (idx - 1);
		return s;
	}

	inline arrRef1(arrRef2<T> const& ar2) {
		data = ar2.data;
		size = ar2.size1 * ar2.size2;
	}


	T* begin() const { return data; }
	T* end() const { return data + size; }

	T* contiguous() {
		return data;
	}
	int getSize() { return size; }
};

template<typename T>
class arrData1 : public arrRef1<T> {
	std::vector<T> mydata;
public:
	arrData1(int s = 0) {
		mydata.resize(s);
		this->data = mydata.data();
		this->size = mydata.size();
	}
	arrData1(std::initializer_list<T> l) {
		mydata.clear();
		mydata.insert(mydata.end(), l.begin(), l.end());
		this->data = mydata.data();
		this->size = mydata.size();
	}
	arrData1(arrRef1<T> const& ar1) {
		// copy from another existing data
		mydata.resize(ar1.size);
		std::copy(ar1.begin(), ar1.end(), mydata.begin());
		this->data = mydata.data();
		this->size = mydata.size();
	}
	void copyFrom(arrRef1<T> const& ar1) {
		// copy from another existing data
		assertm(ar1.size <= mydata.size(), "my size is too low");
		std::copy(ar1.begin(), ar1.end(), mydata.begin());
		this->data = mydata.data();
		/// size remains
	}
};



template<typename T>
class arrRef2 {
	template<typename T> friend class arrRef1;
protected:
	T* data=nullptr;
	int size1=0, size2=0;
	int cap=0;
	arrRef2(){}
public:
	void reshape(int s1, int s2) {
		assertm(s1*s2 <= cap, "Capacity does not allow this reshape!");
		size1 = s1;
		size2 = s2; // might be 1 for unspecified
	}
	void assertDim(int s1, int s2) {
		assertm(s1 * s2 <= cap, "Capacity too low!");

		
		if (size1 == 0) {
			// The reference comes from a one-dimensional and the dimensions are not yet specified.
			// In this case we set the dimensions now.
			size1 = s1;
			size2 = s2;
		}
		//else if (s2 == 1) {
		//	// The last dimension is left unspecified
		//	assertm(size1 == s1, "Dimension 1 does not match!");
		//}
		//else if (size2 == 1) {
		//	// the last dimension was unspecified. We leave it along
		//	assertm(size1 == s1, "Dimension 1 does not match!");
		//}
		//else {
		//	assertm(size1 == s1, "Dimension 1 does not match!");
		//	assertm(size2 == s2, "Dimension 2 does not match!");
		//}
	}
	inline arrRef2(arrRef1<T> const& ar1) {
		data = ar1.data;
		size1 = 0; size2 = 1; // set the matrix size to 0,0. It will be specified later (hopefully)
		cap = ar1.size;
	}
	inline T& operator()(int idx1, int idx2) {
		assertm(idx1 >= 1 && idx1 <= size1, "Dimension 1 index out of range!");
		if(size2 != 1)
			assertm(idx2 >= 1 && idx2 <= size2, "Dimension 2 index out of range!");
		return data[(idx2 - 1)*size1 + (idx1-1)];
	}
	arrRef2 sub(int idx1, int idx2) {
		arrRef2 s;
		s.data = data + (idx2 - 1) * size1 + (idx1 - 1);
		s.size1 = size1 - (idx1 - 1);
		s.size2 = size2 - (idx2 - 1);
		s.cap = s.size1 * s.size2;
		return s;
	}
	T* contiguous() {
		return data;
	}

	void print() {
		for (int i = 1; i <= size1; ++i) {
			for (int j = 1; j <= size2; ++j) {
				std::cout << std::setw(10) << (*this)(i,j);
			}
			std::cout << std::endl;
		}
	}
};

template<typename T>
class arrData2 : public arrRef2<T> {
	std::vector<T> mydata;
public:
	arrData2(){}
	arrData2(int s1, int s2) {
		mydata.resize(s1*s2);
		this->data = mydata.data();
		this->size1 = s1;
		this->size2 = s2;
		this->cap = mydata.size();
	}
};


using dad1 = arrData1<double>;
using iad1 = arrData1<int>;
using dar1 = arrRef1<double>;
using iar1 = arrRef1<int>; 

using dad2 = arrData2<double>;
using iad2 = arrData2<int>;
using dar2 = arrRef2<double>;
using iar2 = arrRef2<int>;

//------------------------------------------------------------------------------------------------------


/* The three functions from Linpack */
void DGEFA(dar2 a, int lda, int n, iar1 ipvt, int info) {
	AutoTimer at(g_timer, "DGEFA"); 
	info = dgefa(a.contiguous(), lda, n, ipvt.contiguous());
}
void DGESL(dar2 a, int lda, int n, iar1 ipvt, dar1 b, int job) {
	AutoTimer at(g_timer, "DGESL"); 
	dgesl(a.contiguous(), lda, n, ipvt.contiguous(),b.contiguous(), job);
}
void DSVDC(dar2 x, int lda, int n, int p, dar1 s, dar1 e,
	dar2 u, int ldu, dar2 v, int ldv, dar1 work, int job, int info) {
	AutoTimer at(g_timer, "DSVDC"); 
	info = dsvdc(x.contiguous(), lda, n, p, s.contiguous(), e.contiguous(),u.contiguous(),
		ldu, v.contiguous(), ldv, work.contiguous(), job);
}


//------------------------------------------------------------------------------------------------------

/* Define callback function types */
using fsub_t  = void (*)(double x, double const z[], double const y[], double f[]);
using dfsub_t = void (*)(double x, double const z[], double const y[], double df[]);
using gsub_t  = void (*)(int i,    double const z[], double& g);
using dgsub_t = void (*)(int i,    double const z[], double dg[]);
using guess_t = void (*)(double x, double const z[], double const y[], double dmval[]);


enum class printMode {
	full = -1,
	selected = 0,
	none = 1
};
enum class meshMode {
	generate = 0, // causes coldae to generate a uniform initial mesh
	custom = 1,   /* if the initial mesh is provided by the user.
					 it is defined in fspace as follows : the mesh
					 aleft = x(1).lt.x(2).lt. ....lt.x(n).lt.x(n + 1) = aright
					 will occupy  fspace(1), ..., fspacen + 1).
					 the user needs to supply only the interior mesh
					 points  fspace(j) = x(j), j = 2, ..., n.*/
	customNoAdaptive = 2 /* if the initial mesh is supplied by the user
						as with ipar(8) = 1, and in addition no adaptive
						mesh selection is to be done.*/
};
enum class guessMode {
	none = 0,  // if no initial guess for the solution is provided
	custom = 1, // if an initial guess is provided by the user in subroutine  guess
	customCoefficients = 2, /*if an initial mesh and approximate solution
							 coefficients are provided by the user in  fspace.
							 (the former and new mesh are the same).*/
	customCoefficientsRefine = 3, /* if a former meshand approximate solution
								coefficients are provided by the user in fspace,
								and the new mesh is to be taken twice as coarse;
								i.e., every second point from the former mesh.*/
	newMesh = 4  /*if in addition to a former initial mesh
					approximate solution coefficients, a new mesh  is provided in fspace as well.
						(see description of output for further details on iguess = 2, 3, and 4.)*/
};
enum class regularControl {
	sensitive = -1, //if the first relax factor is RSTART (use for an extra sensitive nonlinear problem only)
	regular = 0,    // if the problem is regular
	noDamping = 1, /*if the newton iterations are not to be damped
					 (use for initial value problems).*/
	lazy = 2 /* if we are to return immediately upon(a) two
			successive nonconvergences, or (b)after obtaining
			error estimate for the first time.*/
};
enum class indexControl {
	automatic = 0, /* determines the appropriate projection needed at the right end of each
					mesh subinterval using SVD. this is the most expensive and most general option. */
	one = 1,           // if the index of the dae is 1.
	twoHessenberg = 2  //if the index of the dae is 2 and it is in Hessenberg form
};
enum class output_t {
	normal = 1,        // for normal return
	singular = 0,      // if the collocation matrix is singular
	outOfMemory = -1,  // if the expected no. of subintervals exceeds storage specifications
	notConverged = -2, // if the nonlinear iteration has not converged
	inputError = -3    // if there is an input data error.
};

/* System specific properties */
struct systemParams {
	int ncomp;          // number of differential equations (<= 20)
	int ny;             // number of constraints (<= 20)
	iad1 orders ;       // orders of odes
	double left;        // left end of interval
	double right;       // right end of interval
	dad1 bcpoints;      // j-th side condition point (boundary point)

	bool isNonLinear;   // if the problem is nonlinear
	regularControl reg;
	indexControl index; // index of DAE (ignored if ny=0)
};


/* Options for the solver */
struct options {
	int numCollPoints;	   // no. of collocation points per subinterval
	int numSubIntervals;   // no. of subintervals in the initial mesh
	int fdim;              // dimension of fspace
	int idim;              // dimension of ispace
	printMode printLevel;  // output control
	meshMode meshSource;   // mesh control
	guessMode guessSource; // guess control
	
	int numTolerances;     // no. of solution and derivative tolerances
	iad1 ltol;
	dad1 tol;

	int numFixedPoints;    // no. of fixed points in the mesh other than aleft and aright.
	dad1 fixpnt;
};


class cda{

	struct { // COLOUT 
		double PRECIS; int IOUT, IPRINT;
	};
	struct { // COLLOC 
		dad1 RHO = dad1( 7 );
		dad2 COEF = dad2( 7, 7 );
	};
	struct { // COLORD 
		int K; 
		union {
			int NC;
			int NCOMP;
			int NCDUM;
		};
		union {
			int NNY;
			int NYD;
			int NY;
		};
		union {
			int NCY;
			int NCYD;
			int NDM;
		};
		int MSTAR, KD;
		union {
			int KDY;
			int KDUM;
			int KDYM;
		};
		int MMAX;
		iad1 MT = iad1( 20);

		//// aliases
		//auto& M = MT;
	};
	struct { // COLAPR 
		int N, NOLD, NMAX, NZ, NDMZ;
	};
	struct { // COLMSH 
		int MSHFLG, MSHNUM, MSHLMT, MSHALT;
	};
	struct { // COLSID 
		
		dad1 TZETA = dad1( 40);
		union {
			double TLEFT;
			double ALEFT;
		};
		union {
			double TRIGHT;
			double ARIGHT;
		};
		int IZETA;
		union {
			int IDUM;
			int IZSAVE;
		};

		//// aliases
		//auto& ZETA = TZETA;
	};
	struct { // COLNLN 
		int NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX;
	};
	struct { // COLEST 
		dad1 TOL = dad1(40), WGTMSH = dad1(40), WGTERR = dad1(40), TOLIN = dad1(40), ROOT = dad1(40);
		iad1 JTOL = iad1(40), LTOL = iad1(40);
		int NTOL;

		// aliases
		/*auto& LTOL = LTTOL;
		auto& TOL = TTL;*/
	};
	struct { // COLBAS 
		dad2 B = dad2( 7, 4 ), ACOL = dad2( 28, 7 ), ASAVE = dad2( 28, 4 );
	};

	void dumpState() {
		std::ofstream file("state_cpp.txt");
		{
			file << PRECIS << " " << IOUT << " " << IPRINT << std::endl;
		}
		{
			for (int i = 1; i <= 7; ++i)
				file << RHO(i) << " ";
			file << std::endl;
			for (int i = 1; i <= 49; ++i)
				file << COEF.contiguous()[i-1] << " ";
			file << std::endl;
		}
		{
			file << fmt::format("{} {} {} {} {} {} {} {} ", K, NC, NNY, NCY, MSTAR, KD, KDY, MMAX);
			for (int i = 1; i <= 20; ++i)
				file << MT(i) << " ";
			file << std::endl;
		}
		{
			file << fmt::format("{} {} {} {} {} ", N, NOLD, NMAX, NZ, NDMZ);
			file << std::endl;
		}
		{
			file << fmt::format("{} {} {} {} ", MSHFLG, MSHNUM, MSHLMT, MSHALT);
			file << std::endl;
		}
		{
			for (int i = 1; i <= 40; ++i)
				file << TZETA(i) << " ";
			file << std::endl;
			file << fmt::format("{} {} {} {} ", TLEFT, TRIGHT, IZETA, IDUM);
			file << std::endl;
		}
		{
			file << fmt::format("{} {} {} {} {} {} ", NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX);
			file << std::endl;
		}
		{
			file << NTOL << std::endl;
			for (int i = 1; i <= 40; ++i)
				file << TOL(i) << " ";
			file << std::endl;
			for (int i = 1; i <= 40; ++i)
				file << WGTMSH(i) << " ";
			file << std::endl;
			for (int i = 1; i <= 40; ++i)
				file << WGTERR(i) << " ";
			file << std::endl;
			for (int i = 1; i <= 40; ++i)
				file << TOLIN(i) << " ";
			file << std::endl;
			for (int i = 1; i <= 40; ++i)
				file << ROOT(i) << " ";
			file << std::endl;

			for (int i = 1; i <= 40; ++i)
				file << JTOL(i) << " ";
			file << std::endl;
			for (int i = 1; i <= 40; ++i)
				file << LTOL(i) << " ";
			file << std::endl;
		}
		file.close();
	}

public:
void COLDAE(systemParams const& params, options const& opts,
	iar1 ispace, dar1 fspace, output_t& iflag,
	fsub_t fsub, dfsub_t dfsub, gsub_t gsub, dgsub_t dgsub, guess_t guess)
{
	ispace.assertDim(opts.idim);
	fspace.assertDim(opts.fdim);

	this->NCOMP = params.ncomp;
	this->NY = params.ny;
	this->ALEFT = params.left;
	this->ARIGHT = params.right;

	this->LTOL.copyFrom(opts.ltol);
	this->TOL.copyFrom(opts.tol);


	dad1 DUMMY(1), DUMMY2(840);


	//*********************************************************************
	//
	//     the actual subroutine coldae serves as an interface with
	//     the package of subroutines referred to collectively as
	//     coldae. the subroutine serves to test some of the input
	//     parameters, rename some of the parameters (to make under-
	//     standing of the coding easier), to do some initialization,
	//     and to break the work areas fspace and ispace up into the
	//     arrays needed by the program.
	//
	//**********************************************************************

	//  specify machine dependent output unit  iout  and compute machine
	//  dependent constant  precis = 100 * machine unit roundoff


	if (opts.printLevel != printMode::none) 
		fmt::print("VERSION *1* OF COLDAE\n");

	IOUT = 6;
	PRECIS = 1.0;
	double PRECP1;
	do {
		PRECIS = PRECIS / 2.0;
		PRECP1 = PRECIS + 1.0;
	} while (PRECP1 > 1.0);

	PRECIS = PRECIS * 100.0;

	//  in case incorrect input data is detected, the program returns
	//  immediately with iflag=-3.
	iflag = output_t::inputError;
	NCY = NCOMP + NY;
	if (NCOMP < 0 || NCOMP > 20)
		return;
	if (NY < 0 || NY > 20)
		return;
	if (NCY < 1 || NCY > 40)
		return;
	for (int i = 1; i <= NCOMP; ++i)
		if (params.orders(i) < 1 || params.orders(i) > 4)
			return;


	//  rename some of the parameters and set default values.
	NONLIN = params.isNonLinear ? 1 : 0;
	K = opts.numCollPoints;
	N = opts.numSubIntervals;
	if (N == 0)  
		N = 5;
	int IREAD = static_cast<int>(opts.meshSource);
	IGUESS = static_cast<int>(opts.guessSource);
	if (NONLIN == 0 && IGUESS == 1)  IGUESS = 0;
	if (IGUESS >= 2 && IREAD == 0)   IREAD = 1;
	ICARE = static_cast<int>(params.reg);
	NTOL = opts.numTolerances;
	int NDIMF = opts.fdim;
	int NDIMI = opts.idim;
	int NFXPNT = opts.numFixedPoints;
	IPRINT = static_cast<int>(opts.printLevel);
	INDEX = static_cast<int>(params.index);
	if (NY == 0) INDEX = 0;
	MSTAR = 0;
	MMAX = 0;

	for (int i = 1; i <= NCOMP; ++i) {
		MMAX = std::max(MMAX, params.orders(i));
		MSTAR = MSTAR + params.orders(i);
		MT(i) = params.orders(i);
	}
	if (K == 0)   K = std::max(MMAX + 1, 5 - MMAX);
	for (int i = 1; i <= MSTAR; ++i)
		TZETA(i) = params.bcpoints(i);
	for (int i = 1; i <= NTOL; ++i) {
		LTOL(i) = LTOL(i);
		TOLIN(i) = TOL(i);
	}
	TLEFT = ALEFT;
	TRIGHT = ARIGHT;
	NC = NCOMP;
	NNY = NY;
	KD = K * NCOMP;
	KDY = K * NCY;

	//  print the input data for checking.
	if (IPRINT <= -1)
	{
		if (NONLIN == 0) {
			fmt::print("THE NUMBER OF (LINEAR) DIFF EQNS IS {}, THEIR ORDERS ARE {}\n",
				NCOMP, params.orders);
		}
		else {
			fmt::print("THE NUMBER OF (NONLINEAR) DIFF EQNS IS {}, THEIR ORDERS ARE {}\n", 
				NCOMP, params.orders);
		}

		fmt::print("THERE ARE {} ALGEBRAIC CONSTRAINTS\n", NY);
		if (NY > 0 && INDEX == 0) {
			fmt::print("THE PROBLEM HAS MIXED INDEX CONSTRAINTS\n");
		}
		else {
			fmt::print("THE INDEX IS {}\n", INDEX);
		}
		fmt::print("SIDE CONDITION POINTS ZETA: {}\n", params.bcpoints);		
		if (NFXPNT > 0) {
			fmt::print("THERE ARE {} FIXED POINTS IN THE MESH - {}\n",
				NFXPNT, opts.fixpnt);
		}
		fmt::print("NUMBER OF COLLOC PTS PER INTERVAL IS {}\n", K);	
		fmt::print("COMPONENTS OF Z REQUIRING TOLERANCES: {}\n", LTOL);
		fmt::print("CORRESPONDING ERROR TOLERANCES: {}\n", TOL);

		if (IGUESS >= 2) {
			fmt::print("INITIAL MESH(ES) AND Z, DMZ PROVIDED BY USER\n");
		}
		if (IREAD == 2) {
			fmt::print("NO ADAPTIVE MESH SELECTION\n");
		}
	}

	//  check for correctness of data
	if (K < 0 || K > 7)                return;
	if (N < 0)                         return;
	if (IREAD < 0 || IREAD > 2)        return;
	if (IGUESS < 0 || IGUESS > 4)      return;
	if (ICARE < -1 || ICARE > 2)	   return;
	if (INDEX < 0 || INDEX > 2)        return;
	if (NTOL < 0 || NTOL > MSTAR)      return;
	if (NFXPNT < 0)                    return;
	if (IPRINT < (-1) || IPRINT > 1)   return;
	if (MSTAR < 0 || MSTAR > 40)       return;

	int IP = 1;
	for (int i = 1; i <= MSTAR; ++i) {
		if (abs(params.bcpoints(i) - ALEFT) < PRECIS || abs(params.bcpoints(i) - ARIGHT) < PRECIS)
			continue;

		while (true) {
			if (IP > NFXPNT)
				return;
			if (params.bcpoints(i) - PRECIS < opts.fixpnt(IP))
				break;
			IP = IP + 1;
		}

		if (params.bcpoints(i) + PRECIS < opts.fixpnt(IP))
			return;
	}

	//  set limits on iterations and initialize counters.
	//  limit = maximum number of newton iterations per mesh.
	//  see subroutine  newmsh  for the roles of  mshlmt , mshflg , mshnum , and  mshalt .
	MSHLMT = 3;
	MSHFLG = 0;
	MSHNUM = 1;
	MSHALT = 1;
	LIMIT = 40;

	//  compute the maxium possible n for the given sizes of ispace  and  fspace.
	int NREC = 0;
	for (int i = 1; i <= MSTAR; ++i) {
		int IB = MSTAR + 1 - i;
		if (params.bcpoints(IB) >= ARIGHT)  NREC = i;
	}
	int NFIXI = MSTAR;
	int NSIZEI = 3 + KDY + MSTAR;
	int NFIXF = NREC * (2 * MSTAR) + 5 * MSTAR + 3;
	int NSIZEF = 4 + 3 * MSTAR + (KDY + 5) * (KDY + MSTAR) + (2 * MSTAR - NREC) * 2 * MSTAR + (MSTAR + NY + 2) * NCOMP + KDY;
	int NMAXF = (NDIMF - NFIXF) / NSIZEF;
	int NMAXI = (NDIMI - NFIXI) / NSIZEI;
	if (IPRINT < 1) {
		fmt::print("THE MAXIMUM NUMBER OF SUBINTERVALS IS MIN({}(ALLOWED FROM FSPACE),{}(ALLOWED FROM ISPACE))\n",
			NMAXF, NMAXI);
	}
	NMAX = std::min(NMAXF, NMAXI);
	if (NMAX < N)
		return;
	if (NMAX < NFXPNT + 1)   
		return;
	if (NMAX < 2 * NFXPNT + 2 && IPRINT < 1) {
		fmt::print("INSUFFICIENT SPACE TO DOUBLE MESH FOR ERROR ESTIMATE\n");
	}


	//  generate pointers to break up  fspace  and  ispace .
	int LXI = 1;
	int LG = LXI + NMAX + 1;
	int LXIOLD = LG + 2 * MSTAR * (NMAX * (2 * MSTAR - NREC) + NREC);
	int LW = LXIOLD + NMAX + 1;
	int LV = LW + (KDY * KDY) * NMAX;
	int LFC = LV + MSTAR * KDY * NMAX;
	int LZ = LFC + (MSTAR + NY + 2) * NCOMP * NMAX;
	int LDMZ = LZ + MSTAR * (NMAX + 1);
	int LDMV = LDMZ + KDY * NMAX;
	int LDELZ = LDMV + KDY * NMAX;
	int LDELDZ = LDELZ + MSTAR * (NMAX + 1);
	int LDQZ = LDELDZ + KDY * NMAX;
	int LDQDMZ = LDQZ + MSTAR * (NMAX + 1);
	int LRHS = LDQDMZ + KDY * NMAX;
	int LVALST = LRHS + KDY * NMAX + MSTAR;
	int LSLOPE = LVALST + 4 * MSTAR * NMAX;
	int LACCUM = LSLOPE + NMAX;
	int LSCL = LACCUM + NMAX + 1;
	int LDSCL = LSCL + MSTAR * (NMAX + 1);
	int LPVTG = 1;
	int LPVTW = LPVTG + MSTAR * (NMAX + 1);
	int LINTEG = LPVTW + KDY * NMAX;


	//  if  iguess .ge. 2, move  xiold, z, and  dmz  to their proper
	//  locations in  fspace.
	if (IGUESS >= 2) {
		NOLD = N;
		if (IGUESS == 4)
			NOLD = ispace(1);
		NZ = MSTAR * (NOLD + 1);
		NDMZ = KDY * NOLD;
		int NP1 = N + 1;
		if (IGUESS == 4)
			NP1 = NP1 + NOLD + 1;
		for (int i = 1; i <= NZ; ++i)
			fspace(LZ + i - 1) = fspace(NP1 + i);
		int IDMZ = NP1 + NZ;
		for (int i = 1; i <= NDMZ; ++i)
			fspace(LDMZ + i - 1) = fspace(IDMZ + i);
		NP1 = NOLD + 1;
		if (IGUESS == 4) {
			for (int i = 1; i <= NP1; ++i)
				fspace(LXIOLD + i - 1) = fspace(N + 1 + i);
		}
		else {
			for (int i = 1; i <= NP1; ++i)
				fspace(LXIOLD + i - 1) = fspace(LXI + i - 1);
		}
	}

	//dumpState();

	//  initialize collocation points, constants, mesh.
	CONSTS();
	int NYCB;
	if (NY == 0)
		NYCB = 1;
	else
		NYCB = NY;

	int meshmode = 3 + IREAD;
	NEWMSH(meshmode, fspace.sub(LXI), fspace.sub(LXIOLD),
		DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY,
		NFXPNT, opts.fixpnt, DUMMY2, dfsub, DUMMY2, DUMMY2, NYCB);

	//  determine first approximation, if the problem is nonlinear.
	if (IGUESS < 2) {
		for (int i = 1; i <= N + 1; ++i)
			fspace(i + LXIOLD - 1) = fspace(i + LXI - 1);
		NOLD = N;
		if (NONLIN != 0 && IGUESS != 1) {
			//  system provides first approximation of the solution.
			//  choose z(j) = 0  for j=1,...,mstar.
			for (int i = 1; i <= NZ; ++i)
				fspace(LZ - 1 + i) = 0.0;
			for (int i = 1; i <= NDMZ; ++i)
				fspace(LDMZ - 1 + i) = 0.0;
		}
	}
	if (IGUESS >= 2)  
		IGUESS = 0;

	
	CONTRL(fspace.sub(LXI), fspace.sub(LXIOLD), fspace.sub(LZ), fspace.sub(LDMZ), fspace.sub(LDMV),
		fspace.sub(LRHS), fspace.sub(LDELZ), fspace.sub(LDELDZ), fspace.sub(LDQZ),
		fspace.sub(LDQDMZ), fspace.sub(LG), fspace.sub(LW), fspace.sub(LV), fspace.sub(LFC),
		fspace.sub(LVALST), fspace.sub(LSLOPE), fspace.sub(LSCL), fspace.sub(LDSCL),
		fspace.sub(LACCUM), ispace.sub(LPVTG), ispace.sub(LINTEG), ispace.sub(LPVTW),
		NFXPNT, opts.fixpnt, iflag, fsub, dfsub, gsub, dgsub, guess);

	//  prepare output
	ispace(1) = N;
	ispace(2) = K;
	ispace(3) = NCOMP;
	ispace(4) = NY;
	ispace(5) = MSTAR;
	ispace(6) = MMAX;
	ispace(7) = NZ + NDMZ + N + 2;
	int K2 = K * K;
	ispace(8) = ispace(7) + K2 - 1;
	for (int i = 1; i <= NCOMP; ++i)
		ispace(8 + i) = params.orders(i);
	for (int i = 1; i <= NZ; ++i)
		fspace(N + 1 + i) = fspace(LZ - 1 + i);
	int IDMZ = N + 1 + NZ;

	for (int i = 1; i <= NDMZ; ++i)
		fspace(IDMZ + i) = fspace(LDMZ - 1 + i);
	int IC = IDMZ + NDMZ;
	for (int i = 1; i <= K2; ++i)
		fspace(IC + i) = COEF.contiguous()[i-1];
}

// * ****************************************************************
//
//     purpose
//
//           set up a standard call to  approx  to evaluate the
//           approximate solution  z = z(u(x)), y = y(x)  at a
//           point x(it has been computed by a call to  coldae).
//           the parameters needed for  approx  are retrieved
//           from the work arrays  ispaceand fspace .
//
//*****************************************************************
void APPSLN(double& X, dar1 Z, dar1 Y, dar1 FSPACE, iar1 ISPACE) {
	Z.assertDim(1);
	Y.assertDim(1);
	FSPACE.assertDim(1);
	ISPACE.assertDim(1);

	dad1 A(28);

	int IS6 = ISPACE(7);
	int IS5 = ISPACE(1) + 2;
	int IS4 = IS5 + ISPACE(5) * (ISPACE(1) + 1);
	int i = 1;
	APPROX(i, X, Z.contiguous(), Y.contiguous(), A.contiguous(), FSPACE.sub(IS6).contiguous(),
		FSPACE.sub(1).contiguous(), ISPACE(1),
		FSPACE.sub(IS5).contiguous(), FSPACE.sub(IS4).contiguous(), ISPACE(2), ISPACE(3),
		ISPACE(4), ISPACE(6), ISPACE.sub(9).contiguous(), ISPACE(5), 2, nullptr, 1);
}



private:

//**********************************************************************
//
//   purpose
//     this subroutine is the actual driver.  the nonlinear iteration
//     strategy is controlled here ( see [6] ). upon convergence, errchk
//     is called to test for satisfaction of the requested tolerances.
//
//   variables
//
//     check  - maximum tolerance value, used as part of criteria for
//              checking for nonlinear iteration convergence
//     relax  - the relaxation factor for damped newton iteration
//     relmin - minimum allowable value for relax  (otherwise the
//              jacobian is considered singular).
//     rlxold - previous relax
//     rstart - initial value for relax when problem is sensitive
//     ifrz   - number of fixed jacobian iterations
//     lmtfrz - maximum value for ifrz before performing a reinversion
//     iter   - number of iterations (counted only when jacobian
//              reinversions are performed).
//     xi     - current mesh
//     xiold  - previous mesh
//     ipred  = 0  if relax is determined by a correction
//            = 1  if relax is determined by a prediction
//     ifreez = 0  if the jacobian is to be updated
//            = 1  if the jacobian is currently fixed (frozen)
//     iconv  = 0  if no previous convergence has been obtained
//            = 1  if convergence on a previous mesh has been obtained
//     icare  =-1  no convergence occurred (used for regular problems)
//            = 0  a regular problem
//            = 1  no damped newton
//            = 2  used for continuation (see description of ipar(10)
//                 in coldae).
//     rnorm  - norm of rhs (right hand side) for current iteration
//     rnold  - norm of rhs for previous iteration
//     anscl  - scaled norm of newton correction
//     anfix  - scaled norm of newton correction at next step
//     anorm  - scaled norm of a correction obtained with jacobian fixed
//     nz     - number of components of  z  (see subroutine approx)
//     ndmz   - number of components of  dmz  (see subroutine approx)
//     imesh  - a control variable for subroutines newmsh and errchk
//            = 1  the current mesh resulted from mesh selection
//                 or is the initial mesh.
//            = 2  the current mesh resulted from doubling the
//                 previous mesh
//
//**********************************************************************
void CONTRL(dar1 XI, dar1 XIOLD, dar1 Z, dar1 DMZ, dar1 DMV, dar1 RHS, dar1 DELZ, dar1 DELDMZ,
	dar1 DQZ, dar1 DQDMZ, dar1 G, dar1 W, dar1 V, dar1 FC, dar1 VALSTR, dar1 SLOPE, dar1 SCALE, dar1 DSCALE,
	dar1 ACCUM, iar1 IPVTG, iar1 INTEGS, iar1 IPVTW, const int NFXPNT, dar1 FIXPNT, output_t& iflag,
	fsub_t fsub, dfsub_t dfsub, gsub_t gsub, dgsub_t dgsub, guess_t guess)
{
	XI.assertDim(1);
	XIOLD.assertDim(1);
	Z.assertDim(1);
	DMZ.assertDim(1);
	RHS.assertDim(1);
	DMV.assertDim(1);

	G.assertDim(1);
	W.assertDim(1);
	V.assertDim(1);
	VALSTR.assertDim(1);
	SLOPE.assertDim(1);
	ACCUM.assertDim(1);

	DELZ.assertDim(1);
	DELDMZ.assertDim(1);
	DQZ.assertDim(1);
	DQDMZ.assertDim(1);
	FIXPNT.assertDim(1);

	SCALE.assertDim(1);
	DSCALE.assertDim(1);
	FC.assertDim(1);
	INTEGS.assertDim(1);
	IPVTG.assertDim(1);
	IPVTW.assertDim(1);
	
	
	//using namespace COLOUT; // double PRECIS; int IOUT, IPRINT; 
	//using namespace COLORD; // int K, NC, NNY, NCY, MSTAR, KD, KDY, MMAX; iad1 MT(20);
	//using namespace COLAPR; // int N, NOLD, NMAX, NZ, NDMZ;
	//using namespace COLMSH; // int MSHFLG, MSHNUM, MSHLMT, MSHALT;
	//using namespace COLSID; // dad1 TZETA(40);  double TLEFT, TRIGHT;  int IZETA, IDUM;
	//using namespace COLNLN; // int NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX;
	//using namespace COLEST; // dad1 TTL(40), WGTMSH(40), WGTERR(40), TOLIN(40), ROOT(40);  iad1 JTOL(40), LTTOL(40); int NTOL;
	//
	//dumpState();

	dad1 DUMMY(1), DF(800);
	dad2 FCSP(NCOMP, 60), CBSP(20, 20);
	double RNORM, RNOLD;

	// constants for control of nonlinear iteration
	double RELMIN = 1.e-3;
	double RSTART = 1.e-2;
	double LMTFRZ = 4;

	// compute the maximum tolerance
	double CHECK = 0.0;
	for (int i = 1; i <= NTOL; ++i)
		CHECK = std::max(TOLIN(i), CHECK);
	int IMESH = 1;
	int ICONV = 0;
	if (NONLIN == 0) ICONV = 1;
	int ICOR = 0;
	int NOCONV = 0;
	int MSING = 0;
	int ISING = 0;

	// TODO: make lokal
	double RLXOLD;
	int IPRED, IFREEZ;
	int IFRZ;
	double RELAX;
	double ANORM = 0.0;
	double ANFIX = 0.0;

	
	//  the main iteration begins here .
	//  loop 20 is executed until error tolerances are satisfied or
	//  the code fails (due to a singular matrix or storage limitations)

	while (true) {
		{
			//fmt::print(fg(fmt::color::cornflower_blue), "while:\n");

			// initialization for a new mesh
			ITER = 0;
			if (NONLIN <= 0) {

				// the linear case.
				// set up and solve equations
				LSYSLV(MSING, XI, XIOLD, DUMMY, DUMMY, Z, DMZ, G,
					W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 0,
					fsub, dfsub, gsub, dgsub, guess, ISING);

				// check for a singular matrix
				if (ISING != 0) {
					if (IPRINT < 1) {
						fmt::print("SINGULAR PROJECTION MATRIX DUE TO INDEX > 2\n"); 
					}
					iflag = output_t::singular;
					return;
				}
				if (MSING == 0)
					goto n400;
			n30:
				//fmt::print(fg(fmt::color::cornflower_blue), "30:\n");

				if (MSING >= 0) {
					if (IPRINT < 1) {
						fmt::print(fg(fmt::color::red), "A LOCAL ELIMINATION MATRIX IS SINGULAR\n"); 
					}
					goto n460;
				}
				if (IPRINT < 1) {
					fmt::print(fg(fmt::color::red), "THE GLOBAL BVP - MATRIX IS SINGULAR\n");
				}
				iflag = output_t::singular;
				return;
			}

			// iteration loop for nonlinear case
			// define the initial relaxation parameter (= relax)
			RELAX = 1.0;

			// check for previous convergence and problem sensitivity
			if (ICARE == -1)
				RELAX = RSTART;
			if (ICARE == 1)
				RELAX = 1.0;
			if (ICONV == 0)
				goto n160;

			// convergence on a previous mesh has been obtained.    thus
			// we have a very good initial approximation for the newton
			// process.    proceed with one full newton and then iterate
			// with a fixed jacobian.
			IFREEZ = 0;

			// evaluate right hand side and its norm  and
			// find the first newton correction
			LSYSLV(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
				W, V, FC, RHS, DQDMZ, INTEGS, IPVTG, IPVTW, RNOLD, 1,
				fsub, dfsub, gsub, dgsub, guess, ISING);

			if (IPRINT < 0) {
				fmt::print("FIXED JACOBIAN ITERATIONS\n");
			}
			if (IPRINT < 0) {
				fmt::print("ITERATION = {}, NORM(RHS) = {}\n", ITER, RNOLD);
			}
			goto n70;
		}
	n60:
		{
			//fmt::print(fg(fmt::color::cornflower_blue), "60:\n");
			// solve for the next iterate .
			// the value of ifreez determines whether this is a full
			// newton step (=0) or a fixed jacobian iteration (=1).
			if (IPRINT < 0) {
				fmt::print("ITERATION = {}, NORM(RHS) = {}\n", ITER, RNORM);
			}
			RNOLD = RNORM;
			LSYSLV(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
				W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM,
				3 + IFREEZ, fsub, dfsub, gsub, dgsub, guess, ISING);
		}
	n70:
		{
			//fmt::print(fg(fmt::color::cornflower_blue), "70:\n");

			// check for a singular matrix
			if (MSING != 0)
				goto n30;
			if (ISING != 0) {
				if (IPRINT < 1) {
					fmt::print(fg(fmt::color::red), "SINGULAR PROJECTION MATRIX DUE TO INDEX > 2\n");
				}
				iflag = output_t::singular;
				return;
			}
			if (IFREEZ != 1) {
				// a full newton step
				ITER = ITER + 1;
				IFRZ = 0;
			}

			// update   z and dmz , compute new  rhs  and its norm
			for (int i = 1; i <= NZ; ++i)
				Z(i) = Z(i) + DELZ(i);
			for (int i = 1; i <= NDMZ; ++i)
				DMZ(i) = DMZ(i) + DELDMZ(i);

			LSYSLV(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
				W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 2,
				fsub, dfsub, gsub, dgsub, guess, ISING);

			//       check monotonicity. if the norm of  rhs  gets smaller,
			//       proceed with a fixed jacobian; else proceed cautiously,
			//       as if convergence has not been obtained before (iconv=0).
			if (RNORM < PRECIS)
				goto n390;
			if (RNORM > RNOLD)
				goto n130;
			if (IFREEZ == 1)
				goto n110;
			IFREEZ = 1;
			goto n60;
		}
	n110:
		{
			//fmt::print(fg(fmt::color::cornflower_blue), "110:\n");

			// verify that the linear convergence with fixed jacobian
			// is fast enough.
			IFRZ = IFRZ + 1;
			if (IFRZ >= LMTFRZ)       IFREEZ = 0;
			if (RNOLD < 4.0 * RNORM)  IFREEZ = 0;

			// check convergence (iconv = 1).
			for (int IT = 1; IT <= NTOL; ++IT) {
				int INZ = LTOL(IT);
				for (int IZ = INZ; IZ <= NZ; IZ += MSTAR) {
					if (abs(DELZ(IZ)) > TOLIN(IT) * (abs(Z(IZ)) + 1.0))
						goto n60;
				}
			}

			// convergence obtained
			if (IPRINT < 1) 
				fmt::print("CONVERGENCE AFTER {} ITERATIONS\n", ITER);
			goto n400;
		}
	n130:
		{
			//fmt::print(fg(fmt::color::cornflower_blue), "130:\n");

			// convergence of fixed jacobian iteration failed.
			if (IPRINT < 0) 
				fmt::print("ITERATION = {}, NORM(RHS) = {}\nSWITCH TO DAMPED NEWTON ITERATION\n", ITER, RNORM);
			ICONV = 0;
			if (ICARE != 1)  
				RELAX = RSTART;
			for (int i = 1; i <= NZ; ++i)
				Z(i) = Z(i) - DELZ(i);
			for (int i = 1; i <= NDMZ; ++i)
				DMZ(i) = DMZ(i) - DELDMZ(i);


			// update old mesh
			for (int i = 1; i <= N + 1; ++i)
				XIOLD(i) = XI(i);
			NOLD = N;
			ITER = 0;
		}
	n160:
		{
			//fmt::print(fg(fmt::color::cornflower_blue), "160:\n");
			// no previous convergence has been obtained. proceed
			// with the damped newton method.
			// evaluate rhs and find the first newton correction.
			if (IPRINT < 0) 
				fmt::print("FULL DAMPED NEWTON ITERATION\n");
			LSYSLV(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
				W, V, FC, RHS, DQDMZ, INTEGS, IPVTG, IPVTW, RNOLD, 1,
				fsub, dfsub, gsub, dgsub, guess, ISING);

			// check for a singular matrix
			if (MSING != 0)
				goto n30;
			if (ISING != 0) {
				if (IPRINT < 1)  fmt::print(fg(fmt::color::red), "SINGULAR PROJECTION MATRIX DUE TO INDEX > 2\n");
				iflag = output_t::singular;
				return;
			}

			// bookkeeping for first mesh
			if (IGUESS == 1)
				IGUESS = 0;

			// find initial scaling
			SKALE(Z, DMZ, XI, SCALE, DSCALE);
			
			RLXOLD = RELAX;
			IPRED = 1;

			goto n220;
		}
	n170:
		{
			//fmt::print(fg(fmt::color::cornflower_blue), "170:\n");
			// main iteration loop
			RNOLD = RNORM;
			if (ITER >= LIMIT)
				goto n430;

			// update scaling
			SKALE(Z, DMZ, XI, SCALE, DSCALE);

			// compute norm of newton correction with new scaling
			double ANSCL = 0.0;
			for (int i = 1; i <= NZ; ++i)
				ANSCL = ANSCL + pow(DELZ(i) * SCALE(i), 2);

			for (int i = 1; i <= NDMZ; ++i)
				ANSCL = ANSCL + pow(DELDMZ(i) * DSCALE(i), 2);

			ANSCL = sqrt(ANSCL / double(NZ + NDMZ));

			// find a newton direction
			LSYSLV(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
				W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 3,
				fsub, dfsub, gsub, dgsub, guess, ISING);

			// check for a singular matrix
			if (MSING != 0)
				goto n30;

			if (ISING != 0) {
				if (IPRINT < 1)  fmt::print(fg(fmt::color::red), "SINGULAR PROJECTION MATRIX DUE TO INDEX > 2\n");
				iflag = output_t::singular;
				return;
			}

			// predict relaxation factor for newton step.
			if (ICARE != 1) {
				double ANDIF = 0.0;
				for (int i = 1; i <= NZ; ++i)
					ANDIF = ANDIF + pow((DQZ(i) - DELZ(i)) * SCALE(i), 2);

				for (int i = 1; i <= NDMZ; ++i)
					ANDIF = ANDIF + pow((DQDMZ(i) - DELDMZ(i)) * DSCALE(i), 2);

				ANDIF = sqrt(ANDIF / double(NZ + NDMZ) + PRECIS);
				RELAX = RELAX * ANSCL / ANDIF;
				if (RELAX > 1.0)  RELAX = 1.0;
				RLXOLD = RELAX;
				IPRED = 1;
			}
		}
	n220:
		{
			//fmt::print(fg(fmt::color::cornflower_blue), "220:\n");
			ITER = ITER + 1;

			// determine a new  z and dmz  and find new  rhs  and its norm
			for (int i = 1; i <= NZ; ++i)
				Z(i) = Z(i) + RELAX * DELZ(i);

			for (int i = 1; i <= NDMZ; ++i)
				DMZ(i) = DMZ(i) + RELAX * DELDMZ(i);
		}
	n250:
		{
			//fmt::print(fg(fmt::color::cornflower_blue), "250:\n");

			/*for (int i = 1; i <= NZ; ++i) {
				fmt::print("+ DQZ(i)={:10.6f}\n", DQZ(i));
			}*/

			LSYSLV(MSING, XI, XIOLD, Z, DMZ, DQZ, DQDMZ, G,
				W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 2,
				fsub, dfsub, gsub, dgsub, guess, ISING);

			/*for (int i = 1; i <= NZ; ++i) {
				fmt::print("++ DQZ(i)={:10.6f}\n", DQZ(i));
			}*/

			// compute a fixed jacobian iterate (used to control relax)
			LSYSLV(MSING, XI, XIOLD, Z, DMZ, DQZ, DQDMZ, G,
				W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 4,
				fsub, dfsub, gsub, dgsub, guess, ISING);

			//for (int i = 1; i <= NZ; ++i) {
			//	fmt::print("+++ DQZ(i)={:10.6f}\n", DQZ(i));
			//}

			// find scaled norms of various terms used to correct relax
			ANORM = 0.0;
			ANFIX = 0.0;
			for (int i = 1; i <= NZ; ++i) {
				//fmt::print("---- DQZ(i)={:10.6f}, SCALE(i)={:10.6f}\n", DQZ(i), SCALE(i));
				ANORM = ANORM + pow(DELZ(i) * SCALE(i), 2);
				ANFIX = ANFIX + pow(DQZ(i) * SCALE(i), 2);
			}
			for (int i = 1; i <= NDMZ; ++i) {
				//fmt::print("---- DQDMZ(i)={:10.6f}, DSCALE(i)={:10.6f}\n", DQDMZ(i), DSCALE(i));
				ANORM = ANORM + pow(DELDMZ(i) * DSCALE(i), 2);
				ANFIX = ANFIX + pow(DQDMZ(i) * DSCALE(i), 2);
			}
			ANORM = sqrt(ANORM / double(NZ + NDMZ));
			ANFIX = sqrt(ANFIX / double(NZ + NDMZ));
			if (ICOR == 1) {
				if (IPRINT < 0)
					fmt::print( "RELAXATION FACTOR CORRECTED TO RELAX = {}\n"
								" NORM OF SCALED RHS CHANGES FROM {} TO {}\n"
								" NORM OF        RHS CHANGES FROM {} TO {}\n",
						RELAX, ANORM, ANFIX, RNOLD, RNORM);
			}
			else {
				if (IPRINT < 0)
					fmt::print("ITERATION = {}  RELAXATION FACTOR = {}\n"
						" NORM OF SCALED RHS CHANGES FROM {} TO {}\n"
						" NORM OF        RHS CHANGES FROM {} TO {}\n",
						ITER, RELAX, ANORM, ANFIX, RNOLD, RNORM);
			}

			ICOR = 0;

			// check for monotonic decrease in  delz and deldmz.
			if (ANFIX < PRECIS || RNORM < PRECIS)
				goto n390;

			{
				if (ANFIX <= ANORM || ICARE == 1) {
					// we have a decrease.
					// if  dqz  and dqdmz  small, check for convergence
					if (ANFIX <= CHECK)
						goto n350;

					// correct the predicted  relax  unless the corrected
					// value is within 10 percent of the predicted one.
					if (IPRED != 1)
						goto n170;
				}
				if (ITER >= LIMIT)
					goto n430;
				if (ICARE == 1)
					goto n170;
			}
			{
				// correct the relaxation factor.
				IPRED = 0;
				double ARG = (ANFIX / ANORM - 1.0) / RELAX + 1.0;
				if (ARG < 0.0)
					goto n170;
				if (ARG > .25 * RELAX + .125 * RELAX * RELAX) {
					double FACTOR = -1.0 + sqrt(1.0 + 8.0 * ARG);
					if (abs(FACTOR - 1.0) < .10 * FACTOR)
						goto n170;
					if (FACTOR < 0.50)
						FACTOR = 0.5;
					RELAX = RELAX / FACTOR;
				}
				else {
					if (RELAX >= .9)
						goto n170;
					RELAX = 1.0;
				}
			}
			ICOR = 1;
			if (RELAX >= RELMIN) {
				{
					double FACT = RELAX - RLXOLD;
					for (int i = 1; i <= NZ; ++i)
						Z(i) = Z(i) + FACT * DELZ(i);

					for (int i = 1; i <= NDMZ; ++i) {
						DMZ(i) = DMZ(i) + FACT * DELDMZ(i);
					}
					RLXOLD = RELAX;
					goto n250;
				}			
			n350: 
				{
					//fmt::print(fg(fmt::color::cornflower_blue), "350:\n");
					// check convergence (iconv = 0).
					for (int IT = 1; IT <= NTOL; ++IT) {
						int INZ = LTOL(IT);
						for (int IZ = INZ; IZ <= NZ; IZ += MSTAR) {
							if (abs(DQZ(IZ)) > TOLIN(IT) * (abs(Z(IZ)) + 1.0))
								goto n170;
						}
					}

					// convergence obtained
					if (IPRINT < 1) 
						fmt::print("CONVERGENCE AFTER {} ITERATIONS\n", ITER);

					// since convergence obtained, update  z and dmz  with term
					// from the fixed jacobian iteration.
					for (int i = 1; i <= NZ; ++i)
						Z(i) = Z(i) + DQZ(i);

					for (int i = 1; i <= NDMZ; ++i)
						DMZ(i) = DMZ(i) + DQDMZ(i);
				}
			n390:
				{
					//fmt::print(fg(fmt::color::cornflower_blue), "390:\n");
					if ((ANFIX < PRECIS || RNORM < PRECIS) && IPRINT < 1) 
						fmt::print("CONVERGENCE AFTER {} ITERATIONS\n", ITER);
					ICONV = 1;
					if (ICARE == -1)  ICARE = 0;
				}
			n400:
				{
					//fmt::print(fg(fmt::color::cornflower_blue), "400:\n");
					// if full output has been requested, print values of the
					// solution components   z  at the meshpoints and  y  at
					// collocation points.
					if (IPRINT >= 0)
						goto n420;
					for (int j = 1; j <= MSTAR; ++j) {
						fmt::print("MESH VALUES FOR Z({}): ", j);
						for (int LJ = j; LJ <= NZ; LJ += MSTAR)
							fmt::print("{:.4}, ", Z(LJ));
						fmt::print("\n");
					}
					for (int j = 1; j <= NY; ++j) {	
						fmt::print("VALUES AT 1st COLLOCATION POINTS FOR Y({}): ", j);
						for (int LJ = j + NCOMP; LJ <= NDMZ; LJ += KDY)
							fmt::print("{:.4}, ", DMZ(LJ));
						fmt::print("\n");
					}
				}
			n420:
				{
					//fmt::print(fg(fmt::color::cornflower_blue), "420:\n");
					// check for error tolerance satisfaction
					int IFIN = 1;
					if (IMESH == 2)
						ERRCHK(XI, Z, DMZ, VALSTR, IFIN);
					if (IMESH == 1 || IFIN == 0 && ICARE != 2)
						goto n460;
					iflag = output_t::normal;
					return;
				}
			n430: ;
				//fmt::print(fg(fmt::color::cornflower_blue), "430:\n");

				// diagnostics for failure of nonlinear iteration.
				if (IPRINT < 1)  
					fmt::print("NO CONVERGENCE AFTER {} ITERATIONS\n", ITER);
				
			}
			else {
				if (IPRINT < 1) {
					fmt::print(fg(fmt::color::orange),
						"NO CONVERGENCE. RELAXATION FACTOR = {} IS TOO SMALL (LESS THAN {})\n", RELAX, RELMIN);
				}
			}

			iflag = output_t::notConverged;
			NOCONV = NOCONV + 1;
			if (ICARE == 2 && NOCONV > 1)
				return;
			if (ICARE == 0)
				ICARE = -1;
		}
	n460:
		{
		//fmt::print(fg(fmt::color::cornflower_blue), "460:\n");
			// update old mesh
			for (int i = 1; i <= N + 1; ++i)
				XIOLD(i) = XI(i);
			NOLD = N;

			// pick a new mesh
			// check safeguards for mesh refinement
			IMESH = 1;
			if (ICONV == 0 || MSHNUM >= MSHLMT || MSHALT >= MSHLMT)
				IMESH = 2;
			if (MSHALT >= MSHLMT && MSHNUM < MSHLMT)
				MSHALT = 1;
			
			int NYCB;
			if (NY == 0)
				NYCB = 1;
			else
				NYCB = NY;


			NEWMSH(IMESH, XI, XIOLD, Z, DMZ,
				DMV, VALSTR,
				SLOPE, ACCUM, NFXPNT,
				FIXPNT, DF, dfsub,
				FCSP, CBSP, NYCB);

			// exit if expected n is too large (but may try n=nmax once)
			if (N > NMAX){
				N = N / 2;
				iflag = output_t::outOfMemory;
				if (ICONV == 0 && IPRINT < 1)
					fmt::print("NO CONVERGENCE\n"); 
				if (ICONV == 1 && IPRINT < 1)  
					fmt::print("PROBABLY TOLERANCES TOO STRINGENT OR NMAX TOO SMALL\n");
				return;
			}
			if (ICONV == 0) 
				IMESH = 1;
			if (ICARE == 1)  ICONV = 0;
		}
	}
}


//**********************************************************************
//
//   purpose
//            provide a proper scaling of the state variables, used
//            to control the damping factor for a newton iteration [4].
//
//   variables
//
//            n      = number of mesh subintervals
//            mstar  = number of unknomns in z(u(x))
//            kdy     = number of unknowns in dmz per mesh subinterval
//            z      = the global current solution vector
//            dmz    = the global current highest derivs vector
//            xi     = the current mesh
//            scale  = scaling vector for z
//            dscale = scaling vector for dmz
//
//**********************************************************************
void SKALE(dar2 Z, dar2 DMZ, dar1 XI, dar2 SCALE, dar2 DSCALE)
{
	Z.reshape(MSTAR, 1);
	Z.assertDim(MSTAR, 1);
	SCALE.reshape(MSTAR, 1); 
	SCALE.assertDim(MSTAR, 1);
	DMZ.reshape(KDY, N); 
	DMZ.assertDim(KDY, N);
	DSCALE.reshape(KDY, N); 
	DSCALE.assertDim(KDY, N);
	XI.assertDim(1);
	
	//using namespace COLORD; // int K, NC, NNY, NCY, MSTAR, KD, KDY, MMAX; iad1 MT(20);

	//N = N;
	//MSTAR = MSTAR;
	//KDY = KDY;

	dad1 BASM(5);

	BASM(1) = 1.0;
	for (int j = 1; j <= N; ++j) {
		int IZ = 1;
		double H = XI(j + 1) - XI(j);
		for (int l = 1; l <= MMAX; ++l)
			BASM(l + 1) = BASM(l) * H / double(l);

		for (int ICOMP = 1; ICOMP <= NCOMP; ++ICOMP) {
			double SCAL = (abs(Z(IZ, j)) + abs(Z(IZ, j + 1))) * .5 + 1.0;
			int MJ = MT(ICOMP);
			for (int l = 1; l <= MJ; ++l) {
				SCALE(IZ, j) = BASM(l) / SCAL;
				IZ = IZ + 1;
			}
			SCAL = BASM(MJ + 1) / SCAL;
			for (int IDMZ = ICOMP; IDMZ <= KDY; IDMZ += NCY)
				DSCALE(IDMZ, j) = SCAL;

		}
		for (int ICOMP = 1 + NCOMP; ICOMP <= NCY; ++ICOMP) {
			double SCAL = 1.0 / (abs(DMZ(ICOMP, j)) + 1.0);
			for (int IDMZ = ICOMP; IDMZ <= KDY; IDMZ += NCY) {
				DSCALE(IDMZ, j) = SCAL;
			}
		}
	}
	int NP1 = N + 1;
	for (int IZ = 1; IZ <= MSTAR; ++IZ)
		SCALE(IZ, NP1) = SCALE(IZ, N);
}




//----------------------------------------------------------------------
//                            p a r t  2
//          mesh selection, error estimation, (and related
//          constant assignment) routines -- see [5], [6]
//----------------------------------------------------------------------




//**********************************************************************
//
//   purpose
//            select a mesh on which a collocation solution is to be
//            determined
//
//                           there are 5 possible modes of action:
//            mode = 5,4,3 - deal mainly with definition of an initial
//                           mesh for the current boundary value problem
//                 = 2,1   - deal with definition of a new mesh, either
//                           by simple mesh halving or by mesh selection
//            more specifically, for
//            mode = 5  an initial (generally nonuniform) mesh is
//                      defined by the user and no mesh selection is to
//                      be performed
//                 = 4  an initial (generally nonuniform) mesh is
//                      defined by the user
//                 = 3  a simple uniform mesh (except possibly for some
//                      fixed points) is defined; n= no. of subintervals
//                 = 1  the automatic mesh selection procedure is used
//                      (see [5] for details)
//                 = 2  a simple mesh halving is performed
//
//**********************************************************************
//
//   variables
//
//            n      = number of mesh subintervals
//            nold   = number of subintervals for former mesh
//            xi     - mesh point array
//            xiold  - former mesh point array
//            mshlmt - maximum no. of mesh selections which are permitted
//                     for a given n before mesh halving
//            mshnum - no. of mesh selections which have actually been
//                     performed for the given n
//            mshalt - no. of consecutive times ( plus 1 ) the mesh
//                     selection has alternately halved and doubled n.
//                     if mshalt .ge. mshlmt then  contrl  requires
//                     that the current mesh be halved.
//            mshflg = 1  the mesh is a halving of its former mesh
//                       (so an error estimate has been calculated)
//                   = 0  otherwise
//            iguess - ipar(9) in subroutine coldae.  it is used
//                     here only for mode=5 and 4, where
//                   = 2 the subroutine sets xi=xiold.  this is
//                       used e.g. if continuation is being per-
//                       formed, and a mesh for the old differen-
//                       tial equation is being used
//                   = 3 same as for =2, except xi uses every other
//                       point of xiold (so mesh xiold is mesh xi
//                       halved)
//                   = 4 xi has been defined by the user, and an old
//                       mesh xiold is also available
//                       otherwise, xi has been defined by the user
//                       and we set xiold=xi in this subroutine
//            slope  - an approximate quantity to be equidistributed for
//                     mesh selection (see [5]), viz,
//                             .                        (k+mj)
//                     slope(i)=     max   (weight(l) *u      (xi(i)))
//                               1.le.l.le.ntol         j
//
//                     where j=jtol(l)
//            slphmx - maximum of slope(i)*(xiold(i+1)-xiold(i)) for
//                     i = 1 ,..., nold.
//            accum  - accum(i) is the integral of  slope  from  aleft
//                     to  xiold(i).
//            valstr - is assigned values needed in  errchk  for the
//                     error estimate.
//            fc     - you know
//**********************************************************************
void NEWMSH(int& MODE, dar1 XI, dar1 XIOLD, dar1 Z, dar1 DMZ, dar1 DMV,
 dar1 VALSTR, dar1 SLOPE, dar1 ACCUM, const int NFXPNT, dar1 FIXPNT,
 dar2 DF, dfsub_t dfsub, dar2 FC, dar2 CB, const int NYCB)
{
	//using namespace COLOUT; // double PRECIS; int IOUT, IPRINT; 
	//using namespace COLLOC; // dad1 RHO(7), COEF(49);
	//using namespace COLORD; // int K, NC, NNY, NCY, MSTAR, KD, KDY, MMAX; iad1 MT(20);
	//using namespace COLAPR; // int N, NOLD, NMAX, NZ, NDMZ;
	//using namespace COLMSH; // int MSHFLG, MSHNUM, MSHLMT, MSHALT;
	//using namespace COLSID; // dad1 TZETA(40);  double TLEFT, TRIGHT;  int IZETA, IDUM;
	//using namespace COLNLN; // int NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX;
	//using namespace COLEST; // dad1 TTL(40), WGTMSH(40), WGTERR(40), TOLIN(40), ROOT(40);  iad1 JTOL(40), LTTOL(40); int NTOL;
	//using namespace COLBAS; // dad2 B(7, 4), ACOL(28, 7), ASAVE(28, 4);

	SLOPE.assertDim(1);
	ACCUM.assertDim(1);
	XI.assertDim(1);
	XIOLD.assertDim(NOLD + 1); // simon
	Z.assertDim(1);
	DMZ.assertDim(1);
	DMV.assertDim(1);
	FIXPNT.assertDim(1);
	VALSTR.assertDim(1);
	FC.assertDim(NCOMP, 60);
	DF.assertDim(NCY, 1);
	CB.assertDim(NYCB, NYCB);

	
	dad1 D1(40), D2(40), DUMMY(1), ZVAL(40), YVAL(40), A(28), BCOL(40), U(400), V(400);
	iad1 IPVTCB(40);

	int ISING = 0;

	int NFXP1 = NFXPNT + 1;
	switch (MODE) {
	case 5:
		// mode=5   set mshlmt=1 so that no mesh selection is performed
		MSHLMT = 1;

		[[fallthrough]];
	case 4: {
		//  mode=4   the user-specified initial mesh is already in place.
		if (IGUESS >= 2) {
			//  iguess=2, 3 or 4.	
			if (IPRINT < 1)
				fmt::print("THE FORMER MESH (OF {} SUBINTERVALS): {}\n", NOLD, XIOLD);
			if (IGUESS == 3) {
				//  if iread ( ipar(8) ) .ge. 1 and iguess ( ipar(9) ) .eq. 3
				//  then the first mesh is every second point of the
				//  mesh in  xiold .
				N = NOLD / 2;
				for (int j = 1, i = 0; j <= NOLD; j += 2) {
					i = i + 1;
					XI(i) = XIOLD(j);
				}
			}
		}
		XI(1) = TLEFT;
		XI(N + 1) = TRIGHT;
		break;
	}
	case 3: {
		//  mode=3   generate a (piecewise) uniform mesh. if there are
		//  fixed points then ensure that the n being used is large enough
		if (N < NFXP1)
			N = NFXP1;
		int NP1 = N + 1;
		XI(1) = TLEFT;
		int ILEFT = 1;
		double XLEFT = TLEFT;

		//  loop over the subregions between fixed points.
		for (int j = 1; j <= NFXP1; ++j) {
			double XRIGHT = TRIGHT;
			int IRIGHT = NP1;
			if (j != NFXP1) {
				XRIGHT = FIXPNT(j);

				// determine where the j-th fixed point should fall in the
				// new mesh - this is xi(iright) and the (j-1)st fixed
				// point is in xi(ileft)
				int NMIN = int((XRIGHT - TLEFT) / (TRIGHT - TLEFT) * double(N) + 1.5);
				if (NMIN > N - NFXPNT + j)
					NMIN = N - NFXPNT + j;
				IRIGHT = std::max(ILEFT + 1, NMIN);
			}
			XI(IRIGHT) = XRIGHT;

			// generate equally spaced points between the j-1st and the
			// j-th fixed points.
			int NREGN = IRIGHT - ILEFT - 1;
			if (NREGN != 0) {
				double DX = (XRIGHT - XLEFT) / double(NREGN + 1);
				for (int i = 1; i <= NREGN; ++i)
					XI(ILEFT + i) = XLEFT + double(i) * DX;
				ILEFT = IRIGHT;
			}
			XLEFT = XRIGHT;
		}
		break;
	}
n100:
	case 2: {
		// mode=2  halve the current mesh (i.e. double its size)
		int N2 = 2 * N;

		// check that n does not exceed storage limitations
		if (N2 > NMAX) {
			//  if possible, try with n=nmax. redistribute first.
			if (MODE != 2) {
				N = NMAX / 2;
				goto n220;
			}
			if (IPRINT < 1)  
				fmt::print("EXPECTED N TOO LARGE\n");
			N = N2;
			return;
		}
		//  calculate the old approximate solution values at
		//  points to be used in  errchk  for error estimates.
		//  if  mshflg  =1 an error estimate was obtained for
		//  for the old approximation so half the needed values
		//  will already be in  valstr .
		if (MSHFLG != 0) {
			//  save in  valstr  the values of the old solution
			//  at the relative positions 1/6 and 5/6 in each subinterval.
			int KSTORE = 1;
			for (int i = 1; i <= NOLD; ++i) {
				double HD6 = (XIOLD(i + 1) - XIOLD(i)) / 6.0;
				double X = XIOLD(i) + HD6;
				APPROX(i, X, VALSTR.sub(KSTORE).contiguous(), nullptr,
					ASAVE.sub(1, 1).contiguous(), nullptr, XIOLD.contiguous(),
					NOLD, Z.contiguous(), DMZ.contiguous(), K, NCOMP, NY, MMAX,
					MT.contiguous(), MSTAR, 4, nullptr, 0);
				X = X + 4.0 * HD6;
				KSTORE = KSTORE + 3 * MSTAR;
				APPROX(i, X, VALSTR.sub(KSTORE).contiguous(), nullptr,
					ASAVE.sub(1, 4).contiguous(), nullptr, XIOLD.contiguous(),
					NOLD, Z.contiguous(), DMZ.contiguous(), K, NCOMP, NY, MMAX, 
					MT.contiguous(), MSTAR, 4, nullptr, 0);
				KSTORE += MSTAR;
			}
		}
		//  save in  valstr  the values of the old solution
		//  at the relative positions 1/6, 2/6, 4/6 and 5/6 in
		//  each subinterval.
		else {
			int KSTORE = 1;
			for (int i = 1; i <= N; ++i) {
				double X = XI(i);
				double HD6 = (XI(i + 1) - XI(i)) / 6.0;
				for (int j = 1; j <= 4; ++j) {
					X = X + HD6;
					if (j == 3)
						X = X + HD6;
					APPROX(i, X, VALSTR.sub(KSTORE).contiguous(), nullptr, 
						ASAVE.sub(1, j).contiguous(), nullptr,
						XIOLD.contiguous(), NOLD, Z.contiguous(), DMZ.contiguous(),
						K, NCOMP, NY, MMAX, MT.contiguous(), MSTAR, 4, nullptr, 0);
					KSTORE = KSTORE + MSTAR;
				}
			}
		}
		MSHFLG = 0;
		MSHNUM = 1;
		MODE = 2;

		//  generate the halved mesh.
		for (int i = 1, j = 2; i <= N; ++i) {
			XI(j) = (XIOLD(i) + XIOLD(i + 1)) / 2.0;
			XI(j + 1) = XIOLD(i + 1);
			j = j + 2;
		}
		N = N2;
		break;
	}
	case 1: {
		{
			// mode=1  we do mesh selection if it is deemed worthwhile

			if (NOLD == 1)
				goto n100;
			if (NOLD <= 2 * NFXPNT)
				goto n100;

			//  we now project DMZ for mesh selection strategy, if required
			//  but set DMV = DMZ in case it is not
			int IDMZ = 1;
			for (int i = 1; i <= NOLD; ++i) {
				for (int KK = 1; KK <= K; ++KK) {
					for (int j = 1; j <= NCY; ++j) {
						DMV(IDMZ) = DMZ(IDMZ);
						IDMZ = IDMZ + 1;
					}
				}
			}
			if (INDEX != 1 && NY > 0) {
				IDMZ = 1;
				for (int i = 1; i <= NOLD; ++i) {
					double XI1 = XIOLD(i + 1);
					APPROX(i, XI1, ZVAL.contiguous(), YVAL.contiguous(), A.contiguous(),
						COEF.contiguous(), XIOLD.contiguous(), NOLD, Z.contiguous(), DMZ.contiguous(),
						K, NCOMP, NY, MMAX, MT.contiguous(), MSTAR, 3, nullptr, 1);
					dfsub(XI1, ZVAL.contiguous(), YVAL.contiguous(), DF.contiguous());

					// if index=2, form projection matrices directly
					// otherwise use svd to define appropriate projection
					if (INDEX == 0) {
						PRJSVD(FC.contiguous(), DF.contiguous(), CB.contiguous(),
							U.contiguous(), V.contiguous(), IPVTCB.contiguous(), ISING, 2);
					}
					else {
						// form cb
						for (int j = 1; j <= NY; ++j) {
							for (int J1 = 1; J1 <= NY; ++J1) {
								double FACT = 0.0;
								for (int l = 1, ML = 0; l <= NCOMP; ++l) {
									ML = ML + MT(l);
									FACT = FACT + DF(j + NCOMP, ML) * DF(l, MSTAR + J1);
								}
								CB(j, J1) = FACT;
							}
						}

						// decompose cb
						DGEFA(CB, NY, NY, IPVTCB, ISING);
						if (ISING != 0)
							return;

						// form columns of fc
						int ML = 0;
						for (int l = 1; l <= NCOMP; ++l) {
							ML += MT(l);
							for (int J1 = 1; J1 <= NY; ++J1)
								BCOL(J1) = DF(J1 + NCOMP, ML);

							DGESL(CB, NY, NY, IPVTCB, BCOL, 0);

							for (int J1 = 1; J1 <= NCOMP; ++J1) {
								double FACT = 0.0;
								for (int j = 1; j <= NY; ++j)
									FACT += DF(J1, j + MSTAR) * BCOL(j);
								FC(J1, l) = FACT;
							}
						}
					}

					// finally, replace fc with the true projection SR = i - fc
					for (int j = 1; j <= NCOMP; ++j) {
						for (int l = 1; l <= NCOMP; ++l) {
							FC(j, l) = -FC(j, l);
							if (j == l)
								FC(j, l) = FC(j, l) + 1.0;
						}
					}

					// project DMZ for the k collocation points, store in DMV
					for (int KK = 1; KK <= K; ++KK) {
						for (int j = 1; j <= NCOMP; ++j) {
							double FACT = 0.0;
							for (int l = 1; l <= NCOMP; ++l)
								FACT = FACT + FC(j, l) * DMZ(IDMZ + l - 1);
							DMV(IDMZ + j - 1) = FACT;
						}
						IDMZ = IDMZ + NCY;
					}
				}
			}

			//  the first interval has to be treated separately from the
			//  other intervals (generally the solution on the (i-1)st and ith
			//  intervals will be used to approximate the needed derivative, but
			//  here the 1st and second intervals are used.)
			double HIOLD = XIOLD(2) - XIOLD(1);
			HORDER(1, D1.contiguous(), HIOLD, DMV.contiguous());
			IDMZ = IDMZ + (NCOMP + NY) * K;
			HIOLD = XIOLD(3) - XIOLD(2);
			HORDER(2, D2.contiguous(), HIOLD, DMV.contiguous());
			ACCUM(1) = 0.0;
			SLOPE(1) = 0.0;
			double ONEOVH = 2.0 / (XIOLD(3) - XIOLD(1));
			for (int j = 1; j <= NTOL; ++j) {
				int JJ = JTOL(j);
				int JZ = LTOL(j);
				SLOPE(1) = std::max(SLOPE(1),
					pow(abs(D2(JJ) - D1(JJ)) * WGTMSH(j) * ONEOVH / (1.0 + abs(Z(JZ))),
						ROOT(j)));
			}
			double SLPHMX = SLOPE(1) * (XIOLD(2) - XIOLD(1));
			ACCUM(2) = SLPHMX;
			int IFLIP = 1;
			
			//  go through the remaining intervals generating  slope
			//  and  accum .
			for (int i = 2; i <= NOLD; ++i) {
				HIOLD = XIOLD(i + 1) - XIOLD(i);
				if (IFLIP == -1)
					HORDER(i, D1.contiguous(), HIOLD, DMV.contiguous());
				if (IFLIP == 1)
					HORDER(i, D2.contiguous(), HIOLD, DMV.contiguous());
				ONEOVH = 2.0 / (XIOLD(i + 1) - XIOLD(i - 1));
				SLOPE(i) = 0.0;


				// evaluate function to be equidistributed
				for (int j = 1; j <= NTOL; ++j) {
					int JJ = JTOL(j);
					int JZ = LTOL(j) + (i - 1) * MSTAR;
					auto temp = abs(D2(JJ) - D1(JJ)) * WGTMSH(j) * ONEOVH / (1.0 + abs(Z(JZ)));
					SLOPE(i) = std::max(SLOPE(i), pow(temp, ROOT(j)));
					/*	fmt::print(fg(fmt::color::orange), "abs(D2(JJ) - D1(JJ))  = {}\n", abs(D2(JJ) - D1(JJ)));
						fmt::print(fg(fmt::color::orange), "WGTMSH(j)  = {}\n", WGTMSH(j));
						fmt::print(fg(fmt::color::orange), "ONEOVH  = {}\n", ONEOVH);
						fmt::print(fg(fmt::color::orange), "1/ (1.0 + abs(Z(JZ)))  = {}\n", 1. / (1.0 + abs(Z(JZ))));
						fmt::print(fg(fmt::color::orange), "temp  = {}\n", temp);
						fmt::print(fg(fmt::color::orange), "ROOT(j)  = {}\n", ROOT(j));
						fmt::print(fg(fmt::color::orange_red), "SLOPE(i)  = {}\n", SLOPE(i));*/
				}

				// accumulate approximate integral of function to be equidistributed
				double TEMP = SLOPE(i) * (XIOLD(i + 1) - XIOLD(i));
				SLPHMX = std::max(SLPHMX, TEMP);
				ACCUM(i + 1) = ACCUM(i) + TEMP;
				IFLIP = -IFLIP;
				/*fmt::print(fg(fmt::color::green), "SLOPE({}) = {}\n", i, SLOPE(i)); 
				fmt::print(fg(fmt::color::green), "ACCUM({} + 1) = {}\n", i, ACCUM(i + 1));*/

				/*fmt::print(fg(fmt::color::green), "SLOPE = {}\n", SLOPE(i));
				fmt::print(fg(fmt::color::green), "diff  = {}\n", (XIOLD(i + 1) - XIOLD(i)));
				fmt::print(fg(fmt::color::green), "TEMP  = {}\n", TEMP);*/
			}

			
			double AVRG = ACCUM(NOLD + 1) / double(NOLD);
			double DEGEQU = AVRG / std::max(SLPHMX, PRECIS);

			//  naccum=expected n to achieve .1x user requested tolerances
			int NACCUM = int(ACCUM(NOLD + 1) + 1.0);
			if (IPRINT < 0)  
				fmt::print("MESH SELECTION INFO: DEGREE OF EQUIDISTRIBUTION = {}, "
							"PREDICTION FOR REQUIRED N = {}\n", DEGEQU, NACCUM);

			//  decide if mesh selection is worthwhile (otherwise, halve)
			if (AVRG < PRECIS)
				goto n100;
			if (DEGEQU >= .5)
				goto n100;

			//  nmx assures mesh has at least half as many subintervals as the
			//  previous mesh
			int NMX = std::max(NOLD + 1, NACCUM) / 2;

			//  this assures that halving will be possible later (for error est)
			int NMAX2 = NMAX / 2;

			//  the mesh is at most halved
			N = std::min({ NMAX2, NOLD, NMX });
		}
	n220:
		{
			int NOLDP1 = NOLD + 1;
			if (N < NFXP1)
				N = NFXP1;
			MSHNUM = MSHNUM + 1;

			//  if the new mesh is smaller than the old mesh set mshnum
			//  so that the next call to  newmsh  will produce a halved
			//  mesh. if n .eq. nold / 2 increment mshalt so there can not
			//  be an infinite loop alternating between n and n/2 points.
			if (N < NOLD)
				MSHNUM = MSHLMT;
			if (N > NOLD / 2)
				MSHALT = 1;
			if (N == NOLD / 2)
				MSHALT = MSHALT + 1;
			MSHFLG = 0;

			//  having decided to generate a new mesh with n subintervals we now
			//  do so, taking into account that the nfxpnt points in the array
			//  fixpnt must be included in the new mesh.
			int IN = 1;
			double ACCL = 0.0;
			int LOLD = 2;
			XI(1) = TLEFT;
			XI(N + 1) = TRIGHT;
			for (int i = 1; i <= NFXP1; ++i) {
				double ACCR; int LNEW, NREGN;
				if (i != NFXP1) {
					for (int j = LOLD; j <= NOLDP1; ++j) {
						LNEW = j;
						if (FIXPNT(i) <= XIOLD(j))
							break;
					}
					ACCR = ACCUM(LNEW) + (FIXPNT(i) - XIOLD(LNEW)) * SLOPE(LNEW - 1);
					NREGN = int((ACCR - ACCL) / ACCUM(NOLDP1) * double(N) - .5);
					NREGN = std::min(NREGN, N - IN - NFXP1 + i);
					XI(IN + NREGN + 1) = FIXPNT(i);
				}
				else {
					ACCR = ACCUM(NOLDP1);
					LNEW = NOLDP1;
					NREGN = N - IN;
				}
				if (NREGN != 0) {
					double TEMP = ACCL;
					double TSUM = (ACCR - ACCL) / double(NREGN + 1);
					for (int j = 1; j <= NREGN; ++j) {
						IN = IN + 1;
						TEMP = TEMP + TSUM;
						int LCARRY;
						for (int l = LOLD; l <= LNEW; ++l) {
							LCARRY = l;
							if (TEMP <= ACCUM(l))
								break;
						}
						LOLD = LCARRY;
						XI(IN) = XIOLD(LOLD - 1) + (TEMP - ACCUM(LOLD - 1)) / SLOPE(LOLD - 1);
					}
				}
				IN = IN + 1;;
				ACCL = ACCR;
				LOLD = LNEW;
			}
			MODE = 1;
			break;
		}
	}
	} // end of switch


	if (IPRINT < 1) {
		//assert(XI.getSize() == N + 1);
		fmt::print("THE NEW MESH (OF {} SUBINTERVALS): ", N);
		for (int i = 1;i<=N+1;++i)
			fmt::print("{:.2}, ", XI(i));
		fmt::print("\n");
	}
	NZ = MSTAR * (N + 1);
	NDMZ = KDY * N;
}


//**********************************************************************
//
//   purpose
//            assign (once) values to various array constants.
//
//   arrays assigned during compilation:
//     cnsts1 - weights for extrapolation error estimate
//     cnsts2 - weights for mesh selection
//              (the above weights come from the theoretical form for
//              the collocation error -- see [5])
//
//   arrays assigned during execution:
//     wgterr - the particular values of cnsts1 used for current run
//              (depending on k, m)
//     wgtmsh - gotten from the values of cnsts2 which in turn are
//              the constants in the theoretical expression for the
//              errors. the quantities in wgtmsh are 10x the values
//              in cnsts2 so that the mesh selection algorithm
//              is aiming for errors .1x as large as the user
//              requested tolerances.
//     jtol   - components of differential system to which tolerances
//              refer (viz, if ltol(i) refers to a derivative of u(j),
//              then jtol(i)=j)
//     root   - reciprocals of expected rates of convergence of compo-
//              nents of z(j) for which tolerances are specified
//     rho    - the k collocation points on (0,1)
//     coef   -
//     acol  -  the runge-kutta coefficients values at collocation
//              points
//
//**********************************************************************
void CONSTS()
{
	RHO.assertDim(7);
	COEF.reshape(K, K);
	COEF.assertDim(K, K);
	
	//using namespace COLORD; // int K, NC, NNY, NCY, MSTAR, KD, KDY, MMAX; iad1 MT(20);
	//using namespace COLEST; // dad1 TTL(40), WGTMSH(40), WGTERR(40), TOLIN(40), ROOT(40);  iad1 JTOL(40), LTTOL(40); int NTOL;
	//using namespace COLBAS; // dad2 B(7, 4), ACOL(28, 7), ASAVE(28, 4);

	dad1 CNSTS1(28), CNSTS2(28);

	CNSTS1 = { 0.25e0, 0.625e-1,  7.2169e-2, 1.8342e-2,
	      1.9065e-2, 5.8190e-2, 5.4658e-3, 5.3370e-3, 1.8890e-2,
	      2.7792e-2, 1.6095e-3, 1.4964e-3, 7.5938e-3, 5.7573e-3,
	      1.8342e-2, 4.673e-3,  4.150e-4,  1.919e-3,  1.468e-3,
	      6.371e-3,  4.610e-3,  1.342e-4,  1.138e-4,  4.889e-4,
	      4.177e-4,  1.374e-3,  1.654e-3,  2.863e-3 };
	CNSTS2 = { 1.25e-1,   2.604e-3,  8.019e-3,  2.170e-5,
	     7.453e-5,  5.208e-4,  9.689e-8,  3.689e-7,  3.100e-6,
	    2.451e-5,  2.691e-10, 1.120e-9,  1.076e-8,  9.405e-8,
	      1.033e-6,  5.097e-13, 2.290e-12, 2.446e-11, 2.331e-10,
	      2.936e-9,  3.593e-8,  7.001e-16, 3.363e-15, 3.921e-14,
	      4.028e-13, 5.646e-12, 7.531e-11, 1.129e-9 };

	// assign weights for error estimate
	int KOFF = K * (K + 1) / 2;
	int IZ = 1;
	for (int j = 1; j <= NCOMP; ++j) {
		int MJ = MT(j);
		for (int l = 1; l <= MJ; ++l) {
			WGTERR(IZ) = CNSTS1(KOFF - MJ + l);
			IZ = IZ + 1;
		}
	}

	// assign array values for mesh selection: wgtmsh, jtol, and root
	int JCOMP = 1;
	int MTOT = MT(1);
	for (int i = 1; i <= NTOL; ++i) {
		int LTOLI = LTOL(i);
		while (true) {
			if (LTOLI <= MTOT)
				break;
			JCOMP = JCOMP + 1;
			MTOT = MTOT + MT(JCOMP);
		}
		JTOL(i) = JCOMP;
		WGTMSH(i) = 1.e1 * CNSTS2(KOFF + LTOLI - MTOT) / TOLIN(i);
		ROOT(i) = 1.0 / double(K + MTOT - LTOLI + 1);
	}

	// specify collocation points
	switch (K) {
	case 1: RHO(1) = 0.0;
		break;

	case 2:
		RHO(2) = .577350269189625764510;
		RHO(1) = -RHO(2);
		break;

	case 3:
		RHO(3) = .774596669241483377040;
		RHO(2) = .00;
		RHO(1) = -RHO(3);
		break;

	case 4:
		RHO(4) = .861136311594052575230;
		RHO(3) = .339981043584856264800;
		RHO(2) = -RHO(3);
		RHO(1) = -RHO(4);
		break;

	case 5:
		RHO(5) = .906179845938663992800;
		RHO(4) = .538469310105683091040;
		RHO(3) = .00;
		RHO(2) = -RHO(4);
		RHO(1) = -RHO(5);
		break;

	case 6:
		RHO(6) = .932469514203152027810;
		RHO(5) = .661209386466264513660;
		RHO(4) = .238619186083196908630;
		RHO(3) = -RHO(4);
		RHO(2) = -RHO(5);
		RHO(1) = -RHO(6);
		break;

	case 7:
		RHO(7) = .9491079912342758524520;
		RHO(6) = .741531185599394439860;
		RHO(5) = .405845151377397166900;
		RHO(4) = 0.0;
		RHO(3) = -RHO(5);
		RHO(2) = -RHO(6);
		RHO(1) = -RHO(7);
		break;
	}

	// map (-1,1) to (0,1) by  t = .5 * (1. + x)
	for (int j = 1; j<=K; ++j)
		RHO(j) = .50 * (1.0 + RHO(j));


	// now find runge-kutta coeffitients b, acol and asave
	// the values of asave are to be used in  newmsh  and errchk .
	for (int j = 1; j <= K; ++j) {
		for (int i = 1; i <= K; ++i)
			COEF(i, j) = 0.0;
		COEF(j, j) = 1.0;
		VMONDE(COEF.sub(1, j).contiguous(), K);
	}
	RKBAS(1.0, K, MMAX, B.contiguous(), nullptr, 0);
	for (int i = 1; i <= K; ++i)
		RKBAS(RHO(i), K, MMAX, ACOL.sub(1, i).contiguous(), nullptr, 0);

	RKBAS(1.0 / 6.0, K, MMAX, ASAVE.sub(1, 1).contiguous(), nullptr, 0);
	RKBAS(1.0 / 3.0, K, MMAX, ASAVE.sub(1, 2).contiguous(), nullptr, 0);
	RKBAS(2.0 / 3.0, K, MMAX, ASAVE.sub(1, 3).contiguous(), nullptr, 0);
	RKBAS(5.0 / 6.0, K, MMAX, ASAVE.sub(1, 4).contiguous(), nullptr, 0);

}



//**********************************************************************
//
//      purpose
//               determine the error estimates and test to see if the
//               error tolerances are satisfied.
//
//      variables
//        xi     - current mesh points
//        valstr - values of the previous solution which are needed
//                 for the extrapolation- like error estimate.
//        wgterr - weights used in the extrapolation-like error
//                 estimate. the array values are assigned in
//                 subroutine  consts.
//        errest - storage for error estimates
//        err    - temporary storage used for error estimates
//        z      - approximate solution on mesh xi
//        ifin   - a 0-1 variable. on return it indicates whether
//                 the error tolerances were satisfied
//        mshflg - is set by errchk to indicate to newmsh whether
//                 any values of the current solution are stored in
//                 the array valstr. (0 for no, 1 for yes)
//
//**********************************************************************
void ERRCHK(dar1 XI, dar1 Z, dar1 DMZ, dar1 VALSTR, int& IFIN)
{
	XI.assertDim(1);
	Z.assertDim(1);
	DMZ.assertDim(1);
	VALSTR.assertDim(1);

	//using namespace COLOUT; // double PRECIS; int IOUT, IPRINT; 
	//using namespace COLORD; // int K, NC, NNY, NCY, MSTAR, KD, KDY, MMAX; iad1 MT(20);
	//using namespace COLAPR; // int N, NOLD, NMAX, NZ, NDMZ;
	//using namespace COLMSH; // int MSHFLG, MSHNUM, MSHLMT, MSHALT;
	//using namespace COLEST; // dad1 TTL(40), WGTMSH(40), WGTERR(40), TOLIN(40), ROOT(40);  iad1 JTOL(40), LTTOL(40); int NTOL;
	//using namespace COLBAS; // dad2 B(7, 4), ACOL(28, 7), ASAVE(28, 4);

	MSHFLG = 1;

	dad1 ERR(40), ERREST(40), DUMMY(1);


	//  error estimates are to be generated and tested
	//  to see if the tolerance requirements are satisfied.
	IFIN = 1;
	for (int j = 1; j <= MSTAR; ++j)
		ERREST(j) = 0.0;
	for (int IBACK = 1; IBACK <= N; ++IBACK) {
		int i = N + 1 - IBACK;

		// the error estimates are obtained by combining values of
		// the numerical solutions for two meshes.
		// for each value of iback we will consider the two approximations at 2 points in each of
		// the new subintervals. 
		// we work backwards through the subinterval so that new values can be stored
		// in valstr in case they prove to be needed later for an error estimate.
		// The routine  newmsh filled in the needed values of the old solution in valstr.
		int KNEW = (4 * (i - 1) + 2) * MSTAR + 1;
		int KSTORE = (2 * (i - 1) + 1) * MSTAR + 1;
		double X = XI(i) + (XI(i + 1) - XI(i)) * 2.0 / 3.0;
		APPROX(i, X, VALSTR.sub(KNEW).contiguous(), nullptr, ASAVE.sub(1, 3).contiguous(),
			nullptr, XI.contiguous(), N, Z.contiguous(), DMZ.contiguous(), K, NCOMP, NY, 
			MMAX, MT.contiguous(), MSTAR, 4, nullptr, 0);
		for (int l = 1; l <= MSTAR; ++l) {
			ERR(l) = WGTERR(l) * abs(VALSTR(KNEW) - VALSTR(KSTORE));
			KNEW = KNEW + 1;
			KSTORE = KSTORE + 1;
		}
		KNEW = (4 * (i - 1) + 1) * MSTAR + 1;
		KSTORE = 2 * (i - 1) * MSTAR + 1;
		X = XI(i) + (XI(i + 1) - XI(i)) / 3.0;
		APPROX(i, X, VALSTR.sub(KNEW).contiguous(), nullptr, ASAVE.sub(1, 2).contiguous(), 
			nullptr, XI.contiguous(), N, Z.contiguous(), DMZ.contiguous(), K, NCOMP, NY,
			MMAX, MT.contiguous(), MSTAR, 4, nullptr, 0);
		for (int l = 1; l <= MSTAR; ++l) {
			ERR(l) = ERR(l) + WGTERR(l) * abs(VALSTR(KNEW) - VALSTR(KSTORE));
			KNEW = KNEW + 1;
			KSTORE = KSTORE + 1;
		}

		// find component-wise maximum error estimate
		for (int l = 1; l <= MSTAR; ++l)
			ERREST(l) = std::max(ERREST(l), ERR(l));

		// test whether the tolerance requirements are satisfied
		// in the i-th interval.
		if (IFIN == 0)
			continue;
		for (int j = 1; j <= NTOL; ++j) {
			int LTOLJ = LTOL(j);
			int LTJZ = LTOLJ + (i - 1) * MSTAR;
			if (ERR(LTOLJ) > TOLIN(j) * (abs(Z(LTJZ)) + 1.0))
				IFIN = 0;
		}
	}

	if (IPRINT >= 0)
		return;
	fmt::print("THE ESTIMATED ERRORS ARE\n");
	int LJ = 1;
	for (int j = 1; j <= NCOMP; ++j) {
		int MJ = LJ - 1 + MT(j);
		fmt::print("{}:  ", j);
		for (int l = LJ; l <= MJ; ++l)
			fmt::print("{}, ", ERREST(l));
		LJ = MJ + 1;
	}
	fmt::print("\n");
}




//---------------------------------------------------------------------
//                            p a r t  3
//          collocation system setup routines
//---------------------------------------------------------------------


//*********************************************************************
//
//   purpose
//         this routine controls the set up and solution of a linear
//      system of collocation equations.
//         the matrix  g  is cast into an almost block diagonal
//      form by an appropriate ordering of the columns and solved
//      using the package of de boor-weiss [7] modified.
//      the matrix is composed of n blocks. the i-th block has the size
//                  integs(1,i) * integs(2,i).
//      it contains in its last rows the linearized collocation
//      equations, condensed as described in [4],
//      and the linearized side conditions corresponding to
//      the i-th subinterval.  integs(3,i)  steps of gaussian
//      elimination are applied to it to achieve a  partial plu
//      decomposition.  the right hand side vector is put into  rhs
//      and the solution vector is returned in  delz and deldmz.
//      note that the presence of algebraic solution components
//      does not affect the structure (size) of g -- only the contents
//      of the blocks (and the size of deldmz) changes.
//
//         lsyslv operates according to one of 5 modes:
//      mode = 0 - set up the collocation matrices  v , w , g
//                 and the right hand side  rhs ,  and solve.
//                 (for linear problems only.)
//      mode = 1 - set up the collocation matrices  v , w , g
//                 and the right hand sides  rhs  and  dmzo ,
//                 and solve. also set up  integs .
//                 (first iteration of nonlinear problems only).
//      mode = 2 - set up  rhs  only and compute its norm.
//      mode = 3 - set up  v, w, g  only and solve system.
//      mode = 4 - perform forward and backward substitution only
//                 (do not set up the matrices nor form the rhs).
//
//   variables
//
//      ig,izeta  - pointers to g,zeta respectively
//                       (necessary to keep track of blocks of g
//                       during matrix manipulations)
//      idmz,irhs,iv,iw - pointers to  rhs,v,w rspectively
//      df    - partial derivatives of f from dfsub
//      rnorm - euclidean norm of rhs
//      lside - number of side conditions in current and previous blocks
//      iguess = 1 when current soln is user specified via  guess
//             = 0 otherwise
//      dmzo  - an array used to project the initial solution into
//              the current pp-space, when mode=1.
//              in the case mode=1 the current solution iterate may not
//              be in the right space, being defined by an arbitrary
//              users guess or as a pp on a different mesh.
//              when forming collocation equations we are using values
//              of z, y and dmval at collocation points and of z at
//              boundary points. at the end of lsyslv (with mode=1)
//              a similar projection used to obtain the corrections
//              delz and deldmz is used to obtain the projected initial
//              iterate z and dmz.
//      fc    - an array used to store projection matrices for
//              the case of projected collocation
//
//
//*********************************************************************

double* observer = nullptr;

void LSYSLV(int& MSING, dar1 XI, dar1 XIOLD, dar1 Z, dar1 DMZ, dar1 DELZ, dar1 DELDMZ,
	dar1 G, dar1 W, dar1 V, dar1 FC, dar1 RHS, dar1 DMZO,
	iar2 INTEGS, iar1 IPVTG, iar1 IPVTW, double& RNORM,
	const int MODE, fsub_t fsub, dfsub_t dfsub, gsub_t gsub,
	dgsub_t dgsub, guess_t guess, int& ISING)
{
	Z.assertDim(1);
	DMZ.assertDim(1);
	DELZ.assertDim(1);
	DELDMZ.assertDim(1);
	XI.assertDim(1);
	XIOLD.assertDim(1);
	G.assertDim(1);
	W.assertDim(1);
	V.assertDim(1);
	RHS.assertDim(1);
	DMZO.assertDim(1);
	INTEGS.reshape(3, 1);
	INTEGS.assertDim(3, 1);
	IPVTG.assertDim(1);
	IPVTW.assertDim(1);
	FC.assertDim(1);

	//using namespace COLOUT; // double PRECIS; int IOUT, IPRINT; 
	//using namespace COLLOC; // dad1 RHO(7), COEF(49);
	//using namespace COLORD; // int K, NC, NNY, NCY, MSTAR, KD, KDY, MMAX; iad1 MT(20);
	//using namespace COLAPR; // int N, NOLD, NMAX, NZ, NDMZ;
	//using namespace COLSID; // dad1 TZETA(40);  double TLEFT, TRIGHT;  int IZETA, IDUM;
	//using namespace COLNLN; // int NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX;
	//using namespace COLBAS; // dad2 B(7, 4), ACOL(28, 7), ASAVE(28, 4);

	//fmt::print(fg(fmt::color::orange_red),"Enter LSYSLV and RHS(41) = {}\n", RHS(41));
	//fmt::print(fg(fmt::color::orange_red), "Address is {}\n", 
		//(void*)(RHS.contiguous() + 40));
	//observer = RHS.contiguous() + 40;


	//fmt::print(fg(fmt::color::medium_purple), "observer = {}\n", *observer);
	
	dad1 YVAL(20), ZVAL(40), F(40), DGZ(40), DMVAL(20), DF(800);
	dad1 DUMMY(1), Y(1), AT(28), CB(400);
	iad1 IPVTCB(20);

	int NYCB;
	if (NY == 0)
		NYCB = 1;
	else
		NYCB = NY;

	int INFC = (MSTAR + NY) * NCOMP;
	int M1 = MODE + 1;
	double XI1 = 0;


	switch (M1) {
	case 1: {
		//  linear problem initialization
		for (int i = 1; i <= MSTAR; ++i)
			ZVAL(i) = 0.0;
		for (int i = 1; i <= NY; ++i)
			YVAL(i) = 0.0;

		[[fallthrough]];
	}
	case 2:
	case 3:
	case 4: {
		//  initialization
		int IDMZ = 1;
		int IDMZO = 1;
		int IRHS = 1;
		int IG = 1;
		int IW = 1;
		int IV = 1;
		int IFC = 1;
		IZETA = 1;
		int LSIDE = 0;
		int IOLD = 1;
		int NCOL = 2 * MSTAR;
		RNORM = 0.0;
		if (MODE <= 1) {
			//  build integs (describing block structure of matrix)
			for (int i = 1; i <= N; ++i) {
				INTEGS(2, i) = NCOL;
				if (i >= N) {
					INTEGS(3, N) = NCOL;
					LSIDE = MSTAR;
				}
				else {
					INTEGS(3, i) = MSTAR;
					while (true) {
						if (LSIDE == MSTAR)
							break;
						if (TZETA(LSIDE + 1) >= XI(i) + PRECIS)
							break;
						LSIDE = LSIDE + 1;
					}
				}
				int NROW = MSTAR + LSIDE;
				INTEGS(1, i) = NROW;
			}
		}
		if (MODE != 2) {
			//  zero the matrices to be computed
			//fmt::print(fg(fmt::color::yellow_green), "Zero the W matrix\n");
			int LW = KDY * KDY * N;
			for (int l = 1; l <= LW; ++l)
				W(l) = 0.0;
		}
		// set up the linear system of equations
		for (int i = 1; i <= N; ++i) {
			// construct a block of  a  and a corresponding piece of  rhs.
			double XII = XI(i);
			double H = XI(i + 1) - XI(i);
			int NROW = INTEGS(1, i);

			// go thru the ncomp collocation equations and side conditions
			// in the i-th subinterval
			while (true) {
				if (IZETA > MSTAR)
					break;
				if (TZETA(IZETA) > XII + PRECIS)
					break;

				// build equation for a side condition.
				if (MODE != 0) {
					if (IGUESS == 1) {
						// case where user provided current approximation
						guess(XII, ZVAL.contiguous(), YVAL.contiguous(), DMVAL.contiguous());
					}
					else {
						// other nonlinear case
						if (MODE == 1) {
							APPROX(IOLD, XII, ZVAL.contiguous(), Y.contiguous(), AT.contiguous(),
								COEF.contiguous(), XIOLD.contiguous(), NOLD, Z.contiguous(),
								DMZ.contiguous(),K, NCOMP, NY, MMAX, MT.contiguous(),
								MSTAR, 2, nullptr, 0);
						}
						else {
							APPROX(i, XII, ZVAL.contiguous(), Y.contiguous(), AT.contiguous(),
								nullptr, XI.contiguous(), N, Z.contiguous(), DMZ.contiguous(),
								K, NCOMP, NY, MMAX, MT.contiguous(), MSTAR, 1, nullptr, 0);
							if (MODE == 3)
								goto n120;
						}
					}
				}
				// find  rhs  boundary value.
				double GVAL;
				gsub(IZETA, ZVAL.contiguous(), GVAL);
				RHS(NDMZ + IZETA) = -GVAL;
				RNORM = RNORM + GVAL * GVAL;
				if (MODE != 2) {
				n120:
					// build a row of  a  corresponding to a boundary point
					GDERIV(G.sub(IG).contiguous(), NROW, IZETA, ZVAL.contiguous(),
						DGZ.contiguous(), 1, dgsub);
				}
				IZETA = IZETA + 1;
			}

			// assemble collocation equations
			for (int j = 1; j <= K; ++j) {
				double HRHO = H * RHO(j);
				double XCOL = XII + HRHO;

				// this value corresponds to a collocation (interior)
				// point. build the corresponding  ncy  equations.
				if (MODE == 0)
					goto n200;
				if (IGUESS != 1)
					goto n160;

				// use initial approximation provided by the user.
				guess(XCOL, ZVAL.contiguous(), YVAL.contiguous(), DMZO.sub(IRHS).contiguous());
				goto n170;


			n160:
				if (MODE == 1) {
					// find  rhs  values
					APPROX(IOLD, XCOL, ZVAL.contiguous(), YVAL.contiguous(), AT.contiguous(),
						COEF.contiguous(), XIOLD.contiguous(),
						NOLD, Z.contiguous(), DMZ.contiguous(), K, NCOMP, NY, MMAX,
						MT.contiguous(), MSTAR, 2,
						DMZO.sub(IRHS).contiguous(), 2);

				n170:
					fsub(XCOL, ZVAL.contiguous(), YVAL.contiguous(), F.contiguous());
					for (int JJ = NCOMP + 1; JJ <= NCY; ++JJ)
						DMZO(IRHS + JJ - 1) = 0.0;

					for (int JJ = 1; JJ <= NCY; ++JJ) {
						double VALUE = DMZO(IRHS) - F(JJ);
						RHS(IRHS) = -VALUE;
						RNORM = RNORM + VALUE * VALUE;
						IRHS = IRHS + 1;
					}
				}
				else {
					// evaluate former collocation solution
					APPROX(i, XCOL, ZVAL.contiguous(), Y.contiguous(), ACOL.sub(1, j).contiguous(),
						COEF.contiguous(), XI.contiguous(), N, Z.contiguous(),
						DMZ.contiguous(), K, NCOMP, NY, MMAX, MT.contiguous(), MSTAR, 
						4, nullptr, 0);
					if (MODE == 3)
						goto n210;

					// fill in  rhs  values (and accumulate its norm).
					fsub(XCOL, ZVAL.contiguous(), DMZ.sub(IRHS + NCOMP).contiguous(), F.contiguous());
					for (int JJ = 1; JJ <= NCY; ++JJ) {
						double VALUE = F(JJ);
						if (JJ <= NCOMP)
							VALUE = VALUE - DMZ(IRHS);
						RHS(IRHS) = VALUE;
						RNORM = RNORM + VALUE * VALUE;
						IRHS = IRHS + 1;
					}
					continue;

					// the linear case
				n200:
					fsub(XCOL, ZVAL.contiguous(), YVAL.contiguous(), RHS.sub(IRHS).contiguous());
					IRHS = IRHS + NCY;
				}
			n210:
				auto test1 = W(16);

				//fmt::print(fg(fmt::color::light_blue), "vor VWBLOK W(9) = {}\n", W(9));

				/*static int counter = 0;
				counter ++;
				if (counter == 4)
					logall = true;

				for (int ii = 1; ii <= 28; ++ii) {
					for (int jj = 1; jj <= 7; ++jj)
						fmt::print(fg(jj>=j?fmt::color::orange:fmt::color::white),
							"{:6.4f} ", ACOL(ii, jj));
					std::cout << std::endl;
				}
			*/

				// fill in ncy rows of  w and v
				VWBLOK(XCOL, HRHO, j, W.sub(IW).contiguous(), V.sub(IV).contiguous(),
					IPVTW.sub(IDMZ).contiguous(),
					ZVAL.contiguous(), YVAL.contiguous(), DF.contiguous(), 
					ACOL.sub(1, j).contiguous(), DMZO.sub(IDMZO).contiguous(), dfsub, MSING);

				//logall = false;

				//auto test2 = W(16);
				//fmt::print(fg(fmt::color::light_blue), "nach VWBLOK W(9) = {}\n", W(9));

				if (MSING != 0)
					return;
			}

			// build global bvp matrix  g
			if (INDEX != 1 && NY > 0) {
				// projected collocation: find solution at xi(i+1)
				XI1 = XI(i + 1);
				if (MODE != 0) {
					if (IGUESS == 1) {
						guess(XI1, ZVAL.contiguous(), YVAL.contiguous(), DMVAL.contiguous());
					}
					else {
						if (MODE == 1) {
							APPROX(IOLD, XI1, ZVAL.contiguous(), YVAL.contiguous(), AT.contiguous(),
								COEF.contiguous(), // here IOLD gets changed! simon
								XIOLD.contiguous(), NOLD, Z.contiguous(), DMZ.contiguous(),
								K, NCOMP, NY, MMAX,
								MT.contiguous(), MSTAR, 2, nullptr, 1);
							if (i == N) {
								auto temp = NOLD + 1;
								APPROX(temp, XI1, ZVAL.contiguous(), YVAL.contiguous(), 
									AT.contiguous(), COEF.contiguous(),
									XIOLD.contiguous(), NOLD, Z.contiguous(), 
									DMZ.contiguous(), K, NCOMP, NY, MMAX,
									MT.contiguous(), MSTAR, 1, nullptr, 0);
							}
						}
						else {
							APPROX(i, XI1, ZVAL.contiguous(), YVAL.contiguous(), AT.contiguous(),
								COEF.contiguous(),
								XI.contiguous(), N, Z.contiguous(), DMZ.contiguous(), K, NCOMP, NY, MMAX,
								MT.contiguous(), MSTAR, 3, nullptr, 1);
							auto temp = i + 1;
							APPROX(temp, XI1, ZVAL.contiguous(), YVAL.contiguous(), AT.contiguous(),
								COEF.contiguous(),XI.contiguous(), N, Z.contiguous(),
								DMZ.contiguous(), K, NCOMP, NY, MMAX,
								MT.contiguous(), MSTAR, 1, nullptr, 0);
						}
					}
				}

				// find rhs at next mesh point (also for linear case)
				fsub(XI1, ZVAL.contiguous(), YVAL.contiguous(), F.contiguous());
			}

			GBLOCK(H, G.sub(IG), NROW, IZETA, W.sub(IW).contiguous(), V.sub(IV),
				nullptr, DELDMZ.sub(IDMZ).contiguous(),
				IPVTW.sub(IDMZ), 1, MODE, XI1, ZVAL.contiguous(), YVAL.contiguous(),
				F.contiguous(), DF,
				CB, IPVTCB, FC.sub(IFC), dfsub, ISING, NYCB);
			
			if (ISING != 0)
				return;
			if (i >= N) {
				IZSAVE = IZETA;
				while (true) {
					if (IZETA > MSTAR)
						break;

					// build equation for a side condition.
					if (MODE != 0) {
						if (IGUESS == 1) {
							// case where user provided current approximation
							guess(ARIGHT, ZVAL.contiguous(), YVAL.contiguous(), DMVAL.contiguous());
						}
						else {
							// other nonlinear case
							if (MODE == 1) {
								auto temp = NOLD + 1;
								APPROX(temp, ARIGHT, ZVAL.contiguous(), Y.contiguous(), AT.contiguous(),
									COEF.contiguous(), XIOLD.contiguous(), NOLD, Z.contiguous(), 
									DMZ.contiguous(), K,
									NCOMP, NY, MMAX, MT.contiguous(), MSTAR, 1, nullptr, 0);
							}
							else {
								auto temp = N + 1;
								APPROX(temp, ARIGHT, ZVAL.contiguous(), Y.contiguous(), AT.contiguous(),
									COEF.contiguous(), XI.contiguous(), N, Z.contiguous(), 
									DMZ.contiguous(), K,
									NCOMP, NY, MMAX, MT.contiguous(), MSTAR, 1, nullptr, 0);
								if (MODE == 3)
									goto n260;
							}
						}
					}

					// find  rhs  boundary value.
					double GVAL;
					gsub(IZETA, ZVAL.contiguous(), GVAL);
					RHS(NDMZ + IZETA) = -GVAL;
					RNORM = RNORM + GVAL * GVAL;
					if (MODE != 2) {
					n260:
						// build a row of  a  corresponding to a boundary point
						GDERIV(G.sub(IG).contiguous(), NROW, IZETA + MSTAR, ZVAL.contiguous(),
							DGZ.contiguous(), 2, dgsub);
					}
					IZETA = IZETA + 1;
				}
			}
			else {
				// update counters -- i-th block completed
				IG = IG + NROW * NCOL;
				IV = IV + KDY * MSTAR;
				IW = IW + KDY * KDY;
				IDMZ = IDMZ + KDY;
				if (MODE == 1)
					IDMZO = IDMZO + KDY;
				IFC = IFC + INFC + 2 * NCOMP;
			}
		}

		// assembly process completed
		if (MODE != 0 && MODE != 3) {
			RNORM = sqrt(RNORM / double(NZ + NDMZ));
			if (MODE == 2) {
				//fmt::print(fg(fmt::color::orange_red), "Leave LSYSLV and RHS(41) = {}\n", RHS(41));
				return;
			}
		}
		//  solve the linear system.
		//  matrix decomposition
		FCBLOK(G.contiguous(), INTEGS.contiguous(), N, IPVTG.contiguous(), 
			DF.contiguous(), MSING);

		//  check for singular matrix
		MSING = -MSING;
		if (MSING != 0)
			return;
	} [[fallthrough]];
	case 5: {
		//  perform forward and backward substitution .
		for (int l = 1; l <= NDMZ; ++l)
			DELDMZ(l) = RHS(l);

		int IZ = 1;
		int IDMZ = 1;
		int IW = 1;
		int IFC = 1;
		int IZET = 1;
		for (int i = 1; i <= N; ++i) {
			int NROW = INTEGS(1, i);
			IZETA = NROW + 1 - MSTAR; 
			if (i == N)
				IZETA = IZSAVE;
			while (true) {
				if (IZET == IZETA)
					break;
				/*fmt::print("---------------------------------------\n");
				fmt::print("Set DELZ({}) = {}\n", IZ - 1 + IZET, RHS(NDMZ + IZET));*/
				DELZ(IZ - 1 + IZET) = RHS(NDMZ + IZET);
				IZET = IZET + 1;
			}
			double H = XI(i + 1) - XI(i);
//			logall = true;
			
			/*std::ofstream file("dgesl_prev.txt");
			for (int ii = 0; ii < KDY * KDY; ++ii)
				file << W.sub(IW).contiguous()[ii] << std::endl;
			file.close();

			fmt::print(fg(fmt::color::gold), "vorher DELZ = \n");
			for(int ii = 1; ii<=4;++ii)
				fmt::print(fg(fmt::color::gold), "{}  ", DELZ(ii));
			std::cout << std::endl;*/

			GBLOCK(H, G.sub(1), NROW, IZETA, W.sub(IW).contiguous(), V.sub(1), DELZ.sub(IZ).contiguous(),
				DELDMZ.sub(IDMZ).contiguous(),
				IPVTW.sub(IDMZ), 2, MODE, XI1, ZVAL.contiguous(), YVAL.contiguous(),
				FC.sub(IFC + INFC).contiguous(),
				DF, CB, IPVTCB, FC.sub(IFC), dfsub, ISING, NYCB);
			
		/*	fmt::print(fg(fmt::color::golden_rod), "nachher DELZ = \n");
			for (int ii = 1; ii <= 4; ++ii)
				fmt::print(fg(fmt::color::golden_rod), "{}  ", DELZ(ii));
			std::cout << std::endl;*/

			//logall = false;
			IZ = IZ + MSTAR;
			IDMZ = IDMZ + KDY;
			IW = IW + KDY * KDY;
			IFC = IFC + INFC + 2 * NCOMP;
			if (i < N)
				continue;
			while (true) {
				if (IZET > MSTAR)
					break;
				//fmt::print("Set DELZ({}) = {}\n", IZ - 1 + IZET, RHS(NDMZ + IZET));
				DELZ(IZ - 1 + IZET) = RHS(NDMZ + IZET);
				IZET = IZET + 1;
			}
		}

		for (int i = 1; i <= NZ; ++i) {
			//fmt::print("# DQZ(i)={:10.6f}\n", DELZ(i));
		}

		//  perform forward and backward substitution for mode=0,2, or 3.
		SBBLOK(G.contiguous(), INTEGS.contiguous(), N, IPVTG.contiguous(), DELZ.contiguous());

		//  finally find deldmz
		DMZSOL(V.contiguous(), DELZ.contiguous(), DELDMZ.contiguous());

		if (MODE != 1) {
			//fmt::print(fg(fmt::color::orange_red), "Leave LSYSLV and RHS(41) = {}\n", RHS(41));
			return;
		}

		//  project current iterate into current pp-space
		for (int l = 1; l <= NDMZ; ++l)
			DMZ(l) = DMZO(l);
		IZ = 1;
		IDMZ = 1;
		IW = 1;
		IFC = 1;
		IZET = 1;
		for (int i = 1; i <= N; ++i) {
			int NROW = INTEGS(1, i);
			IZETA = NROW + 1 - MSTAR;
			if (i == N)
				IZETA = IZSAVE;
			while (true) {
				if (IZET == IZETA)
					break;
				Z(IZ - 1 + IZET) = DGZ(IZET);
				IZET = IZET + 1;
			}
			double H = XI(i + 1) - XI(i);
			GBLOCK(H, G.sub(1), NROW, IZETA, W.sub(IW).contiguous(),
				DF, Z.sub(IZ).contiguous(), DMZ.sub(IDMZ).contiguous(),
				IPVTW.sub(IDMZ), 2, MODE, XI1,
				ZVAL.contiguous(), YVAL.contiguous(),
				FC.sub(IFC + INFC + NCOMP).contiguous(),
				DF, CB, IPVTCB, FC.sub(IFC), dfsub, ISING, NYCB);
			IZ = IZ + MSTAR;
			IDMZ = IDMZ + KDY;
			IW = IW + KDY * KDY;
			IFC = IFC + INFC + 2 * NCOMP;
			if (i < N)
				continue;

			while (true) {
				if (IZET > MSTAR)
					break;
				Z(IZ - 1 + IZET) = DGZ(IZET);
				IZET = IZET + 1;
			}
		}
		SBBLOK(G.contiguous(), INTEGS.contiguous(), N, IPVTG.contiguous(), Z.contiguous());

		//  finally find dmz
		DMZSOL(V.contiguous(), Z.contiguous(), DMZ.contiguous());
	}
	}
	//fmt::print(fg(fmt::color::orange_red), "Leave LSYSLV and RHS(41) = {}\n", RHS(41));
}



//**********************************************************************
//
//   purpose:
//
//      construct a collocation matrix row according to mode:
//      mode = 1  -  a row corresponding to a initial condition
//                   (i.e. at the left end of the subinterval).
//      mode = 2  -  a row corresponding to a condition at  aright.
//
//   variables:
//
//      gi     - the sub-block of the global bvp matrix in
//               which the equations are to be formed.
//      nrow   - no. of rows in gi.
//      irow   - the row in gi to be used for equations.
//      zval   - z(xi)
//      dg     - the derivatives of the side condition.
//
//**********************************************************************
void GDERIV(double* const GI, const int NROW, const int IROW,
	double const* ZVAL, double* DGZ, const int MODE, dgsub_t dgsub)
{
	//GI(NROW, 1);

	double DG[40];

	//  zero jacobian dg
	for (int j = 0; j < MSTAR; ++j)
		DG[j] = 0.0;

	//  evaluate jacobian dg
	dgsub(IZETA, ZVAL, DG);

	//  evaluate  dgz = dg * zval  once for a new mesh
	if (NONLIN != 0 && ITER <= 0) {
		double DOT = 0.0;
		for (int j = 0; j < MSTAR; ++j)
			DOT += DG[j] * ZVAL[j];
		DGZ[IZETA-1] = DOT;
	}

	//  branch according to  m o d e
	if (MODE != 2) {
		//  provide coefficients of the j-th linearized side condition.
		//  specifically, at x=zeta(j) the j-th side condition reads
		//  dg(1)*z(1) + ... +dg(mstar)*z(mstar) + g = 0

		//  handle an initial condition
		for (int j = 0; j < MSTAR; ++j) {
			GI[IROW - 1 + j * NROW] = DG[j];
			GI[IROW - 1 + (MSTAR + j) * NROW] = 0.0;
		}
	}
	else {
		//  handle a final condition
		for (int j = 0; j < MSTAR; ++j) {
			GI[IROW - 1 + j * NROW] = 0.0;
			GI[IROW - 1 + (MSTAR + j) * NROW] = DG[j];
		}
	}
}



//**********************************************************************
//
//   purpose:
//
//      construct a group of  ncomp  rows of the matrices  wi  and  vi.
//      corresponding to an interior collocation point.
//
//
//   variables:
//
//      xcol   - the location of the collocation point.
//      jj     - xcol is the jj-th of k collocation points
//               in the i-th subinterval.
//      wi,vi  - the i-th block of the collocation matrix
//               before parameter condensation.
//      kdy    - no. of rows in vi and wi .
//      zval   - z(xcol)
//      yval   - y(xcol)
//      df     - the jacobian at xcol .
//      jcomp  - counter for the component being dealt with.
//
//**********************************************************************
void VWBLOK(const double XCOL, const double HRHO, const int JJ, 
	double* const WI, double* const VI, int* const IPVTW,
	double const * const ZVAL, double const* const YVAL, 
	double* const DF, double const* const acol, double* const DMZO,
	dfsub_t dfsub, int& MSING)
{
	//WI(KDY, 1);
	//VI(KDY, 1);
	//DF(NCY, 1);
	//acol(7, 4);
	
	double BASM[5];
	dad2 HA(7, 4);
	
	// initialize  wi
	int I1 = (JJ - 1) * NCY;
	for (int ID = I1; ID < NCOMP + I1; ++ID)
		WI[ID + ID*KDY] = 1.0;

	//  calculate local basis
	double FACT = 1.0;
	for (int l = 0; l < MMAX; ++l) {
		FACT = FACT * HRHO / double(l+1);
		BASM[l] = FACT;
		for (int j = 0; j < K; ++j) {
			HA(j+1, l+1) = FACT * acol[j + l*7];
		}
	}

	// zero jacobian
	for (int JCOL = 0; JCOL < MSTAR + NY; ++JCOL) {
		for (int IR = 0; IR < NCY; ++IR) {
			DF[IR+ JCOL*NCY] = 0.0;
		}
	}
	
	//  build ncy rows for interior collocation point x.
	//  the linear expressions to be constructed are:
	//   (m(id))
	//  u     -  df(id,1)*z(1) - ... - df(id,mstar)*z(mstar) -
	//   id
	//        -  df(id,mstar+1)*u(1) - ... - df(id,mstar+ny)*y(ny)
	//  for id = 1 to ncy  (m(id)=0 for id > ncomp).
	dfsub(XCOL, ZVAL, YVAL, DF);
	int I0 = (JJ - 1) * NCY;
	I1 = I0 + 1;
	int I2 = I0 + NCY;

	// evaluate  dmzo = dmzo - df * (zval,yval)  once for a new mesh
	if (NONLIN != 0 && ITER <= 0) {
		for (int j = 0; j < MSTAR + NY; ++j) {
			if (j + 1 <= MSTAR)
				FACT = -ZVAL[j];
			else
				FACT = -YVAL[j - MSTAR];

			for (int ID = 0; ID < NCY; ++ID) {
				DMZO[I0 + ID] += FACT * DF[ID + j * NCY];
			}
		}
	}

	//  loop over the  ncomp  expressions to be set up for the
	//  current collocation point.
	for (int j = 0; j < MSTAR; ++j) {
		for (int ID = 0; ID < NCY; ++ID) {
			VI[I0 + ID + j * KDY] = DF[ID + j * NCY];
		}
	}
	int JN = 1;
	for (int JCOMP = 1; JCOMP <= NCOMP; ++JCOMP) {
		int MJ = MT(JCOMP);
		JN = JN + MJ;
		for (int l = 1; l <= MJ; ++l) {
			int JV = JN - l - 1;
			int JW = JCOMP - 1;
			for (int j = 1; j <= K; ++j) {
				double AJL = -HA(j, l);
				for (int IW = I1 - 1; IW < I2; ++IW) {
					WI[IW + JW * KDY] += AJL * VI[IW + JV * KDY];
				}
				JW = JW + NCY;
			}
			int LP1 = l + 1;
			if (l == MJ)
				continue;
			for (int LL = LP1; LL <= MJ; ++LL) {
				int JDF = JN - LL - 1;
				double BL = BASM[LL - l - 1];
				for (int IW = I1 - 1; IW < I2; ++IW)
					VI[IW + JV * KDY] += BL * VI[IW + JDF * KDY];
			}
		}
	}
	//  loop for the algebraic solution components
	for (int JCOMP = 0; JCOMP < NY; ++JCOMP) {
		int JD = NCOMP + JCOMP;
		for (int ID = 0; ID<NCY; ++ID)
			WI[(I0 + ID) + (I0 + JD) * KDY] = -DF[ID + (MSTAR + JCOMP) * NCY];
	}

	if (JJ < K)
		return;

	//...decompose the wi block and solve for the mstar columns of vi


	//  do parameter condensation
	MSING = dgefa(WI, KDY, KDY, IPVTW);

	//   check for singularity
	if (MSING != 0)
		return;
	for (int j = 0; j < MSTAR; ++j)
		dgesl(WI, KDY, KDY, IPVTW, VI+(j*KDY), 0);
}



//**********************************************************************
//
//   purpose:
//
//      construct projection matrices in  fc  in the general case
//      where the problem may have a higher index but is not in
//      a Hessenberg index-2 form
//
//      calls the linpack routine  dsvdc  for a singular value
//      decomposition.
//
//      mode = 1 - called from gblock
//           = 2 - called from newmsh
//                 (then fc consists of only ncomp columns)
//
//**********************************************************************
void PRJSVD(double* const FC, double const* const DF, double* const D, 
	double* const U, double* const V,
	int* const IPVTCB, int& ISING, const int MODE)
{
	//FC(NCOMP, 1);
	//DF(NCY, 1);
	//D(NY, NY);
	//U(NY, NY);
	//V(NY, NY);
	
	double WORK[20], S[21], E[20];


	//  compute the maximum tolerance
	double CHECK = 0.0;
	for (int i = 1; i <= NTOL; ++i)
		CHECK = std::max(TOLIN(i), CHECK);

	//  construct d and find its svd
	for (int i = 0; i < NY; ++i)
		for (int j = 0; j < NY; ++j)
			D[i + j*NY] = DF[i + NCOMP + (j + MSTAR)*NCY];

	int JOB = 11;
	int INFO = dsvdc(D, NY, NY, NY, S, E, U, NY, V, NY, WORK, JOB);

	//  determine rank of d
	S[NY] = 0;
	int IRANK = 0;

	while (S[IRANK] >= CHECK){
		IRANK++;
	}

	//  if d has full rank then no projection is needed
	if (IRANK == NY) {
		for (int i = 0; i < NCOMP; ++i)
			for (int j = 0; j < MSTAR + NY; ++j)
				FC[i + j*NCY] = 0.0;
		return;
	}
	else{
		//  form projected cb
		int IR = NY - IRANK;
		for (int i = 0; i < NY; ++i) {
			for (int j = 0; j < NY; ++j) {
				double FACT = 0;
				int ML = 0;
				for (int l = 0; l < NCOMP; ++l) {
					ML += MT(l + 1);
					FACT += DF[i + NCOMP + (ML - 1) * NCY] * DF[l + (MSTAR + j) * NCY];
				}
				D[i + j * NY] = FACT;
			}
		}
		for (int i = 0; i < NY; ++i) {
			for (int j = 0; j < IR; ++j) {
				WORK[j] = 0;
				for (int l = 0; l < NY; ++l) {
					WORK[j] += D[i + l * NY] * V[l + (j + IRANK) * NY];
				}
			}
			for (int j = 0; j < NCOMP; ++j) {
				D[i + j * NY] = WORK[j];
			}
		}
		for (int i = 0; i < IR; ++i) {
			for (int j = 0; j < IR; ++j) {
				WORK[j] = 0;
				for (int l = 0; l < NY; ++l) {
					WORK[j] += U[l + (i + IRANK) * NY] * D[l + j * NY];
				}
			}
			for (int j = 0; j < IR; ++j) {
				D[i + j * NY] = WORK[j];
			}
		}
		//  decompose projected cb
		ISING = dgefa(D, NY, IR, IPVTCB);
		if (ISING != 0)
			return;

		//  form columns of fc
		for (int j = MSTAR; j < MSTAR + NY; ++j) {
			for (int i = 0; i < IR; ++i)
				WORK[i] = U[(j - MSTAR) + (i + IRANK)*NY];
			dgesl(D, NY, IR, IPVTCB, WORK, 0);
			for (int i = 0; i < NY; ++i) {
				U[(j - MSTAR) + i*NY] = 0;
				for (int l = 0; l < IR; ++l)
					U[(j - MSTAR) + i * NY] += V[i + (l + IRANK) * NY] * WORK[l];
			}
			for (int i = 0; i < NCOMP; ++i) {
				double FACT = 0;
				for (int l = 0; l < NY; ++l)
					FACT += DF[i + (MSTAR + l) * NCY] * U[j - MSTAR + l * NY];
				FC[i + j * NCOMP] = FACT;
			}
		}
		if (MODE == 1) {
			for (int i = 0; i < NCOMP; ++i) {
				for (int j = 0; j < MSTAR; ++j) {
					double FACT = 0;
					for (int l = 0; l < NY; ++l)
						FACT += FC[i + (l + MSTAR) * NCOMP] * DF[l + NCOMP + j * NCY];
					FC[i + j * NCOMP] = FACT;
				}
			}
		}
		else {
			for (int i = 0; i < NCOMP; ++i){
				int MJ = 0;
				for (int j = 0; j < NCOMP; ++j) {
					MJ += MT(j+1);
					double FACT = 0;
					for (int l = 0; l < NY; ++l)
						FACT += FC[i + (l + MSTAR) * NCOMP] * DF[l + NCOMP + (MJ-1) * NCY];
					FC[i + j * NCOMP] = FACT;
				}
			}
		}

	}
}


//**********************************************************************
//
//purpose:
//
//construct collocation matrix rows according to mode :
//mode = 1 - a group of  mstar    rows corresponding
//to a mesh interval.
//= 2 - compute condensed form of rhs
//modl = mode of lsyslv
//
//variables :
//
//h - the  local stepsize.
//gi - the sub - block of the collocation matrix in
//which the equations are to be formed.
//wi - the sub - block of noncondensed collocation equations,
//left - hand side part.
//vi - the sub - block of noncondensed collocation equations,
//right - hand side part.
//rhsdmz - the inhomogenous term of the uncondensed collocation
//equations.
//rhsz - the inhomogenous term of the condensed collocation
//equations.
//nrow - no.of rows in gi.
//irow - the first row in gi to be used for equations.
//xi1 - next mesh point(right end)
//zval, yval - current solution at xi1
//fc - projection matrices
//
//**********************************************************************
// double const * WI (dgesl)
void GBLOCK(const double H, dar2 GI, const int NROW, const int IROW,
	double* const  WI, dar2 VI,
	double* const RHSZ, double* const RHSDMZ,
	iar1 IPVTW, const int MODE,
	const int MODL, const double XI1, double const* const  ZVAL,
	double const* const  YVAL, double* const  F,
	dar2 DF, dar2 CB, iar1 IPVTCB, dar2 FC, dfsub_t dfsub,
	int& ISING, const int NYCB)
{	
	AutoTimer at(g_timer, _FUNC_);

	GI.reshape(NROW, 1); 
	VI.reshape(KDY, 1);
	DF.reshape(NCY, 1); 
	CB.reshape(NYCB, NYCB);
	FC.reshape(NCOMP, 1);
	
	dad2 HB(7, 4);
	double BASM[5], BCOL[40], U[400], V[400];

	//  compute local basis
	double FACT = 1.0;
	BASM[0] = 1.0;
	for (int l = 1; l <= MMAX; ++l) {
		FACT = FACT * H / double(l);
		BASM[l] = FACT;
		for (int j = 1; j <= K; ++j)
			HB(j, l) = FACT * B(j, l);
	}

	//  branch according to  m o d e
	switch (MODE) {
	case 1:
		//  set right gi-block to identity
		if (MODL != 2)
		{
			for (int j = 1; j <= MSTAR; ++j) {
				for (int IR = 1; IR <= MSTAR; ++IR) {
					GI(IROW - 1 + IR, j) = 0.0;
					GI(IROW - 1 + IR, MSTAR + j) = 0.0;
				}
				GI(IROW - 1 + j, MSTAR + j) = 1.0;
			}
			//  compute the block gi
			int IR = IROW;
			for (int ICOMP = 1; ICOMP <= NCOMP; ++ICOMP) {
				int MJ = MT(ICOMP);
				IR = IR + MJ;
				for (int l = 1; l <= MJ; ++l) {
					int ID = IR - l;
					for (int JCOL = 1; JCOL <= MSTAR; ++JCOL) {
						int IND = ICOMP;
						double RSUM = 0.0;
						for (int j = 1; j <= K; ++j) {
							RSUM = RSUM - HB(j, l) * VI(IND, JCOL);
							IND = IND + NCY;
						}
						GI(ID, JCOL) = RSUM;
					}
					int JD = ID - IROW;
					for (int LL = 1; LL <= l; ++LL)
						GI(ID, JD + LL) = GI(ID, JD + LL) - BASM[LL-1];

				}
			}

			if (INDEX == 1 || NY == 0)
				return;

			//  projected collocation
			//  set up projection matrix and update gi-block
			dfsub(XI1, ZVAL, YVAL, DF.contiguous());

			//  if index=2 then form projection matrices directly
			//  otherwise use svd to define appropriate projection
			if (INDEX == 0) {
				PRJSVD(FC.contiguous(), DF.contiguous(), CB.contiguous(),
					U, V, IPVTCB.contiguous(), ISING, 1); // TODO warum haben wir U,V?
				if (ISING != 0)
					return;
			}
			else {
				//  form  cb
				for (int i = 1; i <= NY; ++i) {
					for (int j = 1; j <= NY; ++j) {
						FACT = 0;
						int ML = 0;
						for (int l = 1; l <= NCOMP; ++l) {
							ML = ML + MT(l);
							FACT = FACT + DF(i + NCOMP, ML) * DF(l, MSTAR + j);
						}
						CB(i, j) = FACT;
					}
				}

				//  decompose cb
				DGEFA(CB, NY, NY, IPVTCB, ISING);
				if (ISING != 0)
					return;

				//  form columns of fc
				for (int j = 1; j <= MSTAR + NY; ++j) {

					if (j <= MSTAR)
						for (int i = 1; i<= NY;++i)
							BCOL[i-1] = DF(i + NCOMP, j);
					else {
						for (int i = 0; i < NY; ++i)
							BCOL[i] = 0.0;
						BCOL[j - MSTAR - 1] = 1.0;
					}

					dgesl(CB.contiguous(), NY, NY, IPVTCB.contiguous(), BCOL, 0);
					
					for (int i = 1; i <= NCOMP; ++i) {
						FACT = 0.0;
						for (int l = 1; l <= NY; ++l)
							FACT = FACT + DF(i, l + MSTAR) * BCOL[l-1];

						FC(i, j) = FACT;
					}
				}
			}

			//  update gi
			for (int j = 1; j <= MSTAR; ++j) {
				for (int i = 1; i <= NCOMP; ++i) {
					FACT = 0;
					for (int l = 1; l <= MSTAR; ++l)
						FACT = FACT + FC(i, l) * GI(IROW - 1 + l, j);

					BCOL[i-1] = FACT;
				}
				int ML = 0;
				for (int i = 1; i <= NCOMP; ++i) {
					ML = ML + MT(i);
					GI(IROW - 1 + ML, j) = GI(IROW - 1 + ML, j) - BCOL[i-1];
				}
			}
		}

		//  prepare extra rhs piece; two if new mesh
		if (INDEX == 1 || NY == 0)
			return;
		for (int JCOL = 1; JCOL <= 2; ++JCOL) {
			for (int i = 1; i <= NCOMP; ++i) {
				FACT = 0;
				for (int l = 1; l <= NY; ++l)
					FACT = FACT + FC(i, l + MSTAR) * F[l + NCOMP-1];
				FC(i, JCOL + MSTAR + NY) = FACT;
			}

			if (MODL != 1 || JCOL == 2)
				return;
			for (int i = 1 + NCOMP; i <= NY + NCOMP; ++i)
				F[i - 1] = 0;
			for (int j = 1; j <= MSTAR;++j) {
				FACT = -ZVAL[j-1];
				for (int i = 1 + NCOMP; i <= NY + NCOMP; ++i)
					F[i - 1] += DF(i, j) * FACT;
			}
		}
		return;

	case 2:
		//  compute the appropriate piece of  rhsz
		dgesl(WI, KDY, KDY, IPVTW.contiguous(), RHSDMZ, 0);
	
		int IR = IROW;
		for (int JCOMP = 1; JCOMP <= NCOMP; ++JCOMP) {
			int MJ = MT(JCOMP);
			IR = IR + MJ;
			for (int l = 1; l <= MJ; ++l) {
				int	IND = JCOMP;
				double RSUM = 0.0;
				for (int j = 1; j <= K; ++j) {
					RSUM += HB(j, l) * RHSDMZ[IND-1];				
					IND += NCY;
				}
				RHSZ[IR - l-1] = RSUM;
			}
		}

		if (INDEX == 1 || NY == 0)
			return;

		//  projected collocation
		//  calculate projected rhsz
		for (int i = 1; i <= NCOMP; ++i) {
			FACT = 0;
			for (int l = 1; l <= MSTAR; ++l)
				FACT += FC(i, l) * RHSZ[l + IROW - 2];
			BCOL[i-1] = FACT;
		}
		int ML = 0;
		for (int i = 1; i <= NCOMP; ++i) {
			ML += MT(i);
			RHSZ[IROW - 1 + ML-1] -= BCOL[i-1] + F[i-1];
		}	

	}
}




//----------------------------------------------------------------------
//                             p a r t  4
//               polynomial and service routines
//----------------------------------------------------------------------


/***********************************************************************

purpose
		evaluate mesh independent runge - kutta basis for given s

variables
	s - argument, i.e.the relative position for which
			the basis is to be evaluated(0..le.s.le. 1.).
	coef - precomputed derivatives of the basis
	k - number of collocatin points per subinterval
	m - maximal order of the differential equation
	rkb - the runge - kutta basis(0 - th to(m - 1) - th derivatives)
	dm - basis elements for m - th derivative
			these are evaluated if mode > 0.

* **********************************************************************/
void RKBAS(double const S, int const k, int const M, double* const RKB,
	double* const DM, int const MODE)
{
	//COEF(k, k);
	//RKB(7, 1);

	AutoTimer at(g_timer, _FUNC_);

	if (k != 1) {
		double T[10];
		for (int i = 0; i < k + M - 1; ++i)
			T[i] = S / double(i + 1);
		for (int l = 0; l < M; ++l) {
			int LB = k + l + 2;
			for (int i = 0; i < k; ++i) {
				double P = COEF(1, i+1);
				for (int j = 2; j <= k; ++j)
					P = P * T[LB - j - 1] + COEF(j, i+1);
				RKB[i + l*7] = P;
			}
		}
		if (MODE == 0)
			return;
		for (int i = 0; i < k; ++i) {
			double P = COEF(1, i+1);
			for (int j = 2; j <= k; ++j)
				P = P * T[k - j] + COEF(j, i + 1);
			DM[i] = P;
		}
		return;
	}
	RKB[0] = 1.0;
	DM[0] = 1.0;
}




//**********************************************************************
//
// purpose
//                                  (1)       (m1-1)     (mncomp-1)
//         evaluate z(u(x))=(u (x),u (x),...,u  (x),...,u  (x)      )
//                            1     1         1          mncomp
//         as well as optionally y(x) and dmval(x) at one point x.
//
// variables
//   a      - array of mesh independent rk-basis coefficients
//   xi     - the current mesh (having n subintervals)
//   z      - the current solution vector (differential components).
//            it is convenient to imagine z as a two-dimensional
//            array with dimensions mstar x (n+1). then
//            z(j,i) = the jth component of z at the ith mesh point
//   dmz    - the array of mj-th derivatives of the current solution
//            plus algebraic solution components at collocation points
//            it is convenient to imagine dmz as a 3-dimensional
//            array with dimensions ncy x k x n. then
//            dmz(l,j,i) = a solution value at the jth collocation
//            point in the ith mesh subinterval: if l <= ncomp then
//            dmz(l,j,i) is the ml-th derivative of ul, while if
//            l > ncomp then dmz(l,j,i) is the value of the current
//            (l-ncomp)th component of y at this collocation point
//   mode   - determines the amount of initialization needed
//          = 4  forms z(u(x)) using z, dmz and ha
//          = 3  as in =4, but computes local rk-basis
//          = 2  as in =3, but determines i such that
//                     xi(i) .le. x .lt. xi(i+1) (unless x=xi(n+1))
//          = 1  retrieve  z=z(u(x(i)))  directly
//   modm   = 0  evaluate only zval
//          = 1  evaluate also yval
//          = 2  evaluate in addition dmval
// output
//   zval   - the solution vector z(u(x)) (differential components)
//   yval   - the solution vector y(x)  (algebraic components)
//   dmval  - the mth derivatives of u(x)
//
//**********************************************************************
void APPROX(int& i, double& X, double* const ZVAL, double* const YVAL,
	double* const A, double const* const coef, double const* const XI,
	const int n, double const* const Z, double const * const DMZ, const int k, const int ncomp,
	const int ny, const int mmax, int const * const M,
	const int mstar, const int MODE, double* const DMVAL, const int MODM)
{
	//A(7, 1);

	AutoTimer at(g_timer, _FUNC_);

	dad1 BM(4), DM(7);
	int IZ, ILEFT, IRIGHT;

	switch (MODE) {
	case 1:
		//  mode = 1, retrieve  z(u(x))  directly for x = xi(i).
		X = XI[i-1];
		IZ = (i - 1) * mstar;
		for (int j = 0; j < mstar; ++j) {
			IZ++;
			ZVAL[j] = Z[IZ-1];
		}
		return;

	case 2:
		//  mode = 2, locate i so  xi(i).le.x.lt.xi(i + 1)
		if (X < XI[0] - PRECIS || X > XI[n] + PRECIS)
		{
			if (IPRINT < 1)
				fmt::print(fg(fmt::color::red),"****** DOMAIN ERROR IN APPROX ******\n"
					" X = {}, ALEFT = {}, ARIGHT = {}\n",
					X, XI[0], XI[n]);
			if (X < XI[0])
				X = XI[0];
			if (X > XI[n])
				X = XI[n];
		}
		if (i > n || i < 1)
			i = (n + 1) / 2;
		ILEFT = i;
		if (X >= XI[ILEFT-1]) {
			for (int l = ILEFT; l <= n; ++l) {
				i = l;
				if (X < XI[l])
					break;
			}
		}
		else {
			IRIGHT = ILEFT - 1;
			for (int l = 1; l <= IRIGHT; ++l) {
				i = IRIGHT + 1 - l;
				if (X >= XI[i-1])
					break;
			}
		}
		[[fallthrough]];
	case 3:	 {
		//  mode = 2 or 3, compute mesh independent rk - basis.
		double S = (X - XI[i-1]) / (XI[i] - XI[i-1]);
		RKBAS(S, k, mmax, A, DM.contiguous(), MODM);
		}
		[[fallthrough]];
	case 4:
		//  mode = 2, 3, or 4, compute mesh dependent rk - basis.
		BM(1) = X - XI[i-1];

		for (int l = 2; l <= mmax; ++l)
			BM(l) = BM(1) / double(l);

		//  evaluate  z(u(x)).
		int IR = 1;
		NCY = ncomp + ny;
		IZ = (i - 1) * mstar + 1;
		int IDMZ = (i - 1) * k * NCY;
		for (int JCOMP = 1; JCOMP <= ncomp; ++JCOMP) {
			int MJ = M[JCOMP-1];
			IR = IR + MJ;
			IZ = IZ + MJ;
			for (int l = 1; l <= MJ; ++l) {
				int IND = IDMZ + JCOMP;
				double ZSUM = 0.0;
				for (int j = 1; j <= k; ++j) {
					ZSUM = ZSUM + A[(j - 1) + (l-1)*7] * DMZ[IND - 1];
					IND = IND + NCY;
				}
				for (int LL = 1; LL <= l; ++LL) {
					int LB = l + 1 - LL;
					ZSUM = ZSUM * BM(LB) + Z[IZ - LL-1];
				}
				ZVAL[IR - l-1] = ZSUM;
			}
		}
		if (MODM == 0)
			return;

		//  for modm = 1 evaluate  y(j) = j - th component of y.
		for (int JCOMP = 1; JCOMP <= ny; ++JCOMP)
			YVAL[JCOMP-1] = 0.0;
		for (int j = 1; j <= k; ++j) {
			int IND = IDMZ + (j - 1) * NCY + ncomp + 1;
			double FACT = DM(j);
			for (int JCOMP = 1; JCOMP <= ny; ++JCOMP) {
				YVAL[JCOMP - 1] += FACT * DMZ[IND-1];
				IND = IND + 1;
			}
		}
		if (MODM == 1)
			return;

		//  for modm = 2 evaluate  dmval(j) = mj - th derivative of uj.
		for (int JCOMP = 1; JCOMP <= ncomp; ++JCOMP)
			DMVAL[JCOMP - 1] = 0.0;
		for (int j = 1; j <= k; ++j) {
			int IND = IDMZ + (j - 1) * NCY + 1;
			double FACT = DM(j);
			for (int JCOMP = 1; JCOMP <= ncomp; ++JCOMP) {
				DMVAL[JCOMP - 1] += FACT * DMZ[IND-1];
				IND = IND + 1;
			}
		}
	}
}




/* * *********************************************************************

purpose
		solve vandermonde system v * x = e
		with  v(i, j) = rho(j) * *(i - 1) / (i - 1)!.

* **********************************************************************/
void VMONDE(double* const coef, int const k)
{
	RHO.assertDim(k);
	//coef(k);
	
	AutoTimer at(g_timer, _FUNC_);

	if (k == 1)
		return;
	int KM1 = k - 1;
	for (int i = 1; i <= KM1; ++i) {
		int KMI = k - i;
		for (int j = 1; j <= KMI; ++j) {
			coef[j-1] = (coef[j] - coef[j-1]) / (RHO(j + i) - RHO(j));
		}
	}

	int IFAC = 1;
	for (int i = 1; i <= KM1; ++i) {
		int KMI = k + 1 - i;
		for (int j = 2; j <= KMI; ++j)
			coef[j-1] -= RHO(j + i - 1) * coef[j - 2];
		coef[KMI-1] *= double(IFAC);
		IFAC *= i;
	}
	coef[0] *= double(IFAC);
}



/* * *********************************************************************

purpose
		determine highest order(piecewise constant) derivatives
		of the current collocation solution

variables
	hi - the stepsize, hi = xi(i + 1) - xi(i)
	dmz - vector of mj - th derivative of the solution
	uhigh - the array of highest order(piecewise constant)
			derivatives of the approximate solution on
(xi(i), xi(i + 1)), viz,
(k + mj - 1)
			uhigh(j) = u(x)    on(xi(i), xi(i + 1))
						j

* **********************************************************************/
void HORDER(int const i, double* const UHIGH, double const HI, double const * const DMZ)
{
	AutoTimer at(g_timer, _FUNC_);
	
	double DN = 1.0 / pow(HI, (K - 1));

	//  loop over the ncomp solution components
	for (int ID = 0; ID < NCOMP; ++ID)
		UHIGH[ID] = 0.0;

	int KIN = 1;
	int IDMZ = (i - 1) * K * NCY;
	for (int j = 1; j <= K;++j) {
		double FACT = DN * COEF.contiguous()[KIN-1];
		for (int ID = 0; ID < NCOMP;++ID) {
			UHIGH[ID] += FACT * DMZ[IDMZ];
			IDMZ++;
		}
		KIN += K;
	}
}



/* * *********************************************************************

purpose
		compute dmz in a blockwise manner
		dmz(i) = dmz(i) + v(i) * z(i), i = 1, ..., n

* **********************************************************************/
void DMZSOL(double* V, double* Z, double* DMZ)
{
	//V(KDY, 1);
	//DMZ(KDY, 1);
	
	AutoTimer at(g_timer, _FUNC_);
	
	int JZ = 0;
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < MSTAR; ++j) {
			double FACT = Z[JZ];
			for (int l = 0; l < KDY; ++l) {
				DMZ[l + i*KDY] += FACT * V[l + JZ*KDY];
			}
			JZ++;
		}
	}
}




/*----------------------------------------------------------------------
	                        p a r t  5
	        we list here a modified(column oriented, faster)
	        version of the package solveblok of de boor - weiss[7].
	        we also give a listing of the linpack
	        routines dgefa und dgesl used by coldae.
----------------------------------------------------------------------*/




//calls subroutines  factrband shiftb .
//	
//	     fcblok  supervises the plu factorization with pivoting of
//	     scaled rows of the almost block diagonal matrix stored in the
//	     arrays  bloks and integs .
//	
//	     factrb = subprogram which carries out steps 1, ..., last of gauss
//	            elimination(with pivoting) for an individual block.
//	     shiftb = subprogram which shifts the remaining rows to the top of
//	            the next block
//	
//	     parameters
//	      bloks   an array that initially contains the almost block diago -
//	            nal matrix  a  to be factored, and on return contains the
//	            computed factorization of  a .
//	      integs  an integer array describing the block structure of  a .
//	      nbloks  the number of blocks in  a .
//	      ipivot  an integer array of dimension   sum(integs(3, n); n = 1,
//		            ..., nbloks) which, on return, contains the pivoting stra -
//	            tegy used.
//	      scrtch  work area required, of length  max(integs(1, n); n = 1,
//		            ..., nbloks).
//	      info    output parameter;
// = 0  in case matrix was found to be nonsingular.
//            otherwise,
// = n if the pivot element in the nth gauss step is zero.
//
//**********************************************************************
void FCBLOK(double* const BLOKS, int const * const INTEGS, 
	int const NBLOKS, int* const IPIVOT, double* const SCRTCH, int& INFO)
{
	AutoTimer at(g_timer, _FUNC_);
	//INTEGS.assertDim(3, NBLOKS);
	/*IPIVOT.assertDim(1);
	BLOKS.assertDim(1);
	SCRTCH.assertDim(1);*/

	INFO = 0;
	int INDEXX = 1;
	int INDEXN = 1;
	
	//  loop over the blocks. i is loop index
	int i = 1;
	while (true) {
		int INDEX = INDEXN;  // soll ausblenden
		int NROW = INTEGS[0 + 3 * (i - 1)];
		int NCOL = INTEGS[1 + 3 * (i - 1)];
		int LAST = INTEGS[2 + 3 * (i - 1)];

		// carry out elimination on the i - th block until next block
		// enters, i.e., for columns 1, ..., last  of i - th block.
		FACTRB(BLOKS+(INDEX-1), IPIVOT+(INDEXX-1), SCRTCH, NROW, NCOL, LAST, INFO);

		// check for having reached a singular block or the last block
		if (INFO != 0)
			break;
		if (i == NBLOKS)
			return;

		i = i + 1;
		INDEXN = NROW * NCOL + INDEX;
		INDEXX = INDEXX + LAST;

		// put the rest of the i - th block onto the next block
		SHIFTB(BLOKS + (INDEX - 1), NROW, NCOL, LAST, BLOKS + (INDEXN - 1),
			INTEGS[0 + (i - 1) * 3], INTEGS[1 + (i - 1) * 3]);
	}
	INFO = INFO + INDEXX - 1;
}


//
//********************************************************************
//
//     adapted from p.132 of  element.numer.analysis  by conte - de boor
//
//     constructs a partial plu factorization, corresponding to steps
//      1, ..., last   in gauss elimination, for the matrix  w  of
//      order(nrow, ncol), using pivoting of scaled rows.
//
//     parameters
//       w       contains the(nrow, ncol) matrix to be partially factored
//               on input, and the partial factorization on output.
//       ipivot  an integer array of length last containing a record of
//               the pivoting strategy used; explicit interchanges
//               are used for pivoting.
//       d       a work array of length nrow used to store row sizes
//               temporarily.
//       nrow    number of rows of w.
//       ncol    number of columns of w.
//       last    number of elimination steps to be carried out.
//       info    on output, zero if the matrix is found to be non -
//               singular, in case a zero pivot was encountered in row
//               n, info = n on output.
//
// * *********************************************************************
//
void FACTRB(double* const W, int* const IPIVOT, double* const D,
            const int NROW, const int NCOL, const int LAST, int& INFO)
{
	//IPIVOT(NROW);
	//W(NROW, NCOL);
	//D(NROW);
	
	AutoTimer at(g_timer, _FUNC_);
	
	double COLMAX, T, S;
	int k, l, KP1;


	//  initialize  d
	for (int i = 0; i < NROW; ++i)
		D[i] = 0.0;

	for (int j = 0; j < NCOL; ++j)
		for (int i = 0; i < NROW; ++i)
			D[i] = std::max(D[i], abs(W[i + j*NROW]));

	//  gauss elimination with pivoting of scaled rows, loop over
	//  k = 1, ., last
	//  as pivot row for k - th step, pick among the rows not yet used,
	//  i.e., from rows  k, ..., nrow, the one whose k - th entry
	//  (compared to the row size) is largest.then, if this row
	//  does not turn out to be row k, interchange row k with this
	//  particular rowand redefine ipivot(k).

	k = 1;
	while (true) {

		if (D[k-1] == 0.0) {
			INFO = k;
			return;
		}
		if (k == NROW) {
			// if  last.eq.nrow, check now that pivot element in last row is nonzero.
			if (abs(W[NROW-1 + (NROW-1)* NROW]) <= 0)
				INFO = k;
		}

		l = k;
		KP1 = k + 1;
		COLMAX = abs(W[(k-1)+(k-1)*NROW]) / D[k-1];
		// find the (relatively) largest pivot
		for (int i = KP1; i <= NROW; ++i) {
			if (abs(W[(i - 1) + (k - 1) * NROW]) <= COLMAX * D[i - 1])
				continue;
			COLMAX = abs(W[(i - 1) + (k - 1) * NROW]) / D[i - 1];
			l = i;
		}

		IPIVOT[k - 1] = l;
		T = W[(l - 1) + (k - 1) * NROW];
		S = D[l-1];
		if (l != k) {
			W[(l - 1) + (k - 1) * NROW] = W[(k - 1) + (k - 1) * NROW];
			W[(k - 1) + (k - 1) * NROW] = T;
			D[l - 1] = D[k - 1];
			D[k - 1] = S;
		}

		// if pivot element is too small in absolute value, declare
		// matrix to be noninvertible and quit.
		if (abs(T) <= 0)
		{
			INFO = k; //  singularity flag set
			return;
		}

		// otherwise, subtract the appropriate multiple of the pivot
		// row from remaining rows, i.e., the rows(k + 1), ..., (nrow)
		// to make k - th entry zero. save the multiplier in its place.
		// for high performance do this operations column oriented.
		T = -1.0 / T;
		for (int i = KP1; i <= NROW; ++i)
			W[(i - 1) + (k - 1) * NROW] *= T;

		for (int j = KP1; j <= NCOL; ++j) {
			T = W[(l - 1) + (j - 1) * NROW];
			if (l != k) {
				W[(l - 1) + (j - 1) * NROW] = W[(k - 1) + (j - 1) * NROW];
				W[(k - 1) + (j - 1) * NROW] = T;
			}
			if (T == 0.0)
				continue;
			for (int i = KP1; i <= NROW; ++i)
				W[(i - 1) + (j - 1) * NROW] += W[(i - 1) + (k - 1) * NROW] * T;
		}
		k = KP1;

		// check for having reached the next block.
		if (k > LAST)
			break;
	}
}



//
//*********************************************************************
//
//     shifts the rows in current block, ai, not used as pivot rows, if
//     any, i.e., rows(last + 1), ..., (nrowi), onto the first mmax =
// = nrow - last  rows of the next block, ai1, with column last + j of
//      ai  going to column j, j = 1, ..., jmax = ncoli - last.the remaining
//     columns of these rows of ai1 are zeroed out.
//
//                                picture
//
//          original situation after         results in a new block i + 1
//          last = 2 columns have been       created and ready to be
//          done in factrb(assuming no      factored by next factrb
//          interchanges of rows)            call.
//                      1
//                 x  x 1x  x  x           x  x  x  x  x
//                      1
//                 0  x 1x  x  x           0  x  x  x  x
//     block i          1                       -------------- -
//     nrowi = 4   0  0 1x  x  x           0  0 1x  x  x  0  01
//     ncoli = 5        1                       1             1
//     last = 2    0  0 1x  x  x           0  0 1x  x  x  0  01
//------------------------------ - 1             1   new
//                      1x  x  x  x  x          1x  x  x  x  x1  block
//                      1                       1             1   i + 1
//     block i + 1        1x  x  x  x  x          1x  x  x  x  x1
//     nrowi1 = 5        1                       1             1
//     ncoli1 = 5        1x  x  x  x  x          1x  x  x  x  x1
//------------------------------ - 1------------ - 1
//                      1
//
// * ********************************************************************
//
void SHIFTB(double* const AI, const int NROWI, const int NCOLI, const int LAST,
            double* const AI1, const int NROWI1, const int NCOLI1)
{
	//AI(NROWI, NCOLI);
	//AI1(NROWI1, NCOLI1)
	
	AutoTimer at(g_timer, _FUNC_);
	
	int MMAX = NROWI - LAST; // soll ausblenden
	int JMAX = NCOLI - LAST;
	if (MMAX < 1 || JMAX < 1)
		return;

	// put the remainder of block i into ai1
	for (int j = 0; j < JMAX; ++j)
		for (int m = 0; m < MMAX; ++m)
			AI1[m + j * NROWI1] = AI[(LAST + m) + (LAST + j) * NROWI];

	if (JMAX == NCOLI1)
		return;

	// zero out the upper right corner of ai1
	for (int j = JMAX; j < NCOLI1; ++j)
		for (int m = 0; m < MMAX; ++m)
			AI1[m + j * NROWI1] = 0.0;
}




//
//**********************************************************************
//
//     calls subroutines  subforand subbak .
//
//     supervises the solution(by forward and backward substitution) of
//     the linear system  a* x = b  for x, with the plu factorization of
//     a  already generated in  fcblok.individual blocks of
//     equations are solved via  subforand subbak .
//
//    parameters
//       bloks, integs, nbloks, ipivot    are as on return from fcblok.
//       x       on input : the right hand side, in dense storage
//               on output : the solution vector
//
//*********************************************************************
//
void SBBLOK(double* const BLOKS, int const * const INTEGS, 
	int const NBLOKS, int const * const IPIVOT, double* const X)
{
	// double BLOKS[?]
	// int INTEGS[3][NBLOKS]
	// int IPIVOT[?]
	// double X[?]

	AutoTimer at(g_timer, _FUNC_);
	
	
	int NCOL, NROW, LAST;

	//  forward substitution pass
	int INDEX = 1;  // soll ausblenden
	int INDEXX = 1;
	for (int i = 1; i <= NBLOKS; ++i) {
		NROW = INTEGS[0 + 3 * (i - 1)];
		NCOL = INTEGS[1 + 3 * (i - 1)];
		LAST = INTEGS[2 + 3 * (i - 1)];
		SUBFOR(BLOKS+(INDEX-1), IPIVOT+(INDEXX-1), NROW, LAST, X+(INDEXX-1));
		INDEX += NROW * NCOL;
		INDEXX += LAST;
	}

	//  back substitution pass
	int NBP1 = NBLOKS + 1;
	for (int j = 1; j <= NBLOKS; ++j) {
		int i = NBP1 - j;
		NROW = INTEGS[0 + 3 * (i - 1)];
		NCOL = INTEGS[1 + 3 * (i - 1)];
		LAST = INTEGS[2 + 3 * (i - 1)];
		INDEX -= NROW * NCOL;
		INDEXX -= LAST;
		SUBBAK(BLOKS+(INDEX-1), NROW, NCOL, LAST, X+(INDEXX-1));
	}
}



//
//**********************************************************************
//
//     carries out the forward pass of substitution for the current
//     block, i.e., the action on the right side corresponding to the
//     elimination carried out in  factrb  for this block.
//
//    parameters
//       w, ipivot, nrow, last  are as on return from factrb.
//       x(j)  is expected to contain, on input, the right side of j - th
//             equation for this block, j = 1, ..., nrow.
//       x(j)  contains, on output, the appropriately modified right
//             side of equation(j) in this block, j = 1, ..., lastand
//             for j = last + 1, ..., nrow.
//
//*********************************************************************
//
void SUBFOR(double const * const W, int const * const IPIVOT, int const NROW,
	int const LAST, double* const X)
{
	// double W[NROW, LAST]
	// double X[NROW]
	// int IPIVOT[LAST]

	AutoTimer at(g_timer, _FUNC_);

	if (NROW == 1)
		return;

	int LSTEP = std::min(NROW - 1, LAST);
	for (int k = 1; k <= LSTEP; ++k)
	{
		int IP = IPIVOT[k - 1];
		double T = X[IP - 1];
		X[IP - 1] = X[k - 1];
		X[k - 1] = T;
		if (T == 0.0)
			continue;
		for (int i = k + 1; i <= NROW; ++i)
			X[i - 1] = X[i - 1] + W[(i - 1) + (k - 1) * NROW] * T;
	}
}



//
//*********************************************************************
//
//     carries out backsubstitution for current block.
//
//    parameters
//       w, ipivot, nrow, ncol, last  are as on return from factrb.
//       x(1), ..., x(ncol)  contains, on input, the right side for the
//               equations in this block after backsubstitution has been
//               carried up to but not including equation(last).
//               means that x(j) contains the right side of equation(j)
//               as modified during elimination, j = 1, ..., last, while
//               for j.gt.last, x(j) is already a component of the
//               solution vector.
//       x(1), ..., x(ncol) contains, on output, the components of the
//               solution corresponding to the present block.
//
//*********************************************************************
//

void SUBBAK(double* const W, int const NROW, int const NCOL, int const LAST, double* const X)
{
	// W(NROW, NCOL), X(NCOL)
	AutoTimer at(g_timer, _FUNC_);
	
	int LP1 = LAST + 1;
	if (LP1 <= NCOL) {
		for (int j = LP1; j <= NCOL; ++j) {
			double T = -X[j - 1];
			if (T == 0.0)
				continue;
			for (int i = 1; i <= LAST; ++i)
				X[i - 1] += W[(i - 1) + (j - 1) * NROW] * T;
		}
	}

	if (LAST != 1) {
		int LM1 = LAST - 1;
		for (int KB = 1; KB <= LM1; ++KB) {
			int KM1 = LAST - KB;
			int k = KM1 + 1;
			X[k - 1] /= W[(k - 1) + (k - 1) * NROW];
			double T = -X[k - 1];
			if (T == 0.0)
				continue;
			for (int i = 1; i <= KM1; ++i)
				X[i-1] += W[(i - 1) + (k - 1) * NROW] * T;
		}
	}
	X[0] /= W[0];
}

};

