#pragma once

#include "test/linpack/linpack_d.hpp"
#include "test/timer.h"

#define FMT_HEADER_ONLY 
#include "fmt/include/fmt/format.h"
#include "fmt/include/fmt/ranges.h"
#include "fmt/include/fmt/color.h"

//------------------------------------------------------------------------------------------------------

namespace coldae
{

/* Define callback function types */
using fsub_t  = void (*)(double x, double const z[], double const y[], double f[]);
using dfsub_t = void (*)(double x, double const z[], double const y[], double df[]);
using gsub_t  = void (*)(int i,    double const z[],                   double& g);
using dgsub_t = void (*)(int i,    double const z[],                   double dg[]);
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
	customCoefficientsRefine = 3, /* if a former mesh and approximate solution
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
enum class result_t {
	normal = 1,        // for normal return
	singular = 0,      // if the collocation matrix is singular
	outOfMemory = -1,  // if the expected no. of subintervals exceeds storage specifications
	notConverged = -2, // if the nonlinear iteration has not converged
	inputError = -3    // if there is an input data error.
};

/* System specific properties */
struct systemParams {
	int ncomp;                    // number of differential equations (<= 20)
	int ny;                       // number of constraints (<= 20)
	std::vector<int> orders ;     // orders of odes
	double left;                  // left end of interval
	double right;                 // right end of interval
	std::vector<double> bcpoints; // j-th side condition point (boundary point)

	bool isNonLinear;             // if the problem is nonlinear
	regularControl regularity;
	indexControl index;           // index of DAE (ignored if ny=0)

    // Returns the total order of the system
    int getMstar() const{
        int mstar{0};
        for(auto const& o: orders) mstar+=o;
        return mstar;
    }
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
	
	std::vector<int> ltol;      // what component of z the tolerances are for
	std::vector<double> tol;    // the tolerances tol

	int numFixedPoints;         // no. of fixed points in the mesh other than aleft and aright.
	std::vector<double> fixpnt; // the fixed points
};


/* Short notation for pointers. 
   Matrices are stores one-dimensional, column-wise.
   Use restrict for 4% performance increase (g++). */
using dvec = double* const __restrict;
using ivec = int* const __restrict;
using dmat = double* const __restrict;
using imat = int* const __restrict;

using cdvec = double const * const __restrict;
using civec = int const * const __restrict;
using cdmat = double const * const __restrict;
using cimat = int const * const __restrict;



/* Class that holds the solver state. */
class cda{

	// COLOUT 
	double PRECIS;
    int IPRINT;
	
	// COLLOC 
	double RHO[7];
	double COEF[7*7];
	
	// COLORD 
	int K; 
	int NCOMP;
	int NY;
	int NCY;
	int MSTAR;
	int KDY;
	int MMAX;
	int MT[20];
	
	// COLAPR 
	int N, NOLD, NMAX, NZ, NDMZ;
	
	// COLMSH 
	int MSHFLG, MSHNUM, MSHLMT, MSHALT;
	
	// COLSID 	
	double TZETA[40]; // boundary condition points
	double TLEFT;     // left domain end
	double TRIGHT;    // right domain end
	int    IZETA;
	int    IZSAVE;
	
	// COLNLN 
	int NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX;
	
	// COLEST 
	double TOL[40], WGTMSH[40], WGTERR[40], TOLIN[40], ROOT[40];
	int    JTOL[40], LTOL[40];
	int    NTOL; // number of tolerances
	
	// COLBAS 
	double B[7*4];
	double ACOL[28*7];
	double ASAVE[28*4];
	

public:
	result_t COLDAE(systemParams const& params, options const& opts,
                    ivec ispace, dvec fspace,
                    fsub_t fsub, dfsub_t dfsub, gsub_t gsub, dgsub_t dgsub, guess_t guess)
	{
		this->NCOMP = params.ncomp;
		this->NY = params.ny;
		this->TLEFT = params.left;
		this->TRIGHT = params.right;

		if (opts.ltol.size() != opts.tol.size()) {
			return result_t::inputError;
		}

		std::copy(opts.ltol.begin(), opts.ltol.end(), this->LTOL);
		std::copy(opts.tol.begin(), opts.tol.end(), this->TOL);
	

	
		//*********************************************************************
		//
		//     the actual subroutine coldae serves as an interface with
		//     the package of subroutines referred to collectively as
		//     coldae. the subroutine serves to test some of the input
		//     parameters, rename some of the parameters (to make understanding
		//     of the coding easier), to do some initialization,
		//     and to break the work areas fspace and ispace up into the
		//     arrays needed by the program.
		//
		//**********************************************************************


		if (opts.printLevel != printMode::none) 
			fmt::print("VERSION *1* OF COLDAE\n");

		//  specify machine dependent constant  precis = 100 * machine unit roundoff
		PRECIS = std::numeric_limits<double>::epsilon();
		PRECIS = PRECIS * 100.0;

		//  in case incorrect input data is detected, the program returns
		//  immediately with iflag=-3.
		NCY = NCOMP + NY;
		if (NCOMP < 0 || NCOMP > 20){
		    fmt::print(fg(fmt::color::red), "Violated (0 <= NCOMP <= 20).\n");
			return result_t::inputError;
        }
		if (NY < 0 || NY > 20){
            fmt::print(fg(fmt::color::red), "Violated (0 <= NY <= 20).\n");
            return result_t::inputError;
        }
		if (NCY < 1 || NCY > 40) {
            fmt::print(fg(fmt::color::red), "Violated (1 <= NCOMP + NY <= 40).\n");
            return result_t::inputError;
        }
		for (int i = 0; i < NCOMP; ++i)
			if (params.orders[i] < 1 || params.orders[i] > 4) {
                fmt::print(fg(fmt::color::red),
                           "Violated (1 <= orders[{}] <= 4).\n", i);
                return result_t::inputError;
            }

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
		ICARE = static_cast<int>(params.regularity);
		NTOL = static_cast<int>(opts.tol.size());
		int NDIMF = opts.fdim;
		int NDIMI = opts.idim;
		int NFXPNT = opts.numFixedPoints;
		IPRINT = static_cast<int>(opts.printLevel);
		INDEX = static_cast<int>(params.index);
		if (NY == 0) INDEX = 0;
		MSTAR = 0;
		MMAX = 0;

		for (int i = 1; i <= NCOMP; ++i) {
			MMAX = std::max(MMAX, params.orders[i-1]);
			MSTAR = MSTAR + params.orders[i - 1];
			MT[i-1] = params.orders[i - 1];
		}
		if (K == 0)   K = std::max(MMAX + 1, 5 - MMAX);
		for (int i = 1; i <= MSTAR; ++i)
			TZETA[i-1] = params.bcpoints[i - 1];
		for (int i = 1; i <= NTOL; ++i) {
			LTOL[i-1] = LTOL[i-1];
			TOLIN[i-1] = TOL[i-1];
		}
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
			fmt::print("COMPONENTS OF Z REQUIRING TOLERANCES: {}\n", 
				std::vector<int>(LTOL, LTOL + NTOL));
			fmt::print("CORRESPONDING ERROR TOLERANCES: {}\n",
				std::vector<double>(TOL, TOL + NTOL));

			if (IGUESS >= 2) {
				fmt::print("INITIAL MESH(ES) AND Z, DMZ PROVIDED BY USER\n");
			}
			if (IREAD == 2) {
				fmt::print("NO ADAPTIVE MESH SELECTION\n");
			}
		}

		//  check for correctness of data
		if (K < 0 || K > 7){
            fmt::print(fg(fmt::color::red), "Violated (0 <= K <= 7).\n");
            return result_t::inputError;
        }
		if (N < 0){
            fmt::print(fg(fmt::color::red), "Violated (0 <= N).\n");
            return result_t::inputError;
        }
		if (NTOL < 0 || NTOL > MSTAR){
            fmt::print(fg(fmt::color::red), "Violated (0 <= NTOL <= MSTAR).\n");
            return result_t::inputError;
        }
		if (NFXPNT < 0) {
            fmt::print(fg(fmt::color::red), "Violated (0 <= NFXPNT).\n");
            return result_t::inputError;
        }
		if (MSTAR < 0 || MSTAR > 40){
            fmt::print(fg(fmt::color::red), "Violated (0 <= MSTAR <= 40).\n");
            return result_t::inputError;
        }

		int IP = 1;
		for (int i = 1; i <= MSTAR; ++i) {
			if (std::abs(params.bcpoints[i-1] - TLEFT) < PRECIS || std::abs(params.bcpoints[i - 1] - TRIGHT) < PRECIS)
				continue;

			while (true) {
				if (IP > NFXPNT)
					return result_t::inputError;
				if (params.bcpoints[i - 1] - PRECIS < opts.fixpnt[IP])
					break;
				IP = IP + 1;
			}

			if (params.bcpoints[i - 1] + PRECIS < opts.fixpnt[IP])
				return result_t::inputError;
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
			if (params.bcpoints[IB - 1] >= TRIGHT)
				NREC = i;
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
			return result_t::inputError;
		if (NMAX < NFXPNT + 1)   
			return result_t::inputError;
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
				NOLD = ispace[1-1];
			NZ = MSTAR * (NOLD + 1);
			NDMZ = KDY * NOLD;
			int NP1 = N + 1;
			if (IGUESS == 4)
				NP1 = NP1 + NOLD + 1;
			for (int i = 1; i <= NZ; ++i)
				fspace[LZ + i - 1-1] = fspace[NP1 + i-1];
			int IDMZ = NP1 + NZ;
			for (int i = 1; i <= NDMZ; ++i)
				fspace[LDMZ + i - 1-1] = fspace[IDMZ + i-1];
			NP1 = NOLD + 1;
			if (IGUESS == 4) {
				for (int i = 1; i <= NP1; ++i)
					fspace[LXIOLD + i - 1-1] = fspace[N + 1 + i-1];
			}
			else {
				for (int i = 1; i <= NP1; ++i)
					fspace[LXIOLD + i - 1-1] = fspace[LXI + i - 1-1];
			}
		}

		//  initialize collocation points, constants, mesh.
		CONSTS();
		int NYCB = (NY == 0 ? 1 : NY);

		int meshmode = 3 + IREAD;
		NEWMSH(meshmode, fspace+(LXI-1), fspace+(LXIOLD-1),
			nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
			NFXPNT, opts.fixpnt.data(), nullptr, dfsub, nullptr, nullptr, NYCB);

		//  determine first approximation, if the problem is nonlinear.
		if (IGUESS < 2) {
			for (int i = 1; i <= N + 1; ++i)
				fspace[i + LXIOLD - 1-1] = fspace[i + LXI - 1-1];
			NOLD = N;
			if (NONLIN != 0 && IGUESS != 1) {
				//  system provides first approximation of the solution.
				//  choose z(j) = 0  for j=1,...,mstar.
				for (int i = 1; i <= NZ; ++i)
					fspace[LZ - 1 + i-1] = 0.0;
				for (int i = 1; i <= NDMZ; ++i)
					fspace[LDMZ - 1 + i-1] = 0.0;
			}
		}
		if (IGUESS >= 2)  
			IGUESS = 0;

		result_t iflag;
		CONTRL(fspace+(LXI-1), fspace+(LXIOLD-1), fspace+(LZ-1), fspace+(LDMZ-1), fspace+(LDMV-1),
			fspace+(LRHS-1), fspace+(LDELZ-1), fspace+(LDELDZ-1), fspace+(LDQZ-1),
			fspace+(LDQDMZ-1), fspace+(LG-1), fspace+(LW-1), fspace+(LV-1), fspace+(LFC-1),
			fspace+(LVALST-1), fspace+(LSLOPE-1), fspace+(LSCL-1), fspace+(LDSCL-1),
			fspace+(LACCUM-1), ispace+(LPVTG-1), ispace+(LINTEG-1), ispace+(LPVTW-1),
			NFXPNT, opts.fixpnt.data(), iflag, fsub, dfsub, gsub, dgsub, guess);

		//  prepare output
		ispace[1-1] = N;
		ispace[2-1] = K;
		ispace[3-1] = NCOMP;
		ispace[4-1] = NY;
		ispace[5-1] = MSTAR;
		ispace[6-1] = MMAX;
		ispace[7-1] = NZ + NDMZ + N + 2;
		int K2 = K * K;
		ispace[8-1] = ispace[7-1] + K2 - 1;
		for (int i = 1; i <= NCOMP; ++i)
			ispace[8 + i-1] = params.orders[i - 1];
		for (int i = 1; i <= NZ; ++i)
			fspace[N + 1 + i-1] = fspace[LZ - 1 + i-1];
		int IDMZ = N + 1 + NZ;

		for (int i = 1; i <= NDMZ; ++i)
			fspace[IDMZ + i-1] = fspace[LDMZ - 1 + i-1];
		int IC = IDMZ + NDMZ;
		for (int i = 1; i <= K2; ++i)
			fspace[IC + i-1] = COEF[i-1];

		return iflag;
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
	void APPSLN(double& X, dvec Z, dvec Y, cdvec FSPACE, civec ISPACE)
	{
		double A[28];

		int IS6 = ISPACE[6];
		int IS5 = ISPACE[0] + 2;
		int IS4 = IS5 + ISPACE[4] * (ISPACE[0] + 1);
		int i = 1;

		APPROX(i, X, Z, Y, A, FSPACE+(IS6-1),
			FSPACE, ISPACE[0],
			FSPACE+(IS5-1), FSPACE+(IS4-1), ISPACE[1], ISPACE[2],
			ISPACE[3], ISPACE[5], ISPACE+8, ISPACE[4], 2, nullptr, 1);
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
	//     rnorm  - norm of rhs (right-hand side) for current iteration
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
	void CONTRL(dvec XI, dvec XIOLD, dvec Z, dvec DMZ, dvec DMV, dvec RHS, dvec DELZ, dvec DELDMZ,
                dvec DQZ, dvec DQDMZ, dvec G, dvec W, dvec V, dvec FC, dvec VALSTR, dvec SLOPE, dvec SCALE, dvec DSCALE,
                dvec ACCUM, ivec IPVTG, ivec INTEGS, ivec IPVTW, const int NFXPNT, cdvec FIXPNT, result_t& iflag,
                fsub_t fsub, dfsub_t dfsub, gsub_t gsub, dgsub_t dgsub, guess_t guess)
	{
		double DF[800];
		std::vector<double> FCSP(NCOMP * 60);

        double CBSP[20*20];
		double RNORM, RNOLD;

		// constants for control of nonlinear iteration
		double RELMIN = 1.e-3;
		double RSTART = 1.e-2;
		double LMTFRZ = 4;

		// compute the maximum tolerance
		double CHECK = 0.0;
		for (int i = 1; i <= NTOL; ++i)
			CHECK = std::max(TOLIN[i-1], CHECK);
		int IMESH = 1;
		int ICONV = (NONLIN == 0);
		int ICOR = 0;
		int NOCONV = 0;
		int MSING = 0;
		int ISING = 0;

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
				// initialization for a new mesh
				ITER = 0;
				if (NONLIN <= 0) {

					// the linear case.
					// set up and solve equations
					LSYSLV(MSING, XI, XIOLD, nullptr, nullptr,
						Z, DMZ, G,
						W, V, FC, RHS,
						nullptr, INTEGS, IPVTG, IPVTW,
						RNORM, 0, fsub, dfsub, gsub, dgsub, guess, ISING);

					// check for a singular matrix
					if (ISING != 0) {
						if (IPRINT < 1) {
							fmt::print("SINGULAR PROJECTION MATRIX DUE TO INDEX > 2\n");
						}
						iflag = result_t::singular;
						return;
					}
					if (MSING == 0)
						goto n400;
				n30:
					if (MSING >= 0) {
						if (IPRINT < 1) {
							fmt::print(fg(fmt::color::red), "A LOCAL ELIMINATION MATRIX IS SINGULAR\n");
						}
						goto n460;
					}
					if (IPRINT < 1) {
						fmt::print(fg(fmt::color::red), "THE GLOBAL BVP - MATRIX IS SINGULAR\n");
					}
					iflag = result_t::singular;
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
				LSYSLV(MSING, XI, XIOLD, Z, DMZ,
					DELZ, DELDMZ, G,
					W, V, FC, RHS,
					DQDMZ, INTEGS, IPVTG, IPVTW,
					RNOLD, 1, fsub, dfsub, gsub, dgsub, guess, ISING);

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
				// solve for the next iterate .
				// the value of ifreez determines whether this is a full
				// newton step (=0) or a fixed jacobian iteration (=1).
				if (IPRINT < 0) {
					fmt::print("ITERATION = {}, NORM(RHS) = {}\n", ITER, RNORM);
				}
				RNOLD = RNORM;
				LSYSLV(MSING, XI, XIOLD, Z, DMZ,
					DELZ, DELDMZ, G,
					W, V, FC, RHS, nullptr,
					INTEGS, IPVTG, IPVTW, RNORM,
					3 + IFREEZ, fsub, dfsub, gsub, dgsub, guess, ISING);
			}
		n70:
			{
				// check for a singular matrix
				if (MSING != 0)
					goto n30;
				if (ISING != 0) {
					if (IPRINT < 1) {
						fmt::print(fg(fmt::color::red), "SINGULAR PROJECTION MATRIX DUE TO INDEX > 2\n");
					}
					iflag = result_t::singular;
					return;
				}
				if (IFREEZ != 1) {
					// a full newton step
					ITER = ITER + 1;
					IFRZ = 0;
				}

				// update   z and dmz , compute new  rhs  and its norm
				for (int i = 1; i <= NZ; ++i)
					Z[i-1] +=  DELZ[i-1];
				for (int i = 1; i <= NDMZ; ++i)
					DMZ[i-1] +=  DELDMZ[i-1];

				LSYSLV(MSING, XI, XIOLD, Z, DMZ,
					DELZ, DELDMZ, G,
					W, V, FC, RHS, nullptr,
					INTEGS, IPVTG, IPVTW, RNORM, 2,
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
				// verify that the linear convergence with fixed jacobian
				// is fast enough.
				IFRZ = IFRZ + 1;
				if (IFRZ >= LMTFRZ)       IFREEZ = 0;
				if (RNOLD < 4.0 * RNORM)  IFREEZ = 0;

				// check convergence (iconv = 1).
				for (int IT = 1; IT <= NTOL; ++IT) {
					int INZ = LTOL[IT-1];
					for (int IZ = INZ; IZ <= NZ; IZ += MSTAR) {
						if (std::abs(DELZ[IZ-1]) > TOLIN[IT-1] * (std::abs(Z[IZ-1]) + 1.0))
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
				// convergence of fixed jacobian iteration failed.
				if (IPRINT < 0)
					fmt::print("ITERATION = {}, NORM(RHS) = {}\nSWITCH TO DAMPED NEWTON ITERATION\n", ITER, RNORM);
				ICONV = 0;
				if (ICARE != 1)
					RELAX = RSTART;
				for (int i = 1; i <= NZ; ++i)
					Z[i-1] -= DELZ[i-1];
				for (int i = 1; i <= NDMZ; ++i)
					DMZ[i-1] -=  DELDMZ[i-1];


				// update old mesh
				for (int i = 1; i <= N + 1; ++i)
					XIOLD[i-1] = XI[i-1];
				NOLD = N;
				ITER = 0;
			}
		n160:
			{
				// no previous convergence has been obtained. proceed
				// with the damped newton method.
				// evaluate rhs and find the first newton correction.
				if (IPRINT < 0)
					fmt::print("FULL DAMPED NEWTON ITERATION\n");
				LSYSLV(MSING, XI, XIOLD, Z, DMZ,
					DELZ, DELDMZ, G,
					W, V, FC, RHS,
					DQDMZ, INTEGS, IPVTG, IPVTW,
					RNOLD, 1, fsub, dfsub, gsub, dgsub, guess, ISING);

				// check for a singular matrix
				if (MSING != 0)
					goto n30;
				if (ISING != 0) {
					if (IPRINT < 1)  fmt::print(fg(fmt::color::red), "SINGULAR PROJECTION MATRIX DUE TO INDEX > 2\n");
					iflag = result_t::singular;
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
				// main iteration loop
				RNOLD = RNORM;
				if (ITER >= LIMIT)
					goto n430;

				// update scaling
				SKALE(Z, DMZ, XI, SCALE, DSCALE);

				// compute norm of newton correction with new scaling
				double ANSCL = 0.0;
				for (int i = 1; i <= NZ; ++i)
					ANSCL = ANSCL + pow(DELZ[i-1] * SCALE[i-1], 2);

				for (int i = 1; i <= NDMZ; ++i)
					ANSCL = ANSCL + pow(DELDMZ[i-1] * DSCALE[i-1], 2);

				ANSCL = sqrt(ANSCL / double(NZ + NDMZ));

				// find a newton direction
				LSYSLV(MSING, XI, XIOLD, Z, DMZ,
					DELZ, DELDMZ, G,
					W, V, FC, RHS, nullptr,
					INTEGS, IPVTG, IPVTW, RNORM, 3,
					fsub, dfsub, gsub, dgsub, guess, ISING);

				// check for a singular matrix
				if (MSING != 0)
					goto n30;

				if (ISING != 0) {
					if (IPRINT < 1)  fmt::print(fg(fmt::color::red), "SINGULAR PROJECTION MATRIX DUE TO INDEX > 2\n");
					iflag = result_t::singular;
					return;
				}

				// predict relaxation factor for newton step.
				if (ICARE != 1) {
					double ANDIF = 0.0;
					for (int i = 1; i <= NZ; ++i)
						ANDIF = ANDIF + pow((DQZ[i-1] - DELZ[i-1]) * SCALE[i-1], 2);

					for (int i = 1; i <= NDMZ; ++i)
						ANDIF = ANDIF + pow((DQDMZ[i-1] - DELDMZ[i-1]) * DSCALE[i-1], 2);

					ANDIF = sqrt(ANDIF / double(NZ + NDMZ) + PRECIS);
					RELAX *= ANSCL / ANDIF;
					if (RELAX > 1.0)  RELAX = 1.0;
					RLXOLD = RELAX;
					IPRED = 1;
				}
			}
		n220:
			{
				ITER++;

				// determine a new  z and dmz  and find new  rhs  and its norm
				for (int i = 1; i <= NZ; ++i)
					Z[i-1] += RELAX * DELZ[i-1];

				for (int i = 1; i <= NDMZ; ++i)
					DMZ[i-1] +=  RELAX * DELDMZ[i-1];
			}
		n250:
			{
				LSYSLV(MSING, XI, XIOLD, Z, DMZ,
					DQZ, DQDMZ, G,
					W, V, FC, RHS, nullptr,
					INTEGS, IPVTG, IPVTW, RNORM, 2,
					fsub, dfsub, gsub, dgsub, guess, ISING);

				// compute a fixed jacobian iterate (used to control relax)
				LSYSLV(MSING, XI, XIOLD, Z, DMZ,
					DQZ, DQDMZ, G,
					W, V, FC, RHS, nullptr,
					INTEGS, IPVTG, IPVTW, RNORM, 4,
					fsub, dfsub, gsub, dgsub, guess, ISING);

				// find scaled norms of various terms used to correct relax
				ANORM = 0.0;
				ANFIX = 0.0;
				for (int i = 1; i <= NZ; ++i) {
					ANORM += pow(DELZ[i-1] * SCALE[i-1], 2);
					ANFIX += pow(DQZ[i-1] * SCALE[i-1], 2);
				}
				for (int i = 1; i <= NDMZ; ++i) {
					ANORM += pow(DELDMZ[i-1] * DSCALE[i-1], 2);
					ANFIX += pow(DQDMZ[i-1] * DSCALE[i-1], 2);
				}
				ANORM = sqrt(ANORM / double(NZ + NDMZ));
				ANFIX = sqrt(ANFIX / double(NZ + NDMZ));
				if (ICOR == 1) {
					if (IPRINT < 0)
						fmt::print("RELAXATION FACTOR CORRECTED TO RELAX = {}\n"
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
						if (std::abs(FACTOR - 1.0) < .10 * FACTOR)
							goto n170;
						if (FACTOR < 0.50)
							FACTOR = 0.5;
						RELAX /= FACTOR;
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
							Z[i-1] += FACT * DELZ[i-1];

						for (int i = 1; i <= NDMZ; ++i) {
							DMZ[i-1] += FACT * DELDMZ[i-1];
						}
						RLXOLD = RELAX;
						goto n250;
					}
				n350:
					{
						// check convergence (iconv = 0).
						for (int IT = 1; IT <= NTOL; ++IT) {
							int INZ = LTOL[IT-1];
							for (int IZ = INZ; IZ <= NZ; IZ += MSTAR) {
								if (std::abs(DQZ[IZ-1]) > TOLIN[IT-1] * (std::abs(Z[IZ-1]) + 1.0))
									goto n170;
							}
						}

						// convergence obtained
						if (IPRINT < 1)
							fmt::print("CONVERGENCE AFTER {} ITERATIONS\n", ITER);

						// since convergence obtained, update  z and dmz  with term
						// from the fixed jacobian iteration.
						for (int i = 1; i <= NZ; ++i)
							Z[i-1] += DQZ[i-1];

						for (int i = 1; i <= NDMZ; ++i)
							DMZ[i-1] += DQDMZ[i-1];
					}
				n390:
					{
						if ((ANFIX < PRECIS || RNORM < PRECIS) && IPRINT < 1)
							fmt::print("CONVERGENCE AFTER {} ITERATIONS\n", ITER);
						ICONV = 1;
						if (ICARE == -1)  ICARE = 0;
					}
				n400:
					{
						// if full output has been requested, print values of the
						// solution components   z  at the meshpoints and  y  at
						// collocation points.
						if (IPRINT >= 0)
							goto n420;
						for (int j = 1; j <= MSTAR; ++j) {
							fmt::print("MESH VALUES FOR Z[{}-1]: ", j);
							for (int LJ = j; LJ <= NZ; LJ += MSTAR)
								fmt::print("{:.4}, ", Z[LJ-1]);
							fmt::print("\n");
						}
						for (int j = 1; j <= NY; ++j) {
							fmt::print("VALUES AT 1st COLLOCATION POINTS FOR Y({}): ", j);
							for (int LJ = j + NCOMP; LJ <= NDMZ; LJ += KDY)
								fmt::print("{:.4}, ", DMZ[LJ-1]);
							fmt::print("\n");
						}
					}
				n420:
					{
						// check for error tolerance satisfaction
						int IFIN = 1;
						if (IMESH == 2)
							ERRCHK(XI, Z, DMZ,VALSTR, IFIN);
						if (IMESH == 1 || (IFIN == 0 && ICARE != 2))
							goto n460;
						iflag = result_t::normal;
						return;
					}
				n430:
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

				iflag = result_t::notConverged;
				NOCONV = NOCONV + 1;
				if (ICARE == 2 && NOCONV > 1)
					return;
				if (ICARE == 0)
					ICARE = -1;
			}
		n460:
			{
				// update old mesh
				for (int i = 1; i <= N + 1; ++i)
					XIOLD[i-1] = XI[i-1];
				NOLD = N;

				// pick a new mesh
				// check safeguards for mesh refinement
				IMESH = 1;
				if (ICONV == 0 || MSHNUM >= MSHLMT || MSHALT >= MSHLMT)
					IMESH = 2;
				if (MSHALT >= MSHLMT && MSHNUM < MSHLMT)
					MSHALT = 1;

                int NYCB = (NY == 0 ? 1 : NY);


				NEWMSH(IMESH, XI, XIOLD, Z, DMZ,
					DMV, VALSTR,
					SLOPE, ACCUM, NFXPNT,
					FIXPNT, DF, dfsub,
                    FCSP.data(), CBSP, NYCB);

				// exit if expected n is too large (but may try n=nmax once)
				if (N > NMAX) {
					N = N / 2;
					iflag = result_t::outOfMemory;
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
	void SKALE(dmat Z, dmat DMZ, cdvec XI, dmat SCALE, dmat DSCALE)
	{
		//Z(MSTAR, 1);
		//SCALE(MSTAR, 1); 
		//DMZ(KDY, N); 
		//DSCALE(KDY, N); 
	
		double BASM[5];

		BASM[0] = 1.0;
		for (int j = 1; j <= N; ++j) {
			int IZ = 1;
			double H = XI[j] - XI[j-1];
			for (int l = 1; l <= MMAX; ++l)
				BASM[l] = BASM[l-1] * H / double(l);

			for (int ICOMP = 1; ICOMP <= NCOMP; ++ICOMP) {
				double SCAL = (std::abs(Z[IZ-1+(j-1)* MSTAR]) + std::abs(Z[IZ - 1 + (j) * MSTAR])) * .5 + 1.0;
				int MJ = MT[ICOMP-1];
				for (int l = 1; l <= MJ; ++l) {
					SCALE[IZ - 1 + (j - 1) * MSTAR] = BASM[l-1] / SCAL;
					IZ = IZ + 1;
				}
				SCAL = BASM[MJ] / SCAL;
				for (int IDMZ = ICOMP; IDMZ <= KDY; IDMZ += NCY)
					DSCALE[IDMZ-1+ (j-1)*KDY] = SCAL;

			}
			for (int ICOMP = 1 + NCOMP; ICOMP <= NCY; ++ICOMP) {
				double SCAL = 1.0 / (std::abs(DMZ[ICOMP-1+(j-1)*KDY]) + 1.0);
				for (int IDMZ = ICOMP; IDMZ <= KDY; IDMZ += NCY) {
					DSCALE[IDMZ - 1 + (j - 1) * KDY] = SCAL;
				}
			}
		}
		for (int IZ = 1; IZ <= MSTAR; ++IZ)
			SCALE[IZ - 1 + (N) * MSTAR] = SCALE[IZ - 1 + (N - 1) * MSTAR];

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
	//                     where j=JTOL[l-1]
	//            slphmx - maximum of slope(i)*(xiold(i+1)-xiold(i)) for
	//                     i = 1 ,..., nold.
	//            accum  - accum(i) is the integral of  slope  from  aleft
	//                     to  xiold(i).
	//            valstr - is assigned values needed in  errchk  for the
	//                     error estimate.
	//            fc     - you know
	//**********************************************************************
	void NEWMSH(int& MODE, dvec XI, cdvec XIOLD, 
		cdvec Z, cdvec DMZ, dvec DMV,
		dvec VALSTR, dvec SLOPE, dvec ACCUM, 
		const int NFXPNT, cdvec FIXPNT,
		dmat DF, dfsub_t dfsub, dmat FC,
		dmat CB, const int NYCB)
	{
		//XIOLD(NOLD + 1);
		//FC(NCOMP, 60);
		//DF(NCY, 1);
		//CB(NYCB, NYCB);

	
		double D1[40], D2[40], ZVAL[40], YVAL[40], A[28], BCOL[40], U[400], V[400];
		int IPVTCB[40];

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
					fmt::print("THE FORMER MESH (OF {} SUBINTERVALS): {}\n",
						NOLD, std::vector<double>(XIOLD, XIOLD + NOLD));
				if (IGUESS == 3) {
					//  if iread ( ipar(8) ) .ge. 1 and iguess ( ipar(9) ) .eq. 3
					//  then the first mesh is every second point of the mesh in  xiold .
					N = NOLD / 2;
					for (int j = 1, i = 0; j <= NOLD; j += 2) {
						XI[i] = XIOLD[j-1];
						i = i + 1;
					}
				}
			}
			XI[0] = TLEFT;
			XI[N] = TRIGHT;
			break;
		}
		case 3: {
			//  mode=3   generate a (piecewise) uniform mesh. if there are
			//  fixed points then ensure that the n being used is large enough
			if (N < NFXP1)
				N = NFXP1;
			int NP1 = N + 1;
			XI[0] = TLEFT;
			int ILEFT = 1;
			double XLEFT = TLEFT;

			//  loop over the subregions between fixed points.
			for (int j = 1; j <= NFXP1; ++j) {
				double XRIGHT = TRIGHT;
				int IRIGHT = NP1;
				if (j != NFXP1) {
					XRIGHT = FIXPNT[j-1];

					// determine where the j-th fixed point should fall in the
					// new mesh - this is xi(iright) and the (j-1)st fixed
					// point is in xi(ileft)
					int NMIN = int((XRIGHT - TLEFT) / (TRIGHT - TLEFT) * double(N) + 1.5);
					if (NMIN > N - NFXPNT + j)
						NMIN = N - NFXPNT + j;
					IRIGHT = std::max(ILEFT + 1, NMIN);
				}
				XI[IRIGHT-1] = XRIGHT;

				// generate equally spaced points between the j-1st and the
				// j-th fixed points.
				int NREGN = IRIGHT - ILEFT - 1;
				if (NREGN != 0) {
					double DX = (XRIGHT - XLEFT) / double(NREGN + 1);
					for (int i = 1; i <= NREGN; ++i)
						XI[ILEFT + i-1] = XLEFT + double(i) * DX;
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
					double HD6 = (XIOLD[i ] - XIOLD[i-1]) / 6.0;
					double X = XIOLD[i-1] + HD6;
					APPROX(i, X, VALSTR+(KSTORE-1), nullptr,
						ASAVE, nullptr, XIOLD,
						NOLD, Z, DMZ, K, NCOMP, NY, MMAX,
						MT, MSTAR, 4, nullptr, 0);
					X = X + 4.0 * HD6;
					KSTORE = KSTORE + 3 * MSTAR;
					APPROX(i, X, VALSTR+(KSTORE-1), nullptr,
						ASAVE+(3*28), nullptr, XIOLD,
						NOLD, Z, DMZ, K, NCOMP, NY, MMAX, 
						MT, MSTAR, 4, nullptr, 0);
					KSTORE += MSTAR;
				}
			}
			//  save in  valstr  the values of the old solution
			//  at the relative positions 1/6, 2/6, 4/6 and 5/6 in
			//  each subinterval.
			else {
				int KSTORE = 1;
				for (int i = 1; i <= N; ++i) {
					double X = XI[i-1];
					double HD6 = (XI[i] - XI[i - 1]) / 6.0;
					for (int j = 1; j <= 4; ++j) {
						X = X + HD6;
						if (j == 3)
							X = X + HD6;
						APPROX(i, X, VALSTR+(KSTORE-1), nullptr, 
							ASAVE+(j-1)*28, nullptr,
							XIOLD, NOLD, Z, DMZ,
							K, NCOMP, NY, MMAX, MT, MSTAR, 4, nullptr, 0);
						KSTORE = KSTORE + MSTAR;
					}
				}
			}
			MSHFLG = 0;
			MSHNUM = 1;
			MODE = 2;

			//  generate the halved mesh.
			for (int i = 1, j = 2; i <= N; ++i) {
				XI[j-1] = (XIOLD[i - 1] + XIOLD[i]) / 2.0;
				XI[j] = XIOLD[i];
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
							DMV[IDMZ-1] = DMZ[IDMZ-1];
							IDMZ = IDMZ + 1;
						}
					}
				}
				if (INDEX != 1 && NY > 0) {
					IDMZ = 1;
					for (int i = 1; i <= NOLD; ++i) {
						double XI1 = XIOLD[i];
						APPROX(i, XI1, ZVAL, YVAL, A,
							COEF, XIOLD, NOLD, Z, DMZ,
							K, NCOMP, NY, MMAX, MT, MSTAR, 3, nullptr, 1);
						dfsub(XI1, ZVAL, YVAL, DF);

						// if index=2, form projection matrices directly
						// otherwise use svd to define appropriate projection
						if (INDEX == 0) {
							PRJSVD(FC, DF, CB,
								U, V, IPVTCB, ISING, 2);
						}
						else {
							// form cb
							for (int j = 1; j <= NY; ++j) {
								for (int J1 = 1; J1 <= NY; ++J1) {
									double FACT = 0.0;
									for (int l = 1, ML = 0; l <= NCOMP; ++l) {
										ML += MT[l-1];
										FACT += DF[j + NCOMP-1 + (ML-1)*NCY] * DF[l-1+ (MSTAR + J1-1) * NCY];
									}
									CB[j-1+(J1-1)*NYCB] = FACT;
								}
							}

							// decompose cb
							ISING = dgefa(CB, NY, NY, IPVTCB);
							if (ISING != 0)
								return;

							// form columns of fc
							int ML = 0;
							for (int l = 1; l <= NCOMP; ++l) {
								ML += MT[l-1];
								for (int J1 = 1; J1 <= NY; ++J1)
									BCOL[J1-1] = DF[J1 + NCOMP-1+ (ML-1) * NCY];

								dgesl(CB, NY, NY, IPVTCB, BCOL, 0);

								for (int J1 = 1; J1 <= NCOMP; ++J1) {
									double FACT = 0.0;
									for (int j = 1; j <= NY; ++j)
										FACT += DF[J1-1+(j + MSTAR-1) * NCY] * BCOL[j-1];
									FC[J1-1+ (l-1)*NCOMP] = FACT;
								}
							}
						}

						// finally, replace fc with the true projection SR = i - fc
						for (int j = 1; j <= NCOMP; ++j) {
							for (int l = 1; l <= NCOMP; ++l) {
								FC[j-1+ (l-1)* NCOMP] *= -1;
								if (j == l)
									FC[j - 1 + (l - 1) * NCOMP] += 1.0;
							}
						}

						// project DMZ for the k collocation points, store in DMV
						for (int KK = 1; KK <= K; ++KK) {
							for (int j = 1; j <= NCOMP; ++j) {
								double FACT = 0.0;
								for (int l = 1; l <= NCOMP; ++l)
									FACT += FC[j - 1 + (l - 1) * NCOMP] * DMZ[IDMZ + l - 1-1];
								DMV[IDMZ + j - 1-1] = FACT;
							}
							IDMZ += NCY;
						}
					}
				}

				//  the first interval has to be treated separately from the
				//  other intervals (generally the solution on the (i-1)st and ith
				//  intervals will be used to approximate the needed derivative, but
				//  here the 1st and second intervals are used.)
				double HIOLD = XIOLD[1] - XIOLD[0];
				HORDER(1, D1, HIOLD, DMV);
				IDMZ = IDMZ + (NCOMP + NY) * K;
				HIOLD = XIOLD[2] - XIOLD[1];
				HORDER(2, D2, HIOLD, DMV);
				ACCUM[0] = 0.0;
				SLOPE[0] = 0.0;
				double ONEOVH = 2.0 / (XIOLD[2] - XIOLD[0]);
				for (int j = 1; j <= NTOL; ++j) {
					int JJ = JTOL[j-1];
					int JZ = LTOL[j-1];
					SLOPE[0] = std::max(SLOPE[0],
						pow(std::abs(D2[JJ-1] - D1[JJ-1]) * WGTMSH[j-1] * ONEOVH / (1.0 + std::abs(Z[JZ-1])),
							ROOT[j-1]));
				}
				double SLPHMX = SLOPE[0] * (XIOLD[1] - XIOLD[0]);
				ACCUM[1] = SLPHMX;
				int IFLIP = 1;
			
				//  go through the remaining intervals generating  slope
				//  and  accum .
				for (int i = 2; i <= NOLD; ++i) {
					HIOLD = XIOLD[i] - XIOLD[i-1];
					if (IFLIP == -1)
						HORDER(i, D1, HIOLD, DMV);
					if (IFLIP == 1)
						HORDER(i, D2, HIOLD, DMV);
					ONEOVH = 2.0 / (XIOLD[i] - XIOLD[i - 1-1]);
					SLOPE[i-1] = 0.0;


					// evaluate function to be equidistributed
					for (int j = 1; j <= NTOL; ++j) {
						int JJ = JTOL[j-1];
						int JZ = LTOL[j-1] + (i - 1) * MSTAR;
						auto temp = std::abs(D2[JJ-1] - D1[JJ-1]) * WGTMSH[j-1] * ONEOVH / (1.0 + std::abs(Z[JZ-1]));
						SLOPE[i-1] = std::max(SLOPE[i-1], pow(temp, ROOT[j-1]));
					}

					// accumulate approximate integral of function to be equidistributed
					double TEMP = SLOPE[i-1] * (XIOLD[i ] - XIOLD[i - 1]);
					SLPHMX = std::max(SLPHMX, TEMP);
					ACCUM[i ] = ACCUM[i - 1] + TEMP;
					IFLIP = -IFLIP;
				}

			
				double AVRG = ACCUM[NOLD] / double(NOLD);
				double DEGEQU = AVRG / std::max(SLPHMX, PRECIS);

				//  naccum=expected n to achieve .1x user requested tolerances
				int NACCUM = int(ACCUM[NOLD] + 1.0);
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
				XI[0] = TLEFT;
				XI[N] = TRIGHT;
				for (int i = 1; i <= NFXP1; ++i) {
					double ACCR; int LNEW, NREGN;
					if (i != NFXP1) {
						for (int j = LOLD; j <= NOLDP1; ++j) {
							LNEW = j;
							if (FIXPNT[i - 1] <= XIOLD[j - 1])
								break;
						}
						ACCR = ACCUM[LNEW-1] + (FIXPNT[i - 1] - XIOLD[LNEW - 1]) * SLOPE[LNEW - 1 - 1];
						NREGN = int((ACCR - ACCL) / ACCUM[NOLDP1 - 1] * double(N) - .5);
						NREGN = std::min(NREGN, N - IN - NFXP1 + i);
						XI[IN + NREGN] = FIXPNT[i - 1];
					}
					else {
						ACCR = ACCUM[NOLDP1 - 1];
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
								if (TEMP <= ACCUM[l - 1])
									break;
							}
							LOLD = LCARRY;
							XI[IN - 1] = XIOLD[LOLD - 1 - 1] + (TEMP - ACCUM[LOLD - 1 - 1]) / SLOPE[LOLD - 1 - 1];
						}
					}
					IN++;
					ACCL = ACCR;
					LOLD = LNEW;
				}
				MODE = 1;
				break;
			}
		}
        default: throw std::invalid_argument("NEWMSH: Invalid case.\n");
		} // end of switch


		if (IPRINT < 1) {
			//assert(XI.getSize() == N + 1);
			fmt::print("THE NEW MESH (OF {} SUBINTERVALS): ", N);
			for (int i = 1;i<=N+1;++i)
				fmt::print("{:.2}, ", XI[i - 1]);
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
	//              refer (viz, if LTOL[i-1] refers to a derivative of u(j),
	//              then JTOL[i-1]=j)
	//     root   - reciprocals of expected rates of convergence of components
    //              of z(j) for which tolerances are specified
	//     rho    - the k collocation points on (0,1)
	//     coef   -
	//     acol  -  the runge-kutta coefficients values at collocation
	//              points
	//
	//**********************************************************************
	void CONSTS()
	{
		//COEF(K, K);
	
		double CNSTS1[28] = {0.25e0, 0.625e-1, 7.2169e-2, 1.8342e-2,
			  1.9065e-2, 5.8190e-2, 5.4658e-3, 5.3370e-3, 1.8890e-2,
			  2.7792e-2, 1.6095e-3, 1.4964e-3, 7.5938e-3, 5.7573e-3,
			  1.8342e-2, 4.673e-3,  4.150e-4,  1.919e-3,  1.468e-3,
			  6.371e-3,  4.610e-3,  1.342e-4,  1.138e-4,  4.889e-4,
			  4.177e-4,  1.374e-3,  1.654e-3,  2.863e-3 };
		double CNSTS2[28] = {1.25e-1, 2.604e-3,  8.019e-3, 2.170e-5,
			 7.453e-5,  5.208e-4,  9.689e-8,  3.689e-7,  3.100e-6,
			 2.451e-5,  2.691e-10, 1.120e-9,  1.076e-8,  9.405e-8,
			 1.033e-6,  5.097e-13, 2.290e-12, 2.446e-11, 2.331e-10,
			 2.936e-9,  3.593e-8,  7.001e-16, 3.363e-15, 3.921e-14,
			 4.028e-13, 5.646e-12, 7.531e-11, 1.129e-9 };

		// assign weights for error estimate
		int KOFF = K * (K + 1) / 2;
		int IZ = 1;
		for (int j = 1; j <= NCOMP; ++j) {
			int MJ = MT[j-1];
			for (int l = 1; l <= MJ; ++l) {
				WGTERR[IZ-1] = CNSTS1[KOFF - MJ + l-1];
				IZ = IZ + 1;
			}
		}

		// assign array values for mesh selection: wgtmsh, jtol, and root
		int JCOMP = 1;
		int MTOT = MT[1-1];
		for (int i = 1; i <= NTOL; ++i) {
			int LTOLI = LTOL[i-1];
			while (true) {
				if (LTOLI <= MTOT)
					break;
				JCOMP = JCOMP + 1;
				MTOT = MTOT + MT[JCOMP-1];
			}
			JTOL[i-1] = JCOMP;
			WGTMSH[i-1] = 10 * CNSTS2[KOFF + LTOLI - MTOT-1] / TOLIN[i-1];
			ROOT[i-1] = 1.0 / double(K + MTOT - LTOLI + 1);
		}

		// specify collocation points
		switch (K) {
		case 1: RHO[1-1] = 0.0;
			break;

		case 2:
			RHO[2-1] = .577350269189625764510;
			RHO[1-1] = -RHO[2-1];
			break;

		case 3:
			RHO[3-1] = .774596669241483377040;
			RHO[2-1] = .00;
			RHO[1-1] = -RHO[3-1];
			break;

		case 4:
			RHO[4-1] = .861136311594052575230;
			RHO[3-1] = .339981043584856264800;
			RHO[2-1] = -RHO[3-1];
			RHO[1-1] = -RHO[4-1];
			break;

		case 5:
			RHO[5-1] = .906179845938663992800;
			RHO[4-1] = .538469310105683091040;
			RHO[3-1] = .00;
			RHO[2-1] = -RHO[4-1];
			RHO[1-1] = -RHO[5-1];
			break;

		case 6:
			RHO[6-1] = .932469514203152027810;
			RHO[5-1] = .661209386466264513660;
			RHO[4-1] = .238619186083196908630;
			RHO[3-1] = -RHO[4-1];
			RHO[2-1] = -RHO[5-1];
			RHO[1-1] = -RHO[6-1];
			break;

		case 7:
			RHO[7-1] = .9491079912342758524520;
			RHO[6-1] = .741531185599394439860;
			RHO[5-1] = .405845151377397166900;
			RHO[4-1] = 0.0;
			RHO[3-1] = -RHO[5-1];
			RHO[2-1] = -RHO[6-1];
			RHO[1-1] = -RHO[7-1];
			break;
		}

		// map (-1,1) to (0,1) by  t = .5 * (1. + x)
		for (int j = 1; j<=K; ++j)
			RHO[j-1] = .50 * (1.0 + RHO[j-1]);


		// now find runge-kutta coefficients b, acol and asave
		// the values of asave are to be used in  newmsh  and errchk .
		for (int j = 1; j <= K; ++j) {
			for (int i = 1; i <= K; ++i)
				COEF[i-1+(j-1)*K] = 0.0;
			COEF[j - 1 + (j - 1) * K] = 1.0;
			VMONDE(COEF+(j - 1) * K, K);
		}
		RKBAS(1.0, COEF, K, MMAX, B, nullptr, 0);
		for (int i = 1; i <= K; ++i)
			RKBAS(RHO[i-1], COEF, K, MMAX, ACOL+(i-1)*28, nullptr, 0);

		RKBAS(1.0 / 6.0, COEF, K, MMAX, ASAVE + (0) * 28, nullptr, 0);
		RKBAS(1.0 / 3.0, COEF, K, MMAX, ASAVE + (1) * 28, nullptr, 0);
		RKBAS(2.0 / 3.0, COEF, K, MMAX, ASAVE + (2) * 28, nullptr, 0);
		RKBAS(5.0 / 6.0, COEF, K, MMAX, ASAVE + (3) * 28, nullptr, 0);

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
	void ERRCHK(cdvec XI, cdvec Z, cdvec DMZ, dvec VALSTR, int& IFIN)
	{
		MSHFLG = 1;

		double ERR[40], ERREST[40];


		//  error estimates are to be generated and tested
		//  to see if the tolerance requirements are satisfied.
		IFIN = 1;
		for (int j = 0; j < MSTAR; ++j)
			ERREST[j] = 0.0;
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
			double X = XI[i-1] + (XI[i ] - XI[i-1]) * 2.0 / 3.0;
			APPROX(i, X, VALSTR+(KNEW-1), nullptr, ASAVE + (2) * 28,
				nullptr, XI, N, Z, DMZ, K, NCOMP, NY, 
				MMAX, MT, MSTAR, 4, nullptr, 0);
			for (int l = 0; l < MSTAR; ++l) {
				ERR[l] = WGTERR[l+1-1] * std::abs(VALSTR[KNEW-1] - VALSTR[KSTORE-1]);
				KNEW++;
				KSTORE++;
			}
			KNEW = (4 * (i - 1) + 1) * MSTAR + 1;
			KSTORE = 2 * (i - 1) * MSTAR + 1;
			X = XI[i-1] + (XI[i] - XI[i - 1]) / 3.0;
			APPROX(i, X, VALSTR+(KNEW-1), nullptr, ASAVE + (1) * 28,
				nullptr, XI, N, Z, DMZ, K, NCOMP, NY,
				MMAX, MT, MSTAR, 4, nullptr, 0);
			for (int l = 0; l < MSTAR; ++l) {
				ERR[l] += WGTERR[l+1-1] * std::abs(VALSTR[KNEW - 1] - VALSTR[KSTORE-1]);
				KNEW++;
				KSTORE++;
			}

			// find component-wise maximum error estimate
			for (int l = 0; l < MSTAR; ++l)
				ERREST[l] = std::max(ERREST[l], ERR[l]);

			// test whether the tolerance requirements are satisfied
			// in the i-th interval.
			if (IFIN == 0)
				continue;
			for (int j = 1; j <= NTOL; ++j) {
				int LTOLJ = LTOL[j-1];
				int LTJZ = LTOLJ + (i - 1) * MSTAR;
				if (ERR[LTOLJ-1] > TOLIN[j-1] * (std::abs(Z[LTJZ-1]) + 1.0))
					IFIN = 0;
			}
		}

		if (IPRINT >= 0)
			return;
		fmt::print("THE ESTIMATED ERRORS ARE\n");
		int LJ = 1;
		for (int j = 1; j <= NCOMP; ++j) {
			int MJ = LJ - 1 + MT[j-1];
			fmt::print("{}:  ", j);
			for (int l = LJ; l <= MJ; ++l)
				fmt::print("{}, ", ERREST[l-1]);
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
	//         this routine controls the set-up and solution of a linear
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
	//                 and the right-hand side  rhs ,  and solve.
	//                 (for linear problems only.)
	//      mode = 1 - set up the collocation matrices  v , w , g
	//                 and the right-hand sides  rhs  and  dmzo ,
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
	//      idmz,irhs,iv,iw - pointers to  rhs,v,w respectively
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
	void LSYSLV(int& MSING, cdvec XI, cdvec XIOLD, dvec Z, dvec DMZ, dvec DELZ, dvec DELDMZ,
		dvec G, dvec W, dvec V, dvec FC, dvec RHS, dvec DMZO,
		imat INTEGS, ivec IPVTG, ivec IPVTW, double& RNORM,
		const int MODE, fsub_t fsub, dfsub_t dfsub, gsub_t gsub, dgsub_t dgsub, guess_t guess, int& ISING)
	{
		//INTEGS(3, 1);

		double YVAL[20], ZVAL[40], F[40], DGZ[40], DMVAL[20], DF[800];
		double AT[28], CB[400];
		int IPVTCB[20];

		int NYCB = (NY == 0 ? 1 : NY);

		int INFC = (MSTAR + NY) * NCOMP;
		int M1 = MODE + 1;
		double XI1 = 0;


		switch (M1) {
		case 1: {
			//  linear problem initialization
			for (int i = 0; i < MSTAR; ++i)
				ZVAL[i] = 0.0;
			for (int i = 0; i < NY; ++i)
				YVAL[i] = 0.0;

			[[fallthrough]];
		}
		case 2: [[fallthrough]];
		case 3: [[fallthrough]];
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
				for (int i = 0; i < N; ++i) {
					INTEGS[1 + i * 3] = NCOL;
					if (i+1 >= N) {
						INTEGS[2 + (N-1) * 3] = NCOL;
						LSIDE = MSTAR;
					}
					else {
						INTEGS[2 + i * 3] = MSTAR;
						while (true) {
							if (LSIDE == MSTAR)
								break;
							if (TZETA[LSIDE] >= XI[i] + PRECIS)
								break;
							LSIDE++;
						}
					}
					int NROW = MSTAR + LSIDE;
					INTEGS[0 + i * 3] = NROW;
				}
			}
			if (MODE != 2) {
				//  zero the matrices to be computed
				int LW = KDY * KDY * N;
				for (int l = 0; l < LW; ++l)
					W[l] = 0.0;
			}
			// set up the linear system of equations
			for (int i = 0; i < N; ++i) {
				// construct a block of  a  and a corresponding piece of  rhs.
				double XII = XI[i];
				double H = XI[i + 1] - XI[i];
				int NROW = INTEGS[0 + i * 3];

				// go through the ncomp collocation equations and side conditions
				// in the i-th subinterval
				while (true) {
					if (IZETA > MSTAR)
						break;
					if (TZETA[IZETA-1] > XII + PRECIS)
						break;

					// build equation for a side condition.
					if (MODE != 0) {
						if (IGUESS == 1) {
							// case where user provided current approximation
							guess(XII, ZVAL, YVAL, DMVAL);
						}
						else {
							// other nonlinear case
							if (MODE == 1) {
								APPROX(IOLD, XII, ZVAL, nullptr, AT,
									COEF, XIOLD, NOLD, Z,
									DMZ, K, NCOMP, NY, MMAX, MT,
									MSTAR, 2, nullptr, 0);
							}
							else {
								auto temp = i + 1;
								APPROX(temp, XII, ZVAL, nullptr, AT,
									nullptr, XI, N, Z, DMZ,
									K, NCOMP, NY, MMAX, MT, MSTAR, 1, nullptr, 0);
								i = temp - 1;
								if (MODE == 3)
									goto n120;
							}
						}
					}
					// find  rhs  boundary value.
					double GVAL;
					gsub(IZETA, ZVAL, GVAL);
					RHS[NDMZ + IZETA - 1] = -GVAL;
					RNORM += GVAL * GVAL;
					if (MODE != 2) {
					n120:
						// build a row of  a  corresponding to a boundary point
						GDERIV(G + (IG - 1), NROW, IZETA, ZVAL, DGZ, 1, dgsub);
					}
					IZETA++;
				}

				// assemble collocation equations
				for (int j = 1; j <= K; ++j) {
					double HRHO = H * RHO[j-1];
					double XCOL = XII + HRHO;

					// this value corresponds to a collocation (interior)
					// point. build the corresponding  ncy  equations.
					if (MODE == 0)
						goto n200;
					if (IGUESS != 1)
						goto n160;

					// use initial approximation provided by the user.
					guess(XCOL, ZVAL, YVAL, DMZO + (IRHS - 1));
					goto n170;


				n160:
					if (MODE == 1) {
						// find  rhs  values
						APPROX(IOLD, XCOL, ZVAL, YVAL, AT,
							COEF, XIOLD,
							NOLD, Z, DMZ, K, NCOMP, NY, MMAX,
							MT, MSTAR, 2,
							DMZO + (IRHS - 1), 2);

					n170:
						fsub(XCOL, ZVAL, YVAL, F);
						for (int JJ = NCOMP; JJ < NCY; ++JJ)
							DMZO[IRHS + JJ - 1] = 0.0;

						for (int JJ = 0; JJ < NCY; ++JJ) {
							double VALUE = DMZO[IRHS - 1] - F[JJ];
							RHS[IRHS - 1] = -VALUE;
							RNORM += VALUE * VALUE;
							IRHS++;
						}
					}
					else {
						// evaluate former collocation solution
						{auto temp = i + 1;
						APPROX(temp, XCOL, ZVAL, nullptr, ACOL + (j-1) * 28,
							COEF, XI, N, Z,
							DMZ, K, NCOMP, NY, MMAX, MT, MSTAR,
							4, nullptr, 0);
						i = temp - 1;
						}
						if (MODE == 3)
							goto n210;

						// fill in  rhs  values (and accumulate its norm).
						fsub(XCOL, ZVAL, DMZ + (IRHS + NCOMP - 1), F);
						for (int JJ = 0; JJ < NCY; ++JJ) {
							double VALUE = F[JJ];
							if (JJ + 1 <= NCOMP)
								VALUE -= DMZ[IRHS - 1];
							RHS[IRHS - 1] = VALUE;
							RNORM += VALUE * VALUE;
							IRHS++;
						}
						continue;

					n200:
						// the linear case
						fsub(XCOL, ZVAL, YVAL, RHS + (IRHS - 1));
						IRHS += NCY;
					}
				n210:

					// fill in ncy rows of  w and v
					VWBLOK(XCOL, HRHO, j, W + (IW - 1), V + (IV - 1), IPVTW + (IDMZ - 1), ZVAL, YVAL, DF,
						ACOL + (j-1) * 28, DMZO + (IDMZO - 1), dfsub, MSING);

					if (MSING != 0)
						return;
				}

				// build global bvp matrix  g
				if (INDEX != 1 && NY > 0) {
					// projected collocation: find solution at xi(i+1)
					XI1 = XI[i + 1];
					if (MODE != 0) {
						if (IGUESS == 1) {
							guess(XI1, ZVAL, YVAL, DMVAL);
						}
						else {
							if (MODE == 1) {
								APPROX(IOLD, XI1, ZVAL, YVAL, AT,
									COEF,
									XIOLD, NOLD, Z, DMZ,
									K, NCOMP, NY, MMAX,
									MT, MSTAR, 2, nullptr, 1);
								if (i + 1 == N) {
									auto temp = NOLD + 1;
									APPROX(temp, XI1, ZVAL, YVAL,
										AT, COEF,
										XIOLD, NOLD, Z,
										DMZ, K, NCOMP, NY, MMAX,
										MT, MSTAR, 1, nullptr, 0);
								}
							}
							else {
								auto temp = i + 1;
								APPROX(temp, XI1, ZVAL, YVAL, AT,
									COEF,
									XI, N, Z, DMZ, K, NCOMP, NY, MMAX,
									MT, MSTAR, 3, nullptr, 1);
								i = temp - 1;

								temp = i + 1 + 1;
								APPROX(temp, XI1, ZVAL, YVAL, AT,
									COEF, XI, N, Z,
									DMZ, K, NCOMP, NY, MMAX,
									MT, MSTAR, 1, nullptr, 0);

							}
						}
					}

					// find rhs at next mesh point (also for linear case)
					fsub(XI1, ZVAL, YVAL, F);
				}

				GBLOCK(H, G + (IG - 1), NROW, IZETA, W + (IW - 1), V + (IV - 1),
					nullptr, DELDMZ + (IDMZ - 1),
					IPVTW + (IDMZ - 1), 1, MODE, XI1, ZVAL, YVAL,
					F, DF,
					CB, IPVTCB, FC + (IFC - 1), dfsub, ISING, NYCB);

				if (ISING != 0)
					return;
				if (i + 1 >= N) {
					IZSAVE = IZETA;
					while (true) {
						if (IZETA > MSTAR)
							break;

						// build equation for a side condition.
						if (MODE != 0) {
							if (IGUESS == 1) {
								// case where user provided current approximation
								guess(TRIGHT, ZVAL, YVAL, DMVAL);
							}
							else {
								// other nonlinear case
								if (MODE == 1) {
									auto temp = NOLD + 1;
									APPROX(temp, TRIGHT, ZVAL, nullptr, AT,
										COEF, XIOLD, NOLD, Z,
										DMZ, K,
										NCOMP, NY, MMAX, MT, MSTAR, 1, nullptr, 0);
								}
								else {
									auto temp = N + 1;
									APPROX(temp, TRIGHT, ZVAL, nullptr, AT,
										COEF, XI, N, Z,
										DMZ, K,
										NCOMP, NY, MMAX, MT, MSTAR, 1, nullptr, 0);
									if (MODE == 3)
										goto n260;
								}
							}
						}

						// find  rhs  boundary value.
						double GVAL;
						gsub(IZETA, ZVAL, GVAL);
						RHS[NDMZ + IZETA - 1] = -GVAL;
						RNORM = RNORM + GVAL * GVAL;
						if (MODE != 2) {
						n260:
							// build a row of  a  corresponding to a boundary point
							GDERIV(G + (IG - 1), NROW, IZETA + MSTAR, ZVAL,
								DGZ, 2, dgsub);
						}
						IZETA++;
					}
				}
				else {
					// update counters -- i-th block completed
					IG = IG + NROW * NCOL;
					IV = IV + KDY * MSTAR;
					IW = IW + KDY * KDY;
					IDMZ = IDMZ + KDY;
					if (MODE == 1)
						IDMZO += KDY;
					IFC += INFC + 2 * NCOMP;
				}
			}

			// assembly process completed
			if (MODE != 0 && MODE != 3) {
				RNORM = sqrt(RNORM / double(NZ + NDMZ));
				if (MODE == 2) {
					return;
				}
			}
			//  solve the linear system.
			//  matrix decomposition
			FCBLOK(G, INTEGS, N, IPVTG, DF, MSING);

			//  check for singular matrix
			MSING = -MSING;
			if (MSING != 0)
				return;
		} [[fallthrough]];
		case 5: {
			//  perform forward and backward substitution .
			for (int l = 0; l < NDMZ; ++l)
				DELDMZ[l] = RHS[l];

			int IZ = 1;
			int IDMZ = 1;
			int IW = 1;
			int IFC = 1;
			int IZET = 1;
			for (int i = 0; i < N; ++i) {
				int NROW = INTEGS[0 + i * 3];
				IZETA = NROW + 1 - MSTAR;
				if (i + 1 == N)
					IZETA = IZSAVE;
				while (true) {
					if (IZET == IZETA)
						break;
					DELZ[IZ - 1 + IZET - 1] = RHS[NDMZ + IZET - 1];
					IZET++;
				}

				double H = XI[i + 1] - XI[i];
				GBLOCK(H, G, NROW, IZETA, W + (IW - 1), V, DELZ + (IZ - 1),
					DELDMZ + (IDMZ - 1),
					IPVTW + (IDMZ - 1), 2, MODE, XI1, ZVAL, YVAL,
					FC + (IFC + INFC - 1),
					DF, CB, IPVTCB, FC + (IFC - 1), dfsub, ISING, NYCB);

				IZ = IZ + MSTAR;
				IDMZ = IDMZ + KDY;
				IW = IW + KDY * KDY;
				IFC = IFC + INFC + 2 * NCOMP;
				if (i + 1 < N)
					continue;
				while (true) {
					if (IZET > MSTAR)
						break;
					DELZ[IZ - 1 + IZET - 1] = RHS[NDMZ + IZET - 1];
					IZET++;
				}
			}


			//  perform forward and backward substitution for mode=0,2, or 3.
			SBBLOK(G, INTEGS, N, IPVTG, DELZ);

			//  finally find deldmz
			DMZSOL(V, DELZ, DELDMZ);

			if (MODE != 1) {
				return;
			}

			//  project current iterate into current pp-space
			for (int l = 0; l < NDMZ; ++l)
				DMZ[l] = DMZO[l];
			IZ = 1;
			IDMZ = 1;
			IW = 1;
			IFC = 1;
			IZET = 1;
			for (int i = 0; i < N; ++i) {
				int NROW = INTEGS[0+ i*3];
				IZETA = NROW + 1 - MSTAR;
				if (i+1 == N)
					IZETA = IZSAVE;
				while (true) {
					if (IZET == IZETA)
						break;
					Z[IZ - 1 + IZET-1] = DGZ[IZET-1];
					IZET = IZET + 1;
				}
				double H = XI[i + 1] - XI[i];
				GBLOCK(H, G, NROW, IZETA, W+(IW-1),
					DF, Z+(IZ-1), DMZ+(IDMZ-1), IPVTW+(IDMZ-1), 2, MODE, XI1,
					ZVAL, YVAL,	FC+(IFC + INFC + NCOMP-1),
					DF, CB, IPVTCB, FC+(IFC-1), dfsub, ISING, NYCB);
				IZ = IZ + MSTAR;
				IDMZ = IDMZ + KDY;
				IW = IW + KDY * KDY;
				IFC = IFC + INFC + 2 * NCOMP;
				if (i+1 < N)
					continue;

				while (true) {
					if (IZET > MSTAR)
						break;
					Z[IZ - 1 + IZET-1] = DGZ[IZET-1];
					IZET = IZET + 1;
				}
			}
			SBBLOK(G, INTEGS, N, IPVTG, Z);

			//  finally find dmz
			DMZSOL(V, Z, DMZ);
            break;
		}
        default: throw std::invalid_argument("LSYSLV: Invalid case.\n");
		}
	}



	//**********************************************************************
	//
	//   purpose:
	//
	//      construct a collocation matrix row according to mode:
	//      mode = 1  -  a row corresponding to an initial condition
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
	void GDERIV(dmat GI, const int NROW, const int IROW, dvec ZVAL, dvec DGZ, const int MODE, dgsub_t dgsub) const
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
	void VWBLOK(const double XCOL, const double HRHO, const int JJ, dmat WI, dmat VI, ivec IPVTW,
		cdvec ZVAL, cdvec YVAL, dmat DF, cdmat acol, dvec DMZO, dfsub_t dfsub, int& MSING)
	{
		//WI(KDY, 1);
		//VI(KDY, 1);
		//DF(NCY, 1);
		//acol(7, 4);
	
		double BASM[5];
		double HA[7*4];
	
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
				HA[j+ l*7] = FACT * acol[j + l*7];
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
			int MJ = MT[JCOMP-1];
			JN = JN + MJ;
			for (int l = 1; l <= MJ; ++l) {
				int JV = JN - l - 1;
				int JW = JCOMP - 1;
				for (int j = 1; j <= K; ++j) {
					double AJL = -HA[j-1+(l-1)*7];
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
	void PRJSVD(dmat FC, cdmat DF, dmat D, dmat U, dmat V, ivec IPVTCB, int& ISING, const int MODE)
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
			CHECK = std::max(TOLIN[i-1], CHECK);

		//  construct d and find its svd
		for (int i = 0; i < NY; ++i)
			for (int j = 0; j < NY; ++j)
				D[i + j*NY] = DF[i + NCOMP + (j + MSTAR)*NCY];

		int JOB = 11;
		dsvdc(D, NY, NY, NY, S, E, U, NY, V, NY, WORK, JOB);

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
		}
		else{
			//  form projected cb
			int IR = NY - IRANK;
			for (int i = 0; i < NY; ++i) {
				for (int j = 0; j < NY; ++j) {
					double FACT = 0;
					int ML = 0;
					for (int l = 0; l < NCOMP; ++l) {
						ML += MT[l + 1-1];
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
						MJ += MT[j+1-1];
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
	//wi - the sub - block of non-condensed collocation equations,
	//left - hand side part.
	//vi - the sub - block of non-condensed collocation equations,
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
	void GBLOCK(const double H, dmat GI, const int NROW, const int IROW,
		dvec  WI, cdmat VI, dvec RHSZ, dvec RHSDMZ, ivec IPVTW, const int MODE,
		const int MODL, const double XI1, cdvec ZVAL, cdvec YVAL, dvec F,
		dmat DF, dmat CB, ivec IPVTCB, dmat FC, dfsub_t dfsub, int& ISING, const int NYCB)
    {
        // TODO RHSZ aliases with CB
        // double const * WI (dgesl)

        AutoTimer at(g_timer, _FUNC_);

        //GI(NROW, 1);
        //VI(KDY, 1);
        //DF(NCY, 1);
        //CB(NYCB, NYCB);
        //FC(NCOMP, 1);


        double HB[7 * 4]; //HB(7, 4);
        double BASM[5], BCOL[40], U[400], V[400];

        //  compute local basis
        double FACT = 1.0;
        BASM[0] = 1.0;
        for (int l = 0; l < MMAX; ++l) {
            FACT = FACT * H / double(l + 1);
            BASM[l + 1] = FACT;
            for (int j = 0; j < K; ++j)
                HB[j + l * 7] = FACT * B[j + l * 7];
        }

        //  branch according to  m o d e
        switch (MODE) {
            case 1: {
                //  set right gi-block to identity
                if (MODL != 2) {
                    for (int j = 0; j < MSTAR; ++j) {
                        for (int IR = 0; IR < MSTAR; ++IR) {
                            GI[IROW - 1 + IR + j * NROW] = 0.0;
                            GI[IROW - 1 + IR + (MSTAR + j) * NROW] = 0.0;
                        }
                        GI[IROW - 1 + j + (MSTAR + j) * NROW] = 1.0;
                    }
                    //  compute the block gi
                    int IR = IROW;
                    for (int ICOMP = 1; ICOMP <= NCOMP; ++ICOMP) {
                        int MJ = MT[ICOMP - 1];
                        IR = IR + MJ;
                        for (int l = 0; l < MJ; ++l) {
                            int ID = IR - l - 1;
                            for (int JCOL = 0; JCOL < MSTAR; ++JCOL) {
                                int IND = ICOMP - 1;
                                double RSUM = 0.0;
                                for (int j = 0; j < K; ++j) {
                                    RSUM -= HB[j + l * 7] * VI[IND + JCOL * KDY];
                                    IND += NCY;
                                }
                                GI[ID - 1 + JCOL * NROW] = RSUM;
                            }
                            int JD = ID - IROW;
                            for (int LL = 1; LL <= l + 1; ++LL)
                                GI[ID - 1 + (JD + LL - 1) * NROW] -= BASM[LL - 1];

                        }
                    }

                    if (INDEX == 1 || NY == 0)
                        return;

                    //  projected collocation
                    //  set up projection matrix and update gi-block
                    dfsub(XI1, ZVAL, YVAL, DF);

                    //  if index=2 then form projection matrices directly
                    //  otherwise use svd to define appropriate projection
                    if (INDEX == 0) {
                        PRJSVD(FC, DF, CB, U, V, IPVTCB, ISING, 1);
                        if (ISING != 0)
                            return;
                    } else {
                        //  form  cb
                        for (int i = 0; i < NY; ++i) {
                            for (int j = 0; j < NY; ++j) {
                                FACT = 0;
                                int ML = -1;
                                for (int l = 0; l < NCOMP; ++l) {
                                    ML += MT[l + 1 - 1];
                                    FACT += DF[i + NCOMP + ML * NCY] * DF[l + (MSTAR + j) * NCY];
                                }
                                CB[i + j * NYCB] = FACT;
                            }
                        }

                        //  decompose cb
                        ISING = dgefa(CB, NY, NY, IPVTCB);
                        if (ISING != 0)
                            return;

                        //  form columns of fc
                        for (int j = 0; j < MSTAR + NY; ++j) {
                            if (j + 1 <= MSTAR) {
                                for (int i = 0; i < NY; ++i)
                                    BCOL[i] = DF[i + NCOMP + j * NCY];
                            } else {
                                for (int i = 0; i < NY; ++i)
                                    BCOL[i] = 0.0;
                                BCOL[j - MSTAR] = 1.0;
                            }

                            dgesl(CB, NY, NY, IPVTCB, BCOL, 0);

                            for (int i = 0; i < NCOMP; ++i) {
                                FACT = 0.0;
                                for (int l = 0; l < NY; ++l)
                                    FACT += DF[i + (l + MSTAR) * NCY] * BCOL[l];

                                FC[i + j * NCOMP] = FACT;
                            }
                        }
                    }

                    //  update gi
                    for (int j = 0; j < MSTAR; ++j) {
                        for (int i = 0; i < NCOMP; ++i) {
                            FACT = 0;
                            for (int l = 0; l < MSTAR; ++l)
                                FACT += FC[i + l * NCOMP] * GI[IROW - 1 + l + j * NCY];

                            BCOL[i] = FACT;
                        }
                        int ML = -1;
                        for (int i = 0; i < NCOMP; ++i) {
                            ML += MT[i + 1 - 1];
                            GI[IROW - 1 + ML + j * NCY] -= BCOL[i];
                        }
                    }
                }

                //  prepare extra rhs piece; two if new mesh
                if (INDEX == 1 || NY == 0)
                    return;
                for (int JCOL = 1; JCOL <= 2; ++JCOL) {
                    for (int i = 0; i < NCOMP; ++i) {
                        FACT = 0;
                        for (int l = 0; l < NY; ++l)
                            FACT += FC[i + (l + MSTAR) * NCOMP] * F[l + NCOMP];
                        FC[i + (JCOL + MSTAR + NY - 1) * NCOMP] = FACT;
                    }

                    if (MODL != 1 || JCOL == 2)
                        return;
                    for (int i = NCOMP; i < NY + NCOMP; ++i)
                        F[i] = 0;
                    for (int j = 0; j < MSTAR; ++j) {
                        FACT = -ZVAL[j];
                        for (int i = NCOMP; i < NY + NCOMP; ++i)
                            F[i] += DF[i + j * NCY] * FACT;
                    }
                }
                return;
            }
            case 2: {
                //  compute the appropriate piece of  rhsz
                dgesl(WI, KDY, KDY, IPVTW, RHSDMZ, 0);

                int IR = IROW;
                for (int JCOMP = 1; JCOMP <= NCOMP; ++JCOMP) {
                    int MJ = MT[JCOMP - 1];
                    IR += MJ;
                    for (int l = 1; l <= MJ; ++l) {
                        int IND = JCOMP;
                        double RSUM = 0.0;
                        for (int j = 1; j <= K; ++j) {
                            RSUM += HB[(j - 1) + (l - 1) * 7] * RHSDMZ[IND - 1];
                            IND += NCY;
                        }
                        RHSZ[IR - l - 1] = RSUM;
                    }
                }

                if (INDEX == 1 || NY == 0)
                    return;

                //  projected collocation
                //  calculate projected rhsz
                for (int i = 0; i < NCOMP; ++i) {
                    FACT = 0;
                    for (int l = 0; l < MSTAR; ++l)
                        FACT += FC[i + l * NCOMP] * RHSZ[l + IROW - 1];
                    BCOL[i] = FACT;
                }
                int ML = 0;
                for (int i = 1; i <= NCOMP; ++i) {
                    ML += MT[i - 1];
                    RHSZ[IROW - 1 + ML - 1] -= BCOL[i - 1] + F[i - 1];
                }
                break;
            }
            default: throw std::invalid_argument("GBLOCK: Invalid case.\n");
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
	void RKBAS(const double S, cdmat coef, const int k, const int M, dmat RKB, dvec DM, const int MODE)
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
					double P = coef[0+ i*K];
					for (int j = 2; j <= k; ++j)
						P = P * T[LB - j - 1] + coef[j-1 + i * K];
					RKB[i + l*7] = P;
				}
			}
			if (MODE == 0)
				return;
			for (int i = 0; i < k; ++i) {
				double P = coef[0 + i * K];
				for (int j = 2; j <= k; ++j)
					P = P * T[k - j] + coef[j-1 + i * K];
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
	void APPROX(int& i, double& X, dvec ZVAL, dvec YVAL, dmat A, cdvec coef, cdvec XI,
			const int n, cdvec Z, cdvec DMZ, const int k, const int ncomp, const int ny, const int mmax, civec M,
			const int mstar, const int MODE, dvec DMVAL, const int MODM) 
	{
		//A(7, 1);

		AutoTimer at(g_timer, _FUNC_);

		double BM[4], DM[7];
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
			RKBAS(S, coef, k, mmax, A, DM, MODM);
			}
			[[fallthrough]];
		case 4: {
            //  mode = 2, 3, or 4, compute mesh dependent rk - basis.
            BM[0] = X - XI[i - 1];

            for (int l = 2; l <= mmax; ++l)
                BM[l - 1] = BM[0] / double(l);

            //  evaluate  z(u(x)).
            int IR = 1;
            NCY = ncomp + ny;
            IZ = (i - 1) * mstar + 1;
            int IDMZ = (i - 1) * k * NCY;
            for (int JCOMP = 1; JCOMP <= ncomp; ++JCOMP) {
                int MJ = M[JCOMP - 1];
                IR = IR + MJ;
                IZ = IZ + MJ;
                for (int l = 1; l <= MJ; ++l) {
                    int IND = IDMZ + JCOMP;
                    double ZSUM = 0.0;
                    for (int j = 1; j <= k; ++j) {
                        ZSUM = ZSUM + A[(j - 1) + (l - 1) * 7] * DMZ[IND - 1];
                        IND = IND + NCY;
                    }
                    for (int LL = 1; LL <= l; ++LL) {
                        int LB = l + 1 - LL;
                        ZSUM = ZSUM * BM[LB - 1] + Z[IZ - LL - 1];
                    }
                    ZVAL[IR - l - 1] = ZSUM;
                }
            }
            if (MODM == 0)
                return;

            //  for modm = 1 evaluate  y(j) = j - th component of y.
            for (int JCOMP = 1; JCOMP <= ny; ++JCOMP)
                YVAL[JCOMP - 1] = 0.0;
            for (int j = 1; j <= k; ++j) {
                int IND = IDMZ + (j - 1) * NCY + ncomp + 1;
                double FACT = DM[j - 1];
                for (int JCOMP = 1; JCOMP <= ny; ++JCOMP) {
                    YVAL[JCOMP - 1] += FACT * DMZ[IND - 1];
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
                double FACT = DM[j - 1];
                for (int JCOMP = 1; JCOMP <= ncomp; ++JCOMP) {
                    DMVAL[JCOMP - 1] += FACT * DMZ[IND - 1];
                    IND = IND + 1;
                }
            }
            break;
        }
        default: throw std::invalid_argument("APPROX: Invalid case.\n");
		}
	}




	/* * *********************************************************************

	purpose
			solve vandermonde system v * x = e
			with  v(i, j) = rho(j) * *(i - 1) / (i - 1)!.

	* **********************************************************************/
	void VMONDE(dvec coef, int const k)
	{
		//coef(k);
	
		AutoTimer at(g_timer, _FUNC_);

		if (k == 1)
			return;
		int KM1 = k - 1;
		for (int i = 1; i <= KM1; ++i) {
			int KMI = k - i;
			for (int j = 1; j <= KMI; ++j) {
				coef[j-1] = (coef[j] - coef[j-1]) / (RHO[j + i-1] - RHO[j-1]);
			}
		}

		int IFAC = 1;
		for (int i = 1; i <= KM1; ++i) {
			int KMI = k + 1 - i;
			for (int j = 2; j <= KMI; ++j)
				coef[j-1] -= RHO[j + i - 2] * coef[j - 2];
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
	void HORDER(int const i, dvec UHIGH, double const HI, cdvec DMZ)
	{
		AutoTimer at(g_timer, _FUNC_);
	
		double DN = 1.0 / pow(HI, (K - 1));

		//  loop over the ncomp solution components
		for (int ID = 0; ID < NCOMP; ++ID)
			UHIGH[ID] = 0.0;

		int KIN = 1;
		int IDMZ = (i - 1) * K * NCY;
		for (int j = 1; j <= K;++j) {
			double FACT = DN * COEF[KIN-1];
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
	void DMZSOL(cdmat V, cdvec Z, dmat DMZ) const
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
	//	      bloks   an array that initially contains the almost block diagonal
	//	            matrix  a  to be factored, and on return contains the
	//	            computed factorization of  a .
	//	      integs  an integer array describing the block structure of  a .
	//	      nbloks  the number of blocks in  a .
	//	      ipivot  an integer array of dimension   sum(integs(3, n); n = 1,
	//		            ..., nbloks) which, on return, contains the pivoting strategy used.
	//	      scrtch  work area required, of length  max(integs(1, n); n = 1,
	//		            ..., nbloks).
	//	      info    output parameter;
	// = 0  in case matrix was found to be nonsingular.
	//            otherwise,
	// = n if the pivot element in the nth gauss step is zero.
	//
	//**********************************************************************
	static void FCBLOK(dvec BLOKS, cimat INTEGS, const int NBLOKS, ivec IPIVOT, dvec SCRTCH, int& INFO)
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
	//       info    on output, zero if the matrix is found to be non-
	//               singular, in case a zero pivot was encountered in row
	//               n, info = n on output.
	//
	// * *********************************************************************
	//
	static void FACTRB(dmat W, ivec IPIVOT, dvec D, const int NROW, const int NCOL, const int LAST, int& INFO)
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
				D[i] = std::max(D[i], std::abs(W[i + j*NROW]));

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
				if (std::abs(W[NROW-1 + (NROW-1)* NROW]) <= 0)
					INFO = k;
			}

			l = k;
			KP1 = k + 1;
			COLMAX = std::abs(W[(k-1)+(k-1)*NROW]) / D[k-1];
			// find the (relatively) largest pivot
			for (int i = KP1; i <= NROW; ++i) {
				if (std::abs(W[(i - 1) + (k - 1) * NROW]) <= COLMAX * D[i - 1])
					continue;
				COLMAX = std::abs(W[(i - 1) + (k - 1) * NROW]) / D[i - 1];
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
			if (std::abs(T) <= 0)
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
	static void SHIFTB(cdmat AI, const int NROWI, const int NCOLI, const int LAST, dmat AI1, const int NROWI1, const int NCOLI1)
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
	//     an  already generated in  fcblok. individual blocks of
	//     equations are solved via  subforand subbak .
	//
	//    parameters
	//       bloks, integs, nbloks, ipivot    are as on return from fcblok.
	//       x       on input : the right-hand side, in dense storage
	//               on output : the solution vector
	//
	//*********************************************************************
	//
	static void SBBLOK(dvec BLOKS, cimat INTEGS, const int NBLOKS, civec IPIVOT, dvec X)
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
	static void SUBFOR(cdmat W, civec IPIVOT, const int NROW, const int LAST, dvec X)
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
	//     carries out back-substitution for current block.
	//
	//    parameters
	//       w, ipivot, nrow, ncol, last  are as on return from factrb.
	//       x(1), ..., x(ncol)  contains, on input, the right side for the
	//               equations in this block after back-substitution has been
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

	static void SUBBAK(cdmat W, const int NROW, const int NCOL, const int LAST, dvec X)
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
					X[i - 1] += W[(i - 1) + (k - 1) * NROW] * T;
			}
		}
		X[0] /= W[0];
	}


};

}
