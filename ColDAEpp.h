// ColDAEpp.h: Includedatei für Include-Standardsystemdateien
// oder projektspezifische Includedateien.

#pragma once

#include <iostream>


#include <vector>
#include <cassert>
#define assertm(exp, msg) assert(((void)msg, exp))



template<typename T>
class arr2 {
public:
	int offx = 0, offy = 0;
	std::vector<std::vector<T>> data;
public:
	void assertDim(int dim0, int dim1) {
		assertm(data.size() - offx == dim0, "Dimension 0 does not match!");
		for (auto& a : data)
			assertm(a.size() - offy == dim1, "Dimension 1 does not match!");
	}
	T& operator()(int x, int y) { return data[x - 1 + offx][y - 1 + offy]; }
	arr2 sub(int ox, int oy) {
		arr2 s = *this;
		s.offx = offx + ox - 1;
		s.offy = offy + oy - 1;
		return s;
	}
};


template<typename T>
class arr1 {
public:
	int offx = 0;
	std::vector<T> data;
public:
	arr1(int size=0) {
		data.resize(size);
	}
	void assertDim(int dim0) { assertm(data.size()-offx == dim0, "Dimension 0 does not match!"); }
	T& operator()(int idx) { return data[idx - 1 + offx]; }
	arr1 sub(int ox) {
		arr1 s = *this;
		s.offx = offx+ox-1;
		return s;
	}
	void mergeSub(arr1 sub) {
		std::copy(sub.data().begin() + sub.offx, sub.data.end(), data.begin() + sub.offx);
	}
};

template<typename T>
arr1<T> wrap(T elem) {
	arr2<T> vec;
	vec.data = { elem };
	return vec;
}

template<typename T>
arr2<T> wrap(arr1<T> vec) {
	arr2<T> mat;
	mat.data = { vec.data };
	mat.offy = vec.offx;
	return mat;
}

template<typename T>
arr1<T> unwrap(arr2<T> mat) {
	arr1<T> vec;
	vec.data = mat.data[0];
	vec.offx = mat.offy;
	return vec;
}

using darr1 = arr1<double>;
using darr2 = arr2<double>;
using iarr1 = arr1<int>;
using iarr2 = arr2<int>;


//------------------------------------------------------------------------------------------------------

double DABS(double x) { return abs(x); }
double DMAX1(double x, double y) { return std::max(x,y); }
int MIN0(int x, int y) { return std::min(x, y); }

double PRECIS;
int IOUT;
int IPRINT;


double  RHO[7];
double COEF[49];

//------------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------
//                             p a r t  4
//               polynomial and service routines
//----------------------------------------------------------------------


/* * *********************************************************************

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
void RKBAS(double S, darr2 COEF, int K, int M, darr2 RKB, darr1 DM, int MODE)
{
	COEF.assertDim(K, 1);
	RKB.assertDim(7, 1);
	DM.assertDim(1);
	darr1 T(10);

	if (K == 1)
		goto n70;

	int KPM1 = K + M - 1;
	for (int I = 1; I <= KPM1; ++I)
		T(I) = S / float(I);

	for (int L = 1; L <= M; ++L) {
		int LB = K + L + 1;
		for (int I = 1; I <= K; ++I) {
			double P = COEF(1, I);
			for (int J = 2; J <= K; ++J)
				P = P * T(LB - J) + COEF(J, I);

			RKB(I, L) = P;
		}
	}
	if (MODE == 0)
		return;
	for (int I = 1; I <= K; ++I) {
		double P = COEF(1, I);
		for (int J = 2; J <= K; ++J)
			P = P * T(K + 1 - J) + COEF(J, I);
		DM(I) = P;
	}
	return;
n70:
	RKB(1, 1) = 1.00;
	DM(1) = 1.00;
}




// * *********************************************************************
//
//   purpose
//(1)       (m1 - 1)     (mncomp - 1)
//           evaluate z(u(x)) = (u(x), u(x), ..., u(x), ..., u(x))
//                              1     1         1          mncomp
//           as well as optionally y(x) and dmval(x) at one point x.
//
//   variables
//     a - array of mesh independent rk - basis coefficients
//     xi - the current mesh(having n subintervals)
//     z - the current solution vector(differential components).
//              it is convenient to imagine z as a two - dimensional
//              array with dimensions mstar x(n + 1).then
//              z(j, i) = the jth component of z at the ith mesh point
//     dmz - the array of mj - th derivatives of the current solution
//              plus algebraic solution components at collocation points
//              it is convenient to imagine dmz as a 3 - dimensional
//              array with dimensions ncy x k x n.then
//              dmz(l, j, i) = a solution value at the jth collocation
//              point in the ith mesh subinterval : if l <= ncomp then
//              dmz(l, j, i) is the ml - th derivative of ul, while if
//              l > ncomp then dmz(l, j, i) is the value of the current
//(l - ncomp)th component of y at this collocation point
//     mode - determines the amount of initialization needed
// = 4  forms z(u(x)) using z, dmzand ha
// = 3  as in = 4, but computes local rk - basis
// = 2  as in = 3, but determines i such that
//                       xi(i).le.x.lt.xi(i + 1) (unless x = xi(n + 1))
// = 1  retrieve  z = z(u(x(i)))  directly
//     modm = 0  evaluate only zval
// = 1  evaluate also yval
// = 2  evaluate in addition dmval
//   output
//     zval - the solution vector z(u(x)) (differential components)
//     yval - the solution vector y(x)  (algebraic components)
//     dmval - the mth derivatives of u(x)
//
// * *********************************************************************
void APPROX(int I, double X, darr1 ZVAL, darr1 YVAL, darr2 A, darr1 COEF, darr1 XI,
	int N, darr1 Z, darr1 DMZ, double K, int NCOMP, int NY, int MMAX, iarr1 M,
	int MSTAR, int MODE, darr1 DMVAL, int MODM)
{
	ZVAL.assertDim(1);
	DMVAL.assertDim(1);
	XI.assertDim(1);
	M.assertDim(1);
	A.assertDim(7, 1);
	darr1 DM(7);

	Z.assertDim(1);
	DMZ.assertDim(1);
	darr1 BM(4);
	COEF.assertDim(1);
	YVAL.assertDim(1);

	int IZ, ILEFT, L, IRIGHT;

	switch (MODE) {//10, 30, 80, 90), MODE
	case 1:
		//  mode = 1, retrieve  z(u(x))  directly for x = xi(i).


	n10:
		X = XI(I);
		IZ = (I - 1) * MSTAR;
		for (int J = 1; J <= MSTAR; ++J) {
			IZ = IZ + 1;
			ZVAL(J) = Z(IZ);
		}
		return;

	case 2:
		//  mode = 2, locate i so  xi(i).le.x.lt.xi(i + 1)

	n30:

		if (X >= XI(1) - PRECIS && X <= XI(N + 1) + PRECIS)
			goto n40;
		/*if (IPRINT < 1)
			WRITE(IOUT, 900) X, XI(1), XI(N + 1)*/
		if (X < XI(1))
			X = XI(1);
		if (X > XI(N + 1))
			X = XI(N + 1);
	n40:
		if (I > N || I < 1)
			I = (N + 1) / 2;
		ILEFT = I;
		if (X < XI(ILEFT))
			goto n60;
		for (L = ILEFT; L <= N; ++L) {
			I = L;
			if (X < XI(L + 1))
				goto n80;
		}
		goto n80;
	n60:
		IRIGHT = ILEFT - 1;
		for (L = 1; L <= IRIGHT; ++L) {
			I = IRIGHT + 1 - L;
			if (X >= XI(I))
				goto n80;
		}



	case 3:
		//  mode = 2 or 3, compute mesh independent rk - basis.
	n80: {


		double S = (X - XI(I)) / (XI(I + 1) - XI(I));
		RKBAS(S, wrap(COEF), K, MMAX, A, DM, MODM);
	}



	case 4:
		//  mode = 2, 3, or 4, compute mesh dependent rk - basis.
	n90:
		BM(1) = X - XI(I);

		for (L = 2; L <= MMAX; ++L)
			BM(L) = BM(1) / float(L);

		//  evaluate  z(u(x)).

		int IR = 1;
		int NCY = NCOMP + NY;
		IZ = (I - 1) * MSTAR + 1;
		int IDMZ = (I - 1) * K * NCY;
		for (int JCOMP = 1; JCOMP <= NCOMP; ++JCOMP) {
			int MJ = M(JCOMP);
			IR = IR + MJ;
			IZ = IZ + MJ;
			for (L = 1; L <= MJ; ++L) {
				int IND = IDMZ + JCOMP;
				double ZSUM = 0.0;
				for (int J = 1; J < K; ++J) {
					ZSUM = ZSUM + A(J, L) * DMZ(IND);
					IND = IND + NCY;
				}
				for (int LL = 1; LL <= L; ++L) {
					int LB = L + 1 - LL;
					ZSUM = ZSUM * BM(LB) + Z(IZ - LL);
				}
				ZVAL(IR - L) = ZSUM;
			}
		}
		if (MODM == 0)
			return;

		//  for modm = 1 evaluate  y(j) = j - th component of y.

		for (int JCOMP = 1; JCOMP <= NY; JCOMP++)
			YVAL(JCOMP) = 0.0;
		for (int J = 1; J <= K; ++K) {
			int IND = IDMZ + (J - 1) * NCY + NCOMP + 1;
			double FACT = DM(J);
			for (int JCOMP = 1; JCOMP <= NY; ++JCOMP) {
				YVAL(JCOMP) = YVAL(JCOMP) + FACT * DMZ(IND);
				IND = IND + 1;
			}
		}
		if (MODM == 1)
			return;

		//  for modm = 2 evaluate  dmval(j) = mj - th derivative of uj.

		for (int JCOMP = 1; JCOMP <= NCOMP; ++JCOMP)
			DMVAL(JCOMP) = 0.0;
		for (int J = 1; J <= K; ++J) {
			int IND = IDMZ + (J - 1) * NCY + 1;
			double FACT = DM(J);
			for (int JCOMP = 1; JCOMP <= NCOMP; ++JCOMP) {
				DMVAL(JCOMP) = DMVAL(JCOMP) + FACT * DMZ(IND);
				IND = IND + 1;
			}
		}
	}

	//	--------------------------------------------------------------------
//n900:
	//FORMAT(37H * *****DOMAIN ERROR IN APPROX * *****
	//	1 / 4H X =, D20.10, 10H   ALEFT =, D20.10,
	//	2       11H   ARIGHT =, D20.10)
	//	END
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
void APPSLN(double X, darr1 Z, darr1 Y, darr1 FSPACE, iarr1 ISPACE) {
	Z.assertDim(1);
	Y.assertDim(1);
	FSPACE.assertDim(1);
	ISPACE.assertDim(1);

	darr1 A(28);
	darr1 DUMMY(1);

	int IS6 = ISPACE(7);
	int IS5 = ISPACE(1) + 2;
	int IS4 = IS5 + ISPACE(5) * (ISPACE(1) + 1);
	int I = 1;
	APPROX(I, X, Z, Y, wrap(A), FSPACE(IS6), FSPACE(1), ISPACE(1),
		FSPACE(IS5), FSPACE(IS4), ISPACE(2), ISPACE(3),
		ISPACE(4), ISPACE(6), ISPACE(9), ISPACE(5), 2,
		DUMMY, 1);
}



/* * *********************************************************************

purpose
		solve vandermonde system v * x = e
		with  v(i, j) = rho(j) * *(i - 1) / (i - 1)!.

* **********************************************************************/
void VMONDE(darr1 RHO, darr1 COEF, int K)
{
	int I, IFAC, J, KM1, KMI;
	RHO.assertDim(K);
	COEF.assertDim(K);

	if (K == 1)
		return;
	KM1 = K - 1;
	for (int I = 1; I <= KM1; ++I) {
		KMI = K - I;
		for (int J = 1; J <= KMI; ++J) {
			COEF(J) = (COEF(J + 1) - COEF(J)) / (RHO(J + I) - RHO(J));
		}
	}

	IFAC = 1;
	for (I = 1; I <= KM1; ++J) {
		KMI = K + 1 - I;
		for (J = 2; J <= KMI; ++J)
			COEF(J) = COEF(J) - RHO(J + I - 1) * COEF(J - 1);
		COEF(KMI) = float(IFAC) * COEF(KMI);
		IFAC = IFAC * I;
	}
	COEF(1) = float(IFAC) * COEF(1);
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
void HORDER(int I, darr1 UHIGH, double HI, darr1 DMZ, int NCOMP, int NCY, int K)
{
	UHIGH.assertDim(1);
	DMZ.assertDim(1);


	double DN = 1.0 / pow(HI, (K - 1));

	//  loop over the ncomp solution components

	for (int ID = 1; ID <= NCOMP; ++ID)
		UHIGH(ID) = 0.0;

	int KIN = 1;
	int IDMZ = (I - 1) * K * NCY + 1;
	for (int J = 1; J <= K;++J) {
		double FACT = DN * COEF[KIN];
		for (int ID = 1; ID <= NCOMP;++ID) {
			UHIGH(ID) = UHIGH(ID) + FACT * DMZ(IDMZ);
			IDMZ = IDMZ + 1;
		}
		KIN = KIN + K;
	}
}



/* * *********************************************************************

purpose
		compute dmz in a blockwise manner
		dmz(i) = dmz(i) + v(i) * z(i), i = 1, ..., n

* **********************************************************************/
void  DMZSOL(int KDY, int MSTAR, int N, darr2 V, darr1 Z, darr2 DMZ)
{
	V.assertDim(KDY, 1);
	DMZ.assertDim(KDY, 1);
	Z.assertDim(1);

	int JZ = 1;
	for (int I = 1; I <= N; ++I) {
		for (int J = 1; J <= MSTAR; ++J) {
			double FACT = Z(JZ);
			for (int L = 1; L <= KDY; ++L) {
				DMZ(L, I) = DMZ(L, I) + FACT * V(L, JZ);
			}
			JZ = JZ + 1;
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
void FACTRB(darr2 W, iarr1 IPIVOT, darr1 D, int NROW, int NCOL, int LAST, int INFO)
{
	IPIVOT.assertDim(NROW);
	int I, J, K, L, KP1;
	W.assertDim(NROW, NCOL);
	D.assertDim(NROW);
	double COLMAX, T, S;// , DABS, DMAX1;


	////  initialize  d
	for (I = 1; I <= NROW; ++I)
		D(I) = 0.0;

	for (J = 1; J <= NCOL; ++J)
		for (I = 1; I <= NROW; ++I)
			D(I) = DMAX1(D(I), DABS(W(I, J)));

	////  gauss elimination with pivoting of scaled rows, loop over
	////  k = 1, ., last

	K = 1;
	////  as pivot row for k - th step, pick among the rows not yet used,
	////  i.e., from rows  k, ..., nrow, the one whose k - th entry
	////  (compared to the row size) is largest.then, if this row
	////  does not turn out to be row k, interchange row k with this
	////  particular rowand redefine ipivot(k).

n30:

	if (D(K) == 0.0) {
		INFO = K;
		return;
	}
	if (K == NROW) {
		////  if  last.eq.nrow, check now that pivot element in last row
		////  is nonzero.
		if (DABS(W(NROW, NROW)) + D(NROW) <= D(NROW))
			INFO = K;
	}

	L = K;
	KP1 = K + 1;
	COLMAX = DABS(W(K, K)) / D(K);
	////       find the(relatively) largest pivot
	for (I = KP1; I <= NROW; ++I) {
		if (DABS(W(I, K)) <= COLMAX * D(I))
			continue;
		COLMAX = DABS(W(I, K)) / D(I);
		L = I;
	}

	IPIVOT(K) = L;
	T = W(L, K);
	S = D(L);
	if (L != K) {
		W(L, K) = W(K, K);
		W(K, K) = T;
		D(L) = D(K);
		D(K) = S;
	}

	////       if pivot element is too small in absolute value, declare
	////       matrix to be noninvertible and quit.
	if (DABS(T) + D(K) <= D(K))
	{
		INFO = K; ////  singularity flag set
		return;
	}

	////       otherwise, subtract the appropriate multiple of the pivot
	////       row from remaining rows, i.e., the rows(k + 1), ..., (nrow)
	////       to make k - th entry zero.save the multiplier in its place.
	////       for high performance do this operations column oriented.
	T = -1.00 / T;
	for (I = KP1; I <= NROW; ++I)
		W(I, K) = W(I, K) * T;

	for (J = KP1; J <= NCOL; ++J) {
		T = W(L, J);
		if (L != K)
		{
			W(L, J) = W(K, J);
			W(K, J) = T;
		}
		if (T == 0.0)
			continue;
		for (I = KP1; I <= NROW; ++I)
			W(I, J) = W(I, J) + W(I, K) * T;
	}
	K = KP1;

	////       check for having reached the next block.
	if (K <= LAST)
		goto n30;
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
void SHIFTB(darr2 AI, int NROWI, int NCOLI, int LAST, darr2 AI1, int NROWI1, int NCOLI1)
{
	int J, JMAX, JMAXP1, M, MMAX;
	AI.assertDim(NROWI, NCOLI);
	AI1.assertDim(NROWI1, NCOLI1);


	MMAX = NROWI - LAST;
	JMAX = NCOLI - LAST;
	if (MMAX < 1 || JMAX < 1)
		return;

	////  put the remainder of block i into ai1
	for (J = 1; J <= JMAX; ++J)
		for (M = 1; M <= MMAX; ++M)
			AI1(M, J) = AI(LAST + M, LAST + J);


	if (JMAX == NCOLI1)
		return;

	////  zero out the upper right corner of ai1

	JMAXP1 = JMAX + 1;
	for (J = JMAXP1; J <= NCOLI1; ++J)
		for (M = 1; MMAX; ++M)
			AI1(M, J) = 0.0;
}



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
void FCBLOK(darr1 BLOKS, iarr2 INTEGS, int NBLOKS, iarr1 IPIVOT, darr1 SCRTCH, int INFO)
{
	INTEGS.assertDim(3, NBLOKS);
	IPIVOT.assertDim(1);
	int I, INDEX, INDEXN, LAST, NCOL, NROW;
	BLOKS.assertDim(1);
	SCRTCH.assertDim(1);

	INFO = 0;
	int INDEXX = 1;
	INDEXN = 1;
	I = 1;

	//  loop over the blocks. i is loop index

n10:
	INDEX = INDEXN;
	NROW = INTEGS(1, I);
	NCOL = INTEGS(2, I);
	LAST = INTEGS(3, I);

	//       carry out elimination on the i - th block until next block
	//       enters, i.e., for columns 1, ..., last  of i - th block.

	FACTRB(wrap(BLOKS.sub(INDEX)), IPIVOT(INDEXX), SCRTCH, NROW, NCOL, LAST, INFO);

	//       check for having reached a singular block or the last block

	if (INFO != 0)                       goto n20;
	if (I == NBLOKS)                     return;
	I = I + 1;
	INDEXN = NROW * NCOL + INDEX;
	INDEXX = INDEXX + LAST;

	//       put the rest of the i - th block onto the next block

	SHIFTB(wrap(BLOKS.sub(INDEX)), NROW, NCOL, LAST, wrap(BLOKS.sub(INDEXN)), INTEGS(1, I), INTEGS(2, I));

	goto n10;


n20:
	INFO = INFO + INDEXX - 1;
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
// W(NROW, NCOL), X(NCOL)
void SUBBAK(darr2 W, int NROW, int NCOL, int LAST, darr1 X)
{
	W.assertDim(NROW, NCOL);
	X.assertDim(NCOL);
	int I, J, K, KM1, LM1, LP1, KB;
	double T;

	LP1 = LAST + 1;
	if (LP1 <= NCOL)
		for (J = LP1; J <= NCOL; ++J) {
			T = -X(J);
			if (T == 0.0)
				continue;
			for (I = 1; I <= LAST; ++I)
				X(I) = X(I) + W(I, J) * T;
		}

	if (LAST == 1)
		goto n60;
	LM1 = LAST - 1;

	for (KB = 1; KB <= LM1; ++KB) {
		KM1 = LAST - KB;
		K = KM1 + 1;
		X(K) = X(K) / W(K, K);
		T = -X(K);
		if (T == 0.0)
			continue;
		for (I = 1; I <= KM1; ++I)
			X(I) = X(I) + W(I, K) * T;
	}

n60:
	X(1) = X(1) / W(1, 1);

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
//int IPIVOT(LAST), W(NROW, LAST), X(NROW),
void SUBFOR(darr2 W, iarr1 IPIVOT, int NROW, int LAST, darr1 X)
{
	IPIVOT.assertDim(LAST);
	W.assertDim(NROW, LAST);
	X.assertDim(NROW);
	int IP, K, KP1, LSTEP;
	double T;

	if (NROW == 1)
		return;
	LSTEP = MIN0(NROW - 1, LAST);
	for (K = 1; K <= LSTEP; ++K)
	{
		KP1 = K + 1;
		IP = IPIVOT(K);
		T = X(IP);
		X(IP) = X(K);
		X(K) = T;
		if (T == 0.0)
			continue;
		for (int i = KP1; i <= NROW; ++i)
			X(i) = X(i) + W(i, K) * T;
	}
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
void SBBLOK(darr1 BLOKS, iarr2 INTEGS, int NBLOKS, iarr1 IPIVOT, darr1 X)
{
	INTEGS.assertDim(3, NBLOKS);
	IPIVOT.assertDim(1);
	BLOKS.assertDim(1);
	X.assertDim(1);
	int I, INDEX, INDEXX, J, LAST, NBP1, NCOL, NROW;

	////  forward substitution pass
	INDEX = 1;
	INDEXX = 1;
	for (I = 1; I <= NBLOKS; ++I) {
		NROW = INTEGS(1, I);
		LAST = INTEGS(3, I);
		SUBFOR(wrap(BLOKS.sub(INDEX)), IPIVOT.sub(INDEXX), NROW, LAST, X.sub(INDEXX));
		INDEX = NROW * INTEGS(2, I) + INDEX;
		INDEXX = INDEXX + LAST;
	}

	////  back substitution pass
	NBP1 = NBLOKS + 1;
	for (J = 1; J <= NBLOKS; ++J) {
		I = NBP1 - J;
		NROW = INTEGS(1, I);
		NCOL = INTEGS(2, I);
		LAST = INTEGS(3, I);
		INDEX = INDEX - NROW * NCOL;
		INDEXX = INDEXX - LAST;
	}
	SUBBAK(wrap(BLOKS.sub(INDEX)), NROW, NCOL, LAST, X.sub(INDEXX));
}


