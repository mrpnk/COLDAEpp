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
    arr1(int size = 0) {
        data.resize(size);
    }
    void assertDim(int dim0) { assertm(data.size() - offx == dim0, "Dimension 0 does not match!"); }
    T& operator()(int idx) { return data[idx - 1 + offx]; }
    arr1 sub(int ox) {
        arr1 s = *this;
        s.offx = offx + ox - 1;
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
double DMAX1(double x, double y) { return std::max(x, y); }
int MIN0(int x, int y) { return std::min(x, y); }

double PRECIS;
int IOUT;
int IPRINT;


double  RHO[7];
double COEF[49];

//------------------------------------------------------------------------------------------------------






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
//
//**********************************************************************
//
//     ***************     input to coldae     ***************
//
//     variables
//
//     ncomp - no. of differential equations   (ncomp .le. 20)
//
//     ny   - no. of constraints     (ny .le. 20)
//
//     m(j) - order of the j-th differential equation
//            ( mstar = m(1) + ... + m(ncomp) .le. 40 )
//
//     aleft - left end of interval
//
//     aright - right end of interval
//
//     zeta(j) - j-th side condition point (boundary point). must
//               have  zeta(j) .le. zeta(j+1). all side condition
//               points must be mesh points in all meshes used,
//               see description of ipar(11) and fixpnt below.
//
//     ipar - an integer array dimensioned at least 11.
//            a list of the parameters in ipar and their meaning follows
//            some parameters are renamed in coldae; these new names are
//            given in parentheses.
//
//     ipar(1)     ( = nonlin )
//             = 0 if the problem is linear
//             = 1 if the problem is nonlinear
//
//     ipar(2) = no. of collocation points per subinterval  (= k )
//               where max m(i) .le.  k .le. 7 . if ipar(2)=0 then
//               coldae sets  k = max ( max m(i)+1, 5-max m(i) )
//
//     ipar(3) = no. of subintervals in the initial mesh  ( = n ).
//               if ipar(3) = 0 then coldae arbitrarily sets n = 5.
//
//     ipar(4) = no. of solution and derivative tolerances.  ( = ntol )
//               we require  0 .lt. ntol .le. mstar.
//
//     ipar(5) = dimension of fspace.     ( = ndimf )
//
//     ipar(6) = dimension of ispace.     ( = ndimi )
//
//     ipar(7) -  output control ( = iprint )
//              = -1 for full diagnostic printout
//              = 0 for selected printout
//              = 1 for no printout
//
//     ipar(8)     ( = iread )
//             = 0 causes coldae to generate a uniform initial mesh.
//             = 1 if the initial mesh is provided by the user.  it
//                 is defined in fspace as follows:  the mesh
//                 aleft=x(1).lt.x(2).lt. ... .lt.x(n).lt.x(n+1)=aright
//                 will occupy  fspace(1), ..., fspace(n+1). the
//                 user needs to supply only the interior mesh
//                 points  fspace(j) = x(j), j = 2, ..., n.
//             = 2 if the initial mesh is supplied by the user
//                 as with ipar(8)=1, and in addition no adaptive
//                 mesh selection is to be done.
//
//     ipar(9)     ( = iguess )
//             = 0 if no initial guess for the solution is
//                 provided.
//             = 1 if an initial guess is provided by the user
//                 in subroutine  guess.
//             = 2 if an initial mesh and approximate solution
//                 coefficients are provided by the user in  fspace.
//                 (the former and new mesh are the same).
//             = 3 if a former mesh and approximate solution
//                 coefficients are provided by the user in fspace,
//                 and the new mesh is to be taken twice as coarse;
//                 i.e.,every second point from the former mesh.
//             = 4 if in addition to a former initial mesh and
//                 approximate solution coefficients, a new mesh
//                 is provided in fspace as well.
//                 (see description of output for further details
//                 on iguess = 2, 3, and 4.)
//
//     ipar(10)= -1 if the first relax factor is RSTART
//		(use for an extra sensitive nonlinear problem only)
//	      =  0 if the problem is regular
//	      =  1 if the newton iterations are not to be damped
//                 (use for initial value problems).
//	      =  2 if we are to return immediately upon  (a) two
//                 successive nonconvergences, or  (b) after obtaining
//                 error estimate for the first time.
//
//     ipar(11)= no. of fixed points in the mesh other than aleft
//               and aright. ( = nfxpnt , the dimension of fixpnt)
//               the code requires that all side condition points
//               other than aleft and aright (see description of
//               zeta ) be included as fixed points in fixpnt.
//
//     ipar(12)    ( = index )
//                 this parameter is ignored if ny=0.
//             = 0 if the index of the dae is not as per one of the
//                 following cases
//             = 1 if the index of the dae is 1. in this case the
//                 ny x ny jacobian matrix of the constraints with
//                 respect to the algebraic unknowns, i.e.
//                 df(i,j), i=ncomp+1,...,ncomp+ny,
//                          j=mstar+1,...,mstar+ny
//                 (see description of dfsub below)
//                 is nonsingular wherever it is evaluated. this
//                 allows usual collocation to be safely used.
//             = 2 if the index of the dae is 2 and it is in Hessenberg
//                 form. in this case the
//                 ny x ny jacobian matrix of the constraints with
//                 respect to the algebraic unknowns is 0, and the
//                 ny x ny matrix  CB  is nonsingular, where
//                 (i,j) = df(i+ncomp, m(j)), i=1,...,ny,
//                                             j=1,...,ncomp
//                 B(i,j) = df(i,j+mstar),     i=1,...,ncomp,
//                                             j=1,...,ny
//                 the projected collocation method described in [2]
//                 is then used.
//             in case of ipar(12)=0 and ny > 0, coldae determines the
//             appropriate projection needed at the right end of each
//             mesh subinterval using SVD. this is the most expensive
//             and most general option.
//
//     ltol  -  an array of dimension ipar(4). ltol(j) = l  specifies
//              that the j-th tolerance in  tol  controls the error
//              in the l-th component of z(u).   also require that
//              1.le.ltol(1).lt.ltol(2).lt. ... .lt.ltol(ntol).le.mstar
//
//     tol    - an array of dimension ipar(4). tol(j) is the
//              error tolerance on the ltol(j) -th component
//              of z(u). thus, the code attempts to satisfy
//              for j=1,...,ntol  on each subinterval
//              abs(z(v)-z(u))       .le. tol(j)*abs(z(u))       +tol(j)
//                            ltol(j)                     ltol(j)
//
//              if v(x) is the approximate solution vector.
//
//     fixpnt - an array of dimension ipar(11).   it contains
//              the points, other than aleft and aright, which
//              are to be included in every mesh.
//
//     ispace - an integer work array of dimension ipar(6).
//              its size provides a constraint on nmax,
//              the maximum number of subintervals. choose
//              ipar(6) according to the formula
//                      ipar(6)  .ge.  nmax*nsizei
//                where
//                      nsizei = 3 + kdm
//                with
//                      kdm = kdy + mstar  ;  kdy = k * (ncomp+ny) ;
//                      nrec = no. of right end boundary conditions.
//
//
//     fspace - a real work array of dimension ipar(5).
//              its size provides a constraint on nmax.
//              choose ipar(5) according to the formula
//                      ipar(5)  .ge.  nmax*nsizef
//                where
//                      nsizef = 4 + 3 * mstar + (5+kdy) * kdm +
//                              (2*mstar-nrec) * 2*mstar +
//                              ncomp*(mstar+ny+2) + kdy.
//
//
//     iflag - the mode of return from coldae.
//           = 1 for normal return
//           = 0 if the collocation matrix is singular.
//           =-1 if the expected no. of subintervals exceeds storage
//               specifications.
//           =-2 if the nonlinear iteration has not converged.
//           =-3 if there is an input data error.
//
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
void COLDAE(NCOMP, NY, M, ALEFT, ARIGHT, ZETA, IPAR, LTOL,
    1                   TOL, FIXPNT, ISPACE, FSPACE, IFLAG,
    2                   FSUB, DFSUB, GSUB, DGSUB, GUESS)
{
    IMPLICIT DOUBLE PRECISION(A - H, O - Z)
        DIMENSION M(1), ZETA(1), IPAR(1), LTOL(1), TOL(1), DUMMY(1),
        1          FIXPNT(1), ISPACE(1), FSPACE(1)
        DIMENSION DUMMY2(840)

        COMMON / COLOUT / PRECIS, IOUT, IPRINT
        COMMON / COLLOC / RHO(7), COEF(49)
        COMMON / COLORD / K, NC, NNY, NCY, MSTAR, KD, KDY, MMAX, MT(20)
        COMMON / COLAPR / N, NOLD, NMAX, NZ, NDMZ
        COMMON / COLMSH / MSHFLG, MSHNUM, MSHLMT, MSHALT
        COMMON / COLSID / TZETA(40), TLEFT, TRIGHT, IZETA, IDUM
        COMMON / COLNLN / NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX
        COMMON / COLEST / TTL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
        1                ROOT(40), JTOL(40), LTTOL(40), NTOL

        EXTERNAL FSUB, DFSUB, GSUB, DGSUB, GUESS

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

        if (IPAR(7) <= 0)  WRITE(6, 99)
            99  FORMAT(//,' VERSION *1* OF COLDAE .'    ,//)

                IOUT = 6
                PRECIS = 1.D0
                10 PRECIS = PRECIS / 2.D0
                PRECP1 = PRECIS + 1.D0
                if (PRECP1 > 1.D0)                       goto 10
                    PRECIS = PRECIS * 100.D0

                    //  in case incorrect input data is detected, the program returns
                    //  immediately with iflag=-3.

                    IFLAG = -3
                    NCY = NCOMP + NY
                    if (NCOMP < 0 || NCOMP > 20)        return;
    if (NY < 0 || NY > 20)              return;
    if (NCY < 1 || NCY > 40)            return;
    for (20 I = 1, NCOMP
        if (M(I) < 1 || M(I) > 4)        return;
        20 CONTINUE

            //  rename some of the parameters and set default values.

            NONLIN = IPAR(1)
            K = IPAR(2)
            N = IPAR(3)
            if (N == 0)  N = 5
                IREAD = IPAR(8)
                IGUESS = IPAR(9)
                if (NONLIN == 0 && IGUESS == 1)  IGUESS = 0
                    if (IGUESS >= 2 && IREAD == 0)   IREAD = 1
                        ICARE = IPAR(10)
                        NTOL = IPAR(4)
                        NDIMF = IPAR(5)
                        NDIMI = IPAR(6)
                        NFXPNT = IPAR(11)
                        IPRINT = IPAR(7)
                        INDEX = IPAR(12)
                        if (NY == 0) INDEX = 0
                            MSTAR = 0
                            MMAX = 0
                            for (30 I = 1, NCOMP
                                MMAX = MAX0(MMAX, M(I))
                                MSTAR = MSTAR + M(I)
                                MT(I) = M(I)
                                30 CONTINUE
                                if (K == 0)   K = MAX0(MMAX + 1, 5 - MMAX)
                                    for (40 I = 1, MSTAR
                                        40 TZETA(I) = ZETA(I)
                                        for (50 I = 1, NTOL
                                            LTTOL(I) = LTOL(I)
                                            50 TOLIN(I) = TOL(I)
                                            TLEFT = ALEFT
                                            TRIGHT = ARIGHT
                                            NC = NCOMP
                                            NNY = NY
                                            KD = K * NCOMP
                                            KDY = K * NCY

                                            //  print the input data for checking.

                                            if (IPRINT > -1)                         goto 80
                                                if (NONLIN == 0) THEN
                                                    WRITE(IOUT, 260) NCOMP, (M(IP), IP = 1, NCOMP)
                                                else
                                                    WRITE(IOUT, 270) NCOMP, (M(IP), IP = 1, NCOMP)
                                                    END if
                                                    WRITE(IOUT, 275) NY
                                                    if (NY > 0 && INDEX == 0) THEN
                                                        WRITE(IOUT, 276)
                                                    else
                                                        WRITE(IOUT, 278) INDEX
                                                        END if
                                                        WRITE(IOUT, 279)
                                                        WRITE(IOUT, 280) (ZETA(IP), IP = 1, MSTAR)
                                                        if (NFXPNT > 0)
                                                            1   WRITE(IOUT, 340) NFXPNT, (FIXPNT(IP), IP = 1, NFXPNT)
                                                            WRITE(IOUT, 290) K
                                                            WRITE(IOUT, 299)
                                                            WRITE(IOUT, 300) (LTOL(IP), IP = 1, NTOL)
                                                            WRITE(IOUT, 309)
                                                            WRITE(IOUT, 310) (TOL(IP), IP = 1, NTOL)
                                                            if (IGUESS >= 2) WRITE(IOUT, 320)
                                                                if (IREAD == 2) WRITE(IOUT, 330)
                                                                    80 CONTINUE

                                                                    //  check for correctness of data

                                                                    if (K < 0 || K > 7)                 return;
    if (N < 0)                               return;
    if (IREAD < 0 || IREAD > 2)         return;
    if (IGUESS < 0 || IGUESS > 4)       return;
    if (ICARE < -1 || ICARE > 2)	   return;
    if (INDEX < 0 || INDEX > 2)         return;
    if (NTOL < 0 || NTOL > MSTAR)       return;
    if (NFXPNT < 0)                          return;
    if (IPRINT < (-1) || IPRINT > 1)    return;
    if (MSTAR < 0 || MSTAR > 40)        return;
    IP = 1
        for (100 I = 1, MSTAR
            if (DABS(ZETA(I) - ALEFT) < PRECIS ||
                1     DABS(ZETA(I) - ARIGHT) < PRECIS)     goto 100
                90 if (IP > NFXPNT)                         return;
            if (ZETA(I) - PRECIS < FIXPNT(IP))     goto 95
                IP = IP + 1
                goto 90
                95 if (ZETA(I) + PRECIS < FIXPNT(IP))       return;
            100 CONTINUE

                //  set limits on iterations and initialize counters.
                //  limit = maximum number of newton iterations per mesh.
                //  see subroutine  newmsh  for the roles of  mshlmt , mshflg ,
                //  mshnum , and  mshalt .

                MSHLMT = 3
                MSHFLG = 0
                MSHNUM = 1
                MSHALT = 1
                LIMIT = 40

                //  compute the maxium possible n for the given sizes of
                //  ispace  and  fspace.

                NREC = 0
                for (110 I = 1, MSTAR
                    IB = MSTAR + 1 - I
                    if (ZETA(IB) >= ARIGHT)  NREC = I
                        110 CONTINUE
                        NFIXI = MSTAR
                        NSIZEI = 3 + KDY + MSTAR
                        NFIXF = NREC * (2 * MSTAR) + 5 * MSTAR + 3
                        NSIZEF = 4 + 3 * MSTAR + (KDY + 5) * (KDY + MSTAR) +
                        1(2 * MSTAR - NREC) * 2 * MSTAR + (MSTAR + NY + 2) * NCOMP + KDY
                        NMAXF = (NDIMF - NFIXF) / NSIZEF
                        NMAXI = (NDIMI - NFIXI) / NSIZEI
                        if (IPRINT < 1)  THEN
                            WRITE(IOUT, 350)
                            WRITE(IOUT, 351) NMAXF, NMAXI
                            ENDIF
                            NMAX = MIN0(NMAXF, NMAXI)
                            if (NMAX < N)                            return;
    if (NMAX < NFXPNT + 1)                     return;
    if (NMAX < 2 * NFXPNT + 2 && IPRINT < 1)  WRITE(IOUT, 360)

        //  generate pointers to break up  fspace  and  ispace .

        LXI = 1
        LG = LXI + NMAX + 1
        LXIOLD = LG + 2 * MSTAR * (NMAX * (2 * MSTAR - NREC) + NREC)
        LW = LXIOLD + NMAX + 1
        LV = LW + KDY * *2 * NMAX
        LFC = LV + MSTAR * KDY * NMAX
        LZ = LFC + (MSTAR + NY + 2) * NCOMP * NMAX
        LDMZ = LZ + MSTAR * (NMAX + 1)
        LDMV = LDMZ + KDY * NMAX
        LDELZ = LDMV + KDY * NMAX
        LDELDZ = LDELZ + MSTAR * (NMAX + 1)
        LDQZ = LDELDZ + KDY * NMAX
        LDQDMZ = LDQZ + MSTAR * (NMAX + 1)
        LRHS = LDQDMZ + KDY * NMAX
        LVALST = LRHS + KDY * NMAX + MSTAR
        LSLOPE = LVALST + 4 * MSTAR * NMAX
        LACCUM = LSLOPE + NMAX
        LSCL = LACCUM + NMAX + 1
        LDSCL = LSCL + MSTAR * (NMAX + 1)
        LPVTG = 1
        LPVTW = LPVTG + MSTAR * (NMAX + 1)
        LINTEG = LPVTW + KDY * NMAX

        //  if  iguess .ge. 2, move  xiold, z, and  dmz  to their proper
        //  locations in  fspace.

        if (IGUESS < 2)                          goto 160
            NOLD = N
            if (IGUESS == 4)  NOLD = ISPACE(1)
                NZ = MSTAR * (NOLD + 1)
                NDMZ = KDY * NOLD
                NP1 = N + 1
                if (IGUESS == 4)  NP1 = NP1 + NOLD + 1
                    for (120 I = 1, NZ
                        120 FSPACE(LZ + I - 1) = FSPACE(NP1 + I)
                        IDMZ = NP1 + NZ
                        for (125 I = 1, NDMZ
                            125 FSPACE(LDMZ + I - 1) = FSPACE(IDMZ + I)
                            NP1 = NOLD + 1
                            if (IGUESS == 4)                          goto 140
                                for (130 I = 1, NP1
                                    130 FSPACE(LXIOLD + I - 1) = FSPACE(LXI + I - 1)
                                    goto 160
                                    140 for (150 I = 1, NP1
                                        150 FSPACE(LXIOLD + I - 1) = FSPACE(N + 1 + I)
                                        160 CONTINUE

                                        //  initialize collocation points, constants, mesh.

                                        CALL CONSTS(K, RHO, COEF)
                                        if (NY == 0) THEN
                                            NYCB = 1
                                        else
                                            NYCB = NY
                                            ENDIF
                                            CALL NEWMSH(3 + IREAD, FSPACE(LXI), FSPACE(LXIOLD), DUMMY,
                                                1             DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, NFXPNT, FIXPNT,
                                                2		   DUMMY2, DFSUB, DUMMY2, DUMMY2, NCOMP, NYCB)

                                            //  determine first approximation, if the problem is nonlinear.

                                            if (IGUESS >= 2)                            goto 230
                                                NP1 = N + 1
                                                for (210 I = 1, NP1
                                                    210 FSPACE(I + LXIOLD - 1) = FSPACE(I + LXI - 1)
                                                    NOLD = N
                                                    if (NONLIN == 0 || IGUESS == 1)      goto 230

                                                        //  system provides first approximation of the solution.
                                                        //  choose z(j) = 0  for j=1,...,mstar.

                                                        for (220 I = 1, NZ
                                                            220 FSPACE(LZ - 1 + I) = 0.D0
                                                            for (225 I = 1, NDMZ
                                                                225 FSPACE(LDMZ - 1 + I) = 0.D0
                                                                230 CONTINUE
                                                                if (IGUESS >= 2)  IGUESS = 0
                                                                    CALL CONTRL(FSPACE(LXI), FSPACE(LXIOLD), FSPACE(LZ), FSPACE(LDMZ),
                                                                        +FSPACE(LDMV),
                                                                        1     FSPACE(LRHS), FSPACE(LDELZ), FSPACE(LDELDZ), FSPACE(LDQZ),
                                                                        2     FSPACE(LDQDMZ), FSPACE(LG), FSPACE(LW), FSPACE(LV), FSPACE(LFC),
                                                                        3     FSPACE(LVALST), FSPACE(LSLOPE), FSPACE(LSCL), FSPACE(LDSCL),
                                                                        4     FSPACE(LACCUM), ISPACE(LPVTG), ISPACE(LINTEG), ISPACE(LPVTW),
                                                                        5     NFXPNT, FIXPNT, IFLAG, FSUB, DFSUB, GSUB, DGSUB, GUESS)

                                                                    //  prepare output

                                                                    ISPACE(1) = N
                                                                    ISPACE(2) = K
                                                                    ISPACE(3) = NCOMP
                                                                    ISPACE(4) = NY
                                                                    ISPACE(5) = MSTAR
                                                                    ISPACE(6) = MMAX
                                                                    ISPACE(7) = NZ + NDMZ + N + 2
                                                                    K2 = K * K
                                                                    ISPACE(8) = ISPACE(7) + K2 - 1
                                                                    for (240 I = 1, NCOMP
                                                                        240 ISPACE(8 + I) = M(I)
                                                                        for (250 I = 1, NZ
                                                                            250 FSPACE(N + 1 + I) = FSPACE(LZ - 1 + I)
                                                                            IDMZ = N + 1 + NZ
                                                                            for (255 I = 1, NDMZ
                                                                                255 FSPACE(IDMZ + I) = FSPACE(LDMZ - 1 + I)
                                                                                IC = IDMZ + NDMZ
                                                                                for (258 I = 1, K2
                                                                                    258 FSPACE(IC + I) = COEF(I)
                                                                                    return;
                                                                                    ----------------------------------------------------------------------
                                                                                    260 FORMAT(/// ' THE NUMBER OF (LINEAR) DIFF EQNS IS' , I3/ 1X,
                                                                                        1       16HTHEIR ORDERS ARE, 20I3)
                                                                                    270 FORMAT(/// 40H THE NUMBER OF (NONLINEAR) DIFF EQNS IS , I3/ 1X,
                                                                                        1       16HTHEIR ORDERS ARE, 20I3)
                                                                                    275 FORMAT(' THERE ARE', I3, ' ALGEBRAIC CONSTRAINTS')
                                                                                    276 FORMAT(' THE PROBLEM HAS MIXED INDEX CONSTRAINTS')
                                                                                    278 FORMAT(' THE INDEX IS', I2)
                                                                                    279 FORMAT(27H SIDE CONDITION POINTS ZETA)
                                                                                    280 FORMAT(8F10.6, 4(/ 8F10.6))
                                                                                    290 FORMAT(37H NUMBER OF COLLOC PTS PER INTERVAL IS, I3)
                                                                                    299 FORMAT(40H COMPONENTS OF Z REQUIRING TOLERANCES - )
                                                                                    300 FORMAT(100(/ 8I10))
                                                                                    309 FORMAT(33H CORRESPONDING ERROR TOLERANCES - )
                                                                                    310 FORMAT(8D10.2)
                                                                                    320 FORMAT(44H INITIAL MESH(ES) AND Z, DMZ PROVIDED BY USER)
                                                                                    330 FORMAT(27H NO ADAPTIVE MESH SELECTION)
                                                                                    340 FORMAT(10H THERE ARE, I5, 27H FIXED POINTS IN THE MESH - ,
                                                                                        1       10(6F10.6 / ))
                                                                                    350 FORMAT(35H THE MAXIMUM NUMBER OF SUBINTERVALS)
                                                                                    351 FORMAT(9H IS MIN(, I4, 23H(ALLOWED FROM FSPACE), , I4,
                                                                                        1       24H(ALLOWED FROM ISPACE)))
                                                                                    360 FORMAT(/ 53H INSUFFICIENT SPACE TO DOUBLE MESH FOR ERROR ESTIMATE)
}


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
void CONTRL(XI, XIOLD, Z, DMZ, DMV, RHS, DELZ, DELDMZ,
    1           DQZ, DQDMZ, G, W, V, FC, VALSTR, SLOPE, SCALE, DSCALE,
    2           ACCUM, IPVTG, INTEGS, IPVTW, NFXPNT, FIXPNT, IFLAG,
    3           FSUB, DFSUB, GSUB, DGSUB, GUESS)
{


    IMPLICIT DOUBLE PRECISION(A - H, O - Z)
        DIMENSION XI(1), XIOLD(1), Z(1), DMZ(1), RHS(1), DMV(1)
        DIMENSION G(1), W(1), V(1), VALSTR(1), SLOPE(1), ACCUM(1)
        DIMENSION DELZ(1), DELDMZ(1), DQZ(1), DQDMZ(1), FIXPNT(1)
        DIMENSION DUMMY(1), SCALE(1), DSCALE(1), FC(1), DF(800)
        DIMENSION FCSP(40, 60), CBSP(20, 20)
        DIMENSION INTEGS(1), IPVTG(1), IPVTW(1)

        COMMON / COLOUT / PRECIS, IOUT, IPRINT
        COMMON / COLORD / K, NCOMP, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
        COMMON / COLAPR / N, NOLD, NMAX, NZ, NDMZ
        COMMON / COLMSH / MSHFLG, MSHNUM, MSHLMT, MSHALT
        COMMON / COLSID / ZETA(40), ALEFT, ARIGHT, IZETA, IDUM
        COMMON / COLNLN / NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX
        COMMON / COLEST / TOL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
        1                ROOT(40), JTOL(40), LTOL(40), NTOL

        EXTERNAL FSUB, DFSUB, GSUB, DGSUB, GUESS

        //  constants for control of nonlinear iteration

        RELMIN = 1.D - 3
        RSTART = 1.D - 2
        LMTFRZ = 4

        //  compute the maximum tolerance

        CHECK = 0.D0
        for (10 I = 1, NTOL
            10   CHECK = DMAX1(TOLIN(I), CHECK)
            IMESH = 1
            ICONV = 0
            if (NONLIN == 0) ICONV = 1
                ICOR = 0
                NOCONV = 0
                MSING = 0
                ISING = 0

                //  the main iteration begins here .
                //  loop 20 is executed until error tolerances are satisfied or
                //  the code fails (due to a singular matrix or storage limitations)

                20      CONTINUE

                //       initialization for a new mesh

                ITER = 0
                if (NONLIN > 0)                     goto 50

                    //       the linear case.
                    //       set up and solve equations

                    CALL LSYSLV(MSING, XI, XIOLD, DUMMY, DUMMY, Z, DMZ, G,
                        1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 0,
                        2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING)

                    //       check for a singular matrix

                    if (ISING != 0) THEN
                        if (IPRINT < 1)  WRITE(IOUT, 497)
                            IFLAG = 0
                            return;
    END if
        if (MSING == 0)                      goto 400
            30      if (MSING < 0)                      goto 40
            if (IPRINT < 1)  WRITE(IOUT, 495)
                goto 460
                40      if (IPRINT < 1)  WRITE(IOUT, 490)
                IFLAG = 0
                return;

    //       iteration loop for nonlinear case
    //       define the initial relaxation parameter (= relax)

    50      RELAX = 1.D0

        //       check for previous convergence and problem sensitivity

        if (ICARE == (-1))  RELAX = RSTART
            if (ICARE == 1)     RELAX = 1.D0
                if (ICONV == 0)                      goto 160

                    //       convergence on a previous mesh has been obtained.    thus
                    //       we have a very good initial approximation for the newton
                    //       process.    proceed with one full newton and then iterate
                    //       with a fixed jacobian.

                    IFREEZ = 0

                    //       evaluate right hand side and its norm  and
                    //       find the first newton correction

                    CALL LSYSLV(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
                        1          W, V, FC, RHS, DQDMZ, INTEGS, IPVTG, IPVTW, RNOLD, 1,
                        2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING)

                    if (IPRINT < 0)  WRITE(IOUT, 530)
                        if (IPRINT < 0)  WRITE(IOUT, 510) ITER, RNOLD
                            goto 70

                            //       solve for the next iterate .
                            //       the value of ifreez determines whether this is a full
                            //       newton step (=0) or a fixed jacobian iteration (=1).

                            60      if (IPRINT < 0)  WRITE(IOUT, 510) ITER, RNORM
                            RNOLD = RNORM
                            CALL LSYSLV(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
                                1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM,
                                2          3 + IFREEZ, FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING)

                            //       check for a singular matrix

                            70      if (MSING != 0)                      goto 30
                            if (ISING != 0) THEN
                                if (IPRINT < 1)  WRITE(IOUT, 497)
                                    IFLAG = 0
                                    return;
    END if
        if (IFREEZ == 1)                     goto 80

            //       a full newton step

            ITER = ITER + 1
            IFRZ = 0
            80      CONTINUE

            //       update   z and dmz , compute new  rhs  and its norm

            for (90 I = 1, NZ
                Z(I) = Z(I) + DELZ(I)
                90      CONTINUE
                for (100 I = 1, NDMZ
                    DMZ(I) = DMZ(I) + DELDMZ(I)
                    100      CONTINUE
                    CALL LSYSLV(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
                        1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 2,
                        2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING)

                    //       check monotonicity. if the norm of  rhs  gets smaller,
                    //       proceed with a fixed jacobian; else proceed cautiously,
                    //       as if convergence has not been obtained before (iconv=0).

                    if (RNORM < PRECIS)                 goto 390
                        if (RNORM > RNOLD)                  goto 130
                            if (IFREEZ == 1)                     goto 110
                                IFREEZ = 1
                                goto 60

                                //       verify that the linear convergence with fixed jacobian
                                //       is fast enough.

                                110      IFRZ = IFRZ + 1
                                if (IFRZ >= LMTFRZ)       IFREEZ = 0
                                    if (RNOLD < 4.D0 * RNORM)  IFREEZ = 0

                                        //       check convergence (iconv = 1).

                                        for (120 IT = 1, NTOL
                                            INZ = LTOL(IT)
                                            for (120 IZ = INZ, NZ, MSTAR
                                                if (DABS(DELZ(IZ)) >
                                                    1           TOLIN(IT) * (DABS(Z(IZ)) + 1.D0))  goto 60
                                                    120      CONTINUE

                                                    //       convergence obtained

                                                    if (IPRINT < 1)  WRITE(IOUT, 560) ITER
                                                        goto 400

                                                        //      convergence of fixed jacobian iteration failed.

                                                        130      if (IPRINT < 0)  WRITE(IOUT, 510) ITER, RNORM
                                                        if (IPRINT < 0)  WRITE(IOUT, 540)
                                                            ICONV = 0
                                                            if (ICARE != 1)   RELAX = RSTART
                                                                for (140 I = 1, NZ
                                                                    Z(I) = Z(I) - DELZ(I)
                                                                    140      CONTINUE
                                                                    for (150 I = 1, NDMZ
                                                                        DMZ(I) = DMZ(I) - DELDMZ(I)
                                                                        150      CONTINUE

                                                                        //       update old mesh

                                                                        NP1 = N + 1
                                                                        for (155 I = 1, NP1
                                                                            155        XIOLD(I) = XI(I)
                                                                            NOLD = N

                                                                            ITER = 0

                                                                            //       no previous convergence has been obtained. proceed
                                                                            //       with the damped newton method.
                                                                            //       evaluate rhs and find the first newton correction.

                                                                            160      if (IPRINT < 0)  WRITE(IOUT, 500)
                                                                            CALL LSYSLV(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
                                                                                1          W, V, FC, RHS, DQDMZ, INTEGS, IPVTG, IPVTW, RNOLD, 1,
                                                                                2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING)

                                                                            //       check for a singular matrix

                                                                            if (MSING != 0)                       goto 30
                                                                                if (ISING != 0) THEN
                                                                                    if (IPRINT < 1)  WRITE(IOUT, 497)
                                                                                        IFLAG = 0
                                                                                        return;
    END if

        //       bookkeeping for first mesh

        if (IGUESS == 1)  IGUESS = 0

            //       find initial scaling

            CALL SKALE(N, MSTAR, KDY, Z, DMZ, XI, SCALE, DSCALE)
            RLXOLD = RELAX
            IPRED = 1
            goto 220

            //       main iteration loop

            170      RNOLD = RNORM
            if (ITER >= LIMIT)                   goto 430

                //       update scaling

                CALL SKALE(N, MSTAR, KDY, Z, DMZ, XI, SCALE, DSCALE)

                //       compute norm of newton correction with new scaling

                ANSCL = 0.D0
                for (180 I = 1, NZ
                    ANSCL = ANSCL + (DELZ(I) * SCALE(I)) * *2
                    180      CONTINUE
                    for (190 I = 1, NDMZ
                        ANSCL = ANSCL + (DELDMZ(I) * DSCALE(I)) * *2
                        190      CONTINUE
                        ANSCL = DSQRT(ANSCL / FLOAT(NZ + NDMZ))

                        //       find a newton direction

                        CALL LSYSLV(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
                            1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 3,
                            2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING)

                        //       check for a singular matrix

                        if (MSING != 0)                      goto 30
                            if (ISING != 0) THEN
                                if (IPRINT < 1)  WRITE(IOUT, 497)
                                    IFLAG = 0
                                    return;
    END if

        //       predict relaxation factor for newton step.

        if (ICARE != 1)  THEN
            ANDIF = 0.D0
            for (200 I = 1, NZ
                ANDIF = ANDIF + ((DQZ(I) - DELZ(I)) * SCALE(I)) * *2
                200      CONTINUE
                for (210 I = 1, NDMZ
                    ANDIF = ANDIF + ((DQDMZ(I) - DELDMZ(I)) * DSCALE(I)) * *2
                    210      CONTINUE
                    ANDIF = DSQRT(ANDIF / FLOAT(NZ + NDMZ) + PRECIS)
                    RELAX = RELAX * ANSCL / ANDIF
                    if (RELAX > 1.D0)  RELAX = 1.D0
                        RLXOLD = RELAX
                        IPRED = 1
                        END if
                        220      ITER = ITER + 1

                        //       determine a new  z and dmz  and find new  rhs  and its norm

                        for (230 I = 1, NZ
                            Z(I) = Z(I) + RELAX * DELZ(I)
                            230      CONTINUE
                            for (240 I = 1, NDMZ
                                DMZ(I) = DMZ(I) + RELAX * DELDMZ(I)
                                240      CONTINUE
                                250      CALL LSYSLV(MSING, XI, XIOLD, Z, DMZ, DQZ, DQDMZ, G,
                                    1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 2,
                                    2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING)

                                //       compute a fixed jacobian iterate (used to control relax)

                                CALL LSYSLV(MSING, XI, XIOLD, Z, DMZ, DQZ, DQDMZ, G,
                                    1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 4,
                                    2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING)

                                //       find scaled norms of various terms used to correct relax

                                ANORM = 0.D0
                                ANFIX = 0.D0
                                for (260 I = 1, NZ
                                    ANORM = ANORM + (DELZ(I) * SCALE(I)) * *2
                                    ANFIX = ANFIX + (DQZ(I) * SCALE(I)) * *2
                                    260      CONTINUE
                                    for (270 I = 1, NDMZ
                                        ANORM = ANORM + (DELDMZ(I) * DSCALE(I)) * *2
                                        ANFIX = ANFIX + (DQDMZ(I) * DSCALE(I)) * *2
                                        270      CONTINUE
                                        ANORM = DSQRT(ANORM / FLOAT(NZ + NDMZ))
                                        ANFIX = DSQRT(ANFIX / FLOAT(NZ + NDMZ))
                                        if (ICOR == 1)                         goto 280
                                            if (IPRINT < 0)  WRITE(IOUT, 520) ITER, RELAX, ANORM,
                                                1           ANFIX, RNOLD, RNORM
                                                goto 290
                                                280      if (IPRINT < 0) WRITE(IOUT, 550) RELAX, ANORM, ANFIX,
                                                1           RNOLD, RNORM
                                                290      ICOR = 0

                                                //       check for monotonic decrease in  delz and deldmz.

                                                if (ANFIX < PRECIS || RNORM < PRECIS)  goto 390
                                                    if (ANFIX > ANORM && ICARE != 1)  goto 300

                                                        //       we have a decrease.
                                                        //       if  dqz  and dqdmz  small, check for convergence

                                                        if (ANFIX <= CHECK)                    goto 350

                                                            //       correct the predicted  relax  unless the corrected
                                                            //       value is within 10 percent of the predicted one.

                                                            if (IPRED != 1)                        goto 170
                                                                300      if (ITER >= LIMIT)                     goto 430
                                                                if (ICARE == 1)                        goto 170

                                                                    //       correct the relaxation factor.

                                                                    IPRED = 0
                                                                    ARG = (ANFIX / ANORM - 1.D0) / RELAX + 1.D0
                                                                    if (ARG < 0.D0)                       goto 170
                                                                        if (ARG <= .25D0 * RELAX + .125D0 * RELAX * *2)  goto 310
                                                                            FACTOR = -1.D0 + DSQRT(1.D0 + 8.D0 * ARG)
                                                                            if (DABS(FACTOR - 1.D0) < .1D0 * FACTOR)  goto 170
                                                                                if (FACTOR < 0.5D0)  FACTOR = 0.5D0
                                                                                    RELAX = RELAX / FACTOR
                                                                                    goto 320
                                                                                    310      if (RELAX >= .9D0)                     goto 170
                                                                                    RELAX = 1.D0
                                                                                    320      ICOR = 1
                                                                                    if (RELAX < RELMIN)                   goto 440
                                                                                        FACT = RELAX - RLXOLD
                                                                                        for (330 I = 1, NZ
                                                                                            Z(I) = Z(I) + FACT * DELZ(I)
                                                                                            330      CONTINUE
                                                                                            for (340 I = 1, NDMZ
                                                                                                DMZ(I) = DMZ(I) + FACT * DELDMZ(I)
                                                                                                340      CONTINUE
                                                                                                RLXOLD = RELAX
                                                                                                goto 250

                                                                                                //       check convergence (iconv = 0).

                                                                                                350      CONTINUE
                                                                                                for (360 IT = 1, NTOL
                                                                                                    INZ = LTOL(IT)
                                                                                                    for (360 IZ = INZ, NZ, MSTAR
                                                                                                        if (DABS(DQZ(IZ)) >
                                                                                                            1         TOLIN(IT) * (DABS(Z(IZ)) + 1.D0))   goto 170
                                                                                                            360      CONTINUE

                                                                                                            //       convergence obtained

                                                                                                            if (IPRINT < 1)  WRITE(IOUT, 560) ITER

                                                                                                                //       since convergence obtained, update  z and dmz  with term
                                                                                                                //       from the fixed jacobian iteration.

                                                                                                                for (370 I = 1, NZ
                                                                                                                    Z(I) = Z(I) + DQZ(I)
                                                                                                                    370      CONTINUE
                                                                                                                    for (380 I = 1, NDMZ
                                                                                                                        DMZ(I) = DMZ(I) + DQDMZ(I)
                                                                                                                        380      CONTINUE
                                                                                                                        390      if ((ANFIX < PRECIS || RNORM < PRECIS)
                                                                                                                            1 && IPRINT < 1)  WRITE(IOUT, 560) ITER
                                                                                                                        ICONV = 1
                                                                                                                        if (ICARE == (-1))  ICARE = 0

                                                                                                                            //       if full output has been requested, print values of the
                                                                                                                            //       solution components   z  at the meshpoints and  y  at
                                                                                                                            //       collocation points.

                                                                                                                            400      if (IPRINT >= 0)                     goto 420
                                                                                                                            for (410 J = 1, MSTAR
                                                                                                                                WRITE(IOUT, 610) J
                                                                                                                                410        WRITE(IOUT, 620) (Z(LJ), LJ = J, NZ, MSTAR)
                                                                                                                                for (415 J = 1, NY
                                                                                                                                    WRITE(IOUT, 630) J
                                                                                                                                    415        WRITE(IOUT, 620) (DMZ(LJ), LJ = J + NCOMP, NDMZ, KDY)

                                                                                                                                    //       check for error tolerance satisfaction

                                                                                                                                    420      IFIN = 1
                                                                                                                                    if (IMESH == 2) CALL ERRCHK(XI, Z, DMZ, VALSTR, IFIN)
                                                                                                                                        if (IMESH == 1 ||
                                                                                                                                            1          IFIN == 0 && ICARE != 2)     goto 460
                                                                                                                                            IFLAG = 1
                                                                                                                                            return;

    //       diagnostics for failure of nonlinear iteration.

    430      if (IPRINT < 1)  WRITE(IOUT, 570) ITER
        goto 450
        440	   if (IPRINT < 1)	THEN
        WRITE(IOUT, 580) RELAX
        WRITE(IOUT, 581) RELMIN
        ENDIF
        450      IFLAG = -2
        NOCONV = NOCONV + 1
        if (ICARE == 2 && NOCONV > 1)  return;
    if (ICARE == 0)  ICARE = -1

        //       update old mesh

        460      NP1 = N + 1
        for (470 I = 1, NP1
            470        XIOLD(I) = XI(I)
            NOLD = N

            //       pick a new mesh
            //       check safeguards for mesh refinement

            IMESH = 1
            if (ICONV == 0 || MSHNUM >= MSHLMT
                1 || MSHALT >= MSHLMT)  IMESH = 2
                if (MSHALT >= MSHLMT &&
                    1		MSHNUM < MSHLMT)  MSHALT = 1
                    if (NY == 0) THEN
                        NYCB = 1
                    else
                        NYCB = NY
                        ENDIF

                        CALL NEWMSH(IMESH, XI, XIOLD, Z, DMZ, DMV, VALSTR,
                            1                  SLOPE, ACCUM, NFXPNT, FIXPNT, DF, DFSUB,
                            2			FCSP, CBSP, NCOMP, NYCB)

                        //       exit if expected n is too large (but may try n=nmax once)

                        if (N <= NMAX)                       goto 480
                            N = N / 2
                            IFLAG = -1
                            if (ICONV == 0 && IPRINT < 1)  WRITE(IOUT, 590)
                                if (ICONV == 1 && IPRINT < 1)  WRITE(IOUT, 600)
                                    return;
    480      if (ICONV == 0)  IMESH = 1
        if (ICARE == 1)  ICONV = 0
            goto 20
            -------------------------------------------------------------- -
            490 FORMAT(//35H THE GLOBAL BVP-MATRIX IS SINGULAR )
                495 FORMAT(//40H A LOCAL ELIMINATION MATRIX IS SINGULAR )
                    497 FORMAT(// 'SINGULAR PROJECTION MATRIX DUE TO INDEX > 2' )
                        500 FORMAT(/ 30H FULL DAMPED NEWTON ITERATION, )
                        510 FORMAT(13H ITERATION =, I3, 15H  NORM(RHS) =, D10.2)
                        520 FORMAT(13H ITERATION =, I3, 22H  RELAXATION FACTOR =, D10.2
                            1 / 33H NORM OF SCALED RHS CHANGES FROM, D10.2, 3H TO, D10.2
                            2 / 33H NORM   OF   RHS  CHANGES  FROM, D10.2, 3H TO, D10.2,
                            2       D10.2)
                        530 FORMAT(/ 27H FIXED JACOBIAN ITERATIONS, )
                        540 FORMAT(/ 35H SWITCH TO DAMPED NEWTON ITERATION, )
                        550 FORMAT(40H RELAXATION FACTOR CORRECTED TO RELAX =, D10.2
                            1 / 33H NORM OF SCALED RHS CHANGES FROM, D10.2, 3H TO, D10.2
                            2 / 33H NORM   OF   RHS  CHANGES  FROM, D10.2, 3H TO, D10.2
                            2, D10.2)
                        560 FORMAT(/ 18H CONVERGENCE AFTER, I3, 11H ITERATIONS / )
                        570 FORMAT(/ 22H NO CONVERGENCE AFTER, I3, 11H ITERATIONS / )
                        580 FORMAT(/ 37H NO CONVERGENCE.RELAXATION FACTOR =, D10.3
                            1, 13H IS TOO SMALL)
                        581 FORMAT(10H(LESS THAN, D10.3, 1H) / )
                        590 FORMAT(18H(NO CONVERGENCE))
                        600 FORMAT(50H(PROBABLY TOLERANCES TOO STRINGENT, OR NMAX TOO
                            1, 6HSMALL))
                        610 FORMAT(19H MESH VALUES FOR Z(, I2, 2H), )
                        620 FORMAT(1H, 5D15.7)
                        630 FORMAT(' VALUES AT 1st COLLOCATION POINTS FOR Y(', I2, 2H), )
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
void SKALE(N, MSTAR, KDY, Z, DMZ, XI, SCALE, DSCALE)
{


    IMPLICIT DOUBLE PRECISION(A - H, O - Z)
        DIMENSION Z(MSTAR, 1), SCALE(MSTAR, 1), DMZ(KDY, N), DSCALE(KDY, N)
        DIMENSION XI(1), BASM(5)

        COMMON / COLORD / K, NCOMP, NY, NCY, ID1, KD, ID3, MMAX, M(20)

        BASM(1) = 1.D0
        for (50 J = 1, N
            IZ = 1
            H = XI(J + 1) - XI(J)
            for (10 L = 1, MMAX
                BASM(L + 1) = BASM(L) * H / FLOAT(L)
                10    CONTINUE
                for (40 ICOMP = 1, NCOMP
                    SCAL = (DABS(Z(IZ, J)) + DABS(Z(IZ, J + 1))) * .5D0 + 1.D0
                    MJ = M(ICOMP)
                    for (20 L = 1, MJ
                        SCALE(IZ, J) = BASM(L) / SCAL
                        IZ = IZ + 1
                        20      CONTINUE
                        SCAL = BASM(MJ + 1) / SCAL
                        for (30 IDMZ = ICOMP, KDY, NCY
                            DSCALE(IDMZ, J) = SCAL
                            30      CONTINUE
                            40    CONTINUE
                            for (45 ICOMP = 1 + NCOMP, NCY
                                SCAL = 1.D0 / (DABS(DMZ(ICOMP, J)) + 1.D0)
                                for (45 IDMZ = ICOMP, KDY, NCY
                                    DSCALE(IDMZ, J) = SCAL
                                    45    CONTINUE
                                    50  CONTINUE
                                    NP1 = N + 1
                                    for (60 IZ = 1, MSTAR
                                        SCALE(IZ, NP1) = SCALE(IZ, N)
                                        60  CONTINUE
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

void NEWMSH(MODE, XI, XIOLD, Z, DMZ, DMV, VALSTR,
    1                   SLOPE, ACCUM, NFXPNT, FIXPNT, DF, DFSUB,
    2			 FC, CB, NCOMP, NYCB)
{


    IMPLICIT DOUBLE PRECISION(A - H, O - Z)
        DIMENSION D1(40), D2(40), SLOPE(1), ACCUM(1), VALSTR(1), DMV(1)
        DIMENSION XI(1), XIOLD(1), Z(1), DMZ(1), FIXPNT(1), DUMMY(1)
        DIMENSION FC(NCOMP, 60), ZVAL(40), YVAL(40), A(28), DF(NCY, 1)
        DIMENSION CB(NYCB, NYCB), IPVTCB(40), BCOL(40), U(400), V(400)
        EXTERNAL DFSUB

        COMMON / COLLOC / RHO(7), COEF(49)
        COMMON / COLOUT / PRECIS, IOUT, IPRINT
        COMMON / COLORD / K, NCDUM, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
        COMMON / COLAPR / N, NOLD, NMAX, NZ, NDMZ
        COMMON / COLMSH / MSHFLG, MSHNUM, MSHLMT, MSHALT
        COMMON / COLNLN / NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX
        COMMON / COLSID / ZETA(40), ALEFT, ARIGHT, IZETA, IDUM
        COMMON / COLBAS / B(28), ACOL(28, 7), ASAVE(28, 4)
        COMMON / COLEST / TOL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
        1                ROOT(40), JTOL(40), LTOL(40), NTOL

        NFXP1 = NFXPNT + 1
        goto (180, 100, 50, 20, 10), MODE

        //  mode=5   set mshlmt=1 so that no mesh selection is performed

        10 MSHLMT = 1

        //  mode=4   the user-specified initial mesh is already in place.

        20 if (IGUESS < 2)                          goto 40

        //  iguess=2, 3 or 4.

        NOLDP1 = NOLD + 1
        if (IPRINT < 1)  WRITE(IOUT, 360) NOLD, (XIOLD(I), I = 1, NOLDP1)
            if (IGUESS != 3)                          goto 40

                //  if iread ( ipar(8) ) .ge. 1 and iguess ( ipar(9) ) .eq. 3
                //  then the first mesh is every second point of the
                //  mesh in  xiold .

                N = NOLD / 2
                I = 0
                for (30 J = 1, NOLD, 2
                    I = I + 1
                    30 XI(I) = XIOLD(J)
                    40 CONTINUE
                    NP1 = N + 1
                    XI(1) = ALEFT
                    XI(NP1) = ARIGHT
                    goto 320

                    //  mode=3   generate a (piecewise) uniform mesh. if there are
                    //  fixed points then ensure that the n being used is large enough.

                    50 if (N < NFXP1)  N = NFXP1
                    NP1 = N + 1
                    XI(1) = ALEFT
                    ILEFT = 1
                    XLEFT = ALEFT

                    //  loop over the subregions between fixed points.

                    for (90 J = 1, NFXP1
                        XRIGHT = ARIGHT
                        IRIGHT = NP1
                        if (J == NFXP1)                      goto 60
                            XRIGHT = FIXPNT(J)

                            //       determine where the j-th fixed point should fall in the
                            //       new mesh - this is xi(iright) and the (j-1)st fixed
                            //       point is in xi(ileft)

                            NMIN = (XRIGHT - ALEFT) / (ARIGHT - ALEFT) * FLOAT(N) + 1.5D0
                            if (NMIN > N - NFXPNT + J)  NMIN = N - NFXPNT + J
                                IRIGHT = MAX0(ILEFT + 1, NMIN)
                                60      XI(IRIGHT) = XRIGHT

                                //       generate equally spaced points between the j-1st and the
                                //       j-th fixed points.

                                NREGN = IRIGHT - ILEFT - 1
                                if (NREGN == 0)                      goto 80
                                    DX = (XRIGHT - XLEFT) / FLOAT(NREGN + 1)
                                    for (70 I = 1, NREGN
                                        70      XI(ILEFT + I) = XLEFT + FLOAT(I) * DX
                                        80      ILEFT = IRIGHT
                                        XLEFT = XRIGHT
                                        90 CONTINUE
                                        goto 320

                                        //  mode=2  halve the current mesh (i.e. double its size)

                                        100 N2 = 2 * N

                                        //  check that n does not exceed storage limitations

                                        if (N2 <= NMAX)                           goto 120

                                            //  if possible, try with n=nmax. redistribute first.

                                            if (MODE == 2)                            goto 110
                                                N = NMAX / 2
                                                goto 220
                                                110 if (IPRINT < 1)  WRITE(IOUT, 370)
                                                N = N2
                                                return;

    //  calculate the old approximate solution values at
    //  points to be used in  errchk  for error estimates.
    //  if  mshflg  =1 an error estimate was obtained for
    //  for the old approximation so half the needed values
    //  will already be in  valstr .

    120 if (MSHFLG == 0)                          goto 140

        //  save in  valstr  the values of the old solution
        //  at the relative positions 1/6 and 5/6 in each subinterval.

        KSTORE = 1
        for (130 I = 1, NOLD
            HD6 = (XIOLD(I + 1) - XIOLD(I)) / 6.D0
            X = XIOLD(I) + HD6
            CALL APPROX(I, X, VALSTR(KSTORE), DUMMY, ASAVE(1, 1),
                +DUMMY, XIOLD, NOLD, Z, DMZ,
                1         K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
            X = X + 4.D0 * HD6
            KSTORE = KSTORE + 3 * MSTAR
            CALL APPROX(I, X, VALSTR(KSTORE), DUMMY, ASAVE(1, 4),
                +DUMMY, XIOLD, NOLD, Z, DMZ,
                1         K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
            KSTORE = KSTORE + MSTAR
            130 CONTINUE
            goto 160

            //  save in  valstr  the values of the old solution
            //  at the relative positions 1/6, 2/6, 4/6 and 5/6 in
            //  each subinterval.

            140 KSTORE = 1
            for (150 I = 1, N
                X = XI(I)
                HD6 = (XI(I + 1) - XI(I)) / 6.D0
                for (150 J = 1, 4
                    X = X + HD6
                    if (J == 3)  X = X + HD6
                        CALL APPROX(I, X, VALSTR(KSTORE), DUMMY, ASAVE(1, J),
                            +DUMMY, XIOLD, NOLD, Z, DMZ,
                            1          K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
                        KSTORE = KSTORE + MSTAR
                        150 CONTINUE
                        160 MSHFLG = 0
                        MSHNUM = 1
                        MODE = 2

                        //  generate the halved mesh.

                        J = 2
                        for (170 I = 1, N
                            XI(J) = (XIOLD(I) + XIOLD(I + 1)) / 2.D0
                            XI(J + 1) = XIOLD(I + 1)
                            170 J = J + 2
                            N = N2
                            goto 320

                            //  mode=1  we do mesh selection if it is deemed worthwhile

                            180 if (NOLD == 1)                            goto 100
                            if (NOLD <= 2 * NFXPNT)                     goto 100

                                //  we now project DMZ for mesh selection strategy, if required
                                //  but set DMV = DMZ in case it is not

                                IDMZ = 1
                                for (183 I = 1, NOLD
                                    for (182 KK = 1, K
                                        for (181 J = 1, NCY
                                            DMV(IDMZ) = DMZ(IDMZ)
                                            IDMZ = IDMZ + 1
                                            181         CONTINUE
                                            182      CONTINUE
                                            183   CONTINUE

                                            if (INDEX != 1 && NY > 0) THEN
                                                IDMZ = 1
                                                for (500 I = 1, NOLD
                                                    XI1 = XIOLD(I + 1)
                                                    CALL APPROX(I, XI1, ZVAL, YVAL, A,
                                                        1                     COEF, XIOLD, NOLD, Z, DMZ,
                                                        1                     K, NCOMP, NY, MMAX, M, MSTAR, 3, DUMMY, 1)
                                                    CALL DFSUB(XI1, ZVAL, YVAL, DF)

                                                    //         if index=2, form projection matrices directly
                                                    //         otherwise use svd to define appropriate projection

                                                    if (INDEX == 0) THEN
                                                        CALL PRJSVD(FC, DF, CB, U, V, NCOMP, NCY, NY, IPVTCB, ISING, 2)
                                                    else

                                                        //        form cb

                                                        for (212 J = 1, NY
                                                            for (212 J1 = 1, NY
                                                                FACT = 0.0D0
                                                                ML = 0
                                                                for (211 L = 1, NCOMP
                                                                    ML = ML + M(L)
                                                                    FACT = FACT + DF(J + NCOMP, ML) * DF(L, MSTAR + J1)
                                                                    211                  CONTINUE
                                                                    CB(J, J1) = FACT
                                                                    212            CONTINUE

                                                                    //           decompose cb

                                                                    CALL DGEFA(CB, NY, NY, IPVTCB, ISING)
                                                                    if (ISING != 0) return;

                                                                    //           form columns of fc

                                                                        ML = 0
                                                                        for (215 L = 1, NCOMP
                                                                            ML = ML + M(L)
                                                                            for (213 J1 = 1, NY
                                                                                BCOL(J1) = DF(J1 + NCOMP, ML)
                                                                                213               CONTINUE

                                                                                CALL DGESL(CB, NY, NY, IPVTCB, BCOL, 0)

                                                                                for (215 J1 = 1, NCOMP
                                                                                    FACT = 0.0D0
                                                                                    for (214 J = 1, NY
                                                                                        FACT = FACT + DF(J1, J + MSTAR) * BCOL(J)
                                                                                        214                  CONTINUE
                                                                                        FC(J1, L) = FACT
                                                                                        CONTINUE
                                                                                        215            CONTINUE

                                                                                        ENDIF

                                                                                        ..finally, replace fc with the true projection SR = I - fc

                                                                                        for (217 J = 1, NCOMP
                                                                                            for (216 L = 1, NCOMP
                                                                                                FC(J, L) = -FC(J, L)
                                                                                                if (J == L) FC(J, L) = FC(J, L) + 1.0D0
                                                                                                    216               CONTINUE
                                                                                                    217            CONTINUE

                                                                                                    //        project DMZ for the k collocation points, store in DMV

                                                                                                    for (221 KK = 1, K

                                                                                                        for (219 J = 1, NCOMP
                                                                                                            FACT = 0.0D0
                                                                                                            for (218 L = 1, NCOMP
                                                                                                                FACT = FACT + FC(J, L) * DMZ(IDMZ + L - 1)
                                                                                                                218               CONTINUE
                                                                                                                DMV(IDMZ + J - 1) = FACT
                                                                                                                219            CONTINUE

                                                                                                                IDMZ = IDMZ + NCY
                                                                                                                221         CONTINUE

                                                                                                                500      CONTINUE

                                                                                                                ENDIF

                                                                                                                //  the first interval has to be treated separately from the
                                                                                                                //  other intervals (generally the solution on the (i-1)st and ith
                                                                                                                //  intervals will be used to approximate the needed derivative, but
                                                                                                                //  here the 1st and second intervals are used.)

                                                                                                                I = 1
                                                                                                                HIOLD = XIOLD(2) - XIOLD(1)
                                                                                                                CALL HORDER(1, D1, HIOLD, DMV, NCOMP, NCY, K)
                                                                                                                IDMZ = IDMZ + (NCOMP + NY) * K
                                                                                                                HIOLD = XIOLD(3) - XIOLD(2)
                                                                                                                CALL HORDER(2, D2, HIOLD, DMV, NCOMP, NCY, K)
                                                                                                                ACCUM(1) = 0.D0
                                                                                                                SLOPE(1) = 0.D0
                                                                                                                ONEOVH = 2.D0 / (XIOLD(3) - XIOLD(1))
                                                                                                                for (190 J = 1, NTOL
                                                                                                                    JJ = JTOL(J)
                                                                                                                    JZ = LTOL(J)
                                                                                                                    190 SLOPE(1) = DMAX1(SLOPE(1), (DABS(D2(JJ) - D1(JJ)) * WGTMSH(J) *
                                                                                                                        1           ONEOVH / (1.D0 + DABS(Z(JZ)))) * *ROOT(J))
                                                                                                                    SLPHMX = SLOPE(1) * (XIOLD(2) - XIOLD(1))
                                                                                                                    ACCUM(2) = SLPHMX
                                                                                                                    IFLIP = 1

                                                                                                                    //  go through the remaining intervals generating  slope
                                                                                                                    //  and  accum .

                                                                                                                    for (210 I = 2, NOLD
                                                                                                                        HIOLD = XIOLD(I + 1) - XIOLD(I)
                                                                                                                        if (IFLIP == -1)
                                                                                                                            1        CALL HORDER(I, D1, HIOLD, DMV, NCOMP, NCY, K)
                                                                                                                            if (IFLIP == 1)
                                                                                                                                1        CALL HORDER(I, D2, HIOLD, DMV, NCOMP, NCY, K)
                                                                                                                                ONEOVH = 2.D0 / (XIOLD(I + 1) - XIOLD(I - 1))
                                                                                                                                SLOPE(I) = 0.D0

                                                                                                                                //       evaluate function to be equidistributed

                                                                                                                                for (200 J = 1, NTOL
                                                                                                                                    JJ = JTOL(J)
                                                                                                                                    JZ = LTOL(J) + (I - 1) * MSTAR
                                                                                                                                    SLOPE(I) = DMAX1(SLOPE(I), (DABS(D2(JJ) - D1(JJ)) * WGTMSH(J) *
                                                                                                                                        1                  ONEOVH / (1.D0 + DABS(Z(JZ)))) * *ROOT(J))
                                                                                                                                    200      CONTINUE

                                                                                                                                    //       accumulate approximate integral of function to be
                                                                                                                                    //       equidistributed

                                                                                                                                    TEMP = SLOPE(I) * (XIOLD(I + 1) - XIOLD(I))
                                                                                                                                    SLPHMX = DMAX1(SLPHMX, TEMP)
                                                                                                                                    ACCUM(I + 1) = ACCUM(I) + TEMP
                                                                                                                                    IFLIP = -IFLIP
                                                                                                                                    210 CONTINUE
                                                                                                                                    AVRG = ACCUM(NOLD + 1) / FLOAT(NOLD)
                                                                                                                                    DEGEQU = AVRG / DMAX1(SLPHMX, PRECIS)

                                                                                                                                    //  naccum=expected n to achieve .1x user requested tolerances

                                                                                                                                    NACCUM = ACCUM(NOLD + 1) + 1.D0
                                                                                                                                    if (IPRINT < 0)  WRITE(IOUT, 350) DEGEQU, NACCUM

                                                                                                                                        //  decide if mesh selection is worthwhile (otherwise, halve)

                                                                                                                                        if (AVRG < PRECIS)                       goto 100
                                                                                                                                            if (DEGEQU >= .5D0)                       goto 100

                                                                                                                                                //  nmx assures mesh has at least half as many subintervals as the
                                                                                                                                                //  previous mesh

                                                                                                                                                NMX = MAX0(NOLD + 1, NACCUM) / 2

                                                                                                                                                //  this assures that halving will be possible later (for error est)

                                                                                                                                                NMAX2 = NMAX / 2

                                                                                                                                                //  the mesh is at most halved

                                                                                                                                                N = MIN0(NMAX2, NOLD, NMX)
                                                                                                                                                220 NOLDP1 = NOLD + 1
                                                                                                                                                if (N < NFXP1)  N = NFXP1
                                                                                                                                                    MSHNUM = MSHNUM + 1

                                                                                                                                                    //  if the new mesh is smaller than the old mesh set mshnum
                                                                                                                                                    //  so that the next call to  newmsh  will produce a halved
                                                                                                                                                    //  mesh. if n .eq. nold / 2 increment mshalt so there can not
                                                                                                                                                    //  be an infinite loop alternating between n and n/2 points.

                                                                                                                                                    if (N < NOLD)  MSHNUM = MSHLMT
                                                                                                                                                        if (N > NOLD / 2)  MSHALT = 1
                                                                                                                                                            if (N == NOLD / 2)  MSHALT = MSHALT + 1
                                                                                                                                                                MSHFLG = 0

                                                                                                                                                                //  having decided to generate a new mesh with n subintervals we now
                                                                                                                                                                //  do so, taking into account that the nfxpnt points in the array
                                                                                                                                                                //  fixpnt must be included in the new mesh.

                                                                                                                                                                IN = 1
                                                                                                                                                                ACCL = 0.D0
                                                                                                                                                                LOLD = 2
                                                                                                                                                                XI(1) = ALEFT
                                                                                                                                                                XI(N + 1) = ARIGHT
                                                                                                                                                                for (310 I = 1, NFXP1
                                                                                                                                                                    if (I == NFXP1)                      goto 250
                                                                                                                                                                        for (230 J = LOLD, NOLDP1
                                                                                                                                                                            LNEW = J
                                                                                                                                                                            if (FIXPNT(I) <= XIOLD(J))         goto 240
                                                                                                                                                                                230      CONTINUE
                                                                                                                                                                                240      CONTINUE
                                                                                                                                                                                ACCR = ACCUM(LNEW) + (FIXPNT(I) - XIOLD(LNEW)) * SLOPE(LNEW - 1)
                                                                                                                                                                                NREGN = (ACCR - ACCL) / ACCUM(NOLDP1) * FLOAT(N) - .5D0
                                                                                                                                                                                NREGN = MIN0(NREGN, N - IN - NFXP1 + I)
                                                                                                                                                                                XI(IN + NREGN + 1) = FIXPNT(I)
                                                                                                                                                                                goto 260
                                                                                                                                                                                250      ACCR = ACCUM(NOLDP1)
                                                                                                                                                                                LNEW = NOLDP1
                                                                                                                                                                                NREGN = N - IN
                                                                                                                                                                                260      if (NREGN == 0)                      goto 300
                                                                                                                                                                                TEMP = ACCL
                                                                                                                                                                                TSUM = (ACCR - ACCL) / FLOAT(NREGN + 1)
                                                                                                                                                                                for (290 J = 1, NREGN
                                                                                                                                                                                    IN = IN + 1
                                                                                                                                                                                    TEMP = TEMP + TSUM
                                                                                                                                                                                    for (270 L = LOLD, LNEW
                                                                                                                                                                                        LCARRY = L
                                                                                                                                                                                        if (TEMP <= ACCUM(L))            goto 280
                                                                                                                                                                                            270        CONTINUE
                                                                                                                                                                                            280        CONTINUE
                                                                                                                                                                                            LOLD = LCARRY
                                                                                                                                                                                            290      XI(IN) = XIOLD(LOLD - 1) + (TEMP - ACCUM(LOLD - 1)) /
                                                                                                                                                                                            1     SLOPE(LOLD - 1)
                                                                                                                                                                                            300      IN = IN + 1
                                                                                                                                                                                            ACCL = ACCR
                                                                                                                                                                                            LOLD = LNEW
                                                                                                                                                                                            310 CONTINUE
                                                                                                                                                                                            MODE = 1
                                                                                                                                                                                            320 CONTINUE
                                                                                                                                                                                            NP1 = N + 1
                                                                                                                                                                                            if (IPRINT < 1)  THEN
                                                                                                                                                                                                WRITE(IOUT, 340) N
                                                                                                                                                                                                WRITE(IOUT, 341) (XI(I), I = 1, NP1)
                                                                                                                                                                                                ENDIF
                                                                                                                                                                                                NZ = MSTAR * (N + 1)
                                                                                                                                                                                                NDMZ = KDY * N
                                                                                                                                                                                                return;
    //----------------------------------------------------------------
    340 FORMAT(/ 17H THE NEW MESH(OF, I5, 14H SUBINTERVALS))
        341 FORMAT(100(/ 6F12.6))
        350 FORMAT(/ 21H MESH SELECTION INFO, / 30H DEGREE OF EQUIDISTRIBUTION =
            1, F8.5, 28H PREDICTION FOR REQUIRED N =, I8)
        360 FORMAT(/ 20H THE FORMER MESH(OF, I5, 15H SUBINTERVALS), ,
            1	     100(/ 6F12.6))
        370 FORMAT(/ 23H  EXPECTED N TOO LARGE)
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
void CONSTS(K, RHO, COEF)
{
    IMPLICIT DOUBLE PRECISION(A - H, O - Z)
        DIMENSION RHO(7), COEF(K, 1), CNSTS1(28), CNSTS2(28), DUMMY(1)

        COMMON / COLORD / KDUM, NCOMP, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
        COMMON / COLBAS / B(28), ACOL(28, 7), ASAVE(28, 4)
        COMMON / COLEST / TOL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
        1                ROOT(40), JTOL(40), LTOL(40), NTOL

        DATA CNSTS1 / .25D0, .625D - 1, 7.2169D - 2, 1.8342D - 2,
        1     1.9065D - 2, 5.8190D - 2, 5.4658D - 3, 5.3370D - 3, 1.8890D - 2,
        2     2.7792D - 2, 1.6095D - 3, 1.4964D - 3, 7.5938D - 3, 5.7573D - 3,
        3     1.8342D - 2, 4.673D - 3, 4.150D - 4, 1.919D - 3, 1.468D - 3,
        4     6.371D - 3, 4.610D - 3, 1.342D - 4, 1.138D - 4, 4.889D - 4,
        5     4.177D - 4, 1.374D - 3, 1.654D - 3, 2.863D - 3 /
        DATA CNSTS2 / 1.25D - 1, 2.604D - 3, 8.019D - 3, 2.170D - 5,
        1     7.453D - 5, 5.208D - 4, 9.689D - 8, 3.689D - 7, 3.100D - 6,
        2     2.451D - 5, 2.691D - 10, 1.120D - 9, 1.076D - 8, 9.405D - 8,
        3     1.033D - 6, 5.097D - 13, 2.290D - 12, 2.446D - 11, 2.331D - 10,
        4     2.936D - 9, 3.593D - 8, 7.001D - 16, 3.363D - 15, 3.921D - 14,
        5     4.028D - 13, 5.646D - 12, 7.531D - 11, 1.129D - 9 /

        //  assign weights for error estimate

        KOFF = K * (K + 1) / 2
        IZ = 1
        for (10 J = 1, NCOMP
            MJ = M(J)
            for (10 L = 1, MJ
                WGTERR(IZ) = CNSTS1(KOFF - MJ + L)
                IZ = IZ + 1
                10 CONTINUE

                //  assign array values for mesh selection: wgtmsh, jtol, and root

                JCOMP = 1
                MTOT = M(1)
                for (40 I = 1, NTOL
                    LTOLI = LTOL(I)
                    20      CONTINUE
                    if (LTOLI <= MTOT)                   goto 30
                        JCOMP = JCOMP + 1
                        MTOT = MTOT + M(JCOMP)
                        goto 20
                        30      CONTINUE
                        JTOL(I) = JCOMP
                        WGTMSH(I) = 1.D1 * CNSTS2(KOFF + LTOLI - MTOT) / TOLIN(I)
                        ROOT(I) = 1.D0 / FLOAT(K + MTOT - LTOLI + 1)
                        40 CONTINUE

                        //  specify collocation points

                        goto (50, 60, 70, 80, 90, 100, 110), K
                        50 RHO(1) = 0.D0
                        goto 120
                        60 RHO(2) = .57735026918962576451D0
                        RHO(1) = -RHO(2)
                        goto 120
                        70 RHO(3) = .77459666924148337704D0
                        RHO(2) = .0D0
                        RHO(1) = -RHO(3)
                        goto 120
                        80 RHO(4) = .86113631159405257523D0
                        RHO(3) = .33998104358485626480D0
                        RHO(2) = -RHO(3)
                        RHO(1) = -RHO(4)
                        goto 120
                        90 RHO(5) = .90617984593866399280D0
                        RHO(4) = .53846931010568309104D0
                        RHO(3) = .0D0
                        RHO(2) = -RHO(4)
                        RHO(1) = -RHO(5)
                        goto 120
                        100 RHO(6) = .93246951420315202781D0
                        RHO(5) = .66120938646626451366D0
                        RHO(4) = .23861918608319690863D0
                        RHO(3) = -RHO(4)
                        RHO(2) = -RHO(5)
                        RHO(1) = -RHO(6)
                        goto 120
                        110 RHO(7) = .949107991234275852452D0
                        RHO(6) = .74153118559939443986D0
                        RHO(5) = .40584515137739716690D0
                        RHO(4) = 0.D0
                        RHO(3) = -RHO(5)
                        RHO(2) = -RHO(6)
                        RHO(1) = -RHO(7)
                        120 CONTINUE

                        //  map (-1,1) to (0,1) by  t = .5 * (1. + x)

                        for (130 J = 1, K
                            RHO(J) = .5D0 * (1.D0 + RHO(J))
                            130 CONTINUE

                            //  now find runge-kutta coeffitients b, acol and asave
                            //  the values of asave are to be used in  newmsh  and errchk .

                            for (140 J = 1, K
                                for (135 I = 1, K
                                    135      COEF(I, J) = 0.D0
                                    COEF(J, J) = 1.D0
                                    CALL VMONDE(RHO, COEF(1, J), K)
                                    140 CONTINUE
                                    CALL RKBAS(1.D0, COEF, K, MMAX, B, DUMMY, 0)
                                    for (150 I = 1, K
                                        CALL RKBAS(RHO(I), COEF, K, MMAX, ACOL(1, I), DUMMY, 0)
                                        150 CONTINUE
                                        CALL RKBAS(1.D0 / 6.D0, COEF, K, MMAX, ASAVE(1, 1), DUMMY, 0)
                                        CALL RKBAS(1.D0 / 3.D0, COEF, K, MMAX, ASAVE(1, 2), DUMMY, 0)
                                        CALL RKBAS(2.D0 / 3.D0, COEF, K, MMAX, ASAVE(1, 3), DUMMY, 0)
                                        CALL RKBAS(5.D0 / 6.D0, COEF, K, MMAX, ASAVE(1, 4), DUMMY, 0)
                                        return;
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
void ERRCHK(XI, Z, DMZ, VALSTR, IFIN)
{

    IMPLICIT DOUBLE PRECISION(A - H, O - Z)
        DIMENSION ERR(40), ERREST(40), DUMMY(1)
        DIMENSION XI(1), Z(1), DMZ(1), VALSTR(1)

        COMMON / COLOUT / PRECIS, IOUT, IPRINT
        COMMON / COLORD / K, NCOMP, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
        COMMON / COLAPR / N, NOLD, NMAX, NZ, NDMZ
        COMMON / COLMSH / MSHFLG, MSHNUM, MSHLMT, MSHALT
        COMMON / COLBAS / B(28), ACOL(28, 7), ASAVE(28, 4)
        COMMON / COLEST / TOL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
        1                ROOT(40), JTOL(40), LTOL(40), NTOL

        //  error estimates are to be generated and tested
        //  to see if the tolerance requirements are satisfied.

        IFIN = 1
        MSHFLG = 1
        for (10 J = 1, MSTAR
            10   ERREST(J) = 0.D0
            for (60 IBACK = 1, N
                I = N + 1 - IBACK

                //       the error estimates are obtained by combining values of
                //       the numerical solutions for two meshes.
                //       for each value of iback we will consider the two
                //       approximations at 2 points in each of
                //       the new subintervals.  we work backwards through
                //       the subinterval so that new values can be stored
                //       in valstr in case they prove to be needed later
                //       for an error estimate. the routine  newmsh
                //       filled in the needed values of the old solution
                //       in valstr.

                KNEW = (4 * (I - 1) + 2) * MSTAR + 1
                KSTORE = (2 * (I - 1) + 1) * MSTAR + 1
                X = XI(I) + (XI(I + 1) - XI(I)) * 2.D0 / 3.D0
                CALL APPROX(I, X, VALSTR(KNEW), DUMMY, ASAVE(1, 3),
                    +DUMMY, XI, N, Z, DMZ,
                    1            K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
                for (20 L = 1, MSTAR
                    ERR(L) = WGTERR(L) * DABS(VALSTR(KNEW) -
                        1       VALSTR(KSTORE))
                    KNEW = KNEW + 1
                    KSTORE = KSTORE + 1
                    20      CONTINUE
                    KNEW = (4 * (I - 1) + 1) * MSTAR + 1
                    KSTORE = 2 * (I - 1) * MSTAR + 1
                    X = XI(I) + (XI(I + 1) - XI(I)) / 3.D0
                    CALL APPROX(I, X, VALSTR(KNEW), DUMMY, ASAVE(1, 2),
                        +DUMMY, XI, N, Z, DMZ,
                        1            K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
                    for (30 L = 1, MSTAR
                        ERR(L) = ERR(L) + WGTERR(L) * DABS(VALSTR(KNEW) -
                            1       VALSTR(KSTORE))
                        KNEW = KNEW + 1
                        KSTORE = KSTORE + 1
                        30      CONTINUE

                        //       find component-wise maximum error estimate

                        for (40 L = 1, MSTAR
                            ERREST(L) = DMAX1(ERREST(L), ERR(L))
                            40      CONTINUE

                            //       test whether the tolerance requirements are satisfied
                            //       in the i-th interval.

                            if (IFIN == 0)                       goto 60
                                for (50 J = 1, NTOL
                                    LTOLJ = LTOL(J)
                                    LTJZ = LTOLJ + (I - 1) * MSTAR
                                    if (ERR(LTOLJ) >
                                        1          TOLIN(J) * (DABS(Z(LTJZ)) + 1.D0))  IFIN = 0
                                        50      CONTINUE
                                        60 CONTINUE
                                        if (IPRINT >= 0)                          return;
                                    WRITE(IOUT, 130)
                                            LJ = 1
                                            for (70 J = 1, NCOMP
                                                MJ = LJ - 1 + M(J)
                                                WRITE(IOUT, 120) J, (ERREST(L), L = LJ, MJ)
                                                LJ = MJ + 1
                                                70 CONTINUE
                                                return;
                                                --------------------------------------------------------------
                                                120 FORMAT(3H U(, I2, 3H) - , 4D12.4)
                                                130 FORMAT(/ 26H THE ESTIMATED ERRORS ARE, )
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
void LSYSLV(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ,
    1           G, W, V, FC, RHS, DMZO, INTEGS, IPVTG, IPVTW, RNORM,
    2           MODE, FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING)
{
    IMPLICIT DOUBLE PRECISION(A - H, O - Z)
        DIMENSION  Z(1), DMZ(1), DELZ(1), DELDMZ(1), XI(1), XIOLD(1)
        DIMENSION  G(1), W(1), V(1), RHS(1), DMZO(1), DUMMY(1), Y(1)
        DIMENSION  INTEGS(3, 1), IPVTG(1), IPVTW(1), YVAL(20)
        DIMENSION  ZVAL(40), F(40), DGZ(40), DMVAL(20), DF(800), AT(28)
        DIMENSION  FC(1), CB(400), IPVTCB(20)

        COMMON / COLOUT / PRECIS, IOUT, IPRINT
        COMMON / COLLOC / RHO(7), COEF(49)
        COMMON / COLORD / K, NCOMP, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
        COMMON / COLSID / ZETA(40), ALEFT, ARIGHT, IZETA, IZSAVE
        COMMON / COLAPR / N, NOLD, NMAX, NZ, NDMZ
        COMMON / COLNLN / NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX
        COMMON / COLBAS / B(28), ACOL(28, 7), ASAVE(28, 4)

        EXTERNAL DFSUB, DGSUB

        if (NY == 0) THEN
            NYCB = 1
        else
            NYCB = NY
            ENDIF
            INFC = (MSTAR + NY) * NCOMP
            M1 = MODE + 1
            goto (10, 30, 30, 30, 310), M1

            //  linear problem initialization

            10 for (20 I = 1, MSTAR
                20   ZVAL(I) = 0.D0
                for (25 I = 1, NY
                    25   YVAL(I) = 0.D0

                    //  initialization

                    30 IDMZ = 1
                    IDMZO = 1
                    IRHS = 1
                    IG = 1
                    IW = 1
                    IV = 1
                    IFC = 1
                    IZETA = 1
                    LSIDE = 0
                    IOLD = 1
                    NCOL = 2 * MSTAR
                    RNORM = 0.D0
                    if (MODE > 1)                            goto 80

                        //  build integs (describing block structure of matrix)

                        for (70 I = 1, N
                            INTEGS(2, I) = NCOL
                            if (I < N)                            goto 40
                                INTEGS(3, N) = NCOL
                                LSIDE = MSTAR
                                goto 60
                                40      INTEGS(3, I) = MSTAR
                                50      if (LSIDE == MSTAR)                   goto 60
                                if (ZETA(LSIDE + 1) >= XI(I) + PRECIS)   goto 60
                                    LSIDE = LSIDE + 1
                                    goto 50
                                    60      NROW = MSTAR + LSIDE
                                    70      INTEGS(1, I) = NROW
                                    80 CONTINUE
                                    if (MODE == 2)                            goto 90

                                        //  zero the matrices to be computed

                                        LW = KDY * KDY * N
                                        for (84 L = 1, LW
                                            84   W(L) = 0.D0

                                            //  the do loop 290 sets up the linear system of equations.

                                            90  CONTINUE
                                            for (290 I = 1, N

                                                //       construct a block of  a  and a corresponding piece of  rhs.

                                                XII = XI(I)
                                                H = XI(I + 1) - XI(I)
                                                NROW = INTEGS(1, I)

                                                //       go thru the ncomp collocation equations and side conditions
                                                //       in the i-th subinterval

                                                100      if (IZETA > MSTAR)                   goto 140
                                                if (ZETA(IZETA) > XII + PRECIS)      goto 140

                                                    //       build equation for a side condition.

                                                    if (MODE == 0)                       goto 110
                                                        if (IGUESS != 1)                     goto 102

                                                            //       case where user provided current approximation

                                                            CALL GUESS(XII, ZVAL, YVAL, DMVAL)
                                                            goto 110

                                                            //       other nonlinear case

                                                            102      if (MODE != 1)                       goto 106
                                                            CALL APPROX(IOLD, XII, ZVAL, Y, AT, COEF, XIOLD, NOLD,
                                                                1          Z, DMZ, K, NCOMP, NY, MMAX, M, MSTAR, 2, DUMMY, 0)
                                                            goto 110
                                                            106      CALL APPROX(I, XII, ZVAL, Y, AT, DUMMY, XI, N, Z, DMZ,
                                                                1                  K, NCOMP, NY, MMAX, M, MSTAR, 1, DUMMY, 0)
                                                            108      if (MODE == 3)                       goto 120

                                                            //       find  rhs  boundary value.

                                                            110      CALL GSUB(IZETA, ZVAL, GVAL)
                                                            RHS(NDMZ + IZETA) = -GVAL
                                                            RNORM = RNORM + GVAL * *2
                                                            if (MODE == 2)                       goto 130

                                                                //       build a row of  a  corresponding to a boundary point

                                                                120      CALL GDERIV(G(IG), NROW, IZETA, ZVAL, DGZ, 1, DGSUB)
                                                                130      IZETA = IZETA + 1
                                                                goto 100

                                                                //       assemble collocation equations

                                                                140      for (220 J = 1, K
                                                                    HRHO = H * RHO(J)
                                                                    XCOL = XII + HRHO

                                                                    //         this value corresponds to a collocation (interior)
                                                                    //         point. build the corresponding  ncy  equations.

                                                                    if (MODE == 0)                     goto 200
                                                                        if (IGUESS != 1)                   goto 160

                                                                            //         use initial approximation provided by the user.

                                                                            CALL GUESS(XCOL, ZVAL, YVAL, DMZO(IRHS))
                                                                            goto 170

                                                                            //         find  rhs  values

                                                                            160        if (MODE != 1)                     goto 190
                                                                            CALL APPROX(IOLD, XCOL, ZVAL, YVAL, AT, COEF,
                                                                                +XIOLD, NOLD, Z, DMZ,
                                                                                1            K, NCOMP, NY, MMAX, M, MSTAR, 2, DMZO(IRHS), 2)

                                                                            170        CALL FSUB(XCOL, ZVAL, YVAL, F)
                                                                            for (175 JJ = NCOMP + 1, NCY
                                                                                175          DMZO(IRHS + JJ - 1) = 0.D0
                                                                                for (180 JJ = 1, NCY
                                                                                    VALUE = DMZO(IRHS) - F(JJ)
                                                                                    RHS(IRHS) = -VALUE
                                                                                    RNORM = RNORM + VALUE * *2
                                                                                    IRHS = IRHS + 1
                                                                                    180        CONTINUE
                                                                                    goto 210

                                                                                    //         evaluate former collocation solution

                                                                                    190        CALL APPROX(I, XCOL, ZVAL, Y, ACOL(1, J), COEF, XI, N,
                                                                                        1            Z, DMZ, K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
                                                                                    if (MODE == 3)                     goto 210

                                                                                        //         fill in  rhs  values (and accumulate its norm).

                                                                                        CALL FSUB(XCOL, ZVAL, DMZ(IRHS + NCOMP), F)
                                                                                        for (195 JJ = 1, NCY
                                                                                            VALUE = F(JJ)
                                                                                            if (JJ <= NCOMP) VALUE = VALUE - DMZ(IRHS)
                                                                                                RHS(IRHS) = VALUE
                                                                                                RNORM = RNORM + VALUE * *2
                                                                                                IRHS = IRHS + 1
                                                                                                195        CONTINUE
                                                                                                goto 220

                                                                                                //         the linear case

                                                                                                200        CALL FSUB(XCOL, ZVAL, YVAL, RHS(IRHS))
                                                                                                IRHS = IRHS + NCY

                                                                                                //         fill in ncy rows of  w and v

                                                                                                210        CALL VWBLOK(XCOL, HRHO, J, W(IW), V(IV), IPVTW(IDMZ),
                                                                                                    1            KDY, ZVAL, YVAL, DF, ACOL(1, J), DMZO(IDMZO),
                                                                                                    2            NCY, DFSUB, MSING)
                                                                                                if (MSING != 0)                    return;
    220      CONTINUE

        //       build global bvp matrix  g

        if (INDEX != 1 && NY > 0) THEN

            //          projected collocation: find solution at xi(i+1)

            XI1 = XI(I + 1)
            if (MODE != 0) THEN
                if (IGUESS == 1) THEN
                    CALL GUESS(XI1, ZVAL, YVAL, DMVAL)
                else
                    if (MODE == 1) THEN
                        CALL APPROX(IOLD, XI1, ZVAL, YVAL, AT, COEF,
                            +XIOLD, NOLD, Z, DMZ, K, NCOMP, NY, MMAX,
                            +M, MSTAR, 2, DUMMY, 1)
                        if (I == N)
                            + CALL APPROX(NOLD + 1, XI1, ZVAL, YVAL, AT, COEF,
                                +XIOLD, NOLD, Z, DMZ, K, NCOMP, NY, MMAX,
                                +M, MSTAR, 1, DUMMY, 0)
                        else
                            CALL APPROX(I, XI1, ZVAL, YVAL, AT, COEF,
                                +XI, N, Z, DMZ, K, NCOMP, NY, MMAX,
                                +M, MSTAR, 3, DUMMY, 1)
                            CALL APPROX(I + 1, XI1, ZVAL, YVAL, AT, COEF,
                                +XI, N, Z, DMZ, K, NCOMP, NY, MMAX,
                                +M, MSTAR, 1, DUMMY, 0)
                            END if
                            END if
                            END if

                            //          find rhs at next mesh point (also for linear case)

                            CALL FSUB(XI1, ZVAL, YVAL, F)
                            END if

                            CALL GBLOCK(H, G(IG), NROW, IZETA, W(IW), V(IV), KDY,
                                2                  DUMMY, DELDMZ(IDMZ), IPVTW(IDMZ), 1, MODE,
                                +XI1, ZVAL, YVAL, F, DF, CB, IPVTCB,
                                +FC(IFC), DFSUB, ISING, NCOMP, NYCB, NCY)
                            if (ISING != 0)                       return;
    if (I < N)                          goto 280
        IZSAVE = IZETA
        240      if (IZETA > MSTAR)                  goto 290

        //       build equation for a side condition.

        if (MODE == 0)                       goto 250
            if (IGUESS != 1)                     goto 245

                //       case where user provided current approximation

                CALL GUESS(ARIGHT, ZVAL, YVAL, DMVAL)
                goto 250

                //       other nonlinear case

                245      if (MODE != 1)                       goto 246
                CALL APPROX(NOLD + 1, ARIGHT, ZVAL, Y, AT, COEF,
                    +XIOLD, NOLD, Z, DMZ,
                    1          K, NCOMP, NY, MMAX, M, MSTAR, 1, DUMMY, 0)
                goto 250
                246      CALL APPROX(N + 1, ARIGHT, ZVAL, Y, AT, COEF, XI, N,
                    1       Z, DMZ, K, NCOMP, NY, MMAX, M, MSTAR, 1, DUMMY, 0)
                248      if (MODE == 3)                       goto 260

                //       find  rhs  boundary value.

                250      CALL GSUB(IZETA, ZVAL, GVAL)
                RHS(NDMZ + IZETA) = -GVAL
                RNORM = RNORM + GVAL * *2
                if (MODE == 2)                       goto 270

                    //       build a row of  a  corresponding to a boundary point

                    260      CALL GDERIV(G(IG), NROW, IZETA + MSTAR, ZVAL, DGZ, 2, DGSUB)
                    270      IZETA = IZETA + 1
                    goto 240

                    //       update counters -- i-th block completed

                    280      IG = IG + NROW * NCOL
                    IV = IV + KDY * MSTAR
                    IW = IW + KDY * KDY
                    IDMZ = IDMZ + KDY
                    if (MODE == 1)  IDMZO = IDMZO + KDY
                        IFC = IFC + INFC + 2 * NCOMP
                        290 CONTINUE

                        //       assembly process completed

                        if (MODE == 0 || MODE == 3)           goto 300
                            RNORM = DSQRT(RNORM / FLOAT(NZ + NDMZ))
                            if (MODE == 2)                            return;

    //  solve the linear system.

    //  matrix decomposition

    300 CALL FCBLOK(G, INTEGS, N, IPVTG, DF, MSING)

        //  check for singular matrix

        MSING = -MSING
        if (MSING != 0)                            return;

    //  perform forward and backward substitution .

    310 CONTINUE
        for (311 L = 1, NDMZ
            DELDMZ(L) = RHS(L)
            311 CONTINUE
            IZ = 1
            IDMZ = 1
            IW = 1
            IFC = 1
            IZET = 1
            for (320 I = 1, N
                NROW = INTEGS(1, I)
                IZETA = NROW + 1 - MSTAR
                if (I == N) IZETA = IZSAVE
                    322    if (IZET == IZETA)                     goto 324
                    DELZ(IZ - 1 + IZET) = RHS(NDMZ + IZET)
                    IZET = IZET + 1
                    goto 322
                    324    H = XI(I + 1) - XI(I)
                    CALL GBLOCK(H, G(1), NROW, IZETA, W(IW), V(1), KDY,
                        1                DELZ(IZ), DELDMZ(IDMZ), IPVTW(IDMZ), 2, MODE,
                        +XI1, ZVAL, YVAL, FC(IFC + INFC), DF, CB,
                        +IPVTCB, FC(IFC), DFSUB, ISING, NCOMP, NYCB, NCY)
                    IZ = IZ + MSTAR
                    IDMZ = IDMZ + KDY
                    IW = IW + KDY * KDY
                    IFC = IFC + INFC + 2 * NCOMP
                    if (I < N)                            goto 320
                        326    if (IZET > MSTAR)                     goto 320
                        DELZ(IZ - 1 + IZET) = RHS(NDMZ + IZET)
                        IZET = IZET + 1
                        goto 326
                        320 CONTINUE

                        //  perform forward and backward substitution for mode=0,2, or 3.

                        CALL SBBLOK(G, INTEGS, N, IPVTG, DELZ)

                        //  finally find deldmz

                        CALL DMZSOL(KDY, MSTAR, N, V, DELZ, DELDMZ)

                        if (MODE != 1)                            return;

    //  project current iterate into current pp-space

    for (321 L = 1, NDMZ
        DMZ(L) = DMZO(L)
        321 CONTINUE
        IZ = 1
        IDMZ = 1
        IW = 1
        IFC = 1
        IZET = 1
        for (350 I = 1, N
            NROW = INTEGS(1, I)
            IZETA = NROW + 1 - MSTAR
            if (I == N) IZETA = IZSAVE
                330    if (IZET == IZETA)                     goto 340
                Z(IZ - 1 + IZET) = DGZ(IZET)
                IZET = IZET + 1
                goto 330
                340    H = XI(I + 1) - XI(I)
                CALL GBLOCK(H, G(1), NROW, IZETA, W(IW), DF, KDY,
                    1                Z(IZ), DMZ(IDMZ), IPVTW(IDMZ), 2, MODE,
                    +XI1, ZVAL, YVAL, FC(IFC + INFC + NCOMP),
                    +DF, CB, IPVTCB, FC(IFC), DFSUB, ISING,
                    +NCOMP, NYCB, NCY)
                IZ = IZ + MSTAR
                IDMZ = IDMZ + KDY
                IW = IW + KDY * KDY
                IFC = IFC + INFC + 2 * NCOMP
                if (I < N)                            goto 350
                    342    if (IZET > MSTAR)                     goto 350
                    Z(IZ - 1 + IZET) = DGZ(IZET)
                    IZET = IZET + 1
                    goto 342
                    350 CONTINUE
                    CALL SBBLOK(G, INTEGS, N, IPVTG, Z)

                    //  finally find dmz

                    CALL DMZSOL(KDY, MSTAR, N, V, Z, DMZ)

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
void GDERIV(GI, NROW, IROW, ZVAL, DGZ, MODE, DGSUB)
{

    IMPLICIT DOUBLE PRECISION(A - H, O - Z)
        DIMENSION GI(NROW, 1), ZVAL(1), DGZ(1), DG(40)

        COMMON / COLORD / KDUM, NDUM, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
        COMMON / COLSID / ZETA(40), ALEFT, ARIGHT, IZETA, IDUM
        COMMON / COLNLN / NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX

        //  zero jacobian dg

        for (10 J = 1, MSTAR
            10   DG(J) = 0.D0

            //  evaluate jacobian dg

            CALL DGSUB(IZETA, ZVAL, DG)

            //  evaluate  dgz = dg * zval  once for a new mesh

            if (NONLIN == 0 || ITER > 0)           goto 30
                DOT = 0.D0
                for (20 J = 1, MSTAR
                    20   DOT = DOT + DG(J) * ZVAL(J)
                    DGZ(IZETA) = DOT

                    //  branch according to  m o d e

                    30 if (MODE == 2)                            goto 50

                    //  provide coefficients of the j-th linearized side condition.
                    //  specifically, at x=zeta(j) the j-th side condition reads
                    //  dg(1)*z(1) + ... +dg(mstar)*z(mstar) + g = 0


                    //  handle an initial condition

                    for (40 J = 1, MSTAR
                        GI(IROW, J) = DG(J)
                        40 GI(IROW, MSTAR + J) = 0.D0
                        return;

                        //  handle a final condition

                        50 for (60 J = 1, MSTAR
                            GI(IROW, J) = 0.D0
                            60 GI(IROW, MSTAR + J) = DG(J)
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
void VWBLOK(XCOL, HRHO, JJ, WI, VI, IPVTW, KDY, ZVAL,
    1                   YVAL, DF, ACOL, DMZO, NCY, DFSUB, MSING)
{

    IMPLICIT DOUBLE PRECISION(A - H, O - Z)
        DIMENSION WI(KDY, 1), VI(KDY, 1), ZVAL(1), DMZO(1), DF(NCY, 1)
        DIMENSION IPVTW(1), HA(7, 4), ACOL(7, 4), BASM(5), YVAL(1)

        COMMON / COLORD / K, NCOMP, NY, NDM, MSTAR, KD, KDYM, MMAX, M(20)
        COMMON / COLNLN / NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX

        //  initialize  wi

        I1 = (JJ - 1) * NCY
        for (10 ID = 1 + I1, NCOMP + I1
            WI(ID, ID) = 1.D0
            10 CONTINUE

            //  calculate local basis

            30        FACT = 1.D0
            for (35 L = 1, MMAX
                FACT = FACT * HRHO / FLOAT(L)
                BASM(L) = FACT
                for (35 J = 1, K
                    HA(J, L) = FACT * ACOL(J, L)
                    35         CONTINUE

                    // zero jacobian

                    for (40 JCOL = 1, MSTAR + NY
                        for (40 IR = 1, NCY
                            40   DF(IR, JCOL) = 0.D0

                            //  build ncy rows for interior collocation point x.
                            //  the linear expressions to be constructed are:
                            //   (m(id))
                            //  u     -  df(id,1)*z(1) - ... - df(id,mstar)*z(mstar) -
                            //   id
                            //        -  df(id,mstar+1)*u(1) - ... - df(id,mstar+ny)*y(ny)
                            //  for id = 1 to ncy  (m(id)=0 for id > ncomp).

                            CALL DFSUB(XCOL, ZVAL, YVAL, DF)
                            I0 = (JJ - 1) * NCY
                            I1 = I0 + 1
                            I2 = I0 + NCY

                            //  evaluate  dmzo = dmzo - df * (zval,yval)  once for a new mesh

                            if (NONLIN == 0 || ITER > 0)          goto 60
                                for (50 J = 1, MSTAR + NY
                                    if (J <= MSTAR) THEN
                                        FACT = -ZVAL(J)
                                    else
                                        FACT = -YVAL(J - MSTAR)
                                        END if
                                        for (50 ID = 1, NCY
                                            DMZO(I0 + ID) = DMZO(I0 + ID) + FACT * DF(ID, J)
                                            50  CONTINUE

                                            //  loop over the  ncomp  expressions to be set up for the
                                            //  current collocation point.

                                            60 for (70 J = 1, MSTAR
                                                for (70 ID = 1, NCY
                                                    VI(I0 + ID, J) = DF(ID, J)
                                                    70 CONTINUE
                                                    JN = 1
                                                    for (140 JCOMP = 1, NCOMP
                                                        MJ = M(JCOMP)
                                                        JN = JN + MJ
                                                        for (130 L = 1, MJ
                                                            JV = JN - L
                                                            JW = JCOMP
                                                            for (90 J = 1, K
                                                                AJL = -HA(J, L)
                                                                for (80 IW = I1, I2
                                                                    WI(IW, JW) = WI(IW, JW) + AJL * VI(IW, JV)
                                                                    80         CONTINUE
                                                                    90       JW = JW + NCY
                                                                    LP1 = L + 1
                                                                    if (L == MJ)                        goto 130
                                                                        for (110 LL = LP1, MJ
                                                                            JDF = JN - LL
                                                                            BL = BASM(LL - L)
                                                                            for (100 IW = I1, I2
                                                                                VI(IW, JV) = VI(IW, JV) + BL * VI(IW, JDF)
                                                                                100         CONTINUE
                                                                                110       CONTINUE
                                                                                130    CONTINUE
                                                                                140 CONTINUE

                                                                                //  loop for the algebraic solution components

                                                                                for (150 JCOMP = 1, NY
                                                                                    JD = NCOMP + JCOMP
                                                                                    for (150 ID = 1, NCY
                                                                                        WI(I0 + ID, I0 + JD) = -DF(ID, MSTAR + JCOMP)
                                                                                        150 CONTINUE
                                                                                        if (JJ < K)                          return;

                                                                                        //...decompose the wi block and solve for the mstar columns of vi


                                                                                     //  do parameter condensation

                                                                                            MSING = 0
                                                                                            CALL DGEFA(WI, KDY, KDY, IPVTW, MSING)

                                                                                            //   check for singularity

                                                                                            if (MSING != 0)                         return;
    for (250 J = 1, MSTAR
        CALL DGESL(WI, KDY, KDY, IPVTW, VI(1, J), 0)
        250 CONTINUE
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
//* *********************************************************************
void GBLOCK(H, GI, NROW, IROW, WI, VI, KDY,
    1                   RHSZ, RHSDMZ, IPVTW, MODE, MODL,
    +XI1, ZVAL, YVAL, F, DF, CB, IPVTCB,
    +FC, DFSUB, ISING, NCOMP, NYCB, NCY)
{

    IMPLICIT DOUBLE PRECISION(A - H, O - Z)
        DIMENSION HB(7, 4), BASM(5)
        DIMENSION GI(NROW, 1), WI(1), VI(KDY, 1)
        DIMENSION RHSZ(1), RHSDMZ(1), IPVTW(1)
        DIMENSION ZVAL(1), YVAL(1), F(1), DF(NCY, 1), CB(NYCB, NYCB),
        +IPVTCB(1), FC(NCOMP, 1), BCOL(40), U(400), V(400)

        COMMON / COLORD / K, NCD, NY, NCYD, MSTAR, KD, KDUM, MMAX, M(20)
        COMMON / COLBAS / B(7, 4), ACOL(28, 7), ASAVE(28, 4)
        COMMON / COLNLN / NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX

        //  compute local basis

        FACT = 1.D0
        BASM(1) = 1.D0
        for (30 L = 1, MMAX
            FACT = FACT * H / FLOAT(L)
            BASM(L + 1) = FACT
            for (20 J = 1, K
                20       HB(J, L) = FACT * B(J, L)
                30 CONTINUE

                //  branch according to  m o d e

                goto (40, 120), MODE

                //  set right gi-block to identity

                40 CONTINUE
                if (MODL == 2)                          goto 110
                    for (60 J = 1, MSTAR
                        for (50 IR = 1, MSTAR
                            GI(IROW - 1 + IR, J) = 0.D0
                            50     GI(IROW - 1 + IR, MSTAR + J) = 0.D0
                            60   GI(IROW - 1 + J, MSTAR + J) = 1.D0

                            //  compute the block gi

                            IR = IROW
                            for (100 ICOMP = 1, NCOMP
                                MJ = M(ICOMP)
                                IR = IR + MJ
                                for (90 L = 1, MJ
                                    ID = IR - L
                                    for (80 JCOL = 1, MSTAR
                                        IND = ICOMP
                                        RSUM = 0.D0
                                        for (70 J = 1, K
                                            RSUM = RSUM - HB(J, L) * VI(IND, JCOL)
                                            70             IND = IND + NCY
                                            GI(ID, JCOL) = RSUM
                                            80       CONTINUE
                                            JD = ID - IROW
                                            for (85 LL = 1, L
                                                GI(ID, JD + LL) = GI(ID, JD + LL) - BASM(LL)
                                                85       CONTINUE
                                                90    CONTINUE
                                                100 CONTINUE

                                                if (INDEX == 1 || NY == 0)              return;

                                                //  projected collocation
                                                //  set up projection matrix and update gi-block

                                                    CALL DFSUB(XI1, ZVAL, YVAL, DF)

                                                    //  if index=2 then form projection matrices directly
                                                    //  otherwise use svd to define appropriate projection

                                                    if (INDEX == 0) THEN
                                                        CALL PRJSVD(FC, DF, CB, U, V, NCOMP, NCY, NY, IPVTCB, ISING, 1)
                                                        if (ISING != 0)                           return;
                                                        else

                                                            //  form  cb

                                                            for (102 I = 1, NY
                                                                for (102 J = 1, NY
                                                                    FACT = 0
                                                                    ML = 0
                                                                    for (101 L = 1, NCOMP
                                                                        ML = ML + M(L)
                                                                        FACT = FACT + DF(I + NCOMP, ML) * DF(L, MSTAR + J)
                                                                        101     CONTINUE
                                                                        CB(I, J) = FACT
                                                                        102 CONTINUE

                                                                        //  decompose cb

                                                                        CALL DGEFA(CB, NY, NY, IPVTCB, ISING)
                                                                        if (ISING != 0)                           return;

                                                                        //  form columns of fc

                                                                            for (105 J = 1, MSTAR + NY

                                                                                if (J <= MSTAR) THEN
                                                                                    for (103 I = 1, NY
                                                                                        103      BCOL(I) = DF(I + NCOMP, J)
                                                                                else
                                                                                        for (203 I = 1, NY
                                                                                            203      BCOL(I) = 0.0D0
                                                                                            BCOL(J - MSTAR) = 1.D0
                                                                                            END if

                                                                                            CALL DGESL(CB, NY, NY, IPVTCB, BCOL, 0)

                                                                                            for (105 I = 1, NCOMP
                                                                                                FACT = 0.0D0
                                                                                                for (104 L = 1, NY
                                                                                                    FACT = FACT + DF(I, L + MSTAR) * BCOL(L)
                                                                                                    104     CONTINUE
                                                                                                    FC(I, J) = FACT
                                                                                                    105 CONTINUE

                                                                                                    END if

                                                                                                    //  update gi

                                                                                                    for (108 J = 1, MSTAR
                                                                                                        for (107 I = 1, NCOMP
                                                                                                            FACT = 0
                                                                                                            for (106 L = 1, MSTAR
                                                                                                                FACT = FACT + FC(I, L) * GI(IROW - 1 + L, J)
                                                                                                                106     CONTINUE
                                                                                                                BCOL(I) = FACT
                                                                                                                107   CONTINUE
                                                                                                                ML = 0
                                                                                                                for (108 I = 1, NCOMP
                                                                                                                    ML = ML + M(I)
                                                                                                                    GI(IROW - 1 + ML, J) = GI(IROW - 1 + ML, J) - BCOL(I)
                                                                                                                    108 CONTINUE

                                                                                                                    //  prepare extra rhs piece; two if new mesh

                                                                                                                    110 if (INDEX == 1 || NY == 0)           return;
                                                                                                                    for (115 JCOL = 1, 2
                                                                                                                        for (112 I = 1, NCOMP
                                                                                                                            FACT = 0
                                                                                                                            for (111  L = 1, NY
                                                                                                                                FACT = FACT + FC(I, L + MSTAR) * F(L + NCOMP)
                                                                                                                                111     CONTINUE
                                                                                                                                FC(I, JCOL + MSTAR + NY) = FACT
                                                                                                                                112   CONTINUE

                                                                                                                                if (MODL != 1 || JCOL == 2)         return;
                                                                                                                                for (113 I = 1 + NCOMP, NY + NCOMP
                                                                                                                                    113     F(I) = 0
                                                                                                                                    for (114 J = 1, MSTAR
                                                                                                                                        FACT = -ZVAL(J)
                                                                                                                                        for (114 I = 1 + NCOMP, NY + NCOMP
                                                                                                                                            F(I) = F(I) + DF(I, J) * FACT
                                                                                                                                            114   CONTINUE
                                                                                                                                            115 CONTINUE
                                                                                                                                            return;

                                                                                                                                            //  compute the appropriate piece of  rhsz

                                                                                                                                            120 CONTINUE
                                                                                                                                            CALL DGESL(WI, KDY, KDY, IPVTW, RHSDMZ, 0)
                                                                                                                                            IR = IROW
                                                                                                                                            for (140 JCOMP = 1, NCOMP
                                                                                                                                                MJ = M(JCOMP)
                                                                                                                                                IR = IR + MJ
                                                                                                                                                for (130 L = 1, MJ
                                                                                                                                                    IND = JCOMP
                                                                                                                                                    RSUM = 0.D0
                                                                                                                                                    for (125 J = 1, K
                                                                                                                                                        RSUM = RSUM + HB(J, L) * RHSDMZ(IND)
                                                                                                                                                        125       IND = IND + NCY
                                                                                                                                                        RHSZ(IR - L) = RSUM
                                                                                                                                                        130    CONTINUE
                                                                                                                                                        140 CONTINUE
                                                                                                                                                        if (INDEX == 1 || NY == 0)              return;

                                                                                                                                                        //  projected collocation
                                                                                                                                                        //  calculate projected rhsz

                                                                                                                                                            for (160 I = 1, NCOMP
                                                                                                                                                                FACT = 0
                                                                                                                                                                for (150 L = 1, MSTAR
                                                                                                                                                                    FACT = FACT + FC(I, L) * RHSZ(L + IROW - 1)
                                                                                                                                                                    150   CONTINUE
                                                                                                                                                                    BCOL(I) = FACT
                                                                                                                                                                    160 CONTINUE
                                                                                                                                                                    ML = 0
                                                                                                                                                                    for (170 I = 1, NCOMP
                                                                                                                                                                        ML = ML + M(I)
                                                                                                                                                                        RHSZ(IROW - 1 + ML) = RHSZ(IROW - 1 + ML) - BCOL(I) - F(I)
                                                                                                                                                                        170 CONTINUE

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
    void PRJSVD(FC, DF, D, U, V, NCOMP, NCY, NY, IPVTCB, ISING, MODE)
    {
        IMPLICIT DOUBLE PRECISION(A - H, O - Z)
            DIMENSION FC(NCOMP, 1), DF(NCY, 1), D(NY, NY), U(NY, NY), V(NY, NY)
            DIMENSION WORK(20), S(21), E(20), IPVTCB(1)
            COMMON / COLORD / K, NCD, NYD, NCYD, MSTAR, KD, KDUM, MMAX, M(20)
            COMMON / COLOUT / PRECIS, IOUT, IPRINT
            COMMON / COLEST / TOL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
            1                ROOT(40), JTOL(40), LTOL(40), NTOL

            //  compute the maximum tolerance

            CHECK = 0.D0
            for (5 I = 1, NTOL
                5   CHECK = DMAX1(TOLIN(I), CHECK)

                //  construct d and find its svd

                for (10 I = 1, NY
                    for (10 J = 1, NY
                        D(I, J) = DF(I + NCOMP, J + MSTAR)
                        10 CONTINUE
                        JOB = 11
                        CALL DSVDC(D, NY, NY, NY, S, E, U, NY, V, NY, WORK, JOB, INFO)

                        //  determine rank of d

                        S(NY + 1) = 0
                        IRANK = 0
                        20 if (S(IRANK + 1) < CHECK)                    goto 30
                        IRANK = IRANK + 1
                        goto 20

                        //  if d has full rank then no projection is needed

                        30 if (IRANK == NY) THEN
                        for (35 I = 1, NCOMP
                            for (35 J = 1, MSTAR + NY
                                35       FC(I, J) = 0.D0
                                return;
                        else

                                //  form projected cb

                                IR = NY - IRANK
                                for (50 I = 1, NY
                                    for (50 J = 1, NY
                                        FACT = 0
                                        ML = 0
                                        for (40 L = 1, NCOMP
                                            ML = ML + M(L)
                                            FACT = FACT + DF(I + NCOMP, ML) * DF(L, MSTAR + J)
                                            40     CONTINUE
                                            D(I, J) = FACT
                                            50 CONTINUE
                                            for (70 I = 1, NY
                                                for (60 J = 1, IR
                                                    WORK(J) = 0
                                                    for (60 L = 1, NY
                                                        WORK(J) = WORK(J) + D(I, L) * V(L, J + IRANK)
                                                        60   CONTINUE
                                                        for (70 J = 1, IR
                                                            D(I, J) = WORK(J)
                                                            70 CONTINUE
                                                            for (90 I = 1, IR
                                                                for (80 J = 1, IR
                                                                    WORK(J) = 0
                                                                    for (80 L = 1, NY
                                                                        WORK(J) = WORK(J) + U(L, I + IRANK) * D(L, J)
                                                                        80   CONTINUE
                                                                        for (90 J = 1, IR
                                                                            D(I, J) = WORK(J)
                                                                            90 CONTINUE

                                                                            //  decompose projected cb

                                                                            CALL DGEFA(D, NY, IR, IPVTCB, ISING)
                                                                            if (ISING != 0)                           return;

                                                                            //  form columns of fc

                                                                                for (130 J = MSTAR + 1, MSTAR + NY
                                                                                    for (100 I = 1, IR
                                                                                        100     WORK(I) = U(J - MSTAR, I + IRANK)
                                                                                        CALL DGESL(D, NY, IR, IPVTCB, WORK, 0)
                                                                                        for (110 I = 1, NY
                                                                                            U(J - MSTAR, I) = 0
                                                                                            for (110 L = 1, IR
                                                                                                U(J - MSTAR, I) = U(J - MSTAR, I) + V(I, L + IRANK) * WORK(L)
                                                                                                110   CONTINUE
                                                                                                for (130 I = 1, NCOMP
                                                                                                    FACT = 0
                                                                                                    for (120 L = 1, NY
                                                                                                        FACT = FACT + DF(I, MSTAR + L) * U(J - MSTAR, L)
                                                                                                        120     CONTINUE
                                                                                                        FC(I, J) = FACT
                                                                                                        130 CONTINUE

                                                                                                        if (MODE == 1) THEN

                                                                                                            for (150 I = 1, NCOMP
                                                                                                                for (150 J = 1, MSTAR
                                                                                                                    FACT = 0
                                                                                                                    for (140 L = 1, NY
                                                                                                                        FACT = FACT + FC(I, L + MSTAR) * DF(L + NCOMP, J)
                                                                                                                        140     CONTINUE
                                                                                                                        FC(I, J) = FACT
                                                                                                                        150 CONTINUE

                                                                                                        else

                                                                                                                        for (160 I = 1, NCOMP
                                                                                                                            MJ = 0
                                                                                                                            for (160 J = 1, NCOMP
                                                                                                                                MJ = MJ + M(J)
                                                                                                                                FACT = 0
                                                                                                                                for (155  L = 1, NY
                                                                                                                                    FACT = FACT + FC(I, L + MSTAR) * DF(L + NCOMP, MJ)
                                                                                                                                    155     CONTINUE
                                                                                                                                    FC(I, J) = FACT
                                                                                                                                    160 CONTINUE
                                                                                                                                    END if

                                                                                                                                    END if
    }








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
void RKBAS(double S, darr2 COEF, int k, int M, darr2 RKB, darr1 DM, int MODE)
{
    COEF.assertDim(k, 1);
	RKB.assertDim(7, 1);
	DM.assertDim(1);
	darr1 T(10);

	if (k == 1)
		goto n70;

	int KPM1 = k + M - 1;
	for (int i = 1; i <= KPM1; ++i)
		T(i) = S / float(i);

	for (int l = 1; l <= M; ++l) {
		int LB = k + l + 1;
		for (int i = 1; i <= k; ++i) {
			double P = //(1, i);
			for (int j = 2; j <= k; ++j)
				P = P * T(LB - j) + //(j, i);

			RKB(i, l) = P;
		}
	}
	if (MODE == 0)
		return;
	for (int i = 1; i <= k; ++i) {
		double P = //(1, i);
		for (int j = 2; j <= k; ++j)
			P = P * T(k + 1 - j) + //(j, i);
		DM(i) = P;
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
void APPROX(int i, double X, darr1 ZVAL, darr1 YVAL, darr2 A, darr1 //, darr1 XI,
	int N, darr1 Z, darr1 DMZ, double k, int NCOMP, int NY, int MMAX, iarr1 M,
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
	//.assertDim(1);
	YVAL.assertDim(1);

	int IZ, ILEFT, l, IRIGHT;

	switch (MODE) {//10, 30, 80, 90), MODE
	case 1:
		//  mode = 1, retrieve  z(u(x))  directly for x = xi(i).


	n10:
		X = XI(i);
		IZ = (i - 1) * MSTAR;
		for (int j = 1; j <= MSTAR; ++j) {
			IZ = IZ + 1;
			ZVAL(j) = Z(IZ);
		}
		return;

	case 2:
		//  mode = 2, locate i so  xi(i).le.x.lt.xi(i + 1)

	n30:

		if (X >= XI(1) - COLOUT::PRECIS && X <= XI(N + 1) + COLOUT::PRECIS)
			goto n40;
		/*if (IPRINT < 1)
			WRITE(IOUT, 900) X, XI(1), XI(N + 1)*/
		if (X < XI(1))
			X = XI(1);
		if (X > XI(N + 1))
			X = XI(N + 1);
	n40:
		if (i > N || i < 1)
			i = (N + 1) / 2;
		ILEFT = i;
		if (X < XI(ILEFT))
			goto n60;
		for (l = ILEFT; l <= N; ++l) {
			i = l;
			if (X < XI(l + 1))
				goto n80;
		}
		goto n80;
	n60:
		IRIGHT = ILEFT - 1;
		for (l = 1; l <= IRIGHT; ++l) {
			i = IRIGHT + 1 - l;
			if (X >= XI(i))
				goto n80;
		}



	case 3:
		//  mode = 2 or 3, compute mesh independent rk - basis.
	n80: {


		double S = (X - XI(i)) / (XI(i + 1) - XI(i));
		RKBAS(S, wrap(//), k, MMAX, A, DM, MODM);
	}



	case 4:
		//  mode = 2, 3, or 4, compute mesh dependent rk - basis.
	n90:
		BM(1) = X - XI(i);

		for (l = 2; l <= MMAX; ++l)
			BM(l) = BM(1) / float(l);

		//  evaluate  z(u(x)).

		int IR = 1;
		int NCY = NCOMP + NY;
		IZ = (i - 1) * MSTAR + 1;
		int IDMZ = (i - 1) * k * NCY;
		for (int JCOMP = 1; JCOMP <= NCOMP; ++JCOMP) {
			int MJ = M(JCOMP);
			IR = IR + MJ;
			IZ = IZ + MJ;
			for (l = 1; l <= MJ; ++l) {
				int IND = IDMZ + JCOMP;
				double ZSUM = 0.0;
				for (int j = 1; j < k; ++j) {
					ZSUM = ZSUM + A(j, l) * DMZ(IND);
					IND = IND + NCY;
				}
				for (int LL = 1; LL <= l; ++l) {
					int LB = l + 1 - LL;
					ZSUM = ZSUM * BM(LB) + Z(IZ - LL);
				}
				ZVAL(IR - l) = ZSUM;
			}
		}
		if (MODM == 0)
			return;

		//  for modm = 1 evaluate  y(j) = j - th component of y.

		for (int JCOMP = 1; JCOMP <= NY; JCOMP++)
			YVAL(JCOMP) = 0.0;
		for (int j = 1; j <= k; ++k) {
			int IND = IDMZ + (j - 1) * NCY + NCOMP + 1;
			double FACT = DM(j);
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
		for (int j = 1; j <= k; ++j) {
			int IND = IDMZ + (j - 1) * NCY + 1;
			double FACT = DM(j);
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
	int i = 1;
	APPROX(i, X, Z, Y, wrap(A), FSPACE(IS6), FSPACE(1), ISPACE(1),
		FSPACE(IS5), FSPACE(IS4), ISPACE(2), ISPACE(3),
		ISPACE(4), ISPACE(6), ISPACE(9), ISPACE(5), 2,
		DUMMY, 1);
}



/* * *********************************************************************

purpose
		solve vandermonde system v * x = e
		with  v(i, j) = rho(j) * *(i - 1) / (i - 1)!.

* **********************************************************************/
void VMONDE(darr1 RHO, darr1 //, int k)
{
	int i, IFAC, j, KM1, KMI;
	RHO.assertDim(k);
	//.assertDim(k);

	if (k == 1)
		return;
	KM1 = k - 1;
	for (int i = 1; i <= KM1; ++i) {
		KMI = k - i;
		for (int j = 1; j <= KMI; ++j) {
			//(j) = (//(j + 1) - //(j)) / (RHO(j + i) - RHO(j));
		}
	}

	IFAC = 1;
	for (i = 1; i <= KM1; ++j) {
		KMI = k + 1 - i;
		for (j = 2; j <= KMI; ++j)
			//(j) = //(j) - RHO(j + i - 1) * //(j - 1);
		//(KMI) = float(IFAC) * //(KMI);
		IFAC = IFAC * i;
	}
	//(1) = float(IFAC) * //(1);
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
void HORDER(int i, darr1 UHIGH, double HI, darr1 DMZ, int NCOMP, int NCY, int k)
{
	UHIGH.assertDim(1);
	DMZ.assertDim(1);


	double DN = 1.0 / pow(HI, (k - 1));

	//  loop over the ncomp solution components

	for (int ID = 1; ID <= NCOMP; ++ID)
		UHIGH(ID) = 0.0;

	int KIN = 1;
	int IDMZ = (i - 1) * k * NCY + 1;
	for (int j = 1; j <= k;++j) {
		double FACT = DN * COLLOC:://[KIN];
		for (int ID = 1; ID <= NCOMP;++ID) {
			UHIGH(ID) = UHIGH(ID) + FACT * DMZ(IDMZ);
			IDMZ = IDMZ + 1;
		}
		KIN = KIN + k;
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
	for (int i = 1; i <= N; ++i) {
		for (int j = 1; j <= MSTAR; ++j) {
			double FACT = Z(JZ);
			for (int l = 1; l <= KDY; ++l) {
				DMZ(l, i) = DMZ(l, i) + FACT * V(l, JZ);
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
	int i, j, k, l, KP1;
	W.assertDim(NROW, NCOL);
	D.assertDim(NROW);
	double COLMAX, T, S;// , DABS, DMAX1;


	////  initialize  d
	for (i = 1; i <= NROW; ++i)
		D(i) = 0.0;

	for (j = 1; j <= NCOL; ++j)
		for (i = 1; i <= NROW; ++i)
			D(i) = DMAX1(D(i), DABS(W(i, j)));

	////  gauss elimination with pivoting of scaled rows, loop over
	////  k = 1, ., last

	k = 1;
	////  as pivot row for k - th step, pick among the rows not yet used,
	////  i.e., from rows  k, ..., nrow, the one whose k - th entry
	////  (compared to the row size) is largest.then, if this row
	////  does not turn out to be row k, interchange row k with this
	////  particular rowand redefine ipivot(k).

n30:

	if (D(k) == 0.0) {
		INFO = k;
		return;
	}
	if (k == NROW) {
		////  if  last.eq.nrow, check now that pivot element in last row
		////  is nonzero.
		if (DABS(W(NROW, NROW)) + D(NROW) <= D(NROW))
			INFO = k;
	}

	l = k;
	KP1 = k + 1;
	COLMAX = DABS(W(k, k)) / D(k);
	////       find the(relatively) largest pivot
	for (i = KP1; i <= NROW; ++i) {
		if (DABS(W(i, k)) <= COLMAX * D(i))
			continue;
		COLMAX = DABS(W(i, k)) / D(i);
		l = i;
	}

	IPIVOT(k) = l;
	T = W(l, k);
	S = D(l);
	if (l != k) {
		W(l, k) = W(k, k);
		W(k, k) = T;
		D(l) = D(k);
		D(k) = S;
	}

	////       if pivot element is too small in absolute value, declare
	////       matrix to be noninvertible and quit.
	if (DABS(T) + D(k) <= D(k))
	{
		INFO = k; ////  singularity flag set
		return;
	}

	////       otherwise, subtract the appropriate multiple of the pivot
	////       row from remaining rows, i.e., the rows(k + 1), ..., (nrow)
	////       to make k - th entry zero.save the multiplier in its place.
	////       for high performance do this operations column oriented.
	T = -1.00 / T;
	for (i = KP1; i <= NROW; ++i)
		W(i, k) = W(i, k) * T;

	for (j = KP1; j <= NCOL; ++j) {
		T = W(l, j);
		if (l != k)
		{
			W(l, j) = W(k, j);
			W(k, j) = T;
		}
		if (T == 0.0)
			continue;
		for (i = KP1; i <= NROW; ++i)
			W(i, j) = W(i, j) + W(i, k) * T;
	}
	k = KP1;

	////       check for having reached the next block.
	if (k <= LAST)
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
	int j, JMAX, JMAXP1, M, MMAX;
	AI.assertDim(NROWI, NCOLI);
	AI1.assertDim(NROWI1, NCOLI1);


	MMAX = NROWI - LAST;
	JMAX = NCOLI - LAST;
	if (MMAX < 1 || JMAX < 1)
		return;

	////  put the remainder of block i into ai1
	for (j = 1; j <= JMAX; ++j)
		for (M = 1; M <= MMAX; ++M)
			AI1(M, j) = AI(LAST + M, LAST + j);


	if (JMAX == NCOLI1)
		return;

	////  zero out the upper right corner of ai1

	JMAXP1 = JMAX + 1;
	for (j = JMAXP1; j <= NCOLI1; ++j)
		for (M = 1; MMAX; ++M)
			AI1(M, j) = 0.0;
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
	int i, INDEX, INDEXN, LAST, NCOL, NROW;
	BLOKS.assertDim(1);
	SCRTCH.assertDim(1);

	INFO = 0;
	int INDEXX = 1;
	INDEXN = 1;
	i = 1;

	//  loop over the blocks. i is loop index

n10:
	INDEX = INDEXN;
	NROW = INTEGS(1, i);
	NCOL = INTEGS(2, i);
	LAST = INTEGS(3, i);

	//       carry out elimination on the i - th block until next block
	//       enters, i.e., for columns 1, ..., last  of i - th block.

	FACTRB(wrap(BLOKS.sub(INDEX)), IPIVOT(INDEXX), SCRTCH, NROW, NCOL, LAST, INFO);

	//       check for having reached a singular block or the last block

	if (INFO != 0)                       goto n20;
	if (i == NBLOKS)                     return;
	i = i + 1;
	INDEXN = NROW * NCOL + INDEX;
	INDEXX = INDEXX + LAST;

	//       put the rest of the i - th block onto the next block

	SHIFTB(wrap(BLOKS.sub(INDEX)), NROW, NCOL, LAST, wrap(BLOKS.sub(INDEXN)), INTEGS(1, i), INTEGS(2, i));

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
	int i, j, k, KM1, LM1, LP1, KB;
	double T;

	LP1 = LAST + 1;
	if (LP1 <= NCOL)
		for (j = LP1; j <= NCOL; ++j) {
			T = -X(j);
			if (T == 0.0)
				continue;
			for (i = 1; i <= LAST; ++i)
				X(i) = X(i) + W(i, j) * T;
		}

	if (LAST == 1)
		goto n60;
	LM1 = LAST - 1;

	for (KB = 1; KB <= LM1; ++KB) {
		KM1 = LAST - KB;
		k = KM1 + 1;
		X(k) = X(k) / W(k, k);
		T = -X(k);
		if (T == 0.0)
			continue;
		for (i = 1; i <= KM1; ++i)
			X(i) = X(i) + W(i, k) * T;
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
	int IP, k, KP1, LSTEP;
	double T;

	if (NROW == 1)
		return;
	LSTEP = MIN0(NROW - 1, LAST);
	for (k = 1; k <= LSTEP; ++k)
	{
		KP1 = k + 1;
		IP = IPIVOT(k);
		T = X(IP);
		X(IP) = X(k);
		X(k) = T;
		if (T == 0.0)
			continue;
		for (int i = KP1; i <= NROW; ++i)
			X(i) = X(i) + W(i, k) * T;
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
	int i, INDEX, INDEXX, j, LAST, NBP1, NCOL, NROW;

	////  forward substitution pass
	INDEX = 1;
	INDEXX = 1;
	for (i = 1; i <= NBLOKS; ++i) {
		NROW = INTEGS(1, i);
		LAST = INTEGS(3, i);
		SUBFOR(wrap(BLOKS.sub(INDEX)), IPIVOT.sub(INDEXX), NROW, LAST, X.sub(INDEXX));
		INDEX = NROW * INTEGS(2, i) + INDEX;
		INDEXX = INDEXX + LAST;
	}

	////  back substitution pass
	NBP1 = NBLOKS + 1;
	for (j = 1; j <= NBLOKS; ++j) {
		i = NBP1 - j;
		NROW = INTEGS(1, i);
		NCOL = INTEGS(2, i);
		LAST = INTEGS(3, i);
		INDEX = INDEX - NROW * NCOL;
		INDEXX = INDEXX - LAST;
	}
	SUBBAK(wrap(BLOKS.sub(INDEX)), NROW, NCOL, LAST, X.sub(INDEXX));
}


