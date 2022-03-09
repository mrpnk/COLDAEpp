
void COLDAE(systemParams const& params, options const& opts,
	ivec ispace, dvec fspace, output_t& iflag,
	fsub_t fsub, dfsub_t dfsub, gsub_t gsub, dgsub_t dgsub, guess_t guess)

	void APPSLN(double& X, dvec Z, dvec Y, dvec FSPACE, ivec ISPACE)


void CONTRL(dvec XI, dvec XIOLD, dvec Z, dvec DMZ, dvec DMV, dvec RHS, dvec DELZ, dvec DELDMZ,
	dvec DQZ, dvec DQDMZ, dvec G, dvec W, dvec V, dvec FC, dvec VALSTR, dvec SLOPE, dvec SCALE, dvec DSCALE,
	dvec ACCUM, ivec IPVTG, ivec INTEGS, ivec IPVTW, const int NFXPNT, dvec FIXPNT, output_t& iflag,
	fsub_t fsub, dfsub_t dfsub, gsub_t gsub, dgsub_t dgsub, guess_t guess)

void SKALE(dmat Z, dmat DMZ, dvec XI, dmat SCALE, dmat DSCALE)


void CONSTS()
void ERRCHK(dvec XI, dvec Z, dvec DMZ, dvec VALSTR, int& IFIN)

				

void LSYSLV(int& MSING, dvec XI, dvec XIOLD, dvec Z, dvec DMZ, dvec DELZ, dvec DELDMZ,
	dvec G, dvec W, dvec V, dvec FC, dvec RHS, dvec DMZO,
	imat INTEGS, ivec IPVTG, ivec IPVTW, double& RNORM, 
	const int MODE, fsub_t fsub, dfsub_t dfsub, gsub_t gsub,dgsub_t dgsub, guess_t guess, int& ISING)


void GDERIV(dmat GI, const int NROW, const int IROW, dvec ZVAL, dvec DGZ, const int MODE, dgsub_t dgsub)
void VWBLOK(const double XCOL, const double HRHO, const int JJ, dmat WI, dmat VI, ivec IPVTW,
dvec ZVAL, dvec YVAL, dmat DF, dmat acol, dvec DMZO, dfsub_t dfsub, int& MSING)


void PRJSVD(dmat FC, dmat DF, dmat D, dmat U, dmat V, ivec IPVTCB, int& ISING, const int MODE)


				
void RKBAS(const double S, const int k, const int M, dmat RKB, dvec DM, const int MODE)

																																								
void APPROX(int& i, double& X, dvec ZVAL, dvec YVAL, dmat A, dvec coef, dvec XI,
const int n, dvec Z, dvec DMZ, const int k, const int ncomp, const int ny, const int mmax, ivec M,
const int mstar, const int MODE, dvec DMVAL, const int MODM)

	
void VMONDE(dvec coef, int k)
	
void HORDER(const int i, dvec UHIGH, const double HI, dvec DMZ)
	
void DMZSOL(dmat V, dvec Z, dmat DMZ)

	

void FACTRB(dmat W, ivec IPIVOT, dvec D, const int NROW, const int NCOL, const int LAST, int& INFO)
void SHIFTB(dmat AI, const int NROWI, const int NCOLI, const int LAST, dmat AI1, const int NROWI1, const int NCOLI1)
void FCBLOK(dvec BLOKS, imat INTEGS, const int NBLOKS, ivec IPIVOT, dvec SCRTCH, int& INFO)
void SUBBAK(dmat W, const int NROW, const int NCOL, const int LAST, dvec X)
void SUBFOR(dmat W, ivec IPIVOT, const int NROW, const int LAST, dvec X)
void SBBLOK(dvec BLOKS, imat INTEGS, const int NBLOKS, ivec IPIVOT, dvec X)
	
		