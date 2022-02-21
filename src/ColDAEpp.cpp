#include "ColDAEpp.hpp"


int ncomp = 1;
int ny = 0;
iad1 orders = { 1 };
double left = 0, right = 1;
dad1 bcpoints = { 0 };
iad1 ltol = { 1 };
dad1 tol = { 0.1 };
dad1 fixpnt = {};
iad1 ispace(1000);
dad1 fspace(1000);

void fsub(double x, dar1 z, dar1 y, dar1 f) {
	f(1) = x;
}
void dfsub(double x, dar1 z, dar1 y, dar2 df) {
	df(1, 1) = 0;
}
void gsub(int i, dar1 z, double& g) {
	g = z(1) - 1;
}
void dgsub(int i, dar1 z, dar2 dg) {
	dg(1, 1) = 1;
}


int main()
{
	iad1 ipar = {0,0,0,1,1000,1000,0,0,0,0,0,0};

	int iflag;
	cda solver;
	solver.COLDAE(ncomp,ny,orders,left,right,bcpoints,
		ipar,ltol,tol,fixpnt,ispace,fspace,iflag,
		fsub,dfsub,gsub,dgsub,nullptr);

	std::cout << "iflag = " << iflag << std::endl;
	std::cin.get();
	return 0;
}
