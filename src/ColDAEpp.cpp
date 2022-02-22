#pragma warning( disable : 5045 )

#include "ColDAEpp.hpp"


int ncomp = 1;
int ny = 0;
iad1 orders = { 1 };
double left = 0, right = 1;
dad1 bcpoints = { 0 };
iad1 ltol = { 1 };
dad1 tol = { 0.1 };
dad1 fixpnt = {};
iad1 ispace(10000);
dad1 fspace(10000);

void fsub(double x, dar1 z, dar1 y, dar1 f) {
	f(1) = x;
}
void dfsub(double x, dar1 z, dar1 y, dar2 df) {
	df(1, 1) = 0;
}
void gsub(int i, dar1 z, double& g) {
	g = z(1) - 1;
}
void dgsub(int i, dar1 z, dar1 dg) {
	dg(1) = 1;
}


int main()
{
	iad1 ipar = {0,0,0,1,10000,10000,-1,0,0,0,0,0};
	
	int iflag;
	cda solver;
	solver.COLDAE(ncomp,ny,orders,left,right,bcpoints,
		ipar,ltol,tol,fixpnt,ispace,fspace,iflag,
		fsub,dfsub,gsub,dgsub,nullptr);

	
	fmt::print("iflag = {}", iflag);

	if (iflag == 1) {
		std::ofstream file("result_cpp.txt");
		for (int i = 1; i <= ispace(1) + 1; ++i) {
			double x = fspace(i);
			dad1 z(1), y(0);
			solver.APPSLN(x, z, y, fspace, ispace);
			file << x << " " << z(1) << std::endl;
		}
		file.close();
	}

	std::cin.get();
	return 0;
}
