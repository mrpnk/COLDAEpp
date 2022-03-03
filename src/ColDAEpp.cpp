#pragma warning( disable : 5045 )

#include <float.h>
unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);


#include "ColDAEpp.hpp"
#include "timer.h"



// f'(x) = x, f(0)=1
struct sys1 {
	int ncomp = 1;
	int ny = 0;
	iad1 orders = { 1 };
	double left = 0, right = 1;
	dad1 bcpoints = { 0 };
	iad1 ltol = { 1 };
	dad1 tol = { 0.1 };
	dad1 fixpnt = {};


	iad1 ipar = { 0,0,0,1,10000,10000,-1,0,0,0,0,0 };

	static void fsub(double x, dar1 z, dar1 y, dar1 f) {
		f(1) = x;
	}
	static void dfsub(double x, dar1 z, dar1 y, dar2 df) {
		df(1, 1) = 0;
	}
	static void gsub(int i, dar1 z, double& g) {
		g = z(1) - 1;
	}
	static void dgsub(int i, dar1 z, dar1 dg) {
		dg(1) = 1;
	}
};

// f''(x) = sin(0.6*f'(x)) + x, f(0)=1, f(1)=-0.1
struct sys2 {
	
	static void fsub(double x, dar1 z, dar1 y, dar1 f) {
		f(1) = z(2); // z'(x)
		f(2) = sin(0.6 * z(2)) + x; // z''(x)
	}
	static void dfsub(double x, dar1 z, dar1 y, dar2 df) {
		df(1, 1) = 0; // dz'/dz
		df(1, 2) = 1; // dz'/dz'
		df(2, 1) = 0; // dz'' / dz
		df(2, 2) = 0.6 * cos(0.6 * z(2)); // dz'' / dz'
	}
	static void gsub(int i, dar1 z, double& g) {
		if (i == 1)
			g = z(1) - 1; // 0=z(0)-1 left
		else if (i == 2)
			g = z(2) + 0.1; // 0=z(0)-1 right
		else assert(false);
	}
	static void dgsub(int i, dar1 z, dar1 dg) {
		if (i == 1) {
			dg(1) = 1;  // d/dz  g_left
			dg(2) = 0;  // d/dz' g_left
		}
		else if (i == 2) {
			dg(1) = 0;  // d/dz  g_right
			dg(2) = 1;  // d/dz' g_right
		}
		else assert(false);
	}
};

int main()
{
	sys2 sys;
	

	int iflag;
	cda solver;

	startParams sp;

	for(int i = 0;i<10;++i)
	{
		AutoTimer at(g_timer, "COLDAE "+std::to_string(i));
		solver.COLDAE(sp.ncomp, sp.ny, sp.orders, sp.left, sp.right, sp.bcpoints,
			sp.opts, sp.ltol, sp.tol, sp.fixpnt, sp.ispace, sp.fspace, iflag,
			decltype(sys)::fsub, decltype(sys)::dfsub, decltype(sys)::gsub, decltype(sys)::dgsub, nullptr);
	}
	
	fmt::print("iflag = {}", iflag);

	if (iflag == 1) {
		std::ofstream file("result_cpp.txt");
		for (int i = 1; i <= sp.ispace(1) + 1; ++i) {
			double x = sp.fspace(i);
			dad1 z(2), y(0);
			solver.APPSLN(x, z, y, sp.fspace, sp.ispace);
			file << x << " " << z(1) << std::endl;
		}
		file.close();
	}

	g_timer.print();
	std::cin.get();
	return 0;
}
