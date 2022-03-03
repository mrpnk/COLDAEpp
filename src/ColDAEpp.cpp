#pragma warning( disable : 5045 )

// break on nan
#include <float.h>
unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);

#include "timer.h"
#include "ColDAEpp.hpp"



// f'(x) = x, f(0)=1
struct sys1 {
	systemParams params;
	sys1() {
		params.ncomp = 1;
		params.ny = 0;
		params.orders = { 1 };
		params.left = 0;
		params.right = 1;
		params.bcpoints = { 0 };
		params.isNonLinear = 1;
		params.reg = regularControl::regular;
		params.index = indexControl::automatic; // not used
	}

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
	systemParams params;
	sys2() {
		params.ncomp = 2;
		params.ny = 0;
		params.orders = { 1,1 };
		params.left = 0;
		params.right = 1;
		params.bcpoints = { 0, 1 };
		params.isNonLinear = 1;
		params.reg = regularControl::regular;
		params.index = indexControl::automatic; // not used
	}

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

	options opts;
	{
		opts.numCollPoints = 0;
		opts.numSubIntervals = 0;
		opts.numTolerances = 1;

		opts.fdim = 10000;
		opts.idim = 10000;

		opts.printLevel = printMode::full;
		opts.meshSource = meshMode::generate;
		opts.guessSource = guessMode::none;

		opts.numFixedPoints = 0;
		opts.ltol = { 1 };
		opts.tol = { 0.0001 };
	}

	iad1 ispace(opts.idim);
	dad1 fspace(opts.fdim);

	cda solver;
	output_t iflag;
	{
		AutoTimer at(g_timer, "COLDAE");
		solver.COLDAE(sys.params, opts, ispace, fspace, iflag,
			decltype(sys)::fsub, decltype(sys)::dfsub, decltype(sys)::gsub,
			decltype(sys)::dgsub, nullptr);
	}
	
	
	if (iflag == output_t::normal) {
		fmt::print(fg(fmt::color::green_yellow), "Successful return!\n");
		std::ofstream file("result_cpp.txt");
		for (int i = 1; i <= ispace(1) + 1; ++i) {
			double x = fspace(i);
			dad1 z(2), y(0);
			solver.APPSLN(x, z, y, fspace, ispace);
			file << x << " " << z(1) << std::endl;
		}
		file.close();
	}
	else
		fmt::print(fg(fmt::color::red), "Error return!\n");

	g_timer.print();
	std::cin.get();
	return 0;
}
