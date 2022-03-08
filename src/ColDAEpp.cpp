#pragma warning( disable : 5045 )

// break on nan
#include <float.h>
unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);

#include "timer.h"
#include "ColDAEpp.hpp"



//
//// f'(x) = x, f(0)=1
//struct sys1 {
//	systemParams params;
//	sys1() {
//		params.ncomp = 1;
//		params.ny = 0;
//		params.orders = { 1 };
//		params.left = 0;
//		params.right = 1;
//		params.bcpoints = { 0 };
//		params.isNonLinear = 1;
//		params.reg = regularControl::regular;
//		params.index = indexControl::automatic; // not used
//	}
//
//	static void fsub(double x, dar1 z, dar1 y, dar1 f) {
//		f(1) = x;
//	}
//	static void dfsub(double x, dar1 z, dar1 y, dar2 df) {
//		df(1, 1) = 0;
//	}
//	static void gsub(int i, dar1 z, double& g) {
//		g = z(1) - 1;
//	}
//	static void dgsub(int i, dar1 z, dar1 dg) {
//		dg(1) = 1;
//	}
//};
//
//// f''(x) = sin(0.6*f'(x)) + x, f(0)=1, f(1)=-0.1
//struct sys2 {
//	systemParams params;
//	sys2() {
//		params.ncomp = 2;
//		params.ny = 0;
//		params.orders = { 1,1 };
//		params.left = 0;
//		params.right = 1;
//		params.bcpoints = { 0, 1 };
//		params.isNonLinear = 1;
//		params.reg = regularControl::regular;
//		params.index = indexControl::automatic; // not used
//	}
//
//	static void fsub(double x, dar1 z, dar1 y, dar1 f) {
//		f(1) = z(2); // z'(x)
//		f(2) = sin(0.6 * z(2)) + x; // z''(x)
//	}
//	static void dfsub(double x, dar1 z, dar1 y, dar2 df) {
//		df(1, 1) = 0; // dz'/dz
//		df(1, 2) = 1; // dz'/dz'
//		df(2, 1) = 0; // dz'' / dz
//		df(2, 2) = 0.6 * cos(0.6 * z(2)); // dz'' / dz'
//	}
//	static void gsub(int i, dar1 z, double& g) {
//		if (i == 1)
//			g = z(1) - 1; // 0=z(0)-1 left
//		else if (i == 2)
//			g = z(2) + 0.1; // 0=z(0)-1 right
//		else assert(false);
//	}
//	static void dgsub(int i, dar1 z, dar1 dg) {
//		if (i == 1) {
//			dg(1) = 1;  // d/dz  g_left
//			dg(2) = 0;  // d/dz' g_left
//		}
//		else if (i == 2) {
//			dg(1) = 0;  // d/dz  g_right
//			dg(2) = 1;  // d/dz' g_right
//		}
//		else assert(false);
//	}
//};
//


// f1''(x) = f2+f1', 0 = f2 + f1*f1' - x, f1(0)=1, f1'(0)=0.1
struct sys3 {
	systemParams params;
	sys3() {
		params.ncomp = 1;
		params.ny = 1;
		params.orders = { 2 };
		params.left = 0;
		params.right = 5;
		params.bcpoints = { 0, 5 };
		params.isNonLinear = true;
		params.reg = regularControl::regular;
		params.index = indexControl::automatic;
	}

	static void fsub(double x, double const z[2], double const y[1], double f[2]) {
		f[0] = y[0] + z[1]; // f_1 = f1''(x)
		f[1] = y[0] + z[0] * z[1] - x; // f_2 = 0 = algebraic
	}
	static void dfsub(double x, double const z[2], double const y[1], double df[6]) {
		const int ms = 2;
		df[0 + 0 * ms] = 0; // df_1 / df1
		df[0 + 1 * ms] = 1; // df_1 / df1'
		df[0 + 2 * ms] = 1; // df_1 / df2

		df[1 + 0 * ms] = z[1]; // df_2 / df1
		df[1 + 1 * ms] = z[0]; // df_2 / df1'
		df[1 + 2 * ms] = 1;    // df_2 / df2
	}
	static void gsub(int i, double const z[2], double& g) {
		if (i == 1)
			g = z[0] + z[1] + 0.15;  // g_1 = 0 = z(0)-1   left
		else if (i == 2)
			g = z[0] - 6.5;          // g_2 = 0 = z(0)-1   right
		else assert(false);
	}
	static void dgsub(int i, double const z[2], double dg[2]) {
		if (i == 1) {
			dg[0] = 1; // dg_1 / df1
			dg[1] = 1; // dg_1 / df1'
		}
		else if (i == 2) {
			dg[0] = 1;  // dg_2/df1
			dg[1] = 0;  // dg_2/df1'
		}
		else assert(false);
	}
};


int main()
{
	sys3 sys;

	options opts;
	{
		opts.numCollPoints = 0;
		opts.numSubIntervals = 0;
		
		opts.fdim = 1000000;
		opts.idim = 100000;

		opts.printLevel = printMode::full;
		opts.meshSource = meshMode::generate;
		opts.guessSource = guessMode::none;	
		
		opts.ltol = { 1 };
		opts.tol = { 0.0001 };

		opts.numFixedPoints = 0;
	}

	std::vector<int> ispace(opts.idim);
	std::vector<double> fspace(opts.fdim);

	cda solver;
	output_t iflag;
	for (int i = 0; i < 1; ++i) 
	{
		AutoTimer at(g_timer, "COLDAE");
	
		iflag = solver.COLDAE(sys.params, opts, ispace.data(), fspace.data(),
			decltype(sys)::fsub, decltype(sys)::dfsub, decltype(sys)::gsub,
			decltype(sys)::dgsub, nullptr);
	}
	
	
	if (iflag == output_t::normal) {
		fmt::print(fg(fmt::color::green_yellow), "Successful return!\n");
		std::ofstream file("result_cpp.txt");
		for (int i = 1; i <= ispace[0] + 1; ++i) {
			double x = fspace[i-1];
			double z[2], y[1];
			solver.APPSLN(x, z, y, fspace.data(), ispace.data());
			file << x << " " << z[0] << " " << y[0] << std::endl;
		}
		file.close();



		// Compare result with FORTRAN version
		std::ifstream infile("../../original/result_f77.txt");
		double maxError = 0, x, z1, y1;
		while (infile >> x >> z1 >> y1) {
			double z[2], y[1];
			solver.APPSLN(x, z, y, fspace.data(), ispace.data());
			maxError = std::max({ maxError,z1 - z[0],y1 - y[0] });
		}
		fmt::print(fg(fmt::color::cornflower_blue), "Max deviation from FORTRAN is {}.\n", maxError);
		infile.close();
	}
	else
		fmt::print(fg(fmt::color::red), "Error return!\n");

	g_timer.print();
	std::cin.get();
	return 0;
}
