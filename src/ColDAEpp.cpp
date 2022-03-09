#pragma warning( disable : 5045 )


#if defined _MSC_VER and _DEBUG
// break on nan
#include <float.h>
unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);
#endif

#include "timer.h"
#include "ColDAEpp.hpp"

#include <iostream>
#include <fstream>

#include "test.h"

int main()
{
	sys3 sys;

	options opts;
	{
		opts.numCollPoints = 0;
		opts.numSubIntervals = 0;
		
		opts.fdim = 1000000;
		opts.idim = 100000;

		opts.printLevel = printMode::none;
		opts.meshSource = meshMode::generate;
		opts.guessSource = guessMode::none;	
		
		opts.ltol = { 1 };
		opts.tol = { 0.0001 };

		opts.numFixedPoints = 0;
	}

	std::vector<int> ispace(opts.idim);
	std::vector<double> fspace(opts.fdim);

    cda solver{};
	output_t iflag;
	auto t0 = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 50000; ++i)
	{
		AutoTimer at(g_timer, "COLDAE");
	
		iflag = solver.COLDAE(sys.params, opts, ispace.data(), fspace.data(),
			decltype(sys)::fsub, decltype(sys)::dfsub, decltype(sys)::gsub,
			decltype(sys)::dgsub, nullptr);
	}
	auto t1 = std::chrono::high_resolution_clock::now();
	std::cout << "time for 5000: " << std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()*1e-6 << std::endl;
	
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
#ifdef _MSC_VER
	std::cin.get();
#endif // _MSC_VER
	return 0;
}
