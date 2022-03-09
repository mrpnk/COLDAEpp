#include "test1.hpp"
#include "timer.h"
#include "ColDAEpp.hpp"

#include <iostream>
#include <fstream>
#include <gtest/gtest.h>

#if defined _MSC_VER
// break on nan
#include <float.h>
unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);
#endif

template<typename T>
void compareSolution(T&& sys, std::string comparisonFile) {
    options opts;
    {
        opts.numCollPoints = 0;
        opts.numSubIntervals = 0;

        opts.fdim = 1000000;
        opts.idim = 100000;

        opts.printLevel = printMode::none;
        opts.meshSource = meshMode::generate;
        opts.guessSource = guessMode::none;

        opts.ltol = {1};
        opts.tol = {0.0001};

        opts.numFixedPoints = 0;
    }

    std::vector<int> ispace(opts.idim);
    std::vector<double> fspace(opts.fdim);

    cda solver{};
    output_t iflag;

    iflag = solver.COLDAE(sys.params, opts, ispace.data(), fspace.data(),
                          T::fsub, T::dfsub, T::gsub,
                          T::dgsub, nullptr);

    ASSERT_EQ(iflag, output_t::normal) << "COLDAE did not return normally.\n";


    // Compare result with FORTRAN version
    std::ifstream infile(comparisonFile);
    EXPECT_EQ(infile.fail(), false) << "The Fortran file for comparison is not found.\n";
    double maxError = 0, x, z1, y1;
    int count = 0;
    while (infile >> x >> z1 >> y1) {
        double z[2], y[1];
        solver.APPSLN(x, z, y, fspace.data(), ispace.data());
        maxError = std::max({maxError, z1 - z[0], y1 - y[0]});
        count++;
    }
    infile.close();
    EXPECT_EQ(ispace[0] + 1, count) << "We have a wrong number of mesh points.\n";
    EXPECT_LE(maxError, 1e-6) << "The deviation of the solution is too big.\n";
}


TEST(BasicUses, firstOrder) {
    compareSolution(sys1{}, "../../original/result_f77.txt");
}
TEST(BasicUses, twoComponents) {
    compareSolution(sys2{}, "../../original/result_f77.txt");
}
TEST(BasicUses, algebraic) {
    compareSolution(sys3{}, "../../original/result_f77.txt");
}
