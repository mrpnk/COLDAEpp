#include "test.hpp"
#include "timer.h"
#include "ColDAEpp.hpp"

#include <iostream>
#include <fstream>
#include <gtest/gtest.h>
#include <filesystem>

#if defined _MSC_VER
// break on nan
#include <float.h>
unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);
#endif

template<typename T>
void compareSolution(T const& sys, options const& opts, std::string comparisonFile) {
    std::vector<int> ispace(opts.idim);
    std::vector<double> fspace(opts.fdim);

    cda solver{};
    output_t iflag;

    iflag = solver.COLDAE(sys.params, opts, ispace.data(), fspace.data(),
                          T::fsub, T::dfsub, T::gsub,
                          T::dgsub, nullptr);

    ASSERT_EQ(iflag, output_t::normal) << "COLDAE did not return normally.\n";

//    std::ofstream outfile("out.txt");
//    for(int i=1;i<=ispace[0]+1;++i) {
//        double x = fspace[i-1];
//        double z[2], y[1];
//        solver.APPSLN(x, z, y, fspace.data(), ispace.data());
//        outfile << x<<" " <<z[0]<< std::endl;
//    }
//    outfile.close();

    // Compare result with FORTRAN version
    std::ifstream infile(comparisonFile);



    std::filesystem::path cwd = std::filesystem::current_path();
    //EXPECT_EQ(cwd.string(),"");

    ASSERT_EQ(infile.fail(), false) << "The Fortran output file for comparison was not found.\n";
    double maxError = 0;
    int count = 0;
    int ncomp=sys.params.getMstar(), ny=sys.params.ny;
    while (true) {
        double x;
        std::vector<double> z(ncomp), y(ny);

        infile >> x;
        for(int i = 0;i<ncomp;++i) infile>>z[i];
        for(int i = 0;i<ny;++i) infile>>y[i];
        if(!infile) break;

        std::vector<double> myz(ncomp), myy(ny);
        solver.APPSLN(x, myz.data(), myy.data(),fspace.data(), ispace.data());
        double errorNorm = 0;
        for(int i = 0;i<ncomp;++i) errorNorm += std::pow(myz[i]-z[i],2);
        for(int i = 0;i<ny;++i) errorNorm += std::pow(myy[i]-y[i],2);
        maxError = std::max(maxError, errorNorm);
        count++;
    }
    infile.close();
    EXPECT_EQ(ispace[0] + 1, count) << "We have a wrong number of mesh points.\n";
    EXPECT_LE(maxError, 1e-12) << "The deviation of the solution is too big.\n";
}


TEST(BasicUses, firstOrder) {
    options opts;
    {
        opts.numCollPoints = 0;
        opts.numSubIntervals = 0;

        opts.fdim = 10000;
        opts.idim = 10000;

        opts.printLevel = printMode::none;
        opts.meshSource = meshMode::generate;
        opts.guessSource = guessMode::none;

        opts.ltol = {1};
        opts.tol = {0.1};

        opts.numFixedPoints = 0;
    }
    compareSolution(sys1{}, opts, "original/result1.txt");
}
TEST(BasicUses, twoComponents) {
    options opts;
    {
        opts.numCollPoints = 0;
        opts.numSubIntervals = 0;

        opts.fdim = 10000;
        opts.idim = 10000;

        opts.printLevel = printMode::none;
        opts.meshSource = meshMode::generate;
        opts.guessSource = guessMode::none;

        opts.ltol = {1};
        opts.tol = {0.0001};

        opts.numFixedPoints = 0;
    }
    compareSolution(sys2{}, opts, "original/result2.txt");
}
TEST(BasicUses, algebraic) {
    options opts;
    {
        opts.numCollPoints = 0;
        opts.numSubIntervals = 0;

        opts.fdim = 1000000;
        opts.idim = 100000;

        opts.printLevel = printMode::none;
        opts.meshSource = meshMode::generate;
        opts.guessSource = guessMode::none;

        opts.ltol = {1,2};
        opts.tol = {0.0001, 0.0001};

        opts.numFixedPoints = 0;
    }
    compareSolution(sys3{}, opts, "original/result3.txt");
}
