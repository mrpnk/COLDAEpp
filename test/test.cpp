#include "systems.hpp"

#include <iostream>
#include <fstream>
#include <gtest/gtest.h>


template<typename T>
void compareSolution(T const& sys, coldae::options const& opts, std::string comparisonFile) {
    std::vector<int> ispace(opts.idim);
    std::vector<double> fspace(opts.fdim);

    coldae::cda solver{};
    coldae::result_t iflag;

    iflag = solver.COLDAE(sys.params, opts, ispace.data(), fspace.data(),
                          sys, nullptr);

    ASSERT_EQ(iflag, coldae::result_t::normal) << "COLDAE did not return normally.\n";

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
    coldae::options opts;
    {
        opts.fdim = 10000;
        opts.idim = 10000;

        opts.printLevel = coldae::printMode::none;
        opts.meshSource = coldae::meshMode::generate;
        opts.guessSource = coldae::guessMode::none;

        opts.ltol = {1};
        opts.tol = {0.1};
    }
    compareSolution(sys1{}, opts, "original/result1.txt");
}
TEST(BasicUses, twoComponents) {
    coldae::options opts;
    {
        opts.fdim = 10000;
        opts.idim = 10000;

        opts.printLevel = coldae::printMode::none;
        opts.meshSource = coldae::meshMode::generate;
        opts.guessSource = coldae::guessMode::none;

        opts.ltol = {1};
        opts.tol = {0.0001};
    }
    compareSolution(sys2{}, opts, "original/result2.txt");
}
TEST(BasicUses, algebraic) {
    coldae::options opts;
    {
        opts.fdim = 1000000;
        opts.idim = 100000;

        opts.printLevel = coldae::printMode::none;
        opts.meshSource = coldae::meshMode::generate;
        opts.guessSource = coldae::guessMode::none;

        opts.ltol = {1,2};
        opts.tol = {0.0001, 0.0001};
    }
    compareSolution(sys3{}, opts, "original/result3.txt");
}

