#include "systems.hpp"
#include "timer.h"

#include <iostream>
#include <gtest/gtest.h>

void benchmarkSolution(coldae::system const& sys, coldae::options const& opts, int nRuns) {
    std::vector<int> ispace(opts.idim);
    std::vector<double> fspace(opts.fdim);

    coldae::cda solver{};

    using namespace std::chrono;
    auto t0 = high_resolution_clock::now();
    for(int i = 0; i<nRuns; ++i) {
        solver.COLDAE(sys.params, opts, ispace.data(), fspace.data(),
                      sys, nullptr);
    }
    auto t1 = high_resolution_clock::now();
    fmt::print(fg(fmt::color::cornflower_blue),
               "Time for {} runs: {}ms\n", nRuns, duration_cast<milliseconds>(t1-t0).count());
}


TEST(Benchmarks, algebraic) {
    coldae::options opts;
    {
    opts.numCollPoints = 0;
    opts.numSubIntervals = 0;

    opts.fdim = 1000000;
    opts.idim = 100000;

    opts.printLevel = coldae::printMode::none;
    opts.meshSource = coldae::meshMode::generate;
    opts.guessSource = coldae::guessMode::none;

    opts.ltol = {1,2};
    opts.tol = {0.0001, 0.0001};

    opts.numFixedPoints = 0;
    }
    benchmarkSolution(sys3{}, opts,50000);
}
