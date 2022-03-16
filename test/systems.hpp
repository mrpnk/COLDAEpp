#pragma once
#include "../ColDAEpp.hpp"


// f'(x) = x, f(0)=1
struct sys1 {
    coldae::systemParams params;
    sys1() {
        params.ncomp = 1;
        params.ny = 0;
        params.orders = { 1 };
        params.left = 0;
        params.right = 1;
        params.bcpoints = { 0 };
        params.isNonLinear = false;
        params.regularity = coldae::regularControl::regular;
    }

    static void fsub(double x, const double z[1], const double y[0], double f[1]) {
        f[0] = x;
    }
    static void dfsub(double x, const double z[1], const double y[0], double df[1*1]) {
        df[0 + 0*1] = 0;
    }
    static void gsub(int i, const double z[1], double& g) {
        g = z[0] - 1;
    }
    static void dgsub(int i, const double z[1], double dg[1*1]) {
        dg[0] = 1;
    }
};


// f''(x) = sin(0.6*f'(x)) + x, f(0)=1, f(1)=-0.1
struct sys2 {
    coldae::systemParams params;
    sys2() {
        params.ncomp = 2;
        params.ny = 0;
        params.orders = { 1,1 };
        params.left = 0;
        params.right = 1;
        params.bcpoints = { 0, 1 };
        params.isNonLinear = true;
        params.regularity = coldae::regularControl::regular;
        params.index = coldae::indexControl::automatic; // not used
    }

    static void fsub(double x, const double z[2], const double y[], double f[2]) {
        f[0] = z[1]; // z'(x)
        f[1] = sin(0.6 * z[1]) + x; // z''(x)
    }
    static void dfsub(double x, const double z[2], const double y[], double df[2*2]) {
        df[0 + 0*2] = 0; // dz'/dz
        df[0 + 1*2] = 1; // dz'/dz'
        df[1 + 0*2] = 0; // dz'' / dz
        df[1 + 1*2] = 0.6 * cos(0.6 * z[1]); // dz'' / dz'
    }
    static void gsub(int i, const double z[2], double& g) {
        if (i == 1)
            g = z[0] - 1; // 0=z(0)-1  left
        else if (i == 2)
            g = z[1] + 0.1; // 0=z(0)+0.1  right
        else throw(std::invalid_argument("gsub: index out of range."));
    }
    static void dgsub(int i, const double z[2], double dg[2]) {
        if (i == 1) {
            dg[0] = 1;  // d/dz  g_left
            dg[1] = 0;  // d/dz' g_left
        }
        else if (i == 2) {
            dg[0] = 0;  // d/dz  g_right
            dg[1] = 1;  // d/dz' g_right
        }
        else throw(std::invalid_argument("dgsub: index out of range."));
    }
};


// f1''(x) = f2+f1', 0 = f2 + f1*f1' - x, f1(0)=1, f1'(0)=0.1
struct sys3 {
    coldae::systemParams params;
    sys3() {
        params.ncomp = 1;
        params.ny = 1;
        params.orders = { 2 };
        params.left = 0;
        params.right = 5;
        params.bcpoints = { 0, 5 };
        params.isNonLinear = true;
        params.regularity = coldae::regularControl::regular;
        params.index = coldae::indexControl::automatic;
    }

    static void fsub(double x, double const z[2], double const y[1], double f[2])  {
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
        else throw(std::invalid_argument("gsub: index out of range."));
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
        else throw(std::invalid_argument("dgsub: index out of range."));
    }
};

