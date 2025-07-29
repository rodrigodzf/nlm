#pragma once
#include <Eigen/Sparse>
#include <vector>
#include <tuple>
#include <cmath> // For std::pow

template <typename T>
std::tuple<T, T, T, T, T, T> 
D00_coeffs(T K0y, T R0y, T Kx0, T Rx0, T h, T D, T nu);

template <typename T>
std::tuple<T, T, T, T, T, T, T, T> 
D01_coeffs(T K0y, T R0y, T Rx0, T h, T D, T nu);

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T> 
D02_coeffs(T K0y, T R0y, T h, T D, T nu);

template <typename T>
std::tuple<T, T, T, T, T, T, T, T> 
D10_coeffs(T R0y, T Kx0, T Rx0, T h, T D, T nu);

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T, T, T> 
D11_coeffs(T R0y, T Rx0, T h, T D, T nu);

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T, T, T, T> 
D12_coeffs(T R0y, T h, T D, T nu);

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T> 
D20_coeffs(T Kx0, T Rx0, T h, T D, T nu);

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T, T, T, T> 
D21_coeffs(T Rx0, T h, T D, T nu);

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T, T, T, T, T> 
D22_coeffs(); 

// --- Implementations ---

template <typename T>
std::tuple<T, T, T, T, T, T>
D00_coeffs(T K0y, T R0y, T Kx0, T Rx0, T h, T D, T nu) {
    T D2 = D * D;
    T D3 = D2 * D;
    T D4 = D3 * D;
    T h2 = h * h;
    T h3 = h2 * h;
    T h4 = h3 * h;
    T h5 = h4 * h;
    T h6 = h5 * h;
    T nu2 = nu * nu;
    T nu3 = nu2 * nu;
    T nu4 = nu3 * nu;

    T R0y2 = R0y * R0y;
    T Rx02 = Rx0 * Rx0;

    T term1_num = (static_cast<T>(2) * K0y * Rx0 * R0y2 * h6 + static_cast<T>(4) * K0y * D * R0y2 * h5 + static_cast<T>(8) * K0y * Rx0 * D * R0y * h5 - static_cast<T>(8) * K0y * D2 * R0y * h4 * nu2 + static_cast<T>(16) * K0y * D2 * R0y * h4 - static_cast<T>(12) * Rx0 * D2 * R0y * h2 * nu2 + static_cast<T>(24) * Rx0 * D2 * R0y * h2 * nu + static_cast<T>(24) * Rx0 * D2 * R0y * h2 + static_cast<T>(16) * D3 * R0y * h * nu3 - static_cast<T>(56) * D3 * R0y * h * nu2 + static_cast<T>(48) * D3 * R0y * h + static_cast<T>(8) * K0y * Rx0 * D2 * h4 - static_cast<T>(16) * K0y * D3 * h3 * nu2 + static_cast<T>(16) * K0y * D3 * h3 - static_cast<T>(24) * Rx0 * D3 * h * nu2 + static_cast<T>(48) * Rx0 * D3 * h * nu + static_cast<T>(48) * Rx0 * D3 * h + static_cast<T>(16) * D4 * nu4 - static_cast<T>(112) * D4 * nu2 + static_cast<T>(96) * D4);
    T term1_den = (D * (static_cast<T>(2) * D + R0y * h) * (static_cast<T>(4) * D2 - static_cast<T>(4) * D2 * nu2 + static_cast<T>(2) * D * R0y * h + static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2));
    T term1 = term1_num / term1_den;

    T term2_num = (static_cast<T>(8) * (- static_cast<T>(8) * D2 * nu2 + static_cast<T>(4) * Rx0 * h * D * nu + static_cast<T>(8) * D2 + static_cast<T>(4) * Rx0 * h * D));
    T term2_den = (static_cast<T>(4) * D2 - static_cast<T>(4) * D2 * nu2 + static_cast<T>(2) * D * R0y * h + static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2);
    T term2 = term2_num / term2_den;

    T term3 = (static_cast<T>(4) * D * nu) / (static_cast<T>(2) * D + R0y * h);
    T term4 = (static_cast<T>(4) * D * nu) / (static_cast<T>(2) * D + Rx0 * h);

    T term5_num = (static_cast<T>(2) * (static_cast<T>(8) * D2 * nu + static_cast<T>(2) * D * R0y * h * nu + static_cast<T>(2) * D * Rx0 * h * nu));
    T term5_den = ((static_cast<T>(2) * D + R0y * h) * (static_cast<T>(2) * D + Rx0 * h));
    T term5 = term5_num / term5_den;

    T term6_num = (static_cast<T>(8) * (- static_cast<T>(8) * D2 * nu2 + static_cast<T>(4) * R0y * h * D * nu + static_cast<T>(8) * D2 + static_cast<T>(4) * R0y * h * D));
    T term6_den = (static_cast<T>(4) * D2 - static_cast<T>(4) * D2 * nu2 + static_cast<T>(2) * D * R0y * h + static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2);
    T term6 = term6_num / term6_den;

    T term7_num = (static_cast<T>(2) * Kx0 * R0y * Rx02 * h6 + static_cast<T>(4) * Kx0 * D * Rx02 * h5 + static_cast<T>(8) * Kx0 * R0y * D * Rx0 * h5 - static_cast<T>(8) * Kx0 * D2 * Rx0 * h4 * nu2 + static_cast<T>(16) * Kx0 * D2 * Rx0 * h4 - static_cast<T>(12) * R0y * D2 * Rx0 * h2 * nu2 + static_cast<T>(24) * R0y * D2 * Rx0 * h2 * nu + static_cast<T>(24) * R0y * D2 * Rx0 * h2 + static_cast<T>(16) * D3 * Rx0 * h * nu3 - static_cast<T>(56) * D3 * Rx0 * h * nu2 + static_cast<T>(48) * D3 * Rx0 * h + static_cast<T>(8) * Kx0 * R0y * D2 * h4 - static_cast<T>(16) * Kx0 * D3 * h3 * nu2 + static_cast<T>(16) * Kx0 * D3 * h3 - static_cast<T>(24) * R0y * D3 * h * nu2 + static_cast<T>(48) * R0y * D3 * h * nu + static_cast<T>(48) * R0y * D3 * h + static_cast<T>(16) * D4 * nu4 - static_cast<T>(112) * D4 * nu2 + static_cast<T>(96) * D4);
    T term7_den = (D * (static_cast<T>(2) * D + Rx0 * h) * (static_cast<T>(4) * D2 - static_cast<T>(4) * D2 * nu2 + static_cast<T>(2) * D * R0y * h + static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2));
    T term7 = term7_num / term7_den;

    T u00_c = term1 - term2 - term3 - term4 - term5 - term6 + term7 + static_cast<T>(20);

    T den_common = (static_cast<T>(4) * D2 - static_cast<T>(4) * D2 * nu2 + static_cast<T>(2) * D * R0y * h + static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2);
    T den_prod = ((static_cast<T>(2) * D + R0y * h) * (static_cast<T>(2) * D + Rx0 * h));
    T den_d_term1 = D * (static_cast<T>(2) * D + R0y * h) * den_common;
    T den_d_term2 = D * (static_cast<T>(2) * D + Rx0 * h) * den_common;

    T u10_c_term1 = (static_cast<T>(2) * (static_cast<T>(4) * D + static_cast<T>(4) * D * nu)) / (static_cast<T>(2) * D + Rx0 * h);
    T u10_c_term2 = (static_cast<T>(8) * (static_cast<T>(4) * D2 * nu2 - static_cast<T>(4) * D2 + static_cast<T>(2) * D * R0y * h - static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2)) / den_common;
    T u10_c_term3 = (static_cast<T>(2) * (static_cast<T>(8) * D2 * nu + static_cast<T>(8) * D2 + static_cast<T>(4) * D * R0y * h + static_cast<T>(4) * D * R0y * h * nu)) / den_prod;
    T u10_c_term4_num = (static_cast<T>(32) * D4 * nu - static_cast<T>(96) * D4 + static_cast<T>(96) * D4 * nu2 - static_cast<T>(32) * D4 * nu3 - static_cast<T>(48) * D3 * R0y * h - static_cast<T>(48) * D3 * Rx0 * h + static_cast<T>(16) * D3 * R0y * h * nu + static_cast<T>(16) * D3 * Rx0 * h * nu - static_cast<T>(24) * D2 * R0y * Rx0 * h2 + static_cast<T>(48) * D3 * R0y * h * nu2 - static_cast<T>(16) * D3 * R0y * h * nu3 + static_cast<T>(8) * D2 * R0y * Rx0 * h2 * nu);
    T u10_c_term4 = u10_c_term4_num / den_d_term1;
    T u10_c_term5_num = (static_cast<T>(32) * D4 * nu + static_cast<T>(64) * D4 - static_cast<T>(96) * D4 * nu2 - static_cast<T>(32) * D4 * nu3 + static_cast<T>(32) * D4 * nu4 + static_cast<T>(32) * D3 * R0y * h + static_cast<T>(32) * D3 * Rx0 * h + static_cast<T>(64) * D3 * R0y * h * nu + static_cast<T>(16) * D3 * Rx0 * h * nu + static_cast<T>(16) * D2 * R0y * Rx0 * h2 - static_cast<T>(32) * D3 * R0y * h * nu2 - static_cast<T>(16) * D3 * Rx0 * h * nu2 - static_cast<T>(16) * D2 * R0y * Rx0 * h2 * nu2 + static_cast<T>(32) * D2 * R0y * Rx0 * h2 * nu);
    T u10_c_term5 = u10_c_term5_num / den_d_term2;
    T u10_c_term6 = (static_cast<T>(32) * D * R0y * h * nu) / den_common;
    T u10_c = u10_c_term1 - u10_c_term2 + u10_c_term3 + u10_c_term4 - u10_c_term5 + u10_c_term6 - static_cast<T>(8);

    T u20_c_term1_num = (Rx0 * D * R0y2 * h3 + static_cast<T>(2) * D2 * R0y2 * h2 + static_cast<T>(4) * Rx0 * D2 * R0y * h2 - static_cast<T>(4) * D3 * R0y * h * nu2 + static_cast<T>(8) * D3 * R0y * h + static_cast<T>(4) * Rx0 * D3 * h - static_cast<T>(8) * D4 * nu2 + static_cast<T>(8) * D4);
    T u20_c_term1 = u20_c_term1_num / den_d_term1;
    T u20_c_term2 = (static_cast<T>(2) * (static_cast<T>(4) * nu * D2 + static_cast<T>(2) * R0y * h * nu * D)) / den_prod;
    T u20_c_term3 = (static_cast<T>(4) * D * nu) / (static_cast<T>(2) * D + Rx0 * h);
    T u20_c_term4_num = (static_cast<T>(32) * D4 * nu - static_cast<T>(16) * D4 * nu2 - static_cast<T>(32) * D4 * nu3 + static_cast<T>(16) * D4 * nu4 + static_cast<T>(16) * D3 * R0y * h * nu + static_cast<T>(16) * D3 * Rx0 * h * nu - static_cast<T>(8) * D3 * R0y * h * nu2 - static_cast<T>(8) * D3 * Rx0 * h * nu2 - static_cast<T>(4) * D2 * R0y * Rx0 * h2 * nu2 + static_cast<T>(8) * D2 * R0y * Rx0 * h2 * nu);
    T u20_c_term4 = u20_c_term4_num / den_d_term2;
    T u20_c = u20_c_term1 - u20_c_term2 - u20_c_term3 + u20_c_term4 + static_cast<T>(1);

    T u01_c_term1 = (static_cast<T>(2) * (static_cast<T>(4) * D + static_cast<T>(4) * D * nu)) / (static_cast<T>(2) * D + R0y * h);
    T u01_c_term2 = (static_cast<T>(8) * (static_cast<T>(4) * D2 * nu2 - static_cast<T>(4) * D2 - static_cast<T>(2) * D * R0y * h + static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2)) / den_common;
    T u01_c_term3 = (static_cast<T>(2) * (static_cast<T>(8) * D2 * nu + static_cast<T>(8) * D2 + static_cast<T>(4) * D * Rx0 * h + static_cast<T>(4) * D * Rx0 * h * nu)) / den_prod;
    T u01_c_term4_num = (static_cast<T>(32) * D4 * nu - static_cast<T>(96) * D4 + static_cast<T>(96) * D4 * nu2 - static_cast<T>(32) * D4 * nu3 - static_cast<T>(48) * D3 * R0y * h - static_cast<T>(48) * D3 * Rx0 * h + static_cast<T>(16) * D3 * R0y * h * nu + static_cast<T>(16) * D3 * Rx0 * h * nu - static_cast<T>(24) * D2 * R0y * Rx0 * h2 + static_cast<T>(48) * D3 * Rx0 * h * nu2 - static_cast<T>(16) * D3 * Rx0 * h * nu3 + static_cast<T>(8) * D2 * R0y * Rx0 * h2 * nu);
    T u01_c_term4 = u01_c_term4_num / den_d_term2;
    T u01_c_term5_num = (static_cast<T>(32) * D4 * nu + static_cast<T>(64) * D4 - static_cast<T>(96) * D4 * nu2 - static_cast<T>(32) * D4 * nu3 + static_cast<T>(32) * D4 * nu4 + static_cast<T>(32) * D3 * R0y * h + static_cast<T>(32) * D3 * Rx0 * h + static_cast<T>(16) * D3 * R0y * h * nu + static_cast<T>(64) * D3 * Rx0 * h * nu + static_cast<T>(16) * D2 * R0y * Rx0 * h2 - static_cast<T>(16) * D3 * R0y * h * nu2 - static_cast<T>(32) * D3 * Rx0 * h * nu2 - static_cast<T>(16) * D2 * R0y * Rx0 * h2 * nu2 + static_cast<T>(32) * D2 * R0y * Rx0 * h2 * nu);
    T u01_c_term5 = u01_c_term5_num / den_d_term1;
    T u01_c_term6 = (static_cast<T>(32) * D * Rx0 * h * nu) / den_common;
    T u01_c = u01_c_term1 - u01_c_term2 + u01_c_term3 + u01_c_term4 - u01_c_term5 + u01_c_term6 - static_cast<T>(8);

    T u02_c_term1_num = (R0y * D * Rx02 * h3 + static_cast<T>(2) * D2 * Rx02 * h2 + static_cast<T>(4) * R0y * D2 * Rx0 * h2 - static_cast<T>(4) * D3 * Rx0 * h * nu2 + static_cast<T>(8) * D3 * Rx0 * h + static_cast<T>(4) * R0y * D3 * h - static_cast<T>(8) * D4 * nu2 + static_cast<T>(8) * D4);
    T u02_c_term1 = u02_c_term1_num / den_d_term2;
    T u02_c_term2 = (static_cast<T>(2) * (static_cast<T>(4) * nu * D2 + static_cast<T>(2) * Rx0 * h * nu * D)) / den_prod;
    T u02_c_term3 = (static_cast<T>(4) * D * nu) / (static_cast<T>(2) * D + R0y * h);
    T u02_c_term4 = u20_c_term4_num / den_d_term1; // Reusing calculation from u20_c
    T u02_c = u02_c_term1 - u02_c_term2 - u02_c_term3 + u02_c_term4 + static_cast<T>(1);

    T u11_c_term1 = (static_cast<T>(2) * (static_cast<T>(2) * D - Rx0 * h)) / (static_cast<T>(2) * D + Rx0 * h);
    T u11_c_term2 = (static_cast<T>(2) * (static_cast<T>(12) * D2 + static_cast<T>(2) * D * R0y * h + static_cast<T>(2) * D * Rx0 * h - R0y * Rx0 * h2)) / den_prod;
    T u11_c_term3_num = (static_cast<T>(32) * D4 * nu - static_cast<T>(64) * D4 + static_cast<T>(64) * D4 * nu2 - static_cast<T>(32) * D4 * nu3 - static_cast<T>(32) * D3 * R0y * h - static_cast<T>(32) * D3 * Rx0 * h + static_cast<T>(16) * D3 * R0y * h * nu + static_cast<T>(16) * D3 * Rx0 * h * nu - static_cast<T>(16) * D2 * R0y * Rx0 * h2 + static_cast<T>(8) * D2 * R0y * Rx0 * h2 * nu);
    T u11_c_term3 = u11_c_term3_num / den_d_term1;
    T u11_c_term4 = u11_c_term3_num / den_d_term2; // Reusing numerator from term3
    T u11_c_term5 = (static_cast<T>(2) * (static_cast<T>(2) * D - R0y * h)) / (static_cast<T>(2) * D + R0y * h);
    T u11_c = static_cast<T>(2) - u11_c_term1 - u11_c_term2 - u11_c_term3 - u11_c_term4 - u11_c_term5;

    return std::make_tuple(u00_c, u10_c, u20_c, u01_c, u02_c, u11_c);
}

template <typename T>
std::tuple<T, T, T, T, T, T, T, T>
D01_coeffs(T K0y, T R0y, T Rx0, T h, T D, T nu) {
    T D2 = D * D;
    T D3 = D2 * D;
    T D4 = D3 * D;
    T h2 = h * h;
    T h3 = h2 * h;
    T h4 = h3 * h;
    T h5 = h4 * h;
    T h6 = h5 * h;
    T nu2 = nu * nu;
    T nu3 = nu2 * nu;
    T nu4 = nu3 * nu;
    T R0y2 = R0y * R0y;

    T den_common = (static_cast<T>(4) * D2 - static_cast<T>(4) * D2 * nu2 + static_cast<T>(2) * D * R0y * h + static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2);
    T den_d_term1 = D * (static_cast<T>(2) * D + R0y * h) * den_common;

    T u01_c_term1 = (static_cast<T>(4) * D2 * nu2 - static_cast<T>(4) * D2 - static_cast<T>(2) * D * R0y * h + static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2) / den_common;
    T u01_c_term2 = (static_cast<T>(8) * (static_cast<T>(4) * D + static_cast<T>(4) * D * nu)) / (static_cast<T>(2) * D + R0y * h);
    T u01_c_term3 = (static_cast<T>(4) * D * nu) / (static_cast<T>(2) * D + R0y * h);
    T u01_c_term4_num = (static_cast<T>(2) * K0y * Rx0 * R0y2 * h6 + static_cast<T>(4) * K0y * D * R0y2 * h5 + static_cast<T>(8) * K0y * Rx0 * D * R0y * h5 - static_cast<T>(8) * K0y * D2 * R0y * h4 * nu2 + static_cast<T>(16) * K0y * D2 * R0y * h4 - static_cast<T>(14) * Rx0 * D2 * R0y * h2 * nu2 + static_cast<T>(28) * Rx0 * D2 * R0y * h2 * nu + static_cast<T>(24) * Rx0 * D2 * R0y * h2 - static_cast<T>(20) * D3 * R0y * h * nu2 + static_cast<T>(40) * D3 * R0y * h * nu + static_cast<T>(48) * D3 * R0y * h + static_cast<T>(8) * K0y * Rx0 * D2 * h4 - static_cast<T>(16) * K0y * D3 * h3 * nu2 + static_cast<T>(16) * K0y * D3 * h3 - static_cast<T>(28) * Rx0 * D3 * h * nu2 + static_cast<T>(56) * Rx0 * D3 * h * nu + static_cast<T>(48) * Rx0 * D3 * h + static_cast<T>(40) * D4 * nu4 - static_cast<T>(80) * D4 * nu3 - static_cast<T>(136) * D4 * nu2 + static_cast<T>(80) * D4 * nu + static_cast<T>(96) * D4);
    T u01_c_term4 = u01_c_term4_num / den_d_term1;
    T u01_c_term5 = (static_cast<T>(8) * D * Rx0 * h * nu) / den_common;
    T u01_c = u01_c_term1 - u01_c_term2 - u01_c_term3 + u01_c_term4 - u01_c_term5 + static_cast<T>(20);

    T u11_c_term1 = (static_cast<T>(8) * (static_cast<T>(2) * D - R0y * h)) / (static_cast<T>(2) * D + R0y * h);
    T u11_c_term2_num = (static_cast<T>(32) * D4 * nu - static_cast<T>(96) * D4 + static_cast<T>(96) * D4 * nu2 - static_cast<T>(32) * D4 * nu3 - static_cast<T>(48) * D3 * R0y * h - static_cast<T>(48) * D3 * Rx0 * h + static_cast<T>(16) * D3 * R0y * h * nu + static_cast<T>(16) * D3 * Rx0 * h * nu - static_cast<T>(24) * D2 * R0y * Rx0 * h2 + static_cast<T>(8) * D2 * R0y * Rx0 * h2 * nu);
    T u11_c_term2 = u11_c_term2_num / den_d_term1;
    T u11_c = u11_c_term1 + u11_c_term2 - static_cast<T>(8);

    T u21_c_num = (Rx0 * D * R0y2 * h3 + static_cast<T>(2) * D2 * R0y2 * h2 + static_cast<T>(4) * Rx0 * D2 * R0y * h2 - static_cast<T>(4) * D3 * R0y * h * nu2 + static_cast<T>(8) * D3 * R0y * h + static_cast<T>(4) * Rx0 * D3 * h - static_cast<T>(8) * D4 * nu2 + static_cast<T>(8) * D4);
    T u21_c = u21_c_num / den_d_term1 + static_cast<T>(1);

    T u00_c_term1 = (- static_cast<T>(8) * D2 * nu2 + static_cast<T>(4) * R0y * h * D * nu + static_cast<T>(8) * D2 + static_cast<T>(4) * R0y * h * D) / den_common;
    T u00_c_term2 = (static_cast<T>(2) * (- static_cast<T>(8) * D2 * nu2 + static_cast<T>(4) * Rx0 * h * D * nu + static_cast<T>(8) * D2 + static_cast<T>(4) * Rx0 * h * D)) / den_common;
    T u00_c_term3 = (static_cast<T>(16) * D * nu) / (static_cast<T>(2) * D + R0y * h);
    T u00_c_term4_num = (static_cast<T>(32) * D4 * nu + static_cast<T>(32) * D4 - static_cast<T>(48) * D4 * nu2 - static_cast<T>(32) * D4 * nu3 + static_cast<T>(16) * D4 * nu4 + static_cast<T>(16) * D3 * R0y * h + static_cast<T>(16) * D3 * Rx0 * h + static_cast<T>(16) * D3 * R0y * h * nu + static_cast<T>(32) * D3 * Rx0 * h * nu + static_cast<T>(8) * D2 * R0y * Rx0 * h2 - static_cast<T>(24) * D3 * R0y * h * nu2 + static_cast<T>(8) * D3 * R0y * h * nu3 - static_cast<T>(16) * D3 * Rx0 * h * nu2 - static_cast<T>(8) * D2 * R0y * Rx0 * h2 * nu2 + static_cast<T>(16) * D2 * R0y * Rx0 * h2 * nu);
    T u00_c_term4 = u00_c_term4_num / den_d_term1;
    T u00_c = u00_c_term1 + u00_c_term2 + u00_c_term3 - u00_c_term4 - static_cast<T>(8);

    T u02_c_term1 = (static_cast<T>(2) * (static_cast<T>(4) * D + static_cast<T>(4) * D * nu)) / (static_cast<T>(2) * D + R0y * h);
    T u02_c_term2 = (static_cast<T>(16) * D * nu) / (static_cast<T>(2) * D + R0y * h);
    T u02_c_term3_num = (static_cast<T>(64) * D4 * nu + static_cast<T>(32) * D4 - static_cast<T>(64) * D4 * nu2 - static_cast<T>(64) * D4 * nu3 + static_cast<T>(32) * D4 * nu4 + static_cast<T>(16) * D3 * R0y * h + static_cast<T>(16) * D3 * Rx0 * h + static_cast<T>(32) * D3 * R0y * h * nu + static_cast<T>(32) * D3 * Rx0 * h * nu + static_cast<T>(8) * D2 * R0y * Rx0 * h2 - static_cast<T>(16) * D3 * R0y * h * nu2 - static_cast<T>(16) * D3 * Rx0 * h * nu2 - static_cast<T>(8) * D2 * R0y * Rx0 * h2 * nu2 + static_cast<T>(16) * D2 * R0y * Rx0 * h2 * nu);
    T u02_c_term3 = u02_c_term3_num / den_d_term1;
    T u02_c = u02_c_term1 + u02_c_term2 - u02_c_term3 - static_cast<T>(8);

    T u03_c_num = (static_cast<T>(16) * D4 * nu - static_cast<T>(8) * D4 * nu2 - static_cast<T>(16) * D4 * nu3 + static_cast<T>(8) * D4 * nu4 + static_cast<T>(8) * D3 * R0y * h * nu + static_cast<T>(8) * D3 * Rx0 * h * nu - static_cast<T>(4) * D3 * R0y * h * nu2 - static_cast<T>(4) * D3 * Rx0 * h * nu2 - static_cast<T>(2) * D2 * R0y * Rx0 * h2 * nu2 + static_cast<T>(4) * D2 * R0y * Rx0 * h2 * nu);
    T u03_c_den = D * (static_cast<T>(2) * D + R0y * h) * den_common;
    T u03_c = u03_c_num / u03_c_den - (static_cast<T>(4) * D * nu) / (static_cast<T>(2) * D + R0y * h) + static_cast<T>(1);

    T u12_c_num = (static_cast<T>(16) * D4 * nu - static_cast<T>(32) * D4 + static_cast<T>(32) * D4 * nu2 - static_cast<T>(16) * D4 * nu3 - static_cast<T>(16) * D3 * R0y * h - static_cast<T>(16) * D3 * Rx0 * h + static_cast<T>(8) * D3 * R0y * h * nu + static_cast<T>(8) * D3 * Rx0 * h * nu - static_cast<T>(8) * D2 * R0y * Rx0 * h2 + static_cast<T>(4) * D2 * R0y * Rx0 * h2 * nu);
    T u12_c_den = D * (static_cast<T>(2) * D + R0y * h) * den_common;
    T u12_c = static_cast<T>(2) - u12_c_num / u12_c_den - (static_cast<T>(2) * (static_cast<T>(2) * D - R0y * h)) / (static_cast<T>(2) * D + R0y * h);

    T u10_c_term1 = (static_cast<T>(2) * (static_cast<T>(4) * D2 * nu2 - static_cast<T>(4) * D2 + static_cast<T>(2) * D * R0y * h - static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2)) / den_common;
    T u10_c_term2_num = (static_cast<T>(16) * D4 * nu - static_cast<T>(32) * D4 + static_cast<T>(32) * D4 * nu2 - static_cast<T>(16) * D4 * nu3 - static_cast<T>(16) * D3 * R0y * h - static_cast<T>(16) * D3 * Rx0 * h + static_cast<T>(8) * D3 * R0y * h * nu + static_cast<T>(8) * D3 * Rx0 * h * nu - static_cast<T>(8) * D2 * R0y * Rx0 * h2 + static_cast<T>(16) * D3 * R0y * h * nu2 - static_cast<T>(8) * D3 * R0y * h * nu3 + static_cast<T>(4) * D2 * R0y * Rx0 * h2 * nu);
    T u10_c_term2 = u10_c_term2_num / den_d_term1;
    T u10_c_term3 = (static_cast<T>(4) * D * R0y * h * nu) / den_common;
    T u10_c = u10_c_term1 - u10_c_term2 - u10_c_term3 + static_cast<T>(2);

    return std::make_tuple(u01_c, u11_c, u21_c, u00_c, u02_c, u03_c, u12_c, u10_c);
}

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T>
D02_coeffs(T K0y, T R0y, T h, T D, T nu) {
    T D2 = D * D;
    T h3 = h * h * h;
    T h4 = h3 * h;
    T nu2 = nu * nu;
    T den = D * (static_cast<T>(2) * D + R0y * h);

    T u02_c_term1 = (static_cast<T>(2) * K0y * R0y * h4 + static_cast<T>(4) * K0y * D * h3 - static_cast<T>(12) * D2 * nu2 + static_cast<T>(24) * D2 * nu + static_cast<T>(24) * D2) / den;
    T u02_c_term2 = (static_cast<T>(8) * D * nu) / (static_cast<T>(2) * D + R0y * h);
    T u02_c_term3 = (static_cast<T>(8) * (static_cast<T>(4) * D + static_cast<T>(4) * D * nu)) / (static_cast<T>(2) * D + R0y * h);
    T u02_c = u02_c_term1 - u02_c_term2 - u02_c_term3 + static_cast<T>(20);

    T u12_c_term1 = (static_cast<T>(8) * (static_cast<T>(2) * D - R0y * h)) / (static_cast<T>(2) * D + R0y * h);
    T u12_c_term2 = (static_cast<T>(8) * D2 * nu - static_cast<T>(24) * D2) / den;
    T u12_c = u12_c_term1 + u12_c_term2 - static_cast<T>(8);

    T u22_c = (static_cast<T>(2) * D2 + R0y * h * D) / den + static_cast<T>(1);

    T u01_c_term1 = (static_cast<T>(2) * (static_cast<T>(4) * D + static_cast<T>(4) * D * nu)) / (static_cast<T>(2) * D + R0y * h);
    T u01_c_term2 = (- static_cast<T>(8) * D2 * nu2 + static_cast<T>(16) * D2 * nu + static_cast<T>(8) * D2) / den;
    T u01_c_term3 = (static_cast<T>(16) * D * nu) / (static_cast<T>(2) * D + R0y * h);
    T u01_c = u01_c_term1 - u01_c_term2 + u01_c_term3 - static_cast<T>(8);
    T u03_c = u01_c;

    T u04_c_term1 = (- static_cast<T>(2) * D2 * nu2 + static_cast<T>(4) * D2 * nu) / den;
    T u04_c_term2 = (static_cast<T>(4) * D * nu) / (static_cast<T>(2) * D + R0y * h);
    T u04_c = u04_c_term1 - u04_c_term2 + static_cast<T>(1);
    T u00_c = u04_c;

    T u13_c_term1 = (static_cast<T>(4) * D2 * nu - static_cast<T>(8) * D2) / den;
    T u13_c_term2 = (static_cast<T>(2) * (static_cast<T>(2) * D - R0y * h)) / (static_cast<T>(2) * D + R0y * h);
    T u13_c = static_cast<T>(2) - u13_c_term1 - u13_c_term2;
    T u11_c = u13_c;

    return std::make_tuple(u02_c, u12_c, u22_c, u01_c, u03_c, u04_c, u00_c, u13_c, u11_c);
}

template <typename T>
std::tuple<T, T, T, T, T, T, T, T>
D10_coeffs(T R0y, T Kx0, T Rx0, T h, T D, T nu) {
    // Direct implementation from Python D10_coeffs
    T D2 = D * D; T D3 = D2 * D; T D4 = D3 * D;
    T h2 = h * h; T h3 = h2 * h; T h4 = h3 * h; T h5 = h4 * h; T h6 = h5 * h;
    T nu2 = nu * nu; T nu3 = nu2 * nu; T nu4 = nu3 * nu;
    T Rx02 = Rx0 * Rx0;

    T den_common = (static_cast<T>(4) * D2 - static_cast<T>(4) * D2 * nu2 + static_cast<T>(2) * D * R0y * h + static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2);
    T den_d_term2 = D * (static_cast<T>(2) * D + Rx0 * h) * den_common;

    T u10_c = ((static_cast<T>(4) * D2 * nu2 - static_cast<T>(4) * D2 + static_cast<T>(2) * D * R0y * h - static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2) / den_common - 
               (static_cast<T>(8) * (static_cast<T>(4) * D + static_cast<T>(4) * D * nu)) / (static_cast<T>(2) * D + Rx0 * h) - 
               (static_cast<T>(4) * D * nu) / (static_cast<T>(2) * D + Rx0 * h) + 
               (static_cast<T>(2) * Kx0 * R0y * Rx02 * h6 + static_cast<T>(4) * Kx0 * D * Rx02 * h5 + static_cast<T>(8) * Kx0 * R0y * D * Rx0 * h5 - static_cast<T>(8) * Kx0 * D2 * Rx0 * h4 * nu2 + static_cast<T>(16) * Kx0 * D2 * Rx0 * h4 - static_cast<T>(14) * R0y * D2 * Rx0 * h2 * nu2 + static_cast<T>(28) * R0y * D2 * Rx0 * h2 * nu + static_cast<T>(24) * R0y * D2 * Rx0 * h2 - static_cast<T>(20) * D3 * Rx0 * h * nu2 + static_cast<T>(40) * D3 * Rx0 * h * nu + static_cast<T>(48) * D3 * Rx0 * h + static_cast<T>(8) * Kx0 * R0y * D2 * h4 - static_cast<T>(16) * Kx0 * D3 * h3 * nu2 + static_cast<T>(16) * Kx0 * D3 * h3 - static_cast<T>(28) * R0y * D3 * h * nu2 + static_cast<T>(56) * R0y * D3 * h * nu + static_cast<T>(48) * R0y * D3 * h + static_cast<T>(40) * D4 * nu4 - static_cast<T>(80) * D4 * nu3 - static_cast<T>(136) * D4 * nu2 + static_cast<T>(80) * D4 * nu + static_cast<T>(96) * D4) / den_d_term2 - 
               (static_cast<T>(8) * D * R0y * h * nu) / den_common + static_cast<T>(20));

    T u20_c = ((static_cast<T>(2) * (static_cast<T>(4) * D + static_cast<T>(4) * D * nu)) / (static_cast<T>(2) * D + Rx0 * h) + 
               (static_cast<T>(16) * D * nu) / (static_cast<T>(2) * D + Rx0 * h) - 
               (static_cast<T>(64) * D4 * nu + static_cast<T>(32) * D4 - static_cast<T>(64) * D4 * nu2 - static_cast<T>(64) * D4 * nu3 + static_cast<T>(32) * D4 * nu4 + static_cast<T>(16) * D3 * R0y * h + static_cast<T>(16) * D3 * Rx0 * h + static_cast<T>(32) * D3 * R0y * h * nu + static_cast<T>(32) * D3 * Rx0 * h * nu + static_cast<T>(8) * D2 * R0y * Rx0 * h2 - static_cast<T>(16) * D3 * R0y * h * nu2 - static_cast<T>(16) * D3 * Rx0 * h * nu2 - static_cast<T>(8) * D2 * R0y * Rx0 * h2 * nu2 + static_cast<T>(16) * D2 * R0y * Rx0 * h2 * nu) / den_d_term2 - static_cast<T>(8));

    T u30_c = ((static_cast<T>(16) * D4 * nu - static_cast<T>(8) * D4 * nu2 - static_cast<T>(16) * D4 * nu3 + static_cast<T>(8) * D4 * nu4 + static_cast<T>(8) * D3 * R0y * h * nu + static_cast<T>(8) * D3 * Rx0 * h * nu - static_cast<T>(4) * D3 * R0y * h * nu2 - static_cast<T>(4) * D3 * Rx0 * h * nu2 - static_cast<T>(2) * D2 * R0y * Rx0 * h2 * nu2 + static_cast<T>(4) * D2 * R0y * Rx0 * h2 * nu) / den_d_term2 - 
               (static_cast<T>(4) * D * nu) / (static_cast<T>(2) * D + Rx0 * h) + static_cast<T>(1));

    T u00_c = ((static_cast<T>(2) * (- static_cast<T>(8) * D2 * nu2 + static_cast<T>(4) * R0y * h * D * nu + static_cast<T>(8) * D2 + static_cast<T>(4) * R0y * h * D)) / den_common + 
               (- static_cast<T>(8) * D2 * nu2 + static_cast<T>(4) * Rx0 * h * D * nu + static_cast<T>(8) * D2 + static_cast<T>(4) * Rx0 * h * D) / den_common + 
               (static_cast<T>(16) * D * nu) / (static_cast<T>(2) * D + Rx0 * h) - 
               (static_cast<T>(32) * D4 * nu + static_cast<T>(32) * D4 - static_cast<T>(48) * D4 * nu2 - static_cast<T>(32) * D4 * nu3 + static_cast<T>(16) * D4 * nu4 + static_cast<T>(16) * D3 * R0y * h + static_cast<T>(16) * D3 * Rx0 * h + static_cast<T>(32) * D3 * R0y * h * nu + static_cast<T>(16) * D3 * Rx0 * h * nu + static_cast<T>(8) * D2 * R0y * Rx0 * h2 - static_cast<T>(16) * D3 * R0y * h * nu2 - static_cast<T>(24) * D3 * Rx0 * h * nu2 + static_cast<T>(8) * D3 * Rx0 * h * nu3 - static_cast<T>(8) * D2 * R0y * Rx0 * h2 * nu2 + static_cast<T>(16) * D2 * R0y * Rx0 * h2 * nu) / den_d_term2 - static_cast<T>(8));

    T u11_c = ((static_cast<T>(8) * (static_cast<T>(2) * D - Rx0 * h)) / (static_cast<T>(2) * D + Rx0 * h) + 
               (static_cast<T>(32) * D4 * nu - static_cast<T>(96) * D4 + static_cast<T>(96) * D4 * nu2 - static_cast<T>(32) * D4 * nu3 - static_cast<T>(48) * D3 * R0y * h - static_cast<T>(48) * D3 * Rx0 * h + static_cast<T>(16) * D3 * R0y * h * nu + static_cast<T>(16) * D3 * Rx0 * h * nu - static_cast<T>(24) * D2 * R0y * Rx0 * h2 + static_cast<T>(8) * D2 * R0y * Rx0 * h2 * nu) / den_d_term2 - static_cast<T>(8));

    T u12_c = ((R0y * D * Rx02 * h3 + static_cast<T>(2) * D2 * Rx02 * h2 + static_cast<T>(4) * R0y * D2 * Rx0 * h2 - static_cast<T>(4) * D3 * Rx0 * h * nu2 + static_cast<T>(8) * D3 * Rx0 * h + static_cast<T>(4) * R0y * D3 * h - static_cast<T>(8) * D4 * nu2 + static_cast<T>(8) * D4) / den_d_term2 + static_cast<T>(1));

    T u21_c = (static_cast<T>(2) - 
               (static_cast<T>(16) * D4 * nu - static_cast<T>(32) * D4 + static_cast<T>(32) * D4 * nu2 - static_cast<T>(16) * D4 * nu3 - static_cast<T>(16) * D3 * R0y * h - static_cast<T>(16) * D3 * Rx0 * h + static_cast<T>(8) * D3 * R0y * h * nu + static_cast<T>(8) * D3 * Rx0 * h * nu - static_cast<T>(8) * D2 * R0y * Rx0 * h2 + static_cast<T>(4) * D2 * R0y * Rx0 * h2 * nu) / den_d_term2 - 
               (static_cast<T>(2) * (static_cast<T>(2) * D - Rx0 * h)) / (static_cast<T>(2) * D + Rx0 * h));

    T u01_c = ((static_cast<T>(2) * (static_cast<T>(4) * D2 * nu2 - static_cast<T>(4) * D2 - static_cast<T>(2) * D * R0y * h + static_cast<T>(2) * D * Rx0 * h + R0y * Rx0 * h2)) / den_common - 
               (static_cast<T>(16) * D4 * nu - static_cast<T>(32) * D4 + static_cast<T>(32) * D4 * nu2 - static_cast<T>(16) * D4 * nu3 - static_cast<T>(16) * D3 * R0y * h - static_cast<T>(16) * D3 * Rx0 * h + static_cast<T>(8) * D3 * R0y * h * nu + static_cast<T>(8) * D3 * Rx0 * h * nu - static_cast<T>(8) * D2 * R0y * Rx0 * h2 + static_cast<T>(16) * D3 * Rx0 * h * nu2 - static_cast<T>(8) * D3 * Rx0 * h * nu3 + static_cast<T>(4) * D2 * R0y * Rx0 * h2 * nu) / den_d_term2 - 
               (static_cast<T>(4) * D * Rx0 * h * nu) / den_common + static_cast<T>(2));

    // Return order matches Python D10_coeffs: [u10_c, u20_c, u30_c, u00_c, u11_c, u12_c, u21_c, u01_c]
    // Corresponds to C++ binding: [D10u10, D10u20, D10u30, D10u00, D10u11, D10u12, D10u21, D10u01]
    return std::make_tuple(u10_c, u20_c, u30_c, u00_c, u11_c, u12_c, u21_c, u01_c);
}

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T, T, T>
D11_coeffs(T R0y, T Rx0, T h, T D, T nu) {
    T den1 = (static_cast<T>(2) * D + Rx0 * h);
    T den2 = (static_cast<T>(2) * D + R0y * h);
    T val1 = static_cast<T>(20) - (static_cast<T>(2) * D - Rx0 * h) / den1 - (static_cast<T>(2) * D - R0y * h) / den2;
    T val2 = static_cast<T>(-8);
    T val3 = static_cast<T>(1);
    T val4 = (static_cast<T>(4) * D + static_cast<T>(4) * D * nu) / den1 - static_cast<T>(8);
    T val5 = (static_cast<T>(4) * D + static_cast<T>(4) * D * nu) / den2 - static_cast<T>(8);
    T val6 = static_cast<T>(-8);
    T val7 = static_cast<T>(1);
    T val8 = static_cast<T>(2);
    T val9 = static_cast<T>(2) - (static_cast<T>(2) * D * nu) / den1;
    T val10 = static_cast<T>(2) - (static_cast<T>(2) * D * nu) / den1 - (static_cast<T>(2) * D * nu) / den2;
    T val11 = static_cast<T>(2) - (static_cast<T>(2) * D * nu) / den2;
    return std::make_tuple(val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11);
}

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T, T, T, T>
D12_coeffs(T R0y, T h, T D, T nu) {
    T den = (static_cast<T>(2) * D + R0y * h);
    T val1 = static_cast<T>(20) - (static_cast<T>(2) * D - R0y * h) / den;
    T val2 = static_cast<T>(-8);
    T val3 = static_cast<T>(1);
    T val4 = static_cast<T>(-8);
    T val5 = static_cast<T>(1);
    T val6 = (static_cast<T>(4) * D + static_cast<T>(4) * D * nu) / den - static_cast<T>(8);
    T val7 = static_cast<T>(-8);
    T val8 = static_cast<T>(1);
    T val9 = static_cast<T>(2);
    T val10 = static_cast<T>(2);
    T val11 = static_cast<T>(2) - (static_cast<T>(2) * D * nu) / den;
    T val12 = val11;
    return std::make_tuple(val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12);
}

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T>
D20_coeffs(T Kx0, T Rx0, T h, T D, T nu) {
    // This is symmetric to D02_coeffs, swap K0y<->Kx0, R0y<->Rx0
    // Assuming this symmetry holds and the return order is correct based on D02_coeffs
    return D02_coeffs<T>(Kx0, Rx0, h, D, nu);
}

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T, T, T, T>
D21_coeffs(T Rx0, T h, T D, T nu) {
    // Direct implementation from Python D21_coeffs
    T den = (static_cast<T>(2) * D + Rx0 * h);
    // Consider adding check for den == 0 if necessary based on expected inputs
    // if (den == static_cast<T>(0)) { throw std::runtime_error("Division by zero in D21_coeffs"); }
    
    T val1 = static_cast<T>(20) - (static_cast<T>(2) * D - Rx0 * h) / den; // D21u21
    T val2 = static_cast<T>(-8);                                        // D21u22
    T val3 = static_cast<T>(1);                                         // D21u23
    T val4 = (static_cast<T>(4) * D + static_cast<T>(4) * D * nu) / den - static_cast<T>(8); // D21u20
    T val5 = static_cast<T>(-8);                                        // D21u11 
    T val6 = static_cast<T>(-8);                                        // D21u31 
    T val7 = static_cast<T>(1);                                         // D21u41 
    T val8 = static_cast<T>(1);                                         // D21u01 
    T val9 = static_cast<T>(2);                                         // D21u32 
    T val10 = static_cast<T>(2) - (static_cast<T>(2) * D * nu) / den;    // D21u30 
    T val11 = static_cast<T>(2) - (static_cast<T>(2) * D * nu) / den;    // D21u10 
    T val12 = static_cast<T>(2);                                        // D21u12 

    // Return order matches Python D21_coeffs: [u21, u22, u23, u20, u11, u31, u41, u01, u32, u30, u10, u12]
    // Corresponds to C++ binding: [D21u21, D21u22, D21u23, D21u20, D21u11, D21u31, D21u41, D21u01, D21u32, D21u30, D21u10, D21u12]
    return std::make_tuple(val1, val2, val3, val4, val5, val6, val7, val8, val9, val10, val11, val12);
}

template <typename T>
std::tuple<T, T, T, T, T, T, T, T, T, T, T, T, T>
D22_coeffs() {
    // Return order matches Python D22_coeffs: [c20_22, c11_22, c21_22, c31_22, c02_22, c12_22, c22_22, c32_22, c42_22, c13_22, c23_22, c33_22, c24_22]
    return std::make_tuple(static_cast<T>(1), static_cast<T>(2), static_cast<T>(-8), static_cast<T>(2), static_cast<T>(1), static_cast<T>(-8), static_cast<T>(20), static_cast<T>(-8), static_cast<T>(1), static_cast<T>(2), static_cast<T>(-8), static_cast<T>(2), static_cast<T>(1));
} 