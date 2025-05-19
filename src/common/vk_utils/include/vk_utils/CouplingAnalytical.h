#pragma once

#include <vector>
#include <complex>
#include <map>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <functional> // For std::reference_wrapper

// Additional includes from coupling.cpp
#include <cmath>
#include <numeric>    // for std::iota
#include <algorithm>  // for std::sort, std::stable_sort
#include <iostream>   // For cerr, cout used in some functions
// <Eigen/src/Core/Matrix.h> is generally covered by <Eigen/Dense>

// Define type aliases for clarity
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using MatrixCD = Eigen::MatrixXcd;
using VectorCD = Eigen::VectorXcd;
using Tensor3D = std::vector<Matrix>; // Representing Npsi x Nphi x Nphi

// Helper function for (-1)^n
inline double pow_m1(int n)
{
    return (n % 2 == 0) ? 1.0 : -1.0;
}

//----------------------------------------------------------------------------
// Basic Integrals
//----------------------------------------------------------------------------

inline double int4(int m, int p, double L)
{
    double y = 0.0;
    double pi = M_PI;
    double pi2 = pi * pi;
    double pi4 = pi2 * pi2;
    double m_dbl = static_cast<double>(m);
    double p_dbl = static_cast<double>(p);
    double m2 = m_dbl * m_dbl;
    double p2 = p_dbl * p_dbl;
    double m1_m = pow_m1(m);
    double m1_p = pow_m1(p);

    if (m == 0 && p == 0) { y = 120.0 / (7.0 * L); }
    else if (m == p && m != 0)
    {
        y = (768.0 * pi2 * m2 - 47040.0 * m1_m + 35.0 * pi4 * m2 * m2 +
             432.0 * m1_m * pi2 * m2 - 53760.0) /
            (70.0 * L * pi2 * m2);
    }
    else if (m == 0)
    {
        double num = 60.0 * (m1_p + 1.0) * (pi2 * p2 - 42.0);
        double den = 7.0 * L * pi2 * p2;
        y = (den == 0) ? 0.0 : num / den;
    }
    else if (p == 0)
    {
        double num = 60.0 * (m1_m + 1.0) * (pi2 * m2 - 42.0);
        double den = 7.0 * L * pi2 * m2;
        y = (den == 0) ? 0.0 : num / den;
    }
    else
    {
        double term1 = 192.0 / 35.0 / L * (1.0 + m1_m * m1_p);
        double term2 =
            -192.0 / (m2 * p2 * L * pi2) * ((p2 + m2) * (1.0 + m1_m * m1_p));
        double term3 =
            -168.0 / (m2 * p2 * L * pi2) * ((p2 + m2) * (m1_m + m1_p));
        double term4 = 108.0 / 35.0 / L * (m1_m + m1_p);
        y = term1 + term2 + term3 + term4;
    }
    return y;
}

inline double int1(int m, int p, double L)
{
    double y = 0.0;
    double L3 = L * L * L;
    double pi4 = M_PI * M_PI * M_PI * M_PI;
    double m_dbl = static_cast<double>(m);
    double m4 = m_dbl * m_dbl * m_dbl * m_dbl;
    double m1_m = pow_m1(m);
    double m1_p = pow_m1(p);

    if (m == 0 && p == 0) { y = 720.0 / L3; }
    else if (m == p && m != 0)
    {
        y = (pi4 * m4 - 672.0 * m1_m - 768.0) / (2.0 * L3);
    }
    else if (m == 0 || p == 0) { y = 0.0; }
    else
    {
        double val = 7.0 * m1_m + 7.0 * m1_p + 8.0 * m1_m * m1_p + 8.0;
        y = -24.0 * val / L3;
    }
    return y;
}

inline double int2(int m, int p, double L)
{
    double y = 0.0;
    double pi = M_PI;
    double pi4 = pi * pi * pi * pi;
    double m_dbl = static_cast<double>(m);
    double p_dbl = static_cast<double>(p);
    double m2 = m_dbl * m_dbl;
    double p2 = p_dbl * p_dbl;
    double m4 = m2 * m2;
    double p4 = p2 * p2;
    double m1_m = pow_m1(m);
    double m1_p = pow_m1(p);

    if (m == 0 && p == 0) { y = (10.0 * L) / 7.0; }
    else if (m == p && m != 0)
    {
        if (pi4 * m4 == 0) return 0.0;
        y = (67.0 * L) / 70.0 - (m1_m * L) / 35.0 - (768.0 * L) / (pi4 * m4) -
            (672.0 * m1_m * L) / (pi4 * m4);
    }
    else if (m == 0)
    {
        if (pi4 * p4 == 0) return 0.0;
        double num = 3.0 * L * (m1_p + 1.0) * (pi4 * p4 - 1680.0);
        double den = 14.0 * pi4 * p4;
        y = num / den;
    }
    else if (p == 0)
    {
        if (pi4 * m4 == 0) return 0.0;
        double num = 3.0 * L * (m1_m + 1.0) * (pi4 * m4 - 1680.0);
        double den = 14.0 * pi4 * m4;
        y = num / den;
    }
    else
    {
        if (pi4 * m4 == 0 || p4 == 0) return 0.0;
        double part1 = 11760.0 * m1_m + 11760.0 * m1_p - 16.0 * pi4 * m4 +
                       13440.0 * m1_m * m1_p + m1_m * pi4 * m4 +
                       m1_p * pi4 * m4 - 16.0 * m1_m * m1_p * pi4 * m4 +
                       13440.0;
        double part2 = 13440.0 * m4 + 11760.0 * m1_m * m4 +
                       11760.0 * m1_p * m4 + 13440.0 * m1_m * m1_p * m4;
        double den1 = 70.0 * pi4 * m4;
        double den2 = 70.0 * pi4 * m4 * p4;
        y = -(L * part1) / den1 - (L * part2) / den2;
    }
    return y;
}

//----------------------------------------------------------------------------
// Integral Matrix Builders
//----------------------------------------------------------------------------

inline Matrix build_I1(int N, double L)
{
    Matrix I = Matrix::Zero(N, N);
    for (int m = 0; m < N; ++m)
    {
        for (int p = 0; p < N; ++p) { I(m, p) = int1(m, p, L); }
    }
    return I;
}

inline Matrix build_I2(int N, double L)
{
    Matrix I = Matrix::Zero(N, N);
    for (int m = 0; m < N; ++m)
    {
        for (int p = 0; p < N; ++p) { I(m, p) = int2(m, p, L); }
    }
    return I;
}

inline Matrix build_I4(int N, double L)
{
    Matrix I = Matrix::Zero(N, N);
    for (int m = 0; m < N; ++m)
    {
        for (int p = 0; p < N; ++p) { I(m, p) = int4(m, p, L); }
    }
    return I;
}

inline Matrix int2_mat(int N, double L)
{
    return build_I2(N, L);
}

//----------------------------------------------------------------------------
// Kronecker Product Helper
//----------------------------------------------------------------------------

inline Matrix kron(const Matrix& A, const Matrix& B)
{
    Matrix C(A.rows() * B.rows(), A.cols() * B.cols());
    for (int i = 0; i < A.rows(); ++i)
    {
        for (int j = 0; j < A.cols(); ++j)
        {
            C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) =
                A(i, j) * B;
        }
    }
    return C;
}

//----------------------------------------------------------------------------
// K and M Matrix Assembly
//----------------------------------------------------------------------------

inline std::pair<Matrix, Matrix> assemble_K_and_M(int Npsi, double Lx, double Ly)
{
    Matrix I1x = build_I1(Npsi, Lx);
    Matrix I2x = build_I2(Npsi, Lx);
    Matrix I4x = build_I4(Npsi, Lx);

    Matrix I1y = build_I1(Npsi, Ly);
    Matrix I2y = build_I2(Npsi, Ly);
    Matrix I4y = build_I4(Npsi, Ly);

    Matrix K = kron(I1x, I2y) + kron(I2x, I1y) + 2.0 * kron(I4x, I4y);
    Matrix M = kron(I2x, I2y);

    return {K, M};
}

//----------------------------------------------------------------------------
// Airy Stress Coefficients
//----------------------------------------------------------------------------

struct AiryCoefficients {
    Matrix coeff0;
    Matrix coeff1;
    Matrix coeff2;
    int S; // Number of valid modes used
    Vector auto_vec; // Eigenvalues corresponding to coeffs
    Matrix coeff_sorted;
};

inline AiryCoefficients airy_stress_coefficients(
    int n_psi,
    double Lx,
    double Ly,
    const Vector& vals,
    const Matrix& vecs
)
{
    int dim_total = vals.size();
    if (vecs.cols() != dim_total)
    {
        throw std::runtime_error("Eigenvalues and Eigenvectors size mismatch."
        );
    }

    std::vector<int> good_indices;
    for (int i = 0; i < dim_total; ++i)
    {
        if (vals(i) >= 0.0) { good_indices.push_back(i); }
    }

    int S_good = good_indices.size();
    Vector good_vals_real(S_good);
    Matrix good_vecs(vecs.rows(), S_good);

    for (int i = 0; i < S_good; ++i)
    {
        good_vals_real(i) = vals(good_indices[i]);
        good_vecs.col(i) = vecs.col(good_indices[i]);
    }

    std::vector<int> idx_sort(S_good);
    std::iota(idx_sort.begin(), idx_sort.end(), 0);

    std::stable_sort(
        idx_sort.begin(),
        idx_sort.end(),
        [&good_vals_real](int i1, int i2)
        { return good_vals_real(i1) < good_vals_real(i2); }
    );

    Vector auto_sorted(S_good);
    Matrix coeff_sorted(vecs.rows(), S_good);
    for (int i = 0; i < S_good; ++i)
    {
        auto_sorted(i) = good_vals_real(idx_sort[i]);
        coeff_sorted.col(i) = good_vecs.col(idx_sort[i]);
    }

    int dim = n_psi * n_psi;
    if (dim != vecs.rows())
    {
        throw std::runtime_error(
            "Dimension mismatch in airy_stress_coefficients."
        );
    }

    Matrix coeff0 = Matrix::Zero(dim, S_good);
    Matrix coeff1 = Matrix::Zero(dim, S_good);
    Matrix coeff2 = Matrix::Zero(dim, S_good);

    const int N2 = n_psi * n_psi;
    const int N4 = N2 * N2; // This was N2 * N2, which is n_psi^4. dim = n_psi^2
    Matrix I2x_norm = int2_mat(n_psi, Lx);
    Matrix I2y_norm = int2_mat(n_psi, Ly);
#if 0 // Original normalization based on M_norm matrix
    Matrix M_norm = kron(I2x_norm, I2y_norm);
    for (int d = 0; d < S_good; ++d)
    {
        Vector temp = coeff_sorted.col(d);
        double norm_sq = temp.transpose() * M_norm * temp;
        double norm_val = (norm_sq > 0) ? std::sqrt(norm_sq) : 1.0;
        if (norm_val == 0) norm_val = 1.0;

        double auto_d = auto_sorted(d);
        double sqrt_auto_d = (auto_d > 0) ? std::sqrt(auto_d) : 1.0;

        coeff0.col(d) = temp / norm_val;
        coeff1.col(d) = temp / (norm_val * sqrt_auto_d);
        coeff2.col(d) = (auto_d != 0) ? (temp / (norm_val * auto_d)).eval()
                                      : Vector::Zero(dim);
    }
#else // Current normalization attempt (from user's cpp)
 
    // This section attempts to replicate a specific MATLAB normalization.
    // Ensure I2x_norm and I2y_norm are n_psi x n_psi.
    // NN should map to I2x_norm (n_psi x n_psi) -> vector of size n_psi*n_psi
    // MM should map to I2y_norm (n_psi x n_psi) -> vector of size n_psi*n_psi
    
    // If I2x_norm and I2y_norm are (n_psi x n_psi)
    // NN should be Eigen::Map<const Vector>(I2x_norm.data(), n_psi * n_psi);
    // MM should be Eigen::Map<const Vector>(I2y_norm.data(), n_psi * n_psi);
    // Then NNrep would be (n_psi*n_psi) x (n_psi*n_psi) if N2 = n_psi*n_psi.
    // This part seems to have dimensions mixed up with N2 (n_psi*n_psi) vs dim (n_psi*n_psi).
    // The original Python code's nmatr was N2xN2 where N2 = n_psi.
    // Let's assume N2 here refers to n_psi (from context of original Python where nmatr was derived from I2x, I2y which were n_psi x n_psi)
    // But the variable N2 is defined as n_psi * n_psi. This leads to confusion.
    // If I2x_norm is (n_psi x n_psi), then .data() gives n_psi*n_psi elements.
    // The Python: NN = I2x_norm.flatten('F'), MM = I2y_norm.flatten('F')
    // nmatr = (NN[:,None] @ MM[None,:]) which means nmatr(i,j) = NN[i]*MM[j]
    // So nmatr is (n_psi*n_psi) x (n_psi*n_psi)
    // Mvec was a reshape of nmatr to (n_psi^4, 1)

    Matrix M_norm_kron = kron(I2x_norm, I2y_norm); // This is (n_psi^2 x n_psi^2) = (dim x dim)

    for (int d = 0; d < S_good; ++d) {
        Vector temp = coeff_sorted.col(d); // temp is (dim x 1) or (n_psi^2 x 1)
        
        // The norm calculation: temp.transpose() * M_norm_kron * temp
        // M_norm_kron is (dim x dim).
        double norm_sq = temp.transpose() * M_norm_kron * temp;
        double norm_val = (norm_sq > 0.0) ? std::sqrt(norm_sq) : 1.0;
        if (norm_val == 0.0) norm_val = 1.0; // Prevent division by zero

        double auto_d = auto_sorted(d); // Corrected: use (d) not [d] for Eigen::Vector
        double sqrt_auto_d = (auto_d > 0.0 ? std::sqrt(auto_d) : 1.0);

        coeff0.col(d) = temp / norm_val;
        coeff1.col(d) = temp / (norm_val * sqrt_auto_d);
        coeff2.col(d) = (auto_d != 0) ? (temp / (norm_val * auto_d)).eval() : Vector::Zero(dim);
    }
#endif
    int S2 = S_good / 2; 
    return {
        coeff0.block(0, 0, S2, S2),
        coeff1.block(0, 0, S2, S2),
        coeff2.block(0, 0, S2, S2),
        S2,
        auto_sorted,
        coeff_sorted,
    };
}


//----------------------------------------------------------------------------
// Partial Integral Tensors (i*_mat) - Helper Compute Functions
//----------------------------------------------------------------------------

inline double compute_i1(int m1, int n, int p, double L, double) // pi unused
{
    if (m1 == 0 && n == p) return L / 2.0;
    if (m1 == (p - n) || m1 == (n - p)) return L / 4.0;
    if (m1 == (-n - p) || m1 == (n + p)) return -L / 4.0;
    return 0.0;
}

inline double compute_i2(int m1, int n, int p, double L, double pi)
{
    double m1_term = pow_m1(m1) + 1.0;
    double L2 = L * L;
    double L4 = L2 * L2;
    double L5 = L4 * L;
    double pi2 = pi * pi;
    double pi3 = pi2 * pi;
    double pi5 = pi3 * pi2;
    double n_dbl = static_cast<double>(n);
    double p_dbl = static_cast<double>(p);
    double n2 = n_dbl * n_dbl;
    double p2 = p_dbl * p_dbl;
    double n3 = n2 * n_dbl;
    double p3 = p2 * p_dbl;
    double n5 = n3 * n2;
    double p5 = p3 * p2;

    if (n == p)
    {
        double num =
            L5 * (4.0 * pi5 * p5 - 20.0 * pi3 * p3 + 30.0 * pi * p_dbl);
        double den = 40.0 * pi5 * p5;
        return (den == 0) ? 0.0 : (15.0 / L4 * m1_term) * (num / den);
    }
    else
    {
        double np_sum = n_dbl + p_dbl;
        double np_diff = n_dbl - p_dbl;
        double np_sum2 = np_sum * np_sum;
        double np_sum4 = np_sum2 * np_sum2;
        double np_sum5 = np_sum4 * np_sum;
        double np_diff2 = np_diff * np_diff;
        double np_diff4 = np_diff2 * np_diff2;
        double np_diff5 = np_diff4 * np_diff;

        if (np_sum5 == 0 || np_diff5 == 0)
            return 0.0;

        double term1_num =
            sin(pi * np_sum) *
            ((1713638851887625.0 / 17592186044416.0) * L4 * np_sum4 -
             (8334140006820045.0 / 70368744177664.0) * L4 * np_sum2 +
             24.0 * L4);
        term1_num +=
            4.0 * pi * L2 * cos(pi * np_sum) * np_sum *
            ((2778046668940015.0 / 281474976710656.0) * L2 * np_sum2 -
             6.0 * L2);
        double term1 = term1_num / np_sum5;

        double term2_num =
            sin(pi * np_diff) *
            ((1713638851887625.0 / 17592186044416.0) * L4 * np_diff4 -
             (8334140006820045.0 / 70368744177664.0) * L4 * np_diff2 +
             24.0 * L4);
        term2_num +=
            4.0 * pi * L2 * cos(pi * np_diff) * np_diff *
            ((2778046668940015.0 / 281474976710656.0) * L2 * np_diff2 -
             6.0 * L2);
        double term2 = term2_num / np_diff5;

        double big_term = 8796093022208.0 * L * (term1 - term2);
        return -(15.0 / L4 * m1_term) * big_term / 5383555227996211.0;
    }
}

inline double compute_i3(int m1, int n, int p, double L, double pi)
{
    double m1_term = 7.0 * pow_m1(m1) + 8.0;
    double L2 = L * L;
    double L3 = L2 * L;
    double L4 = L3 * L;
    double pi2 = pi * pi;
    double pi4 = pi2 * pi2;
    double n_dbl = static_cast<double>(n);
    double p_dbl = static_cast<double>(p);
    double n2 = n_dbl * n_dbl;
    double p2 = p_dbl * p_dbl;
    double n4 = n2 * n2;
    double p4 = p2 * p2;

    if (n == p)
    {
        double num = L4 * (6.0 * pi2 * p2 - 2.0 * pi4 * p4);
        double den = 16.0 * pi4 * p4;
        return (den == 0) ? 0.0 : (4.0 / L3 * m1_term) * (num / den);
    }
    else
    {
        double np_sum = n_dbl + p_dbl;
        double np_sum2 = np_sum * np_sum;
        double np_sum4 = np_sum2 * np_sum2;
        double np_diff = n_dbl - p_dbl;
        double np_diff2 = np_diff * np_diff;
        double np_diff4 = np_diff2 * np_diff2;

        if (np_sum4 == 0 || np_diff4 == 0 || pi4 == 0)
            return 0.0;  // Avoid division by zero

        double common_factor = (-4.0 / L3 * m1_term) / (2.0 * pi4);

        double term1_num = L * ((6.0 * L3) / np_diff4 - (6.0 * L3) / np_sum4);
        double term1 = common_factor * term1_num;

        double term2_cos_sum =
            3.0 * L * cos(pi * np_sum) * (2.0 * L2 - L2 * pi2 * np_sum2);
        double term2_cos_diff =
            3.0 * L * cos(pi * np_diff) * (2.0 * L2 - L2 * pi2 * np_diff2);
        double term2_num =
            L * (term2_cos_sum / np_sum4 - term2_cos_diff / np_diff4);
        double term2 = common_factor * term2_num;

        return term1 + term2;
    }
}

inline double compute_i4(int m1, int n, int p, double L, double pi)
{
    double m1_term = 2.0 * pow_m1(m1) + 3.0;
    double L2 = L * L;
    double L3 = L2 * L;
    double pi2 = pi * pi;
    double pi3 = pi2 * pi;
    double n_dbl = static_cast<double>(n);
    double p_dbl = static_cast<double>(p);
    double n2 = n_dbl * n_dbl;
    double p2 = p_dbl * p_dbl;
    double n3 = n2 * n_dbl;
    double p3 = p2 * p_dbl;

    if (n == p)
    {
        double num = L3 * (6.0 * pi * p_dbl - 4.0 * pi3 * p3);
        double den = 24.0 * pi3 * p3;
        return (den == 0) ? 0.0 : -(6.0 / L2 * m1_term) * (num / den);
    }
    else
    {
        double np_sum = n_dbl + p_dbl;
        double np_sum2 = np_sum * np_sum;
        double np_diff = n_dbl - p_dbl;
        double np_diff2 = np_diff * np_diff;

        if (pi2 == 0 || np_diff2 == 0 || np_sum2 == 0)
            return 0.0;  // Avoid division by zero

        double common_factor = (6.0 / L2 * m1_term) * L3 / pi2;
        double term1 = common_factor * cos(pi * np_diff) / np_diff2;
        double term2 = common_factor * cos(pi * np_sum) / np_sum2;
        return term1 - term2;
    }
}

inline double compute_i5(int, int n, int p, double L, double) // m1, pi unused
{
    return (n == p) ? -L / 2.0 : 0.0;
}

inline double compute_i9(int m1, int n, int p, double L, double) // pi unused
{
    if (m1 == 0 && n == p) return L / 2.0;
    if (m1 == (p - n) || m1 == (n - p)) return L / 4.0;
    if (m1 == (-n - p) || m1 == (n + p)) return L / 4.0;
    return 0.0;
}

inline double compute_i10(int m1, int n, int p, double L, double pi)
{
    double m1_term = pow_m1(m1) + 1.0;
    double L2 = L * L;
    double L4 = L2 * L2;
    double L5 = L4 * L;
    double pi2 = pi * pi;
    double pi3 = pi2 * pi;
    double pi5 = pi3 * pi2;
    double n_dbl = static_cast<double>(n);
    double p_dbl = static_cast<double>(p);
    double n2 = n_dbl * n_dbl;
    double p2 = p_dbl * p_dbl;
    double n3 = n2 * n_dbl;
    double p3 = p2 * p_dbl;
    double n4 = n3 * n_dbl;
    double p4 = p3 * p_dbl;
    double n5 = n4 * n_dbl;
    double p5 = p4 * p_dbl;

    if (n == p && n != 0)
    {
        double num =
            L5 * (4.0 * pi5 * n5 + 20.0 * pi3 * n3 - 30.0 * pi * n_dbl);
        double den = 40.0 * pi5 * n5;
        return (den == 0) ? 0.0 : (15.0 / L4 * m1_term) * (num / den);
    }
    else if (n != p)
    {
        double np_sum = n_dbl + p_dbl;
        double np_sum2 = np_sum * np_sum;
        double np_sum4 = np_sum2 * np_sum2;
        double np_diff = n_dbl - p_dbl;
        double np_diff2 = np_diff * np_diff;
        double np_diff4 = np_diff2 * np_diff2;

        if (np_sum4 == 0 || np_diff4 == 0 || pi5 == 0)
            return 0.0;  // Avoid division by zero

        double term1_num = 4.0 * pi * L2 * cos(pi * np_sum) *
                           (6.0 * L2 - L2 * pi2 * np_sum2);
        double term1 = term1_num / np_sum4;

        double term2_num = 4.0 * pi * L2 * cos(pi * np_diff) *
                           (6.0 * L2 - L2 * pi2 * np_diff2);
        double term2 = term2_num / np_diff4;

        double common_factor = -(15.0 / L4 * m1_term) * L / (2.0 * pi5);
        return common_factor * (term1 + term2);
    }
    else
    {
        return 0.0;  // Case n=p=0 or n!=p covered
    }
}

inline double compute_i11(int m1, int n, int p, double L, double pi)
{
    double m1_term = 7.0 * pow_m1(m1) + 8.0;
    double L2 = L * L;
    double L3 = L2 * L;
    double L4 = L3 * L;
    double pi2 = pi * pi;
    double pi4 = pi2 * pi2;
    double n_dbl = static_cast<double>(n);
    double p_dbl = static_cast<double>(p);
    double n2 = n_dbl * n_dbl;
    double p2 = p_dbl * p_dbl;
    double n4 = n2 * n2;
    double p4 = p2 * p2;

    if (n == p && n != 0)
    {
        if (pi2 * p2 == 0) return 0.0;
        return (-4.0 / L3 * m1_term) * L4 / 8.0 +
               (-4.0 / L3 * m1_term) * (3.0 * L4) / (8.0 * pi2 * p2);
    }
    else if (n != p)
    {
        double np_sum = n_dbl + p_dbl;
        double np_sum2 = np_sum * np_sum;
        double np_sum4 = np_sum2 * np_sum2;
        double np_diff = n_dbl - p_dbl;
        double np_diff2 = np_diff * np_diff;
        double np_diff4 = np_diff2 * np_diff2;

        if (np_sum4 == 0 || np_diff4 == 0 || pi4 == 0)
            return 0.0;  // Avoid division by zero

        double common_factor = (-4.0 / L3 * m1_term) / (2.0 * pi4);

        double term1_num = L * ((6.0 * L3) / np_diff4 + (6.0 * L3) / np_sum4);
        double term1 = common_factor * term1_num;

        double term2_cos_sum =
            3.0 * L * cos(pi * np_sum) * (2.0 * L2 - L2 * pi2 * np_sum2);
        double term2_cos_diff =
            3.0 * L * cos(pi * np_diff) * (2.0 * L2 - L2 * pi2 * np_diff2);
        double term2_num =
            L * (term2_cos_sum / np_sum4 + term2_cos_diff / np_diff4);
        double term2 = common_factor * term2_num;

        return term1 - term2;
    }
    else
    {
        return 0.0;  // Case n=p=0 or n!=p covered
    }
}

inline double compute_i12(int m1, int n, int p, double L, double pi)
{
    double m1_term = 2.0 * pow_m1(m1) + 3.0;
    double L2 = L * L;
    double L3 = L2 * L;
    double pi2 = pi * pi;
    double n_dbl = static_cast<double>(n);
    double p_dbl = static_cast<double>(p);
    double n2 = n_dbl * n_dbl;
    double p2 = p_dbl * p_dbl;

    if (n == p && n != 0)
    {
        if (pi2 * p2 == 0) return 0.0;
        return (6.0 / L2 * m1_term) * L3 / 6.0 +
               (6.0 / L2 * m1_term) * L3 / (4.0 * pi2 * p2);
    }
    else if (n != p)
    {
        double np_sum = n_dbl + p_dbl;
        double np_sum2 = np_sum * np_sum;
        double np_diff = n_dbl - p_dbl;
        double np_diff2 = np_diff * np_diff;

        if (pi2 == 0 || np_diff2 == 0 || np_sum2 == 0)
            return 0.0;  // Avoid division by zero

        double common_factor = (6.0 / L2 * m1_term) * L3 / pi2;
        double term1 = common_factor * cos(pi * np_diff) / np_diff2;
        double term2 = common_factor * cos(pi * np_sum) / np_sum2;
        return term1 + term2;
    }
    else
    {
        return 0.0;  // Case n=p=0 or n!=p covered
    }
}
inline double compute_i13(int, int n, int p, double L, double) // m1, pi unused
{
    return (n == p) ? -L / 2.0 : 0.0;
}

// Base function to implement the common loop structure for i*_mat
inline Tensor3D i_mat_base(
    int Npsi,
    int Nphi,
    double L,
    std::function<double(int, int, int, double, double)> compute_entry
)
{
    Tensor3D s(Npsi, Matrix::Zero(Nphi, Nphi));
    double pi_val = M_PI; // Use M_PI directly

    for (int m = 1; m <= Npsi; ++m)
    {
        int m1 = m - 1;
        for (int n_loop = 1; n_loop <= Nphi; ++n_loop) // Renamed n to n_loop to avoid conflict with param n in compute_i*
        {
            for (int p_loop = 1; p_loop <= Nphi; ++p_loop) // Renamed p to p_loop
            {
                s[m1](n_loop - 1, p_loop - 1) = compute_entry(m1, n_loop, p_loop, L, pi_val);
            }
        }
    }
    return s;
}


//----------------------------------------------------------------------------
// Partial Integral Tensors (i*_mat) - Main Functions
//----------------------------------------------------------------------------

inline Tensor3D i1_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i1);
}

inline Tensor3D i2_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i2);
}

inline Tensor3D i3_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i3);
}

inline Tensor3D i4_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i4);
}

inline Tensor3D i5_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i5);
}

inline Tensor3D i9_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i9);
}

inline Tensor3D i10_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i10);
}

inline Tensor3D i11_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i11);
}

inline Tensor3D i12_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i12);
}

inline Tensor3D i13_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i13);
}

//----------------------------------------------------------------------------
// Partial Integral Cache
//----------------------------------------------------------------------------

using PartialIntegralsCache = std::map<std::string, Tensor3D>;

inline PartialIntegralsCache compute_partial_integrals(
    int Npsi,
    int Nphi,
    double Lx,
    double Ly
)
{
    PartialIntegralsCache partials;
    partials["i1_Lx"] = i1_mat(Npsi, Nphi, Lx);
    partials["i2_Lx"] = i2_mat(Npsi, Nphi, Lx);
    partials["i3_Lx"] = i3_mat(Npsi, Nphi, Lx);
    partials["i4_Lx"] = i4_mat(Npsi, Nphi, Lx);
    partials["i5_Lx"] = i5_mat(Npsi, Nphi, Lx);

    partials["i9_Lx"] = i9_mat(Npsi, Nphi, Lx);
    partials["i10_Lx"] = i10_mat(Npsi, Nphi, Lx);
    partials["i11_Lx"] = i11_mat(Npsi, Nphi, Lx);
    partials["i12_Lx"] = i12_mat(Npsi, Nphi, Lx);
    partials["i13_Lx"] = i13_mat(Npsi, Nphi, Lx);

    partials["i1_Ly"] = i1_mat(Npsi, Nphi, Ly);
    partials["i2_Ly"] = i2_mat(Npsi, Nphi, Ly);
    partials["i3_Ly"] = i3_mat(Npsi, Nphi, Ly);
    partials["i4_Ly"] = i4_mat(Npsi, Nphi, Ly);
    partials["i5_Ly"] = i5_mat(Npsi, Nphi, Ly);

    partials["i9_Ly"] = i9_mat(Npsi, Nphi, Ly);
    partials["i10_Ly"] = i10_mat(Npsi, Nphi, Ly);
    partials["i11_Ly"] = i11_mat(Npsi, Nphi, Ly);
    partials["i12_Ly"] = i12_mat(Npsi, Nphi, Ly);
    partials["i13_Ly"] = i13_mat(Npsi, Nphi, Ly);

    return partials;
}

//----------------------------------------------------------------------------
// S Matrix Builder
//----------------------------------------------------------------------------

enum class FactorMode { NONE, N, P, NP };

inline Tensor3D build_s_matrix(
    int Npsi,
    int Nphi,
    const std::vector<std::reference_wrapper<const Tensor3D>>& partials_refs,
    const Eigen::VectorXi& idx_array,
    FactorMode factor_mode
)
{
    Tensor3D s(Npsi, Matrix::Zero(Nphi, Nphi));

    for (int m = 0; m < Npsi; ++m)
    {
        for (int n_outer = 0; n_outer < Nphi; ++n_outer) // Renamed n to n_outer
        {
            int n_idx = idx_array(n_outer);
            if (n_idx < 1) continue;
            double n_idx_dbl = static_cast<double>(n_idx);

            for (int p_outer = 0; p_outer < Nphi; ++p_outer) // Renamed p to p_outer
            {
                int p_idx = idx_array(p_outer);
                if (p_idx < 1) continue;
                double p_idx_dbl = static_cast<double>(p_idx);

                double val = 0.0;
                for (const auto& tensor_ref : partials_refs)
                {
                    const Tensor3D& tensor = tensor_ref.get();
                    if (m < tensor.size() && (n_idx - 1) < tensor[m].rows() &&
                        (p_idx - 1) < tensor[m].cols())
                    {
                        val += tensor[m](n_idx - 1, p_idx - 1);
                    }
                    else
                    {
                        std::cerr << "Warning: Potential out-of-bounds "
                                     "access in build_s_matrix for m="
                                  << m << ", n_idx=" << n_idx
                                  << ", p_idx=" << p_idx << std::endl;
                    }
                }

                switch (factor_mode)
                {
                    case FactorMode::N: val *= (n_idx_dbl * n_idx_dbl); break;
                    case FactorMode::P: val *= (p_idx_dbl * p_idx_dbl); break;
                    case FactorMode::NP:
                        val *= (n_idx_dbl * p_idx_dbl);
                        break;
                    case FactorMode::NONE:
                    default:
                        break;
                }
                s[m](n_outer, p_outer) = val;
            }
        }
    }
    return s;
}


//----------------------------------------------------------------------------
// G Matrix Builders Helper
//----------------------------------------------------------------------------

inline Matrix build_g_matrix(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& k_indices,
    const PartialIntegralsCache& cache,
    const std::vector<std::string>& i_keys_x, // Changed from i_keys_x
    const std::vector<std::string>& i_keys_y, // Changed from i_keys_y, used if !use_primary_keys
    FactorMode factor_mode,
    bool use_kx // true for g1, g2, g5 (use kx/i_keys_x); false for g3, g4, g6 (use ky/i_keys_y)
)
{
    std::vector<std::reference_wrapper<const Tensor3D>> partials_refs;
    const auto& keys = use_kx ? i_keys_x : i_keys_y;
    for (const auto& key : keys)
    {
        auto it = cache.find(key);
        if (it != cache.end())
        {
            partials_refs.push_back(std::cref(it->second));
        }
        else { throw std::runtime_error("Cache miss for key: " + key); }
    }

    Tensor3D s =
        build_s_matrix(Npsi, Nphi, partials_refs, k_indices, factor_mode);

    // Reshape s (list of Npsi matrices of Nphi x Nphi) into a single matrix
    // Python: s_reshaped = s.reshape((Npsi, Nphi**2), order="F")
    // Python: m_mat = np.repeat(s_reshaped, Npsi, axis=0) OR
    // np.tile(s_reshaped, (Npsi, 1))

    int Nphi2 = Nphi * Nphi;
    Matrix s_reshaped(Npsi, Nphi2);
    for (int m = 0; m < Npsi; ++m)
    {
        // Eigen defaults to ColMajor (Fortran order), so map directly
        s_reshaped.row(m) = Eigen::Map<const Matrix>(s[m].data(), Nphi, Nphi)
                                .transpose()
                                .reshaped(1, Nphi2);
        // Check order: Python's order='F' flattens columns first.
        // Eigen ColMajor stores column-wise. Eigen RowMajor stores row-wise.
        // Eigen's reshape defaults to ColMajor output order.
        // Map<Matrix> assumes ColMajor storage of s[m].data().
        // s[m] is Nphi x Nphi. Reshaping to 1 x Nphi2 needs correct order.
        // Let's flatten explicitly matching Fortran order:
        Matrix temp_flat(1, Nphi2);
        for (int j = 0; j < Nphi; ++j)
        {  // col index
            for (int i = 0; i < Nphi; ++i)
            {  // row index
                temp_flat(0, j * Nphi + i) = s[m](i, j);
            }
        }
        s_reshaped.row(m) = temp_flat;
    }

    Matrix m_mat;
    if (use_kx)
    {  // Corresponds to g1, g2, g5 (repeat rows)
        // Python: np.repeat(s_reshaped, Npsi, axis=0)
        // Eigen equivalent: block-wise repeat
        m_mat.resize(Npsi * Npsi, Nphi2);
        for (int i = 0; i < Npsi; ++i)
        {
            m_mat.block(i * Npsi, 0, Npsi, Nphi2) =
                s_reshaped.row(i).replicate(Npsi, 1);
        }
    }
    else
    {  // Corresponds to g3, g4, g6 (tile matrix)
        // Python: np.tile(s_reshaped, (Npsi, 1))
        m_mat = s_reshaped.replicate(Npsi, 1);
    }

    if (S > m_mat.rows())
    {
        std::cerr << "Warning: S (" << S << ") is greater than m_mat rows ("
                  << m_mat.rows() << ") in g function." << std::endl;
        S = m_mat.rows();  // Adjust S to available rows
    }
    return m_mat.topRows(S);
}


//----------------------------------------------------------------------------
// G Matrix Builders (Specific Functions)
//----------------------------------------------------------------------------

inline Matrix g1(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& kx,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi, Nphi, S, kx, cache,
        {"i1_Lx", "i2_Lx", "i3_Lx", "i4_Lx", "i5_Lx"}, // primary (Lx)
        {}, // secondary (unused)
        FactorMode::N,
        true // use primary keys (Lx related), repeat_rows_style for expansion
    );
}

inline Matrix g2(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& kx,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi, Nphi, S, kx, cache,
        {"i1_Lx", "i2_Lx", "i3_Lx", "i4_Lx", "i5_Lx"}, // primary (Lx)
        {}, // secondary (unused)
        FactorMode::P,
        true // use primary keys (Lx related), repeat_rows_style for expansion
    );
}

inline Matrix g3(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& ky,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi, Nphi, S, ky, cache,
        {}, // primary (unused)
        {"i1_Ly", "i2_Ly", "i3_Ly", "i4_Ly", "i5_Ly"},  // secondary (Ly)
        FactorMode::N,
        false // use secondary keys (Ly related), tile_matrix_style for expansion
    );
}

inline Matrix g4(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& ky,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi, Nphi, S, ky, cache,
        {}, // primary (unused)
        {"i1_Ly", "i2_Ly", "i3_Ly", "i4_Ly", "i5_Ly"},  // secondary (Ly)
        FactorMode::P,
        false // use secondary keys (Ly related), tile_matrix_style for expansion
    );
}

inline Matrix g5(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& kx,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi, Nphi, S, kx, cache,
        {"i9_Lx", "i10_Lx", "i11_Lx", "i12_Lx", "i13_Lx"}, // primary (Lx, 9-13)
        {}, // secondary (unused)
        FactorMode::NP,
        true // use primary keys (Lx related), repeat_rows_style for expansion
    );
}

inline Matrix g6(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& ky,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi, Nphi, S, ky, cache,
        {}, // primary (unused)
        {"i9_Ly", "i10_Ly", "i11_Ly", "i12_Ly", "i13_Ly"},  // secondary (Ly, 9-13)
        FactorMode::NP,
        false // use secondary keys (Ly related), tile_matrix_style for expansion
    );
}


//----------------------------------------------------------------------------
// H Tensor Computation
//----------------------------------------------------------------------------

struct HTensors { // Already defined, but ensure it's before usage if not.
    Tensor3D H0;
    Tensor3D H1;
    Tensor3D H2;
};

inline void zero_small_values(Matrix& matrix, double threshold)
{
    if (matrix.size() == 0) return;
    double max_abs = matrix.cwiseAbs().maxCoeff();
    if (max_abs == 0.0) return;
    double actual_threshold = threshold * max_abs;

    matrix =
        (matrix.array().abs() < actual_threshold).select(0.0, matrix.array());
}


inline void H_tensor_rectangular(
    const AiryCoefficients& coeffs,
    int Nphi,
    int Npsi, // Original Npsi for partial integrals
    double Lx,
    double Ly,
    const Eigen::VectorXi& kx,
    const Eigen::VectorXi& ky,
    Matrix& H0_mat_out, // Renamed to avoid conflict if H0_mat was a local var
    Matrix& H1_mat_out,
    Matrix& H2_mat_out
)
{

    int S = coeffs.S;  // Number of valid modes used for coefficients
    int Nphi2 = Nphi * Nphi;
    int dim = coeffs.coeff0.rows();  // Should be Npsi * Npsi

    PartialIntegralsCache partials =
        compute_partial_integrals(Npsi, Nphi, Lx, Ly); // Npsi here is for basis functions of integrals

    // G-matrices must be (dim_K_M x Nphi^2)
    // The S parameter passed to g functions should be dim_K_M.
    Matrix m1 = g1(Npsi, Nphi, S, kx, partials);
    Matrix m2 = g2(Npsi, Nphi, S, kx, partials);
    Matrix m3 = g3(Npsi, Nphi, S, ky, partials);
    Matrix m4 = g4(Npsi, Nphi, S, ky, partials);
    Matrix m5 = g5(Npsi, Nphi, S, kx, partials);
    Matrix m6 = g6(Npsi, Nphi, S, ky, partials);

    Matrix G_combo = (m1.array() * m4.array() + m2.array() * m3.array() -
                      2.0 * m5.array() * m6.array())
                         .matrix();

    // Compute H tensor rows
    for (int n = 0; n < S; ++n)
    {
        H0_mat_out.row(n) = coeffs.coeff0.col(n).transpose() * G_combo;
        H1_mat_out.row(n) = coeffs.coeff1.col(n).transpose() * G_combo;
        H2_mat_out.row(n) = coeffs.coeff2.col(n).transpose() * G_combo;
    }

    double const_factor =
        4.0 * M_PI * M_PI * M_PI * M_PI / (Lx * Lx * Lx * Ly * Ly * Ly);
    H0_mat_out *= const_factor;
    H1_mat_out *= const_factor;
    H2_mat_out *= const_factor;

    zero_small_values(H0_mat_out, 1e-8);
    zero_small_values(H1_mat_out, 1e-10);
    zero_small_values(H2_mat_out, 1e-8);
}


//----------------------------------------------------------------------------
// Main Coupling Matrix Computation
//----------------------------------------------------------------------------

inline void compute_coupling_matrix(
    int n_psi, // Npsi for K,M matrices
    int n_phi, // Nphi for integrals and H tensor output columns
    double lx,
    double ly,
    const Eigen::VectorXi& kx_indices,
    const Eigen::VectorXi& ky_indices,
    Matrix& H1_out // Output coupling matrix
)
{
    int n_phi2 = n_phi * n_phi;
    auto [K, M] = assemble_K_and_M(n_psi, lx, ly);

    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> ges;
    ges.compute(K, M);

    Vector vals = ges.eigenvalues();
    Matrix vecs = ges.eigenvectors();

    AiryCoefficients coeffs =
        airy_stress_coefficients(n_psi, lx, ly, vals, vecs);

    // H0 is divided by the norm
    // H1 is divided by the (norm * sqrt(varpsi))
    // H2 is divided by the (norm * varpsi)
    // We only need H1 because it is the only one used in the coupling matrix
    // because it is normal
    Matrix H0_mat = Matrix::Zero(coeffs.S, n_phi*n_phi);
    Matrix H1_mat = Matrix::Zero(coeffs.S, n_phi*n_phi);
    Matrix H2_mat = Matrix::Zero(coeffs.S, n_phi*n_phi);

    // 4. Compute H Tensors
    H_tensor_rectangular(
        coeffs,
        n_phi,
        n_psi,
        lx,
        ly,
        kx_indices,
        ky_indices,
        H0_mat,
        H1_mat,
        H2_mat
    );
    H1_out = H1_mat.block(0, 0, n_psi, n_phi*n_phi);

}


// Add PI constant if not defined by <cmath> or Eigen
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 