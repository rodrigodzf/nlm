#include "coupling.h"
#include <cmath>
#include <vector>
#include <complex>
#include <numeric>
#include <algorithm>
#include <iostream>  // For potential debugging
#include "matioCpp/EigenConversions.h"
#include <Eigen/src/Core/Matrix.h>
#include <matioCpp/matioCpp.h>
#include <vector>
#include <numeric>    // for std::iota
#include <algorithm>  // for std::sort, std::stable_sort

// Helper function for (-1)^n
inline double pow_m1(int n)
{
    return (n % 2 == 0) ? 1.0 : -1.0;
}

//----------------------------------------------------------------------------
// Basic Integrals
//----------------------------------------------------------------------------

double int4(int m, int p, double L)
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
        y = (den == 0) ? 0.0 : num / den;  // Avoid division by zero if p=0
                                           // (though loop logic might prevent
                                           // this)
    }
    else if (p == 0)
    {
        double num = 60.0 * (m1_m + 1.0) * (pi2 * m2 - 42.0);
        double den = 7.0 * L * pi2 * m2;
        y = (den == 0) ? 0.0 : num / den;  // Avoid division by zero if m=0
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

double int1(int m, int p, double L)
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

double int2(int m, int p, double L)
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
        if (pi4 * m4 == 0) return 0.0;  // Avoid division by zero
        y = (67.0 * L) / 70.0 - (m1_m * L) / 35.0 - (768.0 * L) / (pi4 * m4) -
            (672.0 * m1_m * L) / (pi4 * m4);
    }
    else if (m == 0)
    {
        if (pi4 * p4 == 0) return 0.0;  // Avoid division by zero
        double num = 3.0 * L * (m1_p + 1.0) * (pi4 * p4 - 1680.0);
        double den = 14.0 * pi4 * p4;
        y = num / den;
    }
    else if (p == 0)
    {
        if (pi4 * m4 == 0) return 0.0;  // Avoid division by zero
        double num = 3.0 * L * (m1_m + 1.0) * (pi4 * m4 - 1680.0);
        double den = 14.0 * pi4 * m4;
        y = num / den;
    }
    else
    {
        if (pi4 * m4 == 0 || p4 == 0) return 0.0;  // Avoid division by zero
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

Matrix build_I1(int N, double L)
{
    Matrix I = Matrix::Zero(N, N);
    for (int m = 0; m < N; ++m)
    {
        for (int p = 0; p < N; ++p) { I(m, p) = int1(m, p, L); }
    }
    return I;
}

Matrix build_I2(int N, double L)
{
    Matrix I = Matrix::Zero(N, N);
    for (int m = 0; m < N; ++m)
    {
        for (int p = 0; p < N; ++p) { I(m, p) = int2(m, p, L); }
    }
    return I;
}

Matrix build_I4(int N, double L)
{
    Matrix I = Matrix::Zero(N, N);
    for (int m = 0; m < N; ++m)
    {
        for (int p = 0; p < N; ++p) { I(m, p) = int4(m, p, L); }
    }
    return I;
}

Matrix int2_mat(int N, double L)
{
    // This is identical to build_I2, kept for potential specific usage if
    // needed In Python, airy_stress_coefficients specifically called int2_mat
    // with different L values
    return build_I2(N, L);
}

//----------------------------------------------------------------------------
// Kronecker Product Helper
//----------------------------------------------------------------------------

Matrix kron(const Matrix& A, const Matrix& B)
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

std::pair<Matrix, Matrix> assemble_K_and_M(int Npsi, double Lx, double Ly)
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

AiryCoefficients airy_stress_coefficients(
    int n_psi,
    double Lx,  // Need Lx, Ly for int2_mat call inside
    double Ly,
    const Vector& vals,  // these should always be real
    const Matrix& vecs
)
{
    int dim_total = vals.size();  // Should be n_psi * n_psi initially
    if (vecs.cols() != dim_total)
    {
        throw std::runtime_error("Eigenvalues and Eigenvectors size mismatch."
        );
    }

    std::vector<int> good_indices;
    const double eps_imag = 1e-12;
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

    // Sort by real part of eigenvalues
    std::vector<int> idx_sort(S_good);
    std::iota(idx_sort.begin(), idx_sort.end(), 0);  // 0, 1, 2, ... S_good-1

    // Sort indices based on corresponding eigenvalue magnitude
    std::stable_sort(
        idx_sort.begin(),
        idx_sort.end(),
        [&good_vals_real](int i1, int i2)
        { return good_vals_real(i1) < good_vals_real(i2); }
    );

    Vector auto_sorted(S_good);
    Matrix coeff_sorted(vecs.rows(), S_good);  // vecs.rows() == dim_total
    for (int i = 0; i < S_good; ++i)
    {
        auto_sorted(i) = good_vals_real(idx_sort[i]);
        coeff_sorted.col(i) = good_vecs.col(idx_sort[i]);
    }

    // Build the factor arrays (coeff0, coeff1, coeff2)
    int dim = n_psi * n_psi;  // This is the size related to K, M matrices
                              // (vecs.rows())
    if (dim != vecs.rows())
    {
        // This check might be redundant if dim_total was checked correctly
        throw std::runtime_error(
            "Dimension mismatch in airy_stress_coefficients."
        );
    }

    Matrix coeff0 = Matrix::Zero(dim, S_good);
    Matrix coeff1 = Matrix::Zero(dim, S_good);
    Matrix coeff2 = Matrix::Zero(dim, S_good);

    // --- nmatr calculation ---
    const int N2 = n_psi * n_psi;
    const int N4 = N2 * N2;
    Matrix I2x_norm = int2_mat(n_psi, Lx);
    Matrix I2y_norm = int2_mat(n_psi, Ly);
#if 0
    Matrix M_norm = kron(I2x_norm, I2y_norm);  // Use the actual mass matrix
    for (int d = 0; d < S_good; ++d)
    {
        Vector temp = coeff_sorted.col(d);

        // Calculate norm: temp^T * M_norm * temp
        double norm_sq = temp.transpose() * M_norm * temp;
        double norm_val =
            (norm_sq > 0) ? std::sqrt(norm_sq) : 1.0;  // Avoid sqrt(0) or
                                                       // negative
        if (norm_val == 0) norm_val = 1.0;  // Prevent division by zero

        double auto_d = auto_sorted(d);
        double sqrt_auto_d =
            (auto_d > 0) ? std::sqrt(auto_d) : 1.0;  // Avoid sqrt of
                                                     // zero/negative

        coeff0.col(d) = temp / norm_val;
        coeff1.col(d) = temp / (norm_val * sqrt_auto_d);
        coeff2.col(d) = (auto_d != 0) ? (temp / (norm_val * auto_d)).eval()
                                      : Vector::Zero(dim);
    }
#else
 
    Vector NN = Eigen::Map<const Vector>(I2x_norm.data(), N2);
    Vector MM = Eigen::Map<const Vector>(I2y_norm.data(), N2);

    Matrix NNrep = NN.replicate(1, N2);         // size N2×N2
    Matrix MMrep = MM.transpose().replicate(N2, 1); // size N2×N2

    // 4) element-wise multiply
    Matrix nmatr = NNrep.array() * MMrep.array(); // size N2×N2

    matioCpp::File nmatr_results = matioCpp::File::Create("cpp_airy.mat");
    if (!nmatr_results.isOpen()) {
        std::cout << "set_couplings_and_eigenvalues: could not open file 'cpp_airy.mat'" << std::endl;
    }
    auto mv = matioCpp::make_variable("nmatr_cpp", nmatr);
    nmatr_results.write(mv);

    Vector Mvec(N4);
    for (int i1 = 0; i1 < n_psi; ++i1) {
        for (int i2 = 0; i2 < n_psi; ++i2) {
            for (int i3 = 0; i3 < n_psi; ++i3) {
                for (int i4 = 0; i4 < n_psi; ++i4) {
                    // Original matrix indices
                    int i_mat = i1 + i2 * n_psi;  // row index in M_mat
                    int j_mat = i3 + i4 * n_psi;  // col index in M_mat
                    
                    // After permute [4 1 3 2], the dimensions are:
                    // new1 = i4, new2 = i1, new3 = i3, new4 = i2
                    // For column-major flattening, the index is:
                    int new_index = i4 + i1 * n_psi + i3 * (n_psi * n_psi) + i2 * (n_psi * n_psi * n_psi);
                    
                    Mvec(new_index) = nmatr(i_mat, j_mat);
                }
            }
        }
    }

    auto mv2 = matioCpp::make_variable("mvec_cpp", Mvec);
    nmatr_results.write(mv2);

    Vector temp = coeff_sorted.col(0);
    Matrix temp_replicated = temp.replicate(1, N2);
    Matrix temp2 = temp_replicated.transpose();
    Matrix temp3_mat = temp_replicated.cwiseProduct(temp2);
    Vector temp3 = Eigen::Map<Vector>(temp3_mat.data(), N4);
    auto temp_var = matioCpp::make_variable("temp_cpp", temp);
    nmatr_results.write(temp_var);

    double norm_sq = Mvec.dot(temp3);
    std::cout << "norm_sq: " << norm_sq << std::endl;

    auto norm_sq_var = matioCpp::make_variable("norm_sq_cpp", norm_sq);
    nmatr_results.write(norm_sq_var);


    for (int d = 0; d < S_good; ++d) {

        Vector temp = coeff_sorted.col(d);
        Matrix temp_replicated = temp.replicate(1, N2);
        Matrix temp2 = temp_replicated.transpose();
        Matrix temp3_mat = temp_replicated.cwiseProduct(temp2);
        Vector temp3 = Eigen::Map<Vector>(temp3_mat.data(), N4);
        double norm_sq = Mvec.dot(temp3);
        double norm_val = std::sqrt(norm_sq);

        double auto_d = auto_sorted[d];
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
// Partial Integral Tensors (i*_mat)
//----------------------------------------------------------------------------

// Base function to implement the common loop structure for i*_mat
Tensor3D i_mat_base(
    int Npsi,
    int Nphi,
    double L,
    std::function<double(int, int, int, double, double)> compute_entry
)
{
    Tensor3D s(Npsi, Matrix::Zero(Nphi, Nphi));
    double pi = M_PI;

    // Loop using 1-based indexing consistent with Python/MATLAB source
    for (int m = 1; m <= Npsi; ++m)
    {
        int m1 = m - 1;  // 0-based index for array access
        for (int n = 1; n <= Nphi; ++n)
        {
            for (int p = 1; p <= Nphi; ++p)
            {
                s[m - 1](n - 1, p - 1) = compute_entry(m1, n, p, L, pi);
            }
        }
    }
    return s;
}

// Implement compute_entry for each i*_mat function

double compute_i1(int m1, int n, int p, double L, double pi)
{
    if (m1 == 0 && n == p) return L / 2.0;
    if (m1 == (p - n) || m1 == (n - p)) return L / 4.0;
    // MATLAB version had (-n-p) case returning -L/4.0, Python version had
    // L/4.0? Checking python... Python i1_mat: `elif m1 == (-n-p) or m1 ==
    // (n+p): s[...] = -L/4.0` -> Corrected below
    if (m1 == (-n - p) || m1 == (n + p)) return -L / 4.0;
    return 0.0;
}
Tensor3D i1_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i1);
}

double compute_i2(int m1, int n, int p, double L, double pi)
{
    // NOTE: This function in Python had hardcoded large constants likely from
    // symbolic math. These might be prone to floating point issues. Direct
    // translation follows. Consider deriving these symbolically or testing
    // thoroughly.
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
        // Extremely complex terms from Python, likely from symbolic
        // integration These constants are suspicious (e.g., 1713638851887625
        // / 17592186044416 approx 97.4) Re-derivation might be needed for
        // robustness. Translating directly for now.
        double np_sum = n_dbl + p_dbl;
        double np_diff = n_dbl - p_dbl;
        double np_sum2 = np_sum * np_sum;
        double np_sum4 = np_sum2 * np_sum2;
        double np_sum5 = np_sum4 * np_sum;
        double np_diff2 = np_diff * np_diff;
        double np_diff4 = np_diff2 * np_diff2;
        double np_diff5 = np_diff4 * np_diff;

        if (np_sum5 == 0 || np_diff5 == 0)
            return 0.0;  // Avoid division by zero

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
Tensor3D i2_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i2);
}

double compute_i3(int m1, int n, int p, double L, double pi)
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
Tensor3D i3_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i3);
}

double compute_i4(int m1, int n, int p, double L, double pi)
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
Tensor3D i4_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i4);
}

double compute_i5(int m1, int n, int p, double L, double pi)
{
    return (n == p) ? -L / 2.0 : 0.0;
}
Tensor3D i5_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i5);
}

double compute_i9(int m1, int n, int p, double L, double pi)
{
    if (m1 == 0 && n == p) return L / 2.0;
    if (m1 == (p - n) || m1 == (n - p)) return L / 4.0;
    // MATLAB/Python source: elif m1 == (-n-p) or m1 == (n+p): s[...] = L/4.0
    // (Positive sign)
    if (m1 == (-n - p) || m1 == (n + p)) return L / 4.0;
    return 0.0;
}
Tensor3D i9_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i9);
}

double compute_i10(int m1, int n, int p, double L, double pi)
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
Tensor3D i10_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i10);
}

double compute_i11(int m1, int n, int p, double L, double pi)
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
Tensor3D i11_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i11);
}

double compute_i12(int m1, int n, int p, double L, double pi)
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
Tensor3D i12_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i12);
}

double compute_i13(int m1, int n, int p, double L, double pi)
{
    return (n == p) ? -L / 2.0 : 0.0;
}
Tensor3D i13_mat(int Npsi, int Nphi, double L)
{
    return i_mat_base(Npsi, Nphi, L, compute_i13);
}

//----------------------------------------------------------------------------
// Partial Integral Cache
//----------------------------------------------------------------------------

PartialIntegralsCache
    compute_partial_integrals(int Npsi, int Nphi, double Lx, double Ly)
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

Tensor3D build_s_matrix(
    int Npsi,
    int Nphi,
    const std::vector<std::reference_wrapper<const Tensor3D>>& partials_refs,
    const Eigen::VectorXi& idx_array,  // Assuming 1-based indices
    FactorMode factor_mode
)
{
    Tensor3D s(Npsi, Matrix::Zero(Nphi, Nphi));

    for (int m = 0; m < Npsi; ++m)
    {
        for (int n = 0; n < Nphi; ++n)
        {
            int n_idx = idx_array(n);  // 1-based index from input array
            if (n_idx < 1) continue;   // Skip if index is invalid
            double n_idx_dbl = static_cast<double>(n_idx);

            for (int p = 0; p < Nphi; ++p)
            {
                int p_idx = idx_array(p);  // 1-based index from input array
                if (p_idx < 1) continue;   // Skip if index is invalid
                double p_idx_dbl = static_cast<double>(p_idx);

                double val = 0.0;
                for (const auto& tensor_ref : partials_refs)
                {
                    const Tensor3D& tensor = tensor_ref.get();
                    // Access using 0-based indices derived from 1-based
                    // n_idx, p_idx
                    if (m < tensor.size() && (n_idx - 1) < tensor[m].rows() &&
                        (p_idx - 1) < tensor[m].cols())
                    {
                        val += tensor[m](n_idx - 1, p_idx - 1);
                    }
                    else
                    {
                        // Handle potential out-of-bounds access if Nphi in
                        // idx_array > Nphi used in i*_mat Or if Npsi > Npsi
                        // used in i*_mat
                        std::cerr << "Warning: Potential out-of-bounds "
                                     "access in build_s_matrix for m="
                                  << m << ", n_idx=" << n_idx
                                  << ", p_idx=" << p_idx << std::endl;
                    }
                }

                // Apply factor
                switch (factor_mode)
                {
                    case FactorMode::N: val *= (n_idx_dbl * n_idx_dbl); break;
                    case FactorMode::P: val *= (p_idx_dbl * p_idx_dbl); break;
                    case FactorMode::NP:
                        val *= (n_idx_dbl * p_idx_dbl);
                        break;
                    case FactorMode::NONE:
                    default:
                        // No factor
                        break;
                }
                s[m](n, p) = val;  // Store in 0-based (n,p) of the output
                                   // matrix
            }
        }
    }
    return s;
}

//----------------------------------------------------------------------------
// G Matrix Builders Helper
//----------------------------------------------------------------------------

Matrix build_g_matrix(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& k_indices,  // kx or ky
    const PartialIntegralsCache& cache,
    const std::vector<std::string>& i_keys_x,  // e.g., {"i1_Lx", "i2_Lx",
                                               // ...}
    const std::vector<std::string>& i_keys_y,  // e.g., {"i1_Ly", "i2_Ly",
                                               // ...}
    FactorMode factor_mode,
    bool use_kx  // true for g1, g2, g5; false for g3, g4, g6
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

    // Python slices result: return m_mat[:S, :]
    // S comes from airy_stress_coefficients, which is the number of valid
    // modes The dimension Npsi*Npsi corresponds to the original K, M matrices
    // before filtering. We need to return a matrix compatible with coeff
    // matrices from airy_stress (dim x S) where dim = Npsi*Npsi. The
    // g-matrices seem to be used as (dim x S)^T @ g_mat? No... Python: f0 =
    // coeff0[:, n].T @ G_combo, where coeff is (dim, S) and G_combo is (dim,
    // Nphi2) ? coeff0[:, n].T is (1, dim) m1 etc are computed here. What
    // should their shape be? Python makes them (S_total?, Nphi*Nphi), then
    // slices to (S, Nphi*Nphi). S is S_good here. Let's return m_mat directly
    // and slice later if needed. Check H_tensor computation. The python code
    // slices m_mat to S rows: `return m_mat[:S, :]` S here is the number of
    // valid modes from airy_stress. The number of rows in m_mat is Npsi*Npsi.
    // This slicing seems inconsistent unless S happens to equal Npsi*Npsi, or
    // maybe the original K,M matrices were already sized based on S? Looking
    // at H_tensor: `coeff0[:, n]` has dimension `dim = Npsi*Npsi`. The dot
    // product requires `m1*m4 + ...` to have dimension `dim x Nphi*Nphi`? No,
    // (dim, 1). Let's assume the Python slicing `[:S, :]` was intended to
    // match the number of modes `S`. But the first dimension of m_mat is
    // Npsi*Npsi (or Npsi repeated). Let's re-read g1-g6 python... `S` is
    // passed in. It's used only for slicing the result. `m_mat =
    // np.repeat(...)` has shape (Npsi*Npsi, Nphi2). Slicing `[:S,:]` takes
    // the first S rows. This assumes S <= Npsi*Npsi.
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

Matrix g1(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& kx,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi,
        Nphi,
        S,
        kx,
        cache,
        {"i1_Lx", "i2_Lx", "i3_Lx", "i4_Lx", "i5_Lx"},
        {},  // Only Lx keys
        FactorMode::N,
        true
    );  // Factor 'n', use kx
}

Matrix g2(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& kx,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi,
        Nphi,
        S,
        kx,
        cache,
        {"i1_Lx", "i2_Lx", "i3_Lx", "i4_Lx", "i5_Lx"},
        {},  // Only Lx keys
        FactorMode::P,
        true
    );  // Factor 'p', use kx
}

Matrix g3(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& ky,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi,
        Nphi,
        S,
        ky,
        cache,
        {},
        {"i1_Ly", "i2_Ly", "i3_Ly", "i4_Ly", "i5_Ly"},  // Only Ly keys
        FactorMode::N,
        false
    );  // Factor 'n', use ky
}

Matrix g4(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& ky,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi,
        Nphi,
        S,
        ky,
        cache,
        {},
        {"i1_Ly", "i2_Ly", "i3_Ly", "i4_Ly", "i5_Ly"},  // Only Ly keys
        FactorMode::P,
        false
    );  // Factor 'p', use ky
}

Matrix g5(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& kx,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi,
        Nphi,
        S,
        kx,
        cache,
        {"i9_Lx", "i10_Lx", "i11_Lx", "i12_Lx", "i13_Lx"},
        {},  // Only Lx keys (9-13)
        FactorMode::NP,
        true
    );  // Factor 'np', use kx
}

Matrix g6(
    int Npsi,
    int Nphi,
    int S,
    const Eigen::VectorXi& ky,
    const PartialIntegralsCache& cache
)
{
    return build_g_matrix(
        Npsi,
        Nphi,
        S,
        ky,
        cache,
        {},
        {"i9_Ly", "i10_Ly", "i11_Ly", "i12_Ly", "i13_Ly"},  // Only Ly keys
                                                            // (9-13)
        FactorMode::NP,
        false
    );  // Factor 'np', use ky
}

//----------------------------------------------------------------------------
// H Tensor Computation
//----------------------------------------------------------------------------

void zero_small_values(Matrix& matrix, double threshold)
{
    if (matrix.size() == 0) return;
    double max_abs = matrix.cwiseAbs().maxCoeff();
    if (max_abs == 0.0) return;
    double actual_threshold = threshold * max_abs;

    matrix =
        (matrix.array().abs() < actual_threshold).select(0.0, matrix.array());
}

void H_tensor_rectangular(
    const AiryCoefficients& coeffs,
    int Nphi,
    int Npsi,
    double Lx,
    double Ly,
    const Eigen::VectorXi& kx,
    const Eigen::VectorXi& ky,
    Matrix& H0_mat,
    Matrix& H1_mat,
    Matrix& H2_mat
)
{
    int S = coeffs.S;  // Number of valid modes used for coefficients
    int Nphi2 = Nphi * Nphi;
    int dim = coeffs.coeff0.rows();  // Should be Npsi * Npsi

    PartialIntegralsCache partials =
        compute_partial_integrals(Npsi, Nphi, Lx, Ly);

    // Compute integral components (g matrices)
    int S_coeffs = coeffs.S;
    Matrix m1 = g1(Npsi, Nphi, S_coeffs, kx, partials);  // Shape (S_coeffs,
                                                         // Nphi2)
    Matrix m2 = g2(Npsi, Nphi, S_coeffs, kx, partials);  // Shape (S_coeffs,
                                                         // Nphi2)
    Matrix m3 = g3(Npsi, Nphi, S_coeffs, ky, partials);  // Shape (S_coeffs,
                                                         // Nphi2)
    Matrix m4 = g4(Npsi, Nphi, S_coeffs, ky, partials);  // Shape (S_coeffs,
                                                         // Nphi2)
    Matrix m5 = g5(Npsi, Nphi, S_coeffs, kx, partials);  // Shape (S_coeffs,
                                                         // Nphi2)
    Matrix m6 = g6(Npsi, Nphi, S_coeffs, ky, partials);  // Shape (S_coeffs,
                                                         // Nphi2)

    // First compute the element-wise product combination
    Matrix G_combo = (m1.array() * m4.array() + m2.array() * m3.array() -
                      2.0 * m5.array() * m6.array())
                         .matrix();

    // matioCpp::File matFile = matioCpp::File::Create(
    //     "cpp_H_tensor_rectangular.mat",
    //     matioCpp::FileVersion::MAT5
    // );
    // auto G_combo_var = matioCpp::make_variable("G_combo_cpp", G_combo);
    // matFile.write(G_combo_var);

    // Initialize H matrices
    // Matrix H0_mat = Matrix::Zero(S, Nphi2);
    // Matrix H1_mat = Matrix::Zero(S, Nphi2);
    // Matrix H2_mat = Matrix::Zero(S, Nphi2);

    // Compute H tensor rows
    for (int n = 0; n < S; ++n)
    {
        H0_mat.row(n) = coeffs.coeff0.col(n).transpose() * G_combo;
        H1_mat.row(n) = coeffs.coeff1.col(n).transpose() * G_combo;
        H2_mat.row(n) = coeffs.coeff2.col(n).transpose() * G_combo;
    }

    // Normalize with constants
    double const_factor =
        4.0 * M_PI * M_PI * M_PI * M_PI / (Lx * Lx * Lx * Ly * Ly * Ly);
    H0_mat *= const_factor;
    H1_mat *= const_factor;
    H2_mat *= const_factor;

    zero_small_values(H0_mat, 1e-8);
    zero_small_values(H1_mat, 1e-10);
    zero_small_values(H2_mat, 1e-8);

    // auto H0_mat_var = matioCpp::make_variable("H0_mat_cpp", H0_mat);
    // auto H1_mat_var = matioCpp::make_variable("H1_mat_cpp", H1_mat);
    // auto H2_mat_var = matioCpp::make_variable("H2_mat_cpp", H2_mat);
    // matFile.write(H0_mat_var);
    // matFile.write(H1_mat_var);
    // matFile.write(H2_mat_var);

    return;
}

//----------------------------------------------------------------------------
// Main Coupling Matrix Computation
//----------------------------------------------------------------------------

void compute_coupling_matrix(
    int n_psi,
    int n_phi,
    double lx,
    double ly,
    const Eigen::VectorXi& kx_indices,  // Assuming 1-based indices
    const Eigen::VectorXi& ky_indices,  // Assuming 1-based indices
    Matrix& H1
)
{
    // 1. Assemble K and M
    auto [K, M] = assemble_K_and_M(n_psi, lx, ly);

    // 2. Solve Generalized Eigenvalue Problem K*v = lambda*M*v
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> ges;
    ges.compute(K, M);

    Vector vals = ges.eigenvalues();
    Matrix vecs = ges.eigenvectors();

    // 3. Compute Airy Stress Coefficients
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
    H1 = H1_mat.block(0, 0, n_psi, n_phi*n_phi);

    return;
}