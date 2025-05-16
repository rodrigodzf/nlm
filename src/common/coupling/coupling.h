#pragma once

#include <vector>
#include <complex>
#include <map>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <functional> // For std::reference_wrapper

// Define type aliases for clarity
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using MatrixCD = Eigen::MatrixXcd;
using VectorCD = Eigen::VectorXcd;
using Tensor3D = std::vector<Matrix>; // Representing Npsi x Nphi x Nphi

//----------------------------------------------------------------------------
// Basic Integrals
//----------------------------------------------------------------------------

double int1(int m, int p, double L);
double int2(int m, int p, double L);
double int4(int m, int p, double L);

//----------------------------------------------------------------------------
// Integral Matrix Builders
//----------------------------------------------------------------------------

Matrix build_I1(int N, double L);
Matrix build_I2(int N, double L);
Matrix build_I4(int N, double L);
Matrix int2_mat(int N, double L); // Note: Python version took Lx=0.2, Ly=0.3, but C++ takes generic L

//----------------------------------------------------------------------------
// Kronecker Product Helper
//----------------------------------------------------------------------------

Matrix kron(const Matrix& A, const Matrix& B);

//----------------------------------------------------------------------------
// K and M Matrix Assembly
//----------------------------------------------------------------------------

std::pair<Matrix, Matrix> assemble_K_and_M(int Npsi, double Lx, double Ly);

//----------------------------------------------------------------------------
// Airy Stress Coefficients
//----------------------------------------------------------------------------

struct AiryCoefficients {
    Matrix coeff0;
    Matrix coeff1;
    Matrix coeff2;
    int S; // Number of valid modes used
    Vector auto_vec; // Eigenvalues corresponding to coeffs
};

AiryCoefficients airy_stress_coefficients(
    int n_psi,
    double Lx, // Need Lx, Ly for int2_mat call inside
    double Ly,
    const VectorCD& vals,
    const MatrixCD& vecs
);


//----------------------------------------------------------------------------
// Partial Integral Tensors (i*_mat)
//----------------------------------------------------------------------------

Tensor3D i1_mat(int Npsi, int Nphi, double L);
Tensor3D i2_mat(int Npsi, int Nphi, double L);
Tensor3D i3_mat(int Npsi, int Nphi, double L);
Tensor3D i4_mat(int Npsi, int Nphi, double L);
Tensor3D i5_mat(int Npsi, int Nphi, double L);
Tensor3D i9_mat(int Npsi, int Nphi, double L);
Tensor3D i10_mat(int Npsi, int Nphi, double L);
Tensor3D i11_mat(int Npsi, int Nphi, double L);
Tensor3D i12_mat(int Npsi, int Nphi, double L);
Tensor3D i13_mat(int Npsi, int Nphi, double L);

//----------------------------------------------------------------------------
// Partial Integral Cache
//----------------------------------------------------------------------------

using PartialIntegralsCache = std::map<std::string, Tensor3D>;

PartialIntegralsCache compute_partial_integrals(
    int Npsi,
    int Nphi,
    double Lx,
    double Ly
);

//----------------------------------------------------------------------------
// S Matrix Builder
//----------------------------------------------------------------------------

// Enum to specify factor mode
enum class FactorMode { NONE, N, P, NP };

Tensor3D build_s_matrix(
    int Npsi,
    int Nphi,
    const std::vector<std::reference_wrapper<const Tensor3D>>& partials_refs, // Pass refs to avoid copy
    const Eigen::VectorXi& idx_array, // Assuming 1-based indices
    FactorMode factor_mode
);


//----------------------------------------------------------------------------
// G Matrix Builders
//----------------------------------------------------------------------------

Matrix g1(int Npsi, int Nphi, int S, const Eigen::VectorXi& kx, const PartialIntegralsCache& cache);
Matrix g2(int Npsi, int Nphi, int S, const Eigen::VectorXi& kx, const PartialIntegralsCache& cache);
Matrix g3(int Npsi, int Nphi, int S, const Eigen::VectorXi& ky, const PartialIntegralsCache& cache);
Matrix g4(int Npsi, int Nphi, int S, const Eigen::VectorXi& ky, const PartialIntegralsCache& cache);
Matrix g5(int Npsi, int Nphi, int S, const Eigen::VectorXi& kx, const PartialIntegralsCache& cache);
Matrix g6(int Npsi, int Nphi, int S, const Eigen::VectorXi& ky, const PartialIntegralsCache& cache);


//----------------------------------------------------------------------------
// H Tensor Computation
//----------------------------------------------------------------------------

struct HTensors {
    Tensor3D H0; // Storing as vector<Matrix> for consistency (size n_psi)
    Tensor3D H1;
    Tensor3D H2;
};

HTensors H_tensor_rectangular(
    const AiryCoefficients& coeffs,
    int Nphi,
    int Npsi, // Note: Npsi here is the original one, coeffs.S might be different
    double Lx,
    double Ly,
    const Eigen::VectorXi& kx, // Assuming 1-based indices
    const Eigen::VectorXi& ky  // Assuming 1-based indices
);


//----------------------------------------------------------------------------
// Main Coupling Matrix Computation
//----------------------------------------------------------------------------

HTensors compute_coupling_matrix(
    int n_psi,
    int n_phi,
    double lx,
    double ly,
    const Eigen::VectorXi& kx_indices, // Assuming 1-based indices
    const Eigen::VectorXi& ky_indices  // Assuming 1-based indices
);

// Add PI constant if not defined by <cmath> or Eigen
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif 