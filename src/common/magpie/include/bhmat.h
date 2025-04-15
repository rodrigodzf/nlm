#pragma once

#include <Eigen/Sparse>
#include <vector>
#include <cassert>
#include <cmath>
#include <tuple> // Required for structured bindings
#include <stdexcept> // For runtime_error
#include "Dxx_coeffs.h" // Assuming this contains the templated Dxx_coeffs functions
#include <matioCpp/matioCpp.h>

// Helper function to create a sparse block matrix from its diagonal vectors
template <typename T>
Eigen::SparseMatrix<T> create_block_from_diagonals(
    int Ny,
    const Eigen::Vector<T, Eigen::Dynamic>& D0,  // Main diagonal (size Ny)
    const Eigen::Vector<T, Eigen::Dynamic>& D1,  // Super-diagonal +1 (size Ny-1)
    const Eigen::Vector<T, Eigen::Dynamic>& Dm1, // Sub-diagonal -1 (size Ny-1)
    const Eigen::Vector<T, Eigen::Dynamic>& D2,  // Super-diagonal +2 (size Ny-2)
    const Eigen::Vector<T, Eigen::Dynamic>& Dm2) // Sub-diagonal -2 (size Ny-2)
{
    Eigen::SparseMatrix<T> block(Ny, Ny);
    std::vector<Eigen::Triplet<T>> triplets;
    // Estimate capacity: Ny for D0, 2*(Ny-1) for D1/Dm1, 2*(Ny-2) for D2/Dm2
    triplets.reserve(5 * Ny); 

    // Main diagonal (0)
    if (D0.size() == Ny) {
        for (int i = 0; i < Ny; ++i) {
            if (D0(i) != static_cast<T>(0.0)) {
                triplets.emplace_back(i, i, D0(i));
            }
        }
    }

    // Super-diagonal (+1)
    if (D1.size() == Ny - 1) {
        for (int i = 0; i < Ny - 1; ++i) {
            if (D1(i) != static_cast<T>(0.0)) {
                triplets.emplace_back(i, i + 1, D1(i));
            }
        }
    }

    // Sub-diagonal (-1)
    if (Dm1.size() == Ny - 1) {
        for (int i = 0; i < Ny - 1; ++i) {
            if (Dm1(i) != static_cast<T>(0.0)) {
                triplets.emplace_back(i + 1, i, Dm1(i));
            }
        }
    }

    // Super-diagonal (+2)
    if (D2.size() == Ny - 2) {
        for (int i = 0; i < Ny - 2; ++i) {
            if (D2(i) != static_cast<T>(0.0)) {
                triplets.emplace_back(i, i + 2, D2(i));
            }
        }
    }

    // Sub-diagonal (-2)
    if (Dm2.size() == Ny - 2) {
        for (int i = 0; i < Ny - 2; ++i) {
            if (Dm2(i) != static_cast<T>(0.0)) {
                triplets.emplace_back(i + 2, i, Dm2(i));
            }
        }
    }

    block.setFromTriplets(triplets.begin(), triplets.end());
    return block;
}


// Helper function to add a block matrix's triplets to a global list with offset
template <typename T>
void add_block_to_triplets(
    int block_row, int block_col,
    const Eigen::SparseMatrix<T>& block,
    std::vector<Eigen::Triplet<T>>& global_triplets,
    int Ny, int MM)
{
    int row_offset = block_row * Ny;
    int col_offset = block_col * Ny;
    for (int k = 0; k < block.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(block, k); it; ++it) {
            int global_row = it.row() + row_offset;
            int global_col = it.col() + col_offset;
            // Ensure triplets are within the bounds of the final matrix
            if (global_row >= 0 && global_row < MM && global_col >= 0 && global_col < MM) {
                 global_triplets.emplace_back(global_row, global_col, it.value());
            }
        }
    }
}


/**
 * Generate a biharmonic matrix for a plate using block construction, mirroring Python logic.
 *
 * @tparam T The numerical type (e.g., float, double)
 * @param BCs Boundary conditions as a 4x2 matrix of type T.
 * @param Nxy Number of grid points [Nx, Ny] (integers)
 * @param h Grid spacing
 * @param Lz Plate thickness
 * @param E Young's modulus [Pa]
 * @param nu Poisson's ratio
 * @param format (unused in this implementation, kept for API compatibility)
 * @return Biharmonic matrix of size (Nx*Ny) x (Nx*Ny) of type T
 */
template <typename T>
Eigen::SparseMatrix<T> bhmat(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& BCs,
    const Eigen::Vector2i& Nxy,
    T h,
    T Lz,
    T E,
    T nu,
    const std::string& format = "compressed") // format is unused
{
    // --- Input Validation ---
    assert(BCs.rows() == 4 && BCs.cols() == 2 && "BCs must be a 4x2 matrix");
    assert(Nxy.size() == 2 && "Nxy must have 2 elements");
    assert(Nxy(0) > 0 && Nxy(1) > 0 && "Grid dimensions must be positive");
    assert(h != static_cast<T>(0.0) && "Grid spacing cannot be zero");
    assert(Lz != static_cast<T>(0.0) && "Plate thickness cannot be zero");
    assert(E != static_cast<T>(0.0) && "Young's modulus cannot be zero");

    // --- Unpack Variables ---
    T K0y = BCs(0, 0), R0y = BCs(0, 1), Kx0 = BCs(1, 0), Rx0 = BCs(1, 1);
    T KLy = BCs(2, 0), RLy = BCs(2, 1), KxL = BCs(3, 0), RxL = BCs(3, 1);
    int Nx = Nxy(0), Ny = Nxy(1);
    int MM = Nx * Ny; // Total size

    T D = E * std::pow(Lz, static_cast<T>(3)) / static_cast<T>(12.0) / (static_cast<T>(1.0) - std::pow(nu, static_cast<T>(2)));

    // --- Coefficient Calculation (Grouped by Block Relevance) ---
    // Note: Assuming Dxx_coeffs<T>(...) functions are correctly implemented and match Python.

    // Coeffs for Blk1X (ix=0 boundary)
    auto [D00u00, D00u10, D00u20, D00u01, D00u02, D00u11] = D00_coeffs<T>(K0y, R0y, Kx0, Rx0, h, D, nu);
    auto [D01u01, D01u11, D01u21, D01u00, D01u02, D01u03, D01u12, D01u10] = D01_coeffs<T>(K0y, R0y, Rx0, h, D, nu);
    auto [D02u02, D02u12, D02u22, D02u01, D02u03, D02u04, D02u00, D02u13, D02u11] = D02_coeffs<T>(K0y, R0y, h, D, nu);
    auto [D0Nu0N, D0Nu1N, D0Nu2N, D0Nu0Nm1, D0Nu0Nm2, D0Nu1Nm1] = D00_coeffs<T>(K0y, R0y, KxL, RxL, h, D, nu);
    auto [D0Nm1u0Nm1, D0Nm1u1Nm1, D0Nm1u2Nm1, D0Nm1u0N, D0Nm1u0Nm2, D0Nm1u0Nm3, D0Nm1u1Nm2, D0Nm1u1N] = D01_coeffs<T>(K0y, R0y, RxL, h, D, nu);

    // Coeffs for Blk2X (ix=1 boundary)
    auto [D10u10, D10u20, D10u30, D10u00, D10u11, D10u12, D10u21, D10u01] = D10_coeffs<T>(R0y, Kx0, Rx0, h, D, nu);
    auto [D11u11, D11u12, D11u13, D11u10, D11u01, D11u21, D11u31, D11u22, D11u20, D11u00, D11u02] = D11_coeffs<T>(R0y, Rx0, h, D, nu);
    auto [D12u12, D12u13, D12u14, D12u11, D12u10, D12u02, D12u22, D12u32, D12u23, D12u21, D12u01, D12u03] = D12_coeffs<T>(R0y, h, D, nu);
    auto [D1Nu1N, D1Nu2N, D1Nu3N, D1Nu0N, D1Nu1Nm1, D1Nu1Nm2, D1Nu2Nm1, D1Nu0Nm1] = D10_coeffs<T>(R0y, KxL, RxL, h, D, nu);
    auto [D1Nm1u1Nm1, D1Nm1u1Nm2, D1Nm1u1Nm3, D1Nm1u1N, D1Nm1u0Nm1, D1Nm1u2Nm1, D1Nm1u3Nm1, D1Nm1u2Nm2, D1Nm1u2N, D1Nm1u0N, D1Nm1u0Nm2] = D11_coeffs<T>(R0y, RxL, h, D, nu);

    // Coeffs for Blk3X (ix=2+ interior)
    auto [D20u20, D20u21, D20u22, D20u10, D20u30, D20u40, D20u00, D20u31, D20u11] = D20_coeffs<T>(Kx0, Rx0, h, D, nu);
    auto [D21u21, D21u22, D21u23, D21u20, D21u11, D21u31, D21u41, D21u01, D21u32, D21u30, D21u10, D21u12] = D21_coeffs<T>(Rx0, h, D, nu);
    auto [c20_22, c11_22, c21_22, c31_22, c02_22, c12_22, c22_22, c32_22, c42_22, c13_22, c23_22, c33_22, c24_22] = D22_coeffs<T>();
    auto [D2Nu2N, D2Nu2Nm1, D2Nu2Nm2, D2Nu1N, D2Nu3N, D2Nu4N, D2Nu0N, D2Nu3Nm1, D2Nu1Nm1] = D20_coeffs<T>(KxL, RxL, h, D, nu);
    auto [D2Nm1u2Nm1, D2Nm1u2Nm2, D2Nm1u2Nm3, D2Nm1u2N, D2Nm1u1Nm1, D2Nm1u3Nm1, D2Nm1u4Nm1, D2Nm1u0Nm1, D2Nm1u3Nm2, D2Nm1u3N, D2Nm1u1N, D2Nm1u1Nm2] = D21_coeffs<T>(RxL, h, D, nu);

    // Coeffs for Blk(M-1)X (ix = N-2 boundary) - Note the different first arg (RLy)
    auto [D10u10_Mm1, D10u20_Mm1, D10u30_Mm1, D10u00_Mm1, D10u11_Mm1, D10u12_Mm1, D10u21_Mm1, D10u01_Mm1] = D10_coeffs<T>(RLy, Kx0, Rx0, h, D, nu);
    auto [D11u11_Mm1, D11u12_Mm1, D11u13_Mm1, D11u10_Mm1, D11u01_Mm1, D11u21_Mm1, D11u31_Mm1, D11u22_Mm1, D11u20_Mm1, D11u00_Mm1, D11u02_Mm1] = D11_coeffs<T>(RLy, Rx0, h, D, nu);
    auto [D12u12_Mm1, D12u13_Mm1, D12u14_Mm1, D12u11_Mm1, D12u10_Mm1, D12u02_Mm1, D12u22_Mm1, D12u32_Mm1, D12u23_Mm1, D12u21_Mm1, D12u01_Mm1, D12u03_Mm1] = D12_coeffs<T>(RLy, h, D, nu);
    auto [D1Nu1N_Mm1, D1Nu2N_Mm1, D1Nu3N_Mm1, D1Nu0N_Mm1, D1Nu1Nm1_Mm1, D1Nu1Nm2_Mm1, D1Nu2Nm1_Mm1, D1Nu0Nm1_Mm1] = D10_coeffs<T>(RLy, KxL, RxL, h, D, nu);
    auto [D1Nm1u1Nm1_Mm1, D1Nm1u1Nm2_Mm1, D1Nm1u1Nm3_Mm1, D1Nm1u1N_Mm1, D1Nm1u0Nm1_Mm1, D1Nm1u2Nm1_Mm1, D1Nm1u3Nm1_Mm1, D1Nm1u2Nm2_Mm1, D1Nm1u2N_Mm1, D1Nm1u0N_Mm1, D1Nm1u0Nm2_Mm1] = D11_coeffs<T>(RLy, RxL, h, D, nu);

    // Coeffs for BlkMX (ix = N-1 boundary) - Note the different first args (KLy, RLy)
    auto [D00u00_M, D00u10_M, D00u20_M, D00u01_M, D00u02_M, D00u11_M] = D00_coeffs<T>(KLy, RLy, Kx0, Rx0, h, D, nu);
    auto [D01u01_M, D01u11_M, D01u21_M, D01u00_M, D01u02_M, D01u03_M, D01u12_M, D01u10_M] = D01_coeffs<T>(KLy, RLy, Rx0, h, D, nu);
    auto [D02u02_M, D02u12_M, D02u22_M, D02u01_M, D02u03_M, D02u04_M, D02u00_M, D02u13_M, D02u11_M] = D02_coeffs<T>(KLy, RLy, h, D, nu);
    auto [D0Nu0N_M, D0Nu1N_M, D0Nu2N_M, D0Nu0Nm1_M, D0Nu0Nm2_M, D0Nu1Nm1_M] = D00_coeffs<T>(KLy, RLy, KxL, RxL, h, D, nu);
    auto [D0Nm1u0Nm1_M, D0Nm1u1Nm1_M, D0Nm1u2Nm1_M, D0Nm1u0N_M, D0Nm1u0Nm2_M, D0Nm1u0Nm3_M, D0Nm1u1Nm2_M, D0Nm1u1N_M] = D01_coeffs<T>(KLy, RLy, RxL, h, D, nu);

  
    // Define diagonal vector sizes (handle Ny=1, Ny=2 cases)
    int d1_size = std::max(0, Ny - 1);
    int d2_size = std::max(0, Ny - 2);

    // Empty vectors for blocks with fewer diagonals
    Eigen::Vector<T, Eigen::Dynamic> empty_D1(d1_size), empty_Dm1(d1_size), empty_D2(d2_size), empty_Dm2(d2_size);
    empty_D1.setZero(); empty_Dm1.setZero(); empty_D2.setZero(); empty_Dm2.setZero();

    // -- Blk11 (Pentadiagonal)
    Eigen::Vector<T, Eigen::Dynamic> D0_11(Ny), D1_11(d1_size), Dm1_11(d1_size), D2_11(d2_size), Dm2_11(d2_size);
    D0_11.fill(D02u02); D1_11.fill(D02u03); Dm1_11.fill(D02u01); D2_11.fill(D02u04); Dm2_11.fill(D02u00);
    // D2[[0, 1]] = [D00u02, D01u03]
    if(d2_size >= 1) D2_11(0) = D00u02;
    if(d2_size >= 2) D2_11(1) = D01u03;
    // D1[[0, 1, -1]] = [D00u01, D01u02, D0Nm1u0N]
    if(d1_size >= 1) D1_11(0) = D00u01;
    if(d1_size >= 2) D1_11(1) = D01u02;
    if(d1_size >= 1) D1_11(d1_size - 1) = D0Nm1u0N; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D00u00, D01u01, D0Nm1u0Nm1, D0Nu0N]
    if(Ny >= 1) D0_11(0) = D00u00;
    if(Ny >= 2) D0_11(1) = D01u01;
    if(Ny >= 2) D0_11(Ny - 2) = D0Nm1u0Nm1; // D0[-2]
    if(Ny >= 1) D0_11(Ny - 1) = D0Nu0N;     // D0[-1]
    // Dm1[[0, -2, -1]]   = [D01u00, D0Nm1u0Nm2, D0Nu0Nm1]
    if(d1_size >= 1) Dm1_11(0) = D01u00;
    if(d1_size >= 2) Dm1_11(d1_size - 2) = D0Nm1u0Nm2; // Dm1[-2]
    if(d1_size >= 1) Dm1_11(d1_size - 1) = D0Nu0Nm1; // Dm1[-1]
    // Dm2[[-2, -1]]      = [D0Nm1u0Nm3, D0Nu0Nm2]
    if(d2_size >= 2) Dm2_11(d2_size - 2) = D0Nm1u0Nm3; // Dm2[-2]
    if(d2_size >= 1) Dm2_11(d2_size - 1) = D0Nu0Nm2; // Dm2[-1]
    Eigen::SparseMatrix<T> Blk11 = create_block_from_diagonals(Ny, D0_11, D1_11, Dm1_11, D2_11, Dm2_11);

    // -- Blk12 (Tridiagonal)
    Eigen::Vector<T, Eigen::Dynamic> D0_12(Ny), D1_12(d1_size), Dm1_12(d1_size);
    D0_12.fill(D02u12); D1_12.fill(D02u13); Dm1_12.fill(D02u11);
    // D1[[0, 1, -1]] = [D00u11, D01u12, D0Nm1u1N]
    if(d1_size >= 1) D1_12(0) = D00u11;
    if(d1_size >= 2) D1_12(1) = D01u12;
    if(d1_size >= 1) D1_12(d1_size - 1) = D0Nm1u1N; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D00u10, D01u11, D0Nm1u1Nm1, D0Nu1N]
    if(Ny >= 1) D0_12(0) = D00u10;
    if(Ny >= 2) D0_12(1) = D01u11;
    if(Ny >= 2) D0_12(Ny - 2) = D0Nm1u1Nm1; // D0[-2]
    if(Ny >= 1) D0_12(Ny - 1) = D0Nu1N;     // D0[-1]
    // Dm1[[0, -2, -1]] = [D01u10, D0Nm1u1Nm2, D0Nu1Nm1]
    if(d1_size >= 1) Dm1_12(0) = D01u10;
    if(d1_size >= 2) Dm1_12(d1_size - 2) = D0Nm1u1Nm2; // Dm1[-2]
    if(d1_size >= 1) Dm1_12(d1_size - 1) = D0Nu1Nm1; // Dm1[-1]
    Eigen::SparseMatrix<T> Blk12 = create_block_from_diagonals(Ny, D0_12, D1_12, Dm1_12, empty_D2, empty_Dm2);

    // -- Blk13 (Diagonal)
    Eigen::Vector<T, Eigen::Dynamic> D0_13(Ny);
    D0_13.fill(D02u22);
    // D0[[0, 1, -1, -2]] = [D00u20, D01u21, D0Nu2N, D0Nm1u2Nm1]
    if(Ny >= 1) D0_13(0) = D00u20;
    if(Ny >= 2) D0_13(1) = D01u21;
    if(Ny >= 1) D0_13(Ny - 1) = D0Nu2N;     // D0[-1]
    if(Ny >= 2) D0_13(Ny - 2) = D0Nm1u2Nm1; // D0[-2]
    Eigen::SparseMatrix<T> Blk13 = create_block_from_diagonals(Ny, D0_13, empty_D1, empty_Dm1, empty_D2, empty_Dm2);

    // -- Blk21 (Tridiagonal)
    Eigen::Vector<T, Eigen::Dynamic> D0_21(Ny), D1_21(d1_size), Dm1_21(d1_size);
    D0_21.fill(D12u02); D1_21.fill(D12u03); Dm1_21.fill(D12u01);
    // D1[[0, 1, -1]] = [D10u01, D11u02, D1Nm1u0N]
    if(d1_size >= 1) D1_21(0) = D10u01;
    if(d1_size >= 2) D1_21(1) = D11u02;
    if(d1_size >= 1) D1_21(d1_size - 1) = D1Nm1u0N; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D10u00, D11u01, D1Nm1u0Nm1, D1Nu0N]
    if(Ny >= 1) D0_21(0) = D10u00;
    if(Ny >= 2) D0_21(1) = D11u01;
    if(Ny >= 2) D0_21(Ny - 2) = D1Nm1u0Nm1; // D0[-2]
    if(Ny >= 1) D0_21(Ny - 1) = D1Nu0N;     // D0[-1]
    // Dm1[[0, -2, -1]] = [D11u00, D1Nm1u0Nm2, D1Nu0Nm1]
    if(d1_size >= 1) Dm1_21(0) = D11u00; // Python uses D11u00 here, not D10u00
    if(d1_size >= 2) Dm1_21(d1_size - 2) = D1Nm1u0Nm2; // Dm1[-2]
    if(d1_size >= 1) Dm1_21(d1_size - 1) = D1Nu0Nm1; // Dm1[-1]
    Eigen::SparseMatrix<T> Blk21 = create_block_from_diagonals(Ny, D0_21, D1_21, Dm1_21, empty_D2, empty_Dm2);

    // -- Blk22 (Pentadiagonal)
    Eigen::Vector<T, Eigen::Dynamic> D0_22(Ny), D1_22(d1_size), Dm1_22(d1_size), D2_22(d2_size), Dm2_22(d2_size);
    D0_22.fill(D12u12); D1_22.fill(D12u13); Dm1_22.fill(D12u11); D2_22.fill(D12u14); Dm2_22.fill(D12u10);
    // D2[[0, 1]] = [D10u12, D11u13]
    if(d2_size >= 1) D2_22(0) = D10u12;
    if(d2_size >= 2) D2_22(1) = D11u13;
    // D1[[0, 1, -1]] = [D10u11, D11u12, D1Nm1u1N]
    if(d1_size >= 1) D1_22(0) = D10u11;
    if(d1_size >= 2) D1_22(1) = D11u12;
    if(d1_size >= 1) D1_22(d1_size - 1) = D1Nm1u1N; // D1[-1]
    // D0[[0, 1, -1, -2]] = [D10u10, D11u11, D1Nu1N, D1Nm1u1Nm1]
    if(Ny >= 1) D0_22(0) = D10u10;
    if(Ny >= 2) D0_22(1) = D11u11;
    if(Ny >= 1) D0_22(Ny - 1) = D1Nu1N;     // D0[-1]
    if(Ny >= 2) D0_22(Ny - 2) = D1Nm1u1Nm1; // D0[-2]
    // Dm1[[0, -2, -1]] = [D11u10, D1Nm1u1Nm2, D1Nu1Nm1]
    if(d1_size >= 1) Dm1_22(0) = D11u10;
    if(d1_size >= 2) Dm1_22(d1_size - 2) = D1Nm1u1Nm2; // Dm1[-2]
    if(d1_size >= 1) Dm1_22(d1_size - 1) = D1Nu1Nm1; // Dm1[-1]
    // Dm2[[-2, -1]] = [D1Nm1u1Nm3, D1Nu1Nm2]
    if(d2_size >= 2) Dm2_22(d2_size - 2) = D1Nm1u1Nm3; // Dm2[-2]
    if(d2_size >= 1) Dm2_22(d2_size - 1) = D1Nu1Nm2; // Dm2[-1]
    Eigen::SparseMatrix<T> Blk22 = create_block_from_diagonals(Ny, D0_22, D1_22, Dm1_22, D2_22, Dm2_22);

    // -- Blk23 (Tridiagonal)
    Eigen::Vector<T, Eigen::Dynamic> D0_23(Ny), D1_23(d1_size), Dm1_23(d1_size);
    D0_23.fill(D12u22); D1_23.fill(D12u23); Dm1_23.fill(D12u21);
    // D1[[0, 1, -1]] = [D10u21, D11u22, D1Nm1u2N]
    if(d1_size >= 1) D1_23(0) = D10u21;
    if(d1_size >= 2) D1_23(1) = D11u22;
    if(d1_size >= 1) D1_23(d1_size - 1) = D1Nm1u2N; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D10u20, D11u21, D1Nm1u2Nm1, D1Nu2N]
    if(Ny >= 1) D0_23(0) = D10u20;
    if(Ny >= 2) D0_23(1) = D11u21;
    if(Ny >= 2) D0_23(Ny - 2) = D1Nm1u2Nm1; // D0[-2]
    if(Ny >= 1) D0_23(Ny - 1) = D1Nu2N;     // D0[-1]
    // Dm1[[0, -2, -1]] = [D11u20, D1Nm1u2Nm2, D1Nu2Nm1]
    if(d1_size >= 1) Dm1_23(0) = D11u20;
    if(d1_size >= 2) Dm1_23(d1_size - 2) = D1Nm1u2Nm2; // Dm1[-2]
    if(d1_size >= 1) Dm1_23(d1_size - 1) = D1Nu2Nm1; // Dm1[-1]
    Eigen::SparseMatrix<T> Blk23 = create_block_from_diagonals(Ny, D0_23, D1_23, Dm1_23, empty_D2, empty_Dm2);

    // -- Blk24 (Diagonal)
    Eigen::Vector<T, Eigen::Dynamic> D0_24(Ny);
    D0_24.fill(D12u32);
    // D0[[0, 1, -1, -2]] = [D10u30, D11u31, D1Nu3N, D1Nm1u3Nm1]
    if(Ny >= 1) D0_24(0) = D10u30;
    if(Ny >= 2) D0_24(1) = D11u31;
    if(Ny >= 1) D0_24(Ny - 1) = D1Nu3N;     // D0[-1]
    if(Ny >= 2) D0_24(Ny - 2) = D1Nm1u3Nm1; // D0[-2]
    Eigen::SparseMatrix<T> Blk24 = create_block_from_diagonals(Ny, D0_24, empty_D1, empty_Dm1, empty_D2, empty_Dm2);

    // -- Blk31 (Diagonal) - Uses D2X coeffs and D22 defaults
    Eigen::Vector<T, Eigen::Dynamic> D0_31(Ny);
    D0_31.fill(c02_22); // Default from D22
    // D0[[0, 1, -1, -2]] = [D20u00, D21u01, D2Nu0N, D2Nm1u0Nm1]
    if(Ny >= 1) D0_31(0) = D20u00;
    if(Ny >= 2) D0_31(1) = D21u01;
    if(Ny >= 1) D0_31(Ny - 1) = D2Nu0N;     // D0[-1]
    if(Ny >= 2) D0_31(Ny - 2) = D2Nm1u0Nm1; // D0[-2]
    Eigen::SparseMatrix<T> Blk31 = create_block_from_diagonals(Ny, D0_31, empty_D1, empty_Dm1, empty_D2, empty_Dm2);

    // -- Blk32 (Tridiagonal) - Uses D2X coeffs and D22 defaults
    Eigen::Vector<T, Eigen::Dynamic> D0_32(Ny), D1_32(d1_size), Dm1_32(d1_size);
    D0_32.fill(c12_22); D1_32.fill(c13_22); Dm1_32.fill(c11_22); // Defaults from D22
    // D1[[0, 1, -1]] = [D20u11, D21u12, D2Nm1u1N]
    if(d1_size >= 1) D1_32(0) = D20u11;
    if(d1_size >= 2) D1_32(1) = D21u12;
    if(d1_size >= 1) D1_32(d1_size - 1) = D2Nm1u1N; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D20u10, D21u11, D2Nm1u1Nm1, D2Nu1N]
    if(Ny >= 1) D0_32(0) = D20u10;
    if(Ny >= 2) D0_32(1) = D21u11;
    if(Ny >= 2) D0_32(Ny - 2) = D2Nm1u1Nm1; // D0[-2]
    if(Ny >= 1) D0_32(Ny - 1) = D2Nu1N;     // D0[-1]
    // Dm1[[0, -2, -1]] = [D21u10, D2Nm1u1Nm2, D2Nu1Nm1]
    if(d1_size >= 1) Dm1_32(0) = D21u10;
    if(d1_size >= 2) Dm1_32(d1_size - 2) = D2Nm1u1Nm2; // Dm1[-2]
    if(d1_size >= 1) Dm1_32(d1_size - 1) = D2Nu1Nm1; // Dm1[-1]
    Eigen::SparseMatrix<T> Blk32 = create_block_from_diagonals(Ny, D0_32, D1_32, Dm1_32, empty_D2, empty_Dm2);

    // -- Blk33 (Pentadiagonal) - Uses D2X coeffs and D22 defaults
    Eigen::Vector<T, Eigen::Dynamic> D0_33(Ny), D1_33(d1_size), Dm1_33(d1_size), D2_33(d2_size), Dm2_33(d2_size);
    D0_33.fill(c22_22); D1_33.fill(c23_22); Dm1_33.fill(c21_22); D2_33.fill(c24_22); Dm2_33.fill(c20_22); // Defaults from D22
    // D2[[0, 1]] = [D20u22, D21u23]
    if(d2_size >= 1) D2_33(0) = D20u22;
    if(d2_size >= 2) D2_33(1) = D21u23;
    // D1[[0, 1, -1]] = [D20u21, D21u22, D2Nm1u2N]
    if(d1_size >= 1) D1_33(0) = D20u21;
    if(d1_size >= 2) D1_33(1) = D21u22;
    if(d1_size >= 1) D1_33(d1_size - 1) = D2Nm1u2N; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D20u20, D21u21, D2Nm1u2Nm1, D2Nu2N]
    if(Ny >= 1) D0_33(0) = D20u20;
    if(Ny >= 2) D0_33(1) = D21u21;
    if(Ny >= 2) D0_33(Ny - 2) = D2Nm1u2Nm1; // D0[-2]
    if(Ny >= 1) D0_33(Ny - 1) = D2Nu2N;     // D0[-1]
    // Dm1[[0, -2, -1]] = [D21u20, D2Nm1u2Nm2, D2Nu2Nm1]
    if(d1_size >= 1) Dm1_33(0) = D21u20;
    if(d1_size >= 2) Dm1_33(d1_size - 2) = D2Nm1u2Nm2; // Dm1[-2]
    if(d1_size >= 1) Dm1_33(d1_size - 1) = D2Nu2Nm1; // Dm1[-1]
    // Dm2[[-2, -1]] = [D2Nm1u2Nm3, D2Nu2Nm2]
    if(d2_size >= 2) Dm2_33(d2_size - 2) = D2Nm1u2Nm3; // Dm2[-2]
    if(d2_size >= 1) Dm2_33(d2_size - 1) = D2Nu2Nm2; // Dm2[-1]
    Eigen::SparseMatrix<T> Blk33 = create_block_from_diagonals(Ny, D0_33, D1_33, Dm1_33, D2_33, Dm2_33);

    // -- Blk34 (Tridiagonal) - Uses D2X coeffs and D22 defaults
    Eigen::Vector<T, Eigen::Dynamic> D0_34(Ny), D1_34(d1_size), Dm1_34(d1_size);
    D0_34.fill(c32_22); D1_34.fill(c33_22); Dm1_34.fill(c31_22); // Defaults from D22
    // D1[[0, 1, -1]] = [D20u31, D21u32, D2Nm1u3N]
    if(d1_size >= 1) D1_34(0) = D20u31;
    if(d1_size >= 2) D1_34(1) = D21u32;
    if(d1_size >= 1) D1_34(d1_size - 1) = D2Nm1u3N; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D20u30, D21u31, D2Nm1u3Nm1, D2Nu3N]
    if(Ny >= 1) D0_34(0) = D20u30;
    if(Ny >= 2) D0_34(1) = D21u31;
    if(Ny >= 2) D0_34(Ny - 2) = D2Nm1u3Nm1; // D0[-2]
    if(Ny >= 1) D0_34(Ny - 1) = D2Nu3N;     // D0[-1]
    // Dm1[[0, -2, -1]] = [D21u30, D2Nm1u3Nm2, D2Nu3Nm1]
    if(d1_size >= 1) Dm1_34(0) = D21u30;
    if(d1_size >= 2) Dm1_34(d1_size - 2) = D2Nm1u3Nm2; // Dm1[-2]
    if(d1_size >= 1) Dm1_34(d1_size - 1) = D2Nu3Nm1; // Dm1[-1]
    Eigen::SparseMatrix<T> Blk34 = create_block_from_diagonals(Ny, D0_34, D1_34, Dm1_34, empty_D2, empty_Dm2);

    // -- Blk35 (Diagonal) - Uses D2X coeffs and D22 defaults
    Eigen::Vector<T, Eigen::Dynamic> D0_35(Ny);
    D0_35.fill(c42_22); // Default from D22
    // D0[[0, 1, -2, -1]] = [D20u40, D21u41, D2Nm1u4Nm1, D2Nu4N]
    if(Ny >= 1) D0_35(0) = D20u40;
    if(Ny >= 2) D0_35(1) = D21u41;
    if(Ny >= 2) D0_35(Ny - 2) = D2Nm1u4Nm1; // D0[-2]
    if(Ny >= 1) D0_35(Ny - 1) = D2Nu4N;     // D0[-1]
    Eigen::SparseMatrix<T> Blk35 = create_block_from_diagonals(Ny, D0_35, empty_D1, empty_Dm1, empty_D2, empty_Dm2);

    // -- BlkMM (Pentadiagonal) - Uses D0X_M coeffs (ix=N-1 boundary)
    Eigen::Vector<T, Eigen::Dynamic> D0_MM(Ny), D1_MM(d1_size), Dm1_MM(d1_size), D2_MM(d2_size), Dm2_MM(d2_size);
    D0_MM.fill(D02u02_M); D1_MM.fill(D02u03_M); Dm1_MM.fill(D02u01_M); D2_MM.fill(D02u04_M); Dm2_MM.fill(D02u00_M);
    // D2[[0, 1]] = [D00u02, D01u03] -> Uses _M coeffs
    if(d2_size >= 1) D2_MM(0) = D00u02_M;
    if(d2_size >= 2) D2_MM(1) = D01u03_M;
    // D1[[0, 1, -1]] = [D00u01, D01u02, D0Nm1u0N] -> Uses _M coeffs
    if(d1_size >= 1) D1_MM(0) = D00u01_M;
    if(d1_size >= 2) D1_MM(1) = D01u02_M;
    if(d1_size >= 1) D1_MM(d1_size - 1) = D0Nm1u0N_M; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D00u00, D01u01, D0Nm1u0Nm1, D0Nu0N] -> Uses _M coeffs
    if(Ny >= 1) D0_MM(0) = D00u00_M;
    if(Ny >= 2) D0_MM(1) = D01u01_M;
    if(Ny >= 2) D0_MM(Ny - 2) = D0Nm1u0Nm1_M; // D0[-2]
    if(Ny >= 1) D0_MM(Ny - 1) = D0Nu0N_M;     // D0[-1]
    // Dm1[[0, -2, -1]] = [D01u00, D0Nm1u0Nm2, D0Nu0Nm1] -> Uses _M coeffs
    if(d1_size >= 1) Dm1_MM(0) = D01u00_M;
    if(d1_size >= 2) Dm1_MM(d1_size - 2) = D0Nm1u0Nm2_M; // Dm1[-2]
    if(d1_size >= 1) Dm1_MM(d1_size - 1) = D0Nu0Nm1_M; // Dm1[-1]
    // Dm2[[-2, -1]] = [D0Nm1u0Nm3, D0Nu0Nm2] -> Uses _M coeffs
    if(d2_size >= 2) Dm2_MM(d2_size - 2) = D0Nm1u0Nm3_M; // Dm2[-2]
    if(d2_size >= 1) Dm2_MM(d2_size - 1) = D0Nu0Nm2_M; // Dm2[-1]
    Eigen::SparseMatrix<T> BlkMM = create_block_from_diagonals(Ny, D0_MM, D1_MM, Dm1_MM, D2_MM, Dm2_MM);

    // -- BlkMMm1 (Tridiagonal) - Uses D0X_M coeffs (ix=N-1 boundary)
    Eigen::Vector<T, Eigen::Dynamic> D0_MMm1(Ny), D1_MMm1(d1_size), Dm1_MMm1(d1_size);
    D0_MMm1.fill(D02u12_M); D1_MMm1.fill(D02u13_M); Dm1_MMm1.fill(D02u11_M);
    // D1[[0, 1, -1]] = [D00u11, D01u12, D0Nm1u1N] -> Uses _M coeffs
    if(d1_size >= 1) D1_MMm1(0) = D00u11_M;
    if(d1_size >= 2) D1_MMm1(1) = D01u12_M;
    if(d1_size >= 1) D1_MMm1(d1_size - 1) = D0Nm1u1N_M; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D00u10, D01u11, D0Nm1u1Nm1, D0Nu1N] -> Uses _M coeffs
    if(Ny >= 1) D0_MMm1(0) = D00u10_M;
    if(Ny >= 2) D0_MMm1(1) = D01u11_M;
    if(Ny >= 2) D0_MMm1(Ny - 2) = D0Nm1u1Nm1_M; // D0[-2]
    if(Ny >= 1) D0_MMm1(Ny - 1) = D0Nu1N_M;     // D0[-1]
    // Dm1[[0, -2, -1]] = [D01u10, D0Nm1u1Nm2, D0Nu1Nm1] -> Uses _M coeffs
    if(d1_size >= 1) Dm1_MMm1(0) = D01u10_M;
    if(d1_size >= 2) Dm1_MMm1(d1_size - 2) = D0Nm1u1Nm2_M; // Dm1[-2]
    if(d1_size >= 1) Dm1_MMm1(d1_size - 1) = D0Nu1Nm1_M; // Dm1[-1]
    Eigen::SparseMatrix<T> BlkMMm1 = create_block_from_diagonals(Ny, D0_MMm1, D1_MMm1, Dm1_MMm1, empty_D2, empty_Dm2);

    // -- BlkMMm2 (Diagonal) - Uses D0X_M coeffs (ix=N-1 boundary)
    Eigen::Vector<T, Eigen::Dynamic> D0_MMm2(Ny);
    D0_MMm2.fill(D02u22_M);
    // D0[[0, 1, -1, -2]] = [D00u20, D01u21, D0Nu2N, D0Nm1u2Nm1] -> Uses _M coeffs
    if(Ny >= 1) D0_MMm2(0) = D00u20_M;
    if(Ny >= 2) D0_MMm2(1) = D01u21_M;
    if(Ny >= 1) D0_MMm2(Ny - 1) = D0Nu2N_M;     // D0[-1]
    if(Ny >= 2) D0_MMm2(Ny - 2) = D0Nm1u2Nm1_M; // D0[-2]
    Eigen::SparseMatrix<T> BlkMMm2 = create_block_from_diagonals(Ny, D0_MMm2, empty_D1, empty_Dm1, empty_D2, empty_Dm2);

    // -- BlkMm1M (Tridiagonal) - Uses D1X_Mm1 coeffs (ix=N-2 boundary)
    Eigen::Vector<T, Eigen::Dynamic> D0_Mm1M(Ny), D1_Mm1M(d1_size), Dm1_Mm1M(d1_size);
    D0_Mm1M.fill(D12u02_Mm1); D1_Mm1M.fill(D12u03_Mm1); Dm1_Mm1M.fill(D12u01_Mm1);
    // D1[[0, 1, -1]] = [D10u01, D11u02, D1Nm1u0N] -> Uses _Mm1 coeffs
    if(d1_size >= 1) D1_Mm1M(0) = D10u01_Mm1;
    if(d1_size >= 2) D1_Mm1M(1) = D11u02_Mm1;
    if(d1_size >= 1) D1_Mm1M(d1_size - 1) = D1Nm1u0N_Mm1; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D10u00, D11u01, D1Nm1u0Nm1, D1Nu0N] -> Uses _Mm1 coeffs
    if(Ny >= 1) D0_Mm1M(0) = D10u00_Mm1;
    if(Ny >= 2) D0_Mm1M(1) = D11u01_Mm1;
    if(Ny >= 2) D0_Mm1M(Ny - 2) = D1Nm1u0Nm1_Mm1; // D0[-2]
    if(Ny >= 1) D0_Mm1M(Ny - 1) = D1Nu0N_Mm1;     // D0[-1]
    // Dm1[[0, -2, -1]] = [D11u00, D1Nm1u0Nm2, D1Nu0Nm1] -> Uses _Mm1 coeffs
    if(d1_size >= 1) Dm1_Mm1M(0) = D11u00_Mm1;
    if(d1_size >= 2) Dm1_Mm1M(d1_size - 2) = D1Nm1u0Nm2_Mm1; // Dm1[-2]
    if(d1_size >= 1) Dm1_Mm1M(d1_size - 1) = D1Nu0Nm1_Mm1; // Dm1[-1]
    Eigen::SparseMatrix<T> BlkMm1M = create_block_from_diagonals(Ny, D0_Mm1M, D1_Mm1M, Dm1_Mm1M, empty_D2, empty_Dm2);

    // -- BlkMm1Mm1 (Pentadiagonal) - Uses D1X_Mm1 coeffs (ix=N-2 boundary)
    Eigen::Vector<T, Eigen::Dynamic> D0_Mm1Mm1(Ny), D1_Mm1Mm1(d1_size), Dm1_Mm1Mm1(d1_size), D2_Mm1Mm1(d2_size), Dm2_Mm1Mm1(d2_size);
    D0_Mm1Mm1.fill(D12u12_Mm1); D1_Mm1Mm1.fill(D12u13_Mm1); Dm1_Mm1Mm1.fill(D12u11_Mm1); D2_Mm1Mm1.fill(D12u14_Mm1); Dm2_Mm1Mm1.fill(D12u10_Mm1);
    // D2[[0, 1]] = [D10u12, D11u13] -> Uses _Mm1 coeffs
    if(d2_size >= 1) D2_Mm1Mm1(0) = D10u12_Mm1;
    if(d2_size >= 2) D2_Mm1Mm1(1) = D11u13_Mm1;
    // D1[[0, 1, -1]] = [D10u11, D11u12, D1Nm1u1N] -> Uses _Mm1 coeffs
    if(d1_size >= 1) D1_Mm1Mm1(0) = D10u11_Mm1;
    if(d1_size >= 2) D1_Mm1Mm1(1) = D11u12_Mm1;
    if(d1_size >= 1) D1_Mm1Mm1(d1_size - 1) = D1Nm1u1N_Mm1; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D10u10, D11u11, D1Nm1u1Nm1, D1Nu1N] -> Uses _Mm1 coeffs
    if(Ny >= 1) D0_Mm1Mm1(0) = D10u10_Mm1;
    if(Ny >= 2) D0_Mm1Mm1(1) = D11u11_Mm1;
    if(Ny >= 2) D0_Mm1Mm1(Ny - 2) = D1Nm1u1Nm1_Mm1; // D0[-2]
    if(Ny >= 1) D0_Mm1Mm1(Ny - 1) = D1Nu1N_Mm1;     // D0[-1]
    // Dm1[[0, -2, -1]] = [D11u10, D1Nm1u1Nm2, D1Nu1Nm1] -> Uses _Mm1 coeffs
    if(d1_size >= 1) Dm1_Mm1Mm1(0) = D11u10_Mm1;
    if(d1_size >= 2) Dm1_Mm1Mm1(d1_size - 2) = D1Nm1u1Nm2_Mm1; // Dm1[-2]
    if(d1_size >= 1) Dm1_Mm1Mm1(d1_size - 1) = D1Nu1Nm1_Mm1; // Dm1[-1]
    // Dm2[[-2, -1]] = [D1Nm1u1Nm3, D1Nu1Nm2] -> Uses _Mm1 coeffs
    if(d2_size >= 2) Dm2_Mm1Mm1(d2_size - 2) = D1Nm1u1Nm3_Mm1; // Dm2[-2]
    if(d2_size >= 1) Dm2_Mm1Mm1(d2_size - 1) = D1Nu1Nm2_Mm1; // Dm2[-1]
    Eigen::SparseMatrix<T> BlkMm1Mm1 = create_block_from_diagonals(Ny, D0_Mm1Mm1, D1_Mm1Mm1, Dm1_Mm1Mm1, D2_Mm1Mm1, Dm2_Mm1Mm1);

    // -- BlkMm1Mm2 (Tridiagonal) - Uses D1X_Mm1 coeffs (ix=N-2 boundary)
    Eigen::Vector<T, Eigen::Dynamic> D0_Mm1Mm2(Ny), D1_Mm1Mm2(d1_size), Dm1_Mm1Mm2(d1_size);
    D0_Mm1Mm2.fill(D12u22_Mm1); D1_Mm1Mm2.fill(D12u23_Mm1); Dm1_Mm1Mm2.fill(D12u21_Mm1);
    // D1[[0, 1, -1]] = [D10u21, D11u22, D1Nm1u2N] -> Uses _Mm1 coeffs
    if(d1_size >= 1) D1_Mm1Mm2(0) = D10u21_Mm1;
    if(d1_size >= 2) D1_Mm1Mm2(1) = D11u22_Mm1;
    if(d1_size >= 1) D1_Mm1Mm2(d1_size - 1) = D1Nm1u2N_Mm1; // D1[-1]
    // D0[[0, 1, -2, -1]] = [D10u20, D11u21, D1Nm1u2Nm1, D1Nu2N] -> Uses _Mm1 coeffs
    if(Ny >= 1) D0_Mm1Mm2(0) = D10u20_Mm1;
    if(Ny >= 2) D0_Mm1Mm2(1) = D11u21_Mm1;
    if(Ny >= 2) D0_Mm1Mm2(Ny - 2) = D1Nm1u2Nm1_Mm1; // D0[-2]
    if(Ny >= 1) D0_Mm1Mm2(Ny - 1) = D1Nu2N_Mm1;     // D0[-1]
    // Dm1[[0, -2, -1]] = [D11u20, D1Nm1u2Nm2, D1Nu2Nm1] -> Uses _Mm1 coeffs
    if(d1_size >= 1) Dm1_Mm1Mm2(0) = D11u20_Mm1;
    if(d1_size >= 2) Dm1_Mm1Mm2(d1_size - 2) = D1Nm1u2Nm2_Mm1; // Dm1[-2]
    if(d1_size >= 1) Dm1_Mm1Mm2(d1_size - 1) = D1Nu2Nm1_Mm1; // Dm1[-1]
    Eigen::SparseMatrix<T> BlkMm1Mm2 = create_block_from_diagonals(Ny, D0_Mm1Mm2, D1_Mm1Mm2, Dm1_Mm1Mm2, empty_D2, empty_Dm2);

    // -- BlkMm1Mm3 (Diagonal) - Uses D1X_Mm1 coeffs (ix=N-2 boundary)
    Eigen::Vector<T, Eigen::Dynamic> D0_Mm1Mm3(Ny);
    D0_Mm1Mm3.fill(D12u32_Mm1);
    // D0[[0, 1, -1, -2]] = [D10u30, D11u31, D1Nu3N, D1Nm1u3Nm1] -> Uses _Mm1 coeffs
    if(Ny >= 1) D0_Mm1Mm3(0) = D10u30_Mm1;
    if(Ny >= 2) D0_Mm1Mm3(1) = D11u31_Mm1;
    if(Ny >= 1) D0_Mm1Mm3(Ny - 1) = D1Nu3N_Mm1;     // D0[-1]
    if(Ny >= 2) D0_Mm1Mm3(Ny - 2) = D1Nm1u3Nm1_Mm1; // D0[-2]
    Eigen::SparseMatrix<T> BlkMm1Mm3 = create_block_from_diagonals(Ny, D0_Mm1Mm3, empty_D1, empty_Dm1, empty_D2, empty_Dm2);


    // --- Assemble Final Matrix from Blocks ---
    std::vector<Eigen::Triplet<T>> final_triplets;
    // Estimate required capacity (can be tuned - 13 is max non-zeros per row in interior)
    final_triplets.reserve(static_cast<size_t>(MM) * 13); 

    // Add first block row (ix=0)
    add_block_to_triplets(0, 0, Blk11, final_triplets, Ny, MM);
    if (Nx > 1) add_block_to_triplets(0, 1, Blk12, final_triplets, Ny, MM);
    if (Nx > 2) add_block_to_triplets(0, 2, Blk13, final_triplets, Ny, MM);

    // Add second block row (ix=1)
    if (Nx > 1) {
        add_block_to_triplets(1, 0, Blk21, final_triplets, Ny, MM);
        add_block_to_triplets(1, 1, Blk22, final_triplets, Ny, MM);
        if (Nx > 2) add_block_to_triplets(1, 2, Blk23, final_triplets, Ny, MM);
        if (Nx > 3) add_block_to_triplets(1, 3, Blk24, final_triplets, Ny, MM);
    }

    // Add interior block rows (ix = 2 to Nx-3)
    for (int ix = 2; ix < Nx - 2; ++ix) {
        add_block_to_triplets(ix, ix - 2, Blk31, final_triplets, Ny, MM);
        add_block_to_triplets(ix, ix - 1, Blk32, final_triplets, Ny, MM);
        add_block_to_triplets(ix, ix    , Blk33, final_triplets, Ny, MM);
        add_block_to_triplets(ix, ix + 1, Blk34, final_triplets, Ny, MM);
        add_block_to_triplets(ix, ix + 2, Blk35, final_triplets, Ny, MM);
    }

     // Add second to last block row (ix=N-2) - Needs Nx >= 4 for BlkMm1Mm3, Nx >= 3 for others
     if (Nx >= 3) { // Check if ix=N-2 exists
         if (Nx >= 4) add_block_to_triplets(Nx - 2, Nx - 4, BlkMm1Mm3, final_triplets, Ny, MM); // ix = N-2, col = N-4
         add_block_to_triplets(Nx - 2, Nx - 3, BlkMm1Mm2, final_triplets, Ny, MM); // ix = N-2, col = N-3
         add_block_to_triplets(Nx - 2, Nx - 2, BlkMm1Mm1, final_triplets, Ny, MM); // ix = N-2, col = N-2
         add_block_to_triplets(Nx - 2, Nx - 1, BlkMm1M,   final_triplets, Ny, MM); // ix = N-2, col = N-1
     }
     
     // Add last block row (ix=N-1) - Needs Nx >= 3 for BlkMMm2, Nx >= 2 for others
     if (Nx >= 2) { // Check if ix=N-1 exists
        if (Nx >= 3) add_block_to_triplets(Nx - 1, Nx - 3, BlkMMm2, final_triplets, Ny, MM); // ix = N-1, col = N-3
        add_block_to_triplets(Nx - 1, Nx - 2, BlkMMm1, final_triplets, Ny, MM); // ix = N-1, col = N-2
        add_block_to_triplets(Nx - 1, Nx - 1, BlkMM,   final_triplets, Ny, MM); // ix = N-1, col = N-1
     } else if (Nx == 1) {
         // Handle Nx=1 case separately if needed, though the assembly loop for ix=0 already added Blk11.
         // If Blk11 is the only block, it should be correct.
     }


    Eigen::SparseMatrix<T> biharm(MM, MM);
    biharm.setFromTriplets(final_triplets.begin(), final_triplets.end());

    // --- Scale and Return ---
    T h_pow_4 = std::pow(h, static_cast<T>(4));
    if (h_pow_4 == static_cast<T>(0.0)) {
          throw std::runtime_error("h^4 is too small, potential division by zero.");
    }

    return biharm / h_pow_4;
}
