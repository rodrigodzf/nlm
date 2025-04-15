#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace ftm {

/**
 * Compute the eigenvalues of a string with fixed ends.
 * 
 * The eigenvalues are given by: λ_μ = (μ * π/L)²
 * where μ is the mode number and L is the length of the string.
 * 
 * @param n_modes Number of modes to compute
 * @param length Length of the string
 * @return Eigen vector containing the eigenvalues
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> string_eigenvalues(int n_modes, T length) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> eigenvalues(n_modes);
    for (int mu = 1; mu <= n_modes; ++mu) {
        eigenvalues(mu-1) = std::pow(mu * M_PI / length, 2);
    }
    return eigenvalues;
}

/**
 * Evaluate string eigenfunctions at a specific position.
 * 
 * @param indices Selected mode indices (1-based)
 * @param position Position to evaluate the eigenfunctions
 * @param length Length of the string
 * @return Vector of eigenfunction values at the given position
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> evaluate_string_eigenfunctions(
    const Eigen::Matrix<int, Eigen::Dynamic, 1>& indices,
    T position,
    T length
) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> values(indices.size());
    for (int i = 0; i < indices.size(); ++i) {
        values(i) = std::sin(indices(i) * M_PI * position / length);
    }
    return values;
}

/**
 * Compute the wavenumbers of a rectangular plate with clamped edges.
 * 
 * The wavenumbers are given by:
 * k_x = μ * π/L₁, k_y = ν * π/L₂
 * where μ and ν are the mode numbers, and L₁ and L₂ are the plate dimensions.
 * 
 * @param n_max_modes_x Number of modes in x direction
 * @param n_max_modes_y Number of modes in y direction
 * @param l1 Width of the plate
 * @param l2 Height of the plate
 * @return Pair of Eigen vectors containing wavenumbers in x and y directions
 */
template <typename T>
std::pair<Eigen::Matrix<T, Eigen::Dynamic, 1>, Eigen::Matrix<T, Eigen::Dynamic, 1>> 
plate_wavenumbers(int n_max_modes_x, int n_max_modes_y, T l1, T l2) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> wavenumbers_x(n_max_modes_x);
    Eigen::Matrix<T, Eigen::Dynamic, 1> wavenumbers_y(n_max_modes_y);
    
    for (int mu = 1; mu <= n_max_modes_x; ++mu) {
        wavenumbers_x(mu-1) = mu * M_PI / l1;
    }
    
    for (int nu = 1; nu <= n_max_modes_y; ++nu) {
        wavenumbers_y(nu-1) = nu * M_PI / l2;
    }
    
    return {wavenumbers_x, wavenumbers_y};
}

/**
 * Compute the eigenvalues of a rectangular plate.
 * 
 * The eigenvalues are given by:
 * λ_{μ,ν} = (μ * π/L₁)² + (ν * π/L₂)²
 * 
 * @param wavenumbers_x Wavenumbers in x direction
 * @param wavenumbers_y Wavenumbers in y direction
 * @return Matrix of eigenvalues
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
plate_eigenvalues(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& wavenumbers_x,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& wavenumbers_y
) {
    int n_modes_x = wavenumbers_x.size();
    int n_modes_y = wavenumbers_y.size();
    
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigenvalues(n_modes_x, n_modes_y);
    
    for (int i = 0; i < n_modes_x; ++i) {
        for (int j = 0; j < n_modes_y; ++j) {
            eigenvalues(i, j) = std::pow(wavenumbers_x(i), 2) + std::pow(wavenumbers_y(j), 2);
        }
    }
    
    return eigenvalues;
}

/**
 * Evaluate rectangular plate eigenfunctions at a specific position.
 * 
 * @param mn_indices Matrix of mode indices (m,n pairs) for each mode
 * @param position Position vector [x, y] to evaluate the eigenfunctions
 * @param l1 Width of the plate
 * @param l2 Height of the plate
 * @return Vector of eigenfunction values at the given position
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> evaluate_rectangular_eigenfunctions(
    const Eigen::Matrix<int, Eigen::Dynamic, 2>& mn_indices,
    const Eigen::Matrix<T, 2, 1>& position,
    T l1, 
    T l2
) {
    int n_modes = mn_indices.rows();
    Eigen::Matrix<T, Eigen::Dynamic, 1> values(n_modes);
    
    for (int i = 0; i < n_modes; ++i) {
        T sin_x = std::sin(mn_indices(i, 0) * M_PI * position(0) / l1);
        T sin_y = std::sin(mn_indices(i, 1) * M_PI * position(1) / l2);
        values(i) = sin_x * sin_y;
    }
    
    return values;
}

} // namespace ftm
