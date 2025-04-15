#pragma once

// #include <Eigen/Dense>
#include <functional>
#include <vector>
#include "vk_utils.h"

namespace time_integrators {

/**
 * Implements the Störmer-Verlet time integration scheme for modal excitation.
 * 
 * @tparam T Floating point type (float or double)
 * @param gamma2_mu Damping term vector (n_modes)
 * @param omega_mu_squared Stiffness term vector (n_modes)
 * @param modal_excitation Modal excitation over time (n_steps x n_modes)
 * @param dt Time step
 * @param nl_fn Optional nonlinear function (can be nullptr if not used)
 * @return Pair containing the final state and the full solution
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
solve_sv_excitation(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& gamma2_mu,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& omega_mu_squared,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& modal_excitation,
    T dt,
    std::function<Eigen::Matrix<T, Eigen::Dynamic, 1>(const Eigen::Matrix<T, Eigen::Dynamic, 1>&)> nl_fn = nullptr
) {
    // Compute integration coefficient vectors
    Eigen::Matrix<T, Eigen::Dynamic, 1> A_inv = A_inv_vector(dt, gamma2_mu);
    Eigen::Matrix<T, Eigen::Dynamic, 1> B = B_vector(dt, omega_mu_squared).array() * A_inv.array();
    Eigen::Matrix<T, Eigen::Dynamic, 1> C = C_vector(dt, gamma2_mu).array() * A_inv.array();
    
    int n_modes = A_inv.size();
    int n_steps = modal_excitation.rows();
    
    // Initialize state vectors
    Eigen::Matrix<T, Eigen::Dynamic, 1> q = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n_modes);
    Eigen::Matrix<T, Eigen::Dynamic, 1> q_prev = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n_modes);
    
    // Initialize solution matrix
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> solution(n_steps, n_modes);
    solution.row(0) = q.transpose();
    
    // If no nonlinear function is provided, use a dummy function that returns zeros
    if (!nl_fn) {
        nl_fn = [n_modes](const Eigen::Matrix<T, Eigen::Dynamic, 1>& q) {
            return Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n_modes);
        };
    }
    
    // Main time stepping loop
    for (int t = 1; t < n_steps; ++t) {
        // Get current excitation
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& x = modal_excitation.row(t).transpose();
        
        // Compute nonlinear term
        Eigen::Matrix<T, Eigen::Dynamic, 1> nl = nl_fn(q);
        
        // Update state
        Eigen::Matrix<T, Eigen::Dynamic, 1> q_next = 
            B.array() * q.array() + 
            C.array() * q_prev.array() - 
            A_inv.array() * nl.array() + 
            A_inv.array() * x.array();
        
        // Store current solution
        solution.row(t) = q_next.transpose();
        
        // Update for next iteration
        q_prev = q;
        q = q_next;
    }
    
    // Return the modal solution
    return solution;
}

/**
 * Variant of the Störmer-Verlet solver for modal simulation without external excitation.
 * 
 * @tparam T Floating point type (float or double)
 * @param gamma2_mu Damping term vector
 * @param omega_mu_squared Stiffness term vector
 * @param initial_displacement Initial modal displacement
 * @param initial_velocity Initial modal velocity (can be zero)
 * @param n_steps Number of time steps to simulate
 * @param dt Time step
 * @param nl_fn Optional nonlinear function (can be nullptr if not used)
 * @return Matrix containing the full solution (n_steps x n_modes)
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
solve_sv_initial_conditions(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& gamma2_mu,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& omega_mu_squared,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& initial_displacement,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& initial_velocity,
    int n_steps,
    T dt,
    std::function<Eigen::Matrix<T, Eigen::Dynamic, 1>(const Eigen::Matrix<T, Eigen::Dynamic, 1>&)> nl_fn = nullptr
) {
    // Compute integration coefficient vectors
    Eigen::Matrix<T, Eigen::Dynamic, 1> A_inv = A_inv_vector(dt, gamma2_mu);
    Eigen::Matrix<T, Eigen::Dynamic, 1> B = B_vector(dt, omega_mu_squared).array() * A_inv.array();
    Eigen::Matrix<T, Eigen::Dynamic, 1> C = C_vector(dt, gamma2_mu).array() * A_inv.array();
    
    int n_modes = initial_displacement.size();
    
    // Initialize state vectors
    Eigen::Matrix<T, Eigen::Dynamic, 1> q = initial_displacement;
    
    // For q_prev, use a backward Euler step from the initial conditions
    Eigen::Matrix<T, Eigen::Dynamic, 1> q_prev = q - dt * initial_velocity;
    
    // Initialize solution matrix
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> solution(n_steps, n_modes);
    solution.row(0) = q.transpose();
    
    // If no nonlinear function is provided, use a dummy function that returns zeros
    if (!nl_fn) {
        nl_fn = [n_modes](const Eigen::Matrix<T, Eigen::Dynamic, 1>& q) {
            return Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n_modes);
        };
    }
    
    // Main time stepping loop
    for (int t = 1; t < n_steps; ++t) {
        // Compute nonlinear term
        Eigen::Matrix<T, Eigen::Dynamic, 1> nl = nl_fn(q);
        
        // Update state (no external excitation)
        Eigen::Matrix<T, Eigen::Dynamic, 1> q_next = 
            B.array() * q.array() + 
            C.array() * q_prev.array() - 
            A_inv.array() * nl.array();
        
        // Store current solution
        solution.row(t) = q_next.transpose();
        
        // Update for next iteration
        q_prev = q;
        q = q_next;
    }
    
    return solution;
}

} // namespace time_integrators
