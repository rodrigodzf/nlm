#pragma once

#include <Eigen/Dense>

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> calculate_nonlinear_term_eigen(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H_reshaped,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& q
)
{
    // Calculate t0 = H_reshaped * q
    Eigen::Matrix<T, Eigen::Dynamic, 1> t0_flat = H_reshaped * q;

    // Reshape t0 to n_modes x n_modes
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> t0 =
        Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(t0_flat.data(), q.size(), q.size());

    // Calculate t2 = t0 * q
    Eigen::Matrix<T, Eigen::Dynamic, 1> t2 = t0 * q;

    // Calculate nl = t0^T * t2
    Eigen::Matrix<T, Eigen::Dynamic, 1> nl = t0.transpose() * t2;

    return nl;
}

// Eigen implementation of tension modulation nonlinear term calculation
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> calculate_tm_nonlinear_term_eigen(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& lambda_mu,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& tau_with_norms,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& q
)
{
    // Calculate q^2
    Eigen::Matrix<T, Eigen::Dynamic, 1> q_squared = q.array().square();
    
    // Calculate tau_with_norms DOT q^2 (dot product resulting in a scalar)
    T tmp_scalar = tau_with_norms.dot(q_squared);
    
    // Calculate lambda_mu .* q .* (scalar result of dot product)
    Eigen::Matrix<T, Eigen::Dynamic, 1> nl = lambda_mu.cwiseProduct(q) * tmp_scalar;
    
    return nl;
}

// Calculate A_inv vector for time integration
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> A_inv_vector(
    T h,  // temporal grid spacing
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& damping  // damping term vector
) {
    // Result = (2.0 * h * h) / (2.0 + damping * h)
    return (2.0 * h * h) / (2.0 + damping.array() * h);
}

// Calculate B vector for time integration
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> B_vector(
    T h,  // temporal grid spacing
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& stiffness  // stiffness term vector
) {
    // Create a constant vector and subtract stiffness
    return Eigen::Matrix<T, Eigen::Dynamic, 1>::Constant(
        stiffness.size(), 2.0 / (h * h)) - stiffness;
}

// Calculate C vector for time integration
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> C_vector(
    T h,  // temporal grid spacing
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& damping  // damping term vector
) {
    // Constant vector for -1.0/(h*h)
    Eigen::Matrix<T, Eigen::Dynamic, 1> result = 
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Constant(damping.size(), -1.0 / (h * h));
    
    // Add scaled damping term
    return result + damping * (1.0 / (2.0 * h));
}
