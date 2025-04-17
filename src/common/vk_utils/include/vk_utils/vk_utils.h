#pragma once

#include <Eigen/Dense>

template <typename T>
Eigen::VectorX<T> calculate_nonlinear_term_eigen(
    const Eigen::MatrixX<T>& H_reshaped,
    const Eigen::VectorX<T>& q
)
{
    // Calculate t0 = H_reshaped * q
    Eigen::VectorX<T> t0_flat = H_reshaped * q;

    // Reshape t0 to n_modes x n_modes
    Eigen::MatrixX<T> t0 =
        Eigen::Map<const Eigen::MatrixX<T>>(t0_flat.data(), q.size(), q.size());

    // Calculate t2 = t0 * q
    Eigen::VectorX<T> t2 = t0 * q;

    // Calculate nl = t0^T * t2
    Eigen::VectorX<T> nl = t0.transpose() * t2;

    return nl;
}

// Eigen implementation of tension modulation nonlinear term calculation
template <typename T>
Eigen::VectorX<T> calculate_tm_nonlinear_term_eigen(
    const Eigen::VectorX<T>& lambda_mu,
    const Eigen::VectorX<T>& tau_with_norms,
    const Eigen::VectorX<T>& q
)
{
    // Calculate q^2
    Eigen::VectorX<T> q_squared = q.array().square();
    
    // Calculate tau_with_norms DOT q^2 (dot product resulting in a scalar)
    T tmp_scalar = tau_with_norms.dot(q_squared);
    
    // Calculate lambda_mu .* q .* (scalar result of dot product)
    Eigen::VectorX<T> nl = lambda_mu.cwiseProduct(q) * tmp_scalar;
    
    return nl;
}

// Calculate coefficients for time integration with transfer function approach
// Calculate coefficients for time integration with transfer function approach
template <typename T>
void calculate_coefficients_tf(
    const Eigen::VectorX<T>& gamma2_mu,
    const Eigen::VectorX<T>& omega_mu_squared,
    Eigen::VectorX<T>& B,
    Eigen::VectorX<T>& C,
    Eigen::VectorX<T>& A_inv,
    T dt
) {
    Eigen::VectorX<T> gamma_mu = gamma2_mu.array() * 0.5;
    Eigen::VectorX<T> omega_mu_damped = (omega_mu_squared.array() - gamma_mu.array().square()).sqrt();
    Eigen::VectorX<T> radii = (-gamma_mu.array() * dt).exp();
    Eigen::VectorX<T> imag = radii.array() * (omega_mu_damped.array() * dt).sin();
    Eigen::VectorX<T> real = radii.array() * (omega_mu_damped.array() * dt).cos();

    B = 2.0 * real;
    C = -radii.array().square().matrix();
    A_inv = ((dt * imag.array()) / omega_mu_damped.array()).matrix();
}


// Calculate coefficients for time integration with state-variable approach
template <typename T>
void calculate_coefficients_sv(
    const Eigen::VectorX<T>& gamma2_mu,
    const Eigen::VectorX<T>& omega_mu_squared,
    Eigen::VectorX<T>& B,
    Eigen::VectorX<T>& C,
    Eigen::VectorX<T>& A_inv,
    T dt  // time step
) {
    T dt_squared = dt * dt;
    // Calculate A_inv
    A_inv = (2.0 * dt_squared) / (2.0 + gamma2_mu.array() * dt);
    
    // Calculate B
    B = (Eigen::VectorX<T>::Constant(
        gamma2_mu.size(), 2.0 / dt_squared) - omega_mu_squared).cwiseProduct(A_inv);
    
    // Calculate C
    C = (Eigen::VectorX<T>::Constant(
        gamma2_mu.size(), -1.0 / dt_squared)
        + gamma2_mu * (1.0 / (2.0 * dt))).cwiseProduct(A_inv);
}


template <typename T>
Eigen::VectorX<T> A_inv_vector(
    T h,  // temporal grid spacing
    const Eigen::VectorX<T>& damping  // damping term vector
) {
    // Result = (2.0 * h * h) / (2.0 + damping * h)
    return (2.0 * h * h) / (2.0 + damping.array() * h);
}

template <typename T>
Eigen::VectorX<T> B_vector(
    T h,  // temporal grid spacing
    const Eigen::VectorX<T>& stiffness  // stiffness term vector
) {
    // Create a constant vector and subtract stiffness
    return Eigen::VectorX<T>::Constant(
        stiffness.size(), 2.0 / (h * h)) - stiffness;
}

template <typename T>
Eigen::VectorX<T> C_vector(
    T h,  // temporal grid spacing
    const Eigen::VectorX<T>& damping  // damping term vector
) {
    // Constant vector for -1.0/(h*h)
    Eigen::VectorX<T> result = 
        Eigen::VectorX<T>::Constant(damping.size(), -1.0 / (h * h));
    
    // Add scaled damping term
    return result + damping * (1.0 / (2.0 * h));
}