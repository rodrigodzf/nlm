#pragma once

#include <Eigen/Dense>
#include <numeric> // For std::iota

template <typename T>
void calculate_nonlinear_berger(
    const Eigen::VectorX<T>& lambda_mu,
    const Eigen::VectorX<T>& tau_with_norms,
    const Eigen::VectorX<T>& q,
    Eigen::VectorX<T>& nl
)
{
    // dot product of tau_with_norms and q^2
    T scalar = (tau_with_norms.array() * q.array().square()).sum();
    
    // Apply scalar to the product
    nl.noalias() = lambda_mu.cwiseProduct(q) * scalar;
}


template <typename T>
void calculate_nonlinear_vk(
    const Eigen::MatrixX<T>& H_scaled,
    const Eigen::VectorX<T>& q,
    Eigen::VectorX<T>& nl,
    const int n_psi,
    const int n_phi
)
{
    auto t0 = Eigen::Map<const Eigen::MatrixX<T>>((H_scaled * q).eval().data(), n_psi, n_phi);
    nl.noalias() = t0.transpose() * (t0 * q);
}

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

template <typename T>
Eigen::MatrixX<T> calculate_plate_eigenvalues(
    int n_max_modes_x,
    int n_max_modes_y,
    T l1,
    T l2
) {
    // Create a matrix to store the eigenvalues
    Eigen::MatrixX<T> lambda_mu_2d(n_max_modes_y, n_max_modes_x);
    
    // Calculate eigenvalues directly
    for (int mu_y = 1; mu_y <= n_max_modes_y; ++mu_y) {
        for (int mu_x = 1; mu_x <= n_max_modes_x; ++mu_x) {
            // λ_μ,ν = (μπ/L1)^2 + (νπ/L2)^2
            T kx = mu_x * M_PI / l1;
            T ky = mu_y * M_PI / l2;
            lambda_mu_2d(mu_y-1, mu_x-1) = kx * kx + ky * ky;
        }
    }
    
    return lambda_mu_2d;
}

template <typename T>
void select_modes_and_eigenvalues(
    const Eigen::MatrixX<T>& lambda_mu_2d,
    int n_modes,
    Eigen::VectorX<T>& lambda_mu,
    Eigen::MatrixX<int>& selected_indices
) {
    // Get dimensions
    int n_modes_y = lambda_mu_2d.rows();
    int n_modes_x = lambda_mu_2d.cols();
    
    // Create a vector with all eigenvalues
    Eigen::VectorX<T> lambda_mu_flat(n_modes_x * n_modes_y);
    
    // Fill the flat vector in row-major order
    for (int i = 0; i < n_modes_y; ++i) {
        for (int j = 0; j < n_modes_x; ++j) {
            // Convert 2D coordinates to flat index
            int flat_idx = i * n_modes_x + j;
            lambda_mu_flat(flat_idx) = lambda_mu_2d(i, j);
        }
    }
    
    // Sort indices
    std::vector<size_t> indices(lambda_mu_flat.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
        [&lambda_mu_flat](size_t i1, size_t i2) {
            return lambda_mu_flat(i1) < lambda_mu_flat(i2);
        }
    );
    
    // Take only the first n_modes
    indices.resize(n_modes);
    
    // Create selected_indices matrix and sorted lambda_mu vector
    lambda_mu.resize(n_modes);
    selected_indices.resize(n_modes, 2);
    
    for (int i = 0; i < n_modes; ++i) {
        int flat_idx = indices[i];
        
        // Convert flat index back to 2D coordinates
        // Since we used row-major order: flat_idx = row * n_cols + col
        int ky_idx = flat_idx / n_modes_x;  // row index (y)
        int kx_idx = flat_idx % n_modes_x;  // column index (x)
        
        // Store sorted eigenvalue
        lambda_mu(i) = lambda_mu_flat(flat_idx);
        
        // Store mode indices (1-indexed as in the Python code)
        selected_indices(i, 0) = kx_idx + 1;  // x index (μ)
        selected_indices(i, 1) = ky_idx + 1;  // y index (ν)
    }
}