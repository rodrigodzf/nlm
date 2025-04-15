#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

namespace excitations {

/**
 * Create a 1D raised cosine excitation with time parameters in seconds.
 *
 * @param duration Total duration of the excitation (in seconds)
 * @param start_time Start time of the excitation (in seconds)
 * @param end_time End time of the excitation (in seconds)
 * @param amplitude Amplitude of the excitation
 * @param sample_rate Sample rate (samples per second)
 * @return Vector containing the excitation signal
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> create_1d_raised_cosine(
    T duration,
    T start_time,
    T end_time,
    T amplitude,
    T sample_rate
) {
    int num_samples = static_cast<int>(duration * sample_rate);
    Eigen::Matrix<T, Eigen::Dynamic, 1> excitation = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(num_samples);

    int start_index = static_cast<int>(start_time * sample_rate);
    int end_index = static_cast<int>(end_time * sample_rate);

    if (start_index < 0 || end_index > num_samples || start_index >= end_index) {
        throw std::invalid_argument("Invalid start_time or end_time range.");
    }

    for (int i = start_index; i < end_index; ++i) {
        T t = static_cast<T>(i - start_index) / (end_index - start_index);
        excitation(i) = amplitude * (1.0 - std::cos(2.0 * M_PI * t)) / 2.0;
    }

    return excitation;
}

/**
 * Create a raised cosine function on a 2D grid.
 *
 * @param Nx Number of grid points in the x-direction
 * @param Ny Number of grid points in the y-direction
 * @param h Grid spacing
 * @param center Center of the raised cosine (x, y)
 * @param epsilon Scaling parameter (unused in this implementation but kept for API consistency)
 * @param width Width of the cosine
 * @return Tuple containing the raised cosine matrix and additional computed values
 */
template <typename T>
struct RaisedCosineResult {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> raised_cosine;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> X;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Y;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dist;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dist_x;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dist_y;
};

template <typename T>
RaisedCosineResult<T> create_raised_cosine(
    int Nx, 
    int Ny, 
    T h, 
    const Eigen::Matrix<T, 2, 1>& center, 
    T epsilon,
    T width
) {
    RaisedCosineResult<T> result;
    
    // Create the grid
    result.X = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(Ny+1, Nx+1);
    result.Y = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(Ny+1, Nx+1);
    
    for (int i = 0; i <= Ny; ++i) {
        for (int j = 0; j <= Nx; ++j) {
            result.X(i, j) = j * h;
            result.Y(i, j) = i * h;
        }
    }
    
    // Compute the distance
    result.dist_x = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(Ny+1, Nx+1);
    result.dist_y = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(Ny+1, Nx+1);
    result.dist = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(Ny+1, Nx+1);
    
    for (int i = 0; i <= Ny; ++i) {
        for (int j = 0; j <= Nx; ++j) {
            result.dist_x(i, j) = std::pow(result.X(i, j) - center(0), 2);
            result.dist_y(i, j) = std::pow(result.Y(i, j) - center(1), 2);
            result.dist(i, j) = std::sqrt(result.dist_x(i, j) + result.dist_y(i, j));
        }
    }
    
    // Compute the indicator function and raised cosine
    result.raised_cosine = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(Ny+1, Nx+1);
    
    for (int i = 0; i <= Ny; ++i) {
        for (int j = 0; j <= Nx; ++j) {
            T indicator = (result.dist(i, j) <= width / 2) ? 1.0 : 0.0;
            result.raised_cosine(i, j) = 0.5 * indicator * (1.0 + std::cos(2.0 * M_PI * result.dist(i, j) / width));
        }
    }
    
    return result;
}

/**
 * Create a pluck excitation for a string with a given length and pluck position.
 * The pluck is modeled in the modal domain.
 *
 * @param lambdas Eigenvalues of the Laplacian operator
 * @param pluck_position Position of pluck on the string in meters
 * @param initial_deflection Initial deflection of the string in meters
 * @param string_length Total length of the string in meters
 * @return Vector containing pluck excitation in the modal domain
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> create_pluck_modal(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& lambdas,
    T pluck_position = 0.28,
    T initial_deflection = 0.03,
    T string_length = 1.0
) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> lambdas_sqrt = lambdas.array().sqrt();
    
    // Scaling factor for the initial deflection
    T deflection_scaling = initial_deflection * (string_length / (string_length - pluck_position));
    
    // Compute the coefficients
    Eigen::Matrix<T, Eigen::Dynamic, 1> coefficients = 
        deflection_scaling * 
        (lambdas_sqrt * pluck_position).array().sin() / 
        (lambdas_sqrt * pluck_position).array();
    
    coefficients = coefficients.array() / lambdas_sqrt.array();
    
    return coefficients;
}

} // namespace excitations
