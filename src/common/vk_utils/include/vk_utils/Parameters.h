#pragma once

#include <cmath>
#include <Eigen/Dense>

// Base class for all parameter types
class PhysicalParameters {
public:
    virtual ~PhysicalParameters() = default;

    // Virtual method for getting density (implemented by derived classes)
    virtual double density() const = 0;

    // Virtual method for getting bending stiffness (implemented by derived classes)
    virtual double bending_stiffness() const = 0;
};

class StringParameters : public PhysicalParameters {
public:
    // Physical parameters
    double A{0.5188e-6};    // m**2    Cross section area
    double I{0.141e-12};    // m**4    Moment of inertia
    double rho{1140};       // kg/m**3 Density
    double E{5.4e9};        // Pa      Young's modulus
    double d1{8e-5};        // kg/(ms) Frequency independent loss
    double d3{1.4e-5};      // kg m/s  Frequency dependent loss
    double Ts0{60.97};      // N       Tension
    double length{0.65};    // m       Length of the string

    // Factory methods for different string types
    static StringParameters piano_string() {
        StringParameters params;
        params.A = 1.54e-6;
        params.I = 4.12e-12;
        params.rho = 57.0e3;
        params.E = 19.5e9;
        params.d1 = 3e-3;
        params.d3 = 2e-5;
        params.Ts0 = 2104;
        params.length = 1.08;
        return params;
    }

    static StringParameters bass_string() {
        StringParameters params;
        params.A = 2.4e-6;
        params.I = 0.916e-12;
        params.rho = 6300;
        params.E = 5e9;
        params.d1 = 6e-3;
        params.d3 = 1e-3;
        params.Ts0 = 114;
        params.length = 1.05;
        return params;
    }

    static StringParameters guitar_string_D() {
        StringParameters params;
        params.A = 7.96e-7;
        params.I = 0.171e-12;
        params.rho = 1140;
        params.E = 5.4e9;
        params.d1 = 8e-5;
        params.d3 = 1.4e-5;
        params.Ts0 = 13.35;
        params.length = 0.65;
        return params;
    }

    static StringParameters guitar_string_B_schafer() {
        StringParameters params;
        params.A = 0.5e-6;
        params.I = 0.17e-12;
        params.rho = 1140;
        params.E = 5.4e9;
        params.d1 = 8e-5;
        params.d3 = 1.4e-5;
        params.Ts0 = 60.97;
        params.length = 0.65;
        return params;
    }

    // Property implementations
    double density() const override {
        return rho * A;
    }

    double bending_stiffness() const override {
        return E * I;
    }
};

class PlateParameters : public PhysicalParameters {
public:
    double h{5e-4};         // m           Thickness
    double l1{0.2};         // m           Width
    double l2{0.3};         // m           Height
    double rho{7.8e3};      // kg/m³       Density
    double E{2e12};         // Pa          Young's modulus
    double nu{0.3};         //             Poisson's ratio
    double d1{4.2e-2};      //             Frequency independent loss
    double d3{2.3e-3};      //             Frequency dependent loss
    double Ts0{100};        // N/m         Surface Tension

    // Property implementations
    double density() const override {
        return rho * h;
    }

    double bending_stiffness() const override {
        return E * std::pow(h, 3) / (12 * (1 - nu * nu));
    }
};

class CircularDrumHeadParameters : public PhysicalParameters {
public:
    double h{1.9e-4};       // m           Thickness
    double r0{0.328};       // m           Radius
    double I{0.57e-12};     // m**4        Moment of inertia
    double rho{1.38e3};     // kg/m**3     Density
    double E{3.5e9};        // Pa          Young's modulus
    double nu{0.35};        //             Poisson's ratio
    double d1{0.14};        // kg/(m**2 s) Frequency independent loss
    double d3{0.32};        // kg/s        Frequency dependent loss
    double Ts0{3990};       // N/m         Surface Tension
    double f0{143.95};      // Hz          Fundamental frequency

    // Property implementations
    double density() const override {
        return rho * h;
    }

    double bending_stiffness() const override {
        return E * std::pow(h, 3) / (12 * (1 - nu * nu));
    }

    static CircularDrumHeadParameters avanzini() {
        CircularDrumHeadParameters params;
        params.h = 2e-4;
        params.r0 = 0.20;
        params.rho = 1350.0;
        params.E = 3.5e9;
        params.nu = 0.2;
        params.d1 = 1.25;
        params.d3 = 5e-4;
        params.Ts0 = 1500;
        return params;
    }
};

/**
 * Calculates the damping term for a physical system.
 * 
 * @param params The physical parameters of the system
 * @param lambda_mu The eigenvalues of the Laplacian operator
 * @return Vector containing damping terms for each mode
 */
template <typename T, typename ParamType>
Eigen::Matrix<T, Eigen::Dynamic, 1> damping_term(
    const ParamType& params,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& lambda_mu
) {
    return (params.d1 + params.d3 * lambda_mu.array()) / params.density();
}

/**
 * Calculates a simplified damping term based on the eigenvalues.
 * 
 * @param lambda_mu The eigenvalues of the Laplacian operator
 * @param factor Scaling factor for the damping term (default: 1e-3)
 * @return Vector containing damping terms for each mode
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> damping_term_simple(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& lambda_mu,
    T factor = 1e-3
) {
    return factor * lambda_mu;
}

/**
 * Calculates the stiffness term for a physical system.
 * 
 * @param params The physical parameters of the system
 * @param lambda_mu The eigenvalues of the Laplacian operator
 * @return Vector containing stiffness terms for each mode
 */
template <typename T, typename ParamType>
Eigen::Matrix<T, Eigen::Dynamic, 1> stiffness_term(
    const ParamType& params,
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& lambda_mu
) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> omega_mu = 
        params.bending_stiffness() * lambda_mu.array().square() + 
        params.Ts0 * lambda_mu.array();
    
    return omega_mu / params.density();
}
