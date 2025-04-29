#include "c74_min.h"
#include "c74_min_api.h"
#include "c74_min_attribute.h"
#include "c74_min_operator_sample.h"
#include "c74_min_queue.h"
#include <Eigen/Dense>
#include <cmath>
#include <mutex>
#include <vk_utils/TimeIntegrators.h>
#include <vk_utils/Parameters.h>
#include <vk_utils/FTM.h>
#include <vk_utils/LinearInterpolator.h>
#include "version.h"

using namespace c74::min;

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

class string_tilde : public object<string_tilde>, public sample_operator<1, 1> {
public:
    MIN_DESCRIPTION	{ "A modal string." };
    MIN_TAGS		{ "audio" };
    MIN_AUTHOR		{ "Rodrigo Diaz" };
    MIN_RELATED		{ "biquad~" };

    string_tilde(const atoms& args = {}) {

    }

    inlet<>  input    { this, "(signal) Input to the modal resonator", "signal" };
    outlet<> output   { this, "(signal) Output from the modal resonator", "signal" };

    queue<> update_queue {this, MIN_FUNCTION {
        update_all_parameters();
        return {};
    }};

    attribute<number, threadsafe::no, limit::clamp> force_position { this, "force_position", 0.5,
        description {"Force position (0-1)"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> readout_position { this, "readout_position", 0.5,
        description {"Readout position (0-1)"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args };
        }}
    };

    attribute<number> n_modes { this, "n_modes", 32,
        description {"Number of modes"},
        range { 1, 1000 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> lx { this, "lx", 0.75,
        description {"Length in meters"},
        range { 0.001, 5.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> cross_sectional_area { this, "cross_sectional_area", 0.5188e-6,
        description {"Cross sectional area"},
        range { 0.00001, 0.01 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };
    attribute<number> density { this, "density", 1140.0,
        description {"Density"},
        range { 1000.0, 10000.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> youngs_modulus { this, "youngs_modulus", 5.4,
        description {"Young's modulus in GPa"},
        range { 1.0, 1000.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_independent_loss { this, "frequency_independent_loss", 8e-5,
        description {"Frequency independent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_dependent_loss { this, "frequency_dependent_loss", 1.4e-5,
        description {"Frequency dependent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> tension { this, "tension", 60.97,
        description {"Tension in N"},
        range { 0.0, 1000.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };  
    
    message<> maxclass_setup { this, "maxclass_setup",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            cout << "string~ - " << VERSION << " - 2025 - Rodrigo Diaz" << endl;
            return {};
        }
    };

    message<> dspsetup {this, "dspsetup",
        MIN_FUNCTION {
            m_dt = 1.0 / samplerate();
            m_initialized = false;
            reset_state();
            update_queue.set(); // Force recalculation of all parameters
            return {};
        }
    };
    
    sample operator()(sample input) {
        if (!m_initialized) {
            return 0.0;
        }

        std::unique_lock<std::mutex> lock {m_coeff_mutex};

        // Use pre-allocated vector for input force calculation
        m_force_input = input * m_force_weights;

        calculate_nonlinear_berger(m_lambda_mu, m_string_tau_with_norms, m_q, m_nl);

        m_q_next.noalias() = m_B.cwiseProduct(m_q) +
                             m_C.cwiseProduct(m_q_prev) +
                             m_A_inv.cwiseProduct(m_force_input - m_nl);
        
        m_q_prev = m_q;
        m_q = m_q_next;

        // Get the interpolated readout weights and use them for output
        m_readout_weights_lerp.process(m_current_readout_weights);
        sample out = m_current_readout_weights.dot(m_q);

        lock.unlock();
        return out;
    }

private:
    // Reset state vectors
    void reset_state() {
        m_q = Vector::Zero(n_modes);
        m_q_prev = Vector::Zero(n_modes);
        m_q_next = Vector::Zero(n_modes);
        m_nl = Vector::Zero(n_modes);
        m_force_input = Vector::Zero(n_modes);
    }

    // Unified function to update all parameters
    void update_all_parameters() {
        std::unique_lock<std::mutex> lock {m_coeff_mutex};
        
        // Step 1: Reset state if needed
        if (!m_initialized || m_prev_n_modes != n_modes) {
            reset_state();
            m_prev_n_modes = n_modes;
        }
        
        // Step 2: Update string parameters
        m_string_parameters.length = lx;
        m_string_parameters.A = cross_sectional_area;
        m_string_parameters.rho = density;
        m_string_parameters.E = youngs_modulus * 1e9;
        m_string_parameters.d1 = frequency_independent_loss;
        m_string_parameters.d3 = frequency_dependent_loss;
        m_string_parameters.Ts0 = tension;
        
        // Step 3: Calculate modal coefficients
        int n_mode_long = static_cast<int>(n_modes);
        m_selected_indices = Eigen::VectorXd::LinSpaced(n_mode_long, 1, n_mode_long);
        m_lambda_mu = ftm::string_eigenvalues(n_mode_long, m_string_parameters.length);
        
        m_omega_mu_squared = stiffness_term<double>(m_string_parameters, m_lambda_mu);
        m_gamma2_mu = damping_term<double>(m_string_parameters, m_lambda_mu);
        
        double string_norm = m_string_parameters.length * 0.5;
        double string_tau = (m_string_parameters.E * m_string_parameters.A) / 
                           (m_string_parameters.length * 2.0);
        string_tau = string_tau / m_string_parameters.density() / string_norm;
        m_string_tau_with_norms = string_tau * m_lambda_mu;
        
        calculate_coefficients_tf(
            m_gamma2_mu,
            m_omega_mu_squared,
            m_B,
            m_C,
            m_A_inv,
            m_dt
        );
        
        // Step 4: Calculate force weights
        m_force_weights = ftm::evaluate_string_eigenfunctions(
            m_selected_indices,
            force_position * m_string_parameters.length,
            m_string_parameters.length
        ) / (string_norm * m_string_parameters.density());
        
        // Step 5: Calculate readout weights
        Vector new_readout_weights = ftm::evaluate_string_eigenfunctions(
            m_selected_indices,
            readout_position * m_string_parameters.length,
            m_string_parameters.length
        );
        
        if (m_readout_weights_lerp.isFinished()) {
            if (m_readout_weights.size() == 0) {
                m_readout_weights = new_readout_weights;
                m_readout_weights_lerp.setValue(new_readout_weights);
            }
            m_readout_weights_lerp.setTarget(new_readout_weights);
        }
        
        // Step 6: Mark as initialized
        m_initialized = true;
        
        lock.unlock();
    }
    
    // Track previous mode count to detect changes
    int m_prev_n_modes = 0;
    
    double m_dt = 1.0 / 44100.0;
    double m_interpolation_time_ms = 10.0;
    bool m_initialized = false;
    
    std::mutex m_coeff_mutex;
    
    StringParameters m_string_parameters;
    Vector m_selected_indices;
    
    Vector m_string_tau_with_norms;

    Vector m_lambda_mu;
    Vector m_gamma2_mu;
    Vector m_omega_mu_squared;
    
    Vector m_A_inv;
    Vector m_B;
    Vector m_C;
    
    Vector m_q;
    Vector m_q_prev;
    
    Vector m_force_weights;
    Vector m_readout_weights;
    
    // Add VectorLerp for interpolating readout weights
    VectorInterpolator<double> m_readout_weights_lerp;

    // Temporary vectors used in sample processing to avoid stack allocations
    Vector m_t0_flat;
    Vector m_t2;
    Vector m_nl;
    Vector m_q_next;
    Vector m_current_readout_weights;
    Vector m_force_input;  // Vector for storing input * m_force_weights
};

MIN_EXTERNAL(string_tilde);