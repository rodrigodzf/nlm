#include "c74_min.h"
#include "c74_min_api.h"
#include "c74_min_attribute.h"
#include "c74_min_operator_sample.h"
#include "c74_min_queue.h"
#include <Eigen/Dense>
#include <mutex>
#include <vk_utils/TimeIntegrators.h>
#include <vk_utils/Parameters.h>
#include <vk_utils/FTM.h>
#include <vk_utils/LinearInterpolator.h>
#include <vk_utils/ParallelFilter.h>
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
        m_n_modes = static_cast<int>(args[0]);
    }

    inlet<>  input    { this, "(signal) Input to the modal resonator", "signal" };
    inlet<>  right    { this, "(dictionary) Dictionary with the parameters", "dictionary" };
    outlet<> output   { this, "(signal) Output from the modal resonator", "signal" };

    argument<number> nmodes { this, "nmodes", "Number of modes", true,
        MIN_ARGUMENT_FUNCTION {
            m_n_modes = static_cast<int>(arg);
        }
    };

    queue<> update_queue {this, MIN_FUNCTION {
        update_all_parameters();
        return {};
    }};

    // attribute<number> n_modes { this, "n_modes", 32,
    //     description {"Number of modes"},
    //     range { 1, 1000 },
    //     setter { MIN_FUNCTION {
    //         update_queue.set();
    //         return { args[0] };
    //     }}
    // };

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

    attribute<number, threadsafe::no, limit::clamp> lx { this, "lx", 0.65,
        description {"Length in meters"},
        range { 0.001, 5.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> cross_sectional_area { this, "cross_sectional_area", 0.5188,
        description {"Cross sectional area in m^2"},
        range { 0.00001, 0.01 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> density { this, "density", 1140.0,
        description {"Mass density in kg/m^3"},
        range { 500.0, 80000.0 },
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

    attribute<number, threadsafe::no, limit::clamp> frequency_independent_loss { this, "findependent_loss", 8e-5,
        description {"Frequency independent loss in kg/(ms)"},
        range { 1e-7, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_dependent_loss { this, "fdependent_loss", 3.4e-5,
        description {"Frequency dependent loss in kg m/s"},
        range { 1e-6, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> moment_of_inertia { this, "moment_of_inertia", 0.171,
        description {"Moment of inertia in mm^4"},
        range { 0.1, 6.0 },
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

    message<> dictionary { this, "dictionary", "dictionary with the parameters",
        MIN_FUNCTION {
            try {
                dict d = {args[0]};

                // set the parameters from the dictionary
                lx = d["lx"];
                cross_sectional_area = d["cross_sectional_area"];
                density = d["density"];
                youngs_modulus = d["youngs_modulus"];
                frequency_independent_loss = d["findependent_loss"];
                frequency_dependent_loss = d["fdependent_loss"];
                tension = d["tension"];
            }
            catch (std::exception& e) {
                cerr << e.what() << endl;
            }
            return {};
        }
    };
    
    message<> maxclass_setup { this, "maxclass_setup",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            cout << "modal.string~ - " << VERSION << " - 2025 - Rodrigo Diaz" << endl;
            return {};
        }
    };

    message<> dspsetup {this, "dspsetup",
        MIN_FUNCTION {
            m_dt = 1.0 / samplerate();
            m_initialized = false;
            update_queue.set(); // Force recalculation of all parameters
            return {};
        }
    };

    message<> reset_coefficients { this, "reset_coefficients",
        MIN_FUNCTION {
            std::unique_lock<std::mutex> lock {m_coeff_mutex};
            m_parallel_filter.reset();
            cout << "Reset parallel filter state" << endl;
            return {};
        }
    };

    sample operator()(sample input) {
        if (!m_initialized) {
            return 0.0;
        }

        std::unique_lock<std::mutex> lock {m_coeff_mutex};

        // Calculate force input
        m_force_input = input * m_force_weights;

        // Update nonlinearity
        m_parallel_filter.update_nonlinearity(m_lambda_mu, m_string_tau_with_norms);
        
        // Process input
        m_parallel_filter(m_force_input);

        // Get the interpolated readout weights and use them for output
        m_readout_weights_lerp.process(m_current_readout_weights);
        sample out = m_current_readout_weights.dot(m_parallel_filter.get_q());

        lock.unlock();
        return out;
    }

private:
    void update_all_parameters() {
        if (!m_initialized || m_prev_n_modes != m_n_modes) {
            m_prev_n_modes = m_n_modes;
            m_readout_weights_lerp.resize(m_n_modes, 441);
            m_current_readout_weights.resize(m_n_modes);
            m_parallel_filter.resize(m_n_modes);
            m_force_input.resize(m_n_modes);
        }

        std::unique_lock<std::mutex> lock {m_coeff_mutex};
        
        // Update string parameters
        m_string_parameters.length = lx; // meters
        m_string_parameters.A = cross_sectional_area * 1e-6; // mm^2 to m^2
        m_string_parameters.E = youngs_modulus * 1e9; // GPa to Pa
        m_string_parameters.I = moment_of_inertia * 1e-12; // mm^4 to m^4
        m_string_parameters.rho = density;
        m_string_parameters.d1 = frequency_independent_loss;
        m_string_parameters.d3 = frequency_dependent_loss;
        m_string_parameters.Ts0 = tension;

        // Calculate modal coefficients
        // int n_mode_long = static_cast<int>(n_modes);
        m_selected_indices = Eigen::VectorXd::LinSpaced(m_n_modes, 1, m_n_modes);
        m_lambda_mu = ftm::string_eigenvalues(m_n_modes, m_string_parameters.length);
        
        m_omega_mu_squared = stiffness_term<double>(m_string_parameters, m_lambda_mu);
        m_gamma2_mu = damping_term<double>(m_string_parameters, m_lambda_mu);
        
        double string_norm = m_string_parameters.length * 0.5;
        double string_tau = (m_string_parameters.E * m_string_parameters.A) / 
                           (m_string_parameters.length * 2.0);
        string_tau = string_tau / m_string_parameters.density() / string_norm;
        m_string_tau_with_norms = string_tau * m_lambda_mu;
        
        // Set coefficients for parallel filter
        m_parallel_filter.set_coefficients(m_gamma2_mu, m_omega_mu_squared, m_dt);
        
        // Calculate force weights
        m_force_weights = ftm::evaluate_string_eigenfunctions(
            m_selected_indices,
            force_position * m_string_parameters.length,
            m_string_parameters.length
        ) / (string_norm * m_string_parameters.density());
        
        // Calculate readout weights
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
        
        m_initialized = true;
        lock.unlock();
    }
    
    double m_dt = 1.0 / 44100.0;
    bool m_initialized = false;
    int m_prev_n_modes = 0;
    int m_n_modes = 32;
    
    std::mutex m_coeff_mutex;
    
    StringParameters m_string_parameters;
    Vector m_selected_indices;
    Vector m_string_tau_with_norms;
    Vector m_lambda_mu;
    Vector m_gamma2_mu;
    Vector m_omega_mu_squared;
    
    Vector m_force_weights;
    Vector m_readout_weights;
    Vector m_current_readout_weights;
    Vector m_force_input;
    
    VectorInterpolator<double> m_readout_weights_lerp;
    ParallelFilter<double> m_parallel_filter;
};

MIN_EXTERNAL(string_tilde);