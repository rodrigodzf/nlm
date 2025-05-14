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
#include <vk_utils/ParallelFilter.h>
#include "version.h"

using namespace c74::min;

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

long simplemc_multichanneloutputs(c74::max::t_object* x, long index, long count);

class mcs_string_tilde : public object<mcs_string_tilde>, public vector_operator<> {
public:
    MIN_DESCRIPTION	{ "A modal string." };
    MIN_TAGS		{ "audio" };
    MIN_AUTHOR		{ "Rodrigo Diaz" };
    MIN_RELATED		{ "biquad~" };

    mcs_string_tilde(const atoms& args = {}) {
        if (args.size() > 1) {
            m_num_inputs = static_cast<int>(args[0]);
            m_num_outputs = static_cast<int>(args[1]);
        }

        m_force_positions = Vector::Constant(m_num_inputs, 0.5);
        m_readout_positions = Vector::Constant(m_num_outputs, 0.5);
    }

    inlet<>  input    { this, "(multichannelsignal) Input to the modal resonator" };
    outlet<> output   { this, "(multichannelsignal) Output from the modal resonator", "multichannelsignal" };

    // Single update queue for all parameter changes
    queue<> update_queue {this, MIN_FUNCTION {
        update_all_parameters();
        return {};
    }};

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
            cout << "mcs.modal.string~ - " << VERSION << " - 2025 - Rodrigo Diaz" << endl;
            c74::max::t_class* c = args[0];
            c74::max::class_addmethod(c, (c74::max::method)simplemc_multichanneloutputs, "multichanneloutputs", c74::max::A_CANT, 0);
            return {};
        }
    };

    message<> force_position { this, "force_position",
        MIN_FUNCTION {
            if (args.size() != m_num_inputs) {
                cerr << "force_position: must provide " << m_num_inputs << " positions" << endl;
                return {};
            }
            
            for (int i = 0; i < m_num_inputs; ++i) {
                // Clamp values between 0 and 1
                m_force_positions[i] = std::clamp(double(args[i]), 0.0, 1.0);
            }
            
            update_queue.set();
            return {};
        }
    };

    message<> readout_position { this, "readout_position",
        MIN_FUNCTION {
            if (args.size() != m_num_outputs) {
                cerr << "readout_position: must provide " << m_num_outputs << " positions" << endl;
                return {};
            }

            for (int i = 0; i < m_num_outputs; ++i) {
                // Clamp values between 0 and 1
                m_readout_positions[i] = std::clamp(double(args[i]), 0.0, 1.0);
            }

            update_queue.set();
            return {};
        }
    };

    void operator()(audio_bundle input, audio_bundle output) {
        if (!m_initialized) {
            for (int ch = 0; ch < output.channel_count(); ++ch)
                std::fill(output.samples(ch), output.samples(ch) + output.frame_count(), 0.0);
            return;
        }

        std::unique_lock<std::mutex> lock {m_coeff_mutex};

        for (int frame = 0; frame < input.frame_count(); ++frame) {
            // For each input channel, calculate force input
            for (int in_ch = 0; in_ch < m_num_inputs; ++in_ch) {
                m_force_input.col(in_ch) = input.samples(in_ch)[frame] * m_force_weights.col(in_ch);
            }

            // Update nonlinearity
            m_parallel_filter.update_nonlinearity(m_lambda_mu, m_string_tau_with_norms);
            
            // Process all inputs
            m_parallel_filter(m_force_input.rowwise().sum());

            // Interpolated readout weights
            bool success = m_readout_weights_lerp.process(m_current_readout_weights);
            if (!success) {
                cout << "Failed to process m_readout_weights_lerp" << endl;
            }

            // Output for each output channel
            for (int out_ch = 0; out_ch < m_num_outputs; ++out_ch) {
                output.samples(out_ch)[frame] = m_current_readout_weights.col(out_ch).dot(m_parallel_filter.get_q());
            }
        }
        lock.unlock();
    }

public:
    int m_num_inputs = 1;
    int m_num_outputs = 2;

private:
    void update_all_parameters() {
        if (!m_initialized || m_prev_n_modes != n_modes) {
            m_prev_n_modes = n_modes;
            m_readout_weights_lerp.resize(n_modes, m_num_outputs, 441);
            m_current_readout_weights.resize(n_modes, m_num_outputs);
            m_parallel_filter.resize(n_modes);
            m_force_input.resize(n_modes, m_num_inputs);
        }

        std::unique_lock<std::mutex> lock {m_coeff_mutex};
        
        // Update string parameters
        m_string_parameters.length = lx;
        m_string_parameters.A = cross_sectional_area;
        m_string_parameters.rho = density;
        m_string_parameters.E = youngs_modulus * 1e9;
        m_string_parameters.d1 = frequency_independent_loss;
        m_string_parameters.d3 = frequency_dependent_loss;
        m_string_parameters.Ts0 = tension;
        
        // Calculate modal coefficients
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
        
        // Set coefficients for parallel filter
        m_parallel_filter.set_coefficients(m_gamma2_mu, m_omega_mu_squared, m_dt);
        
        // Calculate force weights matrix
        m_force_weights.resize(n_modes, m_num_inputs);
        for (int i = 0; i < m_num_inputs; ++i) {
            m_force_weights.col(i) = ftm::evaluate_string_eigenfunctions(
                m_selected_indices,
                m_force_positions[i] * m_string_parameters.length,
                m_string_parameters.length
            ) / (string_norm * m_string_parameters.density());
        }
        
        // Calculate readout weights matrix
        Matrix new_readout_weights(static_cast<int>(n_modes), m_num_outputs);
        for (int i = 0; i < m_num_outputs; ++i) {
            new_readout_weights.col(i) = ftm::evaluate_string_eigenfunctions(
                m_selected_indices,
                m_readout_positions[i] * m_string_parameters.length,
                m_string_parameters.length
            );
        }
        
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
    
    std::mutex m_coeff_mutex;
    
    StringParameters m_string_parameters;
    Vector m_selected_indices;
    Vector m_string_tau_with_norms;
    Vector m_lambda_mu;
    Vector m_gamma2_mu;
    Vector m_omega_mu_squared;
    
    Matrix m_force_weights;
    Matrix m_readout_weights;
    Matrix m_current_readout_weights;
    Matrix m_force_input;
    
    Vector m_force_positions;
    Vector m_readout_positions;
    
    MatrixInterpolator<double> m_readout_weights_lerp;
    ParallelFilter<double> m_parallel_filter;
};

long simplemc_multichanneloutputs(c74::max::t_object* x, long index, long count) {
    minwrap<mcs_string_tilde>* ob = (minwrap<mcs_string_tilde>*)(x);
    return ob->m_min_object.m_num_outputs;
}

MIN_EXTERNAL(mcs_string_tilde);