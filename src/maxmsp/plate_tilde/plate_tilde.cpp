#include "c74_min.h"
#include "c74_min_api.h"
#include "c74_min_attribute.h"
#include "c74_min_operator_sample.h"
#include "c74_min_queue.h"
#include <Eigen/Dense>
#include <matioCpp/matioCpp.h>
#include <matioCpp/MultiDimensionalArray.h>
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

class plate_tilde : public object<plate_tilde>, public sample_operator<1, 1> {
public:
    MIN_DESCRIPTION	{ "A modal plate." };
    MIN_TAGS		{ "audio" };
    MIN_AUTHOR		{ "Rodrigo Diaz" };
    MIN_RELATED		{ "biquad~" };

    plate_tilde(const atoms& args = {}) {
        m_lambda_mu = Vector::Zero(m_n_phi);
        m_gamma2_mu = Vector::Zero(m_n_phi);
        reset_state();
        
        // Initialize model type
        if (args.size() > 0) {
            if (args[0] == "vk") {
                m_model_type = "vk";
            } else {
                m_model_type = "berger";
            }
        } else {
            m_model_type = "berger";  // Default model type
        }
        
        // Update the attribute
        model_type = m_model_type;

        cout << "Model type: " << m_model_type << endl;
    }

    inlet<>  input    { this, "(signal) Input to the modal resonator", "signal" };
    outlet<> output   { this, "(signal) Output from the modal resonator", "signal" };

    attribute<symbol> model_type { this, "model_type", 
        title {"Model Type"},
        description {"Plate model type (vk or berger)"},
        getter { MIN_GETTER_FUNCTION { return { m_model_type }; } },
        readonly { true }
    };

    // Single update queue for all parameter changes
    queue<> update_queue {this, MIN_FUNCTION {
        update_all_parameters();
        return {};
    }};

    attribute<numbers> force_position { this, "force_position", {{0.5, 0.5}},
        description {"Force position (0-1)"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args };
        }}
    };

    attribute<numbers> readout_position { this, "readout_position", {{0.5, 0.5}},
        description {"Readout position (0-1)"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> poisson_ratio { this, "poisson_ratio", 0.3,
        description {"Poisson's ratio"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> thickness { this, "thickness", 5e-4,
        description {"Thickness in meters"},
        range { 1e-4, 1e-2 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };
    
    attribute<number, threadsafe::no, limit::clamp> lx { this, "lx", 0.2,
        description {"Length in meters"},
        range { 0.001, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> ly { this, "ly", 0.3,
        description {"Width in meters"},
        range { 0.001, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> density { this, "density", 7850,
        description {"Density in kg/m^3"},
        range { 1000.0, 10000.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> youngs_modulus { this, "youngs_modulus", 2000,
        description {"Young's modulus in GPa"},
        range { 100.0, 10000.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_independent_loss { this, "frequency_independent_loss", 0.01,
        description {"Frequency independent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_dependent_loss { this, "frequency_dependent_loss", 0.01,
        description {"Frequency dependent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> surface_tension { this, "surface_tension", 0.0,
        description {"Surface tension in N/m"},
        range { 0.0, 1000.0 },
        setter { MIN_FUNCTION {
            update_queue.set();
            return { args[0] };
        }}
    };  
    
    message<> maxclass_setup { this, "maxclass_setup",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            cout << "plate~ - " << VERSION << " - 2025 - Rodrigo Diaz" << endl;
            return {};
        }
    };

    message<> dspsetup {this, "dspsetup",
        MIN_FUNCTION {
            m_dt = 1.0 / samplerate();
            // Don't reset initialized flag if data is loaded from file
            if (m_initialized) {
                allocate_temp_vectors(); // Ensure temp vectors are allocated with correct size
                update_queue.set();       // Update parameters but preserve loaded data
            }
            return {};
        }
    };
    message<> set_couplings_and_eigenvalues { this, "set_couplings_and_eigenvalues",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            std::string filename = args[0];
            matioCpp::File file(filename);

            cout << "Loading data from " << filename << endl;
            
            matioCpp::MultiDimensionalArray<double> H_mat = file.read("H").asMultiDimensionalArray<double>();
            matioCpp::Vector<double> lambda_mu_mat = file.read("lambda_mu").asVector<double>();
            matioCpp::Vector<double> selected_indices_x = file.read("selected_indices_x").asVector<double>();
            matioCpp::Vector<double> selected_indices_y = file.read("selected_indices_y").asVector<double>();

            m_lambda_mu = matioCpp::to_eigen(lambda_mu_mat);
            int size_lambda_mu = m_lambda_mu.size();
            cout << "Size of lambda_mu: " << size_lambda_mu << endl;

            // get the dimensions of H_mat
            m_n_psi = H_mat.dimensions()[0];
            m_n_phi = H_mat.dimensions()[1];
            cout << "n_psi: " << m_n_psi << " n_phi: " << m_n_phi << endl;

            m_H_original = Eigen::Map<Matrix>(H_mat.data(), m_n_psi * m_n_phi, m_n_phi);
            m_selected_indices_x = matioCpp::to_eigen(selected_indices_x);
            m_selected_indices_y = matioCpp::to_eigen(selected_indices_y);

            m_data_loaded = true;
            reset_state();
            allocate_temp_vectors();
            update_queue.set(); // Trigger parameter update after data loading
            return {};
        }
    };
    
    message<> info { this, "info",
        MIN_FUNCTION {
            std::stringstream ss;
            ss << "plate~ parameters:" << '\n'
               << "Model type: " << m_model_type << '\n'
               << "Thickness: " << thickness << " m" << '\n'
               << "Dimensions: " << lx << " x " << ly << " m" << '\n'
               << "Density: " << density << " kg/m³" << '\n'
               << "Young's modulus: " << youngs_modulus << " GPa" << '\n'
               << "Poisson's ratio: " << poisson_ratio << '\n'
               << "Surface tension: " << surface_tension << " N/m" << '\n'
               << "Loss parameters: " << frequency_independent_loss << " (freq. indep.), " 
                 << frequency_dependent_loss << " (freq. dep.)" << '\n'
               << "Force position: (" << force_position[0] << ", " << force_position[1] << ")" << '\n'
               << "Readout position: (" << readout_position[0] << ", " << readout_position[1] << ")" << '\n'
               << "Sampling rate: " << samplerate() << " Hz" << '\n'
               << "Modes: " << m_n_phi;
            cout << ss.str() << endl;
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

        if (m_model_type == "vk") {
            calculate_nonlinear_vk(m_H_scaled, m_q, m_nl, m_n_psi, m_n_phi);
        } else if (m_model_type == "berger") {
            calculate_nonlinear_berger(m_lambda_mu, m_plate_tau_with_norms, m_q, m_nl);
        }

        m_q_next.noalias() = m_B.cwiseProduct(m_q) +
                             m_C.cwiseProduct(m_q_prev) +
                             m_A_inv.cwiseProduct(m_force_input - m_nl);
        
        m_q_prev = m_q;
        m_q = m_q_next;

        // Get the interpolated readout weights and use them for output
        m_current_readout_weights = m_readout_weights_lerp.process();
        sample out = m_current_readout_weights.dot(m_q);

        lock.unlock();
        return out;
    }

private:
    void reset_state() {
        m_q = Vector::Zero(m_n_phi);
        m_q_prev = Vector::Zero(m_n_phi);
        m_nl = Vector::Zero(m_n_phi);
        m_q_next = Vector::Zero(m_n_phi);
    }
    
    void allocate_temp_vectors() {
        // Allocate temporary vectors used in audio processing
        m_t0_flat.resize(m_n_psi * m_n_phi);
        m_t2.resize(m_n_phi);
        m_nl.resize(m_n_phi);
        m_q_next.resize(m_n_phi);
        m_current_readout_weights.resize(m_n_phi);
        m_force_input.resize(m_n_phi);  // Pre-allocate vector for input force
    }
    
    // Unified function to update all parameters
    void update_all_parameters() {
        if (!m_data_loaded) {
            return; // Can't update without data loaded from file
        }

        std::unique_lock<std::mutex> lock {m_coeff_mutex};
        
        // Step 1: Update plate parameters
        m_plate_parameters.nu = poisson_ratio;
        m_plate_parameters.h = thickness;
        m_plate_parameters.l1 = lx;
        m_plate_parameters.l2 = ly;
        m_plate_parameters.rho = density;
        m_plate_parameters.E = youngs_modulus * 1e9;
        m_plate_parameters.d1 = frequency_independent_loss;
        m_plate_parameters.d3 = frequency_dependent_loss;
        m_plate_parameters.Ts0 = surface_tension;
        
        // Step 2: Calculate modal coefficients
        m_omega_mu_squared = stiffness_term<double>(m_plate_parameters, m_lambda_mu);
        m_gamma2_mu = damping_term<double>(m_plate_parameters, m_lambda_mu);
        
        double plate_norm = m_plate_parameters.l1 * m_plate_parameters.l2 * 0.25;
        
        // Step 3: Model-specific calculations
        if (m_model_type == "berger") {
            double plate_tau = (m_plate_parameters.E * m_plate_parameters.h) / (
                2.0 * m_plate_parameters.l1 * m_plate_parameters.l2 * (1.0 - m_plate_parameters.nu * m_plate_parameters.nu)
            );
            plate_tau = plate_tau / m_plate_parameters.density() / plate_norm;
            m_plate_tau_with_norms = plate_tau * m_lambda_mu;
        } else if (m_model_type == "vk") {
            double scale = (m_plate_parameters.E * plate_norm) / (2.0 * m_plate_parameters.density());
            m_H_scaled = m_H_original * std::sqrt(scale);
        }
        
        // Step 4: Calculate time integration coefficients
        calculate_coefficients_tf(
            m_gamma2_mu,
            m_omega_mu_squared,
            m_B,
            m_C,
            m_A_inv,
            m_dt
        );
        
        // Step 5: Calculate force weights
        m_force_weights = ftm::evaluate_rectangular_eigenfunctions(
            m_selected_indices_x, m_selected_indices_y, 
            force_position[0] * m_plate_parameters.l1, force_position[1] * m_plate_parameters.l2, 
            m_plate_parameters.l1, m_plate_parameters.l2
        ) / (plate_norm * m_plate_parameters.density());

        // Step 6: Calculate readout weights
        Vector new_readout_weights = ftm::evaluate_rectangular_eigenfunctions(
            m_selected_indices_x, m_selected_indices_y, 
            readout_position[0] * m_plate_parameters.l1, readout_position[1] * m_plate_parameters.l2, 
            m_plate_parameters.l1, m_plate_parameters.l2
        );
        
        if (m_readout_weights_lerp.isFinished()) {
            if (m_readout_weights.size() == 0) {
                m_readout_weights = new_readout_weights;
                m_readout_weights_lerp.setValue(new_readout_weights);
            }
            m_readout_weights_lerp.setTarget(new_readout_weights);
        }
        
        // Step 7: Mark as initialized
        m_initialized = true;
        
        lock.unlock();
    }
    
    double m_dt = 1.0 / 44100.0;
    double m_interpolation_time_ms = 10.0;
    bool m_initialized = false;
    int m_n_psi = 10;
    int m_n_phi = 10;
    
    std::mutex m_coeff_mutex;
    
    PlateParameters m_plate_parameters;
    Vector m_selected_indices_x;
    Vector m_selected_indices_y;
    
    Matrix m_H_original;
    Matrix m_H_scaled;
    Vector m_plate_tau_with_norms;

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
    
    std::string m_model_type = "berger";  // Default model type

    // Track if data is loaded from file
    bool m_data_loaded = false;
};

MIN_EXTERNAL(plate_tilde);