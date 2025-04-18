#include "c74_min.h"
#include "c74_min_api.h"
#include "c74_min_attribute.h"
#include "c74_min_operator_sample.h"
#include "c74_min_queue.h"
#include "matioCpp/MultiDimensionalArray.h"
#include <cmath>
#include <Eigen/Dense>
#include <matioCpp/matioCpp.h>
#include <vk_utils/TimeIntegrators.h>
#include <vk_utils/Parameters.h>
#include <vk_utils/FTM.h>
#include <mutex>
#include <atomic>
using namespace c74::min;

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
class modal_resonator : public object<modal_resonator>, public sample_operator<1, 1> {
public:
    MIN_DESCRIPTION	{ "Modal resonator." };
    MIN_TAGS		{ "audio" };
    MIN_AUTHOR		{ "Rodrigo Diaz" };
    MIN_RELATED		{ "modal_resonator~" };

    // Constructor
    modal_resonator(const atoms& args = {}) {
        // Initialize any variables or parameters here
        // calculate_coefficients();

        m_lambda_mu = Vector::Zero(m_n_phi);
        m_gamma2_mu = Vector::Zero(m_n_phi);
        reset_state();

    }

    queue<> m_queue {this, MIN_FUNCTION {
        calculate_coefficients();
        return {};
    }};

    queue<> update_modal_weights {this, MIN_FUNCTION {
        calculate_modal_weights();
        return {};
    }};

    // Inlets and outlets
    inlet<>  input    { this, "(signal) Input to the modal resonator", "signal" };
    outlet<> output   { this, "(signal) Output from the modal resonator", "signal" };

    // Attributes
    attribute<number> sample_rate { this, "sample_rate", 44100.0,
        description {"Sample rate in Hz"},
        range { 16000.0, 192000.0 },
        setter { MIN_FUNCTION {
            cout << "Setting sample rate to " << args[0] << endl;
            m_sampleRate = args[0];
            m_dt = 1.0 / m_sampleRate;
            return {
                m_sampleRate
            };
        } }
    };

    attribute<numbers> force_position { this, "force_position", {{0.5, 0.5}},
        description {"Force position (0-1)"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting force position to " << args[0] << " " << args[1] << endl;
            update_modal_weights.set();
            return {
                args
            };
        } }
    };

    attribute<numbers> readout_position { this, "readout_position", {{0.5, 0.5}},
        description {"Readout position (0-1)"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting readout position to " << args[0] << " " << args[1] << endl;
            update_modal_weights.set();
            return {
                args
            };
        } }
    };

    attribute<number> poisson_ratio { this, "poisson_ratio", 0.3,
        description {"Poisson's ratio"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting poisson_ratio to " << args[0] << endl;
            // m_plate_parameters.nu = args[0];
            m_queue.set();
            return {
                args[0]
            };
        } }
    };

    attribute<number> thickness { this, "thickness", 5e-4,
        description {"Thickness in meters"},
        range { 1e-4, 1e-3 },
        setter { MIN_FUNCTION {
            cout << "Setting thickness to " << args[0] << endl;
            // m_plate_parameters.h = args[0];
            m_queue.set();
            return {
                args[0]
            };
        } }
    };
    
    attribute<number> lx { this, "lx", 0.2,
        description {"Length in meters"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting lx to " << args[0] << endl;
            // m_plate_parameters.l1 = args[0];
            return {
                args[0]
            };
        } }
    };

    attribute<number> ly { this, "ly", 0.3,
        description {"Width in meters"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting ly to " << args[0] << endl;
            // m_plate_parameters.l2 = args[0];
            return {
                args[0]
            };
        } }
    };

    attribute<number> density { this, "density", 7850,
        description {"Density in kg/m^3"},
        range { 1000.0, 10000.0 },
        setter { MIN_FUNCTION {
            cout << "Setting density to " << args[0] << endl;
            // m_plate_parameters.rho = args[0];
            m_queue.set();
            return {
                args[0]
            };
        } }
    };

    attribute<number> youngs_modulus { this, "youngs_modulus", 2000,
        description {"Young's modulus in GPa"},
        range { 100.0, 10000.0 },
        setter { MIN_FUNCTION {
            cout << "Setting youngs_modulus to " << args[0] << endl;
            // m_plate_parameters.E = static_cast<double>(args[0]) * 1e9;
            m_queue.set();
            return {
                args[0]
            };
        } }
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_independent_loss { this, "frequency_independent_loss", 0.01,
        description {"Frequency independent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting frequency_independent_loss to " << args[0] << endl;
            // m_plate_parameters.d1 = static_cast<double>(args[0]) * 1e-3;
            m_queue.set();
            return {
                args[0]
            };
        } }
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_dependent_loss { this, "frequency_dependent_loss", 0.01,
        description {"Frequency dependent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting frequency_dependent_loss to " << args[0] << endl;
            // m_plate_parameters.d3 = static_cast<double>(args[0]) * 1e-3;
            m_queue.set();
            return {
                args[0]
            };
        } }
    };

    attribute<number> surface_tension { this, "surface_tension", 0.0,
        description {"Surface tension in N/m"},
        range { 0.0, 1000.0 },
        setter { MIN_FUNCTION {
            cout << "Setting surface tension to " << args[0] << endl;
            // m_plate_parameters.Ts0 = args[0];
            m_queue.set();
            return {
                args[0]
            };
        } }
    };  
    

    // Methods using direct binding approach
    // post to max window == but only when the class is loaded the first time
    message<> maxclass_setup { this, "maxclass_setup",
        MIN_FUNCTION {
            cout << "modal_resonator class loaded" << endl;
            return {};
        }
    };

    message<> set_couplings_and_eigenvalues { this, "set_couplings_and_eigenvalues",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            
            // Load data from MATLAB file
            std::string filename = args[0];
            matioCpp::File file(filename);

            cout << "Loading data from " << filename << endl;
            
            matioCpp::MultiDimensionalArray<double> H_mat = file.read("H").asMultiDimensionalArray<double>();
            matioCpp::Vector<double> lambda_mu_mat = file.read("lambda_mu").asVector<double>();
            matioCpp::Vector<double> selected_indices_x = file.read("selected_indices_x").asVector<double>();
            matioCpp::Vector<double> selected_indices_y = file.read("selected_indices_y").asVector<double>();

            m_lambda_mu = matioCpp::to_eigen(lambda_mu_mat);

            // get the size of lambda_mu
            int size_lambda_mu = m_lambda_mu.size();
            cout << "Size of lambda_mu: " << size_lambda_mu << endl;

            // convert to Eigen::Matrix
            m_H = Eigen::Map<Matrix>(H_mat.data(), m_n_psi * m_n_phi, m_n_phi);
            m_selected_indices_x = matioCpp::to_eigen(selected_indices_x);
            m_selected_indices_y = matioCpp::to_eigen(selected_indices_y);

            m_initialized = true;
            calculate_coefficients();
            calculate_modal_weights();
            return {};
        }
    };
    

    // DSP perform method (for Max/MSP audio processing)
    sample operator()(sample input) {
        if (!m_initialized) {
            return 0.0;
        }

        std::unique_lock<std::mutex> lock {m_coeff_mutex};

        // Each input is a scalar which weights the modal gains at the excitation position
        auto x = input * m_force_weights;

        Vector t0_flat = m_H_scaled * m_q;
        Matrix t0 = Eigen::Map<Matrix>(t0_flat.data(), m_q.size(), m_q.size());
        Vector t2 = t0 * m_q;
        Vector nl = t0.transpose() * t2;

        Vector q_next = m_B.cwiseProduct(m_q) + m_C.cwiseProduct(m_q_prev) + m_A_inv.cwiseProduct(x - nl);
        
        // Update previous state
        m_q_prev = m_q;
        m_q = q_next;

        // Get a scalar output at the readout position
        sample out = m_readout_weights.dot(q_next);

        lock.unlock();
        return out;
    }


private:

    // Private methods
    void setup(const atoms& args, const int vector_size) {
        // Get current sample rate from Max
        if (!args.empty())
            m_sampleRate = args[0];
    }

    void reset_state() {
        m_q = Vector::Zero(m_n_phi);
        m_q_prev = Vector::Zero(m_n_phi);
    }
    


    void calculate_coefficients() {
        if (!m_initialized) {
            return;
        }

        m_plate_parameters.nu = poisson_ratio;
        m_plate_parameters.h = thickness;
        m_plate_parameters.l1 = lx;
        m_plate_parameters.l2 = ly;
        m_plate_parameters.rho = density;
        m_plate_parameters.E = youngs_modulus * 1e9;
        m_plate_parameters.d1 = frequency_independent_loss;
        m_plate_parameters.d3 = frequency_dependent_loss;
        m_plate_parameters.Ts0 = surface_tension;

        
        m_omega_mu_squared = stiffness_term<double>(m_plate_parameters, m_lambda_mu);
        m_gamma2_mu = damping_term<double>(m_plate_parameters, m_lambda_mu);

        // These lines will cause a race condition and crash the program
        double plate_norm = m_plate_parameters.l1 * m_plate_parameters.l2 * 0.25;
        double scale = (m_plate_parameters.E * plate_norm) / (2.0 * m_plate_parameters.density());

        // std::lock_guard<std::mutex> lock(m_coeff_mutex);
        std::unique_lock<std::mutex> lock {m_coeff_mutex};

        // m_H *= std::sqrt(scale);
        m_H_scaled = m_H * std::sqrt(scale);
        calculate_coefficients_tf(
            m_gamma2_mu,
            m_omega_mu_squared,
            m_B,
            m_C,
            m_A_inv,
            m_dt
        );
        lock.unlock();

    }

    void calculate_modal_weights() {
        if (!m_initialized) {
            return;
        }


        // phi has shape (ny * nx, n_modes) but we need to index it with ny and nx
        // int force_nx = static_cast<int>(std::round(force_position[0] * m_nx_elem));
        // int force_ny = static_cast<int>(std::round(force_position[1] * m_ny_elem));
        // int force_row_idx = force_ny * m_nx_elem + force_nx;
        // cout << "force_row_idx: " << force_row_idx << endl;

        // int readout_nx = static_cast<int>(std::round(readout_position[0] * m_nx_elem));
        // int readout_ny = static_cast<int>(std::round(readout_position[1] * m_ny_elem));
        // int readout_row_idx = readout_ny * m_nx_elem + readout_nx;
        // cout << "readout_row_idx: " << readout_row_idx << endl;
        // m_force_weights = m_phi.row(force_row_idx);
        // m_readout_weights = m_phi.row(readout_row_idx);
        double plate_norm = m_plate_parameters.l1 * m_plate_parameters.l2 * 0.25;

        std::unique_lock<std::mutex> lock {m_coeff_mutex};

        m_force_weights = ftm::evaluate_rectangular_eigenfunctions(
            m_selected_indices_x, m_selected_indices_y, 
            force_position[0] * m_plate_parameters.l1, force_position[1] * m_plate_parameters.l2, 
            m_plate_parameters.l1, m_plate_parameters.l2
        ) / (plate_norm * m_plate_parameters.density());

        m_readout_weights = ftm::evaluate_rectangular_eigenfunctions(
            m_selected_indices_x, m_selected_indices_y, 
            readout_position[0] * m_plate_parameters.l1, readout_position[1] * m_plate_parameters.l2, 
            m_plate_parameters.l1, m_plate_parameters.l2
        );

        lock.unlock();

    }

    // Private member variables
    double m_sampleRate = 44100.0;
    double m_dt = 1.0 / m_sampleRate;

    // Thread-safe coefficient management

    std::mutex m_coeff_mutex;
    
    bool m_initialized = false;

    int m_n_psi = 10;
    int m_n_phi = 10;
    PlateParameters m_plate_parameters;
    Vector m_selected_indices_x;
    Vector m_selected_indices_y;
    Matrix m_phi;
    Matrix m_H;
    Matrix m_H_scaled;
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

};

MIN_EXTERNAL(modal_resonator);