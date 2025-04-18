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
    }

    inlet<>  input    { this, "(signal) Input to the modal resonator", "signal" };
    outlet<> output   { this, "(signal) Output from the modal resonator", "signal" };

    queue<> m_queue {this, MIN_FUNCTION {
        calculate_coefficients();
        return {};
    }};

    queue<> update_modal_weights {this, MIN_FUNCTION {
        calculate_modal_weights();
        return {};
    }};

    attribute<number> sample_rate { this, "sample_rate", 44100.0,
        description {"Sample rate in Hz"},
        range { 16000.0, 192000.0 },
        setter { MIN_FUNCTION {
            // cout << "Setting sample rate to " << args[0] << endl;
            m_sampleRate = args[0];
            m_dt = 1.0 / m_sampleRate;
            return { m_sampleRate };
        }}
    };

    attribute<numbers> force_position { this, "force_position", {{0.5, 0.5}},
        description {"Force position (0-1)"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            update_modal_weights.set();
            return { args };
        }}
    };

    attribute<numbers> readout_position { this, "readout_position", {{0.5, 0.5}},
        description {"Readout position (0-1)"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            update_modal_weights.set();
            return { args };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> poisson_ratio { this, "poisson_ratio", 0.3,
        description {"Poisson's ratio"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            m_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> thickness { this, "thickness", 5e-4,
        description {"Thickness in meters"},
        range { 1e-4, 1e-2 },
        setter { MIN_FUNCTION {
            m_queue.set();
            return { args[0] };
        }}
    };
    
    attribute<number> lx { this, "lx", 0.2,
        description {"Length in meters"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Not implemented" << endl;
            return { args[0] };
        }}
    };

    attribute<number> ly { this, "ly", 0.3,
        description {"Width in meters"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Not implemented" << endl;
            return { args[0] };
        }}
    };

    attribute<number> density { this, "density", 7850,
        description {"Density in kg/m^3"},
        range { 1000.0, 10000.0 },
        setter { MIN_FUNCTION {
            m_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> youngs_modulus { this, "youngs_modulus", 2000,
        description {"Young's modulus in GPa"},
        range { 100.0, 10000.0 },
        setter { MIN_FUNCTION {
            m_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_independent_loss { this, "frequency_independent_loss", 0.01,
        description {"Frequency independent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            m_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_dependent_loss { this, "frequency_dependent_loss", 0.01,
        description {"Frequency dependent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            m_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> surface_tension { this, "surface_tension", 0.0,
        description {"Surface tension in N/m"},
        range { 0.0, 1000.0 },
        setter { MIN_FUNCTION {
            m_queue.set();
            return { args[0] };
        }}
    };  
    
    message<> maxclass_setup { this, "maxclass_setup",
        MIN_FUNCTION {
            cout << "plate~ - " << VERSION << " - 2025 - Rodrigo Diaz" << endl;
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

            m_initialized = true;
            calculate_coefficients();
            calculate_modal_weights();
            return {};
        }
    };
    
    sample operator()(sample input) {
        if (!m_initialized) {
            return 0.0;
        }

        std::unique_lock<std::mutex> lock {m_coeff_mutex};

        auto x = input * m_force_weights;

        Vector t0_flat = m_H_scaled * m_q;
        Matrix t0 = Eigen::Map<Matrix>(t0_flat.data(), m_n_psi, m_n_phi);
        Vector t2 = t0 * m_q;
        Vector nl = t0.transpose() * t2;

        Vector q_next = m_B.cwiseProduct(m_q) + m_C.cwiseProduct(m_q_prev) + m_A_inv.cwiseProduct(x - nl);
        
        m_q_prev = m_q;
        m_q = q_next;

        sample out = m_readout_weights.dot(q_next);

        lock.unlock();
        return out;
    }

private:
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

        double plate_norm = m_plate_parameters.l1 * m_plate_parameters.l2 * 0.25;
        double scale = (m_plate_parameters.E * plate_norm) / (2.0 * m_plate_parameters.density());

        std::unique_lock<std::mutex> lock {m_coeff_mutex};
        
        m_H_scaled = m_H_original * std::sqrt(scale);
        
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

    double m_sampleRate = 44100.0;
    double m_dt = 1.0 / m_sampleRate;
    bool m_initialized = false;
    int m_n_psi = 10;
    int m_n_phi = 10;
    
    std::mutex m_coeff_mutex;
    
    PlateParameters m_plate_parameters;
    Vector m_selected_indices_x;
    Vector m_selected_indices_y;
    
    Matrix m_H_original;
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

MIN_EXTERNAL(plate_tilde);