#include "c74_min.h"
#include "c74_min_attribute.h"
#include "c74_min_operator_sample.h"
#include "matioCpp/MultiDimensionalArray.h"
#include <cmath>
#include <Eigen/Dense>
#include <matioCpp/matioCpp.h>
#include <vk_utils/TimeIntegrators.h>
#include <vk_utils/Parameters.h>
using namespace c74::min;

class modal_resonator : public object<modal_resonator>, public sample_operator<1, 1> {
public:
    MIN_DESCRIPTION	{ "Modal resonator." };
    MIN_TAGS		{ "audio" };
    MIN_AUTHOR		{ "Rodrigo Diaz" };
    MIN_RELATED		{ "modal_resonator~" };

    // Constructor
    modal_resonator(const atoms& args = {}) {
        // Initialize any variables or parameters here
        calculate_coefficients();
    }

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

    attribute<number> force_position_x { this, "force_position_x", 0.05,
        description {"Force position (0-1)"},
        range { 0.0, 1.0 }
    };

    attribute<number> force_position_y { this, "force_position_y", 0.05,
        description {"Force position (0-1)"},
        range { 0.0, 1.0 }
    };

    attribute<number> readout_position_x { this, "readout_position_x", 0.1,
        description {"Readout position (0-1)"},
        range { 0.0, 1.0 }
    };

    attribute<number> readout_position_y { this, "readout_position_y", 0.1,
        description {"Readout position (0-1)"},
        range { 0.0, 1.0 }
    };

    attribute<number> poisson_ratio { this, "poisson_ratio", 0.3,
        description {"Poisson's ratio"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting poisson_ratio to " << args[0] << endl;
            m_plate_parameters.nu = args[0];
            return {
                m_plate_parameters.nu
            };
        } }
    };

    attribute<number> thickness { this, "thickness", 5e-4,
        description {"Thickness in meters"},
        range { 1e-4, 1e-3 },
        setter { MIN_FUNCTION {
            cout << "Setting thickness to " << args[0] << endl;
            m_plate_parameters.h = args[0];
            return {
                m_plate_parameters.h
            };
        } }
    };
    
    attribute<number> lx { this, "lx", 0.2,
        description {"Length in meters"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting lx to " << args[0] << endl;
            m_plate_parameters.l1 = args[0];
            return {
                m_plate_parameters.l1
            };
        } }
    };

    attribute<number> ly { this, "ly", 0.3,
        description {"Width in meters"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting ly to " << args[0] << endl;
            m_plate_parameters.l2 = args[0];
            return {
                m_plate_parameters.l2
            };
        } }
    };

    attribute<number> density { this, "density", 7850,
        description {"Density in kg/m^3"},
        range { 1000.0, 10000.0 },
        setter { MIN_FUNCTION {
            cout << "Setting density to " << args[0] << endl;
            m_plate_parameters.rho = args[0];
            return {
                m_plate_parameters.rho
            };
        } }
    };

    attribute<number> youngs_modulus { this, "youngs_modulus", 2000,
        description {"Young's modulus in GPa"},
        range { 100.0, 10000.0 },
        setter { MIN_FUNCTION {
            cout << "Setting youngs_modulus to " << args[0] << endl;
            m_plate_parameters.E = static_cast<double>(args[0]) * 1e9;
            return {
                m_plate_parameters.E
            };
        } }
    };

    attribute<number> frequency_independent_loss { this, "frequency_independent_loss", 4.2e-6,
        description {"Frequency independent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting frequency_independent_loss to " << args[0] << endl;
            m_plate_parameters.d1 = args[0];
            return {
                m_plate_parameters.d1
            };
        } }
    };

    attribute<number> frequency_dependent_loss { this, "frequency_dependent_loss", 2.3e-3,
        description {"Frequency dependent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            cout << "Setting frequency_dependent_loss to " << args[0] << endl;
            m_plate_parameters.d3 = args[0];
            return {
                m_plate_parameters.d3
            };
        } }
    };

    attribute<number> surface_tension { this, "surface_tension", 0.0,
        description {"Surface tension in N/m"},
        range { 0.0, 1000.0 },
        setter { MIN_FUNCTION {
            cout << "Setting surface tension to " << args[0] << endl;
            m_plate_parameters.Ts0 = args[0];
            return {
                m_plate_parameters.Ts0
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
    
    message<> bang { this, "bang",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            trigger();
            return {};
        }
    };
    
    message<> reset { this, "reset",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            reset_state();
            return {};
        }
    };

    // Message for setting parameters
    message<> set_force_position { this, "set_force_position",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            force_position_x = args[0];
            force_position_y = args[1];
            return {};
        }
    };

    message<> set_readout_position { this, "set_readout_position",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            readout_position_x = args[0];
            readout_position_y = args[1];
            return {};
        }
    };
    
    message<> set_poisson_ratio { this, "set_poisson_ratio",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            poisson_ratio = args[0];
            return {};
        }
    };

    message<> set_thickness { this, "set_thickness",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            thickness = args[0];
            return {};
        }
    };
    
    message<> set_lx { this, "set_lx",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            lx = args[0];
            return {};
        }
    };

    message<> set_ly { this, "set_ly",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            ly = args[0];
            return {};
        }
    };
    
    message<> set_density { this, "set_density",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            density = args[0];
            return {};
        }
    };  

    message<> set_youngs_modulus { this, "set_youngs_modulus",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            youngs_modulus = args[0];
            return {};
        }
    };

    message<> set_frequency_independent_loss { this, "set_frequency_independent_loss",  
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            frequency_independent_loss = args[0];
            return {};
        }
    };      

    message<> set_frequency_dependent_loss { this, "set_frequency_dependent_loss",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            frequency_dependent_loss = args[0];
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
            matioCpp::Element<double> ny_elem = file.read("ny").asElement<double>();
            matioCpp::Element<double> nx_elem = file.read("nx").asElement<double>();
            matioCpp::MultiDimensionalArray<double> phi_mat = file.read("phi").asMultiDimensionalArray<double>();
            
            cout << "nx_elem: " << nx_elem << endl;
            cout << "ny_elem: " << ny_elem << endl;


            m_lambda_mu = matioCpp::to_eigen(lambda_mu_mat);

            // get the size of lambda_mu
            int size_lambda_mu = m_lambda_mu.size();
            cout << "Size of lambda_mu: " << size_lambda_mu << endl;
            m_n_modes = size_lambda_mu;

            // convert to Eigen::Matrix
            m_H = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(H_mat.data(), m_n_modes * m_n_modes, m_n_modes);
            m_phi = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(phi_mat.data(), phi_mat.dimensions()[0], phi_mat.dimensions()[1]);


            
            // phi has shape (ny * nx, n_modes) but we need to index it with ny and nx
            int force_nx = static_cast<int>(std::round(force_position_x * nx_elem));
            int force_ny = static_cast<int>(std::round(force_position_y * ny_elem));
            int force_row_idx = force_ny * nx_elem + force_nx;
            cout << "force_row_idx: " << force_row_idx << endl;

            int readout_nx = static_cast<int>(std::round(readout_position_x * nx_elem));
            int readout_ny = static_cast<int>(std::round(readout_position_y * ny_elem));
            int readout_row_idx = readout_ny * nx_elem + readout_nx;
            cout << "readout_row_idx: " << readout_row_idx << endl;
            m_force_weights = m_phi.row(force_row_idx);
            m_readout_weights = m_phi.row(readout_row_idx);

            cout << "force_weights: " << m_force_weights << endl;
            cout << "readout_weights: " << m_readout_weights << endl;

            // divide m_force_weights by m_plate_parameters.density
            m_force_weights = m_force_weights / m_plate_parameters.density();

            // get the size of H
            int size = m_H.rows();
            cout << "Size of H: " << size << endl;



            // print m_plate_parameters
            cout << "m_plate_parameters: " << m_plate_parameters << endl;
            m_omega_mu_squared = stiffness_term<double>(m_plate_parameters, m_lambda_mu);
            m_gamma2_mu = damping_term<double>(m_plate_parameters, m_lambda_mu);

            cout << "omega_mu: " << m_omega_mu_squared.cwiseSqrt() << endl;
            cout << "gamma2_mu: " << m_gamma2_mu << endl;

            m_A_inv = A_inv_vector(m_dt, m_gamma2_mu);
            m_B = B_vector(m_dt, m_omega_mu_squared).array() * m_A_inv.array();
            m_C = C_vector(m_dt, m_gamma2_mu).array() * m_A_inv.array();

            cout << "A_inv: " << m_A_inv << endl;
            cout << "B: " << m_B << endl;
            cout << "C: " << m_C << endl;

            cout << "m_n_modes: " << m_n_modes << endl;

            reset_state();

            // test the operator() with a few samples
            // sample out = operator()(1.0);
            // cout << "out: " << out << endl;
            // sample out2 = operator()(1.0);
            // cout << "out2: " << out2 << endl;
            // sample out3 = operator()(1.0);
            // cout << "out3: " << out3 << endl;
            // sample out4 = operator()(1.0);
            // cout << "out4: " << out4 << endl;

            m_initialized = true;
            return {};
        }
    };
    

    // DSP perform method (for Max/MSP audio processing)
    sample operator()(sample input) {

        if (!m_initialized) {
            return 0.0;
        }

        sample out = 0.0;

        // Each input is a scalar which weights the modal gains at the excitation position
        // That is this is a scalar times a vector
        auto x = input * m_force_weights;

        Eigen::Matrix<double, Eigen::Dynamic, 1> t0_flat = m_H * m_q;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> t0 =
            Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(t0_flat.data(), m_q.size(), m_q.size());
        Eigen::Matrix<double, Eigen::Dynamic, 1> t2 = t0 * m_q;
        Eigen::Matrix<double, Eigen::Dynamic, 1> nl = t0.transpose() * t2;


        Eigen::Matrix<double, Eigen::Dynamic, 1> q_next = m_B.cwiseProduct(m_q) + m_C.cwiseProduct(m_q_prev) + m_A_inv.cwiseProduct(x - nl);
        
        // Update previous state
        m_q_prev = m_q;
        m_q = q_next;

        // Get a scalar output at the readout position
        out = m_readout_weights.dot(q_next);

        return out;

        // const int vector_size = input.frame_count();
        
        // Store current values to check for changes
        // static double last_frequency = 0.0;
        // static double last_decay = 0.0;

        
        // Check if parameters have changed
        // if (!m_initialized || 
        //     last_frequency != static_cast<double>(frequency) || 
        //     last_decay != static_cast<double>(decay)) {
            
        //     calculate_coefficients();
            
        //     // Update stored values
        //     last_frequency = static_cast<double>(frequency);
        //     last_decay = static_cast<double>(decay);
        // }
        
        // auto in = input.samples(0);   // Get input samples
        // auto out = output.samples(0); // Get output buffer
        
        // Process each sample
        // for (int i = 0; i < vector_size; ++i) {
            // Apply biquad filter equation: y[n] = b0*x[n] - a1*y[n-1] - a2*y[n-2]
            // double y0 = m_b0 * in[i] - m_a1 * m_y1 - m_a2 * m_y2;
            
            // Store output with gain applied
            // out[i] = y0 * gain;
            
            // Update delay registers
            // m_y2 = m_y1;
            // m_y1 = y0;
            // out[i] = in[i];
        // }
    }

    // We could potentially make these public methods that take atoms directly
    // Then we wouldn't need wrapper lambdas in our messages
    atoms handle_bang(const atoms& args, int inlet) {
        trigger();
        return {};
    }
    
    atoms handle_reset(const atoms& args, int inlet) {
        reset_state();
        return {};
    }

private:

    double bending_stiffness() const {
        return youngs_modulus * std::pow(thickness, 3) / (12 * (1 - poisson_ratio * poisson_ratio));
    }

    double mass_density() const {
        return density * thickness;
    }

    // Private methods
    void setup(const atoms& args, const int vector_size) {
        // Get current sample rate from Max
        if (!args.empty())
            m_sampleRate = args[0];
        
        // Recalculate filter coefficients based on new sample rate
        calculate_coefficients();
        
        // Reset filter state
        reset_state();
    }
    
    void trigger() {
        // Trigger with an impulse (1.0)
        // m_y1 = 1.0;
        // m_y2 = 0.0;
    }
    
    void reset_state() {
        m_q = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(m_n_modes);
        m_q_prev = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(m_n_modes);
        // m_y1 = 0.0;
        // m_y2 = 0.0;
    }
    
    void calculate_coefficients() {
        // // Convert frequency to normalized frequency (0-1)
        // double omega = 2.0 * M_PI * frequency / m_sampleRate;
        
        // // Calculate coefficients for a resonant filter
        // double r = decay;
        // m_a1 = -2.0 * r * cos(omega);
        // m_a2 = r * r;
        // m_b0 = 1.0 - r * r;  // Scale factor for unity gain at resonance
        
        // m_initialized = true;
    }

    // Private member variables
    double m_y1 = 0.0;
    double m_y2 = 0.0;
    double m_a1 = 0.0;
    double m_a2 = 0.0;
    double m_b0 = 1.0;
    double m_sampleRate = 44100.0;
    double m_dt = 1.0 / m_sampleRate;

    bool m_initialized = false;
    int m_n_modes = 0;

    PlateParameters m_plate_parameters;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_phi;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_H;
    Eigen::Matrix<double, Eigen::Dynamic, 1> m_lambda_mu;
    Eigen::Matrix<double, Eigen::Dynamic, 1> m_gamma2_mu;
    Eigen::Matrix<double, Eigen::Dynamic, 1> m_omega_mu_squared;
    Eigen::Matrix<double, Eigen::Dynamic, 1> m_A_inv;
    Eigen::Matrix<double, Eigen::Dynamic, 1> m_B;
    Eigen::Matrix<double, Eigen::Dynamic, 1> m_C;
    Eigen::Matrix<double, Eigen::Dynamic, 1> m_q;
    Eigen::Matrix<double, Eigen::Dynamic, 1> m_q_prev;
    Eigen::Matrix<double, Eigen::Dynamic, 1> m_force_weights;
    Eigen::Matrix<double, Eigen::Dynamic, 1> m_readout_weights;

};

MIN_EXTERNAL(modal_resonator);