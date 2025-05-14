#include "c74_min.h"
#include "c74_min_api.h"
#include "c74_min_attribute.h"
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
#include <vk_utils/ParallelFilter.h>

using namespace c74::min;

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

long simplemc_multichanneloutputs(c74::max::t_object* x, long index, long count);

class mcs_plate_tilde : public object<mcs_plate_tilde>, public vector_operator<> {
public:
    MIN_DESCRIPTION	{ "A modal plate." };
    MIN_TAGS		{ "audio" };
    MIN_AUTHOR		{ "Rodrigo Diaz" };
    MIN_RELATED		{ "biquad~" };

    mcs_plate_tilde(const atoms& args = {}) {

        
        if (args.size() > 2) {
            m_num_inputs = static_cast<int>(args[1]);
            m_num_outputs = static_cast<int>(args[2]);
        }


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
        cout << "Number of inputs: " << m_num_inputs << endl;
        cout << "Number of outputs: " << m_num_outputs << endl;

        // Default number of modes for berger model
        if (m_model_type == "berger") {
            m_n_phi = 32;
        }

        m_lambda_mu = Vector::Zero(m_n_phi);
        m_gamma2_mu = Vector::Zero(m_n_phi);

        m_force_positions_x = Vector::Constant(m_num_inputs, 0.5);
        m_force_positions_y = Vector::Constant(m_num_inputs, 0.5);
        m_readout_positions_x = Vector::Constant(m_num_outputs, 0.5);
        m_readout_positions_y = Vector::Constant(m_num_outputs, 0.5);
    }

    inlet<>  input    { this, "(multichannelsignal) Input to the modal resonator" };
    outlet<> output   { this, "(multichannelsignal) Output from the modal resonator", "multichannelsignal" };

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
            // Reset eigenmode calculation flag when length changes
            if (m_model_type == "berger" && !m_data_loaded) {
                m_eigenmode_calculation_complete = false;
            }
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> ly { this, "ly", 0.3,
        description {"Width in meters"},
        range { 0.001, 1.0 },
        setter { MIN_FUNCTION {
            // Reset eigenmode calculation flag when width changes
            if (m_model_type == "berger" && !m_data_loaded) {
                m_eigenmode_calculation_complete = false;
            }
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
            cout << "mcs.modal.plate~ - " << VERSION << " - 2025 - Rodrigo Diaz" << endl;
            
            c74::max::t_class* c = args[0];
            c74::max::class_addmethod(c, (c74::max::method)simplemc_multichanneloutputs, "multichanneloutputs", c74::max::A_CANT, 0);
            // c74::max::class_addmethod(c, (c74::max::method)simplemc_inputchanged, "inputchanged", c74::max::A_CANT, 0);
            return {};
        }
    };

    message<> dspsetup {this, "dspsetup",
        MIN_FUNCTION {
            m_dt = 1.0 / samplerate();
            // Don't reset initialized flag if data is loaded from file
            if (m_initialized) {
                // allocate_temp_vectors(); // Ensure temp vectors are allocated with correct size
                update_queue.set();       // Update parameters but preserve loaded data
            }
            return {};
        }
    };
    message<> set_couplings_and_eigenvalues { this, "set_couplings_and_eigenvalues",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            if (args.size() < 1) {
                cout << "set_couplings_and_eigenvalues: no filename provided" << endl;
                return {};
            }

            std::string filename = args[0];
            
            // Check if file exists
            std::ifstream file_check(filename);
            if (!file_check.good()) {
                cout << "set_couplings_and_eigenvalues: file '" << filename << "' does not exist or is not accessible" << endl;
                return {};
            }
            file_check.close();

            matioCpp::File file(filename);
            if (!file.isOpen()) {
                cout << "set_couplings_and_eigenvalues: could not open file '" << filename << "'" << endl;
                return {};
            }

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
            m_readout_weights_lerp.resize(m_n_phi, m_num_outputs, 441);
            m_current_readout_weights.resize(m_n_phi, m_num_outputs);

            // Only resize parallel filter if needed
            if (m_parallel_filter.get_n_modes() != m_n_phi) {
                m_parallel_filter.resize(m_n_phi);
            }

            m_force_input.resize(m_n_phi, m_num_inputs);
            update_queue.set(); // Trigger parameter update after data loading
            return {};
        }
    };
    
    message<> info { this, "info",
        MIN_FUNCTION {
            std::stringstream ss;
            ss << "mcs.modal.plate~ parameters:" << '\n'
               << "Model type: " << m_model_type << '\n'
               << "Thickness: " << thickness << " m" << '\n'
               << "Dimensions: " << lx << " x " << ly << " m" << '\n'
               << "Density: " << density << " kg/m³" << '\n'
               << "Young's modulus: " << youngs_modulus << " GPa" << '\n'
               << "Poisson's ratio: " << poisson_ratio << '\n'
               << "Surface tension: " << surface_tension << " N/m" << '\n'
               << "Loss parameters: " << frequency_independent_loss << " (freq. indep.), " 
                 << frequency_dependent_loss << " (freq. dep.)" << '\n'
               << "Sampling rate: " << samplerate() << " Hz" << '\n'
               << "Modes: " << m_n_phi << '\n'
               << "Force positions: ";
            for (int i = 0; i < m_force_positions_x.size(); ++i) {
                ss << "(" << m_force_positions_x[i] << ", " << m_force_positions_y[i] << ") ";
            }
            ss << '\n'
               << "Readout positions: ";
            for (int i = 0; i < m_readout_positions_x.size(); ++i) {
                ss << "(" << m_readout_positions_x[i] << ", " << m_readout_positions_y[i] << ") ";
            }
            ss << '\n';
            cout << ss.str() << endl;
            return {};
        }
    };
    
    message<> force_position { this, "force_position",
        MIN_FUNCTION {
            if (args.size() % 2 != 0) {
                cerr << "force_position: must provide pairs of x,y coordinates" << endl;
                return {};
            }
            // check if args.size() / 2 is equal to m_num_inputs
            if (args.size() / 2 != m_num_inputs) {
                cerr << "force_position: must provide " << m_num_inputs << " pairs of x,y coordinates" << endl;
                return {};
            }
            
            for (int i = 0; i < args.size() / 2; ++i) {
                // Clamp values between 0 and 1
                m_force_positions_x[i] = std::clamp(double(args[i * 2]), 0.0, 1.0);
                m_force_positions_y[i] = std::clamp(double(args[i * 2 + 1]), 0.0, 1.0);
            }
            
            update_queue.set();
            return {};
        }
    };

    message<> readout_position { this, "readout_position",
        MIN_FUNCTION {
            if (args.size() % 2 != 0) {
                cerr << "readout_position: must provide pairs of x,y coordinates" << endl;
                return {};
            }
            // check if args.size() / 2 is equal to m_num_outputs
            if (args.size() / 2 != m_num_outputs) {
                cerr << "readout_position: must provide " << m_num_outputs << " pairs of x,y coordinates" << endl;
                return {};
            }

            for (int i = 0; i < args.size() / 2; ++i) {
                // Clamp values between 0 and 1
                m_readout_positions_x[i] = std::clamp(double(args[i * 2]), 0.0, 1.0);
                m_readout_positions_y[i] = std::clamp(double(args[i * 2 + 1]), 0.0, 1.0);
            }

            update_queue.set();
            return {};
        }
    };

    message<> set_num_modes { this, "set_num_modes",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            if (args.size() < 1) {
                cout << "set_num_modes: no mode count provided" << endl;
                return {};
            }

            int num_modes = args[0];
            if (num_modes < 1) {
                cout << "set_num_modes: mode count must be positive" << endl;
                return {};
            }

            // Reset eigenmode calculation flag when number of modes changes
            if (m_model_type == "berger" && !m_data_loaded && m_n_phi != num_modes) {
                m_eigenmode_calculation_complete = false;
            }
            
            m_n_phi = num_modes;
            cout << "Number of modes set to " << m_n_phi << endl;
            
            // If using berger model, recalculate modes
            if (m_model_type == "berger") {
                update_queue.set();
            }
            
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

            // Nonlinear update
            if (m_model_type == "vk") {
                m_parallel_filter.update_nonlinearity(m_H_scaled, m_n_psi, m_n_phi);
            } else if (m_model_type == "berger") {
                m_parallel_filter.update_nonlinearity(m_lambda_mu, m_plate_tau_with_norms);
            }
            // add the contribution for each input channel
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
    
    // Unified function to update all parameters
    void update_all_parameters() {
        // Can't update VK model without data loaded from file
        if (m_model_type == "vk" && !m_data_loaded) {
            return;
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
        
        // For berger model, initialize vectors if needed
        if (m_model_type == "berger") {
            // Initialize vectors for all berger model cases
            if (m_readout_weights.size() == 0 || m_readout_weights.rows() != m_n_phi || m_readout_weights.cols() != m_num_outputs) {
                m_readout_weights_lerp.resize(m_n_phi, m_num_outputs, 441);
            }
            
            if (m_current_readout_weights.rows() != m_n_phi || m_current_readout_weights.cols() != m_num_outputs) {
                m_current_readout_weights.resize(m_n_phi, m_num_outputs);
            }
            
            // Only resize parallel filter if needed
            if (m_parallel_filter.get_n_modes() != m_n_phi) {
                m_parallel_filter.resize(m_n_phi);
            }
            
            if (m_force_input.rows() != m_n_phi || m_force_input.cols() != m_num_inputs) {
                m_force_input.resize(m_n_phi, m_num_inputs);
            }
            
            // For berger model, calculate eigenvalues and select modes ONLY if data is not loaded from file
            // and if eigenmode calculation is not already complete
            if (!m_data_loaded && !m_eigenmode_calculation_complete) {
                // Set number of modes to consider in calculation (larger than needed to ensure we have enough)
                int n_max_modes_x = 20;
                int n_max_modes_y = 20;
                
                // Calculate eigenvalues matrix
                Matrix lambda_mu_2d = calculate_plate_eigenvalues<double>(
                    n_max_modes_x, n_max_modes_y, m_plate_parameters.l1, m_plate_parameters.l2);
                
                // Allocate vectors for selected modes
                m_lambda_mu.resize(m_n_phi);
                Eigen::MatrixX<int> selected_indices(m_n_phi, 2);
                
                // Select modes and eigenvalues
                select_modes_and_eigenvalues<double>(lambda_mu_2d, m_n_phi, m_lambda_mu, selected_indices);
                
                // Extract indices into separate x and y vectors
                m_selected_indices_x.resize(m_n_phi);
                m_selected_indices_y.resize(m_n_phi);
                
                for (int i = 0; i < m_n_phi; ++i) {
                    m_selected_indices_x(i) = selected_indices(i, 0);
                    m_selected_indices_y(i) = selected_indices(i, 1);
                }
                
                cout << "Berger model: calculated " << m_n_phi << " modes" << endl;
                
                // Mark eigenmode calculation as complete
                m_eigenmode_calculation_complete = true;
            }
        }
        
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
     
        m_parallel_filter.set_coefficients(m_gamma2_mu, m_omega_mu_squared, m_dt);
        
        // Calculate force weights
        m_force_weights = ftm::evaluate_rectangular_eigenfunctions(
            m_selected_indices_x, m_selected_indices_y,
            (m_force_positions_x * m_plate_parameters.l1).eval(),
            (m_force_positions_y * m_plate_parameters.l2).eval(),
            m_plate_parameters.l1, m_plate_parameters.l2
        ) / (plate_norm * m_plate_parameters.density());

        // Calculate readout weights
        Eigen::MatrixXd new_readout_weights = ftm::evaluate_rectangular_eigenfunctions(
            m_selected_indices_x, m_selected_indices_y,
            (m_readout_positions_x * m_plate_parameters.l1).eval(),
            (m_readout_positions_y * m_plate_parameters.l2).eval(),
            m_plate_parameters.l1, m_plate_parameters.l2
        );
        
        // Comment out the print statements to reduce console output
        // cout << "First 5 rows of new_readout_weights:" << endl;
        // for (int col = 0; col < new_readout_weights.cols(); ++col) {
        //     cout << "Column " << col << ": ";
        //     for (int row = 0; row < 5 && row < new_readout_weights.rows(); ++row) {
        //         cout << new_readout_weights(row, col) << " ";
        //     }
        //     cout << endl;
        // }

        if (m_readout_weights_lerp.isFinished()) {
            if (m_readout_weights.size() == 0) {
                m_readout_weights = new_readout_weights;
                m_readout_weights_lerp.setValue(new_readout_weights);
            }
            m_readout_weights_lerp.setTarget(new_readout_weights);
        }
        bool success = m_readout_weights_lerp.process(m_current_readout_weights);
        if (!success) {
            cout << "Failed to process m_readout_weights_lerp" << endl;
        }

        // Comment out additional debugging prints
        // cout << "m_current_readout_weights: " << m_current_readout_weights.rows() << " x " << m_current_readout_weights.cols() << endl;
        // cout << "First 5 rows of m_current_readout_weights:" << endl;
        // for (int col = 0; col < m_current_readout_weights.cols(); ++col) {
        //     cout << "Column " << col << ": ";
        //     for (int row = 0; row < 5 && row < m_current_readout_weights.rows(); ++row) {
        //         cout << m_current_readout_weights(row, col) << " ";
        //     }
        //     cout << endl;
        // }

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
    
    Eigen::MatrixXd m_force_weights; // (modes, inputs)
    Eigen::MatrixXd m_readout_weights; // (modes, outputs)
    Eigen::MatrixXd m_current_readout_weights; // (modes, outputs)
    Eigen::MatrixXd m_force_input; // (modes, inputs)
    MatrixInterpolator<double> m_readout_weights_lerp;

    std::string m_model_type = "berger";  // Default model type

    // Track if data is loaded from file
    bool m_data_loaded = false;
    
    // Track if eigenmode calculation is complete
    bool m_eigenmode_calculation_complete = false;

    Vector m_force_positions_x;
    Vector m_force_positions_y;
    Vector m_readout_positions_x;
    Vector m_readout_positions_y;

    ParallelFilter<double> m_parallel_filter;
};

long simplemc_multichanneloutputs(c74::max::t_object* x, long index, long count) {
    minwrap<mcs_plate_tilde>* ob = (minwrap<mcs_plate_tilde>*)(x);
    return ob->m_min_object.m_num_outputs;
}

// long simplemc_inputchanged(c74::max::t_object* x, long index, long count) {
//     minwrap<mcs_plate_tilde>* ob = (minwrap<mcs_plate_tilde>*)(x);
//     auto chan_number = ob->m_min_object.m_num_inputs;
//     ob->m_min_object.input_chans[index] = count;
//     if (chan_number != count) {
//       c74::max::object_error(x, (std::string("invalid channel number for input ") + std::to_string(index)
//                                           + std::string("; should be ") + std::to_string(chan_number)).c_str());
//     } 
//     return false;
// }

MIN_EXTERNAL(mcs_plate_tilde);