#include "c74_min.h"
#include "c74_min_api.h"
#include "c74_min_attribute.h"
#include "c74_min_queue.h"
#include <Eigen/Dense>
#include <matioCpp/matioCpp.h>
#include <matioCpp/MultiDimensionalArray.h>
#include <cmath>
#include <mutex>
#include <vk_utils/ModalSynthesizer.h>
#include "version.h"

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

    int m_num_inputs;
    int m_num_outputs;
    std::string m_model_type;
    MultiChannelPlateModalSynthesizer m_synthesizer;

    mcs_plate_tilde(const atoms& args = {}) 
        : m_num_inputs(args.size() > 2 ? static_cast<int>(args[1]) : 1),
          m_num_outputs(args.size() > 2 ? static_cast<int>(args[2]) : 2),
          m_model_type(args.size() > 0 && args[0] == "vk" ? "vk" : "berger"),
          m_synthesizer(64, m_num_inputs, m_num_outputs, m_model_type) {
        
        // Update the attribute
        model_type = m_model_type;


        // Initialize positions
        Vector force_pos_x = Vector::Constant(m_num_inputs, 0.5);
        Vector force_pos_y = Vector::Constant(m_num_inputs, 0.5);
        Vector readout_pos_x = Vector::Constant(m_num_outputs, 0.5);
        Vector readout_pos_y = Vector::Constant(m_num_outputs, 0.5);
        m_synthesizer.set_force_positions(force_pos_x, force_pos_y);
        m_synthesizer.set_readout_positions(readout_pos_x, readout_pos_y);
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
        update_parameters();
        return {};
    }};

    attribute<number, threadsafe::no, limit::clamp> poisson_ratio { this, "poisson_ratio", 0.3,
        description {"Poisson's ratio"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            m_parameters.nu = args[0];
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> thickness { this, "thickness", 1.9e-4,
        description {"Thickness in meters"},
        range { 1e-4, 1e-2 },
        setter { MIN_FUNCTION {
            m_parameters.h = args[0];
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };
    
    attribute<number, threadsafe::no, limit::clamp> lx { this, "lx", 0.4,
        description {"Length in meters"},
        range { 0.001, 1.0 },
        setter { MIN_FUNCTION {
            m_parameters.l1 = args[0];
            m_geometry_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> ly { this, "ly", 0.5,
        description {"Width in meters"},
        range { 0.001, 1.0 },
        setter { MIN_FUNCTION {
            m_parameters.l2 = args[0];
            m_geometry_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> density { this, "density", 3000,
        description {"Density in kg/m^3"},
        range { 1000.0, 10000.0 },
        setter { MIN_FUNCTION {
            m_parameters.rho = args[0];
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> youngs_modulus { this, "youngs_modulus", 3.5,
        description {"Young's modulus in GPa"},
        range { 1.0, 10000.0 },
        setter { MIN_FUNCTION {
            m_parameters.E = double(args[0]) * 1e9; // GPa to Pa
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_independent_loss { this, "findependent_loss", 0.01,
        description {"Frequency independent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            m_parameters.d1 = args[0];
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_dependent_loss { this, "fdependent_loss", 0.008,
        description {"Frequency dependent loss"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            m_parameters.d3 = args[0];
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> tension { this, "tension", 3000.0,
        description {"Tension in N/m"},
        range { 0.0, 50000.0 },
        setter { MIN_FUNCTION {
            m_parameters.Ts0 = args[0];
            m_material_changed = true;
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
            std::unique_lock<std::mutex> lock(m_update_mutex);
            // Set sample rate once during setup
            m_synthesizer.set_sample_rate(samplerate());
            // Mark all parameters for update
            m_material_changed = true;
            m_geometry_changed = true;
            m_position_changed = true;
            update_queue.set();
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
            
            try {
                matioCpp::File file(filename);
                if (!file.isOpen()) {
                    cout << "set_couplings_and_eigenvalues: could not open file '" << filename << "'" << endl;
                    return {};
                }

                
                auto H_mat = file.read("H").asMultiDimensionalArray<double>();
                auto lambda_mu_mat = file.read("lambda_mu").asVector<double>();
                auto selected_indices_x = file.read("selected_indices_x").asVector<double>();
                auto selected_indices_y = file.read("selected_indices_y").asVector<double>();

                Vector lambda_mu = matioCpp::to_eigen(lambda_mu_mat);
                Vector indices_x = matioCpp::to_eigen(selected_indices_x);
                Vector indices_y = matioCpp::to_eigen(selected_indices_y);

                int n_psi = H_mat.dimensions()[0];
                int n_phi = H_mat.dimensions()[1];
                
                Matrix H_original = Eigen::Map<Matrix>(H_mat.data(), n_psi * n_phi, n_phi);
                
                std::unique_lock<std::mutex> lock(m_update_mutex);
                m_synthesizer.set_precomputed_data(H_original, lambda_mu, indices_x, indices_y);
                m_synthesizer.set_n_psi(n_psi);
                m_n_psi = n_psi;
                
                
                // Force all updates after loading data
                m_geometry_changed = true;
                m_material_changed = true;
                m_position_changed = true;
                update_queue.set();
            }
            catch (std::exception& e) {
                cout << "Error loading precomputed data: " << e.what() << endl;
            }
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
               << "Surface tension: " << tension << " N/m" << '\n'
               << "Loss parameters: " << frequency_independent_loss << " (freq. indep.), " 
                 << frequency_dependent_loss << " (freq. dep.)" << '\n'
               << "Sampling rate: " << samplerate() << " Hz" << '\n'
               << "Modes: " << m_synthesizer.get_n_modes() << '\n'
               << "Number of inputs: " << m_num_inputs << '\n'
               << "Number of outputs: " << m_num_outputs << '\n';
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
            
            Vector positions_x(m_num_inputs);
            Vector positions_y(m_num_inputs);
            for (int i = 0; i < args.size() / 2; ++i) {
                // Clamp values between 0 and 1
                positions_x[i] = std::clamp(double(args[i * 2]), 0.0, 1.0);
                positions_y[i] = std::clamp(double(args[i * 2 + 1]), 0.0, 1.0);
            }
            
            m_synthesizer.set_force_positions(positions_x, positions_y);
            m_position_changed = true;
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

            Vector positions_x(m_num_outputs);
            Vector positions_y(m_num_outputs);
            for (int i = 0; i < args.size() / 2; ++i) {
                // Clamp values between 0 and 1
                positions_x[i] = std::clamp(double(args[i * 2]), 0.0, 1.0);
                positions_y[i] = std::clamp(double(args[i * 2 + 1]), 0.0, 1.0);
            }

            m_synthesizer.set_readout_positions(positions_x, positions_y);
            m_position_changed = true;
            update_queue.set();
            return {};
        }
    };

    message<> set_n_phi { this, "set_n_phi",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            if (args.size() < 1) return {};
            
            int n_modes = args[0];
            if (n_modes > 0) {
                std::unique_lock<std::mutex> lock(m_update_mutex);
                m_synthesizer.set_n_modes(n_modes);
                m_geometry_changed = true;
                update_queue.set();
            }
            return {};
        }
    };

    message<> set_n_psi { this, "set_n_psi",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            if (args.size() < 1) return {};
            
            int n_psi = args[0];
            if (n_psi > 0) {
                std::unique_lock<std::mutex> lock(m_update_mutex);
                m_synthesizer.set_n_psi(n_psi);
                m_n_psi = n_psi;
                m_geometry_changed = true;
                update_queue.set();
            }
            return {};
        }
    };

    message<> reset_coefficients { this, "reset_coefficients",
        MIN_FUNCTION {
            std::unique_lock<std::mutex> lock(m_update_mutex);
            m_synthesizer.reset_filter();
            cout << "Reset parallel filter state" << endl;
            return {};
        }
    };

    void operator()(audio_bundle input, audio_bundle output) {
        if (!m_initialized) {
            for (int frame = 0; frame < input.frame_count(); ++frame) {
                for (int out_ch = 0; out_ch < m_num_outputs; ++out_ch) {
                    output.samples(out_ch)[frame] = 0.0;
                }
            }
            return;
        }
        
        std::unique_lock<std::mutex> lock(m_update_mutex);

        for (int frame = 0; frame < input.frame_count(); ++frame) {
            Vector inputs(m_num_inputs);
            for (int in_ch = 0; in_ch < m_num_inputs; ++in_ch) {
                inputs[in_ch] = input.samples(in_ch)[frame];
            }

            m_synthesizer.process_frame_multichannel(inputs);

            for (int out_ch = 0; out_ch < m_num_outputs; ++out_ch) {
                output.samples(out_ch)[frame] = m_synthesizer.get_output_multichannel(out_ch);
            }
        }
    }

private:
    void update_parameters() {
        std::unique_lock<std::mutex> lock(m_update_mutex);
        
        update_geometry_if_needed();
        update_material_if_needed();
        update_position_if_needed();
        m_initialized = true;
    }
    
    void update_geometry_if_needed() {
        if (!m_geometry_changed) return;
        
        m_synthesizer.update_geometry_parameters(m_parameters);
        m_geometry_changed = false;
    }
    
    void update_material_if_needed() {
        if (!m_material_changed) return;
        
        m_synthesizer.update_material_parameters(m_parameters);
        m_material_changed = false;
    }
    
    void update_position_if_needed() {
        if (!m_position_changed) return;
        
        m_synthesizer.update_multichannel_position_parameters();
        m_position_changed = false;
    }

    std::mutex m_update_mutex;
    
    // Parameter storage
    PlateParameters m_parameters;
    
    // Model type and configuration
    int m_n_psi = 10;
    
    // Update flags (handled externally, not by synthesizer)
    bool m_material_changed = true;
    bool m_geometry_changed = true;
    bool m_position_changed = true;
    bool m_initialized = false;
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