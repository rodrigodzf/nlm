#include "c74_min.h"
#include "c74_min_api.h"
#include "c74_min_attribute.h"
#include "c74_min_operator_sample.h"
#include "c74_min_queue.h"
#include <Eigen/Dense>
#include <matioCpp/matioCpp.h>
#include <matioCpp/MultiDimensionalArray.h>
#include <mutex>
#include <vk_utils/ModalSynthesizer.h>
#include "version.h"

using namespace c74::min;

class plate_tilde : public object<plate_tilde>, public sample_operator<1, 1> {
public:
    MIN_DESCRIPTION	{ "A non-linear modal plate." };
    MIN_TAGS		{ "audio" };
    MIN_AUTHOR		{ "Rodrigo Diaz" };
    MIN_RELATED		{ "biquad~" };

    plate_tilde(const atoms& args = {}) {
        // Initialize model type
        if (args.size() > 0) {
            if (args[0] == "vk") {
                m_model_type = "vk";
            } else {
                m_model_type = "berger";
            }
        }

        // Initialize n_phi and n_psi from arguments
        if (args.size() > 1) {
            int n_phi = static_cast<int>(args[1]);
            m_synthesizer.set_n_modes(n_phi);
        }
        if (args.size() > 2) {
            int n_psi = static_cast<int>(args[2]);
            m_synthesizer.set_n_psi(n_psi);
            m_n_psi = n_psi;
        }

        m_synthesizer.set_model_type(m_model_type);
    }

    inlet<>  left    { this, "(signal) Input to the modal resonator", "signal" };
    inlet<>  right   { this, "(dictionary) Dictionary with the parameters", "dictionary" };
    outlet<> output   { this, "(signal) Output from the modal resonator", "signal" };

    argument<symbol> model { this, "model", "Plate model type (vk or berger)", true,
        MIN_ARGUMENT_FUNCTION {
            m_model_type = std::string(arg);
            m_synthesizer.set_model_type(m_model_type);
        }
    };

    queue<> update_queue {this, MIN_FUNCTION {
        update_parameters();
        return {};
    }};

    // Position parameters (2D for plates)
    attribute<numbers> force_position { this, "force_position", {{0.5, 0.5}},
        description {"Force position (0-1)"},
        setter { MIN_FUNCTION {
            atoms result;
            for (auto& a : args)
                result.push_back(clamp<double>(a, 0.0, 1.0));
            
            if (result.size() >= 2) {
                m_force_position_x = result[0];
                m_force_position_y = result[1];
                m_position_changed = true;
                update_queue.set();
            }
            return result;
        }}
    };

    attribute<numbers> readout_position { this, "readout_position", {{0.5, 0.5}},
        description {"Readout position (0-1)"},
        setter { MIN_FUNCTION {
            atoms result;
            for (auto& a : args)
                result.push_back(clamp<double>(a, 0.0, 1.0));
                
            if (result.size() >= 2) {
                m_readout_position_x = result[0];
                m_readout_position_y = result[1];
                m_position_changed = true;
                update_queue.set();
            }
            return result;
        }}
    };

    // Material parameters
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

    // Geometry parameters
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

    // Messages for loading precomputed data
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

    // Dictionary parameter loading
    message<> dictionary { this, "dictionary", "dictionary with the parameters",
        MIN_FUNCTION {
            try {
                dict d = {args[0]};

                poisson_ratio = d["poisson_ratio"];
                thickness = d["thickness"];
                lx = d["lx"];
                ly = d["ly"];
                density = d["density"];
                youngs_modulus = d["youngs_modulus"];
                frequency_independent_loss = d["findependent_loss"];
                frequency_dependent_loss = d["fdependent_loss"];
                tension = d["tension"];
                
                // Manually ensure all flags are set
                m_geometry_changed = true;
                m_material_changed = true;
                m_position_changed = true; // Position weights depend on geometry/material parameters
                
                update_queue.set();
            }
            catch (std::exception& e) {
                cerr << e.what() << endl;
            }
            return {};
        }
    };

    message<> info { this, "info",
        MIN_FUNCTION {
            std::stringstream ss;
            ss << "modal.plate~ parameters:" << '\n'
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
               << "Force position: " << force_position[0] << " " << force_position[1] << '\n'
               << "Readout position: " << readout_position[0] << " " << readout_position[1] << '\n';
            cout << ss.str() << endl;
            return {};
        }
    };

    message<> maxclass_setup { this, "maxclass_setup",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            cout << "modal.plate~ - " << VERSION << " - 2025 - Rodrigo Diaz" << endl;
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

    message<> reset_coefficients { this, "reset_coefficients",
        MIN_FUNCTION {
            std::unique_lock<std::mutex> lock(m_update_mutex);
            m_synthesizer.reset_filter();
            return {};
        }
    };

    sample operator()(sample input) {
        if (!m_initialized) {
            return 0.0;
        }
        
        std::unique_lock<std::mutex> lock(m_update_mutex);
        return m_synthesizer.process_sample(input);
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
        
        m_synthesizer.set_force_position(m_force_position_x, m_force_position_y);
        m_synthesizer.set_readout_position(m_readout_position_x, m_readout_position_y);
        m_synthesizer.update_position_parameters(0.0, 0.0); // Dummy call to trigger weight calculation
        m_position_changed = false;
    }

    PlateModalSynthesizer m_synthesizer;
    std::mutex m_update_mutex;
    
    // Parameter storage
    PlateParameters m_parameters;
    double m_force_position_x = 0.5;
    double m_force_position_y = 0.5;
    double m_readout_position_x = 0.5;
    double m_readout_position_y = 0.5;
    
    // Model type and configuration
    std::string m_model_type = "berger";
    int m_n_psi = 10;
    
    // Update flags (handled externally, not by synthesizer)
    bool m_material_changed = true;
    bool m_geometry_changed = true;
    bool m_position_changed = true;
    bool m_initialized = false;
};

MIN_EXTERNAL(plate_tilde);