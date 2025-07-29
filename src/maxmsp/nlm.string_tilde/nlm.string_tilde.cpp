#include "c74_min.h"
#include "c74_min_api.h"
#include "c74_min_attribute.h"
#include "c74_min_operator_sample.h"
#include "c74_min_queue.h"
#include <mutex>
#include <vk_utils/ModalSynthesizer.h>
#include "version.h"

using namespace c74::min;

class string_tilde : public object<string_tilde>, public sample_operator<1, 1> {
public:
    MIN_DESCRIPTION	{ "A modal string." };
    MIN_TAGS		{ "audio" };
    MIN_AUTHOR		{ "Rodrigo Diaz" };
    MIN_RELATED		{ "biquad~" };

    string_tilde(const atoms& args = {}) {
        if (!args.empty()) {
            m_synthesizer.set_n_modes(static_cast<int>(args[0]));
        }
    }

    inlet<>  input    { this, "(signal) Input to the modal resonator", "signal" };
    inlet<>  right    { this, "(dictionary) Dictionary with the parameters", "dictionary" };
    outlet<> output   { this, "(signal) Output from the modal resonator", "signal" };

    argument<number> nmodes { this, "nmodes", "Number of modes", true,
        MIN_ARGUMENT_FUNCTION {
            std::unique_lock<std::mutex> lock(m_update_mutex);
            m_synthesizer.set_n_modes(static_cast<int>(arg));
            m_geometry_changed = true;
        }
    };

    queue<> update_queue {this, MIN_FUNCTION {
        update_parameters();
        return {};
    }};

    // Position parameters - these update immediately
    attribute<number, threadsafe::no, limit::clamp> force_position { this, "force_position", 0.5,
        description {"Force position (0-1)"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            m_force_position = args[0];
            m_position_changed = true;
            update_queue.set();
            return { args };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> readout_position { this, "readout_position", 0.5,
        description {"Readout position (0-1)"},
        range { 0.0, 1.0 },
        setter { MIN_FUNCTION {
            m_readout_position = args[0];
            m_position_changed = true;
            update_queue.set();
            return { args };
        }}
    };

    // Geometry parameters
    attribute<number, threadsafe::no, limit::clamp> lx { this, "lx", 0.65,
        description {"Length in meters"},
        range { 0.001, 5.0 },
        setter { MIN_FUNCTION {
            m_parameters.length = args[0];
            m_geometry_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    // Material parameters
    attribute<number> cross_sectional_area { this, "cross_sectional_area", 0.5188,
        description {"Cross sectional area in mm^2"},
        range { 0.00001, 10.0 },
        setter { MIN_FUNCTION {
            m_parameters.A = double(args[0]) * 1e-6; // mm^2 to m^2
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> density { this, "density", 1140.0,
        description {"Mass density in kg/m^3"},
        range { 500.0, 80000.0 },
        setter { MIN_FUNCTION {
            m_parameters.rho = args[0];
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> youngs_modulus { this, "youngs_modulus", 5.4,
        description {"Young's modulus in GPa"},
        range { 1.0, 1000.0 },
        setter { MIN_FUNCTION {
            m_parameters.E = double(args[0]) * 1e9; // GPa to Pa
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_independent_loss { this, "findependent_loss", 8e-5,
        description {"Frequency independent loss in kg/(ms)"},
        range { 1e-7, 1.0 },
        setter { MIN_FUNCTION {
            m_parameters.d1 = args[0];
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number, threadsafe::no, limit::clamp> frequency_dependent_loss { this, "fdependent_loss", 3.4e-5,
        description {"Frequency dependent loss in kg m/s"},
        range { 1e-6, 1.0 },
        setter { MIN_FUNCTION {
            m_parameters.d3 = args[0];
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> moment_of_inertia { this, "moment_of_inertia", 0.171,
        description {"Moment of inertia in mm^4"},
        range { 0.1, 6.0 },
        setter { MIN_FUNCTION {
            m_parameters.I = double(args[0]) * 1e-12; // mm^4 to m^4
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    attribute<number> tension { this, "tension", 60.97,
        description {"Tension in N"},
        range { 0.0, 1000.0 },
        setter { MIN_FUNCTION {
            m_parameters.Ts0 = args[0];
            m_material_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

    // Dictionary parameter loading  
    message<> dictionary { this, "dictionary", "dictionary with the parameters",
        MIN_FUNCTION {
            try {
                dict d = {args[0]};

                lx = d["lx"];
                cross_sectional_area = d["cross_sectional_area"];
                density = d["density"];
                youngs_modulus = d["youngs_modulus"];
                frequency_independent_loss = d["findependent_loss"];
                frequency_dependent_loss = d["fdependent_loss"];
                tension = d["tension"];
                
                m_geometry_changed = true;
                m_material_changed = true;
                m_position_changed = true;
                update_queue.set();
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
            std::unique_lock<std::mutex> lock(m_update_mutex);
            m_synthesizer.set_sample_rate(samplerate());
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
        // Note: Geometry update recalculates material coefficients, but we don't clear
        // m_material_changed here in case new material parameters were set after geometry
    }
    
    void update_material_if_needed() {
        if (!m_material_changed) return;
        
        m_synthesizer.update_material_parameters(m_parameters);
        m_material_changed = false;
    }
    
    void update_position_if_needed() {
        if (!m_position_changed) return;
        
        m_synthesizer.update_position_parameters(m_force_position, m_readout_position);
        m_position_changed = false;
    }

    StringModalSynthesizer m_synthesizer;
    std::mutex m_update_mutex;
    
    // Parameter storage
    StringParameters m_parameters;
    double m_force_position = 0.5;
    double m_readout_position = 0.5;
    
    // Update flags (handled externally, not by synthesizer)
    bool m_material_changed = true;
    bool m_geometry_changed = true;
    bool m_position_changed = true;
    bool m_initialized = false;
};

MIN_EXTERNAL(string_tilde);