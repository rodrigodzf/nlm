#include "c74_min.h"
#include "c74_min_api.h"
#include "c74_min_attribute.h"
#include "c74_min_queue.h"
#include <Eigen/Dense>
#include <mutex>
#include <sstream>
#include <vk_utils/ModalSynthesizer.h>
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

    int m_num_inputs;
    int m_num_outputs;
    MultiChannelStringModalSynthesizer m_synthesizer;

    mcs_string_tilde(const atoms& args = {}) 
        : m_num_inputs(args.size() > 0 ? static_cast<int>(args[0]) : 1),
          m_num_outputs(args.size() > 1 ? static_cast<int>(args[1]) : 2),
          m_synthesizer(32, m_num_inputs, m_num_outputs)
    {
        
        // Initialize positions
        Vector force_positions = Vector::Constant(m_num_inputs, 0.5);
        Vector readout_positions = Vector::Constant(m_num_outputs, 0.5);
        m_synthesizer.set_force_positions(force_positions);
        m_synthesizer.set_readout_positions(readout_positions);
    }

    inlet<>  input    { this, "(multichannelsignal) Input to the modal resonator" };
    outlet<> output   { this, "(multichannelsignal) Output from the modal resonator", "multichannelsignal" };

    // Single update queue for all parameter changes
    queue<> update_queue {this, MIN_FUNCTION {
        update_parameters();
        return {};
    }};

    attribute<number> n_modes { this, "n_modes", 32,
        description {"Number of modes"},
        range { 1, 1000 },
        setter { MIN_FUNCTION {
            m_synthesizer.set_n_modes(static_cast<int>(args[0]));
            m_geometry_changed = true;
            update_queue.set();
            return { args[0] };
        }}
    };

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
    
    message<> maxclass_setup { this, "maxclass_setup",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            cout << "mcs.modal.string~ - " << VERSION << " - 2025 - Rodrigo Diaz" << endl;
            
            c74::max::t_class* c = args[0];
            c74::max::class_addmethod(c, (c74::max::method)simplemc_multichanneloutputs, "multichanneloutputs", c74::max::A_CANT, 0);
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
    
    message<> info { this, "info",
        MIN_FUNCTION {
            std::stringstream ss;
            ss << "mcs.modal.string~ parameters:" << '\n'
               << "Length: " << lx << " m" << '\n'
               << "Cross sectional area: " << cross_sectional_area << " mm²" << '\n'
               << "Density: " << density << " kg/m³" << '\n'
               << "Young's modulus: " << youngs_modulus << " GPa" << '\n'
               << "Moment of inertia: " << moment_of_inertia << " mm⁴" << '\n'
               << "Tension: " << tension << " N" << '\n'
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
            if (args.size() != m_num_inputs) {
                cerr << "force_position: must provide " << m_num_inputs << " positions" << endl;
                return {};
            }
            
            Vector positions(m_num_inputs);
            for (int i = 0; i < m_num_inputs; ++i) {
                // Clamp values between 0 and 1
                positions[i] = std::clamp(double(args[i]), 0.0, 1.0);
            }
            
            m_synthesizer.set_force_positions(positions);
            m_position_changed = true;
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

            Vector positions(m_num_outputs);
            for (int i = 0; i < m_num_outputs; ++i) {
                // Clamp values between 0 and 1
                positions[i] = std::clamp(double(args[i]), 0.0, 1.0);
            }

            m_synthesizer.set_readout_positions(positions);
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
    StringParameters m_parameters;
    
    // Update flags (handled externally, not by synthesizer)
    bool m_material_changed = true;
    bool m_geometry_changed = true;
    bool m_position_changed = true;
    bool m_initialized = false;
};

long simplemc_multichanneloutputs(c74::max::t_object* x, long index, long count) {
    minwrap<mcs_string_tilde>* ob = (minwrap<mcs_string_tilde>*)(x);
    return ob->m_min_object.m_num_outputs;
}

MIN_EXTERNAL(mcs_string_tilde);