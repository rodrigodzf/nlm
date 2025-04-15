#include "c74_min.h"
#include "c74_min_attribute.h"
#include <cmath>
#include <Eigen/Dense>
#include <matioCpp/matioCpp.h>
#include <vk_utils/TimeIntegrators.h>

using namespace c74::min;

class modal_resonator : public object<modal_resonator> {
public:
    // Constructor
    modal_resonator(const atoms& args = {}) {
        // Initialize any variables or parameters here
        calculate_coefficients();
    }

    // Inlets and outlets
    inlet<>  input    { this, "(signal) Input to the modal resonator" };
    outlet<> output   { this, "(signal) Output from the modal resonator" };

    // Attributes
    attribute<number> sample_rate { this, "sample_rate", 44100.0,
        description {"Sample rate in Hz"},
        range { 16000.0, 192000.0 }
    };

    attribute<number> excitation_position { this, "excitation_position", 0.2,
        description {"Excitation position (0-1)"},
        range { 0.0, 1.0 }
    };

    attribute<number> readout_position { this, "readout_position", 0.3,
        description {"Readout position (0-1)"},
        range { 0.0, 1.0 }
    };

    attribute<number> poisson_ratio { this, "poisson_ratio", 0.35,
        description {"Poisson's ratio"},
        range { 0.0, 1.0 }
    };

    attribute<number> thickness { this, "thickness", 5e-4,
        description {"Thickness in meters"},
        range { 1e-4, 1e-3 }
    };
    
    attribute<number> lx { this, "lx", 0.2,
        description {"Length in meters"},
        range { 0.0, 1.0 }
    };

    attribute<number> ly { this, "ly", 0.3,
        description {"Width in meters"},
        range { 0.0, 1.0 }
    };

    attribute<number> density { this, "density", 7800,
        description {"Density in kg/m^3"},
        range { 1000.0, 10000.0 }
    };

    attribute<number> youngs_modulus { this, "youngs_modulus", 2000,
        description {"Young's modulus in GPa"},
        range { 100.0, 10000.0 }
    };

    attribute<number> frequency_independent_loss { this, "frequency_independent_loss", 4.2e-2,
        description {"Frequency independent loss"},
        range { 0.0, 1.0 }
    };

    attribute<number> frequency_dependent_loss { this, "frequency_dependent_loss", 2.3e-3,
        description {"Frequency dependent loss"},
        range { 0.0, 1.0 }
    };
    

    // Methods using direct binding approach
    message<> dspsetup { this, "dspsetup",
        MIN_FUNCTION {
            setup(args, 64); // Default vector size if not provided
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
    message<> set_excitation_position { this, "set_excitation_position",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            excitation_position = args[0];
            return {};
        }
    };

    message<> set_readout_position { this, "set_readout_position",
        [this](const c74::min::atoms& args, const int inlet) -> c74::min::atoms {
            readout_position = args[0];
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
            
            matioCpp::MultiDimensionalArray<double> H_mat = file.read("H").asMultiDimensionalArray<double>();
            matioCpp::Vector<double> lambda_mu_mat = file.read("lambda_mu").asVector<double>();
            
            // convert to Eigen::Matrix
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(H_mat.data(), H_mat.dimensions()[0], H_mat.dimensions()[1]);
            Eigen::VectorXd lambda_mu = matioCpp::to_eigen(lambda_mu_mat);
            
            return {};
        }
    };
    

    // DSP perform method (for Max/MSP audio processing)
    void operator()(audio_bundle input, audio_bundle output) {
        const int vector_size = input.frame_count();
        
        // Store current values to check for changes
        static double last_frequency = 0.0;
        static double last_decay = 0.0;
        
        // Check if parameters have changed
        // if (!m_initialized || 
        //     last_frequency != static_cast<double>(frequency) || 
        //     last_decay != static_cast<double>(decay)) {
            
        //     calculate_coefficients();
            
        //     // Update stored values
        //     last_frequency = static_cast<double>(frequency);
        //     last_decay = static_cast<double>(decay);
        // }
        
        auto in = input.samples(0);   // Get input samples
        auto out = output.samples(0); // Get output buffer
        
        // Process each sample
        for (int i = 0; i < vector_size; ++i) {
            // Apply biquad filter equation: y[n] = b0*x[n] - a1*y[n-1] - a2*y[n-2]
            // double y0 = m_b0 * in[i] - m_a1 * m_y1 - m_a2 * m_y2;
            
            // Store output with gain applied
            // out[i] = y0 * gain;
            
            // Update delay registers
            // m_y2 = m_y1;
            // m_y1 = y0;
            out[i] = in[i];
        }
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
    bool m_initialized = false;
};

MIN_EXTERNAL(modal_resonator);