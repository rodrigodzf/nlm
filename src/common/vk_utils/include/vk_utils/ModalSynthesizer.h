#pragma once

#include <Eigen/Dense>
#include <vk_utils/Parameters.h>
#include <vk_utils/FTM.h>
#include <vk_utils/LinearInterpolator.h>
#include <vk_utils/ParallelFilter.h>
#include <vk_utils/vk_utils.h>
#include <vk_utils/CouplingAnalytical.h>

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * Base class for modal synthesis objects with consolidated update patterns.
 * 
 * This class encapsulates the three types of updates needed for modal synthesis:
 * 1. Material parameter updates (affecting damping/stiffness)
 * 2. Geometry updates (affecting eigenvalues/indices) 
 * 3. Position updates (affecting force/readout weights)
 * 
 * Thread safety must be handled externally by the caller.
 */
template<typename ParameterType>
class ModalSynthesizer {
public:
    explicit ModalSynthesizer(int n_modes = 32) 
        : m_n_modes(n_modes) {
        resize_vectors();
    }

    virtual ~ModalSynthesizer() = default;

    // Set sample rate (should be called once during setup)
    void set_sample_rate(double sample_rate) {
        m_dt = 1.0 / sample_rate;
    }

    // Separate update functions for each type of parameter change
    void update_material_parameters(const ParameterType& params) {
        m_parameters = params;
        calculate_material_coefficients();
        m_parallel_filter.set_coefficients(m_gamma2_mu, m_omega_mu_squared, m_dt);
    }

    void update_geometry_parameters(const ParameterType& params) {
        m_parameters = params;
        calculate_geometry();
        // Material coefficients depend on eigenvalues, so recalculate them
        calculate_material_coefficients();
        m_parallel_filter.set_coefficients(m_gamma2_mu, m_omega_mu_squared, m_dt);
    }

    void update_position_parameters(double force_pos, double readout_pos) {
        m_force_position = force_pos;
        m_readout_position = readout_pos;
        calculate_position_weights();
    }

    // Audio processing function
    double process_sample(double input) {
        // Calculate force input
        m_force_input = input * m_force_weights;

        // Update nonlinearity (model-specific)
        update_nonlinearity();
        
        // Process input through parallel filter
        m_parallel_filter(m_force_input);

        // Get interpolated readout weights and compute output
        m_readout_weights_lerp.process(m_current_readout_weights);
        return m_current_readout_weights.dot(m_parallel_filter.get_q());
    }

    void set_n_modes(int n_modes) {
        if (m_n_modes != n_modes) {
            m_n_modes = n_modes;
            resize_vectors();
        }
    }

    void reset_filter() {
        m_parallel_filter.reset();
    }

    // Getters
    int get_n_modes() const { return m_n_modes; }
    const ParameterType& get_parameters() const { return m_parameters; }
    double get_force_position() const { return m_force_position; }
    double get_readout_position() const { return m_readout_position; }

protected:
    // Pure virtual functions to be implemented by derived classes
    virtual void calculate_geometry() = 0;
    virtual void calculate_material_coefficients() = 0;
    virtual void calculate_position_weights() = 0;
    virtual void update_nonlinearity() = 0;

    // Common member variables
    int m_n_modes;
    double m_dt = 1.0 / 44100.0;
    
    // Parameters
    ParameterType m_parameters;
    double m_force_position = 0.5;
    double m_readout_position = 0.5;
    
    // Core synthesis components
    ParallelFilter<double> m_parallel_filter;
    VectorInterpolator<double> m_readout_weights_lerp;
    
    // Modal coefficients (computed from parameters)
    Vector m_lambda_mu;
    Vector m_gamma2_mu;
    Vector m_omega_mu_squared;
    Vector m_selected_indices;
    
    // Force and readout weights
    Vector m_force_weights;
    Vector m_readout_weights;
    Vector m_current_readout_weights;
    Vector m_force_input;

private:
    void resize_vectors() {
        m_readout_weights_lerp.resize(m_n_modes, 441);
        m_current_readout_weights.resize(m_n_modes);
        m_parallel_filter.resize(m_n_modes);
        m_force_input.resize(m_n_modes);
        m_lambda_mu.resize(m_n_modes);
        m_gamma2_mu.resize(m_n_modes);
        m_omega_mu_squared.resize(m_n_modes);
        m_selected_indices.resize(m_n_modes);
        m_force_weights.resize(m_n_modes);
        m_readout_weights.resize(m_n_modes);
    }
};

/**
 * String modal synthesizer implementation
 */
class StringModalSynthesizer : public ModalSynthesizer<StringParameters> {
public:
    explicit StringModalSynthesizer(int n_modes = 32) 
        : ModalSynthesizer<StringParameters>(n_modes) {}

protected:
    void calculate_geometry() override {
        // Calculate modal indices and eigenvalues for string
        m_selected_indices = Eigen::VectorXd::LinSpaced(m_n_modes, 1, m_n_modes);
        m_lambda_mu = ftm::string_eigenvalues(m_n_modes, m_parameters.length);
    }

    void calculate_material_coefficients() override {
        // Calculate damping and stiffness terms
        m_omega_mu_squared = stiffness_term<double>(m_parameters, m_lambda_mu);
        m_gamma2_mu = damping_term<double>(m_parameters, m_lambda_mu);
        
        // Calculate nonlinear coupling coefficients  
        double string_norm = m_parameters.length * 0.5;
        double string_tau = (m_parameters.E * m_parameters.A) / 
                           (m_parameters.length * 2.0);
        string_tau = string_tau / m_parameters.density() / string_norm;
        m_string_tau_with_norms = string_tau * m_lambda_mu;
    }

    void calculate_position_weights() override {
        // Calculate force weights
        double string_norm = m_parameters.length * 0.5;
        m_force_weights = ftm::evaluate_string_eigenfunctions(
            m_selected_indices,
            m_force_position * m_parameters.length,
            m_parameters.length
        ) / (string_norm * m_parameters.density());
        
        // Calculate readout weights with interpolation
        Vector new_readout_weights = ftm::evaluate_string_eigenfunctions(
            m_selected_indices,
            m_readout_position * m_parameters.length,
            m_parameters.length
        );
        
        if (m_readout_weights_lerp.isFinished()) {
            if (m_readout_weights.size() == 0) {
                m_readout_weights = new_readout_weights;
                m_readout_weights_lerp.setValue(new_readout_weights);
            }
            m_readout_weights_lerp.setTarget(new_readout_weights);
        }
    }

    void update_nonlinearity() override {
        // String uses lambda_mu based nonlinearity
        m_parallel_filter.update_nonlinearity(m_lambda_mu, m_string_tau_with_norms);
    }

private:
    Vector m_string_tau_with_norms;
};

/**
 * Plate modal synthesizer implementation
 */
class PlateModalSynthesizer : public ModalSynthesizer<PlateParameters> {
public:
    explicit PlateModalSynthesizer(int n_modes = 64, const std::string& model_type = "berger") 
        : ModalSynthesizer<PlateParameters>(n_modes), m_model_type(model_type), m_n_psi(10) {}

    void set_model_type(const std::string& model_type) {
        m_model_type = model_type;
    }

    void set_n_psi(int n_psi) {
        m_n_psi = n_psi;
    }

    void set_precomputed_data(const Matrix& H_original, const Vector& lambda_mu, 
                             const Vector& selected_indices_x, const Vector& selected_indices_y) {
        m_H_original = H_original;
        m_lambda_mu = lambda_mu;
        m_selected_indices_x = selected_indices_x;
        m_selected_indices_y = selected_indices_y;
        m_data_loaded = true;
        
        // Resize vectors to match precomputed data
        m_n_modes = lambda_mu.size();
        resize_vectors();
    }

protected:
    void calculate_geometry() override {
        if (m_data_loaded) {
            // Use precomputed data from .mat file
            return;
        }
        
        // Calculate eigenvalues matrix for berger model
        int n_max_modes_x = 20;
        int n_max_modes_y = 20;
        
        Matrix lambda_mu_2d = calculate_plate_eigenvalues<double>(
            n_max_modes_x, n_max_modes_y, m_parameters.l1, m_parameters.l2);
        
        // Allocate vectors for selected modes
        m_lambda_mu.resize(m_n_modes);
        Eigen::MatrixX<int> selected_indices(m_n_modes, 2);
        
        // Select modes and eigenvalues
        select_modes_and_eigenvalues<double>(lambda_mu_2d, m_n_modes, m_lambda_mu, selected_indices);
        
        // Extract indices into separate x and y vectors
        m_selected_indices_x.resize(m_n_modes);
        m_selected_indices_y.resize(m_n_modes);
        
        for (int i = 0; i < m_n_modes; ++i) {
            m_selected_indices_x(i) = selected_indices(i, 0);
            m_selected_indices_y(i) = selected_indices(i, 1);
        }
        
        // For von Karman model, compute coupling matrix
        if (m_model_type == "vk") {
            compute_coupling_matrix(
                m_n_psi,
                m_n_modes,
                m_parameters.l1,
                m_parameters.l2,
                m_selected_indices_x.cast<int>(),
                m_selected_indices_y.cast<int>(),
                m_H_original
            );
            
            // Reshape the H tensor to (n_psi * n_phi, n_phi)
            m_H_original = Eigen::Map<const Matrix>(m_H_original.transpose().data(), m_n_psi * m_n_modes, m_n_modes);
        }
    }

    void calculate_material_coefficients() override {
        // Calculate damping and stiffness terms
        m_omega_mu_squared = stiffness_term<double>(m_parameters, m_lambda_mu);
        m_gamma2_mu = damping_term<double>(m_parameters, m_lambda_mu);
        
        double plate_norm = m_parameters.l1 * m_parameters.l2 * 0.25;
        
        // Model-specific calculations
        if (m_model_type == "berger") {
            double plate_tau = (m_parameters.E * m_parameters.h) / (
                2.0 * m_parameters.l1 * m_parameters.l2 * (1.0 - m_parameters.nu * m_parameters.nu)
            );
            plate_tau = plate_tau / m_parameters.density() / plate_norm;
            m_plate_tau_with_norms = plate_tau * m_lambda_mu;
        } else if (m_model_type == "vk") {
            double scale = m_parameters.E / (2.0 * m_parameters.rho * plate_norm);
            m_H_scaled = m_H_original * std::sqrt(scale);
        }
    }

    void calculate_position_weights() override {
        // Calculate force weights
        m_force_weights = ftm::evaluate_rectangular_eigenfunctions(
            m_selected_indices_x, m_selected_indices_y, 
            m_force_position_x * m_parameters.l1, m_force_position_y * m_parameters.l2, 
            m_parameters.l1, m_parameters.l2
        ) / m_parameters.density();

        // Calculate readout weights with interpolation
        Vector new_readout_weights = ftm::evaluate_rectangular_eigenfunctions(
            m_selected_indices_x, m_selected_indices_y, 
            m_readout_position_x * m_parameters.l1,
            m_readout_position_y * m_parameters.l2,
            m_parameters.l1, m_parameters.l2
        ) / (m_parameters.l1 * m_parameters.l2 * 0.25);
        
        if (m_readout_weights_lerp.isFinished()) {
            if (m_readout_weights.size() == 0) {
                m_readout_weights = new_readout_weights;
                m_readout_weights_lerp.setValue(new_readout_weights);
            }
            m_readout_weights_lerp.setTarget(new_readout_weights);
        }
    }

    void update_nonlinearity() override {
        if (m_model_type == "vk") {
            m_parallel_filter.update_nonlinearity(m_H_scaled, m_n_psi, m_n_modes);
        } else if (m_model_type == "berger") {
            m_parallel_filter.update_nonlinearity(m_lambda_mu, m_plate_tau_with_norms);
        }
    }

public:
    // Override the base class update_position_parameters to handle 2D positions
    void update_position_parameters(double force_pos, double readout_pos) {
        // For plates, we use the internal 2D positions, not the parameters
        calculate_position_weights();
    }
    // 2D position setters
    void set_force_position(double x, double y) {
        m_force_position_x = x;
        m_force_position_y = y;
    }

    void set_readout_position(double x, double y) {
        m_readout_position_x = x;
        m_readout_position_y = y;
    }

protected:
    // Plate-specific data (accessible to derived classes)
    Vector m_selected_indices_x;
    Vector m_selected_indices_y;
    Matrix m_H_original;
    Matrix m_H_scaled;
    Vector m_plate_tau_with_norms;

private:
    void resize_vectors() {
        m_readout_weights_lerp.resize(m_n_modes, 441);
        m_current_readout_weights.resize(m_n_modes);
        m_parallel_filter.resize(m_n_modes);
        m_force_input.resize(m_n_modes);
        m_lambda_mu.resize(m_n_modes);
        m_gamma2_mu.resize(m_n_modes);
        m_omega_mu_squared.resize(m_n_modes);
        m_force_weights.resize(m_n_modes);
        m_readout_weights.resize(m_n_modes);
    }

private:
    std::string m_model_type;
    int m_n_psi;
    bool m_data_loaded = false;
    
    // 2D positions
    double m_force_position_x = 0.5;
    double m_force_position_y = 0.5;
    double m_readout_position_x = 0.5;
    double m_readout_position_y = 0.5;
};

/**
 * Multi-channel string modal synthesizer implementation
 */
class MultiChannelStringModalSynthesizer : public StringModalSynthesizer {
public:
    explicit MultiChannelStringModalSynthesizer(int n_modes = 32, int n_inputs = 1, int n_outputs = 2) 
        : StringModalSynthesizer(n_modes), m_n_inputs(n_inputs), m_n_outputs(n_outputs) {
        resize_multichannel_vectors();
    }

    void set_channel_counts(int n_inputs, int n_outputs) {
        if (m_n_inputs != n_inputs || m_n_outputs != n_outputs) {
            m_n_inputs = n_inputs;
            m_n_outputs = n_outputs;
            resize_multichannel_vectors();
        }
    }

    void set_force_positions(const Vector& positions) {
        if (positions.size() == m_n_inputs) {
            m_force_positions = positions;
        }
    }

    void set_readout_positions(const Vector& positions) {
        if (positions.size() == m_n_outputs) {
            m_readout_positions = positions;
        }
    }

    void process_frame_multichannel(const Vector& inputs) {
        // Safety check - ensure vectors are properly sized
        if (m_force_weights_multichannel.cols() != m_n_inputs) {
            return;
        }
        
        // Calculate force input for all channels
        for (int ch = 0; ch < m_n_inputs; ++ch) {
            m_force_input_multichannel.col(ch) = inputs[ch] * m_force_weights_multichannel.col(ch);
        }

        // Update nonlinearity (same for all channels)
        update_nonlinearity();
        
        // Process summed input through parallel filter
        m_parallel_filter(m_force_input_multichannel.rowwise().sum());

        // Get interpolated readout weights (advance interpolator once per frame)
        m_readout_weights_lerp_multichannel.process(m_current_readout_weights_multichannel);
    }

    double get_output_multichannel(int output_channel) {
        // Safety check - ensure vectors are properly sized
        if (m_current_readout_weights_multichannel.cols() != m_n_outputs ||
            output_channel >= m_n_outputs) {
            return 0.0;
        }
        
        // Return output for specific channel
        return m_current_readout_weights_multichannel.col(output_channel).dot(m_parallel_filter.get_q());
    }

    void update_multichannel_position_parameters() {
        calculate_multichannel_position_weights();
    }

protected:
    void calculate_multichannel_position_weights() {
        // Calculate force weights for all input channels
        double string_norm = m_parameters.length * 0.5;
        m_force_weights_multichannel.resize(m_n_modes, m_n_inputs);
        
        for (int ch = 0; ch < m_n_inputs; ++ch) {
            m_force_weights_multichannel.col(ch) = ftm::evaluate_string_eigenfunctions(
                m_selected_indices,
                m_force_positions[ch] * m_parameters.length,
                m_parameters.length
            ) / (string_norm * m_parameters.density());
        }
        
        // Calculate readout weights for all output channels with interpolation
        Matrix new_readout_weights(m_n_modes, m_n_outputs);
        for (int ch = 0; ch < m_n_outputs; ++ch) {
            new_readout_weights.col(ch) = ftm::evaluate_string_eigenfunctions(
                m_selected_indices,
                m_readout_positions[ch] * m_parameters.length,
                m_parameters.length
            );
        }
        
        if (m_readout_weights_lerp_multichannel.isFinished()) {
            if (m_readout_weights_multichannel.size() == 0) {
                m_readout_weights_multichannel = new_readout_weights;
                m_readout_weights_lerp_multichannel.setValue(new_readout_weights);
            }
            m_readout_weights_lerp_multichannel.setTarget(new_readout_weights);
        }
    }

private:
    void resize_multichannel_vectors() {
        m_readout_weights_lerp_multichannel.resize(m_n_modes, m_n_outputs, 441);
        m_current_readout_weights_multichannel.resize(m_n_modes, m_n_outputs);
        m_force_input_multichannel.resize(m_n_modes, m_n_inputs);
        m_force_weights_multichannel.resize(m_n_modes, m_n_inputs);
        m_readout_weights_multichannel.resize(m_n_modes, m_n_outputs);
        m_force_positions.resize(m_n_inputs);
        m_readout_positions.resize(m_n_outputs);
        m_force_positions.setConstant(0.5);
        m_readout_positions.setConstant(0.5);
    }

private:
    int m_n_inputs;
    int m_n_outputs;
    
    // Multi-channel position vectors
    Vector m_force_positions;
    Vector m_readout_positions;
    
    // Multi-channel weight matrices
    Matrix m_force_weights_multichannel;
    Matrix m_readout_weights_multichannel;
    Matrix m_current_readout_weights_multichannel;
    Matrix m_force_input_multichannel;
    
    MatrixInterpolator<double> m_readout_weights_lerp_multichannel;
};

/**
 * Multi-channel plate modal synthesizer implementation
 */
class MultiChannelPlateModalSynthesizer : public PlateModalSynthesizer {
public:
    explicit MultiChannelPlateModalSynthesizer(int n_modes = 64, int n_inputs = 1, int n_outputs = 2, const std::string& model_type = "berger") 
        : PlateModalSynthesizer(n_modes, model_type), m_n_inputs(n_inputs), m_n_outputs(n_outputs) {
        resize_multichannel_vectors();
    }

    void set_channel_counts(int n_inputs, int n_outputs) {
        if (m_n_inputs != n_inputs || m_n_outputs != n_outputs) {
            m_n_inputs = n_inputs;
            m_n_outputs = n_outputs;
            resize_multichannel_vectors();
        }
    }

    void set_force_positions(const Vector& positions_x, const Vector& positions_y) {
        if (positions_x.size() == m_n_inputs && positions_y.size() == m_n_inputs) {
            m_force_positions_x = positions_x;
            m_force_positions_y = positions_y;
        }
    }

    void set_readout_positions(const Vector& positions_x, const Vector& positions_y) {
        if (positions_x.size() == m_n_outputs && positions_y.size() == m_n_outputs) {
            m_readout_positions_x = positions_x;
            m_readout_positions_y = positions_y;
        }
    }

    void process_frame_multichannel(const Vector& inputs) {
        // Safety check - ensure vectors are properly sized
        if (m_force_weights_multichannel.cols() != m_n_inputs) {
            return;
        }
        
        // Calculate force input for all channels
        for (int ch = 0; ch < m_n_inputs; ++ch) {
            m_force_input_multichannel.col(ch) = inputs[ch] * m_force_weights_multichannel.col(ch);
        }

        // Update nonlinearity (same for all channels)
        update_nonlinearity();
        
        // Process summed input through parallel filter
        m_parallel_filter(m_force_input_multichannel.rowwise().sum());

        // Get interpolated readout weights (advance interpolator once per frame)
        m_readout_weights_lerp_multichannel.process(m_current_readout_weights_multichannel);
    }

    double get_output_multichannel(int output_channel) {
        // Safety check - ensure vectors are properly sized
        if (m_current_readout_weights_multichannel.cols() != m_n_outputs ||
            output_channel >= m_n_outputs) {
            return 0.0;
        }
        
        // Return output for specific channel
        return m_current_readout_weights_multichannel.col(output_channel).dot(m_parallel_filter.get_q());
    }

    void update_multichannel_position_parameters() {
        calculate_multichannel_position_weights();
    }

protected:
    void calculate_multichannel_position_weights() {
        // Calculate force weights for all input channels
        m_force_weights_multichannel.resize(m_n_modes, m_n_inputs);
        
        for (int ch = 0; ch < m_n_inputs; ++ch) {
            Vector force_weights = ftm::evaluate_rectangular_eigenfunctions(
                m_selected_indices_x, m_selected_indices_y,
                m_force_positions_x[ch] * m_parameters.l1, m_force_positions_y[ch] * m_parameters.l2,
                m_parameters.l1, m_parameters.l2
            ) / m_parameters.density();
            m_force_weights_multichannel.col(ch) = force_weights;
        }
        
        // Calculate readout weights for all output channels with interpolation
        Matrix new_readout_weights(m_n_modes, m_n_outputs);
        for (int ch = 0; ch < m_n_outputs; ++ch) {
            Vector readout_weights = ftm::evaluate_rectangular_eigenfunctions(
                m_selected_indices_x, m_selected_indices_y,
                m_readout_positions_x[ch] * m_parameters.l1, m_readout_positions_y[ch] * m_parameters.l2,
                m_parameters.l1, m_parameters.l2
            ) / (m_parameters.l1 * m_parameters.l2 * 0.25);
            new_readout_weights.col(ch) = readout_weights;
        }
        
        if (m_readout_weights_lerp_multichannel.isFinished()) {
            if (m_readout_weights_multichannel.size() == 0) {
                m_readout_weights_multichannel = new_readout_weights;
                m_readout_weights_lerp_multichannel.setValue(new_readout_weights);
            }
            m_readout_weights_lerp_multichannel.setTarget(new_readout_weights);
        }
    }

private:
    void resize_multichannel_vectors() {
        m_readout_weights_lerp_multichannel.resize(m_n_modes, m_n_outputs, 441);
        m_current_readout_weights_multichannel.resize(m_n_modes, m_n_outputs);
        m_force_input_multichannel.resize(m_n_modes, m_n_inputs);
        m_force_weights_multichannel.resize(m_n_modes, m_n_inputs);
        m_readout_weights_multichannel.resize(m_n_modes, m_n_outputs);
        m_force_positions_x.resize(m_n_inputs);
        m_force_positions_y.resize(m_n_inputs);
        m_readout_positions_x.resize(m_n_outputs);
        m_readout_positions_y.resize(m_n_outputs);
        m_force_positions_x.setConstant(0.5);
        m_force_positions_y.setConstant(0.5);
        m_readout_positions_x.setConstant(0.5);
        m_readout_positions_y.setConstant(0.5);
    }

private:
    int m_n_inputs;
    int m_n_outputs;
    
    // Multi-channel 2D position vectors
    Vector m_force_positions_x;
    Vector m_force_positions_y;
    Vector m_readout_positions_x;
    Vector m_readout_positions_y;
    
    // Multi-channel weight matrices
    Matrix m_force_weights_multichannel;
    Matrix m_readout_weights_multichannel;
    Matrix m_current_readout_weights_multichannel;
    Matrix m_force_input_multichannel;
    
    MatrixInterpolator<double> m_readout_weights_lerp_multichannel;
};