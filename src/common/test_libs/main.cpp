#include <iostream>
#include <cxxopts.hpp>
#include <vector>
#include <fstream>
#include <Eigen/Dense>

#include "vk_utils/Parameters.h"
#include "vk_utils/FTM.h"
#include "vk_utils/Excitations.h"
#include "vk_utils/vk_utils.h"
#include "vk_utils/TimeIntegrators.h"

// Function to write output to a file
template <typename T>
void write_vector_to_file(const std::string& filename, const Eigen::Matrix<T, Eigen::Dynamic, 1>& vector) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    
    for (int i = 0; i < vector.size(); ++i) {
        outfile << vector(i) << std::endl;
    }
    
    outfile.close();
}

int main(int argc, char* argv[])
{
    // Parse command line arguments
    cxxopts::Options options("test_libs", "Test the vk_utils library");
    
    options.add_options()
        ("i,input", "Input file path", cxxopts::value<std::string>()->default_value("benchmark_input_010.mat"))
        ("o,output", "Output file path", cxxopts::value<std::string>()->default_value("output.csv"))
        ("t,use_tm", "Use tension modulation", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage");
    
    auto result = options.parse(argc, argv);
    
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }
    
    std::string output_file = result["output"].as<std::string>();
    
    std::cout << "Testing vk_utils library..." << std::endl;
    
    // Choose precision based on command line argument
    // 1. Define the setup with float precision
    const int n_modes = 50;
    const int n_steps = 44100;
    const float sample_rate = 44100.0f;
    const float dt = 1.0f / sample_rate;
    const float excitation_position = 0.2f;
    const float readout_position = 0.5f;
    const float initial_deflection = 0.03f;
    const int n_gridpoints = 101;
    
    StringParameters string_params;
    
    // Calculate eigenvalues
    Eigen::Matrix<float, Eigen::Dynamic, 1> lambda_mu = 
        ftm::string_eigenvalues<float>(n_modes, string_params.length);
    
    // 2. Create the excitation (initial conditions)
    Eigen::Matrix<float, Eigen::Dynamic, 1> u0_modal = excitations::create_pluck_modal<float>(
        lambda_mu,
        excitation_position,
        initial_deflection,
        string_params.length
    );
    
    // 3. Get damping and stiffness terms
    Eigen::Matrix<float, Eigen::Dynamic, 1> gamma2_mu = damping_term<float>(
        string_params,
        lambda_mu
    );
    
    Eigen::Matrix<float, Eigen::Dynamic, 1> omega_mu_squared = stiffness_term<float>(
        string_params,
        lambda_mu
    );
    
    // Create a zero excitation signal (since we're using initial conditions)
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> zero_excitation = 
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>::Zero(n_steps, n_modes);
    
    // Set initial displacement as first step in excitation to use solve_sv_excitation
    zero_excitation.row(0) = u0_modal.transpose();
    
    // Solve the modal system using the Störmer-Verlet time integrator
    auto modal_sol = time_integrators::solve_sv_excitation<float>(
        gamma2_mu,
        omega_mu_squared,
        zero_excitation,
        dt,
        nullptr  // No nonlinear function
    );
    
    // 4. Get the readout weights
    Eigen::Matrix<int, Eigen::Dynamic, 1> mu(n_modes);
    for (int i = 0; i < n_modes; ++i) {
        mu(i) = i + 1;  // Mode indices starting from 1
    }
    
    Eigen::Matrix<float, Eigen::Dynamic, 1> readout_weights = 
        ftm::evaluate_string_eigenfunctions<float>(mu, readout_position, string_params.length);
    
    // Compute readout signal
    Eigen::Matrix<float, Eigen::Dynamic, 1> u_readout(n_steps);
    for (int i = 0; i < n_steps; ++i) {
        u_readout(i) = readout_weights.dot(modal_sol.row(i));
    }
    
    // Write output to file
    write_vector_to_file(output_file, u_readout);
    
    std::cout << "Simulation completed with single precision." << std::endl;
    std::cout << "Output written to: " << output_file << std::endl;
  

    return 0;
}