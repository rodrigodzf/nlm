#include "bhmat.h"
#include <matioCpp/matioCpp.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <algorithm>  // For std::sort, std::min

#ifdef EIGEN_USE_ARPACK
#include <unsupported/Eigen/ArpackSupport>
#endif

// Trapezoid integration function (similar to scipy.integrate.trapezoid)
template <typename T>
T trapezoid(const Eigen::Matrix<T, Eigen::Dynamic, 1>& y, T dx) {
    int n = y.size();
    if (n < 2) return T(0);
    
    T result = T(0);
    // First point
    result += y(0) / 2.0;
    // Middle points (weight = 1)
    for (int i = 1; i < n - 1; i++) {
        result += y(i);
    }
    // Last point
    result += y(n - 1) / 2.0;
    
    return result * dx;
}

// Double trapezoid integration on 2D matrix
template <typename T>
T double_trapezoid(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& f, T dx, T dy) {
    int rows = f.rows();
    int cols = f.cols();
    
    // First apply trapezoid along each column
    Eigen::Matrix<T, Eigen::Dynamic, 1> column_integrals(cols);
    for (int j = 0; j < cols; j++) {
        column_integrals(j) = trapezoid<T>(f.col(j), dy);
    }
    
    // Then apply trapezoid to the result
    return trapezoid<T>(column_integrals, dx);
}

// Double trapezoid integration on flattened array
template <typename T>
T double_trapezoid_flat(const Eigen::Matrix<T, Eigen::Dynamic, 1>& f, T dx, T dy, int Ny, int Nx) {
    // Reshape flattened array into 2D matrix (Fortran order - column major)
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> F(Ny, Nx);
    
    for (int j = 0; j < Nx; j++) {
        for (int i = 0; i < Ny; i++) {
            F(i, j) = f(i + j * Ny);
        }
    }
    
    return double_trapezoid<T>(F, dx, dy);
}

int main(int argc, char *argv[]) {
    // Parameters for the biharmonic matrix
    const int Nx = 50;
    const int Ny = 75;
    Eigen::Vector2i Nxy(Nx, Ny);
    
    const double h = 0.004;     // Grid spacing
    const double Lz = 5.0e-4; // Plate thickness
    const double E = 2.0e12;  // Young's modulus (steel)
    const double nu = 0.3;    // Poisson's ratio (steel)
    
    // Boundary conditions (4x2 matrix)
    // Format: [K0y, R0y; Kx0, Rx0; KLy, RLy; KxL, RxL]
    Eigen::Matrix<double, 4, 2> BCs;
    // Simply supported edges (pinned)
    BCs << 1.0e15, 0.0,  // K0y, R0y
           1.0e15, 0.0,  // Kx0, Rx0
           1.0e15, 0.0,  // KLy, RLy
           1.0e15, 0.0;  // KxL, RxL
    
    std::cout << "Generating biharmonic matrix..." << std::endl;
    
    // Generate the biharmonic matrix
    Eigen::SparseMatrix<double> biharmonicMatrix = bhmat<double>(BCs, Nxy, h, Lz, E, nu);
    
    std::cout << "Matrix size: " << biharmonicMatrix.rows() << " x " << biharmonicMatrix.cols() << std::endl;
    std::cout << "Number of non-zeros: " << biharmonicMatrix.nonZeros() << std::endl;
    
    Eigen::MatrixXd denseBiharmonic = biharmonicMatrix;
    
    // Number of modes to compute
    const int n_modes = 10; // Adjust as needed
    
    std::cout << "Computing eigenvalues and eigenvectors..." << std::endl;
    
    // Use Eigen's sparse eigenvalue solver
    Eigen::SparseMatrix<double> A = biharmonicMatrix;
    Eigen::MatrixXd eigenvectors;
    Eigen::VectorXd eigenvalues;

    try {
        // Create a .mat file to save the results
        std::string filename = "biharmonic_matrix.mat";
        std::cout << "Saving matrix to " << filename << std::endl;
        
        matioCpp::File matFile = matioCpp::File::Create(filename, matioCpp::FileVersion::MAT5);
        
        // Create a variable from the dense matrix
        auto matrixVar = matioCpp::make_variable("biharmonic_matrix", denseBiharmonic);
        
        // Save the parameters as well
        auto nxVar = matioCpp::make_variable("Nx", Nx);
        auto nyVar = matioCpp::make_variable("Ny", Ny);
        auto hVar = matioCpp::make_variable("h", h);
        auto lzVar = matioCpp::make_variable("Lz", Lz);
        auto eVar = matioCpp::make_variable("E", E);
        auto nuVar = matioCpp::make_variable("nu", nu);
        
        // Write variables to the file
        matFile.write(matrixVar);
        matFile.write(nxVar);
        matFile.write(nyVar);
        matFile.write(hVar);
        matFile.write(lzVar);
        matFile.write(eVar);
        matFile.write(nuVar);
        
        // Compute eigenvalues and eigenvectors
        // Using SelfAdjointEigenSolver for dense matrix
        // Note: For large sparse matrices, consider using Spectra or ArpackGeneralizedSelfAdjointEigenSolver
        // to match Python's eigs(sigma=0, which="LR")
        #ifdef EIGEN_USE_ARPACK
        // This would be preferable for large sparse matrices
        Eigen::ArpackGeneralizedSelfAdjointEigenSolver<Eigen::SparseMatrix<double>> arpackSolver;
        arpackSolver.compute(biharmonicMatrix, n_modes, "SM"); // "SM" = smallest magnitude (closest to 0)
        eigenvalues = arpackSolver.eigenvalues();
        eigenvectors = arpackSolver.eigenvectors();
        #else
        // Fallback to dense solver
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigSolver(denseBiharmonic);
        if (eigSolver.info() == Eigen::Success) {
            // Get n_modes largest magnitude eigenvalues (from the end)
            eigenvalues = eigSolver.eigenvalues().tail(n_modes);
            eigenvectors = eigSolver.eigenvectors().rightCols(n_modes);
        } else {
            std::cerr << "Eigenvalue computation failed!" << std::endl;
            return 1;
        }
        #endif
        
        // Sort eigenvalues and eigenvectors
        std::vector<std::pair<double, int>> eigenPairs;
        for (int i = 0; i < n_modes; i++) {
            eigenPairs.push_back(std::make_pair(eigenvalues(i), i));
        }
        std::sort(eigenPairs.begin(), eigenPairs.end());
        
        // Create sorted eigenvalues and eigenvectors
        Eigen::VectorXd sortedEigenvalues(n_modes);
        Eigen::MatrixXd sortedEigenvectors(eigenvectors.rows(), n_modes);
        
        for (int i = 0; i < n_modes; i++) {
            sortedEigenvalues(i) = eigenPairs[i].first;
            sortedEigenvectors.col(i) = eigenvectors.col(eigenPairs[i].second);
        }
        
        eigenvalues = sortedEigenvalues;
        eigenvectors = sortedEigenvectors;
        
        // Calculate norms and normalize eigenvectors
        const bool normalise_eigenvectors = true; // Set to false if normalization not required
        Eigen::VectorXd norms(n_modes);
        
        for (int i = 0; i < n_modes; i++) {
            // Implement double trapezoid integration for norm calculation
            Eigen::VectorXd eigenvec_squared = eigenvectors.col(i).array().square();
            
            // Use the double_trapezoid_flat implementation
            double norm = double_trapezoid_flat<double>(eigenvec_squared, h, h, Ny, Nx);
            
            norms(i) = norm;
            
            if (normalise_eigenvectors) {
                eigenvectors.col(i) /= std::sqrt(norm);
            }
        }
        
        std::cout << "First few eigenvalues:" << std::endl;
        for (int i = 0; i < std::min(5, n_modes); i++) {
            std::cout << eigenvalues(i) << std::endl;
        }
        
        // Save eigenvalues and eigenvectors to the .mat file
        auto eigenvaluesVar = matioCpp::make_variable("eigenvalues", eigenvalues);
        auto eigenvectorsVar = matioCpp::make_variable("eigenvectors", eigenvectors);
        auto normsVar = matioCpp::make_variable("norms", norms);
        
        matFile.write(eigenvaluesVar);
        matFile.write(eigenvectorsVar);
        matFile.write(normsVar);
        
        std::cout << "Matrix saved successfully." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error saving matrix: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
