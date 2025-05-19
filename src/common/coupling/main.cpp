#include "coupling.h"
#include <iostream>
#include <Eigen/src/Core/Matrix.h>
#include <matioCpp/matioCpp.h>
#include "vk_utils/vk_utils.h"

int main(int argc, char *argv[])
{
    int n_psi = 10;
    int n_phi = 10;
    double lx = 0.4;
    double ly = 0.7;

    // calculate the eigenvalues for the plate
    Eigen::Matrix lambda_mu_2d = calculate_plate_eigenvalues<double>(
        n_phi,
        n_phi,
        lx,
        ly
    );

    Eigen::VectorX<double> lambda_mu;
    Eigen::MatrixX<int> selected_indices;
    select_modes_and_eigenvalues(
        lambda_mu_2d,
        n_phi,
        lambda_mu,
        selected_indices
    );

    // Extract indices into separate x and y vectors
    Eigen::VectorXi selected_indices_x(n_phi);
    Eigen::VectorXi selected_indices_y(n_phi);
    
    for (int i = 0; i < n_phi; ++i) {
        selected_indices_x(i) = static_cast<int>(selected_indices(i, 0));
        selected_indices_y(i) = static_cast<int>(selected_indices(i, 1));
    }

#if 1
    Matrix H1;
    compute_coupling_matrix(
        n_psi,
        n_phi,
        lx,
        ly,
        selected_indices_x,
        selected_indices_y,
        H1
    );

    // reshape the H1 tensor to (n_psi * n_phi, n_phi)
    Matrix H1_transposed = Eigen::Map<const Matrix>(H1.transpose().data(), n_psi * n_phi, n_phi);

    // save the H1 tensor to a .mat file
    matioCpp::File matFile_H1 = matioCpp::File::Create(
        "cpp_H_tensor_rectangular.mat",
        matioCpp::FileVersion::MAT5
    );
    auto H1_var = matioCpp::make_variable("H1_cpp", H1_transposed);
    matFile_H1.write(H1_var);
#else

    // load the eigenvalues and eigenvectors from a .mat file
    matioCpp::File eigenproblem_results("/home/diaz/projects/VKGong/matlab/VK-Gong/eigenproblem_results.mat");
    if (!eigenproblem_results.isOpen()) {
        std::cout << "set_couplings_and_eigenvalues: could not open file '" << "eigenproblem_results.mat" << "'" << std::endl;
        return 1;
    }
    matioCpp::MultiDimensionalArray<double> vecs_mat = eigenproblem_results.read("VEC_old").asMultiDimensionalArray<double>();
    matioCpp::Vector<double> vals_mat = eigenproblem_results.read("VAL_old").asVector<double>();
    matioCpp::MultiDimensionalArray<double> coeff0 = eigenproblem_results.read("coeff0").asMultiDimensionalArray<double>();
    matioCpp::MultiDimensionalArray<double> coeff1 = eigenproblem_results.read("coeff1").asMultiDimensionalArray<double>();
    matioCpp::MultiDimensionalArray<double> coeff2 = eigenproblem_results.read("coeff2").asMultiDimensionalArray<double>();
    Matrix coeffs0_mat = Eigen::Map<Matrix>(coeff0.data(), coeff0.dimensions()[0], coeff0.dimensions()[1]);
    Matrix coeffs1_mat = Eigen::Map<Matrix>(coeff1.data(), coeff1.dimensions()[0], coeff1.dimensions()[1]);
    Matrix coeffs2_mat = Eigen::Map<Matrix>(coeff2.data(), coeff2.dimensions()[0], coeff2.dimensions()[1]);

    Eigen::VectorX<double> vals_matlab = Eigen::Map<Eigen::VectorX<double>>(vals_mat.data(), vals_mat.size());
    Eigen::MatrixX<double> vecs_matlab = Eigen::Map<Eigen::MatrixX<double>>(vecs_mat.data(), vecs_mat.dimensions()[0], vecs_mat.dimensions()[1]);


    // 1. Assemble K and M
    auto [K, M] = assemble_K_and_M(n_psi, lx, ly);

    // 2. Solve Generalized Eigenvalue Problem K*v = lambda*M*v
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> ges;
    ges.compute(K, M);

    Vector vals = ges.eigenvalues();
    Matrix vecs = ges.eigenvectors();



    AiryCoefficients coeffs = airy_stress_coefficients(n_psi, lx, ly, vals, vecs);
    // AiryCoefficients coeffs_matlab;
    // coeffs_matlab.coeff0 = coeffs0_mat;
    // coeffs_matlab.coeff1 = coeffs1_mat;
    // coeffs_matlab.coeff2 = coeffs2_mat;
    // coeffs_matlab.S = (n_psi * n_psi) / 2;

    std::cout << "coeffs.S: " << coeffs.S << std::endl;
    std::cout << "n_phi: " << n_phi << std::endl;
    std::cout << "n_psi: " << n_psi << std::endl;
    // 4. Compute H Tensors
    Matrix H0_mat = Matrix::Zero(coeffs.S, n_phi*n_phi);
    Matrix H1_mat = Matrix::Zero(coeffs.S, n_phi*n_phi);
    Matrix H2_mat = Matrix::Zero(coeffs.S, n_phi*n_phi);
    H_tensor_rectangular(
        coeffs,
        n_phi,
        n_psi, // Pass original n_psi
        lx,
        ly,
        selected_indices_x,
        selected_indices_y,
        H0_mat,
        H1_mat,
        H2_mat
    );

    // we need to slice the H tensors to keep the shape (n_psi, n_phi * n_phi)
    H0_mat = H0_mat.block(0, 0, n_psi, n_phi*n_phi);
    H1_mat = H1_mat.block(0, 0, n_psi, n_phi*n_phi);
    H2_mat = H2_mat.block(0, 0, n_psi, n_phi*n_phi);

    // reshape them
    Matrix H0 = Eigen::Map<const Matrix>(H0_mat.transpose().data(), n_psi * n_phi, n_phi);
    Matrix H1 = Eigen::Map<const Matrix>(H1_mat.transpose().data(), n_psi * n_phi, n_phi);
    Matrix H2 = Eigen::Map<const Matrix>(H2_mat.transpose().data(), n_psi * n_phi, n_phi);

    matioCpp::File matFile_H = matioCpp::File::Create(
        "cpp_H_tensor_rectangular.mat",
        matioCpp::FileVersion::MAT5
    );
    auto H0_var = matioCpp::make_variable("H0_cpp", H0);
    auto H1_var = matioCpp::make_variable("H1_cpp", H1);
    auto H2_var = matioCpp::make_variable("H2_cpp", H2);
    matFile_H.write(H0_var);
    matFile_H.write(H1_var);
    matFile_H.write(H2_var);

    // H tensors have the shape (n_psi, n_phi, n_phi)

    // save the Airy coefficients to a .mat file
    matioCpp::File matFile = matioCpp::File::Create(
        "cpp_results.mat",
        matioCpp::FileVersion::MAT5
    );

    std::cout << "writing to file" << std::endl;

    auto auto_vec_var = matioCpp::make_variable("auto_vec_cpp", coeffs.auto_vec);
    auto vals_var = matioCpp::make_variable("vals_cpp", vals);
    auto vecs_var = matioCpp::make_variable("vecs_cpp", coeffs.coeff_sorted);
    auto coeff0Var = matioCpp::make_variable("coeff0_cpp", coeffs.coeff0);
    auto coeff1Var = matioCpp::make_variable("coeff1_cpp", coeffs.coeff1);
    auto coeff2Var = matioCpp::make_variable("coeff2_cpp", coeffs.coeff2);
    auto lambda_mu_var = matioCpp::make_variable("lambda_mu_cpp", lambda_mu);
    auto selected_indices_x_var = matioCpp::make_variable("selected_indices_x_cpp", selected_indices_x);
    auto selected_indices_y_var = matioCpp::make_variable("selected_indices_y_cpp", selected_indices_y);
    // auto H0_0_var = matioCpp::make_variable("H0_0_cpp", H.H0[0]);
    
    matFile.write(auto_vec_var);
    matFile.write(vals_var);
    matFile.write(vecs_var);
    matFile.write(coeff0Var);
    matFile.write(coeff1Var);
    matFile.write(coeff2Var);
    matFile.write(lambda_mu_var);
    matFile.write(selected_indices_x_var);
    matFile.write(selected_indices_y_var);
    // matFile.write(H0_0_var);
#endif
    return 0;
}