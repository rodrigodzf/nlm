#pragma once

#include <Eigen/Dense>
#include <vk_utils/vk_utils.h>


template <typename DerivedV>
inline void clip_inplace(
          Eigen::MatrixBase<DerivedV>& v,
          const double min,
          const double max)
{
  v = v.cwiseMin(max).cwiseMax(min);
}

template <typename T>
class ParallelFilter
{
    public:
        ParallelFilter() : m_initialized(false), m_n_modes(0) {}

        void resize(int n_modes) {
            m_n_modes = n_modes;
            m_q = Eigen::VectorX<T>::Zero(n_modes);
            m_q_prev = Eigen::VectorX<T>::Zero(n_modes);
            m_q_next = Eigen::VectorX<T>::Zero(n_modes);
            m_A_inv = Eigen::VectorX<T>::Zero(n_modes);
            m_B = Eigen::VectorX<T>::Zero(n_modes);
            m_C = Eigen::VectorX<T>::Zero(n_modes);
            m_nl = Eigen::VectorX<T>::Zero(n_modes);
            if (!m_initialized) {
                m_initialized = true;
            }
        }

        int get_n_modes() const {
            return m_n_modes;
        }

        void set_coefficients(const Eigen::VectorX<T>& m_gamma2_mu, const Eigen::VectorX<T>& m_omega_mu_squared, double m_dt) {
            if (!m_initialized) {
                return;
            }
            calculate_coefficients_tf(
                m_gamma2_mu,
                m_omega_mu_squared,
                m_B,
                m_C,
                m_A_inv,
                m_dt
            );
        }

        void update_nonlinearity(const Eigen::MatrixX<T>& H_scaled, const int n_psi, const int n_phi) {
            if (!m_initialized) {
                return;
            }
            calculate_nonlinear_vk(H_scaled, m_q, m_nl, n_psi, n_phi);
        }

        void update_nonlinearity(const Eigen::VectorX<T>& lambda_mu, const Eigen::VectorX<T>& plate_tau_with_norms) {
            if (!m_initialized) {
                return;
            }
            calculate_nonlinear_berger(lambda_mu, plate_tau_with_norms, m_q, m_nl);
        }

        const Eigen::VectorX<T>& get_q() const {
            return m_q;
        }

        void operator()(const Eigen::VectorX<T>& input) {
            if (!m_initialized) {
                return;
            }
            
            // Modal update
            // clip_inplace(m_q, -100, 100);

            m_q_next.noalias() = m_B.cwiseProduct(m_q) +
                                 m_C.cwiseProduct(m_q_prev) +
                                 m_A_inv.cwiseProduct(input - m_nl);

            // Clip the output to avoid overflow - using same scale as input

            m_q_prev = m_q;
            m_q = m_q_next;
        }

        void reset() {
            if (!m_initialized) {
                return;
            }
            m_q.setZero();
            m_q_prev.setZero();
            m_q_next.setZero();
            m_nl.setZero();
        }

    private:
        int m_n_modes;
        Eigen::VectorX<T> m_q;
        Eigen::VectorX<T> m_q_prev;
        Eigen::VectorX<T> m_q_next;

        Eigen::VectorX<T> m_A_inv;
        Eigen::VectorX<T> m_B;
        Eigen::VectorX<T> m_C;
        Eigen::VectorX<T> m_nl;
        bool m_initialized;
};