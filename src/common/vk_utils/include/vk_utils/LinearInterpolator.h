#pragma once

#include <vector>
#include <Eigen/Dense>

template <typename T>
class LinearInterpolator {
public:
    LinearInterpolator() 
        : m_val(0)
        , m_increment(0)
        , m_counter(0)
        , m_delta(441)
        , m_target(0)
    {}
    
    LinearInterpolator(unsigned int delta) : m_val(0), m_increment(0), m_counter(0), m_delta(delta), m_target(0) {}
    
    void setup(unsigned int delta) {
        m_val = 0;
        m_increment = 0;
        m_counter = 0;
        m_delta = delta;
        m_target = 0;
    }
    
    void setTarget(T target) {
        m_target = target;
        calcIncrement(m_target, m_delta);
    }
    
    void setValue(T value) {
        m_val = value;
        m_counter = 0;
        m_increment = 0;
    }
    
    void setDelta(unsigned int delta) {
        m_delta = delta;
        // calcIncrement(m_target, m_delta);
    }
    
    T process(bool bound=true) {
        if (m_counter > 0) {
            m_val += m_increment;
            m_counter--;
        }
        if (bound)
            m_val = boundValue();
        return m_val;
    }
    
    T getTarget() const { return m_target; }
    T getValue() const { return m_val; }
    unsigned int getDelta() const { return m_delta; }
    
    bool isFinished() const { return m_counter == 0; }
    
private:
    T m_val;
    T m_increment;
    int m_counter;
    unsigned int m_delta;
    T m_target;
    
    void calcIncrement(T target, unsigned int delta) {
        if(delta == 0)
            m_increment = 0.0;
        else
            m_increment = (target - m_val) / static_cast<T>(delta);
        
        if(m_increment == 0.0) {
            m_counter = 0;
            m_val = target;
        }
        else {
            m_counter = delta;
        }
    }
    
    T boundValue() {
        if (m_increment > 0 && m_val > m_target)
            m_val = m_target;
        else if (m_increment < 0 && m_val < m_target)
            m_val = m_target;
        return m_val;
    }
};

template <typename T>
class VectorInterpolator {
public:
    VectorInterpolator() : m_initialized(false) {}
    
    void resize(int size, unsigned int delta) {
        m_lerps.resize(size);
        for (auto& lerp : m_lerps) {
            lerp.setup(delta);
        }
        m_initialized = true;
    }
    
    void setTarget(const Eigen::VectorX<T>& target) {
        if (!m_initialized) {
            return;
        }
        if (m_lerps.size() != target.size()) {
            m_lerps.resize(target.size());
        }
        
        for (int i = 0; i < target.size(); i++) {
            m_lerps[i].setTarget(target[i]);
        }
    }
    
    void setValue(const Eigen::VectorX<T>& value) {
        if (!m_initialized) {
            return;
        }
        if (m_lerps.size() != value.size()) {
            m_lerps.resize(value.size());
        }
        
        for (int i = 0; i < value.size(); i++) {
            m_lerps[i].setValue(value[i]);
        }
    }
    
    void setDelta(unsigned int delta) {
        for (auto& lerp : m_lerps) {
            lerp.setDelta(delta);
        }
    }
    
    bool process(Eigen::VectorX<T>& result) {
        if (!m_initialized) {
            return false;
        }
        for (int i = 0; i < m_lerps.size(); i++) {
            result[i] = m_lerps[i].process();
        }
        return true;
    }
    
    bool isFinished() const {
        for (const auto& lerp : m_lerps) {
            if (!lerp.isFinished()) {
                return false;
            }
        }
        return true;
    }
    
private:
    bool m_initialized;
    std::vector<LinearInterpolator<T>> m_lerps;
};

template <typename T>
class MatrixInterpolator {
public:
    MatrixInterpolator() : m_initialized(false) {}

    void resize(int rows, int cols, unsigned int delta) {
        m_lerps.resize(rows);
        for (auto& row : m_lerps) {
            row.resize(cols);
            for (auto& lerp : row) {
                lerp.setup(delta);
            }
        }
        m_initialized = true;
    }

    void setTarget(const Eigen::MatrixX<T>& target) {
        if (!m_initialized) {
            return;
        }
        // if (m_lerps.size() != target.rows() || m_lerps[0].size() != target.cols()) {
            // return;  // Ignore if dimensions don't match
        // }
        for (int i = 0; i < target.rows(); ++i) {
            for (int j = 0; j < target.cols(); ++j) {
                m_lerps[i][j].setTarget(target(i, j));
            }
        }
    }

    void setValue(const Eigen::MatrixX<T>& value) {
        if (!m_initialized) {
            return;
        }
        // if (m_lerps.size() != value.rows() || m_lerps[0].size() != value.cols()) {
            // return;  // Ignore if dimensions don't match
        // }
        for (int i = 0; i < value.rows(); ++i) {
            for (int j = 0; j < value.cols(); ++j) {
                m_lerps[i][j].setValue(value(i, j));
            }
        }
    }

    void setDelta(unsigned int delta) {
        for (auto& row : m_lerps) {
            for (auto& lerp : row) {
                lerp.setDelta(delta);
            }
        }
    }

    bool process(Eigen::MatrixX<T>& result) {
        int rows = m_lerps.size();
        int cols = m_lerps[0].size();

        // check the result matrix has the correct dimensions
        if (result.rows() != rows || result.cols() != cols) {
            return false;
        }

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = m_lerps[i][j].process();
            }
        }
        return true;
    }

    bool isFinished() const {
        for (const auto& row : m_lerps) {
            for (const auto& lerp : row) {
                if (!lerp.isFinished()) {
                    return false;
                }
            }
        }
        return true;
    }


private:
    bool m_initialized;
    std::vector<std::vector<LinearInterpolator<T>>> m_lerps;
};