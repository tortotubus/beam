#pragma once

#include <Eigen/Dense>
#include <iomanip>
#include <iostream>

namespace beam {

inline std::ostream& operator<<(std::ostream& os, const Eigen::VectorXd& v) {
    // Save stream state
    std::ios oldFmt(nullptr);
    oldFmt.copyfmt(os);

    // Set format
    os << std::scientific << std::setprecision(6);

    os << "[";
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        os << std::setw(15) << v(i);
        if (i < v.size() - 1) os << ",";
    }
    os << "]";

    // Restore stream state
    os.copyfmt(oldFmt);
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const Eigen::MatrixXd& m) {
    // Save stream state
    std::ios oldFmt(nullptr);
    oldFmt.copyfmt(os);

    // Set format
    os << std::scientific << std::setprecision(6);

    os << "[";
    for (Eigen::Index i = 0; i < m.rows(); ++i) {
        if (i == 0)
          os << "[";
        else 
          os << " [";
        for (Eigen::Index j = 0; j < m.cols(); ++j) {
            os << std::setw(15) << m(i,j);
            if (j < m.cols() - 1) os << ",";
        }
        os << "]";
        if (i < m.rows() - 1) os << ",";
        os << "\n";
    }
    os << "]";

    // Restore stream state
    os.copyfmt(oldFmt);
    return os;
}

} // namespace beam
