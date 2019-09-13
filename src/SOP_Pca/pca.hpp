//
// Created by symek on 9/13/19.
//

#pragma once

namespace pca {

    using StringVec = std::vector<std::string>;
    using Matrix    = Eigen::MatrixXd;
    using Vector    = Eigen::VectorXd;

// Should we just use EIGEN::QRMatrix?
    void orthogonalize_matrix(Matrix & matrix, int c=0);
// Reduced deformation space cutting out columns with eigenvalues bellow variance.
    bool computePCA(Matrix & matrix, Matrix & pcamatrix,
                    double variance, bool shift=false, bool orthogonalize=false);
// Saves matrix to a dummy binary format.
    bool write_matrix(const Matrix & matrix, const char * filename);
// Reads matrix from binary format.
    bool read_matrix(const char * filename, Matrix & matrix);

} // end of pca namespace