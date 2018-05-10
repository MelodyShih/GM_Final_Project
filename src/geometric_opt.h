#ifndef GEOMETRIC_OPT_H
#define GEOMETRIC_OPT_H

#include <Eigen/Core>

Eigen::MatrixXd weighted_laplacian_flow
        (
                const Eigen::MatrixXd& V,          // Vertices of the mesh (3D)
                const Eigen::MatrixXi& F,          // Faces
                const Eigen::MatrixXd& P                 // Vertices of the mesh (2D)
        );
#endif
