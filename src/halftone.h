#ifndef HALFTONE
#define HALFTONE

#include <Eigen/Core>

Eigen::MatrixXd error_diffusion
        (
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& V_uv,       // Adjacency triangle-triangle
        const Eigen::VectorXd& control_map
        );

#endif
