#ifndef GEOMETRY_MAPS
#define GEOMETRY_MAPS

#include <Eigen/Core>

Eigen::VectorXd area_distortion_map
        (
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& V_uv        // Adjacency triangle-triangle
        );

Eigen::VectorXd gaussian_curvature_map
        (
        const Eigen::MatrixXd& V,          // Vertices of the mesh (3D)
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& V_uv       // Vertices of the mesh (2D)
        );

Eigen::VectorXd mean_curvature_map
        (
        const Eigen::MatrixXd& V,          // Vertices of the mesh (3D)
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& V_uv       // Vertices of the mesh (2D)
        );
#endif
