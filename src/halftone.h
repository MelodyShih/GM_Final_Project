#ifndef HALFTONE_H
#define HALFTONE_H

#include <Eigen/Core>

void level0_diffusion
        (
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& V_uv,       // Adjacency triangle-triangle
        Eigen::VectorXd& control_map,
        const Eigen::MatrixXd& bnd_uv,
        std::vector<Eigen::RowVector2d>& p
        );

void level1_diffusion
        (
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& V_uv,       // Adjacency triangle-triangle
        Eigen::VectorXd& control_map,
        const Eigen::MatrixXd& bnd_uv,
        std::vector<Eigen::RowVector2d>& p
        );

void level2_diffusion
        (
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& V_uv,       // Adjacency triangle-triangle
        Eigen::VectorXd& control_map,
        const Eigen::MatrixXd& bnd_uv,
        std::vector<Eigen::RowVector2d>& p
        );

Eigen::MatrixXd error_diffusion
        (
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& V_uv,       // Adjacency triangle-triangle
        Eigen::VectorXd& control_map,
        const Eigen::MatrixXd& bnd_uv
        );

#endif
