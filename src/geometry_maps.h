#ifndef GEOMETRY_MAPS_H
#define GEOMETRY_MAPS_H

#include <Eigen/Core>

Eigen::VectorXd area_distortion_map
        (
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,          
        const Eigen::MatrixXd& V_uv        
        );

Eigen::VectorXd gaussian_curvature_map
        (
        const Eigen::MatrixXd& V,          
        const Eigen::MatrixXi& F,          
        const Eigen::MatrixXd& V_uv       
        );

Eigen::VectorXd mean_curvature_map
        (
        const Eigen::MatrixXd& V,          
        const Eigen::MatrixXi& F,          
        const Eigen::MatrixXd& V_uv       
        );
#endif
