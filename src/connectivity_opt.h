#ifndef CONNECTIVITY_OPT_H
#define CONNECTIVITY_OPT_H

#include <Eigen/Core>

void regularity_opt
    (
        const Eigen::MatrixXd& V,          // Vertices of the mesh (3D)
        Eigen::MatrixXi& F          // Faces
    );

void face_aspect_ratio_opt
    (
        const Eigen::MatrixXd& V,          // Vertices of the mesh (3D)
        Eigen::MatrixXi& F,
        const Eigen::MatrixXd& P          // Faces
    );
void regularity_statistic
	(
		const Eigen::MatrixXd& V,   // Vertices of the mesh (3D)
	    Eigen::MatrixXi& F          // Faces
	);
#endif
