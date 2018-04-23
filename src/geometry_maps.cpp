#include "geometry_maps.h"
#include <igl/doublearea.h>

#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

#include <igl/principal_curvature.h>

using namespace std;
using namespace Eigen;

VectorXd area_distortion_map
        (
        const MatrixXd& V,          // Vertices of the mesh (3D)
        const MatrixXi& F,          // Faces
        const MatrixXd& V_uv        // Vertices of the mesh (2D)
        )
{
	VectorXd area_distortion;
	VectorXd area_3D, area_2D;

	igl::doublearea(V_uv,F,area_2D);
	igl::doublearea(V,F,area_3D);

	area_distortion = area_3D.array()/area_2D.array();
	return area_distortion;
}

VectorXd gaussian_curvature_map
        (
        const MatrixXd& V,          // Vertices of the mesh (3D)
        const MatrixXi& F,          // Faces
        const MatrixXd& V_uv       // Vertices of the mesh (2D)
        )
{
	VectorXd Kg;
	// Compute integral of Gaussian curvature
	igl::gaussian_curvature(V,F,Kg);
	// Compute mass matrix
	SparseMatrix<double> M,Minv;
	igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
	igl::invert_diag(M,Minv);
	// Divide by area to get integral average
	Kg = (Minv*Kg).eval();

	return Kg;
}

VectorXd mean_curvature_map
        (
        const MatrixXd& V,          // Vertices of the mesh (3D)
        const MatrixXi& F,          // Faces
        const MatrixXd& V_uv       // Vertices of the mesh (2D)
        )
{
	VectorXd H;
	
	// Compute curvature directions via quadric fitting
	MatrixXd PD1,PD2;
	VectorXd PV1,PV2;
	igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);
	
	// mean curvature
	H = 0.5*(PV1+PV2);

	return H;
}