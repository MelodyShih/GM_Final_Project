#include "halftone.h"
#include <igl/doublearea.h>
#include <Eigen/Sparse>


using namespace std;
using namespace Eigen;

MatrixXd error_diffusion
        (
        const MatrixXd& V,          // Vertices of the mesh (3D)
        const MatrixXi& F,          // Faces
        const MatrixXd& V_uv,        // Vertices of the mesh (2D)
        const VectorXd& control_map     // Vertices of the mesh (2D)
        )
{
	MatrixXd sampling_aux;
	MatrixXd P;

	double u_min, u_max, v_min, v_max;
	double du, dv;
	double epsilon = 0.01;
	int nu, nv;

	nu = 200;
	nv = 200;

	sampling_aux.resize(nu, nv);
	sampling_aux.setZero();
	cout << sampling_aux.rows() << " " << sampling_aux.cols() << endl;

	u_min = V_uv.col(0).minCoeff();
	u_max = V_uv.col(0).maxCoeff();
	v_min = V_uv.col(1).minCoeff();
	v_max = V_uv.col(1).maxCoeff();

	du = (u_max - u_min + epsilon)/nu;
	dv = (v_max - v_min + epsilon)/nv;

	cout << "u:" << u_min << "," << u_max << endl;
	cout << "v:" << v_min << "," << v_max << endl;
	cout << F.rows() << ", "<< control_map.size() << ", " << V_uv.rows() << endl;

	int u_index, v_index;
	double u,v;

	for (int i = 0; i < F.rows(); ++i)
	{
		u = 0;
		v = 0;
		for (int j = 0; j < F.cols(); ++j)
		{
			u += V_uv(F(i,j),0);
			v += V_uv(F(i,j),1);
		}
		u = u/F.cols();
		v = v/F.cols();

		u_index = floor((u - u_min)/du);
		v_index = floor((v - v_min)/dv);
		assert(u_index < 200 && v_index < 200);
		sampling_aux(u_index,v_index) += control_map[i];
	}

	double threshold = 1;
	std::vector<Eigen::RowVector3d> p;
	for (int i = 0; i < nu; ++i)
	{
		for (int j = 0; j < nv; ++j)
		{
			if (sampling_aux(i,j) >= threshold || i==0 || j==0 || i==nu-1 || j==nv-1)
			{
				p.push_back(Eigen::RowVector3d(u_min+(i+0.5)*du, v_min+(j+0.5)*dv, 0));
			}
		}
	}
	int nrows = p.size();
	cout << nrows << endl;
	P.resize(nrows, 3);
	for (int i = 0; i < nrows; ++i)
	{
		P.row(i) = p[i];
	}
	// cout << "sampling_aux:" << sampling_aux.minCoeff() << "," << sampling_aux.maxCoeff() << endl;
	return P;
}
