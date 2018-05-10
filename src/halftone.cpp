#include "halftone.h"
#include <math.h>
#include <igl/doublearea.h>
#include <igl/median.h>
#include <Eigen/Sparse>
#include <igl/png/readPNG.h>


using namespace std;
using namespace Eigen;

void level0_diffusion
        (
	        const MatrixXd& V,          	// Vertices of the mesh (3D)
	        const MatrixXi& F,          	// Faces
	        const MatrixXd& V_uv,        	// Vertices of the mesh (2D)
	        VectorXd& control_map,    		// Vertices of the mesh (2D)
	        const MatrixXd& bnd_uv,
	        std::vector<Eigen::RowVector2d>& p
        )
{
	MatrixXd sampling_aux;
	MatrixXd P;

	Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;

	// Save it to a PNG
	igl::png::readPNG("../data/area_halftone.png",R,G,B,A);
	sampling_aux = (0.2989 * R.cast<double>() + 0.5870 * G.cast<double>() + 0.1140 * B.cast<double>());
	sampling_aux = 256*MatrixXd::Ones(sampling_aux.rows(), sampling_aux.cols()) - sampling_aux;


	double u_min, u_max, v_min, v_max;
	double du, dv;
	double epsilon = 0.01;
	int nu, nv;

	nu = sampling_aux.cols();
	nv = sampling_aux.rows();

	u_min = V_uv.col(0).minCoeff();
	u_max = V_uv.col(0).maxCoeff();
	v_min = V_uv.col(1).minCoeff();
	v_max = V_uv.col(1).maxCoeff();

	du = (u_max - u_min)/nu;
	dv = (v_max - v_min)/nv;

	// std::vector<Eigen::RowVector2d> p;
	for (int i = 0; i < nv; ++i)
	{
		for (int j = 0; j < nu; ++j)
		{
			double ucoo = u_min+i*du;
			double vcoo = v_min+j*dv;
			double ucenter = 0.5*(u_min + u_max);
			double vcenter = 0.5*(v_min + v_max);
			double r = sqrt((ucoo - ucenter)*(ucoo - ucenter)+ (vcoo - vcenter)*(vcoo - vcenter));

			if ( sampling_aux(i,j) > 127.5 && r < (ucenter - u_min - 1e-6) && r < (vcenter - v_min - 1e-6))
			{
				p.push_back(Eigen::RowVector2d(u_min+i*du, v_min+j*dv));
			}
		}
	}

	// int nrows = p.size();
	// P.resize(nrows, 2);
	// for (int i = 0; i < nrows; ++i)
	// {
	// 	P.row(i) = p[i];
	// }
	// return p;
}

void level1_diffusion
        (
	        const MatrixXd& V,          	// Vertices of the mesh (3D)
	        const MatrixXi& F,          	// Faces
	        const MatrixXd& V_uv,        	// Vertices of the mesh (2D)
	        VectorXd& control_map,    		// Vertices of the mesh (2D)
	        const MatrixXd& bnd_uv,
	        std::vector<Eigen::RowVector2d>& p
        )
{
	MatrixXd sampling_aux;
	MatrixXd P;

	Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;

	// Save it to a PNG
	igl::png::readPNG("../data/curvature_halftone.png",R,G,B,A);
	sampling_aux = (0.2989 * R.cast<double>() + 0.5870 * G.cast<double>() + 0.1140 * B.cast<double>());
	sampling_aux = 256*MatrixXd::Ones(sampling_aux.rows(), sampling_aux.cols()) - sampling_aux;


	double u_min, u_max, v_min, v_max;
	double du, dv;
	double epsilon = 0.01;
	int nu, nv;

	nu = sampling_aux.cols();
	nv = sampling_aux.rows();

	u_min = V_uv.col(0).minCoeff();
	u_max = V_uv.col(0).maxCoeff();
	v_min = V_uv.col(1).minCoeff();
	v_max = V_uv.col(1).maxCoeff();

	du = (u_max - u_min)/nu;
	dv = (v_max - v_min)/nv;

	// std::vector<Eigen::RowVector2d> p;
	for (int i = 0; i < nv; ++i)
	{
		for (int j = 0; j < nu; ++j)
		{
			double ucoo = u_min+i*du;
			double vcoo = v_min+j*dv;
			double ucenter = 0.5*(u_min + u_max);
			double vcenter = 0.5*(v_min + v_max);
			double r = sqrt((ucoo - ucenter)*(ucoo - ucenter)+ (vcoo - vcenter)*(vcoo - vcenter));

			if ( sampling_aux(i,j) > 127.5 && r < (ucenter - u_min - 1e-6) && r < (vcenter - v_min - 1e-6))
			{
				p.push_back(Eigen::RowVector2d(u_min+i*du, v_min+j*dv));
			}
		}
	}
}

void level2_diffusion
        (
	        const MatrixXd& V,          	// Vertices of the mesh (3D)
	        const MatrixXi& F,          	// Faces
	        const MatrixXd& V_uv,        	// Vertices of the mesh (2D)
	        VectorXd& control_map,    		// Vertices of the mesh (2D)
	        const MatrixXd& bnd_uv,
	        std::vector<Eigen::RowVector2d>& p
        )
{
	cout<<p.size()<<endl;
	MatrixXd sampling_aux;
	MatrixXd P;

	Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;

	// Save it to a PNG
	igl::png::readPNG("../data/curvature_fine_halftone.png",R,G,B,A);
	sampling_aux = (0.2989 * R.cast<double>() + 0.5870 * G.cast<double>() + 0.1140 * B.cast<double>());
	sampling_aux = 256*MatrixXd::Ones(sampling_aux.rows(), sampling_aux.cols()) - sampling_aux;


	double u_min, u_max, v_min, v_max;
	double du, dv;
	double epsilon = 0.01;
	int nu, nv;

	nu = sampling_aux.cols();
	nv = sampling_aux.rows();

	u_min = V_uv.col(0).minCoeff();
	u_max = V_uv.col(0).maxCoeff();
	v_min = V_uv.col(1).minCoeff();
	v_max = V_uv.col(1).maxCoeff();

	du = (u_max - u_min)/nu;
	dv = (v_max - v_min)/nv;

	// std::vector<Eigen::RowVector2d> p;
	for (int i = 0; i < nv; ++i)
	{
		for (int j = 0; j < nu; ++j)
		{
			double ucoo = u_min+i*du;
			double vcoo = v_min+j*dv;
			double ucenter = 0.5*(u_min + u_max);
			double vcenter = 0.5*(v_min + v_max);
			double r = sqrt((ucoo - ucenter)*(ucoo - ucenter)+ (vcoo - vcenter)*(vcoo - vcenter));

			if ( sampling_aux(i,j) > 127.5 && r < (ucenter - u_min - 1e-6) && r < (vcenter - v_min - 1e-6))
			{
				p.push_back(Eigen::RowVector2d(u_min+i*du, v_min+j*dv));
			}
		}
	}
}



MatrixXd error_diffusion
        (
	        const MatrixXd& V,          	// Vertices of the mesh (3D)
	        const MatrixXi& F,          	// Faces
	        const MatrixXd& V_uv,        	// Vertices of the mesh (2D)
	        VectorXd& control_map,    // Vertices of the mesh (2D)
	        const MatrixXd& bnd_uv
        )
{
	MatrixXd sampling_aux;
	MatrixXd P;
	// double epsilon = 0.01;

	Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;

	// Save it to a PNG
	igl::png::readPNG("../data/curvature_halftone.png",R,G,B,A);
	sampling_aux = (0.2989 * R.cast<double>() + 0.5870 * G.cast<double>() + 0.1140 * B.cast<double>());
	sampling_aux = 256*MatrixXd::Ones(sampling_aux.rows(), sampling_aux.cols()) - sampling_aux;


	double u_min, u_max, v_min, v_max;
	double du, dv;
	double epsilon = 0.01;
	int nu, nv;

	nu = sampling_aux.cols();
	nv = sampling_aux.rows();

	u_min = V_uv.col(0).minCoeff();
	u_max = V_uv.col(0).maxCoeff();
	v_min = V_uv.col(1).minCoeff();
	v_max = V_uv.col(1).maxCoeff();

	du = (u_max - u_min)/nu;
	dv = (v_max - v_min)/nv;

	std::vector<Eigen::RowVector2d> p;
	for (int i = 0; i < nv; ++i)
	{
		for (int j = 0; j < nu; ++j)
		{
			double ucoo = u_min+i*du;
			double vcoo = v_min+j*dv;
			double ucenter = 0.5*(u_min + u_max);
			double vcenter = 0.5*(v_min + v_max);
			double r = sqrt((ucoo - ucenter)*(ucoo - ucenter)+ (vcoo - vcenter)*(vcoo - vcenter));

			if ( sampling_aux(i,j) > 127.5 && r < (ucenter - u_min - 1e-6) && r < (vcenter - v_min - 1e-6))
			{
				p.push_back(Eigen::RowVector2d(u_min+i*du, v_min+j*dv));
			}
		}
	}

	int nowrow = p.size();

	// Add boundary points
	for (int i = 0; i < bnd_uv.rows(); ++i)
	{
		p.push_back(Eigen::RowVector2d(bnd_uv(i,0),bnd_uv(i,1) ));
	}

	int nrows = p.size();
	P.resize(nrows, 2);
	for (int i = 0; i < nrows; ++i)
	{
		P.row(i) = p[i];
	}
	return P;
}

#if 0
	MatrixXd sampling_aux;
	MatrixXd P;


	double u_min, u_max, v_min, v_max;
	double du, dv;
	double epsilon = 0.01;
	int nu, nv;

	nu = sampling_aux.cols();
	nv = sampling_aux.rows();

	nu = n;
	nv = n;

	sampling_aux.resize(nu, nv);
	sampling_aux.setZero();

	u_min = V_uv.col(0).minCoeff();
	u_max = V_uv.col(0).maxCoeff();
	v_min = V_uv.col(1).minCoeff();
	v_max = V_uv.col(1).maxCoeff();

	du = (u_max - u_min)/nu;
	dv = (v_max - v_min)/nv;

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
		
		if (u_index == nu || v_index == nv)
		{
			continue;
		}
		sampling_aux(u_index,v_index) += control_map[i];
	}

    du = (u_max - u_min)/nu;
	dv = (v_max - v_min)/nv;

    double dt = 0.01*du*dv;
    std::vector<Eigen::RowVector2d> p;
    
    for (int timestep = 0; timestep < step; ++timestep)
    {
		for (int y = 1; y < nv-1; y++){

	        for (int x = 1; x < nu-1; x++){

	            double ucoo = u_min+x*du;
				double vcoo = v_min+y*dv;
				double ucenter = 0.5*(u_min + u_max);
				double vcenter = 0.5*(v_min + v_max);
				double r = (ucoo - ucenter)*(ucoo - ucenter)+ (vcoo - vcenter)*(vcoo - vcenter);
				if ( sampling_aux(x,y) > threshold && r < (ucenter - u_min -epsilon) && r < (vcenter - v_min - epsilon))
				{
					sampling_aux(x,y) += dt/(du*dv) * (sampling_aux(x+1,y)+
		            						   		   sampling_aux(x-1,y)+
		            						           sampling_aux(x,y+1)+
		            						           sampling_aux(x,y-1)-
		            						           4*sampling_aux(x,y));
				}
	        }
	    }
    }

	for (int y = 1; y < nv-1; y++){

        for (int x = 1; x < nu-1; x++){

            double ucoo = u_min+x*du;
			double vcoo = v_min+y*dv;
			double ucenter = 0.5*(u_min + u_max);
			double vcenter = 0.5*(v_min + v_max);
			double r = (ucoo - ucenter)*(ucoo - ucenter)+ (vcoo - vcenter)*(vcoo - vcenter);

			if ( sampling_aux(x,y) > threshold && r < (ucenter - u_min -epsilon) && r < (vcenter - v_min - epsilon))
			{
				p.push_back(Eigen::RowVector2d(u_min+x*du, v_min+y*dv));
			}

        }
    }


	int nrows = p.size();
	P.resize(nrows, 2);
	cout << V.rows() << endl;
	cout << nrows << endl;
	for (int i = 0; i < nrows; ++i)
	{
		P.row(i) = p[i];
	}
#endif
