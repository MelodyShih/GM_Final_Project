#include "geometric_opt.h"
#include <igl/adjacency_list.h>
#include <igl/doublearea.h>
#include <igl/boundary_loop.h>
#include <igl/edges.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/edge_topology.h>

#include <math.h>
#include <iostream>

using namespace std;
using namespace Eigen;

MatrixXd weighted_laplacian_flow
        (
	        const MatrixXd& V,          // Vertices of the mesh (3D)
	        const MatrixXi& F,          // Faces
	        const MatrixXd& P        	// Vertices of the mesh (2D)
        )
{
	MatrixXd dp = MatrixXd::Zero(P.rows(), P.cols());
	double w, totalarea;


	vector<vector<double> > A;


	VectorXd area_3D;
	area_3D.resize(F.rows());
	igl::doublearea(V,F,area_3D);
	totalarea = area_3D.sum();

	// get adjacency list
	igl::adjacency_list(F,A);

	VectorXi bnd;
	igl::boundary_loop(F,bnd);

	Eigen::MatrixXi EV; 
    Eigen::MatrixXi FE; 
    Eigen::MatrixXi EF;

    igl::edge_topology(V,F,EV,FE,EF);

	int v1, v2, v3;
	Vector3d vv1, vv2, vv3;
	Vector2d e1, e2;
	Vector3d ee1, ee2;
	int face1index, face2index;

	double cot1, cot2, cos1, cos2, area1, area2;
	
	for (int i = 0; i < EV.rows(); ++i)
	{
		face1index = EF(i,0);
		face2index = EF(i,1);
		
		if (face1index == -1 || face2index == -1)
		{
			continue;
		}
		
		v1 = EV(i,0);
		v2 = EV(i,1);

		if (F(face1index,0) != v1 && F(face1index,0) != v2){
			v3 = F(face1index,0);
		}
		
		if (F(face1index,1) != v1 && F(face1index,1) != v2){
			v3 = F(face1index,1);
		}
		
		if (F(face1index,2) != v1 && F(face1index,2) != v2){
			v3 = F(face1index,2);
		}
		assert(v3 < P.rows());
		area1 = area_3D(face1index);
		e1 = P.row(v1) - P.row(v3);
		e2 = P.row(v2) - P.row(v3);

		cos1 = e1.dot(e2) / (e1.norm()*e2.norm());
		cot1 = cos1 / sqrt(1-cos1*cos1);
		if (F(face2index,0) != v1 && F(face2index,0) != v2){
			v3 = F(face2index,0);
		}
		
		if (F(face2index,1) != v1 && F(face2index,1) != v2){
			v3 = F(face2index,1);
		}
		
		if (F(face2index,2) != v1 && F(face2index,2) != v2){
			v3 = F(face2index,2);
		}
		e1 = P.row(v1) - P.row(v3);
		e2 = P.row(v2) - P.row(v3);

		cos2 = e1.dot(e2) / (e1.norm()*e2.norm());
		cot2 = cos2 / sqrt(1-cos2*cos2);
		area2 = area_3D(face2index);
		w = area1*cot1 + area2*cot2;
		
		if (isnan(w) || isinf(w))
		{
			w = 0;
			cout << "i = " << i << endl;
			cout << F(face1index,0) << ","<< F(face1index,1) << ","<< F(face1index,2)<<endl;
			cout << F(face2index,0) << ","<< F(face2index,1) << ","<< F(face2index,2)<<endl;
			cout << "index = " << v1 << "," << v2 << ", "<<v3 << endl;
			cout << "a1 = " << area1 << " a2 = " << area2 << " cot1 = " << cot1 << " cot2 = " << cot2<< endl;

		}
		
		dp.row(v1) += w*(P.row(v2) - P.row(v1));
		dp.row(v2) += w*(P.row(v1) - P.row(v2));
	}

	return dp/totalarea;
}
#if 0
void smooth
	(
        const MatrixXd& V,          // Vertices of the mesh (3D)
        const MatrixXi& F,          // Faces
        const MatrixXd& P        	// Vertices of the mesh (2D)
    )
{
	Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V,F,L);

	VectorXd dblA;
    igl::doublearea(V,F,dblA);
    double area = 0.5*dblA.sum();
    MatrixXd BC;
    igl::barycenter(V,F,BC);
    RowVector3d centroid(0,0,0);
    for(int i = 0;i<BC.rows();i++)
    {
      centroid += 0.5*dblA(i)/area*BC.row(i);
    }
    V.rowwise() -= centroid;
}
#endif