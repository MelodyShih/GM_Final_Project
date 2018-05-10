#include "connectivity_opt.h"
#include <igl/edge_topology.h>
#include <igl/unique_edge_map.h>
#include <igl/flip_edge.h>
#include <igl/edge_flaps.h>
#include <igl/adjacency_list.h>

#include <math.h>
#include <iostream>

using namespace std;
using namespace Eigen;

void regularity_opt
    (
        const Eigen::MatrixXd& V,          // Vertices of the mesh (3D)
        Eigen::MatrixXi& F          // Faces
    )
{
    int uei;

    // for edge flip
    Eigen::MatrixXi E; 
    Eigen::MatrixXi uE; 
    Eigen::MatrixXi EMAP;
    std::vector<std::vector<int> > uE2E;

    igl::unique_edge_map(F,E,uE,EMAP,uE2E);

    // get EV,FE,EF relationship
    Eigen::MatrixXi EV; 
    Eigen::MatrixXi FE; 
    Eigen::MatrixXi EF;
    Eigen::MatrixXi EI;
    // igl::edge_topology(V,F,EV,FE,EF);

    // to calculate degree
    vector<vector<double> > A;

   	int v[4];
   	int face1index, face2index;
   	int old_degree[4], new_degree[4];
   	double cosnew, cosold;
   	Vector2d e1, e2;

    for (int ee = 0; ee < uE.rows(); ++ee)
    {
    	igl::edge_flaps(F,uE,EMAP,EF,EI);
    	//igl::edge_topology(V,F,EV,FE,EF);

    	int edge_index = ee;
    	// edge_index = 1166;
	    v[0] = uE(edge_index,0);
	    v[1] = uE(edge_index,1);
	    
	    if (EF(edge_index,0) == -1 || EF(edge_index,1) == -1)
	    {
	      	continue;
	    }

	    v[2] = F(EF(edge_index,0), EI(edge_index,0));
	    v[3] = F(EF(edge_index,1), EI(edge_index,1));

	    igl::adjacency_list(F,A);
	    for (int i = 0; i < 4; ++i)
	    {
	      old_degree[i] = A[v[i]].size();
	      if ( i<2 )
	        new_degree[i] = old_degree[i] - 1;
	      else
	        new_degree[i] = old_degree[i] + 1;
	    }

	    int old_var = 0, new_var = 0;
	    for (int i = 0; i < 4; ++i)
	    {
	      old_var += abs(old_degree[i] - 6);
	      new_var += abs(new_degree[i] - 6);
	    }

	    bool is_valid;
	    
	    e1 = V.row(v[2]) - V.row(v[3]);
	    e2 = V.row(v[1]) - V.row(v[3]); 
	    cosnew = e1.dot(e2)/(e1.norm()*e2.norm());

	    e1 = V.row(v[1]) - V.row(v[3]);
	    e2 = V.row(v[0]) - V.row(v[3]); 
	    cosold = e1.dot(e2)/(e1.norm()*e2.norm());
	    
	    is_valid = (cosold < cosnew);

	    if (new_var < old_var && is_valid)
	    {
	        igl::flip_edge(F,E,uE,EMAP,uE2E,edge_index);
	    }
    }
}

void regularity_statistic
	(
		const Eigen::MatrixXd& V,   // Vertices of the mesh (3D)
	    Eigen::MatrixXi& F          // Faces
	)
{
	vector<int> v(15,0);
	vector<vector<double> > A;

	igl::adjacency_list(F,A);
	for (int i = 0; i < V.rows(); ++i)
	{
		v[A[i].size() - 1] += 1;
	}
	for (int i = 0; i < v.size(); ++i)
	{
		cout << v[i] << endl;
	}
}

void face_aspect_ratio_opt
    (
        const Eigen::MatrixXd& V,          // Vertices of the mesh (3D)
        Eigen::MatrixXi& F,          // Faces
    	const Eigen::MatrixXd& P
    )
{
    int uei;

    // for edge flip
    Eigen::MatrixXi E; 
    Eigen::MatrixXi uE; 
    Eigen::MatrixXi EMAP;
    std::vector<std::vector<int> > uE2E;

    igl::unique_edge_map(F,E,uE,EMAP,uE2E);

    // get EV,FE,EF relationship
    Eigen::MatrixXi EF;
    Eigen::MatrixXi EI;

    // to calculate degree
    vector<vector<double> > A;

   	int v[4];
   	int face1index, face2index;
   	double old_area, old_peri, new_area, new_peri;
   	double old_face_aspect_ratio, new_face_aspect_ratio; 
   	Vector3d e1, e2;
   	double cosnew, cosold;
   	Vector2d ee1, ee2;

    for (int ee = 0; ee < 100; ++ee)
    {
    	igl::edge_flaps(F,uE,EMAP,EF,EI);

    	int edge_index = (rand() % uE.rows());
	    v[0] = uE(edge_index,0);
	    v[1] = uE(edge_index,1);
	    
	    if (EF(edge_index,0) == -1 || EF(edge_index,1) == -1)
	    {
	      continue;
	    }

	    v[2] = F(EF(edge_index,0), EI(edge_index,0));
	    v[3] = F(EF(edge_index,1), EI(edge_index,1));

	    e1 = V.row(v[1]) - V.row(v[0]);
	    e2 = V.row(v[2]) - V.row(v[0]);

	    old_area = 0.5*e1.cross(e2).norm();
	    old_peri = e1.norm() + e2.norm() + (V.row(v[2]) - V.row(v[1])).norm();

	    old_face_aspect_ratio = old_area/(old_peri*old_peri);
	    e2 = V.row(v[3]) - V.row(v[0]);

	    old_area = 0.5*e1.cross(e2).norm();
	    old_peri = e1.norm() + e2.norm() + (V.row(v[3]) - V.row(v[1])).norm();
	    
	    old_face_aspect_ratio += old_area/(old_peri*old_peri);

	    e1 = V.row(v[0]) - V.row(v[2]);
	    e2 = V.row(v[3]) - V.row(v[2]);

	    new_area = 0.5*e1.cross(e2).norm();
	    new_peri = e1.norm() + e2.norm() + (V.row(v[3]) - V.row(v[0])).norm();
	    new_face_aspect_ratio = new_area/(new_peri*new_peri);
	    
	    e1 = V.row(v[1]) - V.row(v[2]);

	    new_area = 0.5*e1.cross(e2).norm();
	    new_peri = e1.norm() + e2.norm() + (V.row(v[3]) - V.row(v[1])).norm();
	    new_face_aspect_ratio += new_area/(new_peri*new_peri);

	    bool is_valid;
	    
	    ee1 = P.row(v[2]) - P.row(v[3]);
	    ee2 = P.row(v[1]) - P.row(v[3]); 

	    cosnew = ee1.dot(ee2)/(ee1.norm()*ee2.norm());

	    ee1 = P.row(v[2]) - P.row(v[3]);
	    ee2 = P.row(v[0]) - P.row(v[3]); 
	    cosnew = min(ee1.dot(ee2)/(ee1.norm()*ee2.norm()), cosnew);
	    
	    ee1 = P.row(v[1]) - P.row(v[3]);
	    ee2 = P.row(v[0]) - P.row(v[3]); 
	    cosold = ee1.dot(ee2)/(ee1.norm()*ee2.norm());
	    
	    is_valid = (cosold < cosnew);

	    if ((old_area/old_peri) < (new_area/new_peri) && is_valid)
	    {
	        igl::flip_edge(F,E,uE,EMAP,uE2E,edge_index);
	    }
    }
}