#include "vertex_to_triangle.h"
#include <igl/doublearea.h>
#include <igl/matlab/matlabinterface.h>
#include <igl/triangle/triangulate.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>

using namespace std;

static
double area(double x1, double y1, 
		   double x2, double y2, 
		   double x3, double y3)
{
   return abs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.0);
}
 

static
bool isInside(double x1, double y1, 
			  double x2, double y2, 
			  double x3, double y3, 
			  double x, double y)
{   
   /* Calculate area of triangle ABC */
   double A = area (x1, y1, x2, y2, x3, y3);
 
   /* Calculate area of triangle PBC */  
   double A1 = area (x, y, x2, y2, x3, y3);
 
   /* Calculate area of triangle PAC */  
   double A2 = area (x1, y1, x, y, x3, y3);
 
   /* Calculate area of triangle PAB */   
   double A3 = area (x1, y1, x2, y2, x, y);
   
   /* Check if sum of A1, A2 and A3 is same as A */
   if (abs(A - A1 - A2 - A3) < 1e-15)
   		return true;
   else
   		return false;

   // return (A == A1 + A2 + A3);
}

// Matlab instance
Engine* engine;

void vertex_to_triangle
        (
        	const Eigen::MatrixXd& P,
        	const Eigen::MatrixXi& F,
	        const Eigen::MatrixXd& V_uv,
	        Eigen::MatrixXi& F_new, 
	        Eigen::MatrixXd& V_new
        )
{
	// Eigen::MatrixXi F_aux;
	// Eigen::MatrixXi F_new;
	Eigen::MatrixXi E, H;
	Eigen::VectorXi bnd;
	Eigen::MatrixXd VV;

	
	igl::boundary_loop(F,bnd);
	E.resize(bnd.size(),2);
	for (int i = 0; i < bnd.size(); ++i)
	{
		E(i,0) = P.rows() - bnd.size() + i;
		E(i,1) = P.rows() - bnd.size() + i + 1;
	}
	E(bnd.size()-1,1) = P.rows() - bnd.size();

	igl::triangle::triangulate(P,E,H,"a0.005q",V_new,F_new);
	VV = Eigen::MatrixXd::Zero(V_new.rows(),3);
	VV.leftCols(2) = V_new;
	igl::writeOBJ("./test.obj",VV,F_new);
}

Eigen::MatrixXi locate_newvertex_to_old_mesh
        (
	        const Eigen::MatrixXd& P,
	        const Eigen::MatrixXi& F,
	        const Eigen::MatrixXd& V_uv
        )
{
	Eigen::MatrixXi Face_index = Eigen::MatrixXi::Zero(P.rows(),1);
	for (int i = 0; i < P.rows(); ++i)
	{
		/* code */
		double x, y;

		x = P(i,0);
		y = P(i,1);

		for (int j = 0; j < F.rows(); ++j)
		{
			double x1, y1, x2, y2, x3, y3;
			
			x1 = V_uv(F(j,0), 0);
			y1 = V_uv(F(j,0), 1);

			x2 = V_uv(F(j,1), 0);
			y2 = V_uv(F(j,1), 1);

			x3 = V_uv(F(j,2), 0);
			y3 = V_uv(F(j,2), 1);
			
			if (isInside(x1,y1,x2,y2,x3,y3,x,y))
			{
				Face_index.row(i) << j;
				break;
			}
			assert(j != F.rows()-1);
		}
	}
  	return Face_index;
}

Eigen::MatrixXd project_back_to_3d
        (
	        const Eigen::MatrixXd& P,
	        const Eigen::MatrixXi& F,
	        const Eigen::MatrixXd& V_uv,
	        const Eigen::MatrixXd& V, 
	        const Eigen::MatrixXi& Face_index
        )
{
	Eigen::MatrixXd V_new = Eigen::MatrixXd::Zero(P.rows(),3);
	Eigen::MatrixXd B   = Eigen::MatrixXd::Zero(3,3);
	Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(3,1);
	B.row(2) << 1,1,1;

	for (int i = 0; i < P.rows(); ++i)
	{
		double x1,y1,x2,y2,x3,y3;
		int faceid = Face_index(i,0);

		rhs.col(0) << P(i,0),P(i,1),1;
		x1 = V_uv( F(faceid,0), 0 );
		y1 = V_uv( F(faceid,0), 1 );

		x2 = V_uv( F(faceid,1), 0 );
		y2 = V_uv( F(faceid,1), 1 );

		x3 = V_uv( F(faceid,2), 0 );
		y3 = V_uv( F(faceid,2), 1 );

		B.row(0) << x1, x2, x3;
		B.row(1) << y1, y2, y3;

		Eigen::FullPivLU<Eigen::MatrixXd> lusolver(B);
        Eigen::MatrixXd lambda = lusolver.solve(rhs);

        V_new.row(i) = lambda(0,0)*V.row(F(faceid,0)) + 
        			   lambda(1,0)*V.row(F(faceid,1)) +
        			   lambda(2,0)*V.row(F(faceid,2));

	}
  	return V_new;
}