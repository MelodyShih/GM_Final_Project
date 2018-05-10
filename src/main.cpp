#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <imgui/imgui.h>
#include "geometry_maps.h"
#include "halftone.h"
#include "vertex_to_triangle.h"
#include "geometric_opt.h"
#include "connectivity_opt.h"

/* key == 1 */
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/lscm.h>

#include <igl/colon.h>

#include <igl/jet.h>
#include <igl/parula.h>
#include <igl/png/writePNG.h>

#include <vector>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/edge_topology.h>

// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXd V_new;
Eigen::MatrixXi F;
Eigen::MatrixXi F_new;
Eigen::MatrixXi F_index;
Eigen::MatrixXd V_uv;
Eigen::VectorXi bnd;
Eigen::MatrixXd bnd_uv;

// Geometric maps
Eigen::VectorXd area_distortion;
Eigen::VectorXd K_gaussian;
Eigen::VectorXd H;

// Halftone
Eigen::MatrixXd test;
// Input: imported points, #P x3
Eigen::MatrixXd P1, P2;
Eigen::MatrixXd P;
std::vector<Eigen::RowVector2d> p;
Eigen::MatrixXd P_new;
Eigen::MatrixXi E_new;

Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> texture_I;
float r_curvature, r_area;

void line_texture() 
{
  int size = 128;              // Texture size
  int w    = 7;                // Line width
  int pos  = size / 2 - w / 2; // Center the line
  texture_I.setConstant(size, size, 255);
  texture_I.block(0, pos, size, w).setZero();
  texture_I.block(pos, 0, w, size).setZero();

}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
	using namespace Eigen;
	using namespace std;

	/* parameterization */
	if (key == 'a')
	{
	    /* Harmonic mapping */ 
	    igl::boundary_loop(F,bnd);

	    igl::map_vertices_to_circle(V,bnd,bnd_uv);
	    igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);

	    /* Scale the uv */
	    // V_uv *= 10;
	    
	    viewer.data().clear();
	    viewer.data().set_mesh(V_uv,F);
	    viewer.data().set_uv(V_uv);
	    viewer.core.align_camera_center(V_uv,F);
	}

	if (key == 'b')
	{
		// Fix two points on the boundary
		VectorXi b(2,1);
		igl::boundary_loop(F,bnd);
		b(0) = bnd(0);
		b(1) = bnd(round(bnd.size()/2));
		MatrixXd bc(2,2);
		bc<<0,0,1,0;

		// LSCM parametrization
		igl::lscm(V,F,b,bc,V_uv);
		/* Scale the uv */
	    // V_uv *= 10;

		viewer.data().clear();
		viewer.data().set_mesh(V_uv, F);
  		viewer.data().set_uv(V_uv);
  		viewer.core.align_camera_center(V_uv,F);
	}

	/* visualize geometric maps */
	if (key == 'c')
	{
		area_distortion = area_distortion_map(V,F,V_uv);
		
		viewer.data().clear();
	    viewer.data().set_mesh(V_uv,F);
	    viewer.data().set_uv(V_uv);
	    viewer.core.align_camera_center(V_uv,F);

		MatrixXd C;
		double max, min, range;
		max = area_distortion.maxCoeff();
		min = area_distortion.minCoeff();
		range = max - min;

		// rescale to [0,1]
		area_distortion = (area_distortion - min*VectorXd::Ones(area_distortion.size()))/(range);
		for (int i = 0; i < area_distortion.size(); ++i)
		{
			area_distortion(i) = pow(area_distortion(i), r_area);
		}
		
    	igl::jet(area_distortion,false,C);
    	viewer.data().set_colors(C);

    	Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);

		// Draw the scene in the buffers
		viewer.data().show_lines = false;
		viewer.core.draw_buffer(viewer.data(),false,R,G,B,A);

		// Save it to a PNG
		igl::png::writePNG(R,G,B,A,"/Users/yuhsuan/Downloads/ImageHalftoningFloyd/ImageHalftoningFloyd/area_map.png");
	}

	if (key == 'd')
	{
		K_gaussian = gaussian_curvature_map(V,F,V_uv);
		
		viewer.data().clear();
	    viewer.data().set_mesh(V_uv,F);
	    viewer.data().set_uv(V_uv);
	    viewer.core.align_camera_center(V_uv,F);

		MatrixXd C;
		double max, min, range;
		max = K_gaussian.maxCoeff();
		min = K_gaussian.minCoeff();
		range = max - min;

		K_gaussian = K_gaussian/(range);

    	igl::jet(K_gaussian,false,C);
    	viewer.data().set_colors(C);

    	Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(640,400);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(640,400);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(640,400);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(640,400);

		// Draw the scene in the buffers
		viewer.data().show_lines = false;
		viewer.core.draw_buffer(viewer.data(),false,R,G,B,A);

		// Save it to a PNG
		igl::png::writePNG(R,G,B,A,"/Users/yuhsuan/Downloads/ImageHalftoningFloyd/ImageHalftoningFloyd/gaussian_curvature_map.png");
	}

	if (key == 'e')
	{
		H = mean_curvature_map(V,F,V_uv);
		
		viewer.data().clear();
	    viewer.data().set_mesh(V_uv,F);
	    viewer.data().set_uv(V_uv);
	    viewer.core.align_camera_center(V_uv,F);

		MatrixXd C;
		double max, min, range, median;
		
		// rescale to [0,1]

		for (int i = 0; i < H.size(); ++i)
		{
			H(i) = abs(H(i));
		}

		max = H.maxCoeff();
		min = H.minCoeff();
		range = max - min;

		H = H/(range);

		for (int i = 0; i < H.size(); ++i)
		{
			H(i) = pow(abs(H(i)),r_curvature);
		}

    	igl::jet(H,false,C);
    	viewer.data().set_colors(C);

    	viewer.data().set_colors(C);

    	Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);

		// Draw the scene in the buffers
		viewer.data().show_lines = false;
		viewer.core.draw_buffer(viewer.data(),false,R,G,B,A);

		// Save it to a PNG
		igl::png::writePNG(R,G,B,A,"/Users/yuhsuan/Downloads/ImageHalftoningFloyd/ImageHalftoningFloyd/mean_curvature_map.png");
	}

	if (key == 'g')
	{
		level0_diffusion(V,F,V_uv,area_distortion,bnd_uv,p);
		
		int nrows = p.size();
		cout << nrows << endl;
		P1.resize(nrows, 2);
		for (int i = 0; i < nrows; ++i)
		{
			P1.row(i) = p[i];
		}
		viewer.data().clear();
        viewer.core.align_camera_center(V_uv);
        viewer.data().point_size = 5;
        viewer.data().add_points(P1, Eigen::RowVector3d(1,0,0));
	}
	
	if (key == 'r')
	{
		level1_diffusion(V,F,V_uv,area_distortion,bnd_uv,p);
		
		int nrows = p.size();
		cout << nrows << endl;
		P1.resize(nrows, 2);
		for (int i = 0; i < nrows; ++i)
		{
			P1.row(i) = p[i];
		}
		viewer.data().clear();
        viewer.core.align_camera_center(V_uv);
        viewer.data().point_size = 5;
        viewer.data().add_points(P1, Eigen::RowVector3d(1,0,0));
	}
	
	if (key == 's')
	{
		level2_diffusion(V,F,V_uv,area_distortion,bnd_uv,p);
		
		int nrows = p.size();
		cout << nrows << endl;
		P1.resize(nrows, 2);
		for (int i = 0; i < nrows; ++i)
		{
			P1.row(i) = p[i];
		}
		viewer.data().clear();
        viewer.core.align_camera_center(V_uv);
        viewer.data().point_size = 5;
        viewer.data().add_points(P1, Eigen::RowVector3d(1,0,0));
	}

	if (key == 'n')
	{
		P2 = error_diffusion(V,F,V_uv,H,bnd_uv);
		viewer.data().clear();
        viewer.core.align_camera_center(V_uv);
        viewer.data().point_size = 5;
        viewer.data().add_points(P2, Eigen::RowVector3d(1,0,0));
	}

	if (key == 'o')
	{
		P.resize(P1.rows()+P2.rows(), P1.cols());
		P << P1, P2;
		
		viewer.data().clear();
        viewer.core.align_camera_center(P);
        viewer.data().point_size = 5;
        viewer.data().add_points(P, Eigen::RowVector3d(1,0,0));
	}


	if (key == 'h')
	{
		vertex_to_triangle(P,F,V_uv,F_new,P_new);
		viewer.data().clear();
		viewer.data().set_mesh(P_new, F_new);
		viewer.core.align_camera_center(P_new,F_new);
	}

	if (key == 'i')
	{
		F_index = locate_newvertex_to_old_mesh(P_new,F,V_uv);
		viewer.data().clear();
		viewer.data().set_mesh(P_new, F_new);
	}

	if (key == 'j')
	{
		V_new = project_back_to_3d(P_new,F,V_uv,V,F_index);
		viewer.data().clear();
		viewer.data().set_mesh(V_new,F_new);
		viewer.core.align_camera_center(V_new,F_new);
	}

	if (key == 'k')
	{
		MatrixXd dp;
		dp = weighted_laplacian_flow(V_new,F_new,P_new);

		P_new = P_new + dp;

		viewer.data().clear();
		viewer.data().set_mesh(P_new, F_new);
		viewer.core.align_camera_center(P_new,F_new);
	}

	if (key == 'l')
	{
		cout << "Before" << endl;
		regularity_statistic(P_new,F_new);
		regularity_opt(P_new,F_new);
		cout << "After" << endl;
		regularity_statistic(P_new,F_new);
		viewer.data().set_mesh(P_new, F_new);
		viewer.core.align_camera_center(P_new,F_new);
	}

	if (key == 'm')
	{
		face_aspect_ratio_opt(V_new,F_new,P_new);
		viewer.data().set_mesh(P_new, F_new);
		viewer.core.align_camera_center(P_new,F_new);
	}


	return false;
}

int main(int argc, char *argv[])
{
	using namespace std;
	using namespace Eigen;

	// Load a mesh in OBJ format
	igl::readOBJ(argv[1], V, F);
	igl::opengl::glfw::Viewer viewer;

	line_texture();

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::InputFloat("r_curvature", &r_curvature, 0, 0, 3);
		ImGui::InputFloat("r_area", &r_area, 0, 0, 3);

		ImGui::Begin(
	    	"Parametrization", nullptr,
	    	ImGuiWindowFlags_NoSavedSettings
		);

		if (ImGui::Button("Harmonic", ImVec2(-1,0)))
			key_down(viewer, 'a', 0);

		ImGui::End();

		ImGui::SetNextWindowPos(ImVec2(400.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
	    	"Geometry maps", nullptr,
	    	ImGuiWindowFlags_NoSavedSettings
		);

		if (ImGui::Button("Area distortion", ImVec2(-1,0))){
			key_down(viewer, 'c', 0);
		}
		if (ImGui::Button("Gaussian curvature", ImVec2(-1,0))){
			key_down(viewer, 'd', 0);
		}
		if (ImGui::Button("Mean curvature", ImVec2(-1,0))){
			key_down(viewer, 'e', 0);
		}
		ImGui::End();

		ImGui::SetNextWindowPos(ImVec2(620.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
	    	"Halftone", nullptr,
	    	ImGuiWindowFlags_NoSavedSettings
		);

		if (ImGui::Button("resampling level 0", ImVec2(-1,0))){
			key_down(viewer, 'g', 0);
		}

		if (ImGui::Button("resampling level 1", ImVec2(-1,0))){
			key_down(viewer, 'r', 0);
		}

		if (ImGui::Button("resampling level 2", ImVec2(-1,0))){
			key_down(viewer, 's', 0);
		}

		if (ImGui::Button("boundary points", ImVec2(-1,0))){
			key_down(viewer, 'n', 0);
		}

		if (ImGui::Button("add points", ImVec2(-1,0))){
			key_down(viewer, 'o', 0);
		}

		if (ImGui::Button("clean points", ImVec2(-1,0))){
			p.clear();
			P1.resize(0,2);
			viewer.data().clear();
		}
		ImGui::End();

		ImGui::SetNextWindowPos(ImVec2(840.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
	    	"Remesh", nullptr,
	    	ImGuiWindowFlags_NoSavedSettings
		);

		if (ImGui::Button("build_triangle", ImVec2(-1,0))){
			key_down(viewer, 'h', 0);
		}

		if (ImGui::Button("locate_vertex_to_old_mesh", ImVec2(-1,0))){
			key_down(viewer, 'i', 0);
		}

		if (ImGui::Button("project_back_to_3d", ImVec2(-1,0))){
			key_down(viewer, 'j', 0);
		}
		ImGui::End();

		ImGui::SetNextWindowPos(ImVec2(1040.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
	    	"Optimization", nullptr,
	    	ImGuiWindowFlags_NoSavedSettings
		);

		if (ImGui::Button("geometric opt", ImVec2(-1,0))){
			key_down(viewer, 'k', 0);
		}

		if (ImGui::Button("regularity opt", ImVec2(-1,0))){
			key_down(viewer, 'l', 0);
		}

		if (ImGui::Button("face aspect ratio opt", ImVec2(-1,0))){
			key_down(viewer, 'm', 0);
		}
		
		ImGui::End();
	};

	viewer.callback_key_down = &key_down;
	// Launch the viewer
	viewer.data().set_mesh(V, F);
	viewer.data().set_texture(texture_I, texture_I, texture_I);
	viewer.launch();
}
