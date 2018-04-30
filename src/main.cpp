#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include <imgui/imgui.h>
#include "geometry_maps.h"
#include "halftone.h"
#include "vertex_to_triangle.h"

/* key == 1 */
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/lscm.h>

#include <igl/colon.h>

#include <igl/jet.h>
#include <igl/parula.h>
#include <igl/png/writePNG.h>


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
Eigen::MatrixXd P;

Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> texture_I;

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
	if (key == '1')
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

	if (key == '2')
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
	if (key == '3')
	{
		area_distortion = area_distortion_map(V,F,V_uv);
		
		viewer.data().clear();
	    viewer.data().set_mesh(V_uv,F);
	    viewer.data().set_uv(V_uv);
	    viewer.core.align_camera_center(V_uv,F);

		MatrixXd C;
		// cout << area_distortion.minCoeff()<< endl;
		// cout << area_distortion.maxCoeff()<< endl;
    	igl::parula(area_distortion,false,C);
    	viewer.data().set_colors(C);
	}

	if (key == '4')
	{
		K_gaussian = gaussian_curvature_map(V,F,V_uv);
		
		viewer.data().clear();
	    viewer.data().set_mesh(V_uv,F);
	    viewer.data().set_uv(V_uv);
	    viewer.core.align_camera_center(V_uv,F);

		MatrixXd C;
		// cout << K_gaussian.minCoeff()<< endl;
		// cout << K_gaussian.maxCoeff()<< endl;

    	igl::parula(K_gaussian,false,C);
    	viewer.data().set_colors(C);
	}

	if (key == '5')
	{
		H = mean_curvature_map(V,F,V_uv);
		
		viewer.data().clear();
	    viewer.data().set_mesh(V_uv,F);
	    viewer.data().set_uv(V_uv);
	    viewer.core.align_camera_center(V_uv,F);

		MatrixXd C;
		// cout << H.minCoeff()<< endl;
		// cout << H.maxCoeff()<< endl;

    	igl::parula(H,false,C);
    	viewer.data().set_colors(C);
	}

	if (key == '6')
	{
		F_index = igl::colon<int,int,int>(0,F.rows()-1);
		cout << F_index << endl;

		viewer.data().clear();
	    viewer.data().set_mesh(V_uv,F);
	    viewer.data().set_uv(V_uv);
	    viewer.core.align_camera_center(V_uv,F);

		MatrixXd C;
    	igl::parula(F_index,true,C);
    	viewer.data().set_colors(C);
/*
    	// Allocate temporary buffers for 1280x800 image
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);

		// Draw the scene in the buffers
		viewer.core.draw_buffer(
      		viewer.data(),false,R,G,B,A);

		// Save it to a PNG
		igl::png::writePNG(R,G,B,A,"out.png");
*/
	}

	if (key == '7')
	{
		P = error_diffusion(V,F,V_uv,area_distortion,bnd_uv);
		
		// viewer.data().clear();
        viewer.core.align_camera_center(P);
        viewer.data().point_size = 5;
        viewer.data().add_points(P, Eigen::RowVector3d(1,0,0));
	}

	if (key == '8')
	{
		F_new = vertex_to_triangle(P);
		viewer.data().clear();
		viewer.data().set_mesh(P, F_new);
  		// viewer.data().set_uv(P);
	}

	if (key == '9')
	{
		F_index = locate_newvertex_to_old_mesh(P,F,V_uv);
		viewer.data().clear();
		viewer.data().set_mesh(P, F_new);

		MatrixXd C;
    	igl::parula(F_index,false,C);
    	viewer.data().set_colors(C);
	}

	if (key == 'a')
	{
		V_new = project_back_to_3d(P,F,V_uv,V,F_index);
		viewer.data().clear();
		viewer.data().set_mesh(V_new, F_new);
	}

	return false;
}

int main(int argc, char *argv[])
{
	using namespace std;
	using namespace Eigen;

	// Load a mesh in OBJ format
	igl::readOFF("../data/camel_head.off", V, F);
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
		ImGui::Begin(
	    	"Parametrization", nullptr,
	    	ImGuiWindowFlags_NoSavedSettings
		);

		if (ImGui::Button("Harmonic", ImVec2(-1,0)))
			key_down(viewer, '1', 0);
		
		if (ImGui::Button("LSCM", ImVec2(-1,0)))
			key_down(viewer, '2', 0);
		ImGui::End();

		ImGui::SetNextWindowPos(ImVec2(400.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
	    	"Geometry maps", nullptr,
	    	ImGuiWindowFlags_NoSavedSettings
		);

		if (ImGui::Button("Area distortion", ImVec2(-1,0))){
			key_down(viewer, '3', 0);
		}
		if (ImGui::Button("Gaussian curvature", ImVec2(-1,0))){
			key_down(viewer, '4', 0);
		}
		if (ImGui::Button("Mean curvature", ImVec2(-1,0))){
			key_down(viewer, '5', 0);
		}
		if (ImGui::Button("Face index", ImVec2(-1,0))){
			key_down(viewer, '6', 0);
		}
		ImGui::End();

		ImGui::SetNextWindowPos(ImVec2(620.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
	    	"Halftone", nullptr,
	    	ImGuiWindowFlags_NoSavedSettings
		);

		if (ImGui::Button("error_diffusion", ImVec2(-1,0))){
			key_down(viewer, '7', 0);
		}
		ImGui::End();

		ImGui::SetNextWindowPos(ImVec2(840.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
	    	"Remesh", nullptr,
	    	ImGuiWindowFlags_NoSavedSettings
		);

		if (ImGui::Button("build_triangle", ImVec2(-1,0))){
			key_down(viewer, '8', 0);
		}

		if (ImGui::Button("locate_vertex_to_old_mesh", ImVec2(-1,0))){
			key_down(viewer, '9', 0);
		}

		if (ImGui::Button("project_back_to_3d", ImVec2(-1,0))){
			key_down(viewer, 'a', 0);
		}
		
		ImGui::End();
	};

	viewer.callback_key_down = &key_down;
	// Launch the viewer
	viewer.data().set_mesh(V, F);
	viewer.data().set_texture(texture_I, texture_I, texture_I);
	viewer.launch();
}
