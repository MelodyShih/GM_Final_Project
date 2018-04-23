#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include "geometry_maps.h"
#include "halftone.h"

/* key == 1 */
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/lscm.h>

#include <igl/jet.h>
#include <igl/parula.h>


// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V_uv;

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
	    Eigen::VectorXi bnd;
	    igl::boundary_loop(F,bnd);

	    Eigen::MatrixXd bnd_uv;
	    igl::map_vertices_to_circle(V,bnd,bnd_uv);
	    igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);

	    /* Scale the uv */
	    V_uv *= 10;
	    
	    viewer.data().clear();
	    viewer.data().set_mesh(V_uv,F);
	    viewer.data().set_uv(V_uv);
	    viewer.core.align_camera_center(V_uv,F);
	}

	if (key == '2')
	{
		// Fix two points on the boundary
		VectorXi bnd,b(2,1);
		igl::boundary_loop(F,bnd);
		b(0) = bnd(0);
		b(1) = bnd(round(bnd.size()/2));
		MatrixXd bc(2,2);
		bc<<0,0,1,0;

		// LSCM parametrization
		igl::lscm(V,F,b,bc,V_uv);

		/* Scale the uv */
	    V_uv *= 10;

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
		cout << area_distortion.minCoeff()<< endl;
		cout << area_distortion.maxCoeff()<< endl;
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
		cout << K_gaussian.minCoeff()<< endl;
		cout << K_gaussian.maxCoeff()<< endl;

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
		cout << H.minCoeff()<< endl;
		cout << H.maxCoeff()<< endl;

    	igl::parula(H,false,C);
    	viewer.data().set_colors(C);
	}

	if (key == '6')
	{
		P = error_diffusion(V,F,V_uv,area_distortion);

		// viewer.data().clear();
        viewer.core.align_camera_center(P);
        viewer.data().point_size = 5;
        viewer.data().add_points(P, Eigen::RowVector3d(1,0,0));
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
		ImGui::End();

		ImGui::SetNextWindowPos(ImVec2(620.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
	    	"Halftone", nullptr,
	    	ImGuiWindowFlags_NoSavedSettings
		);

		if (ImGui::Button("error_diffusion", ImVec2(-1,0))){
			key_down(viewer, '6', 0);
		}
		
		ImGui::End();
	};

	viewer.callback_key_down = &key_down;
	// Launch the viewer
	viewer.data().set_mesh(V, F);
	viewer.data().set_texture(texture_I, texture_I, texture_I);
	viewer.launch();
}
