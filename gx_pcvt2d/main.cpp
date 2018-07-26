/*
 *  _____ _____ ________  _
 * /  __//  __//  __/\  \//
 * | |  _|  \  |  \   \  / 
 * | |_//|  /_ |  /_  /  \ 
 * \____\\____\\____\/__/\\
 *
 * Graphics Environment for EXperimentations.
 *  Copyright (C) 2006 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: 
 *
 *     ALICE Project - INRIA
 *     INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 */

#include <Geex/graphics/geexapp.h>
#include <Geex/basics/file_system.h>
#include <Geex/basics/stopwatch.h>
#include <glut_viewer/tweak_bar.h>
#include <fstream>
#include "cvt.h"
//#include "psa-src\psa.h"
//#include "irreg_editor.h"

void Lloyd() ;
void TW_CALL CB_action(void *clientData) ;
void TW_CALL CB_convert_mesh(void *clientData) ;
void TW_CALL CB_create_grid(void *clientData) ;

namespace Geex {

    class CVTApp : public GeexApp {
    public:
        CVTApp(int argc, char** argv) : GeexApp(argc, argv) { 
            hdr_ = false ;
            boundary_filename_ = get_file_arg("line") ;
            if(boundary_filename_.length() > 0) {
                if(!Geex::FileSystem::is_file(boundary_filename_)) {
                    boundary_filename_ = Geex::FileSystem::get_project_root() + 
                        "/gx_pcvt2d/" + boundary_filename_ ;
                }
            }
            nb_points_ = 100 ;
            get_arg("nb_pts", nb_points_) ;
            nb_iter_ = 50 ;
            get_arg("nb_iter", nb_iter_) ;            
            non_convex_ = GL_FALSE ;
            get_arg("non_convex", non_convex_) ;
            edit_ = GL_FALSE ;
			insert_boundary_ = GL_FALSE ;
			get_arg("insert_boundary", insert_boundary_) ;
			min_max_ = 30 ;
			nb_runs_ = 1 ;
        }

        CVT* cvt() { return static_cast<CVT*>(scene()) ; }

        GLboolean& edit() { return edit_ ; }

		virtual void init_scene() {
            scene_ = new CVT ;
            std::cerr << "Non convex = " 
                      << (non_convex_ ? "true" : "false") 
                      << "(use +non_convex to set)" << std::endl ;
            cvt()->set_non_convex_mode(non_convex_) ;
            if(boundary_filename_.length() > 0) {
                cvt()->load_boundary(boundary_filename_) ;
            }
			//cvt()->insert_boundary() = insert_boundary_ ;
            cvt()->insert_random_vertices(nb_points_) ;
			//insert_copies() ;
        }

		void insert_grid() {
			cvt()->clear() ;
			cvt()->insert_grid(500) ;
			insert_copies() ;
		}

        void Lloyd() {
            cvt()->lloyd(nb_iter_) ;
        }

        void Lloyd(int nb_iter) {
            cvt()->lloyd(nb_iter) ;
        }

        void NewtonLloyd() {
            cvt()->newton_lloyd(nb_iter_) ;
        }

		void Lloyd_fpo() {
			cvt()->lloyd_fpo(nb_iter_) ;
		}

        void reset() {
            cvt()->clear() ;
            cvt()->insert_random_vertices(nb_points_) ;
			insert_copies() ; 
        }

		void insert_copies() {
			//cvt()->insert_copies(cvt()->pvd_mode()==FULL_COPY) ;
		}

		void clear_copies() {
			//cvt()->clear_copies(true) ;
		}

		void save() {
			std::string filename = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/" + "out.pts" ;
			cvt()->save(filename) ;
		}
		void load() {
			std::string filename = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/" + "out.pts" ;
			cvt()->load(filename) ;
		}

		void save(const char* name) {
			std::string filename = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/" + name + ".pts" ;
			cvt()->save(filename) ;
		}
		void load(const char* name) {
			std::string filename = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/" + name + ".pts" ;
			cvt()->load(filename) ;
		}

		void generate_poisson_disk()  {
			std::string filename = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/pointsets/random/" ;
			if(nb_runs_ > 1) 
				cvt()->generate_pointset(cvt()->sample_radius(), nb_runs_, filename) ;
			else {
				cvt()->generate_poisson_disk() ;
			}
		}

		void analysis_pointsets() {
			std::string path = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/pointsets/PVoronoi/r0.005" ;
			cvt()->analyze_pointsets(path) ;
		}

		void smooth() {
			cvt()->smooth(nb_iter_, true) ;
		}

        virtual void init_gui() {
            GeexApp::init_gui() ;

            // New-style GUI =====================================================================================

            TwBar* graphics_bar = TwNewBar("Graphics") ;
			TwDefine(" Graphics position='16 10' size='200 480'") ;          
            TwAddVarRW(graphics_bar, "Vertices", TW_TYPE_FLOAT, &cvt()->vertices_size(), "min=0 max=1 step=0.01") ;
            TwAddVarRW(graphics_bar, "Vertex Color", TW_TYPE_BOOL8, &cvt()->vertices_color(), "") ;
            TwAddVarRW(graphics_bar, "Centers", TW_TYPE_FLOAT, &cvt()->centers_size(), "min=0 max=1 step=0.01") ;
            TwAddVarRW(graphics_bar, "Primal", TW_TYPE_BOOL8, &cvt()->show_primal_mesh(), "") ;
            TwAddVarRW(graphics_bar, "Dual", TW_TYPE_BOOL8, &cvt()->show_dual_mesh(), "") ;
			TwAddVarRW(graphics_bar, "InnerDual", TW_TYPE_BOOL8, &cvt()->show_inner_voronoi(), "") ;
			TwAddVarRW(graphics_bar, "Colorize", TW_TYPE_BOOL8, &cvt()->colorize(), "") ;
			TwAddVarRW(graphics_bar, "Non-hex", TW_TYPE_BOOL8, &cvt()->show_cells(), "") ;
			TwAddVarRW(graphics_bar, "Snap border", TW_TYPE_BOOL8, &cvt()->snap_boundary(), "") ;
			TwAddVarRW(graphics_bar, "Bndcells", TW_TYPE_BOOL8, &cvt()->show_boundary_cells(), "") ;
			TwAddVarRW(graphics_bar, "Euclidean", TW_TYPE_BOOL8, &cvt()->show_pvd_euclidean(), "") ;
			TwAddVarRW(graphics_bar, "Copies", TW_TYPE_BOOL8, &cvt()->show_copies(), "") ;
			TwEnumVal pvd_def[] = { {NO_COPY, "No Copy"}, {FULL_COPY, "Full Copy"}, {MIN_COPY, "Min Copy"}} ;
			TwType tw_pvd = TwDefineEnum("PVDMode", pvd_def, 3) ;
			TwAddVarRW(graphics_bar, "PVD Mode", tw_pvd, &cvt()->pvd_mode(), "") ;
			TwAddVarRW(graphics_bar, "show regularity", TW_TYPE_BOOL8, &cvt()->show_regularity(), "") ;            
            TwAddVarRW(graphics_bar, "Energy", TW_TYPE_BOOL8, &cvt()->show_energy(), "") ;
            TwAddVarRW(graphics_bar, "Bkgnd. field", TW_TYPE_BOOL8, &cvt()->show_field(), "") ;
            TwAddVarRW(graphics_bar, "Quads", TW_TYPE_FLOAT, &cvt()->quad_ratio(), "min=0.0 max=1.5 step=0.01") ;
            TwAddVarRW(graphics_bar, "show Lp", TW_TYPE_BOOL8, &cvt()->show_lp(), "") ;            
            TwAddVarRW(graphics_bar, "Lp view", TW_TYPE_INT32, &cvt()->lp_shader(), "min=0 max=10 step=1") ;
            TwAddVarRW(graphics_bar, "Lp scale", TW_TYPE_FLOAT, &cvt()->new_uniform("scale", 1.0), "min=0.001 max=10 step=0.1") ;
			TwAddVarRW(graphics_bar, "perturb", TW_TYPE_FLOAT, &cvt()->perturb(), "min=0 max=1 step=0.001") ;

			TwAddSeparator(graphics_bar, "Poisson sampling", "");
			TwAddVarRW(graphics_bar, "number of runs", TW_TYPE_INT32, &nb_runs_, "min=1 max=10000 step=1") ;
			TwEnumVal pds_def[] = { {PD_CCVT, "ccvt"}, {PD_VORONOI, "Voronoi"}, {PD_BOUNDARY, "Pure"}, {PD_DARTTHROW, "Dart throw"}, {PD_PENROSE, "Penrose"} } ;
            TwType tw_pds = TwDefineEnum("Sampling mode", pds_def, 5) ;
            TwAddVarRW(graphics_bar, "Sampling mode", tw_pds, &cvt()->sampling_mode(), "") ;
			TwAddVarRW(graphics_bar, "sample radius", TW_TYPE_DOUBLE, &cvt()->sample_radius(), "min=0 max=1 step=0.000001") ;
			TwAddVarRW(graphics_bar, "show radius", TW_TYPE_BOOL8, &cvt()->show_disk(), "") ;
			TwAddVarRW(graphics_bar, "show min/max", TW_TYPE_BOOL8, &cvt()->show_min_max(), "") ;
			TwAddVarRW(graphics_bar, "edge histogram", TW_TYPE_BOOL8, &cvt()->show_edge_hist(), "") ;

			TwBar* numerics_bar = TwNewBar("Numerics") ;
            TwDefine(" Numerics position='16 500' size='200 240'") ;            
			TwAddVarRW(numerics_bar, "nb iter",  TW_TYPE_FLOAT, &nb_iter_, "min=0 max=1000 step=1") ;
            TwAddVarRW(numerics_bar, "Lp order", TW_TYPE_FLOAT, &cvt()->Lp(), "min=0 max=9 step=1") ;
            TwEnumVal aniso_def[] = { {CONSTANT, "const"}, {R_INV, "1/R"}, {BORDER, "Border"} } ;
            TwType tw_aniso = TwDefineEnum("AnisoType", aniso_def, 3) ;
            TwAddVarRW(numerics_bar, "Aniso", tw_aniso, &cvt()->aniso_mode(), "") ;
            TwType tw_center_mode = TwDefineEnum("CenterMode", "centroid, quasi-incenter") ;
            TwAddVarRW(numerics_bar, "CenterMode", tw_center_mode, &cvt()->center_mode(), "") ;
            TwAddVarRW(numerics_bar, "X scale", TW_TYPE_FLOAT, &cvt()->Xscale(), "min=-2 max=2 step=1") ;            
            TwAddVarRW(numerics_bar, "Y scale", TW_TYPE_FLOAT, &cvt()->Yscale(), "min=-2 max=2 step=1") ;            
            TwEnumVal cvt_def[] = { {LLOYD, "Lloyd"}, {NEWTON, "Newton"}, {THETA, "Theta"} } ;
            TwType tw_cvt = TwDefineEnum("CVTType", cvt_def, 3) ;
            TwAddVarRW(numerics_bar, "Mode", tw_cvt, &cvt()->mode(), "") ;
            //TwAddVarRW(numerics_bar, "Bndry.", TW_TYPE_BOOL8, &cvt()->insert_boundary(), "") ;            
            TwAddVarRW(numerics_bar, "Edit", TW_TYPE_BOOL8, &edit_, "") ;    
		
            viewer_properties_->add_separator("Graphics") ;
            viewer_properties_->add_slider("Vertices", cvt()->vertices_size()) ;
            viewer_properties_->add_slider("Centers", cvt()->centers_size()) ;
            viewer_properties_->add_toggle("Primal mesh", cvt()->show_primal_mesh()) ;
            viewer_properties_->add_toggle("Dual mesh", cvt()->show_dual_mesh()) ;
			viewer_properties_->add_toggle("Colorize", cvt()->colorize()) ;
			viewer_properties_->add_toggle("Snap border", cvt()->snap_boundary()) ;
            viewer_properties_->add_toggle("Non-hex", cvt()->show_cells()) ;
            viewer_properties_->add_toggle("Energy", cvt()->show_energy()) ;
            viewer_properties_->add_toggle("Bkgnd. field", cvt()->show_field()) ;
            viewer_properties_->add_slider("Quads", cvt()->quad_ratio(), 0.1, 1.5) ;
			
            viewer_properties_->add_separator("Optimizer") ;
            viewer_properties_->add_slider("Lp order", cvt()->Lp(), 0.0, double(DelaunayCVT::MAX_P))->set_integer(GL_TRUE) ;
            viewer_properties_->add_enum("Aniso", cvt()->aniso_mode(), GlutViewerGUI::LabelList() | AnisoModeNames) ;            
            viewer_properties_->add_slider("X scale", cvt()->Xscale(), -2.0, 2.0)->set_integer(GL_TRUE) ;
            viewer_properties_->add_slider("Y scale", cvt()->Yscale(), -2.0, 2.0)->set_integer(GL_TRUE) ;
            viewer_properties_->add_enum("Mode", cvt()->mode(), GlutViewerGUI::LabelList() | CVTModeNames) ;
            //viewer_properties_->add_toggle("Bndry.", cvt()->insert_boundary()) ;
            viewer_properties_->add_toggle("Edit", edit_) ;

            toggle_skybox_CB() ;

            glut_viewer_add_toggle('B', glut_viewer_is_enabled_ptr(GLUT_VIEWER_BACKGROUND), "switch Color/BW") ;
//            glut_viewer_add_toggle('T', &viewer_properties_->visible(), "viewer properties") ;
            viewer_properties_->hide() ;
            glut_viewer_add_toggle('T', glut_viewer_is_enabled_ptr(GLUT_VIEWER_TWEAKBARS), "Tweak bars") ;

            glut_viewer_disable(GLUT_VIEWER_BACKGROUND) ;
        }

        virtual GLboolean mouse(float x, float y, int button, enum GlutViewerEvent event) {
            static int timestamp = 0 ;
            static int last_timestamp = 0 ;
            timestamp++ ;
            static int mode = 0 ;

            GLdouble p[3] ;
            GLdouble v[3] ;
            if(GeexApp::mouse(x, y, button, event)) { return GL_TRUE ; }
            if(edit_) {
 
				bool change = false ;
                if(event == GLUT_VIEWER_DOWN) { mode = button ; change = true ; }
                
                if(event == GLUT_VIEWER_UP) { return GL_FALSE ; }

                glut_viewer_get_picked_ray(p,v) ;
                vec2 pt(p[0], p[1]) ;

                if(mode == 0 && (change || timestamp - last_timestamp > 3)) {
                    if(cvt()->in_boundary(pt)) {
                        cvt()->begin_insert() ;
                        cvt()->insert(pt) ;
                        cvt()->end_insert() ;
                    }
                    last_timestamp = timestamp ;
                }

                if(mode == 1 && (change || timestamp - last_timestamp > 3)) {
                    if(cvt()->in_boundary(pt)) {
//						if(cvt()->period()) cvt()->clear_copies() ;
                        cvt()->begin_insert() ;
                        cvt()->remove(pt) ;
                        cvt()->end_insert(false) ;
                        cvt()->begin_insert() ;
                        cvt()->insert(pt) ; // ->locked = true ;
                        cvt()->end_insert() ;
//						if(cvt()->period()) cvt()->insert_copies(cvt()->pvd_mode()==FULL_COPY, false) ;
                    }
                    last_timestamp = timestamp ;
                }

                if(mode == 2 && (change || timestamp - last_timestamp > 3)) {
                    if(cvt()->in_boundary(pt)) {
//						if(cvt()->period()) cvt()->clear_copies() ;
                        cvt()->begin_insert() ;
                        cvt()->remove(pt) ;
                        cvt()->end_insert() ;
//						if(cvt()->period()) cvt()->insert_copies(cvt()->pvd_mode()==FULL_COPY) ;
                    }
                    last_timestamp = timestamp ;
                }

                return GL_TRUE ;
			} 
			// pick vertex, edge
			if(cvt()->map_edit_mode()) {
				bool change = false ;
                if(event == GLUT_VIEWER_DOWN) { mode = button ; change = true ; }
                if(event == GLUT_VIEWER_UP) { return GL_FALSE ; }

                glut_viewer_get_picked_ray(p,v) ;
                vec2 pt(p[0], p[1]) ;

			}
            
			return GL_FALSE ;
        }

    private:
        std::string boundary_filename_ ;
        int nb_points_ ;
		int min_max_ ;
        GLfloat nb_iter_ ;
        GLboolean non_convex_ ;
        GLboolean edit_ ;
		GLboolean insert_boundary_ ;
		int nb_runs_ ;
    } ;
}

Geex::CVTApp* cvt_app() { return static_cast<Geex::CVTApp*>(Geex::GeexApp::instance()) ; }

void Lloyd() {
    cvt_app()->Lloyd() ;
}

void Lloyd1() {
    cvt_app()->Lloyd(1) ;
}

void NewtonLloyd() {
    cvt_app()->NewtonLloyd() ;
	cvt_app()->cvt()->lloyd_energy() ; // update regularity
}

void reset() {
    cvt_app()->reset() ;
	cvt_app()->cvt()->lloyd_energy() ;
}

void insert_copies() {
	cvt_app()->insert_copies() ;
	cvt_app()->cvt()->lloyd_energy() ;
}

void clear_copies() {
	cvt_app()->clear_copies() ;
}

void euclidean() {
	cvt_app()->cvt()->show_pvd_euclidean() = !cvt_app()->cvt()->show_pvd_euclidean() ;
}
void inc_Lp() {
    cvt_app()->cvt()->Lp() += 1.0 ;
    cvt_app()->cvt()->Lp() = Geex::gx_min(cvt_app()->cvt()->Lp(), float(Geex::DelaunayCVT::MAX_P)) ;
    cvt_app()->cvt()->lp_shader() = int(cvt_app()->cvt()->Lp()) ;
}

void dec_Lp() {
    cvt_app()->cvt()->Lp() -= 1.0 ;
    cvt_app()->cvt()->Lp() = Geex::gx_max(cvt_app()->cvt()->Lp(), 0.0f) ;
    cvt_app()->cvt()->lp_shader() = int(cvt_app()->cvt()->Lp()) ;
}

void save() {
	cvt_app()->save() ;
}

void load() {
	cvt_app()->load() ;
	cvt_app()->cvt()->lloyd_energy() ; // update regularity
}

void stats() { // statistic 
//	cvt_app()->stats() ;
//	cvt_app()->degree_stats(100) ;
	cvt_app()->analysis_pointsets() ;
}

void show_edge_length() {
	cvt_app()->cvt()->compute_edge_length() ;
}

void do_perturb() {
	cvt_app()->cvt()->do_perturb() ;
}

void compute_hessian() {
	cvt_app()->cvt()->compute_hessian() ;
}

void generate_poisson_disk() {
	cvt_app()->generate_poisson_disk() ;
	//std::string fname = "D:/project/gx_pcvt/gx_pcvt2d/pointsets/PVoronoi" ;
	//Geex::psa(true, false, false, fname) ;
}

void update_gapss() {
	cvt_app()->cvt()->update_gaps(cvt_app()->cvt()->sample_radius()) ;
}

void batch_snapshot() {
	/*std::string dir = Geex::FileSystem::get_project_root() + "/gx_pcvt2d/" ;
	char filename[1024] ;
	int N = 6, nSol = 5 ;
	for(int i=0; i<nSol; ++i) {
		sprintf(filename, "%sdata/6/pcvt%d_%d.pts", dir.c_str(), N, i+1) ;
		cvt_app()->cvt()->load(filename) ;
		cvt_app()->cvt()->insert_copies(true) ;
		glut_viewer_redraw() ;
		sprintf(filename, "%sdata/6/pcvt%d_%d.png", dir.c_str(), N, i+1) ;
		cvt_app()->cvt()->snapshot(filename) ;
	}*/
}

void Lloyd_fpo() {
	cvt_app()->Lloyd_fpo() ;
}

void insert_grid() {
	cvt_app()->insert_grid() ;
}

int main(int argc, char** argv) {
    Geex::initialize() ;
    Geex::CVTApp app(argc, argv) ;
    glut_viewer_add_key_func('k', Lloyd, "Lloyd iterations") ;
    glut_viewer_add_key_func('K', Lloyd1, "Lloyd one iteration") ;
	glut_viewer_add_key_func('f', Lloyd_fpo, "Lloyd fpo iterations") ;
    glut_viewer_add_key_func('m', NewtonLloyd, "Newton-Lloyd iterations") ;
	//glut_viewer_add_key_func('M', smooth, "smooth point set") ;
    glut_viewer_add_key_func('Z', reset, "reset") ;
	glut_viewer_add_key_func('g', insert_grid, "insert_grid") ;
	glut_viewer_add_key_func('j', compute_hessian, "Hessian") ;
	//glut_viewer_add_key_func('i', insert_copies, "Insert copies") ;
	//glut_viewer_add_key_func('c', clear_copies, "clear copies") ;
	glut_viewer_add_key_func('u', euclidean, "show euclidean") ;
	glut_viewer_add_key_func('s', save, "save points") ;
	glut_viewer_add_key_func('o', load, "load points") ;
	glut_viewer_add_key_func('d', stats, "statistic") ;
	glut_viewer_add_key_func('z', generate_poisson_disk, "generate poisson disk") ;
	glut_viewer_add_key_func('v', update_gapss, "update void regions") ;
	//glut_viewer_add_key_func('b', batch_snapshot, "snapshots") ;
	glut_viewer_add_key_func('p', do_perturb, "perturb") ;
    glut_viewer_add_toggle('e', &(cvt_app()->edit()), "edit") ;
    glut_viewer_add_key_func('<', dec_Lp, "decrement Lp") ;
    glut_viewer_add_key_func('>', inc_Lp, "increment Lp") ;
    glut_viewer_disable(GLUT_VIEWER_3D) ;
    app.main_loop() ;
    Geex::terminate() ;
}
