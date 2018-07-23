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

#include "delaunay_cvt.h"
#include "geometry.h"
#include "lloyd_energy.h"
#include <glut_viewer/glut_viewer.h>
#include <TNT/tnt_linalg.h>

#include "generated/L.h"

#include "generated/L2.h"
#include "generated/L4.h"
#include "generated/L6.h"
#include "generated/L8.h"
#include "generated/L10.h"
#include "generated/L12.h"
#include "generated/L14.h"
#include "generated/L16.h"
#include "generated/L18.h"
#include "generated/L20.h"

#include "generated/Ltheta2.h"
#include "generated/Ltheta4.h"
#include "generated/Ltheta6.h"
#include "generated/Ltheta8.h"
#include "generated/Ltheta10.h"
#include "generated/Ltheta12.h"
#include "generated/Ltheta14.h"
#include "generated/Ltheta16.h"
#include "generated/Ltheta18.h"
#include "generated/Ltheta20.h"

#include "generated/Ltheta_rho2.h"
#include "generated/Ltheta_rho4.h"
#include "generated/Ltheta_rho6.h"
#include "generated/Ltheta_rho8.h"
#include "generated/Ltheta_rho10.h"
#include "generated/Ltheta_rho12.h"
#include "generated/Ltheta_rho14.h"
#include "generated/Ltheta_rho16.h"
#include "generated/Ltheta_rho18.h"
#include "generated/Ltheta_rho20.h"

static bool lbfgs_redraw = true ;

void funcgrad_cvt2d(int N, double* x, double& f, double* g);
void newiteration_cvt2d(int N, const double* x, double f, const double* g, double gnorm);


namespace Geex {

    DelaunayCVT* DelaunayCVT::instance_ = nil ;

    DelaunayCVT::DelaunayCVT(Delaunay* delaunay) : delaunay_(delaunay) {
        gx_assert(instance_ == nil) ;
        instance_ = this ;
        symbolic_ = false ;

//      funcs_table_[0] = new GenericLloydFuncs<L> ;

        funcs_table_[0] = new GenericLloydFuncs<L2> ;
        funcs_table_[1] = new GenericLloydFuncs<L4> ;
        funcs_table_[2] = new GenericLloydFuncs<L6> ;
        funcs_table_[3] = new GenericLloydFuncs<L8> ;
        funcs_table_[4] = new GenericLloydFuncs<L10> ;
        funcs_table_[5] = new GenericLloydFuncs<L12> ;
        funcs_table_[6] = new GenericLloydFuncs<L14> ;
        funcs_table_[7] = new GenericLloydFuncs<L16> ;
        funcs_table_[8] = new GenericLloydFuncs<L18> ;
        funcs_table_[9] = new GenericLloydFuncs<L20> ;

        funcs_table_theta_[0] = new GenericLloydFuncs<Ltheta_rho2>(2) ;
        funcs_table_theta_[1] = new GenericLloydFuncs<Ltheta_rho4>(2) ;
        funcs_table_theta_[2] = new GenericLloydFuncs<Ltheta_rho6>(2) ;
        funcs_table_theta_[3] = new GenericLloydFuncs<Ltheta_rho8>(2) ;
        funcs_table_theta_[4] = new GenericLloydFuncs<Ltheta_rho10>(2) ;
        funcs_table_theta_[5] = new GenericLloydFuncs<Ltheta_rho12>(2) ;
        funcs_table_theta_[6] = new GenericLloydFuncs<Ltheta_rho14>(2) ;
        funcs_table_theta_[7] = new GenericLloydFuncs<Ltheta_rho16>(2) ;
        funcs_table_theta_[8] = new GenericLloydFuncs<Ltheta_rho18>(2) ;
        funcs_table_theta_[9] = new GenericLloydFuncs<Ltheta_rho20>(2) ;


        funcs_ = funcs_table_[0] ; 

        Lp_ = 0.0 ; 
        X_scale_ = 0.0 ;
        Y_scale_ = 0.0 ;

        use_theta_ = GL_FALSE ;
        mode_ = NEWTON ;
        aniso_mode_ = CONSTANT ;

        center_mode_ = CENTROID;
		snap_boundary_ = GL_FALSE ;

		pvd_mode_ = FULL_COPY ;
    }

    DelaunayCVT::~DelaunayCVT() {
		for(int i=0; i<MAX_P+1; i++){
			delete funcs_table_[i];
			delete funcs_table_theta_[i];
		}
		instance_ = nil ; 
    }

	void DelaunayCVT::smooth(int nb_iter, bool redraw) {

	}

	void DelaunayCVT::get_vertex_center(Delaunay::Vertex_handle v, vec2& g) {

	}

	static double dxx(double xix, double xkx, double p1x, double p2x) {
		double dxx ;
		dxx = (p1x*p1x + p2x*p2x + p1x*p2x)/3.0 - (xix*p1x+xix*p2x+p1x*xkx+p2x*xkx)/2.0 + xix*xkx;
		return dxx ;
	}
	
	static double dxy(const vec2& xi, const vec2& xk, const vec2& p1, const vec2& p2) {
		double dxy ;
		dxy = (p1[0]*p1[1]+p2[0]*p2[1])/3.0 + (p1[0]*p2[1]+p1[1]*p2[0])/6.0 - (xi[0]*p1[1]+xi[0]*p2[1]+p1[0]*xk[1]+p2[0]*xk[1])/2.0 + xi[0]*xk[1] ;
		return dxy ;
	}

	static double dyx(const vec2& xi, const vec2& xk, const vec2& p1, const vec2& p2) {
		double dyx ;
		dyx = (p1[0]*p1[1]+p2[0]*p2[1])/3.0 + (p1[0]*p2[1]+p1[1]*p2[0])/6.0 - (xi[1]*p1[0]+xi[1]*p2[0]+xk[0]*p1[1]+xk[0]*p2[1])/2.0 + xi[1]*xk[0] ;
		return dyx ;
	}

	void DelaunayCVT::compute_hessian() {

	}

    void DelaunayCVT::lloyd(int nb_iter, bool redraw) {
        for(unsigned int k=0; k<nb_iter; k++) {
            std::vector<vec2> new_points ;
			std::vector<bool> flags ;
			

			//delaunay_->compute_rvd() ;


            FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
/*				if(v->locked) {
					new_points.push_back(to_geex(v->point())) ;
					flags.push_back(true) ;
				}
				else */{ //if(delaunay_->in_boundary(to_geex(v->point()))) {
					double V ;
					vec2 g ;
					if (center_mode_ == CENTROID) {
						if(period()) {
							if(!delaunay_->is_primary(v))
								continue ;
							get_cell_primary_centroid(v, g, V) ;
						}
						else if(snap_boundary() && delaunay_->is_boundary_cell(v))
							get_boundary_cell_centroid(v, g, V) ;
						else {
							get_cell_centroid(v, g, V) ;
						//	gx_assert(delaunay_->in_boundary(g)) ;
						}
					}
					else if(center_mode_ == QUASI_INCENTER)
						get_cell_quasi_incenter(v, g);
					new_points.push_back(g) ;
					flags.push_back(false);
				}
            }
            delaunay_->clear() ;
            delaunay_->begin_insert() ;
            for(unsigned int j=0; j<new_points.size(); j++) {
				Delaunay::Vertex_handle v = delaunay_->insert(new_points[j]) ;
            }
            convert_to_9_sheeted_covering();
            delaunay_->end_insert(false) ;

            double F = lloyd_energy() ;
			//if(period()) {
			//	delaunay_->clear_copies(false) ;
			//}
            
			if(redraw) {
                //std::cerr << "Lloyd energy = " << F << std::endl ;
				glut_viewer_redraw() ;
				std::cout << "Lloyd energy = " << 16*F/12 << std::endl ;
            }
            FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
                vec2 p0 = to_geex(v->point()) ;
                //std::cerr << "p0 = (" <<p0.x << "," << p0.y<< ")" << std::endl ;
                Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
                do {
                    const vec2& p1 = it->dual ;
                    //std::cerr << "p1 = (" <<p1.x << "," << p1.y<< ")" << std::endl ;

                    it++ ;
                    
                } while(it != delaunay_->incident_faces(v)) ;
            }
            
        }
    }

	void DelaunayCVT::lloyd_fpo(int nb_iter, bool redraw) {

	}

    //----------------------------------------------------------------------------------------------------

    void DelaunayCVT::newton_lloyd(int nb_iter, bool redraw) {
        use_theta_ = (mode() == THETA) ;

        if(use_theta_) {
            //FOR_EACH_VERTEX_DT(Delaunay, delaunay_, it) {
            //    if(false && it->dual_intersects_boundary) {
            //        it->theta = default_theta(it) ;
            //    }
            //    it->rho = default_rho(it) ;
            //}
        } 

        lbfgs_redraw = redraw ;

        int iLp = int(Lp_) ;
        iLp = gx_min(iLp, int(MAX_P)) ;
        iLp = gx_max(iLp, 0) ;

        if(use_theta_) {
            symbolic_ = true ;
            funcs_ = funcs_table_theta_[iLp] ;
        } else {
            if(iLp == 0) { 
                symbolic_ = false ; 
            } else {
                symbolic_ = true ;
                iLp-- ;
                funcs_ = funcs_table_[iLp] ;
            }
        }

        int n = use_theta_ ? delaunay_->nb_vertices() * 4 : delaunay_->nb_vertices() * 2 ;
        int m = 7 ;

        double* x = new double[n];
        
		
		if(use_theta_){
			//for(unsigned int i=0; i<delaunay_->all_vertices_.size(); i++) {
		 //       Delaunay::Vertex_handle it = delaunay_->all_vertices_[i] ;
			//	// periodic CVT
			//	if(period() && !delaunay_->is_primary(it))
			//		continue ;
			//	x[4*i  ] = it->point().x() ;
			//	x[4*i+1] = it->point().y() ;
			//	x[4*i+2] = it->theta ;
   //             x[4*i+3] = it->rho ;
			//}
		} else {
			for(unsigned int i=0; i<delaunay_->all_vertices_.size(); i++) {
		        Delaunay::Vertex_handle it = delaunay_->all_vertices_[i] ;
				// periodic CVT
				if(period() && !delaunay_->is_primary(it))
					continue ;
				x[2*i  ] = it->point().x() ;
				x[2*i+1] = it->point().y() ;
			}
		}

        double epsg = 0, epsf=0, epsx=0;
        
		//Optimizer* opt = new LBFGSOptimizer();
		Optimizer* opt = new HLBFGSOptimizer();
		
		opt->set_epsg(epsg);
		opt->set_epsf(epsf);
		opt->set_epsx(epsx);
		
		opt->set_M(m);
		opt->set_N(n);
		opt->set_max_iter(nb_iter);

		opt->set_newiteration_callback(newiteration_cvt2d);
		opt->set_funcgrad_callback(funcgrad_cvt2d);
		
		opt->optimize(x) ;
       

        set_vertices(x) ;
		delete opt;
		delete [] x;
    }

	void DelaunayCVT::move_vertices(vec2 delta) {


	}

    void DelaunayCVT::set_vertices(const double* x) {


    }


    //----------------------------------------------------------------------------------------------------

    void DelaunayCVT::get_cell_centroid(Delaunay::Vertex_handle v, vec2& g, double& V) {

    }

	void DelaunayCVT::get_boundary_cell_centroid(Delaunay::Vertex_handle v, vec2& g, double& V) {

	}

	//----------------------------------------------------------------------------------------------------
	void DelaunayCVT::get_primary_position(vec2& g) {
        if(g[0]<0) {
				while(g[0]<0)
					g[0]+=1 ;
			}
			if(g[0]>=1)
				while(g[0]>=1)
					g[0]-=1 ;
			if(g[1]<0) {
				while(g[1]<0)
					g[1]+=1 ;
			}
			if(g[1]>=1)
				while(g[1]>=1)
					g[1]-=1 ;
	}

	void DelaunayCVT::get_cell_primary_centroid(Delaunay::Vertex_handle v, vec2& g, double& V) {
        vec2 p0 = to_geex(v->point()) ;
        //std::cerr << "p0 = (" <<p0.x << "," << p0.y<< ")" << std::endl ;
		g.x = 0.0 ; g.y = 0.0 ; V = 0.0 ;
		Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
		Delaunay::Face_circulator jt = it ; jt++ ;
		do {
			const vec2& p1 = it->dual ;
			const vec2& p2 = jt->dual ;
			double Vi = triangle_area(p0, p1, p2) ;
			vec2 Gi = triangle_centroid(p0, p1, p2) ;
			V += Vi ;
			g += Vi * Gi ;
			it++ ; jt++ ;
            
		} while(it != delaunay_->incident_faces(v)) ;
        Delaunay::Vertex_circulator vt = delaunay_->adjacent_vertices(v);
		double s = (1.0 / V) ;
		g.x *= s ;
		g.y *= s ;	
		get_primary_position(g) ;
	}

    //----------------------------------------------------------------------------------------------------
    void DelaunayCVT::get_cell_quasi_incenter(Delaunay::Vertex_handle v, vec2& g) {

    }
   //----------------------------------------------------------------------------------------------------
    
    bool DelaunayCVT::get_fg(Delaunay::Vertex_handle v, double& f, vec2& grad_f) {

    }

    bool DelaunayCVT::get_fgv(Delaunay::Vertex_handle v, double& f, vec2& grad_f, double& area) {

    }

    void DelaunayCVT::set_anisotropy(const vec2& X, const vec2& Y) {

    }

   
    void DelaunayCVT::query_anisotropy(const vec2& P, vec2& U, vec2& V) {

    }
   
    void DelaunayCVT::set_anisotropy(const vec2& P) {

    }

    void DelaunayCVT::set_anisotropy(Delaunay::Vertex_handle V) {

    }

    void DelaunayCVT::get_PQR(PolygonVertex& V, int center_index) {

    }

    double DelaunayCVT::lloyd_energy() {
        double result = 0.0 ;

        return result ;
    }
    
}


//------------------------- LBFGS interface ------------------------------------------------

namespace Geex {

    void DelaunayCVT::funcgrad(const double* x, double& f, double* g, bool& valid) {

    }
    
    void DelaunayCVT::funcgrad_simple(const double* x, double& f, double* g, bool& valid) {

    }

    double DelaunayCVT::default_theta(const vec2& P) {
        double result = 0.0 ;

        return result ;
    }

    double DelaunayCVT::default_theta(Delaunay::Vertex_handle v) {
        return 0;
    }

    double DelaunayCVT::default_rho(Delaunay::Vertex_handle v) {

        return 1.0 ; 

    }

    void DelaunayCVT::funcgrad_symbolic(const double* x, double& f, double* g, bool& valid) {
        
    }
    
    
    void DelaunayCVT::add_to_fg(double& f, double* g) {
        
    }

}

//------------ Optimizer interface -------------

void funcgrad_cvt2d(int n, double* x, double& f, double* g) {

}

void newiteration_cvt2d(int n, const double* x, double f, const double* g, double gnorm) {

}

