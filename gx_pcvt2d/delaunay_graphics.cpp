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

#include "delaunay_graphics.h"
#include <delaunay_cvt.h>
#include <delaunay_io.h>
#include <Geex/combinatorics/map.h>
#include <glut_viewer/glut_viewer.h>

namespace Geex {


    static const double c1 = 0.35 ;
    static const double c2 = 0.5 ;
    static const double c3 = 1.0 ;

    static double color_table[12][3] = 
    {
        {c3, c2, c2},
        {c2, c3, c2},
        {c2, c2, c3},
        {c2, c3, c3},
        {c3, c2, c3},
        {c3, c3, c2},

        {c1, c2, c2},
        {c2, c1, c2},
        {c2, c2, c1},
        {c2, c1, c1},
        {c1, c2, c1},
        {c1, c1, c2}

    } ;

    static int random_color_index_ = 0 ;


    static void gl_random_color() {
        glColor3f(
            color_table[random_color_index_][0], 
            color_table[random_color_index_][1], 
            color_table[random_color_index_][2]
        ) ;
        random_color_index_ = (random_color_index_ + 1) % 12 ;
    }
    
    static void gl_random_color(int index) {
        random_color_index_ = index % 12 ;
        gl_random_color() ;
    }

    static void gl_randomize_colors() {
        random_color_index_ = 0 ;
    }

	static void gl_vertex_color(int degree) {
		if(degree==6) {
			glColor3f(0, 1, 0) ;
		}
		else if(degree==7) {
			glColor3f(1.0, 0.5, 0.25) ;
		}
		else if(degree==5) {
			glColor3f(0.42, 0.55, 0.9) ;
		}
		else if(degree>7) {
			glColor3f(0.5, 0, 0) ;
		}
		else {
			glColor3f(0, 0, 0.5) ;
		}
	}

    
    DelaunayGraphics::DelaunayGraphics(DelaunayCVT* CVT, DelaunayIO* IO) 
		: delaunay_(CVT->delaunay())
		, CVT_(CVT)
		, IO_(IO) {
        show_domain_ = GL_TRUE ; 
        vertices_size_ = 0.01 ;
		vertices_color_ = GL_TRUE ;
        centers_size_ = 0.0 ;
        show_primal_mesh_ = GL_FALSE ;
        show_dual_mesh_ = GL_TRUE ;
        colorize_ = GL_FALSE ;
        show_cells_ = GL_FALSE ;
        show_energy_ = GL_FALSE ;
        show_field_ = GL_FALSE ;
        quad_ratio_ = 0.0 ;
		show_boundary_cells_ = GL_FALSE ;
		show_pvd_euclidean_ = GL_FALSE ;
		show_copies_ = GL_FALSE ;
		show_disk_ = GL_FALSE ;
		show_mesh_ = GL_FALSE ;
		show_min_max_ = GL_FALSE ;
		show_inner_voronoi_ = GL_FALSE ;
		show_regularity_ = GL_FALSE ;
    }
    
    void DelaunayGraphics::draw() {
        non_convex_ = delaunay_->non_convex_mode_ ;
		glDisable(GL_LIGHTING) ;
        glUseProgramObjectARB(0) ;

        gl_randomize_colors() ;

        if(show_domain_) {
            draw_domain() ;
        }

        if(vertices_size_ != 0.0) {
            draw_vertices() ;
        }

        if(centers_size_ != 0.0) {
            draw_centers() ;
        }

        if(show_primal_mesh_) {
            draw_primal_mesh() ;
        }

		if(show_inner_voronoi_) {
			draw_inner_voronoi() ;
		}

        if(show_dual_mesh_) {
//            draw_dual_mesh() ;
			draw_dual_edges() ;
        }

        if(show_cells_) {
			if(CVT_->period())
				draw_non_hex_periodic_cells() ;
			else
				draw_non_hex_cells() ;
        }

        if(show_energy_) {
            draw_energy() ;
        }

        if(show_field_) {
            draw_field() ;
        }
		
		if(show_boundary_cells_) {
			draw_boundary_cells() ;
		}

		if(show_disk_) {
			draw_disk() ;
		}

		if(show_min_max_) {
			draw_min_max() ;
		}

		if(show_mesh_) {
			draw_mesh() ;
		}

		if(IO_->show_edge_hist()) {
			draw_edge_hist() ;
		}

		if(show_regularity_) {
            draw_regularity() ;
        }

		draw_selection() ;
    }

    void DelaunayGraphics::draw_domain() {
        glLineWidth(3) ;
        glColor3f(0.1, 0.1, 0.1) ;
        
		double x_min, y_min;
        double x_max, y_max;
		double dx, dy;
		x_min = delaunay_->x_min_;
		x_max = delaunay_->x_max_;
		y_min = delaunay_->y_min_;
		y_max = delaunay_->y_max_;
		dx = x_max - x_min;
		dy = y_max - y_min;
		
		for(unsigned int i=0; i<4; i++) {
			glBegin(GL_LINES) ;
			glVertex(vec2(x_min,y_min+i*dy)) ;
			glVertex(vec2(x_min+3*dx,y_min+i*dy)) ;
			glEnd();
		}
		for(unsigned int i=0; i<4; i++) {
			glBegin(GL_LINES) ;
			glVertex(vec2(x_min+i*dx,y_min)) ;
			glVertex(vec2(x_min+i*dx,y_min+3*dy)) ;
			glEnd();
		}
		/*glBegin(GL_LINES) ;
        for(unsigned int i=0; i<delaunay_->boundary_.size(); i++) {
            glVertex(delaunay_->boundary_[i].vertex[0]) ;
            glVertex(delaunay_->boundary_[i].vertex[1]) ;
        }
        glEnd() ;*/
    }

    void DelaunayGraphics::draw_vertices() {
        glDisable(GL_LIGHTING) ;
        glPointSize(GLfloat(vertices_size_ * 20)) ;
        int w,h ;
        glut_viewer_get_screen_size(&w, &h) ;
        notify_screen_resize(w,h) ;
        glEnable(GL_POINT_SMOOTH) ;
        glBegin(GL_POINTS) ;
//        begin_spheres() ;
		if(!vertices_color_) {
			FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
				if(!show_copies_ && !delaunay_->is_primary(v))
					continue ;
/*				if(v->locked) {
					glColor3f(1.0,0.0,0.0) ;
				} else*/ {
					glColor3f(0.0,0.0,0.0) ;
				}
				glVertex2f(v->point().x(), v->point().y()) ;
			}
		}
		else {
			FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
				if(!show_copies_ && !delaunay_->is_primary(v))
					continue ;
				int deg = delaunay_->degree(v) ;
				gl_vertex_color(deg) ;
				Point p = delaunay_->point(delaunay_->periodic_point(v));
				//std::cerr << pp.x() << std::endl ;
				glVertex2f(p.x(), p.y()) ;
				//std::cerr << "p = (" <<p.x() << "," << p.y()<< ")" << std::endl ;
			}
		}
//        end_spheres() ;
        glEnd() ;
        if(vertices_size_ > 0.6) {
            glBegin(GL_LINES) ;
            FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
                double R = 0.5 * radius(v) ;
                vec2 p = to_geex(v->point()) ;
                
				if(!show_copies_ && !delaunay_->is_primary(v))
					continue ;

                vec2 U,V ;
                CVT_->query_anisotropy(p,U,V) ;
/*
                vec2 U(cos(v->theta), sin(v->theta)) ;
                U *= R ;
                vec2 V(-U.y, U.x) ; 
*/              

                U *= 0.3 ;
                V *= 0.3 ;
                glVertex(p-U) ; glVertex(p+U) ;
                glVertex(p-V) ; glVertex(p+V) ;
            }
            glEnd() ;
        }
    }

    void DelaunayGraphics::draw_centers() {
        glDisable(GL_LIGHTING) ;
        glPointSize(int(centers_size_ * 20)) ;
        glColor3f(0.0,0.5,0.0) ;
        int w,h ;
        glut_viewer_get_screen_size(&w, &h) ;
        notify_screen_resize(w,h) ;
        glEnable(GL_POINT_SMOOTH) ;
        glBegin(GL_POINTS) ;
//        begin_spheres() ;
        FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
            double V ;
            vec2 g ;
			
			if(!show_copies_ && !delaunay_->is_primary(v))
				continue ;

			if(CVT_->period()) {
				if(!delaunay_->is_primary(v))
					continue ;
				CVT_->get_cell_primary_centroid(v, g, V) ;
			}
			else 
				CVT_->get_cell_centroid(v, g, V) ;
            vec2 p = g ;
            glVertex2f(p.x, p.y) ;
        }
        glEnd() ;
//        end_spheres() ;
    }

    void DelaunayGraphics::draw_cells(bool colorize, bool mesh) {
        FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
            if(colorize) {gl_random_color() ; }
			if(delaunay_->dual_cell_intersects_boundary(v) || delaunay_->dimension()==1) { 
                Polygon2* P = delaunay_->dual_convex_clip(v) ;
                if(mesh) {
                    glBegin(GL_LINES) ;
                    for(unsigned int i=0; i<P->size(); i++) {
                        glVertex((*P)[i].vertex[0]) ;
                        glVertex((*P)[i].vertex[1]) ;
                    }
                    glEnd() ;
                } else {
                    glBegin(GL_TRIANGLES) ;
                    for(unsigned int i=0; i<P->size(); i++) {
                        glVertex(to_geex(v->point())) ;
                        glVertex((*P)[i].vertex[0]) ;
                        glVertex((*P)[i].vertex[1]) ;
                    }
                    glEnd() ;
                }
            } else {
                glBegin(GL_POLYGON) ;
                Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
                do {
                    glVertex(it->dual) ;
                    it++ ;
                } while(it != delaunay_->incident_faces(v)) ;
                glEnd() ;
            } 
			//if(!delaunay_->dual_cell_infinite(v)) { 
			//	glBegin(GL_POLYGON) ;
			//	Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
			//	do {
			//		glVertex(it->dual) ;
			//		it++ ;
			//	} while(it != delaunay_->incident_faces(v)) ;
			//	glEnd() ;
			//}
        }
    }

	void DelaunayGraphics::draw_inner_voronoi() {
		glColor3f(0.0f, 1.0f, 0.0f) ;
		glLineWidth(2.0f) ;	
		glBegin(GL_LINES) ;
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(!delaunay_->is_primary(v)) continue ;
			std::vector<vec2> P ;
            delaunay_->compute_inner_voronoi(v, P) ;
			if(P.size()<4) continue ;
            for(unsigned int i=0; i<P.size(); i++) {
                glVertex(P[i]) ;
                glVertex(P[(i+1)%P.size()]) ;
            }
		}
		glEnd() ;

		glPointSize(vertices_size_*20.f) ;
		glColor3f(1.0, 0.5, 0.0) ;
		glBegin(GL_POINTS) ;
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(!delaunay_->is_primary(v)) continue ;
			std::vector<vec2> P ;
            delaunay_->compute_inner_voronoi(v, P) ;
			glVertex(polygon_centroid(P)) ;
		}
		glEnd() ;
	}

	void DelaunayGraphics::draw_dual_edges() {
		glColor3f(0.0f, 0.0f, 0.0f) ;
		glLineWidth(1.0f) ;

		FOR_EACH_EDGE_DT(Delaunay, delaunay_, e) {
			CGAL_Segment seg = delaunay_->baseclass::dual(e) ;
			glBegin(GL_LINES) ;
			/*if(CGAL::assign(ray, o)) {
				glVertex(to_geex(ray.source())) ;
				glVertex(to_geex(ray.source()+1e3*ray.to_vector())) ;
			}
			else if(CGAL::assign(seg, o)) {*/
				glVertex(to_geex(seg.start())) ;
				glVertex(to_geex(seg.end())) ;
			/*}
			else if(CGAL::assign(line, o)) {
				if(line.b()!=0) {
					glVertex(vec2(-1e3, -1e3*(-line.a()/line.b())-line.c()/line.b())) ;
					glVertex(vec2(1e3,  1e3*(-line.a()/line.b())-line.c()/line.b())) ;
				} else {
					glVertex(vec2(-1e3*(-line.b()/line.a())-line.c()/line.a(), -1e3)) ;
					glVertex(vec2( 1e3*(-line.b()/line.a())-line.c()/line.a(),  1e3)) ;
				}
			}*/
			glEnd() ;
		}
	}

    void DelaunayGraphics::draw_non_hex_cells() {
        FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
            if(delaunay_->dual_cell_intersects_boundary(v)) { 
                Polygon2* P = delaunay_->dual_convex_clip(v) ;
                glColor3f(0.0, 0.5, 0.0) ;
                glBegin(GL_POLYGON) ;
                for(unsigned int i=0; i<P->size(); i++) {
//                    glVertex(to_geex(v->point())) ;
                    glVertex((*P)[i].vertex[0]) ;
                    glVertex((*P)[i].vertex[1]) ;
                }
                glEnd() ;
            } else  {
				int degree = delaunay_->dual_facet_degree(v, CVT_->period()) ;
                if(delaunay_->dual_facet_degree(v, CVT_->period()) != 6 && !delaunay_->dual_cell_infinite(v)) {
                    if(degree > 7)
						glColor3f(0.0, 0.0, 0.5) ;
					else if(degree == 7)
						glColor3f(0.0, 0.0, 1) ;
					else if(degree == 5)
						glColor3f(0.0, 0.64, 0.9) ;
					else
						glColor3f(0.0, 0.9, 0.9) ;
                    glBegin(GL_POLYGON) ;
                    Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
                    do {
                        glVertex(it->dual) ;
                        it++ ;
                    } while(it != delaunay_->incident_faces(v)) ;
                    glEnd() ;
                }
            } 
        }
    }

	void DelaunayGraphics::draw_non_hex_periodic_cells() {
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			int degree = delaunay_->dual_facet_degree_period(v) ;
			if( degree != 6) {
				if(degree > 7)
					glColor3f(0.0, 0.0, 0.5) ;
				else if(degree == 7)
					glColor3f(0.0, 0.0, 1) ;
				else if(degree == 5)
					glColor3f(0.0, 0.64, 0.9) ;
				else
					glColor3f(0.0, 0.9, 0.9) ;

				if(delaunay_->dual_cell_intersects_boundary(v) || delaunay_->dimension()==1) { 
					Polygon2* P = delaunay_->dual_convex_clip(v) ;
					glBegin(GL_POLYGON) ;
					for(unsigned int i=0; i<P->size(); i++) {
						glVertex((*P)[i].vertex[0]) ;
						glVertex((*P)[i].vertex[1]) ;
					}
					glEnd() ;
				}
				else if(delaunay_->is_primary(v)) {
					glBegin(GL_POLYGON) ;
					Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
					if(delaunay_->dimension()>1) {
						do {
							glVertex(it->dual) ;
							it++ ;
						} while(it != delaunay_->incident_faces(v)) ;
					}
					glEnd() ;
				}
			}
		}
	}

    inline bool is_quad(Delaunay* del, Delaunay::Face_handle f, int i, double ratio) {
        Delaunay::Face_handle g = f->neighbor(i) ;
        if(del->is_infinite(g)) { return false ; }

        vec2 p1 = to_geex(f->vertex(0)->point()) ;
        vec2 q1 = to_geex(g->vertex(0)->point()) ;

        vec2 c1 = f->dual ;
        vec2 c2 = g->dual ;

        double R1 = (p1 - c1).length() ;
        double R2 = (q1 - c2).length() ;
        double R = 0.5 * (R1 + R2) ;

        double dC = (c2-c1).length() ;
        return (dC < ratio * R) ;
    }

	bool is_obtuse(Delaunay::Face_handle fh) {
	    vec2 v1 = to_geex(fh->vertex(0)->point()) ;
        vec2 v2 = to_geex(fh->vertex(1)->point()) ;
		vec2 v3 = to_geex(fh->vertex(2)->point()) ;
		double d1 = distance2(v2, v3) ;
		double d2 = distance2(v3, v1) ;
		double d3 = distance2(v1, v2) ;
		return d1+d2<d3 || d2+d3<d1 || d3+d1<d2 ;
	}

    void DelaunayGraphics::draw_primal_mesh() {
        glLineWidth(2) ;
        glColor3f(0.0, 0.0, 0.9) ;
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
        glBegin(GL_LINES) ;
        FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
            if(delaunay_->is_infinite(it)) { continue ; }
		//	if(it->dual_outside) continue ;

			/* we don't draw quad here */
            for(unsigned int i=0; i<3; i++) {
                unsigned int j1 = i + 1 ; 
                if(j1 == 3) { j1 = 0 ; }
                unsigned int j2 = (j1 + 1) ;
                if(j2 == 3) { j2 = 0 ; }

                //if(!is_quad(delaunay_, it, i, quad_ratio_)) {
                    Delaunay::Vertex_handle v1 = it->vertex(j1) ;
                    Delaunay::Vertex_handle v2 = it->vertex(j2) ;
                    glVertex2f(v1->point().x(), v1->point().y()) ;
                    glVertex2f(v2->point().x(), v2->point().y()) ;
                //}
            }

        }
        glEnd() ;
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
		glColor3f(0.8, 0.8, 0.8) ;
		glBegin(GL_TRIANGLES) ;
		FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) { continue ; }
			if(is_obtuse(it)) {
                for(int i=0; i<3; ++i) {
                    Delaunay::Vertex_handle vi = it->vertex(i) ;
					glVertex(to_geex(vi->point())) ;
				}
			}
		}
		glEnd() ;
    }

    void DelaunayGraphics::draw_dual_mesh() {

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
        //glLineWidth(3) ;
        //glColor3f(0.0, 0.0, 0.0) ;
        //draw_cells(false, true) ;
        //glLineWidth(9) ;
        //glColor3f(1.0, 1.0, 1.0) ;
        glColor3f(0.0, 0.0, 0.0) ;
		glLineWidth(1) ;
        draw_cells(false, true) ;
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
    }

    void DelaunayGraphics::is_min_max(Delaunay::Vertex_handle v, bool& is_min, bool& is_max) {
        is_min = true ;
        is_max = true ;
        Delaunay::Vertex_circulator it = delaunay_->adjacent_vertices(v) ;
        do {
            if(!delaunay_->is_infinite(it)) {
//                is_min = is_min && v->energy < it->energy ;
//                is_max = is_max && v->energy > it->energy ;
            }
            it++ ;
        } while(it != delaunay_->adjacent_vertices(v)) ;
    }

    void DelaunayGraphics::draw_energy() {

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

        double emin = 1e30 ;
        double emax = -1e30 ;

        Delaunay::Vertex_handle vmin = 0 ;
        Delaunay::Vertex_handle vmax = 0 ;

        FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
//            if(v->energy < emin) { emin = v->energy; vmin = v ; }
//            if(v->energy > emax) { emax = v->energy; vmax = v ; }
        }

//        std::cerr << "emin = " << emin << " emax =" << emax << std::endl ;

        double emean = 0.5 * (emin + emax) ;


        /* if(emax - emean > emean - emin) */ {
            glColor3f(1.0, 0.5, 0.5) ;
            if(vmax != 0) { draw_dual_facet(vmax) ; }
        } 
        /* else */ {
            glColor3f(0.5, 0.5, 1.0) ;
            if(vmin != 0) { draw_dual_facet(vmin) ; }
        }


        double scale = (emax - emin) ;
        if(::fabs(scale) < 1e-30) { scale = 1.0 ; }
        scale = 1.0 / scale ;

        const vec3 c1(0.0, 0.0, 1.0) ;
        const vec3 c2(1.0, 0.0, 0.0) ;

        FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
            //double e = scale * (v->energy - emin) ;
            //glColor(mix(c1, c2, e)) ;
            draw_dual_facet(v) ;
        }
    }

	void DelaunayGraphics::draw_regularity() {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

        double emin = 1e30 ;
        double emax = -1e30 ;

        Delaunay::Vertex_handle vmin = 0 ;
        Delaunay::Vertex_handle vmax = 0 ;

        FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
			if(!delaunay_->is_primary(v))
				continue ;
			if(v->regularity < emin) { emin = v->regularity; vmin = v ; }
            if(v->regularity > emax) { emax = v->regularity; vmax = v ; }
        }

//        std::cerr << "emin = " << emin << " emax =" << emax << std::endl ;

        double emean = 0.5 * (emin + emax) ;

		// 2d optimital measure
		emin = 5.0/(36*sqrt(3.0)) ;
//		emax = 5.0/(36*sqrt(3.0)) * 1.2;

        ///* if(emax - emean > emean - emin) */ {
        //    glColor3f(1.0, 0.5, 0.5) ;
        //    if(vmax != 0) { draw_dual_facet(vmax) ; }
        //} 
        ///* else */ {
        //    glColor3f(0.5, 0.5, 1.0) ;
        //    if(vmin != 0) { draw_dual_facet(vmin) ; }
        //}

        double scale = (emax - emin) ;
        if(::fabs(scale) < 1e-30) { scale = 1.0 ; }
        scale = 1.0 / scale ;

//		std::cout << "max/min ratio: " << emax/emin << std::endl ;

        const vec3 c1(0.0, 0.0, 1.0) ;
        const vec3 c2(1.0, 0.0, 0.0) ;

        FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
			double e ;
			if(delaunay_->is_primary(v)) {
				e = scale * (v->regularity - emin) ;
			}
			else {
				Delaunay::Vertex_handle vp = delaunay_->vertices()[v->index] ;
				e = scale * (vp->regularity - emin) ;
			}
            glColor(mix(c1, c2, e)) ;
//            draw_dual_facet(v) ;

			if(delaunay_->dual_cell_intersects_boundary(v) || delaunay_->dimension()==1) { 
				Polygon2* P = delaunay_->dual_convex_clip(v) ;
				glBegin(GL_POLYGON) ;
				for(unsigned int i=0; i<P->size(); i++) {
					glVertex((*P)[i].vertex[0]) ;
					glVertex((*P)[i].vertex[1]) ;
				}
				glEnd() ;
			}
			else if(delaunay_->is_primary(v)) {
				glBegin(GL_POLYGON) ;
				Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
				if(delaunay_->dimension()>1) {
					do {
						glVertex(it->dual) ;
						it++ ;
					} while(it != delaunay_->incident_faces(v)) ;
				}
				glEnd() ;
			}
        }
	}
	
	void DelaunayGraphics::draw_field() {
        /*double x_min, y_min, z_min ;
        double x_max, y_max, z_max ;
        delaunay_->get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
        double dx = x_max - x_min ;
        double dy = y_max - y_min ;
        glPointSize(4) ;
        glBegin(GL_POINTS) ;

        for(int i=0; i<100; i++) {
            for(int j=0; j<100; j++) {
                vec2 p(x_min + dx * double(i)/100.0, y_min + dy * double(j)/100.0) ;
                gl_random_color(delaunay_->segments().locate(p)) ;
                glVertex(p) ;
            }
        }

        glPointSize(8) ;
        for(SegmentDelaunay::Vertex_iterator it = delaunay_->segments().finite_vertices_begin() ;
            it != delaunay_->segments().finite_vertices_end(); it++
        ) {
            gl_random_color(it->index) ;
            glVertex(to_geex(it->point())) ;
       }
    
        glEnd() ;*/
    }

    void DelaunayGraphics::draw_dual_facet(Delaunay::Vertex_handle v) {
        if(delaunay_->dual_cell_intersects_boundary(v)) { 
            Polygon2* P = delaunay_->dual_convex_clip(v) ;
            glBegin(GL_TRIANGLES) ;
            for(unsigned int i=0; i<P->size(); i++) {
                glVertex(to_geex(v->point())) ;
                glVertex((*P)[i].vertex[0]) ;
                glVertex((*P)[i].vertex[1]) ;
            }
            glEnd() ;
        } else  {
            glBegin(GL_POLYGON) ;
            Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
            do {
                glVertex(it->dual) ;
                it++ ;
            } while(it != delaunay_->incident_faces(v)) ;
            glEnd() ;
        } 
    }


    double DelaunayGraphics::radius(Delaunay::Vertex_handle v) {
        double result = 1e30 ;
        vec2 p1 = to_geex(v->point()) ;
        Delaunay::Vertex_circulator it = delaunay_->adjacent_vertices(v) ;
        do {
            if(!delaunay_->is_infinite(it)) {
                vec2 p2 = to_geex(it->point()) ;
                result = gx_min(result, (p2 - p1).length2()) ;
            }
            it++ ;
        } while(it != delaunay_->adjacent_vertices(v)) ;
        return sqrt(result) ;
    }

	void DelaunayGraphics::draw_boundary_cells() {
		int current_color = 0 ;
		//		delaunay_->compute_rvd() ;

		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(!delaunay_->is_primary(v)) continue ;

			if(v->dual_intersects_boundary || delaunay_->dimension()==1) {
				if(colorize_) {
					current_color = random_color_index_ ;
					gl_random_color() ;
				}
				else glColor3f(1.f, 1.f, 0.0f) ;

				Polygon2 *pl = delaunay_->dual_convex_clip(v) ; //rvd_[v->index] ;
				glBegin(GL_POLYGON) ; // convex only, need to improve
				for(unsigned int j=0; j<pl->size(); ++j) {
					glVertex((*pl)[j].vertex[0]) ;
					glVertex((*pl)[j].vertex[1]) ;
				}
				glEnd() ;

				// draw border
				glColor3f(0.f, 0.f, 0.f) ;
				glBegin(GL_LINES) ;
				for(unsigned int j=0; j<pl->size(); ++j) {
					glVertex((*pl)[j].vertex[0]) ;
					glVertex((*pl)[j].vertex[1]) ;
				}
				glEnd() ;

				// draw mirrors
				std::set<Delaunay::Vertex_handle>& mirrors = delaunay_->mirrors(v->index) ;
				for(std::set<Delaunay::Vertex_handle>::iterator it=mirrors.begin(); it!=mirrors.end(); ++it) {
					if((*it)->dual_intersects_boundary) {
						Polygon2 *pl = delaunay_->dual_convex_clip(*it) ; //rvd_[v->index] ;
						if(colorize_) gl_random_color(current_color) ;

						if(show_pvd_euclidean_) { // offset polygon
							//int mid = delaunay_->domain_idx(to_geex((*it)->point())) ;
							vec2 offset = to_geex(v->point()) - to_geex((*it)->point());
							for(unsigned int i=0; i<pl->size(); ++i) {
								(*pl)[i].vertex[0] = (*pl)[i].vertex[0] + offset ;
								(*pl)[i].vertex[1] = (*pl)[i].vertex[1] + offset ;
							}
						}

						glBegin(GL_POLYGON) ; // convex only, need to improve
						for(unsigned int j=0; j<pl->size(); ++j) {
							glVertex((*pl)[j].vertex[0]) ;
							glVertex((*pl)[j].vertex[1]) ;
						}
						glEnd() ;

						// draw border
						glColor3f(0.f, 0.f, 0.f) ;
						glBegin(GL_LINES) ;
						for(unsigned int j=0; j<pl->size(); ++j) {
							glVertex((*pl)[j].vertex[0]) ;
							glVertex((*pl)[j].vertex[1]) ;
						}
						glEnd() ;
					}
				}
			}
			else {
				glColor3f(0.0, 0.0, 0.0) ;
				glBegin(GL_LINE_LOOP) ;
				Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
				if(delaunay_->dimension()>1) { // to avoid degenerate case: all the points are on a line.
					do {
						glVertex(it->dual) ;
						it++ ;
					} while(it != delaunay_->incident_faces(v)) ;
				}
				glEnd() ;
			}
		}
	}

	void DelaunayGraphics::draw_edge_hist() {
		std::vector<double>& hist = IO_->edge_histogram() ;
		int w, h ;
		glut_viewer_get_screen_size(&w, &h) ;
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0, w, 0, h);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		double maxhist = 0 ;
		for(unsigned int i=0; i<hist.size(); ++i) {
			if(maxhist < hist[i])
				maxhist = hist[i] ;
		}

		glColor4f(0.76, 0.87, 0.85, 0.5) ;
		glBegin(GL_QUADS) ;
		glVertex2i(w-256, h-256) ;
		glVertex2i(w, h-256) ;
		glVertex2i(w, h) ;
		glVertex2i(w-256, h) ;
		glEnd() ;
		
		glLineWidth(3.0) ;
		glColor3f(0.2, 0.2, 0.2) ;
		glBegin(GL_LINES) ;
		for(int i=0; i<256; ++i) {
			glVertex2i(w-256+i,   h-256) ;
			glVertex2i(w-256+i,   h-256+hist[i]/maxhist*256) ;
		}
		glEnd() ;

		glMatrixMode(GL_PROJECTION);
		glPopMatrix() ;

		glMatrixMode(GL_MODELVIEW);
		glPopMatrix() ;
	}

	static void drawCircle(float cx, float cy, float r, int num_segments) 
	{ 
		glBegin(GL_LINE_LOOP); 
		for(int ii = 0; ii < num_segments; ii++) 
		{ 
			float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);//get the current angle 
			float x = r * cosf(theta);//calculate the x component 
			float y = r * sinf(theta);//calculate the y component 
			glVertex2f(x + cx, y + cy);//output vertex 
		} 
		glEnd(); 
	}

	static inline void draw_vertex_disk(GLUquadric *q, const Delaunay::Vertex_handle& v, double radius) {
		glPushMatrix() ;
		glTranslatef(v->point().x(), v->point().y(), 0) ;
		glColor3f(.7, 0.7, 0.7) ;
		gluDisk(q, 0, radius, 100, 100) ;
		glPopMatrix() ;
		glColor3f(0.0, 0.0, 0.0) ;
		drawCircle(v->point().x(), v->point().y(), radius, 100) ;
	}

	void DelaunayGraphics::draw_disk() {
		GLUquadric * q = gluNewQuadric( ) ;
		glLineWidth(2.0) ;

		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(!delaunay_->is_primary(v))
				continue ;
			draw_vertex_disk(q, v, delaunay_->sample_radius()) ;
		}
		gluDeleteQuadric(q) ;

		GLint func ;
		glGetIntegerv(GL_DEPTH_FUNC, &func) ;
		glDepthFunc(GL_LEQUAL) ;

		// draw skeleton of the gaps
		std::vector<std::vector<vec2> >& skeleton = delaunay_->void_skeleton() ;
		glLineWidth(4.0) ;
		glColor3f(0.0, 1.0, 0.0) ;
		for(int i=0; i<skeleton.size(); ++i) {
				glBegin(GL_LINES) ;
				for(unsigned int j=0; j<skeleton[i].size()-1; j+=2)  {
					glVertex(skeleton[i][j]) ;
					glVertex(skeleton[i][j+1]) ;
				}
				glEnd() ;
		}

		// draw gaps
		std::vector<std::vector<vec2> >& regions = delaunay_->void_regions() ;
		for(int i=0; i<regions.size(); ++i) {
				gl_random_color() ;
				glBegin(GL_POLYGON) ;
				for(unsigned int j=0; j<regions[i].size(); ++j)  {
					glVertex(regions[i][j]) ;
				}
				glEnd() ;
		}
		glDepthFunc(func) ;
	}

	void DelaunayGraphics::draw_min_max() {
		GLUquadric * q = gluNewQuadric( ) ;
		glLineWidth(2.0) ;

		// draw triangle with angle < 30
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(!show_copies_ && !(delaunay_->is_primary(v) || delaunay_->neighbor_to_primary(v)))
				continue ;

			vec2 p = to_geex(v->point()) ;
			Delaunay::Vertex_circulator cir = delaunay_->adjacent_vertices(v) ;
			Delaunay::Vertex_circulator cir2 = cir++ ;
			do {
				vec2 e1 = to_geex(cir->point())-p ;
				vec2 e2 = to_geex(cir2->point())-p ;
				double angle = to_degree(acos(dot(e1, e2)/(e1.length()*e2.length()))) ;
				if(angle < 30) {
					glColor3f(1.0, 0, 0) ;
					glBegin(GL_POLYGON) ;
					glVertex(p) ;
					glVertex(to_geex(cir->point())) ;
					glVertex(to_geex(cir2->point())) ;
					glEnd() ;

					draw_vertex_disk(q, v, delaunay_->sample_radius()) ;
					draw_vertex_disk(q, cir, delaunay_->sample_radius()) ;
					draw_vertex_disk(q, cir2, delaunay_->sample_radius()) ;
				}
				++cir  ;
				++cir2 ;
			} while(cir!=delaunay_->adjacent_vertices(v)) ;
		}
		gluDeleteQuadric(q) ;

		// draw faces with max/min area
		Delaunay::Face_handle fmin, fmax ;
		double area_min=1e10, area_max=-1e10 ;

		FOR_EACH_FINITE_FACE_DT(Delaunay, delaunay_, f) {
			if(!f->dual_outside) {
				double area = triangle_area(to_geex(f->vertex(0)->point()), 
											to_geex(f->vertex(1)->point()), 
											to_geex(f->vertex(2)->point())) ;
				if(area < area_min) {
					area_min = area ;
					fmin = f ;
				}
				if(area > area_max) {
					area_max = area ;
					fmax = f ;
				}
			}
		}

		glColor3f(1.0, 0, 0) ;
		glBegin(GL_POLYGON) ;
		glVertex(to_geex(fmin->vertex(0)->point())) ;
		glVertex(to_geex(fmin->vertex(1)->point())) ;
		glVertex(to_geex(fmin->vertex(2)->point())) ;
		glEnd() ;

		glBegin(GL_POLYGON) ;
		glVertex(to_geex(fmax->vertex(0)->point())) ;
		glVertex(to_geex(fmax->vertex(1)->point())) ;
		glVertex(to_geex(fmax->vertex(2)->point())) ;
		glEnd() ;
	} 

	void DelaunayGraphics::draw_mesh() {
		Map* M = delaunay_->map() ;

		GLint func ;
		glGetIntegerv(GL_DEPTH_FUNC, &func) ;
		glDepthFunc(GL_LEQUAL) ;

		glColor3f(0.2, 0.2, 0.2) ;
		glBegin(GL_LINES) ;
		for(Map::Facet_iterator it=M->facets_begin(); it!=M->facets_end(); ++it) {
			Map::Halfedge *h = it->halfedge() ;
			do {
				glVertex(h->prev()->vertex()->point()) ;
				glVertex(h->vertex()->point()) ;
				h = h->next() ;
			} while(h!=it->halfedge()) ;

		}
		glEnd() ;

		if(vertices_size_>0) {
			glBegin(GL_POINTS) ;
			for(Map::Vertex_iterator it=M->vertices_begin(); it!=M->vertices_end(); ++it) {
				gl_vertex_color(it->degree()) ;
				glVertex(it->point()) ;
			}
			glEnd() ;
		}

		//// draw picked vertices and edges 
		//IrregEditor *editor = delaunay_->editor() ;
		//std::vector<Map::Vertex*>&   vsel = editor->picked_verts() ;
		//std::vector<Map::Halfedge*>& esel = editor->picked_edges() ; 
		//GLUquadric *q = gluNewQuadric() ;

		//glEnable(GL_LIGHTING) ;
		//glColor3f(0.1, 0.8, 0.8) ;
		//for(unsigned int i=0; i<esel.size(); ++i) {
		//	Map::Halfedge * h = esel[i] ;
		//	vec3 v0 = h->prev()->vertex()->point() ;
		//	vec3 v1 = h->vertex()->point() ;
		//	vec3 dir = normalize(v1-v0) ;
		//	
		//	glPushMatrix() ;
		//	glTranslatef(v0.x, v0.y, 0) ;
		//	glRotatef(90, -dir.y, dir.x, 0.0) ;
		//	gluCylinder(q, 0.01, 0.01, distance(v0, v1), 20, 20) ;
		//	glPopMatrix() ;
		//}

		//glDisable(GL_LIGHTING) ;
		//for(unsigned int i=0; i<vsel.size(); ++i) {
		//	vec3 c = vsel[i]->point() ;
		//	int  deg = vsel[i]->degree() ;
		//	
		//	gl_vertex_color(deg) ;
		//	glPushMatrix() ;
		//	glTranslatef(c.x, c.y, 0) ;
		//	gluSphere(q, vertices_size_, 20, 20) ;
		//	glPopMatrix() ;
		//}

		//gluDeleteQuadric(q) ;

		//// draw current facet where the picked point locates
		//Map::Facet* f = editor->cur_facet() ;
		//if(f!=nil) {
		//	glColor3f(1.0, 0, 0) ;
		//	glBegin(GL_POLYGON) ;
		//	Map::Halfedge* h = f->halfedge() ;
		//	do{
		//		glVertex(h->vertex()->point()) ;
		//		h = h->next() ;
		//	} while(h!= editor->cur_facet()->halfedge()) ;
		//	glEnd() ;
		//}

		glDepthFunc(func) ;
	}

	void DelaunayGraphics::draw_selection() {
		std::vector<Delaunay::Vertex_handle>&  vsel = delaunay_->vertices_sel() ;
		std::vector<Delaunay::Edge>& esel = delaunay_->edges_sel() ; 
		GLUquadric *q = gluNewQuadric() ;

		glEnable(GL_LIGHTING) ;
		glColor3f(0.1, 0.8, 0.8) ;
		for(unsigned int i=0; i<esel.size(); ++i) {
			Delaunay::Face_handle f = esel[i].first ;
			int  idx = esel[i].second ;
			vec2 v0 = to_geex(f->vertex(f->cw(idx))->point()) ;
			vec2 v1 = to_geex(f->vertex(f->ccw(idx))->point()) ;
			vec2 dir = normalize(v1-v0) ;
			
			glPushMatrix() ;
			glTranslatef(v0.x, v0.y, 0) ;
			glRotatef(90, -dir.y, dir.x, 0.0) ;
			gluCylinder(q, vertices_size_, vertices_size_, distance(v0, v1), 20, 20) ;
			glPopMatrix() ;
		}

//		glDisable(GL_LIGHTING) ;
		for(unsigned int i=0; i<vsel.size(); ++i) {
			vec2 c = to_geex(vsel[i]->point()) ;
			int  deg = delaunay_->degree(vsel[i]) ;
			
			gl_vertex_color(deg) ;
			glPushMatrix() ;
			glTranslatef(c.x, c.y, 0) ;
			gluSphere(q, vertices_size_, 20, 20) ;
			glPopMatrix() ;
		}

		gluDeleteQuadric(q) ;
	}
}
