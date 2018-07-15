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

#ifndef __DELAUNAY__
#define __DELAUNAY__

//#define NOMINMAX 

#include "delaunay_cgal.h"
#include "polygons.h"
#include <GL/gl.h>

namespace Geex {

enum SamplingMode { PD_CCVT, PD_VORONOI, PD_BOUNDARY, PD_DARTTHROW, PD_PENROSE } ;
#define SamplingModeNames "Vornoi | CCVT" | "Pure" | "DartThrow" | "Penrose"

class DelaunayGraphics ;
class DelaunayCVT ;
class IrregEditor ;
class Map ;
class GridCell ;

//______________________________________________________________________________________

    class SegmentDelaunay : public DelaunayBase {
        typedef DelaunayBase baseclass ;
    public:
        void insert(int id, const vec2& p1, const vec2& p2, double step_length) ;
        int locate(const vec2& p) ;
    protected:
        void insert(int id, const vec2& p) ;
    } ;

    class Delaunay : public DelaunayBase {
        typedef DelaunayBase baseclass ;

    public:
        Delaunay() ;
		~Delaunay() ;

        void save(const std::string& filename) ;
        void load(const std::string& filename) ;
        void load_boundary(const std::string& filename) ;
        void set_non_convex_mode(bool x) { non_convex_mode_ = x ; }

        void get_bbox(
            real& x_min, real& y_min, real& z_min,
            real& x_max, real& y_max, real& z_max
        ) ;

        bool in_boundary(const vec2& p) { 
            return non_convex_mode_ ? boundary_.contains(p) : boundary_convex_.contains(p) ; 
        }

        GLboolean& insert_boundary() { return insert_boundary_ ; }
//		GLboolean& insert_full_copy() { return insert_full_copy_ ; }
//		GLboolean& period() { return period_ ; }

		void compute_inner_voronoi(Vertex_handle vh, std::vector<vec2>& poly) ;
        
		// ------------------------------------ Posson disk
		
		vec2 origin() { return 0.5*vec2(x_min_+y_min_, x_max_+y_max_) ; } 
		GLfloat& perturb() { return perturb_ ; }
		void do_perturb() ;
		void digree_hist(double& plt5, double& p5, double& p6, double& p7, double& pgt7) ;
		void edge_hist(double& p66, double& p56, double& p76, double& p55, double& p77, double& p57) ;
		std::vector<double>& edge_hist() { return edge_hist_ ; }
		
		
		void generate_poisson_disk() ;
		void sample_poisson_grid(double radius, int nmaxfail) ;
		void init_sample_grid(int res, double radius, int nmaxfail) ;
		void insert_grid_samples_period(int res) ;
		void insert_grid_samples_clipped(int res) ;
//		bool sample_grid_gaps(int res, double radius) ;
		inline bool is_hit_period(vec2& p, int u, int v, int res, double sqr_radius) ;

		void sample_poisson(double radius=0.001, bool isTiled=false, bool maximize=1, int minMaxThrows=1000, 
	        int multiplier=false, int relax=0, const char* method="Penrose") ;
		void set_vertices(std::vector<vec2>& points) ;
		double& sample_radius() { return sample_radius_ ; }
		SamplingMode& sampling_mode() { return sampling_mode_ ; }

		void fast_poisson_disk(double exclude_radius) ;
		Vertex_handle insert_vertex_periodic(vec2& p, double exclude_radius) ;
		void compute_voronoi_area(Delaunay::Vertex_handle& v, double exclude_radius) ;
		vec2 new_vertex_in_voronoi(Delaunay::Vertex_handle& v, double exclude_radius);
		void compute_gap_area(Delaunay::Vertex_handle& v, double exclude_radius) ;
		vec2 new_vertex_in_gap(Delaunay::Vertex_handle& v, double exclude_radius) ;
		bool on_convex_hull(Delaunay::Vertex_handle& v) ;
		void compute_triangle_voids() ;
		bool is_dual_gap(Delaunay::Face_handle f, double exclude_radius) ;
		void compute_void_region(Delaunay::Face_handle f, double exclude_radius, std::vector<vec2>& poly) ;
		void update_gaps(double exclude_radius) ;
		void compute_void_skeleton(Face_handle f, double exclude_radius) ;
		void compute_void_region(Face_handle f, double exclude_radius) ;
		std::vector<std::vector<vec2> >& void_skeleton() { return void_skeleton_ ; } 
		std::vector<std::vector<vec2> >& void_regions() { return void_regions_ ; } 
		void insert_point_poisson(vec2& p, double radius) ;
		void insert_copy_poisson(vec2& p, int index, double radius) ;
		void update_neighbor_cells(Vertex_handle& vh, double radius) ;

		bool pick_vertex(vec2& pt) ;
		bool pick_edge(vec2& pt) ;

		std::vector<Vertex_handle>&  vertices_sel() { return vertices_sel_ ; }// selected vertices 
		std::vector<Edge>& edges_sel() { return edges_sel_ ; }     // selected edges
		void clear_selection() {
			vertices_sel_.clear() ; 
			edges_sel_.clear() ; 
		}
		Map* map()             { return map_ ;    }
		IrregEditor* editor()  { return editor_ ; }
		bool& map_edit_mode()  { return map_edit_mode_ ; }
		// ------------------------------------ Delaunay 

        int nb_vertices() const { return all_vertices_.size() ; }
        void clear() ;
        void begin_insert() ;
        Vertex_handle insert(const vec2& p) ;
		Vertex_handle insert(const vec2& p, int index) ;
        Vertex_handle nearest(const vec2& p) ;
        void remove(const vec2& p) ;
        void end_insert(bool redraw = true) ;

        void insert_random_vertex() ;
        void insert_random_vertices(int nb) ;
		void insert_grid(int nb) ;

        std::vector<Vertex_handle>& vertices() { return all_vertices_ ; }

		// ------------------------------------ Periodic Delaunay
		void compute_pvd() ;
		void compute_pvd2() ;
		void clear_copies(bool redraw=false) ;
		void insert_copies(bool full_copy=false, bool redraw=false) ;
		void begin_insert_copies() ;
		void end_insert_copies(bool redraw=false) ;
		void insert_copies_full(bool redraw=false) ;
		void insert_copies_ring(bool redraw=false) ;
		void insert_copies_poisson(double radius, bool redraw=false) ;

		bool is_primary(Vertex_handle v) { return v->domain == 0 ; }
		bool neighbor_to_primary(Vertex_handle v) ;
		bool is_full_hex() ; // check pvd configuration
		int  nb_primary() ;
		vec2 copy_point(vec2& p, PolygonEdge& e) ; // 
		vec2 copy_point(vec2& p, PolygonVertex& e) ; // 
		int  domain_idx(vec2& p) ;
		vec2 mirror_vertex_point(vec2& p, int vidx) ; // 
		vec2 mirror_edge_point(vec2& p, int eidx) ; // 
		vec2 translate(int domain_idx, vec2 p) ; // translate periodic coordinate to Euclidean
		std::set<Vertex_handle>& mirrors(int pid) {
			return mirrors_[pid] ;
		}
		void compute_edge_length() ; // test whether the edges of a hex with same length
		std::map<int, std::set<Vertex_handle> >& mirrors() { return mirrors_ ; }

		//void save_vertices(const std::string& filename) ;
		//void load_vertices(const std::string& filename) ;

		// ----------------------------------- RVD 2D
//		typedef std::vector<PolygonEdge> PolyLine ;
		std::vector<Polygon2> rvd_ ;
		void compute_rvd(bool closed=true) ;
		bool is_boundary_cell(Delaunay::Vertex_handle v) { return rvd_[v->index].size()>0 ; }

        // ------------------------------------ Delaunay combinatorics ----------------

        vec2 dual(Face_handle c) { return c->dual ; }
        bool dual_cell_intersects_boundary(Vertex_handle v) {
            return v->dual_intersects_boundary ;
        }
		bool dual_cell_infinite(Vertex_handle v) {
			return v->dual_infinite ;
		}
 
        Line<real> get_dual_line(Edge e) {
			if(dimension()==1)
				return median_line(to_geex(e.first->vertex(1)->point()), to_geex(e.first->vertex(0)->point())) ;
			else
				return median_line(
					to_geex(e.first->vertex(ccw(e.second))->point()),
					to_geex(e.first->vertex( cw(e.second))->point())
				) ;
        }

        Polygon2* dual_convex_clip(Vertex_handle v, bool close = true) ;
		Polygon2* dual_convex_clip(Polygon2& from, Vertex_handle v, bool close = true) ;
        int dual_facet_degree(Vertex_handle v, bool period=false) ;
		int dual_facet_degree_period(Vertex_handle v) ;

        // --------------

        SegmentDelaunay& segments() { return segments_ ; }

    protected:
        std::vector<Vertex_handle> all_vertices_ ;
		std::map<int, std::vector<Vertex_handle> > mirror_vertices_ ;
        bool non_convex_mode_ ;

        Polygon2 boundary_ ;
        Convex boundary_convex_ ;
        bool cached_bbox_ ;
        double x_min_, y_min_, x_max_, y_max_ ;

        Polygon2 ping_ ;
        Polygon2 pong_ ;

        SegmentDelaunay segments_ ;

        friend class DelaunayGraphics ;
        friend class DelaunayCVT ;

        bool opened_ ;
        GLboolean insert_boundary_ ;
	//	GLboolean period_ ;
//		GLboolean insert_full_copy_ ;
		unsigned int nb_master_ ;
		const static vec2 v_offset_[4] ;
		const static vec2 e_offset_[4] ;
		const static vec2 domain_offset_[9] ;

		// used for poisson disk sampling
		GLfloat perturb_ ;
		std::vector<double> edge_hist_ ;
		int bin_size_ ;
		double sample_radius_ ; 
		GridCell **grid_ ;
		std::map<int, std::vector<Vertex_handle> > vertices_ ; 
		std::map<int, std::set<Vertex_handle> > mirrors_ ;
		SamplingMode sampling_mode_ ;
		std::vector<std::vector<vec2> > void_skeleton_ ;
		std::vector<std::vector<vec2> > void_regions_ ;

		// user interaction
		std::vector<Vertex_handle>  vertices_sel_ ; // selected vertices 
		std::vector<Delaunay::Edge> edges_sel_ ;     // selected edges

		bool        map_edit_mode_ ;
		IrregEditor *editor_ ;
		Map         *map_ ;
    } ;

//______________________________________________________________________________________

#define FOR_EACH_VERTEX_DT(TCLASS, INSTANCE, IT)                       \
    for(                                                            \
        TCLASS::Vertex_iterator IT = (INSTANCE)->vertices_begin() ; \
        IT != (INSTANCE)->vertices_end(); IT++                      \
    )

#define FOR_EACH_FACE_DT(TCLASS, INSTANCE, IT)                      \
    for(                                                         \
        TCLASS::Face_iterator IT = (INSTANCE)->faces_begin() ;   \
        IT != (INSTANCE)->faces_end(); IT++                      \
    )

#define FOR_EACH_FINITE_FACE_DT(TCLASS, INSTANCE, IT)                              \
    for(                                                                        \
        TCLASS::Finite_faces_iterator IT = (INSTANCE)->finite_faces_begin() ;   \
        IT != (INSTANCE)->finite_faces_end(); IT++                              \
    )

#define FOR_EACH_EDGE_DT(TCLASS, INSTANCE, IT)                      \
	for(                                                         \
	TCLASS::Edge_iterator IT = (INSTANCE)->edges_begin() ;   \
	IT != (INSTANCE)->edges_end(); IT++                      \
	)
    
}

#endif
