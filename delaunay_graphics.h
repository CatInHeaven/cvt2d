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


#ifndef __DELAUNAY_GRAPHICS__
#define __DELAUNAY_GRAPHICS__

#include "delaunay.h"
#include <Geex/graphics/opengl.h>

namespace Geex {

    class DelaunayCVT ;
	class DelaunayIO ;
    
    class DelaunayGraphics {
    public:
        DelaunayGraphics(DelaunayCVT* CVT, DelaunayIO* IO=nil) ;
        void draw() ;

        GLboolean& show_domain() { return show_domain_ ; }
        GLfloat& vertices_size() { return vertices_size_ ; }
		GLboolean& vertices_color() { return vertices_color_ ; }
        GLfloat& centers_size() { return centers_size_ ; }
        GLboolean& show_primal_mesh() { return show_primal_mesh_ ; }
        GLboolean& show_dual_mesh() { return show_dual_mesh_ ; }
        GLboolean& colorize() { return colorize_ ; }
        GLboolean& show_cells() { return show_cells_ ; }
        GLboolean& show_energy() { return show_energy_ ; }
        GLboolean& show_field() { return show_field_ ; }
        GLfloat& quad_ratio() { return quad_ratio_ ; }
		GLboolean& show_boundary_cells() { return show_boundary_cells_ ; }
		GLboolean& show_pvd_euclidean() { return show_pvd_euclidean_ ; }
		GLboolean& show_copies() { return show_copies_ ; }
		GLboolean& show_disk() { return show_disk_ ; }
		GLboolean& show_mesh() { return show_mesh_ ; }
		GLboolean& show_min_max() { return show_min_max_ ;  } 
		GLboolean& show_inner_voronoi() { return show_inner_voronoi_ ; } 
		GLboolean& show_regularity() { return show_regularity_ ; }

    protected:
        void draw_domain() ;
        void draw_vertices() ;
        void draw_centers() ;
        void draw_cells(bool colorize = true, bool mesh = false) ;
        void draw_edge() ;
        void draw_dual_point() ;
        void draw_primal_mesh() ;
        void draw_dual_mesh() ;
		void draw_dual_edges() ;
        void draw_non_hex_cells() ;
		void draw_non_hex_periodic_cells() ;
        void draw_energy() ;
        void draw_field() ;
		void draw_boundary_cells() ;
		void draw_edge_hist() ;
		void draw_disk() ;
		void draw_mesh() ;
		void draw_selection() ;
		void draw_min_max() ;
		void draw_inner_voronoi() ;
        void is_min_max(Delaunay::Vertex_handle v, bool& is_min, bool& is_max) ;
        void draw_dual_facet(Delaunay::Vertex_handle v) ;
        double radius(Delaunay::Vertex_handle v) ;
		void draw_regularity() ;

    private:
        Delaunay* delaunay_ ;
		DelaunayIO* IO_ ;
        DelaunayCVT* CVT_ ;
        GLboolean show_domain_ ;
        GLfloat vertices_size_ ;
		GLboolean vertices_color_ ;
        GLfloat centers_size_ ;
        GLboolean show_primal_mesh_ ;
        GLboolean show_dual_mesh_ ;
        GLboolean colorize_ ;
        GLboolean show_cells_ ;
        GLboolean show_energy_ ;
        GLboolean show_field_ ;
        bool non_convex_ ;
        GLfloat quad_ratio_ ;
		GLboolean show_boundary_cells_ ;
		GLboolean show_pvd_euclidean_ ;
		GLboolean show_copies_ ;
		//GLboolean show_edge_hist_ ;
		GLboolean show_disk_ ;
		GLboolean show_mesh_ ;
		GLboolean show_min_max_ ;
		GLboolean show_inner_voronoi_ ;
		GLboolean show_regularity_ ; 
        GLboolean show_edge_ ;
        GLboolean show_dual_point_ ;
	} ;

} 

#endif
