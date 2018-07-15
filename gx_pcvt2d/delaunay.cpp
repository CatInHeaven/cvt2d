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
 *  You should have received a copy of the GNU General Public License<
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

#include "delaunay.h"
#include <Geex/basics/stopwatch.h>
#include <Geex/basics/file_system.h>
#include <Geex/combinatorics/map.h>
#include <glut_viewer/glut_viewer.h>
#include <fstream>
#include <queue>

#define printf

#define SQR(x) ((x)*(x))

namespace Geex {

    void SegmentDelaunay::insert(int id, const vec2& p1, const vec2& p2, double step_length) {
        vec2 V = p2 - p1 ;
        double l = V.length() ;
        int count = int(l / step_length) ;
        count = gx_max(count, 2) ;
        V = normalize(V) ;
        double delta = l / double(count) ;
        for(int i=0; i<=count; i++) {
			vec2 p = p1 + (delta * double(i)) * V;
			std::cerr << "i = " << i << " p =(" << p.x <<","<< p.y << ")" << " delta = "<< delta << std::endl ;
            insert(id, p) ;
        }
    } 

    int SegmentDelaunay::locate(const vec2& p) {
        return nearest_vertex(to_cgal(p))->index ;
    }

    void SegmentDelaunay::insert(int id, const vec2& p) {
		Point pp = to_cgal(p);
        baseclass::insert(pp)->index = id ;
    }

//____________________________________________________________________________________

	const vec2 Delaunay::v_offset_[4] = {vec2(-1,-1), vec2(1,-1), vec2(1, 1), vec2(-1, 1) } ;
	const vec2 Delaunay::e_offset_[4] = {vec2( 0,-1), vec2(1, 0), vec2(0, 1), vec2(-1, 0) } ;
	//            7 | 6 | 5
	//           -----------
	//            8 | 0 | 4 
	//           -----------
	//            1 | 2 | 3
	const vec2 Delaunay::domain_offset_[9] = {
		vec2(0, 0), vec2(-1, -1), vec2(0, -1), vec2(1, -1), vec2(1, 0), vec2(1, 1), vec2(0, 1), vec2(-1, 1), vec2(-1, 0)
	} ;

    Delaunay::Delaunay() {
        non_convex_mode_ = false ;
        cached_bbox_ = false ;
        load_boundary(FileSystem::get_project_root() + "/gx_pcvt2d/square.line") ;
		opened_ = false ;
        insert_boundary_ = false ;
		perturb_ = 0 ;
		sample_radius_ = 0.1 ;
		sampling_mode_ = PD_VORONOI ;
		map_ = new Map ;
//		editor_ = new IrregEditor(this, map_) ;
		map_edit_mode_ = false ;
		grid_ = nil ;
    }

	Delaunay::~Delaunay() {
		//delete editor_ ;
		delete map_ ;
	} 

    void Delaunay::save(const std::string& filename) {
		std::ofstream out(filename.c_str()) ;
		out.precision(30) ;
		for(Vertex_iterator it = vertices_begin(); it != vertices_end() ; it++) {
			if(is_primary(it)) { // slave vertex
				out << to_geex(it->point()) << std::endl ;
			}
		}

		out.close() ;
    }
    
    void Delaunay::load(const std::string& filename) {
        std::cerr << "loading " << filename << std::endl ;
        clear() ;
        std::ifstream in(filename.c_str()) ;
        if(!in) {
            std::cerr << "could not open file" << std::endl ;
            return ;
        }
        vec2 p ;
		begin_insert() ;
        while(in) {
            in >> p ;
            if(in) {  // we need to do this additional check else we got the last point twice !
				if(in_boundary(p)) {
					insert(p) ;
//					insert(p+vec2(Numeric::random_float64()*1e-8,Numeric::random_float64()*1e-8)) ; 
				}
            }
        } ;
		end_insert(true) ;
		in.close() ;
    }

    void Delaunay::load_boundary(const std::string& filename) {
        cached_bbox_ = false ;
        boundary_.clear() ;
        boundary_convex_.clear() ;
        boundary_.load(filename) ;
		printf("end loading...\n");
//        boundary_.normalize() ;
        if(!non_convex_mode_) {
            for(unsigned int i=0; i<boundary_.size(); i++) {
                boundary_convex_.push_back(boundary_[i].line()) ;
            }
        }

        double x_min, y_min, z_min ;
        double x_max, y_max, z_max ;
        get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
        double dx = x_max - x_min - 1e-3;
        double dy = y_max - y_min - 1e-3;
        double diag = sqrt(dx*dx + dy*dy) ;

        segments_.clear() ;
        for(unsigned int i=0; i<boundary_.size(); i++) {
            segments_.insert(int(i), boundary_[i].vertex[0], boundary_[i].vertex[1], 0.01 * diag) ;
        }
    }


    void Delaunay::get_bbox(
        real& x_min, real& y_min, real& z_min,
        real& x_max, real& y_max, real& z_max
    ) {
        z_min = 0.0 ; z_max = 1.0 ;
        if(!cached_bbox_) {
            x_min_ = y_min_ =  1e30 ;
            x_max_ = y_max_ = -1e30 ;
            if(boundary_.size() == 0) {
                x_min_ = 0.0 ; y_min_ = 0.0 ;
                x_max_ = 1.0 ; y_max_ = 1.0 ;
            }
            for(unsigned int i=0; i<boundary_.size(); i++) {
                for(unsigned int j=0; j<2; j++) {
                    x_min_ = gx_min(x_min_, boundary_[i].vertex[j].x) ;
                    y_min_ = gx_min(y_min_, boundary_[i].vertex[j].y) ;
                    x_max_ = gx_max(x_max_, boundary_[i].vertex[j].x) ;
                    y_max_ = gx_max(y_max_, boundary_[i].vertex[j].y) ;
                }
            }
            cached_bbox_ = true ;
        }
        x_min = x_min_ ; y_min = y_min_ ;
        x_max = x_max_+1e-3 ; y_max = y_max_+1e-3 ;
    }

    inline vec2 random_v() {
        return vec2( 
            Numeric::random_float64(),
            Numeric::random_float64()
        ) ;
    }

	void Delaunay::do_perturb() {
		std::vector<vec2> new_points ;
		for(unsigned int i=0; i<all_vertices_.size(); ++i) {
			if(is_primary(all_vertices_[i])) {
				double angle = Numeric::random_float64()*2*M_PI ;
				vec2 newp = to_geex(all_vertices_[i]->point()) + perturb_*vec2(cos(angle), sin(angle)) ;
				while(!in_boundary(newp)) {
					angle = Numeric::random_float64()*2*M_PI ;
					newp = to_geex(all_vertices_[i]->point()) + perturb_*vec2(cos(angle), sin(angle)) ;
				}
				new_points.push_back(newp) ;
			}
		}

		this->clear() ;
		begin_insert() ;
		for(unsigned i=0; i<new_points.size(); ++i) {
			insert(new_points[i]) ;
		}
		end_insert(true) ;
		std::cout << "end of perturbation" << std::endl ;
	}

	void Delaunay::compute_inner_voronoi(Vertex_handle vh, std::vector<vec2>& poly) {
		DelaunayBase dt ;		
		
		Delaunay::Face_circulator dfcir = incident_faces(vh) ;
		do {
			vec2 p = dual(dfcir) ;
			dt.insert(to_cgal(p)) ;
			dfcir++ ;
		} while(dfcir!=incident_faces(vh)) ;

		Delaunay::Vertex_circulator vcir = adjacent_vertices(vh) ;
		do {
			dt.insert(vcir->point()) ;
			vcir++ ;
		} while(vcir!=adjacent_vertices(vh)) ;

		Vertex_handle vc = dt.insert(vh->point()) ;

		Delaunay::Face_circulator fcir = incident_faces(vc) ;
		do { 
			if(!dt.is_infinite(fcir)) {
				vec2 dc = to_geex(dt.dual(fcir)) ;
				poly.push_back(dc) ;
			}
			fcir++ ;

		} while(fcir!=incident_faces(vc)) ;
	}

	// -------------------------------------- Generate Poisson disk-------------------------
	struct GridCell {
		GridCell() {
			valid = true ;
//			visited = false ;
		}
		bool valid ;
		vec2 point ;
//		bool visited ;
	} ;

	static int get_uv(int i, int res, double& x) { 
		int ii = i ;
		if(ii>=res) {
			ii -= res ;
			x -= 1.0 ;
		} else if(ii<0) {
			ii += res ;
			x += 1.0 ;
		}
		return ii ;
	}

	void Delaunay::sample_poisson_grid(double radius, int nmaxfail) {
		//int    res = ceil(1.0/(radius/sqrt(2.0))) ;
		//bool   done = false ;

		//this->clear() ;

		//grid_ = new GridCell*[res] ;
		//for(int i=0; i<res; ++i) {
		//	grid_[i] = new GridCell[res] ;
		//}

		//std::cerr << "grid resolution " << res << std::endl ;

		//init_sample_grid(res, radius, nmaxfail) ;
		//insert_grid_samples_period(res) ;
		//while(!sample_grid_gaps(res, radius)) {
		//};		

		//std::cerr << "Done! number of points " << all_vertices_.size() << std::endl ;
		//for(int i=0; i<res; ++i) 
		//	delete [] grid_[i] ;
		//delete [] grid_ ;
	}

	bool Delaunay::is_hit_period(vec2& p, int u, int v, int res, double sqr_radius) {
		if(!grid_[u][v].valid) 
			return false ;

		for(int i=u-2; i<u+3; ++i) {
			double x = p.x ;
			int ii=get_uv(i, res, x) ;
			for(int j=v-2; j<v+3; ++j) {
				double y=p.y ;
				int jj = get_uv(j, res, y) ;
				if(!grid_[ii][jj].valid) {
					
					if(distance2(vec2(x, y), grid_[ii][jj].point) < sqr_radius) {
						return false ;
					}
				}
			}
		}
		return true ;
	}

	void Delaunay::init_sample_grid(int res, double radius, int nmaxfail) {
		double gridlen = 1.0/res ;
		double r2 = SQR(radius) ;
		bool   done = false ;
		int    nmiss = 0 ;
		int    npoint = 0 ;
		int    nmaxthrow = 5*SQR(res) ;
		int    nminthrow = SQR(res)/16 ;
		int    nthrow = 0 ;
		
		SystemStopwatch timer ;
		double t0 = timer.now() ;

		while(nmiss<nmaxfail && nthrow<nmaxthrow || nthrow<nminthrow) {
			vec2 p = random_v() ;
			int  u = p.x/gridlen ;
			int  v = p.y/gridlen ;

			u = min(u, res-1) ;
			v = min(v, res-1) ;

			if(is_hit_period(p, u, v, res, r2)) {
				grid_[u][v].valid = false ;
				grid_[u][v].point = p ;
				nmiss = 0 ;
				npoint ++ ;
			} else {
				nmiss ++ ;
			}

			nthrow ++ ;
		} ;

		double t1 = timer.now() ;
		std::cout << "Initial sampling of " << npoint << " points, " << t1-t0 <<" seconds" << std::endl ;
	}

	void Delaunay::insert_grid_samples_clipped(int res) {
		SystemStopwatch timer ;
		double t0 = timer.now() ;
		int idx = 0 ;

		// insert primary points
		begin_insert() ;
		for(int i=0; i<res; ++i) 
			for(int j=0; j<res; ++j) {
				if(!grid_[i][j].valid) {
					baseclass::insert(to_cgal(grid_[i][j].point)) ;
				}
			}
		double t1 = timer.now() ;
		end_insert(false) ;
		glut_viewer_redraw() ;

		std::cout << "Delaunay of the initial sampling " << all_vertices_.size() << " points, " << t1-t0 <<" seconds" << std::endl ;
	}

	void Delaunay::insert_grid_samples_period(int res) {
		int margin = 4 ;

		SystemStopwatch timer ;
		double t0 = timer.now() ;

		this->clear() ;

		// insert primary points
		begin_insert() ;
		for(int i=0; i<res; ++i) 
			for(int j=0; j<res; ++j) {
				if(!grid_[i][j].valid) {
					baseclass::insert(to_cgal(grid_[i][j].point)) ;
				}
			}
		double t1 = timer.now() ;
		// the indices of each vertex can be encoded on the fly.
		end_insert(false) ;
		std::cout << "Delaunay of the initial sampling " << all_vertices_.size() << " points, " << t1-t0 <<" seconds" << std::endl ;

		t1 = timer.now() ;
		insert_copies_poisson(sample_radius(), false) ;
		// insert copies
		double t2 = timer.now() ;
		std::cerr << "Insertion of copies " << t2-t1 << " seconds." << std::endl ;

		// check the area condition
		double total_area=0 ;
		FOR_EACH_FINITE_FACE_DT(Delaunay, this, f) {
			if(!f->dual_outside) {
				double area = triangle_area(to_geex(f->vertex(0)->point()), 
											to_geex(f->vertex(1)->point()), 
											to_geex(f->vertex(2)->point())) ;
				total_area += area ;
			}
		}
		std::cerr << "total area " << total_area << std::endl ;

		glut_viewer_redraw() ;
	}
	void Delaunay::insert_copies_poisson(double radius, bool redraw) {
		unsigned int nv = all_vertices_.size() ;
		vec2 newp ;
		Vertex_handle vh ;

		mirror_vertices_.clear() ;
		opened_ = true ;
		for(unsigned int i=0; i<nv; ++i) {
			int config ;
			vec2 p = to_geex(all_vertices_[i]->point()) ;
			insert_copy_poisson(p, i, radius) ;
		}

		int nb_gaps = 0 ;
		FOR_EACH_FACE_DT(Delaunay, this, it) {
			if(is_infinite(it)) {
//				it->infinite = true ;
				it->dual = vec2(0.0, 0.0) ;
				it->dual_outside = true ;
				it->vertex(0)->dual_infinite = true ;
				it->vertex(1)->dual_infinite = true ;
				it->vertex(2)->dual_infinite = true ;
				it->has_gap = false ;
			} else {
//				it->infinite = false ;
				it->dual = to_geex(baseclass::dual(it)) ;
				it->dual_outside = !in_boundary(it->dual) ;
				it->has_gap = is_dual_gap(it, radius) ;
				if(!it->dual_outside && it->has_gap) {
					nb_gaps ++ ;
				}
			}
		}

		std::cerr << "number of gaps " << nb_gaps << std::endl ;

		opened_ = false ;
	}

//	bool Delaunay::sample_grid_gaps(int res, double radius) {
//		WeightedTree *WT = nil ;
////		int nb_init = points.size() ;
//		int nb_gap_verts = 0 ;
//		std::vector<vec2> new_points ;
//
//		for(unsigned int i=0; i<all_vertices_.size(); ++i) {
//			Vertex_handle it = all_vertices_[i] ;
//			compute_gap_area(it, radius) ;
//			if(it->ex_area > 0) {
//				if(WT==nil) {
//					WT= new WeightedTree(it->index, it->ex_area) ;	
//				}
//				else {
//					WT->add_to_tree(it->index, it->ex_area) ;
//				}
//
//				nb_gap_verts ++ ;
//			}
//		}
//
//		std::cerr << "number of vertices incident to gaps " << nb_gap_verts << std::endl ;
//		if(nb_gap_verts==0) return true ;
//
//		// sample gap
//		double gridlen = 1.0/res ;
//		double r2 = SQR(radius) ;
//		bool   done = false ;
//		int    nmiss = 0 ;
//		int    npoint = 0 ;
//		int    nmaxthrow = 5*nb_gap_verts ;  //SQR(res) ;
//		int    nminthrow = nb_gap_verts/16 ; //SQR(res)/16 ;
//		int    nthrow = 0 ;
//		int    nmaxfail = 100 ;
//		
//		SystemStopwatch timer ;
//		double t0 = timer.now() ;
//
//		while(nmiss<nmaxfail && nthrow<nmaxthrow || nthrow<nminthrow) {
//			int  vidx = WT->select_vertex_from_tree() ;
//			vec2 p = new_vertex_in_voronoi(all_vertices_[vidx], radius) ;
//			int  domain = domain_idx(p) ;
//			p = p - domain_offset_[domain] ;
//			int  u = p.x/gridlen ;
//			int  v = p.y/gridlen ;
//			bool hit = true ;
//
//			u = min(u, res-1) ;
//			v = min(v, res-1) ;
//
//			if(is_hit_period(p, u, v, res, r2)) {
//				grid_[u][v].valid = false ;
//				grid_[u][v].point = p ;
//				nmiss = 0 ;
//				npoint ++ ;
//				//points.push_back(p) ;
//				new_points.push_back(p) ;
//				//// to do: accelerate insertion by indicate the facet
//				//insert_point_poisson(p, radius) ;
//			} else {
//				nmiss ++ ;
//			}
//
//			nthrow ++ ;
//		} ;
//
//		double t1 = timer.now() ;
//		std::cout << "Sampling gaps " << npoint << " points, " << t1-t0 <<" seconds" << std::endl ;
//
//		// to do: accelerate insertion by indicate the facet
//		// this can be done in constant time since we know the location of new sample	
//		for(unsigned int i=0; i<new_points.size(); ++i) {
//			insert_point_poisson(new_points[i], radius) ;
//		}
//
//		delete WT ;
//
//		return new_points.size() == 0 ;
//	}

	void Delaunay::insert_point_poisson(vec2& p, double radius) {
		opened_ = true ;
		Vertex_handle vh = baseclass::insert(to_cgal(p)) ;
		vh->index = all_vertices_.size() ;
		vh->domain = 0 ;
		all_vertices_.push_back(vh) ;
		update_neighbor_cells(vh, radius) ;
		insert_copy_poisson(p, vh->index, radius) ;
		opened_ = false ;	
	}

	void Delaunay::update_neighbor_cells(Vertex_handle& vh, double radius) {
		Face_circulator fcir = incident_faces(vh) ;

		do {
			if(is_infinite(fcir)) {
//                fcir->infinite = true ;
                fcir->dual = vec2(0.0, 0.0) ;
                fcir->dual_outside = true ;
				fcir->has_gap = false ;
            }
			else {
 //               fcir->infinite = false ;
                fcir->dual = to_geex(baseclass::dual(fcir)) ;
                fcir->dual_outside = !in_boundary(fcir->dual) ;			
				fcir->has_gap = is_dual_gap(fcir, radius) ;
			}
			++fcir ;
		} while(fcir!=incident_faces(vh)) ;
	}

	void Delaunay::insert_copy_poisson(vec2& p, int i, double radius) {
		double margin = radius*2.0 ;
		vec2   newp ;
		Vertex_handle vh ;

		if(p.x < margin) {
			newp = p + vec2(1.0, 0) ;
			vh = insert(newp) ;
			vh->index = i ;
			vh->domain = 4 ;
			mirror_vertices_[i].push_back(vh) ;
			update_neighbor_cells(vh, radius) ;
		}
		if(p.x > 1.0-margin) {
			newp = p - vec2(1.0, 0) ;
			vh = insert(newp) ;
			vh->index = i ;
			vh->domain = 8 ;
			mirror_vertices_[i].push_back(vh) ;	
			update_neighbor_cells(vh, radius) ;
		}
		if(p.y < margin) {
			newp = p + vec2(0, 1) ;
			vh = insert(newp) ;
			vh->index = i ;
			vh->domain = 6 ;
			mirror_vertices_[i].push_back(vh) ;
			update_neighbor_cells(vh, radius) ;
		}
		if(p.y > 1.0-margin) {
			newp = p - vec2(0, 1) ;
			vh = insert(newp) ;
			vh->index = i ;
			vh->domain = 2 ;
			mirror_vertices_[i].push_back(vh) ;
			update_neighbor_cells(vh, radius) ;
		}
		//
		if(p.x < margin && p.y < margin) {
			newp = p + vec2(1.0, 1.0) ;
			vh = insert(newp) ;
			vh->index = i ;
			vh->domain = 5 ;
			mirror_vertices_[i].push_back(vh) ;
			update_neighbor_cells(vh, radius) ;
		}
		if(p.x>1.0-margin && p.y<margin) {
			newp = p + vec2(-1.0, 1.0) ;
			vh = insert(newp) ;
			vh->index = i ;
			vh->domain = 7 ;
			mirror_vertices_[i].push_back(vh) ;
			update_neighbor_cells(vh, radius) ;
		}
		if(p.x>1.0-margin && p.y>1.0-margin) {
			newp = p + vec2(-1.0, -1.0) ;
			vh = insert(newp) ;
			vh->index = i ;
			vh->domain = 1 ;
			mirror_vertices_[i].push_back(vh) ;
			update_neighbor_cells(vh, radius) ;
		}
		if(p.x<margin && p.y>1.0-margin) {
			newp = p + vec2(1.0, -1.0) ;
			vh = insert(newp) ;
			vh->index = i ;
			vh->domain = 3 ;
			mirror_vertices_[i].push_back(vh) ;
			update_neighbor_cells(vh, radius) ;
		}
	}

	void Delaunay::generate_poisson_disk() {
		std::vector<vec2>   points ;

		//switch(sampling_mode_) {
		//case PD_CCVT:
		//	sample_ccvt(points, 1024) ;
		//	set_vertices(points) ;
		//	break ;
		//case PD_VORONOI:
		//	fast_poisson_disk(sample_radius()) ;
		//	break; 
		//case PD_BOUNDARY:
		//	break; 
		//case PD_DARTTHROW: {
		//	sample_poisson_grid(sample_radius(), 400) ;
		//	break ;
		//}			
		//case PD_PENROSE:
		//	break; 
		//default: 
		//	std::cerr << "sampling mode is not implemented...\n" ;
		//	break ;
		//}
	}

	void Delaunay::fast_poisson_disk(double exclude_radius) {
//		vertices_.clear() ;
//		all_vertices_.clear() ;
//		clear() ;
//
//		srand(time(NULL)) ;
//
//		std::cerr << "sampling radius: " << exclude_radius << std::endl ;
//
////		vec2 p0(0.68059327982421336, 0.30582598345896789) ;
////		std::cerr << "first point " << p0 << std::endl ;
//		vec2 p0 = random_v() ;
//		Vertex_handle v0 = insert_vertex_periodic(p0, exclude_radius) ;
//		gx_assert(v0->ex_area > 0) ;
//		WeightedTree *WT = new WeightedTree(v0->index, v0->ex_area);
//		SystemStopwatch timer ;
//		double t0 = timer.now() ;
//
//		while(!WT->no_area_left()) {
//			int  v = WT->select_vertex_from_tree() ;
//			if(all_vertices_[v]->ex_area == 0.0) {
//				std::cerr << "Voronoi cell with 0 area is selected..." << std::endl ;
//			}
//			gx_assert(all_vertices_[v]->ex_area > 0.0) ;
//
//			vec2 newv = new_vertex_in_voronoi(all_vertices_[v], exclude_radius) ;
//			int domain = domain_idx(newv) ;
//			newv = newv - domain_offset_[domain] ;
//			Vertex_handle vh = insert_vertex_periodic(newv, exclude_radius) ;
//			WT->add_to_tree(vh->index, vh->ex_area) ;
//
//			Vertex_circulator cir = incident_vertices(vh) ;
//			do {
//			  printf("checking neighbor %f %f   %f  ", 
//					 to_geex(cir->point()) , 
//					 distance(newv, to_geex(cir->point())));
//			  if(distance(newv, to_geex(cir->point())) < exclude_radius) {
//				  std::cerr << "vertex " << v << " is selected, new vertex " << newv <<" is generated" << std::endl ;
//				  std::cerr << "checking neighbor " << cir->index << " " << to_geex(cir->point()) << std::endl ;
//			  }
//			  gx_assert(distance(newv, to_geex(cir->point())) >= exclude_radius);
//			  if (cir->domain==0 && WT->in_tree(cir->index))
//				  WT->update_tree(cir->index, cir->ex_area);
//			  printf("neighbor %d has been updated\n", cir->index);
//
//			  ++cir ;
//			} while(cir!=incident_vertices(vh)) ;
//
//			// the neighbor affected by periodic copies also need to be updated
//			for(int k=1; k<9; ++k) {
//				Vertex_handle vk = vertices_[k][vh->index] ;
//				cir = incident_vertices(vk) ;
//				do {
//					if(cir->domain==0 && WT->in_tree(cir->index)) {
//						WT->update_tree(cir->index, cir->ex_area);
//					}
//					++cir ;
//				} while(cir!=incident_vertices(vk)) ;
//			}
//
//			if ((all_vertices_.size() % 10000) == 0) {
//				double t1 = timer.now() ;
//				std::cerr << "10000 points added in " << t1-t0 << " seconds (ratio "
//					<< (t1-t0) / log((double)all_vertices_.size()) <<") total " << all_vertices_.size() <<", " 
//					<< WT->area_left() <<" arealeft" << std::endl ;
//				t0 = t1 ;
//			}
//		}
//
//		std::cerr << "Done. " << all_vertices_.size() << " points generated..."
//			<< this->number_of_vertices() << " in DT..." << std::endl ;
//
//		delete WT ;
	}

	bool Delaunay::is_dual_gap(Delaunay::Face_handle f, double exclude_radius) {
		return distance(f->dual, to_geex(f->vertex(0)->point())) > exclude_radius &&
			   distance(f->dual, to_geex(f->vertex(1)->point())) > exclude_radius &&
		       distance(f->dual, to_geex(f->vertex(2)->point())) > exclude_radius ;
	}
	
    static inline int gx_sign(double x) {
        return x > 0 ? 1 : (x < 0 ? -1 : 0) ;
    }


	static bool in_triangle(vec2& v0, vec2& v1, vec2& v2, vec2& p) {
		int s01 = gx_sign(det(v0-p, v1-p)) ;
		int s12 = gx_sign(det(v1-p, v2-p)) ;
		int s20 = gx_sign(det(v2-p, v0-p)) ;
		return (s01*s12>=0 && s12*s20>=0) ;
	}

	void Delaunay::update_gaps(double exclude_radius) {
		void_regions_.clear() ;
		void_skeleton_.clear() ;

		for(Delaunay::Finite_faces_iterator it = finite_faces_begin(); it!=finite_faces_end(); ++it) {
			if(!is_infinite(it)) {
				it->dual = to_geex(baseclass::dual(it)) ;
				it->has_gap = is_dual_gap(it, exclude_radius) ;
				it->dual_outside = !in_boundary(it->dual) ;
			}
		}

		int nvoids = 0 ;
		for(Delaunay::Finite_faces_iterator it = finite_faces_begin(); it!=finite_faces_end(); ++it) {
			if(!it->visited && !is_infinite(it)) {
				if(it->has_gap && !it->dual_outside) {
					compute_void_skeleton(it, exclude_radius) ;
					compute_void_region(it, exclude_radius) ;
					nvoids ++ ;
				}
			}
		}

		std::cerr << "number of gaps " << nvoids << std::endl ;

		// reset visit flag
		for(Delaunay::Finite_faces_iterator it = finite_faces_begin(); it!=finite_faces_end(); ++it) {
			it->visited = false ; 
		}
	}

	static inline void get_opposite_vertex(Delaunay::Face_handle f, int vidx, vec2& ov) {
		Delaunay::Face_handle fopp = f->neighbor(vidx) ;
		for(int i=0; i<3; ++i) {
			if(fopp->neighbor(i)==f) {
				ov = to_geex(fopp->vertex(i)->point()) ;
				return ;
			}
		}
		gx_assert(false) ; // cannot reach
	}

	void Delaunay::compute_void_region(Face_handle f, double exclude_radius) {
		std::vector<vec2> poly ;
		bool same_side = false ;

		for(int i=0; i<3; ++i) {
			Face_handle fopp = f->neighbor(i) ;
			vec2 v0 = to_geex(f->vertex(f->ccw(i))->point()) ;
			vec2 v1 = to_geex(f->vertex(f->cw(i))->point()) ;
			vec2 v2 = to_geex(f->vertex(i)->point()) ;
			vec2 v3 ;
			get_opposite_vertex(f, i, v3) ;
			double elen = distance(v0, v1) ;
			vec2 dir = normalize(v1-v0) ;
			vec2 nor = vec2(-dir.y, dir.x) ;
			Line<double> L(v0, nor) ;

			if(L.side(f->dual)*L.side(fopp->dual)<0) {
				if(elen < 2*exclude_radius) {
					double h = sqrt(SQR(exclude_radius)-SQR(elen/2.0)) ;
					vec2 mp = 0.5*(v0 + v1) ;
					vec2 newp = mp + h*nor ;
					if(in_triangle(v0, v1, v2, newp)) {
						poly.push_back(newp) ;
					}
				} else {
					vec2 newv0 = v0 + exclude_radius*dir ;
					vec2 newv1 = v1 - exclude_radius*dir ;
					poly.push_back(newv0) ;
					poly.push_back(newv1) ;
				}
			}
			else {
				// same_side = true ;
				if(L.side(f->dual) < 0) { // obuse triangle
					vec2 newv0 = v0 + exclude_radius*dir ;
					vec2 newv1 = v1 - exclude_radius*dir ;
					bool out0 = distance(newv0, v2) > exclude_radius ;
					bool out1 = distance(newv1, v2) > exclude_radius ;
					if(out0 && out1) {
						poly.push_back(newv0) ;
						poly.push_back(newv1) ;
					} 
					else if(out0) {
						double t = dot(v2-v1, -dir) ;
						double h = det(v2-v1, -dir) ;
						gx_assert(h>0 && t>0) ;
						newv1 = v1 - (t+sqrt(SQR(exclude_radius)-SQR(h)))*dir ;
						poly.push_back(newv0) ;
						poly.push_back(newv1) ;
					}
					else if(out1) {
						double t = dot(v2-v0, dir) ;
						double h = det(dir, v2-v0) ;
						gx_assert(h>0 && t>0) ;
						newv0 = v0 + (t+sqrt(SQR(exclude_radius)-SQR(h)))*dir ;
						poly.push_back(newv0) ;
						poly.push_back(newv1) ;
					} 
					else {
						// both newv0 and newv1 are enclosed by circle of v2
					}	
				} else { // the neighbor of an obuse triangle
					if(elen < 2*exclude_radius) {
						double h = sqrt(SQR(exclude_radius)-SQR(elen/2.0)) ;
						vec2 mp = 0.5*(v0 + v1) ;
						vec2 newp = mp + h*nor ;
						if(distance(newp, v3) > exclude_radius) {
							poly.push_back(newp) ;
						} else	{
							double len03 = distance(v0, v3) ;
							double h03 = sqrt(SQR(exclude_radius)-SQR(len03/2.0)) ;
							vec2 mp03 = 0.5*(v0 + v3) ;
							vec2 dir03 = normalize(v3 - v0) ;
							vec2 newv03 = mp03 + h03*vec2(-dir03.y, dir03.x) ;

							double len31 = distance(v1, v3) ;
							double h31 = sqrt(SQR(exclude_radius)-SQR(len31/2.0)) ;
							vec2 mp31 = 0.5*(v1 + v3) ;
							vec2 dir31 = normalize(v1 - v3) ;
							vec2 newv31 = mp31 + h31*vec2(-dir31.y, dir31.x) ;

							poly.push_back(newv03) ;
							poly.push_back(newv31) ;
						}

					} else {
						vec2 newv0 = v0 + exclude_radius*dir ;
						vec2 newv1 = v1 - exclude_radius*dir ;
						// test the intersection points with the obuse vertex
						bool out0 = distance(newv0, v3) > exclude_radius ;
						bool out1 = distance(newv1, v3) > exclude_radius ;
						if(out0 && out1) {
							poly.push_back(newv0) ;
							poly.push_back(newv1) ;
						} 
						else if(out0) {
							double t = dot(v3-v1, -dir) ;
							double h = det(v3-v1, dir) ;
							gx_assert(h>0 && t>0) ;
							newv1 = v1 - (t+sqrt(SQR(exclude_radius)-SQR(h)))*dir ;
							// add an additional point
							double len31 = distance(v1, v3) ;
							double h31 = sqrt(SQR(exclude_radius)-SQR(len31/2.0)) ;
							vec2 mp31 = 0.5*(v1 + v3) ;
							vec2 dir31 = normalize(v1 - v3) ;
							vec2 newv31 = mp31 + h31*vec2(-dir31.y, dir31.x) ;

							poly.push_back(newv0) ;
							poly.push_back(newv1) ;
							poly.push_back(newv31) ;
						}
						else if(out1) {
							double t = dot(v3-v0, dir) ;
							double h = det(v3-v0, dir) ;
							gx_assert(h>0 && t>0) ;
							newv0 = v0 + (t+sqrt(SQR(exclude_radius)-SQR(h)))*dir ;
							// add an additional point
							double len03 = distance(v0, v3) ;
							double h03 = sqrt(SQR(exclude_radius)-SQR(len03/2.0)) ;
							vec2 mp03 = 0.5*(v0 + v3) ;
							vec2 dir03 = normalize(v3 - v0) ;
							vec2 newv03 = mp03 + h03*vec2(-dir03.y, dir03.x) ;

							poly.push_back(newv03) ;
							poly.push_back(newv0) ;
							poly.push_back(newv1) ;

						} else {
							double len03 = distance(v0, v3) ;
							double h03 = sqrt(SQR(exclude_radius)-SQR(len03/2.0)) ;
							vec2 mp03 = 0.5*(v0 + v3) ;
							vec2 dir03 = normalize(v3 - v0) ;
							vec2 newv03 = mp03 + h03*vec2(-dir03.y, dir03.x) ;

							double len31 = distance(v1, v3) ;
							double h31 = sqrt(SQR(exclude_radius)-SQR(len31/2.0)) ;
							vec2 mp31 = 0.5*(v1 + v3) ;
							vec2 dir31 = normalize(v1 - v3) ;
							vec2 newv31 = mp31 + h31*vec2(-dir31.y, dir31.x) ;

							poly.push_back(newv03) ;
							poly.push_back(newv31) ;
						}
					}
				}
			}
		}

		if(poly.size() > 2 /*&& !same_side*/) {
			void_regions_.push_back(poly) ;
		}
	}

	void Delaunay::compute_void_skeleton(Delaunay::Face_handle f, double exclude_radius) {
		std::vector<vec2> poly ;
		if(f->visited) return ;
		//Edge estart(f, 0), ecur(f, 0) ;
		//Face_handle curf = f ;

		for(int i=0; i<3; ++i) {
			Edge e(f, i) ;
			Face_handle fopp = f->neighbor(e.second) ;
			vec2 v0 = to_geex(f->vertex(f->ccw(e.second))->point()) ;
			vec2 v1 = to_geex(f->vertex(f->cw(e.second))->point()) ;
			double elen = distance(v0, v1) ;
			if(fopp->has_gap) { 
				vec2 n = normalize(v1-v0) ;
				n = vec2(-n.y, n.x) ;
				Line<double> L(v0, n) ;
				if(elen>2*exclude_radius || (L.side(f->dual)*L.side(fopp->dual)>0))  {
					poly.push_back(f->dual) ;
					poly.push_back(fopp->dual) ;
				} else {
					double h = sqrt(SQR(exclude_radius)-SQR(elen/2.0)) ;
					vec2 mp = 0.5*(v0 + v1) ;
					vec2 dir = normalize(f->dual-mp) ;
					vec2 newp = mp + h*dir ;
					poly.push_back(newp) ;
					poly.push_back(f->dual) ;
				}
			}
			else {
				double h = sqrt(SQR(exclude_radius)-SQR(elen/2.0)) ;
				vec2 mp = 0.5*(v0 + v1) ;
				vec2 dir = normalize(f->dual-mp) ;
				vec2 newp = mp + h*dir ;
				poly.push_back(newp) ;
				poly.push_back(f->dual) ;
			}
		}

		if(poly.size()>0) {
			void_skeleton_.push_back(poly) ;
		}
	}

	void Delaunay::compute_void_region(Delaunay::Face_handle f, double exclude_radius, std::vector<vec2>& poly) {
		vec2 v0 = to_geex(f->vertex(0)->point()) ;
		vec2 v1 = to_geex(f->vertex(1)->point()) ;
		vec2 v2 = to_geex(f->vertex(2)->point()) ;

		gx_assert(det(v1-v0, v2-v0)>0) ; // f is ccw 
		
		poly.clear() ;
		for(int i=0; i<3; ++i) {
			int j = (i+1)%3 ;
			int k = (i+2)%3 ;
			vec2 vi = to_geex(f->vertex(i)->point()) ;
			vec2 vj = to_geex(f->vertex(j)->point()) ;
			vec2 vk = to_geex(f->vertex(k)->point()) ;
			vec2 vo ;
			Face_handle fo = f->neighbor(k) ;
			for(int ii=0; ii<3; ++ii) {
				if(fo->vertex(ii) != f->vertex(i) && fo->vertex(ii) != f->vertex(j)) {
					vo = to_geex(fo->vertex(ii)->point()) ;
					break ;
				}
			}
			double distij = distance(vi, vj) ;
			vec2 dir = vj - vi ;
			dir = normalize(dir) ;

			if(distij > 2*exclude_radius) {
				vec2 newvi = vi + exclude_radius*dir ;
				vec2 newvj = vj - exclude_radius*dir ;

				// test with the other vertex vk
				bool newvi_out = distance(newvi, vk) > exclude_radius ;
				bool newvj_out = distance(newvj, vk) > exclude_radius ;
				if(newvi_out && newvj_out) {
					poly.push_back(newvi) ;
					poly.push_back(newvj) ;
				} else if(newvi_out) {
					double diskjk = distance(vj, vk) ;
					vec2 dirjk = vk-vj ;
					dirjk = normalize(dirjk) ;
					double cosa = dot(dirjk, -dir) ;
					double t = diskjk * cosa ;
					double h = sqrt(SQR(diskjk)-SQR(t)) ;
					newvj = vj - (t + sqrt(SQR(exclude_radius)-SQR(h)))*dir ;

					poly.push_back(newvi) ;
					poly.push_back(newvj) ;
				} else if(newvj_out) {
					double diskik = distance(vi, vk) ;
					vec2 dirik = vk-vi ;
					dirik = normalize(dirik) ;
					double cosa = dot(dirik, dir) ;
					double t = diskik * cosa ;
					double h = sqrt(SQR(diskik)-SQR(t)) ;
					newvi = vi + (t + sqrt(SQR(exclude_radius)-SQR(h)))*dir ;

					poly.push_back(newvi) ;
					poly.push_back(newvj) ;
				} else {
					// the triangle is fully covered by 3 circles
				//	break ;
				}

				// test with the opposite vertex vo


			}
			else {
				double h = sqrt(SQR(exclude_radius)-SQR(distij/2.0)) ;
				vec2 newv = 0.5*(vi+vj) + h*vec2(-dir.y, dir.x) ;
				if(in_triangle(vi, vj, vk, newv)) {
					poly.push_back(newv) ;
				}
			}
		}
//		if(poly.size()<3) poly.push_back(f->dual) ;
	}

	Delaunay::Vertex_handle Delaunay::insert_vertex_periodic(vec2& p, double exclude_radius) {
//		vec2 offset[] = {vec2(0, 0), vec2(-1, -1), vec2(0, -1), vec2(1, -1), vec2(1, 0), vec2(1, 1), vec2(0, 1), vec2(-1, 1), vec2(-1, 0)} ;
		int  domain = domain_idx(p) ;
		vec2 newp = p - domain_offset_[domain] ;
		Vertex_handle v0 = baseclass::insert(to_cgal(newp)) ;
		v0->index = all_vertices_.size() ;
		v0->domain = 0 ; //-1 ;
		
		for(int i=1; i<9; ++i) {
			Vertex_handle vh = baseclass::insert(to_cgal(newp+domain_offset_[i])) ;
			vh->index = v0->index ; //i*nb_sample_ + all_vertices_.size() ;
			vh->domain = i ; //v0->index ;
			vertices_[i].push_back(vh) ;
		}

		// update nieghbor cells of the primary point
		Face_circulator fcir = incident_faces(v0) ;
		do {
			if(is_infinite(fcir)) {
 //               fcir->infinite = true ;
                fcir->dual = vec2(0.0, 0.0) ;
                fcir->dual_outside = true ;
				fcir->has_gap = false ;
            }
			else {
 //               fcir->infinite = false ;
                fcir->dual = to_geex(baseclass::dual(fcir)) ;
                fcir->dual_outside = !in_boundary(fcir->dual) ;			
				fcir->has_gap = is_dual_gap(fcir, exclude_radius) ;
			}
			++fcir ;
		} while(fcir!=incident_faces(v0)) ;

		// update neighbor cells the mirror points
		for(int i=1; i<9; ++i) {
			Vertex_handle vi = vertices_[i][v0->index] ;
			fcir = incident_faces(vi) ;
			do {
				if(is_infinite(fcir)) {
//					fcir->infinite = true ;
					fcir->dual = vec2(0.0, 0.0) ;
					fcir->dual_outside = true ;
					fcir->has_gap = false ;
				}
				else {
//					fcir->infinite = false ;
					fcir->dual = to_geex(baseclass::dual(fcir)) ;
					fcir->dual_outside = !in_boundary(fcir->dual) ;			
					fcir->has_gap = is_dual_gap(fcir, exclude_radius) ;
				}
				++ fcir ;
			} while(fcir!=incident_faces(vi)) ;
		}

		// update excluded area 
//		compute_gap_area(v0, exclude_radius) ;
		compute_voronoi_area(v0, exclude_radius) ;
		all_vertices_.push_back(v0) ;

		// update one ring of v0
		Vertex_circulator cir = adjacent_vertices(v0) ;
		do {
			// update area for primay points
			if(is_primary(cir)) { 
				compute_voronoi_area(all_vertices_[cir->index], exclude_radius) ;
				//compute_gap_area(all_vertices_[cir->index], exclude_radius) ;
			}
			++ cir ;
		} while(cir!=adjacent_vertices(v0)) ;

		// update one ring of copies
		for(int i=1; i<9; ++i) {
			Vertex_handle vi = vertices_[i][v0->index] ;
			cir = adjacent_vertices(vi) ;
			do {
				// update area for primay points
				if(is_primary(cir)) { 
					compute_voronoi_area(all_vertices_[cir->index], exclude_radius) ;
					//compute_gap_area(all_vertices_[cir->index], exclude_radius) ;
				}
				++ cir ;
			} while(cir!=adjacent_vertices(vi)) ;
		}

//		glut_viewer_redraw() ;
		return v0 ;
	}

	void Delaunay::compute_voronoi_area(Delaunay::Vertex_handle& v, double exclude_radius) {
//		Delaunay::Face_circulator it = incident_faces(v) ;
//		Delaunay::Face_circulator jt = it ; jt++ ;
//		vec2 p0 = to_geex(v->point()) ;
//
//		double area =0 ;
//		do {
//			const vec2& p1 = it->dual ;
//			const vec2& p2 = jt->dual ;
////			double Vi = triangle_area(p0, p1, p2) ;
//			double tarea = triangle_area_exclude(p0, p1, p2, exclude_radius) ;
//			area += tarea ;
//			it++ ; jt++ ;
//		} while(it != incident_faces(v)) ;
//
//		v->ex_area = area ;
	}

//	vec2 Delaunay::new_vertex_in_voronoi(Delaunay::Vertex_handle& v, double exclude_radius) {
	 // double area = 0.0 ;
	 // vec2 vout ;
	 // std::vector<vec2> last_positive ;
	 // bool done = false ;
	
		//Delaunay::Face_circulator it = incident_faces(v) ;
		//Delaunay::Face_circulator jt = it ; jt++ ;
		//vec2 p0 = to_geex(v->point()) ;
		//do {
		//	const vec2& p1 = it->dual ;
		//	const vec2& p2 = jt->dual ;
		//	double Vi = triangle_area(p0, p1, p2) ;
		//	double tarea = triangle_area_exclude(p0, p1, p2, exclude_radius) ;
		//	area += tarea ;
		//	if (tarea >0.0) {
		//		last_positive.clear() ;
		//		last_positive.push_back(p0) ;
		//		last_positive.push_back(p1) ;
		//		last_positive.push_back(p2) ;
		//	}
		//	it++ ; jt++ ;
		//} while(it != incident_faces(v)) ;

		//if(last_positive.size()==0) {
		//	std::cerr << "last positive is 0, recomputing..." << std::endl ;
		//	it = incident_faces(v) ;
		//	jt = it ; jt++ ;
		//	vec2 p0 = to_geex(v->point()) ;
		//	do {
		//		const vec2& p1 = to_geex(baseclass::dual(it)) ;//it->dual ;
		//		const vec2& p2 = to_geex(baseclass::dual(jt)) ;//jt->dual ;
		//		double tarea = triangle_area_exclude(p0, p1, p2, exclude_radius) ;
		//		area += tarea ;
		//		if (tarea >0.0) {
		//			last_positive.clear() ;
		//			last_positive.push_back(p0) ;
		//			last_positive.push_back(p1) ;
		//			last_positive.push_back(p2) ;
		//		}
		//		it++ ; jt++ ;
		//	} while(it != incident_faces(v)) ;
		//	std::cerr << "recomputed area = " << area << std::endl ;
		//}
		//gx_assert(last_positive.size() > 0) ;

		//area *= Numeric::random_float64() ; //drand48();

		//it = incident_faces(v) ;
		//jt = it ; jt++ ;
		//while(jt!=incident_faces(v) && (area >0.0)) {
		//	const vec2& p1 = it->dual ;
		//	const vec2& p2 = jt->dual ;
		//	double tarea = triangle_area_exclude(p0, p1, p2, exclude_radius) ;
		//	area -= tarea ;// triangle_area_exclude(p0, p1, p2, exclude_radius) ;
		//	if (area <= 0.0) {
		//	  vout = random_point_in_triangle_exclude(p0, p1, p2, exclude_radius);
		//	  done = true ;
		//	  break ;
		//	}
		//	it++ ; jt++ ;
		//} 

	 //// l = voronoi_region;
	 //// while (l && (area > 0.0)) {
		//////area -= triangle_area_exclude(GTS_TRIANGLE(l->data), v, exclude_radius);
		////  area -= triangle_area_exclude(p0, p1, p2, 
		////if (area <= 0.0) {
		////  vout = random_point_in_triangle_exclude(GTS_TRIANGLE(l->data), v, exclude_radius);
		////  done = true ;
		////}
		////l = l->next;
	 //// }

	 // if (! done) {
		//vout = random_point_in_triangle_exclude(last_positive[0], last_positive[1], last_positive[2], exclude_radius);
	 // }

	 //// while (l) {
		////gts_object_destroy(GTS_OBJECT(l->data));
		////l = l->next;
	 //// }
	 //// g_slist_free(voronoi_region);

	 // return (vout);
//	}


	// new algorithm replace the original one

	static inline double compute_sector_gap_area(vec2& p0, vec2& p1, vec2& p2, double exclude_radius) {
		vec2 v01 = p1-p0 ;
		vec2 v02 = p2-p0 ;
		double d01 = v01.length() ;
		double d02 = v02.length() ;
		double T = det(v01, v02) ;
		//gx_assert(T>0) ;
		double valsin = T/(d01*d02) ;
		if(valsin>1) valsin=1 ;
		if(valsin<0) valsin=0 ;
		double angle = asin(valsin) ;
		double garea = 0.5*T-0.5*SQR(exclude_radius)*angle ;
		return garea ;	
	}

	void Delaunay::compute_gap_area(Delaunay::Vertex_handle& v, double exclude_radius) {
		Face_circulator f1 = incident_faces(v) ;
		Face_circulator f2 = f1 ; f2++ ;
		vec2 p0 = to_geex(v->point()) ;
		double area = 0 ;

		do {
			if(f1->has_gap && f2->has_gap) {
				vec2 p1 = f1->dual ;
				vec2 p2 = f2->dual ;
				Vertex_handle vopp = f1->vertex(f1->cw(f1->index(v))) ;
				vec2 p3 = to_geex(vopp->point()) ;
				vec2 dir = normalize(p3-p0) ;
				vec2 nor = vec2(-dir.y, dir.x) ;
				Line<double> L(p0, nor) ;

				// p1 p2 on same side
				if(L.side(p1) * L.side(p2) > 0) {
					double garea = compute_sector_gap_area(p0, p1, p2, exclude_radius) ;
					area += garea ;
				}
				else {
					double d03 = distance(p0, p3) ;
					if(d03>2*exclude_radius) {
						double garea = compute_sector_gap_area(p0, p1, p2, exclude_radius) ;
						area += garea ;
					} else {
						vec2 mp = 0.5*(p0+p3) ;
						double h = sqrt(SQR(exclude_radius)-SQR(d03*0.5)) ;
						vec2 p12 = mp - h*nor;
						vec2 p21 = mp + h*nor;
						double garea1 = compute_sector_gap_area(p0, p1, p12, exclude_radius) ;
						double garea2 = compute_sector_gap_area(p0, p21, p2, exclude_radius) ;
						area += garea1 ;
						area += garea2 ;
					}
				}								
			} 
			else if(f1->has_gap && !f2->has_gap) {
				vec2 p1 = f1->dual ;
				Vertex_handle vopp = f1->vertex(f1->cw(f1->index(v))) ;
				vec2 p3 = to_geex(vopp->point()) ;
				vec2 dir = normalize(p3-p0) ;
				vec2 nor = vec2(-dir.y, dir.x) ;
				double d03 = distance(p0, p3) ;
				vec2 mp = 0.5*(p0+p3) ;
				double h = sqrt(SQR(exclude_radius)-SQR(d03*0.5)) ;
				vec2 p12 = mp - h*nor;
				double garea = compute_sector_gap_area(p0, p1, p12, exclude_radius) ;
				area += garea ;
			}
			else if(!f1->has_gap && f2->has_gap) {
				vec2 p2 = f2->dual ;
				Vertex_handle vopp = f1->vertex(f1->cw(f1->index(v))) ;
				vec2 p3 = to_geex(vopp->point()) ;
				vec2 dir = normalize(p3-p0) ;
				vec2 nor = vec2(-dir.y, dir.x) ;
				double d03 = distance(p0, p3) ;
				vec2 mp = 0.5*(p0+p3) ;
				double h = sqrt(SQR(exclude_radius)-SQR(d03*0.5)) ;
				vec2 p21 = mp + h*nor;
				double garea = compute_sector_gap_area(p0, p21, p2, exclude_radius) ;
				area += garea ;
			}
			else {
				// no gap
			}
			++f1 ;
			++f2 ;
		} while(f1!=incident_faces(v)) ;

		v->ex_area = area ;
	}

	vec2 Delaunay::new_vertex_in_gap(Delaunay::Vertex_handle& v, double exclude_radius) {
		return vec2(0, 0) ;
	}

	void Delaunay::set_vertices(std::vector<vec2>& points) {
		SystemStopwatch timer ;
		double t0 = timer.now() ;

		clear() ;
		begin_insert() ;
		for(unsigned int i=0; i<points.size(); ++i) {
			insert(points[i]) ;
		}
		end_insert() ;
		double t1 = timer.now() ;
		std::cout << "Delaunay triangulation of " << this->number_of_vertices() << " vertices, " 
			<< t1-t0 << " second. " << std::endl ;

//		insert_copies(true, true) ;
		insert_copies_poisson(this->sample_radius(), true) ;
		double t2 = timer.now() ;
		std::cout << "Insertion  mirror vertices, " << t2-t1 << " second. " << std::endl ;
	}

	void Delaunay::digree_hist(double& plt5, double& p5, double& p6, double& p7, double& pgt7) {
		// out put vertex degrees 
		int nb_lt5=0, nb_gt7=0, nb_5=0, nb_6=0, nb_7=0 ;
		for(unsigned int i=0; i<all_vertices_.size(); ++i) {

			if(on_convex_hull(all_vertices_[i]) )
				continue ;

			if(is_primary(all_vertices_[i])) {
				if(degree(all_vertices_[i]) < 5)
					nb_lt5++ ;
				else if(degree(all_vertices_[i]) == 5)
					nb_5++ ;
				else if(degree(all_vertices_[i]) == 6)
					nb_6++ ;
				else if(degree(all_vertices_[i]) == 7)
					nb_7++ ;
				else //(degree(all_vertices_[i]) > 7)
					nb_gt7++ ;
			}
		}
		double total = nb_lt5 + nb_gt7 + nb_5 + nb_6 + nb_7 ;
		plt5 = nb_lt5/total ;
		p5 = nb_5/total ;
		p6 = nb_6/total ;
		p7 = nb_7/total ;
		pgt7 = nb_gt7/total ;
	}

	bool Delaunay::on_convex_hull(Delaunay::Vertex_handle& v) {
		Vertex_circulator cir = adjacent_vertices(v) ;
		do {
			if(is_infinite(cir))
				return true ;
			++cir ;
		} while(cir !=adjacent_vertices(v)) ;
		return false ;
	}

	void Delaunay::edge_hist(double& p66, double& p56, double& p76, double& p55, double& p77, double& p57) {
		int n66=0, n56=0, n76=0, n55=0, n77=0, n57=0, total=0 ;
		double minlen=1e10, maxlen=-1e10 ;
		for(unsigned int i=0; i<all_vertices_.size(); ++i) {
			// exclude points on convex hull for statics
			if(on_convex_hull(all_vertices_[i]) )
				continue ;

			if(is_primary(all_vertices_[i]) && !is_infinite(all_vertices_[i])) {
				int deg = degree(all_vertices_[i]) ;
				Vertex_circulator cir = adjacent_vertices(all_vertices_[i]) ;
				do {
					int d2 = degree(cir) ;
					if(deg==6) {
						if(d2==5) n56++ ;
						if(d2==6) n66++ ;
						if(d2==7) n76++ ;
					}
					else if(deg==5) {
						if(d2==5) n55++ ;
						if(d2==6) n56++ ;
						if(d2==7) n57++ ;
					}
					else if(deg==7) {
						if(d2==5) n57++ ;
						if(d2==6) n76++ ;
						if(d2==7) n77++ ;
					}

					if(!is_infinite(cir)) {
						double curlen = distance2(to_geex(all_vertices_[i]->point()), to_geex(cir->point())) ;
						if(curlen < minlen) minlen = curlen ;
						if(curlen > maxlen) maxlen = curlen ;
					}

					total ++ ;
					++cir  ;
					
				} while(cir!=adjacent_vertices(all_vertices_[i])) ;
			}
		}
		p66 = (double)n66/total ;
		p76 = (double)n76/total ;
		p56 = (double)n56/total ;
		p55 = (double)n55/total ;
		p57 = (double)n57/total ;
		p77 = (double)n77/total ;

		minlen = sqrt(minlen) ;
		maxlen = sqrt(maxlen) ;
//		maxlen = 2.05*minlen ;
		

		std::cout << "min edge length: " << minlen << std::endl ;
		std::cout << "max edge length: " << maxlen << std::endl << std::endl ;

		bin_size_ = 256 ;
		double step = (maxlen-minlen)/(bin_size_-1) ;
		edge_hist_.resize(bin_size_) ;
		std::fill(edge_hist_.begin(), edge_hist_.end(), 0.0) ;
		for(unsigned int i=0; i<all_vertices_.size(); ++i) {
			if(is_primary(all_vertices_[i])) {
				Vertex_circulator cir = adjacent_vertices(all_vertices_[i]) ;
				do {
					double curlen = distance(to_geex(all_vertices_[i]->point()), to_geex(cir->point())) ;
					if(curlen < maxlen) {  // exclude edges on the vonvex hull
						int bin_id = (int)((curlen-minlen)/step) ;
						edge_hist_[bin_id] ++ ;
					}
					++cir  ;
				} while(cir!=adjacent_vertices(all_vertices_[i])) ;
			}
		}

		for(unsigned int i=0; i<bin_size_; ++i) {
			edge_hist_[i] = edge_hist_[i]/total ;
		}
	}

    // ------------------------------------ Delaunay 

    void Delaunay::clear() {
        baseclass::clear() ;
        all_vertices_.clear() ;
		mirrors_.clear() ;
    }

    Delaunay::Vertex_handle Delaunay::insert(const vec2& p) {
        gx_assert(opened_) ;
        if(Numeric::is_nan(p.x) || Numeric::is_nan(p.y)) {
            std::cerr << "Nan !" << std::endl ;
            return 0 ;
        }
        Vertex_handle result = baseclass::insert(to_cgal(p)) ;
        //result->theta = ::Geex::Numeric::random_float64() * M_PI / 2.0 ;
        //result->energy = 0.0 ;
        //result->rho = 1.0 ; 
        //result->locked = false ;
		result->regularity = 0.0 ;
        return result ;
    }

	Delaunay::Vertex_handle Delaunay::insert(const vec2& p, int index) {
		Vertex_handle v = insert(p) ;
		v->index = index ;
		return v ;
	}

    Delaunay::Vertex_handle Delaunay::nearest(const vec2& p) {
		Delaunay::Vertex_handle result ;
//		result = this->nearest_vertex(to_cgal(p)) ;
		if(dimension() == 1) {
			result = this->nearest_vertex(to_cgal(p)) ;
		}
		else {
			Delaunay::Face_handle f = locate(to_cgal(p)) ;
			double dist = 1e10 ;
			for(unsigned int i=0; i<3; i++) {
				double cur_dist = (to_geex(f->vertex(i)->point()) - p).length2() ;
				if(cur_dist < dist && !is_infinite(f->vertex(i))) {
					dist = cur_dist ;
					result = f->vertex(i) ;
				}
			}
		}
        return result ;
    }

    void Delaunay::remove(const vec2& p) {
        if(number_of_vertices() <= 3) { return ; }
        gx_assert(opened_) ;
        Face_handle f = locate(to_cgal(p)) ;


        double min_d = 1e30 ;
        Vertex_handle v = 0 ;
        for(unsigned int i=0; i<3; i++) {
            if(!is_infinite(f->vertex(i))) {
                double cur_d = (to_geex(f->vertex(i)->point()) - p).length2() ;
                if(cur_d < min_d) {
                    min_d = cur_d ;
                    v = f->vertex(i) ;
                }
            }
        }
        baseclass::remove(v) ;
    }

    void Delaunay::begin_insert() { opened_ = true ; }

    void Delaunay::end_insert(bool redraw) {
        all_vertices_.clear() ;

        for(Vertex_iterator it = vertices_begin() ; it != vertices_end() ; it++) {
            it->dual_intersects_boundary = false ;
			it->dual_infinite = false ;
            it->index = -1 ;
			it->domain = -1 ;
        }

		for(Face_iterator it = faces_begin() ; it != faces_end() ; it++) {
            if(is_infinite(it)) {
//                it->infinite = true ;
                it->dual = vec2(0.0, 0.0) ;
                it->dual_outside = true ;
				it->vertex(0)->dual_infinite = true ;
				it->vertex(1)->dual_infinite = true ;
				it->vertex(2)->dual_infinite = true ;
            } else {
//                it->infinite = false ;
                it->dual = to_geex(baseclass::dual(it)) ;
                it->dual_outside = !in_boundary(it->dual) ;
            }
            if(it->dual_outside) {
				if(in_boundary(to_geex(it->vertex(0)->point())))
					it->vertex(0)->dual_intersects_boundary = true ;
					it->vertex(1)->dual_intersects_boundary = true ;
					it->vertex(2)->dual_intersects_boundary = true ;
			} 
        }        
       
        all_vertices_.clear() ;
        int cur_index = 0 ;
        for(Vertex_iterator it = vertices_begin(); it != vertices_end() ; it++) {
            all_vertices_.push_back(it) ;
            it->index = cur_index ;
			it->domain = 0 ; 
            cur_index++ ;
        }

		/*
		** for the case that Voronoi cell intersects the domain but
		** the Voronoi vertices are outside
		*/
		for(unsigned int i=0; i<boundary_.size(); ++i) {
			PolygonEdge& e = boundary_[i] ;
			// this returns the infinite vertex in release mode
			Vertex_handle vh = nearest(e.vertex[1]) ;
			all_vertices_[vh->index]->dual_intersects_boundary = true ;
		}

		opened_ = false ;
        if(redraw) {
            glut_viewer_redraw() ;            
        }
    }
	
    void Delaunay::insert_random_vertex() {
        double x_min, y_min, z_min, x_max, y_max, z_max ;
        get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
        Geex::vec2 p = random_v() ;
        p.x = x_min + (x_max - x_min) * p.x ;
        p.y = y_min + (y_max - y_min) * p.y ;
        int nb_tries = 0 ;
        while(!in_boundary(p)) {
            nb_tries++ ;
            if(!non_convex_mode_ && nb_tries > 1000) {
                std::cerr << "Could not insert point, probably missing +non_convex flag" << std::endl ;
                exit(-1) ;
            }
            p = random_v() ;
            p.x = x_min + (x_max - x_min) * p.x ;
            p.y = y_min + (y_max - y_min) * p.y ;
        }
        insert(p) ;
    }

    void Delaunay::insert_random_vertices(int nb) {
        begin_insert() ;
        for(unsigned int i=0; i<nb; i++) {
            insert_random_vertex() ;
            std::cerr << (i+1) << '/' << nb << std::endl ;
        }
        if(insert_boundary_ ) {
            for(unsigned int i=0; i<boundary_.size(); i++) {
                Vertex_handle v = insert(boundary_[i].vertex[0]) ;
        //        v->locked = true ;
            }
        }
        end_insert(false) ;
    }

	void Delaunay::insert_grid(int nb) {
        double x_min, y_min, z_min, x_max, y_max, z_max ;
        get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
		int nrow = sqrt(nb) ;
		double dx = (x_max-x_min)/nrow ;
		double dy = (y_max-y_min)/nrow ;

        begin_insert() ;
        for(unsigned int i=0; i<nrow; i++) {
			for(unsigned int j=0; j<nrow; j++) {

				vec2 p ;
				p = vec2(x_min, y_min) + vec2((i+0.5)*dx, (j+0.5)*dy) ;

				if(!in_boundary(p)) {
					if(p[0]<0) {
						while(p[0]<x_min)
							p[0]+= (x_max - x_min) ;
					}
					if(p[0]>x_max)
						while(p[0]>x_max)
							p[0]-=(x_max-x_min) ;
					if(p[1]<y_min) {
						while(p[1]<y_min)
							p[1]+=(y_max-y_min) ;
					}
					if(p[1]>y_max)
						while(p[1]>y_max)
							p[1]-=(y_max-y_min) ;
				}

				insert(p) ;
			}        
        }

		nrow = sqrt(nb - nrow*nrow) ;
		dx = (x_max-x_min)/nrow ;
		dy = (y_max-y_min)/nrow ;

        for(unsigned int i=0; i<nrow; i++) {
			for(unsigned int j=0; j<nrow; j++) {

				vec2 p ;
				p = vec2(x_min, y_min) + vec2((i+0.5)*dx, (j+0.5)*dy) ;

				if(!in_boundary(p)) {
					if(p[0]<0) {
						while(p[0]<x_min)
							p[0]+= (x_max - x_min) ;
					}
					if(p[0]>x_max)
						while(p[0]>x_max)
							p[0]-=(x_max-x_min) ;
					if(p[1]<y_min) {
						while(p[1]<y_min)
							p[1]+=(y_max-y_min) ;
					}
					if(p[1]>y_max)
						while(p[1]>y_max)
							p[1]-=(y_max-y_min) ;
				}

				insert(p) ;
			}        
        }
		end_insert(false) ;

		std::cout << "nb grid points: " << nb_primary() << std::endl ;
	}

    Polygon2* Delaunay::dual_convex_clip(Vertex_handle v, bool close) {
        Polygon2* from = &boundary_ ;
        Polygon2* to = &pong_ ;

        std::vector<Edge> edges ;
        Edge_circulator it = incident_edges(v) ;
        do {
			if(!is_infinite(it))
                edges.push_back(*it) ;
            it++ ;
        } while(it != incident_edges(v)) ;

        // No need to have special cases for infinite vertices: they
        // do not yield any clipping plane !
        for(unsigned int i=0; i<edges.size(); i++) {
            Edge e = edges[i] ;
			int ne = edges.size() ;
            if(is_infinite(e)) { continue ; }
            Geex::Line<real> L = get_dual_line(e) ;
            int E = e.first->vertex(ccw(e.second))->index ;
			Vertex_handle v1 = e.first->vertex(ccw(e.second)) ;
			if(dimension()==1) { // degenerate cases
				int id0 = e.first->vertex(0)->index ;
				int id1 = e.first->vertex(1)->index ;
				if(id0 != v->index) {
					L.inverse() ;
				}
				E = e.first->vertex(0)->index==v->index ? e.first->vertex(1)->index : e.first->vertex(0)->index ; 
			}
			gx_assert(E != v->index || (E==v->index && v1->domain!=v->domain)) ;
			from->convex_clip(*to, L, E, close) ;
			if(from == &boundary_) {
				from = &pong_ ;
				to = &ping_ ;
			} else {
				gx_swap(from, to) ;
			}
        }

		//if(from->size()==0) {
		//	std::cerr << "vertex " << v->index << " has non-clipped cell" << std::endl ;
		//}
        return from ;
		//return dual_convex_clip(boundary_, v, close) ;
    }

	Polygon2* Delaunay::dual_convex_clip(Polygon2& poly, Vertex_handle v, bool close) {
		Polygon2* from = &poly ;
		Polygon2* to = &pong_ ;

		std::vector<Edge> edges ;
		Edge_circulator it = incident_edges(v) ;
		do {
			edges.push_back(*it) ;
			it++ ;
		} while(it != incident_edges(v)) ;

		// No need to have special cases for infinite vertices: they
		// do not yield any clipping plane !
		for(unsigned int i=0; i<edges.size(); i++) {
			Edge e = edges[i] ;
			if(is_infinite(e)) { continue ; }
			Geex::Line<real> L = get_dual_line(e) ;
			int E = e.first->vertex(ccw(e.second))->index ;
			gx_assert(E != v->index) ;
			from->convex_clip(*to, L, E, close) ;
			if(from == &poly) {
				from = &pong_ ;
				to = &ping_ ;
			} else {
				gx_swap(from, to) ;
			}
		}
		return from ;
	}

    int Delaunay::dual_facet_degree(Vertex_handle v, bool period) {
        int result = 0 ;
        if(dual_cell_intersects_boundary(v) && !period) { 
            Polygon2* P = dual_convex_clip(v) ;
            result = P->size() ;
        } else {
            Face_circulator it = incident_faces(v) ;
			if(dimension()>1) {
				do {
					result++ ;
					it++ ;
				} while(it != incident_faces(v)) ;
			}
        } 
        return result ;
    }

	int Delaunay::dual_facet_degree_period(Vertex_handle v) {
		int pid = v->index ; //is_primary(v) ? v->index : v->domain ;
		Vertex_handle pv = all_vertices_[pid] ;
		std::set<Vertex_handle>& mset = mirrors_[pid] ;
		int result = 0 ;

		if(dual_cell_intersects_boundary(pv)) {
			Polygon2* P = dual_convex_clip(pv, true) ;
			for(unsigned int i=0; i<P->size(); ++i) {
				PolygonEdge& e = (*P)[i] ;
				if(e.vertex[0].bisectors.size()==2) {
					result ++ ;
				}
			}

			for(std::set<Vertex_handle>::iterator vi=mset.begin(); vi!=mset.end(); ++vi) {
				P = dual_convex_clip(*vi, true) ;
				for(unsigned int i=0; i<P->size(); ++i) {
					PolygonEdge& e = (*P)[i] ;
					if(e.vertex[0].bisectors.size()==2) {
						result ++ ;
					}
				}
			}
		}
		else {
			Face_circulator it = incident_faces(v) ;
			if(dimension()>1) {
				do {
					result++ ;
					it++ ;
				} while(it != incident_faces(v)) ;
			}
		}

		return result ;
	}

	// -------------------------------------- Periodic Delaunay: necessary copies -------------------------
	void Delaunay::insert_copies(bool full_copy, bool redraw) {
		
		/*clear_copies() ;

		
		if(full_copy) {
			insert_copies_full() ;
		} else {
//			compute_pvd() ;
			insert_copies_ring() ;
		}*/
	}

	void Delaunay::insert_copies_full(bool redraw) {
		int nv = nb_vertices() ;  // master vertices
//		int nb_copies = 8 ;
//		double shift[8][2] = {{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1} } ;

		begin_insert_copies() ;
		mirrors_.clear() ;
		for(int i=0; i<nv; ++i) {
			Vertex_handle mc = all_vertices_[i] ;
			for(int j=1; j<9; ++j) {
				vec2 newp = to_geex(mc->point())+domain_offset_[j] ; //vec2(shift[j][0], shift[j][1]) ;
				Vertex_handle sc = insert(newp) ;
				sc->domain = j ; //mc->index ; 
				sc->index = mc->index ;
				mirrors_[mc->index].insert(sc) ;
			}
		}
		end_insert_copies(redraw) ;
	}

	// the opposite domain the current domain
	// domain 0 maps to itself
	static int opposite_domain(int idx) {
		static int opp_domain[] = {0, 5, 6, 7, 8, 1, 2, 3, 4} ;
		gx_assert(idx>=0 && idx<=8) ;
		return opp_domain[idx] ;
	}

	void Delaunay::insert_copies_ring(bool redraw) {
		static int edge_domain[] = {6, 8, 2, 4} ;
		static int vert_domain[] = {5, 7, 1, 3} ;
		std::vector<std::pair<int, int> > to_insert ;
		unsigned int nv = all_vertices_.size() ;

		mirror_vertices_.clear() ;
		mirrors_.clear() ;

		for(unsigned int i=0; i<nv; ++i) {
			if(dual_cell_intersects_boundary(all_vertices_[i])) {
				Polygon2* P = dual_convex_clip(all_vertices_[i]) ;
				for(unsigned int j=0; j<P->size(); ++j) {
					int bidx = (*P)[j].boundary_index() ;

					if(bidx>=0) {
						to_insert.push_back(std::pair<int, int>(i, edge_domain[bidx])) ;
						if((*P)[j].vertex[0].on_boundary) {
							to_insert.push_back(std::pair<int, int>(i, vert_domain[bidx])) ;
						}
					}
				}
			}
		}

		if(to_insert.size()>0) {
			begin_insert_copies() ;
			for(unsigned int i=0; i<to_insert.size(); ++i) {
				vec2 newp = to_geex(all_vertices_[to_insert[i].first]->point()) + domain_offset_[to_insert[i].second] ;
				Vertex_handle newv = insert(newp) ;
				newv->index = to_insert[i].first ; 	 
				newv->domain = to_insert[i].second ; //to_insert[i].first ; 	 
				gx_assert(!is_primary(newv)) ;
				mirror_vertices_[newv->index].push_back(newv) ;
				mirrors_[to_insert[i].first].insert(newv) ; //
			}
			end_insert_copies(true) ;
		}

		// check copies
		to_insert.clear() ;

		for(unsigned int i=0; i<nv; ++i) {
			std::vector<Vertex_handle>& M = mirror_vertices_[i] ;
			for(unsigned int j=0; j<M.size(); ++j) {
				Vertex_handle vh = M[j] ;
				int opp = opposite_domain(M[j]->domain) ;
				int vid = M[j]->index ;
				gx_assert(vid==i) ;

				Vertex_circulator cir = adjacent_vertices(M[j]) ;
				Vertex_handle vpri = all_vertices_[vid] ;
				do {
					if(is_primary(cir)) {
						Vertex_circulator cir2 = adjacent_vertices(vpri) ;
						bool find = false ;
						do {
							if(cir2->index == cir->index) {
								find = true ;
							}
							++cir2 ;
						} while(cir2!=adjacent_vertices(vpri)) ;

						if(!find) {
							to_insert.push_back(std::pair<int, int>(cir->index, opp)) ;
						}
					}
					++cir ;
				} while(cir!=adjacent_vertices(M[j])) ;
			}
		}

		if(to_insert.size()>0) {
			begin_insert_copies() ;
			for(unsigned int i=0; i<to_insert.size(); ++i) {
				vec2 newp = to_geex(all_vertices_[to_insert[i].first]->point()) + domain_offset_[to_insert[i].second] ;
				Vertex_handle newv = insert(newp) ;
				newv->index = to_insert[i].first ; 	 
				newv->domain = to_insert[i].second ; //to_insert[i].first ; 	 
				gx_assert(!is_primary(newv)) ;
				mirror_vertices_[newv->index].push_back(newv) ;
				mirrors_[to_insert[i].first].insert(newv) ; //
			}
			end_insert_copies(true) ;
		}
	}

	vec2 Delaunay::copy_point(vec2& p, PolygonEdge& e) {
		// find the transform vector. only for square
		vec2 dir = e.vertex[1] - e.vertex[0] ;
		dir /= dir.length() ;
		return p + vec2(-dir[1], dir[0]) ;
	}

	vec2 Delaunay::copy_point(vec2& p, PolygonVertex& e) {
		return p ;	
	}

	void Delaunay::clear_copies(bool redraw) {
		std::vector<vec2> points ;
		for(int i=0; i<all_vertices_.size(); ++i) {
			if(is_primary(all_vertices_[i]))
				points.push_back(to_geex(all_vertices_[i]->point())) ;
		}
		mirrors_.clear() ;
		mirror_vertices_.clear() ;

		clear() ;
		begin_insert() ;
		for(int i=0; i<points.size(); ++i)
			insert(points[i]) ;
		end_insert(redraw) ;
	}

	void Delaunay::begin_insert_copies() {	opened_ = true ; }

	void Delaunay::end_insert_copies(bool redraw) {
		for(Vertex_iterator it = vertices_begin() ; it != vertices_end() ; it++) {
			it->dual_intersects_boundary = false ;
			it->dual_infinite = false ;
		//	it->index = -1 ;
		}

		for(Face_iterator it = faces_begin() ; it != faces_end() ; it++) {
			if(is_infinite(it)) {
//				it->infinite = true ;
				it->dual = vec2(0.0, 0.0) ;
				it->dual_outside = true ;
				it->vertex(0)->dual_infinite = true ;
				it->vertex(1)->dual_infinite = true ;
				it->vertex(2)->dual_infinite = true ;
			} else {
//				it->infinite = false ;
				it->dual = to_geex(baseclass::dual(it)) ;
				it->dual_outside = !in_boundary(it->dual) ;
			}
			if(it->dual_outside) {
				for(int i=0; i<3; ++i) {
					if(in_boundary(to_geex(it->vertex(i)->point())))
						it->vertex(i)->dual_intersects_boundary = true ;
				}
			} else{
				for(int i=0; i<3; ++i) {
					if(!in_boundary(to_geex(it->vertex(i)->point())))
						it->vertex(i)->dual_intersects_boundary = true ;
				}
			}
		}        



		/*
		** for the case that Voronoi cell intersects the domain but
		** the Voronoi vertices are outside
		*/
		for(unsigned int i=0; i<boundary_.size(); ++i) {
			PolygonEdge& e = boundary_[i] ;
			Vertex_handle vh = nearest(e.vertex[0]) ;
			all_vertices_[vh->index]->dual_intersects_boundary = true ;
		}

		opened_ = false ;
		if(redraw) {
			glut_viewer_redraw() ;            
		}
	}

	vec2 Delaunay::mirror_vertex_point(vec2& p, int vidx) {
		gx_assert(vidx>=0 && vidx<4) ;
		return p-v_offset_[vidx] ;
	}

	vec2 Delaunay::mirror_edge_point(vec2& p, int eidx) {
		gx_assert(eidx>=0 && eidx<4) ;
		return p-e_offset_[eidx] ;	
	}

	int  Delaunay::domain_idx(vec2& p) {
		if(in_boundary(p)) 
			return 0 ;
		for(int i=0; i<4; ++i) {
			if(in_boundary(p-v_offset_[i]))
				return i*2+1 ;
			if(in_boundary(p-e_offset_[i]))
				return (i+1)*2 ;
		}

		return -1 ;
	}	

	vec2 Delaunay::translate(int domain_idx, vec2 p) {
		if(domain_idx%2==0) 
			return p-e_offset_[(domain_idx/2-1)];
		else
			return p-v_offset_[((domain_idx-1)/2)];
	}

	int  Delaunay::nb_primary() {
		int nb = 0 ;
		FOR_EACH_VERTEX_DT(Delaunay, this, it) {
			if(is_primary(it))
				nb ++ ;
		}
		return nb ;
	}
	bool Delaunay::neighbor_to_primary(Vertex_handle v) {
		Vertex_circulator cir = adjacent_vertices(v) ;
		do {
			if(is_primary(cir))
				return true ;
			++ cir ;
		} while(cir!=adjacent_vertices(v)) ;
		return false ;
	}

	void Delaunay::compute_pvd() {
		std::map<int, std::set<int> > vvmap ;
		bool suc = false ;
		mirrors_.clear() ;
		mirror_vertices_.clear() ;

		while(!suc) {
			std::vector<std::pair<int,int> > to_insert ;
			
			compute_rvd(false) ;

			// add mirror for domain vertices
			for(unsigned int i=0; i<boundary_.size(); ++i) {
				vec2 bv = boundary_[i].vertex[0] ;
				Vertex_handle v = baseclass::nearest_vertex(to_cgal(bv)) ;
				vec2 temp_v = to_geex(v->point());
				vec2 mp = mirror_vertex_point(temp_v, i) ;
				int  pid = v->index ;//is_primary(v) ? v->index : v->domain ;  // index of the master vertex of the nearest vertex
				int  vid = domain_idx(mp) ; // index of the salve vertex
				if(vid==0) continue ; // itself
				if(vvmap[pid].insert(vid).second) {
					to_insert.push_back(std::pair<int,int>(pid, vid)) ;
				}
			}

			// add mirror for boundary sites
			for(unsigned int i=0; i<rvd_.size(); ++i) {
				if(rvd_[i].size()>0) { // clipped by boundary
					Vertex_handle v = all_vertices_[i] ;
					for(unsigned int j=0; j<rvd_[i].size(); ++j) {
						vec2 temp_v = to_geex(v->point());
						vec2 mp = mirror_edge_point(temp_v, rvd_[i][j].boundary_index()) ;
						int  pid = v->index ;// is_primary(v) ? v->index : v->domain ;  // index of the master vertex of the nearest vertex
						int  vid = domain_idx(mp) ; // index of the salve vertex
						if(vid==0) continue ; // itself
						if(vvmap[pid].insert(vid).second) {
							to_insert.push_back(std::pair<int,int>(pid, vid)) ;
						}
					}
				}
			}

			if(to_insert.size()>0) {
				begin_insert_copies() ;
				for(unsigned int i=0; i<to_insert.size(); ++i) {
					vec2 newp = to_geex(all_vertices_[to_insert[i].first]->point()) + domain_offset_[to_insert[i].second] ;
					Vertex_handle newv = insert(newp) ;
					newv->index  = to_insert[i].first ; 	 
					newv->domain = to_insert[i].second ; 	 
					gx_assert(!is_primary(newv)) ;
					mirror_vertices_[newv->index].push_back(newv) ;
					mirrors_[to_insert[i].first].insert(newv) ; //
				}
				end_insert_copies(true) ;
			}
			else {
				suc = true ;
			}
		} 
	}

	// new implementation
	void Delaunay::compute_pvd2() {
		std::map<int, std::set<int> > vvmap ;
		bool suc = false ;
		mirrors_.clear() ;

		while(!suc) {
			std::vector<std::pair<int,int> > to_insert ;

			for(unsigned int i=0; i<all_vertices_.size(); ++i) {
				Vertex_handle v = all_vertices_[i] ;
				int pid = v->index ; // is_primary(v) ? v->index : v->domain ;  // index of the master vertex of the nearest vertex
				int vid ;
				if(dual_cell_intersects_boundary(v) ) {
					Polygon2* P = dual_convex_clip(v, false) ;
					for(unsigned int j=0; j<P->size(); ++j) {
						PolygonEdge& e = (*P)[j] ;

						for(unsigned int k=0; k<2; ++k) {
							vec2 newp = to_geex(v->point()) ;
							if(e.vertex[k].boundary_edges.size()==1) {
								newp = newp - e_offset_[e.boundary_index()] ;
							}
							if(e.vertex[k].boundary_edges.size()==2) {
								std::set<int>& edges = e.vertex[k].boundary_edges ;
								for(std::set<int>::iterator it=edges.begin(); it!=edges.end(); ++it) {
									newp = newp - e_offset_[*it] ;
								}
							}
							vid = domain_idx(newp) ;
							if(vid==0) continue ; // primary site itself
							if(vvmap[pid].insert(vid).second) {
								to_insert.push_back(std::pair<int,int>(pid, vid)) ;
							}
						}
					}
				}
			}

			if(to_insert.size()>0) {
				begin_insert_copies() ;
				for(unsigned int i=0; i<to_insert.size(); ++i) {
					vec2 newp = to_geex(all_vertices_[to_insert[i].first]->point()) + domain_offset_[to_insert[i].second] ; 
					Vertex_handle newv = insert(newp) ;//to_insert[i].second) ;
					newv->index  = to_insert[i].first ; 	 
					newv->domain = to_insert[i].second ; 	 
					gx_assert(!is_primary(newv)) ;
					mirrors_[to_insert[i].first].insert(newv) ; //
				}
				end_insert_copies(false) ;
			}
			else {
				suc = true ;
			}
		} 
	}

	bool Delaunay::is_full_hex() {
		for(unsigned int i=0; i<all_vertices_.size(); ++i) {
			Vertex_handle v = all_vertices_[i] ;
			if(is_primary(v)) {
				if(dual_facet_degree_period(v)!=6)
					return false ;
			}
		}
		return true ;
	}

	void Delaunay::compute_rvd(bool closed) {
		std::set<int> boundary_cells ;
		std::vector<bool> visited ;
		std::queue<int> Q ;

		rvd_.clear() ;
		rvd_.resize(all_vertices_.size()) ;

		visited.resize(boundary_.size()) ;
		std::fill(visited.begin(), visited.end(), false) ;
		
		for(unsigned int j=0; j<boundary_.size(); ++j) {
			// start from the 1st unvisited boundary edge
			if(visited[j]) continue ;

			PolygonEdge e = boundary_[j] ;
			vec2 mp = 0.5*(e.vertex[0] + e.vertex[1]) ;
			Vertex_handle v = baseclass::nearest_vertex(to_cgal(mp)) ;
			boundary_cells.insert(v->index) ;
			Q.push(v->index) ;

			while(!Q.empty()) {
				int index = Q.front() ;
				Q.pop() ;
				Polygon2 *to ;
				to = dual_convex_clip(all_vertices_[index], closed) ;
				rvd_[index] = *to ;

				// find neighboring boundary cells
				for(int i=0; i<to->size(); ++i) {
					PolygonEdge e = (*to)[i] ;
					if(e.on_boundary()) {
						if(e.vertex[0].bisectors.size()==1) {
							int k = *(e.vertex[0].bisectors.begin()) ;
							if(boundary_cells.insert(k).second) {
								Q.push(k) ;
							}
						}
						if(e.vertex[1].bisectors.size()==1) {
							int k = *(e.vertex[1].bisectors.begin()) ;
							if(boundary_cells.insert(k).second) {
								Q.push(k) ;
							}
						}
						visited[e.boundary_index()] = true ;
					}
				}
			}
		}
	}

	void Delaunay::compute_edge_length() {
		for(unsigned int i=0; i<all_vertices_.size(); ++i) {
			Vertex_handle v = all_vertices_[i] ;
			Vertex_circulator cir = adjacent_vertices(v) ;
			if(is_primary(v)) {
				std::cout << "vertex degree: " << dual_facet_degree_period(v) << std::endl ;
				do {
					double len = (to_geex(v->point()) - to_geex(cir->point())).length() ;
					std::cout << "edge length: " << len << std::endl ;
					cir ++ ;
				} while( cir != adjacent_vertices(v) ) ;
			}
		}

	}

	bool Delaunay::pick_vertex(vec2& pt) {
        Face_handle f = locate(to_cgal(pt)) ;

        double min_d = 1e30 ;
        Vertex_handle v = 0 ;
		bool picked = false ;
        for(unsigned int i=0; i<3; i++) {
            if(!is_infinite(f->vertex(i))) {
                double cur_d = (to_geex(f->vertex(i)->point()) - pt).length2() ;
                if(cur_d < min_d && cur_d < 1e-4) {
                    min_d = cur_d ;
                    v = f->vertex(i) ;
					picked = true ;
                }
            }
        }

		if(picked) {
			vertices_sel_.push_back(v) ;
		}
		return picked ;
	}

	static inline double distance_point_segment(vec2& v0, vec2& v1, vec2& p) {
		if(dot(p-v0, v1-v0) > 0) {
			return det(normalize(v1-v0), p-v0) ;
		} else {
			return gx_min(distance(p, v0), distance(p, v1)) ;
		}
	}

	bool Delaunay::pick_edge(vec2& pt) {
       Face_handle f = locate(to_cgal(pt)) ;
	   bool picked = false ;

        double min_d = 1e30 ;
        Edge e(static_cast<Face_handle>(0), 0) ;
        for(unsigned int i=0; i<3; i++) {
            if(!is_infinite(Edge(f,i))) {
		vec2 temp_v1 = to_geex(f->vertex(f->ccw(i))->point());
		vec2 temp_v2 = to_geex(f->vertex(f->cw(i))->point());
                double cur_d = distance_point_segment(
					temp_v1, 
					temp_v2, 
					pt ) ;
				gx_assert(cur_d>=0) ;
                if(cur_d < min_d && cur_d < 1e-2) {
                    min_d = cur_d ;
                    e = Edge(f, i) ;
					picked = true ;
                }
            }
        }
		if(picked) {
			edges_sel_.push_back(e) ;
		}

		return picked ;
	}
}