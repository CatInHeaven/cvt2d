#include "irreg_editor.h"
#include "delaunay.h"

#include <Geex/combinatorics/map.h>
#include <Geex/combinatorics/map_editor.h>
#include <Geex/combinatorics/map_builder.h>

namespace Geex {

	IrregEditor::IrregEditor(Delaunay *delaunay, Map *map) 
	: map_(map)
	, delaunay_(delaunay) {
		delaunay_to_map() ;
		set_target(map_) ;
		cur_f_ = nil ;
		edit_mode_ = EDGE_FLIP ;
	}

	void IrregEditor::delaunay_to_map() {
		std::vector<Delaunay::Vertex_handle> vertices = delaunay_->vertices() ;
		std::map<int, std::set<Delaunay::Vertex_handle> >& mirrors = delaunay_->mirrors() ;
		std::map<Delaunay::Vertex_handle, int> indices ;

		map_->clear() ;
		MapBuilder builder(map_) ;
		builder.begin_surface() ;
		
		int cur = 0 ;
		for(unsigned int i=0; i<vertices.size(); ++i) {
			vec2 v = to_geex(vertices[i]->point()) ;
			builder.add_vertex(vec3(v.x, v.y, 0)) ;
			indices[vertices[i]] = cur ;
			cur ++ ;
			
			std::set<Delaunay::Vertex_handle>& vmirror = mirrors[i] ;
			if(vmirror.size() > 0) {
				for(std::set<Delaunay::Vertex_handle>::iterator it=vmirror.begin(); it!=vmirror.end(); ++it) {
					vec2 v = to_geex((*it)->point()) ;
					builder.add_vertex(vec3(v.x, v.y, 0)) ;
					indices[*it] = cur ;
					cur ++ ;			
				}
			}
		}

		FOR_EACH_FINITE_FACE_DT(Delaunay, delaunay_, it) {
			builder.begin_facet() ;
			for(int i=0; i<3; ++i) {
				builder.add_vertex_to_facet(indices[it->vertex(i)]) ;
			}
			builder.end_facet() ;
		}

		builder.end_surface() ;	
	}

	void IrregEditor::map_to_delaunay() {
	}

	void IrregEditor::clear_picked() {
		picked_verts_.clear() ;
		picked_edges_.clear() ;
		cur_f_ = nil ;
	}

	void IrregEditor::do_editing() {
		switch(edit_mode_) {
		case EDGE_FLIP:
			edge_flip(picked_edges_[0]) ;
			break ;
		case EDGE_COLLAPSE: {
			Halfedge *h = picked_edges_[0] ;
			Vertex *v = h->vertex() ;
			vec3 p = 0.5*(v->point()+h->opposite()->vertex()->point()) ;
			collapse_edge(picked_edges_[0]) ;
			v->set_point(p) ;
			break ;
		}
		case VERTEX_SPLIT: {
			Vertex *v = picked_verts_[0] ;
			Halfedge *f1 = picked_edges_[0] ;
			Halfedge *g1 = picked_edges_[1] ;
			if(f1->vertex()==v) f1 = f1->opposite() ;
			if(g1->vertex()==v) g1 = g1->opposite() ;
			split_vertex(v, f1, g1) ;
			break ;
		}			
		case V4_SPLIT:
			split_v4(picked_verts_[0], picked_edges_[0]) ;
			break ;
		case V4_CREATE:
			create_v4(picked_edges_[0], picked_edges_[1]) ;
			break ;
		case VL7_SPLIT:
			break ;
		case VL7_CREATE:
			break ;
		default:
			std::cerr << "the editing mode is not implemented..." << std::endl ;
			break ;
		}

		clear_picked() ;		
	}

	void IrregEditor::split_v4(Vertex* v, Halfedge* h) {
		Halfedge* ih = h->vertex()==v ? h : h->opposite() ;
		Halfedge* flip = ih->prev() ;
		collapse_edge(ih->next()) ;
		edge_flip(flip->prev()) ;
	}
       
	void IrregEditor::create_v4(Halfedge* h1, Halfedge* h2) {
		gx_assert(h1->prev()->opposite()==h2->prev()) ;
		split_edge(h1->prev(), true) ;
	}

	bool IrregEditor::edge_flip(Map::Halfedge* h)	{
		if(!is_flippable(h))
			return false;

		Map::Halfedge* hopp = h->opposite();

		Map::Halfedge* h00 = h->prev();
		Map::Halfedge* h01 = h->next();

		Map::Halfedge* h10 = hopp->next();
		Map::Halfedge* h11 = hopp->prev();

		Map::Facet* f0 = h->facet();
		Map::Facet* f1 = hopp->facet();


		link(h, h11, 1);
		link(h11, h01, 1);
		link(h01, h, 1);

		link(hopp, h00, 1);
		link(h00, h10, 1);
		link(h10, hopp, 1);

		set_facet_on_orbit(h, f0);
		make_facet_key(h);


		set_facet_on_orbit(hopp, f1);
		make_facet_key(hopp);

		make_vertex_key(h, h10->vertex());
		make_vertex_key(hopp, h01->vertex());

		make_vertex_key(h00);
		make_vertex_key(h10);
		make_vertex_key(h01);
		make_vertex_key(h11);

		return true;
	}
        
    bool IrregEditor::is_flippable(Map::Halfedge* h)	{
		// the two edges involved in the flip

		vec3 plane_v0 = normalize(h->vertex()->point() - h->opposite()->vertex()->point());
		vec3 plane_v1 = normalize(h->opposite()->next()->vertex()->point() - h->next()->vertex()->point());		

		// the plane defined by the two edges

		vec3 plane_n = normalize(cross(plane_v1, plane_v0));


		// orthogonalize in-plane vectors

		plane_v0 = normalize(plane_v0 - dot(plane_v0, plane_v1) * plane_v1);

		vec3 plane_origin = h->next()->vertex()->point();


		// 2d coordinates in plane

		vec3 local_t = local_plane_coords(h->vertex()->point(), plane_v0, plane_v1, plane_n, plane_origin);
		vec3 local_b = local_plane_coords(h->opposite()->vertex()->point(), plane_v0, plane_v1, plane_n, plane_origin);
		vec3 local_l = local_plane_coords(h->next()->vertex()->point(), plane_v0, plane_v1, plane_n, plane_origin);
		vec3 local_r = local_plane_coords(h->opposite()->next()->vertex()->point(), plane_v0, plane_v1, plane_n, plane_origin);


		// check if edge intersections lies inside triangles pair (in plane)

		vec3 tb = local_t - local_b;
		vec3 lr = local_l - local_r;

		double ntb[2];
		double ctb;

		double nlr[2];
		double clr;

		ntb[0] = - tb[1];
		ntb[1] = tb[0];

		ctb = -(ntb[0] * local_t[0] + ntb[1] * local_t[1]);


		nlr[0] = - lr[1];
		nlr[1] = lr[0];

		clr = -(nlr[0] * local_l[0] + nlr[1] * local_l[1]);


		double det = ntb[0] * nlr[1] - nlr[0] * ntb[1];

		vec3 intersection(- (nlr[1] * ctb - ntb[1] * clr) / det, 
			- (-nlr[0] * ctb + ntb[0] * clr) / det, 
			0.0);


		double l0 = dot(intersection - local_r, lr) / dot(lr, lr);
		double l1 = dot(intersection - local_b, tb) / dot(tb, tb);

		return l0 > 0.0 && l0 < 1.0 && l1 > 0.0 && l1 < 1.0;
	}

    vec3 IrregEditor::local_plane_coords(
        const vec3& p, const vec3& v0, const vec3& v1, const vec3& plane_n, const vec3& plane_origin
    ) {
        vec3 q(p + dot(plane_origin - p, plane_n) * plane_n);
        return vec3( dot(q - plane_origin, v0), dot(q - plane_origin, v1), 0.0);
    }

	void IrregEditor::split_vertex(Vertex* v, Halfedge* f1, Halfedge* g1) {
		Vertex *v2 = new_vertex() ;
		Halfedge *f2 = new_edge() ;
		Halfedge *g2 = new_edge() ;
		Halfedge *f =  new_edge() ;
		Halfedge *g =  f->opposite() ;
		std::vector<Halfedge*> fan ;

		gx_assert(f1->opposite()->vertex()==v && g1->opposite()->vertex()==v) ;

		Halfedge *cir = g1->opposite() ;
		do {
			fan.push_back(cir) ;
			cir = cir->opposite()->prev() ;
		} while (cir != f1->opposite()) ;


		link(f2, f1->opposite()->next(), 1) ;
		link(f1->opposite()->prev(), f2, 1) ;
		link(f2->opposite(), f1->opposite(), 1) ;
		link(f1->opposite(), f, 1) ;
		link(f, f2->opposite(), 1) ;
		set_halfedge_vertex(f2->opposite(), f1->vertex()) ;
		set_halfedge_vertex(f2, v2) ;
		set_halfedge_vertex(f, v2) ;
		gx_assert(f1->opposite()->vertex()==v) ;
		set_facet_on_orbit(f, new_facet()) ;
		make_facet_key(f, f->facet()) ;
		make_vertex_key(f) ;
		make_vertex_key(f2->opposite()) ;

		set_facet_on_orbit(f2, f2->next()->facet()); 
		make_facet_key(f2, f2->facet()) ;

		link(g2, g1->opposite()->next(), 1) ;
		link(g1->opposite()->prev(), g2, 1) ;
		link(g1->opposite(), g, 1) ;
		link(g, g2->opposite(), 1) ;
		link(g2->opposite(), g1->opposite(), 1) ;
		set_halfedge_vertex(g, v) ;
		set_halfedge_vertex(g2, v) ;
		set_halfedge_vertex(g2->opposite(), g1->vertex()) ;
		set_halfedge_vertex(g1->opposite(), v2) ;
		set_facet_on_orbit(g, new_facet()) ;
		make_facet_key(g, g->facet()) ;
		make_vertex_key(g) ;

		set_facet_on_orbit(g2, g2->next()->facet()); 
		make_facet_key(g2, g2->facet()) ;

		for(int i=0; i<fan.size(); ++i) {
			set_halfedge_vertex(fan[i], v2) ;
		}

		// smooth
		v2->set_point(v->point()) ;
		std::vector<Vertex*> to_smth ;
		to_smth.push_back(v) ;
		to_smth.push_back(v2) ;
		smooth_vertices(to_smth, 5) ;
	}

	void IrregEditor::smooth(int iter) {
		std::vector<vec3> new_points ;

		for(int i=0; i<iter; ++i) {
			int cur = 0 ;

			new_points.assign(map_->nb_vertices(), vec3(0, 0, 0)) ;
			FOR_EACH_VERTEX(Map, map_, v) {
				Halfedge *h = v->halfedge() ;
				if(!v->is_on_border()) {
					do {
						new_points[cur] += h->opposite()->vertex()->point() ;
						h = h->next_around_vertex() ;
					} while(h!=v->halfedge()) ;
					new_points[cur] = 1.0/v->degree()*new_points[cur] ;
				}
				else {
					new_points[cur] = v->point() ;
				}
				cur ++ ;
			}

			cur = 0 ;
			FOR_EACH_VERTEX(Map, map_, v) {
				v->set_point(new_points[cur]) ;
				cur ++ ;
			}
		}

	}

	void IrregEditor::smooth_vertices(std::vector<Vertex*>& to_smth, int iter) {
		std::vector<vec3> new_points ;
		new_points.resize(to_smth.size()) ;
		for(int i=0; i<iter; ++i) {
			new_points.assign(new_points.size(), vec3(0, 0, 0)) ;
			
			for(unsigned int j=0; j<to_smth.size(); ++j) {
				Halfedge *h = to_smth[j]->halfedge() ;
				do {
					new_points[j] += h->opposite()->vertex()->point() ;
					h = h->next_around_vertex() ;
				} while(h!=to_smth[j]->halfedge()) ;
				new_points[j] = 1.0/to_smth[j]->degree()*new_points[j] ;
			}

			for(unsigned int j=0; j<to_smth.size(); ++j) {
				to_smth[j]->set_point(new_points[j]) ;
			}
		}
	}

	bool IrregEditor::inside_facet(Facet* f, vec2& pt) {
		Halfedge* h = f->halfedge() ;
		vec3 v0 = h->vertex()->point() ;
		vec3 v1 = h->next()->vertex()->point() ;
		vec3 v2 = h->prev()->vertex()->point() ;
		vec3 p(pt.x, pt.y, 0) ;

		vec3 c01 = cross(v0-p, v1-p) ;
		vec3 c12 = cross(v1-p, v2-p) ;
		vec3 c20 = cross(v2-p, v0-p) ;

		if(dot(c01, c12)>0 && dot(c01, c20)>0) {
			return true ;
		}
		return false ;
	}

	Map::Facet* IrregEditor::locate(vec2& pt) {
		FOR_EACH_FACET(Map, map_, it) {
			if(inside_facet(it, pt)) {
				return it ;
			}
		}
		return nil ;
	}

	void IrregEditor::pick_vertex(vec2& pt) {
		Facet* f = locate(pt) ;
		cur_f_ = f ;
		vec3 p3d(pt.x, pt.y, 0) ;
		double min_dist = 1e10 ;
		Vertex* v = nil ;
		if(f!=nil) {
			Halfedge* h = f->halfedge() ;
			do {
				double cur_dist = distance(h->vertex()->point(), p3d) ;
				if(cur_dist < min_dist && cur_dist < 1e-2) {
					min_dist = cur_dist ;
					v = h->vertex() ;
				}
				h = h->next() ;
			} while (h!=f->halfedge()) ;
		}

		if(v!=nil && std::find(picked_verts_.begin(), picked_verts_.end(), v)==picked_verts_.end()) {
			picked_verts_.push_back(v) ;
		}
	}

	static inline double distance_point_segment(vec3& v0, vec3& v1, vec3& p) {
		if(dot(p-v0, v1-v0) > 0) {
			return cross(normalize(v1-v0), p-v0).length() ;
		} else {
			return gx_min(distance(p, v0), distance(p, v1)) ;
		}
	}

	void IrregEditor::pick_edge(vec2& pt) {
		Facet* f = locate(pt) ;
		cur_f_ = f ;
		if(f==nil) return ;
		vec3 p3d(pt.x, pt.y, 0) ;
		double min_dist = 1e10 ;
		Halfedge *e = nil ;
		
		Halfedge* h = f->halfedge() ;
		do {
			double cur_dist = distance_point_segment(h->vertex()->point(), h->opposite()->vertex()->point(), p3d) ;
			if(cur_dist < min_dist) {
				min_dist = cur_dist ;
				e = h ;
			}
			h = h->next() ;
		} while (h!=f->halfedge()) ;

		if(e!=nil) {
			bool flag1 = std::find(picked_edges_.begin(), picked_edges_.end(), e)==picked_edges_.end() ;
			bool flag2 = std::find(picked_edges_.begin(), picked_edges_.end(), e->opposite())==picked_edges_.end() ;
			if(flag1 && flag2) {
				picked_edges_.push_back(e) ;
			}
		}
	}

	static inline int vindex_grid(int i, int j, int n, int margin) {
		return (n+2*margin)*(i+margin) + j+margin ;
	}

	void IrregEditor::create_grid_map(double radius) {
		int    n = ceil(1.0/radius) ;
		double len = 1.0/n ;
		int    margin = 3 ;
	
		map_->clear() ;
		MapBuilder builder(map_) ;
		builder.begin_surface() ;

		for(int i=-margin; i<n+margin; ++i) {
			for(int j=-margin; j<n+margin; ++j) {
				builder.add_vertex(vec3((i+0.5)*len, (j+0.5)*len, 0)) ;
			}
		}

		for(int i=-margin; i<n+margin-1; ++i) {
			for(int j=-margin; j<n+margin-1; ++j) {
				builder.begin_facet() ;
				builder.add_vertex_to_facet(vindex_grid(i  , j  , n, margin)) ;
				builder.add_vertex_to_facet(vindex_grid(i+1, j  , n, margin)) ;
				builder.add_vertex_to_facet(vindex_grid(i+1, j+1, n, margin)) ;
				builder.end_facet() ;

				builder.begin_facet() ;
				builder.add_vertex_to_facet(vindex_grid(i  , j  , n, margin)) ;
				builder.add_vertex_to_facet(vindex_grid(i+1, j+1, n, margin)) ;
				builder.add_vertex_to_facet(vindex_grid(i  , j+1, n, margin)) ;
				builder.end_facet() ;
			}
		}

		builder.end_surface() ;	
	}

	static inline int vindex_regular(int i, int j, int n, int m, int margin) {
		return (n+2*margin)*(i+margin) + j+margin ;
	}

	void IrregEditor::create_regular_map(double radius) {
		int n = 1.0 / radius ;
		int m = int(2.0/sqrt(3.0)*n+0.5) ;
		double newr = sqrt(2.0/(sqrt(3.0)*m*n)) ;
		double width = n*newr ;
		double height = m*sqrt(3.0)/2*newr ;
		int margin = 3 ;

		map_->clear() ;
		MapBuilder builder(map_) ;
		builder.begin_surface() ;

		for(int i=-margin; i<m+margin; ++i) {
			for(int j=-margin; j<n+margin; ++j) {
				if(i%2==0) {
					builder.add_vertex(vec3(j*newr, i*sqrt(3.0)/2*newr, 0)) ;
				} 
				else {
					builder.add_vertex(vec3((j+0.5)*newr, i*sqrt(3.0)/2*newr, 0)) ;
				}				
			}
		}

		for(int i=-margin; i<m+margin-1; ++i) {
			for(int j=-margin; j<n+margin-1; ++j) {
				builder.begin_facet() ;
				if(i%2==0) {
					builder.add_vertex_to_facet(vindex_regular(i  , j  , n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i  , j+1, n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j  , n, m, margin)) ;
				}else {
					builder.add_vertex_to_facet(vindex_regular(i  , j  , n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i  , j+1, n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j+1, n, m, margin)) ;
				}
				builder.end_facet() ;

				builder.begin_facet() ;
				if(i%2==0) {
					builder.add_vertex_to_facet(vindex_regular(i  , j+1, n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j+1, n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j  , n, m, margin)) ;
				}else {
					builder.add_vertex_to_facet(vindex_regular(i  , j  , n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j+1, n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j  , n, m, margin)) ;
				}
				builder.end_facet() ;
			}
		}

		builder.end_surface() ;	
	}
	
	void IrregEditor::perturb(double esp) {
		FOR_EACH_VERTEX(Map, map_, it) {
			vec3 dir = normalize(vec3(Numeric::random_float64(), Numeric::random_float64(), 0)) ;
			it->set_point(it->point() + esp*dir) ;
		}
	}
} // end of namespace Geex		