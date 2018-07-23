#ifndef _EDITOR_H_ 
#define _EDITOR_H_ 

#include <Geex/combinatorics/map_editor.h>

namespace Geex {

class Delaunay ;
class Map ;

enum IrregEditMode {EDGE_FLIP, EDGE_COLLAPSE, VERTEX_SPLIT, V4_SPLIT, V4_CREATE, VL7_SPLIT, VL7_CREATE} ;

class IrregEditor : public MapEditor{
public:
	IrregEditor(Delaunay *delaunay=nil, Map *map=nil) ;
	void delaunay_to_map() ;
	void map_to_delaunay() ;
	void pick_vertex(vec2& pt) ;
	void pick_edge(vec2& pt) ;
	void clear_picked() ;
	void do_editing() ;
	void smooth(int iter=5) ;
	void perturb(double esp) ;
	void split_v4(Vertex* v, Halfedge* h) ;
	void create_v4(Halfedge* h1, Halfedge* h2) ;
	std::vector<Vertex*>&   picked_verts() { return picked_verts_ ; }
	std::vector<Halfedge*>& picked_edges() { return picked_edges_ ; }
	Facet *cur_facet() { return cur_f_ ; }
	IrregEditMode& edit_mode() { return edit_mode_ ; }

	// generate special input
	void create_grid_map(double radius) ;
	void create_regular_map(double radius) ;

protected:
	void smooth_vertices(std::vector<Vertex*>& to_smth, int iter=5) ;

	Facet* locate(vec2& pt) ;
	bool inside_facet(Facet* f, vec2& pt) ;
	void split_vertex(Vertex* v, Halfedge* f1, Halfedge* g1) ;
       
	inline vec3 local_plane_coords(
            const vec3& p, const vec3& v0, 
            const vec3& v1, const vec3& plane_n,
            const vec3& plane_origin
        ) ;
	 bool is_flippable(Map::Halfedge* h) ;
	 bool edge_flip(Map::Halfedge* h) ;

private:
	Map      *map_ ;
	Delaunay *delaunay_ ;
	IrregEditMode mode_ ;
	std::vector<Vertex*>   picked_verts_ ;
	std::vector<Halfedge*> picked_edges_ ;
	Facet    *cur_f_ ;
	IrregEditMode edit_mode_ ;
} ; 

} // end of namespace Geex		

#endif 