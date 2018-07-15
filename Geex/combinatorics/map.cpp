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
 
 
 
 
 
 
 

#include <Geex/combinatorics/map.h>
#include <algorithm>
#include <stack>

namespace Geex {

//_________________________________________________________

    MapCombelObserver<Map::Vertex>::MapCombelObserver(Map* m) : map_(m) {
        map_->add_vertex_observer(this);
    }
    
    MapCombelObserver<Map::Vertex>::~MapCombelObserver() {
        map_->remove_vertex_observer(this);
    }

    MapCombelObserver<Map::Halfedge>::MapCombelObserver(Map* m) : map_(m) {
        map_->add_halfedge_observer(this);
    }
    
    MapCombelObserver<Map::Halfedge>::~MapCombelObserver() {
        map_->remove_halfedge_observer(this);
    }

    MapCombelObserver<Map::Facet>::MapCombelObserver(Map* m) : map_(m) {
        map_->add_facet_observer(this);
    }
    
    MapCombelObserver<Map::Facet>::~MapCombelObserver() {
        map_->remove_facet_observer(this);
    }

//_________________________________________________________

    Map::Map() { 
    }
    
    Map::~Map() {
        clear() ;
    } 

    // ____________________ Predicates _______________________

    bool Map::is_triangulated() const {
        FOR_EACH_FACET_CONST(Map, this, it) {
            if(!it->is_triangle()) {
                return false ;
            }
        }
        return true ;
    }

    // ____________________ Observers ________________________

    void Map::add_vertex_observer(MapCombelObserver<Vertex>* obs) {
        vertex_observers_.push_back(obs) ;
    }
    
    void Map::remove_vertex_observer(MapCombelObserver<Vertex>* obs) {
        std::vector<MapCombelObserver<Vertex>* >::iterator it = std::find(
            vertex_observers_.begin(), vertex_observers_.end(), obs
        ) ;
        gx_assert(it != vertex_observers_.end()) ;
        vertex_observers_.erase(it) ;
    }
    
    void Map::add_halfedge_observer(MapCombelObserver<Halfedge>* obs) {
        halfedge_observers_.push_back(obs)  ;
    }
    
    void Map::remove_halfedge_observer(MapCombelObserver<Halfedge>* obs) {
        std::vector<MapCombelObserver<Halfedge>* >::iterator it = std::find(
            halfedge_observers_.begin(), halfedge_observers_.end(), obs
        ) ;
        gx_assert(it != halfedge_observers_.end()) ;
        halfedge_observers_.erase(it) ;
    }
    
    void Map::add_facet_observer(MapCombelObserver<Facet>* obs) {
        facet_observers_.push_back(obs) ;
    }
    
    void Map::remove_facet_observer(MapCombelObserver<Facet>* obs) {
        std::vector<MapCombelObserver<Facet>* >::iterator it = std::find(
            facet_observers_.begin(), facet_observers_.end(), obs
        ) ;
        gx_assert(it != facet_observers_.end()) ;
        facet_observers_.erase(it) ;
    }


    void Map::notify_add_vertex(Vertex* v) {
        for(
            std::vector<MapCombelObserver<Vertex>* >::iterator
                it=vertex_observers_.begin(); it!=vertex_observers_.end(); it++
        ) {
            (*it)->add(v) ;
        }
    }
    
    void Map::notify_remove_vertex(Vertex* v) {
        for(
            std::vector<MapCombelObserver<Vertex>* >::iterator
                it=vertex_observers_.begin(); it!=vertex_observers_.end(); it++
        ) {
            (*it)->remove(v) ;
        }
    }
    
    void Map::notify_add_halfedge(Halfedge* h) {
        for(
            std::vector<MapCombelObserver<Halfedge>* >::iterator
                it=halfedge_observers_.begin(); 
            it!=halfedge_observers_.end(); it++
        ) {
            (*it)->add(h) ;
        }
    }
    
    void Map::notify_remove_halfedge(Halfedge* h) {
        for(
            std::vector<MapCombelObserver<Halfedge>* >::iterator
                it=halfedge_observers_.begin(); 
            it!=halfedge_observers_.end(); it++
        ) {
            (*it)->remove(h) ;
        }
    }
    
    void Map::notify_add_facet(Facet* f) {
        for(
            std::vector<MapCombelObserver<Facet>* >::iterator
                it=facet_observers_.begin(); 
            it!=facet_observers_.end(); it++
        ) {
            (*it)->add(f) ;
        }
    }
    
    void Map::notify_remove_facet(Facet* f) {
        for(
            std::vector<MapCombelObserver<Facet>* >::iterator
                it=facet_observers_.begin(); 
            it!=facet_observers_.end(); it++
        ) {
            (*it)->remove(f) ;
        }
    }

    // ____________________ Modification _____________________

    void Map::clear() {
        vertices_.clear() ;
        halfedges_.clear() ;
        facets_.clear() ;
        
        vertex_property_manager_.clear() ;
        halfedge_property_manager_.clear() ;
        facet_property_manager_.clear() ;
    }

    void Map::erase_all() {
        clear() ;
    }

    void Map::clear_inactive_items() {
        // TODO: traverse the inactive items list, 
        //  and remove the properties ...
        vertices_.clear_inactive_items() ;
        halfedges_.clear_inactive_items() ;
        facets_.clear_inactive_items() ;
    }

    // _____________________ Low level ______________________

    
    Map::Halfedge* Map::new_edge() {
        Halfedge* h1 = new_halfedge() ;
        Halfedge* h2 = new_halfedge() ;
        h1-> set_opposite(h2) ;
        h2-> set_opposite(h1) ;
        h1-> set_next(h2) ;
        h2-> set_next(h1) ;
        h1-> set_prev(h2) ;
        h2-> set_prev(h1) ;
        return h1 ;
    }

    void Map::delete_edge(Halfedge* h) {
        delete_halfedge(h-> opposite()) ;
        delete_halfedge(h) ;
    }

    Map::Vertex* Map::new_vertex() {
        Vertex* result = vertices_.create() ;
        vertex_property_manager_.new_record(result) ;
        notify_add_vertex(result) ;
        return result ;
    }
    
    Map::Vertex* Map::new_vertex(const Map::Vertex* rhs) {
        Vertex* result = vertices_.create() ;
        result->set_point(rhs->point()) ;
        vertex_property_manager_.new_record(result, rhs) ;
        notify_add_vertex(result) ;
        return result ;
    }

    Map::Halfedge* Map::new_halfedge() {
        Halfedge* result = halfedges_.create() ;
        halfedge_property_manager_.new_record(result) ;
        notify_add_halfedge(result) ;
        return result ;
    }

    Map::Halfedge* Map::new_halfedge(const Map::Halfedge* rhs) {
        Halfedge* result = halfedges_.create() ;
        halfedge_property_manager_.new_record(result, rhs) ;
        notify_add_halfedge(result) ;
        return result ;
    }
    
    Map::Facet* Map::new_facet() {
        Facet* result = facets_.create() ;
        facet_property_manager_.new_record(result) ;
        notify_add_facet(result) ;
        return result ;
    }

    Map::Facet* Map::new_facet(const Map::Facet* rhs) {
        Facet* result = facets_.create() ;
        facet_property_manager_.new_record(result, rhs) ;
        notify_add_facet(result) ;
        return result ;
    }

	Map::TexVertex* Map::new_tex_vertex() {
		// For the moment, TexVertex does Attribute management
		// (it keeps a pointer to the AttributeManager)
		TexVertex* result = new TexVertex(tex_vertex_property_manager()) ;
		return result ;
	}

	Map::TexVertex* Map::new_tex_vertex(const Map::TexVertex* rhs) {
		// For the moment, TexVertex does Attribute management
		// (it keeps a pointer to the AttributeManager)
		TexVertex* result = new TexVertex(tex_vertex_property_manager(),rhs) ;
		return result ;
	}

    void Map::delete_vertex(Vertex* v) {
        notify_remove_vertex(v) ;
        vertex_property_manager_.delete_record(v) ;
        vertices_.destroy(v) ;
    }
    
    void Map::delete_halfedge(Halfedge* h) {
        notify_remove_halfedge(h) ;
        halfedge_property_manager_.delete_record(h) ;
        halfedges_.destroy(h) ;
    }
    
    void Map::delete_facet(Facet* f) {
        notify_remove_facet(f) ;
        facet_property_manager_.delete_record(f) ;
        facets_.destroy(f) ;
    }
	
	void Map::delete_tex_vertex(TexVertex* tv) {
		// For the moment, TexVertex does Attribute management
		// (it keeps a pointer to the AttributeManager)
		delete tv ;
	}

    void Map::activate_vertex(Vertex* v) {
        notify_add_vertex(v) ;
        vertices_.activate(v) ;
    }
    
    void Map::activate_halfedge(Halfedge* h) {
        notify_add_halfedge(h) ;
        halfedges_.activate(h) ;
    }
    
    void Map::activate_facet(Facet* f) {
        notify_add_facet(f) ;
        facets_.activate(f) ;
    }
    
    void Map::deactivate_vertex(Vertex* v) {
        notify_remove_vertex(v) ;
        vertices_.deactivate(v) ;
    }
    
    void Map::deactivate_halfedge(Halfedge* h) {
        notify_remove_halfedge(h) ;
        halfedges_.deactivate(h) ;
    }
    
    void Map::deactivate_facet(Facet* f) {
        notify_remove_facet(f) ;
        facets_.deactivate(f) ;
    }

    void Map::get_connected_component(
        Vertex* h, std::vector<Vertex*>& l
    ) {

        Property<Vertex, bool> visited(vertex_property_manager()) ;

        std::stack<Vertex*> stack ;
        stack.push(h) ;

        while(!stack.empty()) {
            Vertex* top = stack.top() ;
            stack.pop() ;
            if(!visited[top]) {
                visited[top] = true ;
                l.push_back(top) ;
                Halfedge* cir = top->halfedge() ;
                do {
                    Vertex* cur = cir-> opposite()-> vertex() ;
                    if(!visited[cur]) {
                        stack.push(cur) ;
                    }
                    cir = cir->next_around_vertex() ;
                } while(cir != top->halfedge()) ;
            }
        }
    }

    void Map::assert_is_valid() const {
        FOR_EACH_VERTEX_CONST(Map, this, it) {
            it->assert_is_valid() ;
        }
        FOR_EACH_HALFEDGE_CONST(Map, this, it) {
            it->assert_is_valid() ;
        }
        FOR_EACH_FACET_CONST(Map, this, it) {
            it->assert_is_valid() ;
        }
    }


//_________________________________________________________

    MapMutator::~MapMutator() { target_ = nil ;  }

    void MapMutator::set_target(Map* m) { target_ = m ; }

    void MapMutator::link(Halfedge* h1, Halfedge* h2, int dim) {
        switch(dim) {
        case 1:
            h1-> set_next(h2) ;
            h2-> set_prev(h1) ;
            break ;
        case 2:
            h1-> set_opposite(h2) ;
            h2-> set_opposite(h1) ;
            break ;
        }
    }

    Map::Halfedge* MapMutator::new_edge(Map::Halfedge* rhs) {
        Halfedge* h1 = new_halfedge(rhs) ;
        Halfedge* h2 = new_halfedge(rhs->opposite()) ;
        h1-> set_opposite(h2) ;
        h2-> set_opposite(h1) ;
        h1-> set_next(h2) ;
        h2-> set_next(h1) ;
        h1-> set_prev(h2) ;
        h2-> set_prev(h1) ;
        return h1 ;
    }
        
    void MapMutator::set_vertex_on_orbit(Halfedge* h, Vertex* v) {
        Halfedge* it = h ;
        do {
            it->set_vertex(v) ;
            it = it->next_around_vertex() ;
        } while(it != h) ;
    }

    void MapMutator::set_facet_on_orbit(Halfedge* h, Facet* f) {
        Halfedge* it = h ;
        do {
            it-> set_facet(f) ;
            it = it->next() ;
        } while(it != h) ;
    }

//_________________________________________________________

}

