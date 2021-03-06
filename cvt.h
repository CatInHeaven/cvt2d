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

#ifndef __CVT__
#define __CVT__

#include "delaunay.h"
#include "delaunay_graphics.h"
#include "delaunay_io.h"
#include "delaunay_cvt.h"

#include <Geex/graphics/geexob.h>

namespace Geex {

    class CVT : public Geexob , public Delaunay, public DelaunayCVT, public DelaunayGraphics, public DelaunayIO {
    public:
        CVT() ;

        virtual void get_bbox(
            real& x_min, real& y_min, real& z_min,
            real& x_max, real& y_max, real& z_max
        ) ;

        virtual void pre_draw() ;
        virtual void do_draw() ;
        virtual void set_frame(int x) ;

        bool& show_lp() { return show_lp_ ; }
        int& lp_shader() { return lp_shader_ ; }

    private:
        bool show_lp_ ;
        int lp_shader_ ;
        GLuint colormap_id_ ;
        GLuint center_id_ ;
        GLuint lp_id_ ;
        GLuint U_id_ ;
        GLuint V_id_ ;
        GLuint colormap_texture_ ;
    } ;

}

#endif
