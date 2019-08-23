/***************************************************************************
 *   Copyright (C) 2007 by Daniel Iglesias   *
 *   daniel.iglesias@ciemat.es   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "node.h"

Node::Node()
{
}

Node::Node ( double x_in, double y_in, double z_in )
        : x ( x_in ), y ( y_in ), z ( z_in )
        , scalar ( 0 ), density ( 0 )
        , n_densities3D ( 1 ), n_densities ( 1 )
{
}


Node::~Node()
{
}

void Node::applyShift( double x_in, double y_in, double z_in )
{   shift_x = x_in;
    shift_y = y_in;
    shift_z = z_in;
    x += shift_x;
    y += shift_y;
    z += shift_z;
}
