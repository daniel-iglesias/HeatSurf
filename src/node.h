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
#ifndef NODE_H
#define NODE_H

#include <iostream>

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file node.h

  \brief Node class for graphical representation purposes.

  Simple point definition and manipulation. It has a scalar property for storing the power density.

  \author Daniel Iglesias <daniel.iglesias@ciemat.es>

 */
//////////////////////////////////////////// Doxygen file documentation (end)

/**
 @author Daniel Iglesias <daniel.iglesias@ciemat.es>
 */
class Node
{
public:
    Node();

    Node ( double, double, double );

    ~Node();

    void setScalar ( double scalar_in )
    { scalar = scalar_in; }

    void addScalar ( double scalar_in )
    { scalar += scalar_in; }

    void addDensity3D ( double dens_in )
    {
        density = ( ( n_densities3D-1 ) *density + dens_in ) /n_densities3D;
        ++n_densities3D;
    }

    double getScalar() { return scalar; }

    void addDensity ( double density_in )
    {
//       density += density_in;
        density = ( ( n_densities-1 ) *density + density_in ) /n_densities;
        ++n_densities;
    }

    double getDensity()
    {
//       if(density > 1){
//         std::cout << "density = " << density
//         << ", n_densities3D = " << n_densities3D << std::endl;
//       }

        return density;
    }

    double getX() { return x; }
    double getY() { return y; }
    double getZ() { return z; }

    double getShiftX() { return shift_x; }
    double getShiftY() { return shift_y; }
    double getShiftZ() { return shift_z; }

    void setX(double x_in) { x=x_in; }
    void setY(double y_in) { y=y_in; }
    void setZ(double z_in) { z=z_in; }

    void applyShift (double, double, double);

private:
    double x, y, z, scalar, density;
    double shift_x, shift_y, shift_z;
    int n_densities3D, n_densities;
};

#endif
