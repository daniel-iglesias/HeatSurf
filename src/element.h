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
#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>

class vtkCell;

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file element.h

  \brief Facet planar geometry for graphical representation purposes.

  Uses the Node class for the defining the vertices. Usually, the dimensions are proportional to the number of sections (longitudinal) and sectors (transversal divisions).

  \author Daniel Iglesias <daniel.iglesias@ciemat.es>

 */
//////////////////////////////////////////// Doxygen file documentation (end)

/**
 @author Daniel Iglesias <daniel.iglesias@ciemat.es>
 */
class Element
{
public:
    Element();

    Element ( int type_in ) : type ( type_in )
    {}

    ~Element();

    void setNumber ( int number_in )
    { number = number_in; }

    void setArea ( double area_in )
    { area = area_in; }

    void addNode ( int node_number )
    { connectivity.push_back ( node_number ); }

    void addDensity ( double dens_in )
    { density += dens_in; }

    void addDensity3D ( double dens_in )
    { density3D.push_back ( dens_in ); }

    int getType()
    { return type; }

    int getNumber()
    { return number; }

    double getArea()
    { return area; }

    double getDensity()
    { return density; }

    double getDensity3D();

    int getNumberOfNodes()
    {
        if ( this->type == 0 ) return 4; // 2D - 4 nodes
        if ( this->type == 1 ) return 8; // 3D - 8 nodes
        else return 0;
    }

    std::vector<int>& getConnectivity()
    { return connectivity; }

    void generateGeometry();

public:
    vtkCell* geometry;

private:
    int type, number;
    double area, density;
    std::vector<double> density3D;
    std::vector<int> connectivity;
    std::vector<int>::iterator it;
};

#endif
