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
#include "element.h"

#include <vtkQuad.h>
#include <vtkHexahedron.h>

Element::Element()
{
}

Element::~Element()
{
}

double Element::getDensity3D()
{
    if ( density3D.size() == 0 )
        return 0.;
    else
    {
        double density_mean ( 0. );
        std::vector<double>::iterator it_density3D = density3D.begin();
        std::vector<double>::iterator it_density3D_end = density3D.end();
//     cout << "density3D=";
        for ( it_density3D; it_density3D != it_density3D_end; ++it_density3D )
        {
//       cout << *it_density3D<< ",";
            density_mean += ( *it_density3D );
        }
        return ( density_mean / density3D.size() );
    }
}


void Element::generateGeometry()
{
    int counter = 0;
    switch ( this->type )
    {
    case 0 : // quad
        // should check if all necessary nodes are defined
        geometry = vtkQuad::New();
        for ( it = connectivity.begin();
                it!= connectivity.end();
                ++it
            )
        {
            geometry->GetPointIds()->SetId ( counter, *it );
            ++counter;
        }
        break;

    case 1 : // quad
        // should check if all necessary nodes are defined
        geometry = vtkHexahedron::New();
        for ( it = connectivity.begin();
                it!= connectivity.end();
                ++it
            )
        {
            geometry->GetPointIds()->SetId ( counter, *it );
            ++counter;
        }
//         cout << "Quad: " << number << endl;

    }
}


