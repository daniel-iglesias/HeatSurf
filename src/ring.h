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
#ifndef RING_H
#define RING_H

#include "geometry.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
 *\file ring.h
 *
 *\brief Anular geometry.
 *
 *\author Daniel Iglesias <daniel.iglesias@ciemat.es>
 *
 */
//////////////////////////////////////////// Doxygen file documentation (end)


/**
 @author Daniel Iglesias <daniel.iglesias@ciemat.es>
 */
class Ring : public Geometry
{
public:
    Ring();

    Ring ( std::string, double, double, double, int );

    ~Ring();

    void setSections ( double );

    void computeGeometry( );

    void residue ( lmx::Vector<double>&, lmx::Vector<double>& );
    void jacobian ( lmx::Matrix<double>&, lmx::Vector<double>& );
    void computeIntersection ( Particle* );

    void computeNodalPower ( Particle* );

    void outputTable( );

    void outputPowerFile ( int );

    void outputPowerDensityFile( );


private:
    double intDiam, extDiam;
    double x,y,z,vx,vy,vz;
};

#endif
