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
#ifndef CONE3_H
#define CONE3_H

#include "geometry.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
 *\file Cone2.h
 *
 *\brief Truncated conical geometry.
 *
 * Geometry with origin in s=0 defined by s-length and initial diameter and final diameter
 * (constant slope).
 *
 *\author Daniel Iglesias <daniel.iglesias@ciemat.es>
 *
 */
//////////////////////////////////////////// Doxygen file documentation (end)


/**
 @author Daniel Iglesias <daniel.iglesias@ciemat.es>
 */
class Cone3 : public Geometry
{
public:
    Cone3();

    Cone3 ( std::string, double, double, double, double, int );

    ~Cone3();

    void setSections ( double );

    void computeGeometry( );

    void residue ( lmx::Vector<double>&, lmx::Vector<double>& );
    void jacobian ( lmx::Matrix<double>&, lmx::Vector<double>& );
    void computeIntersection ( Particle* );

    void computeNodalPower ( Particle* );

    void outputTable( );

    void outputPowerFile ( int );

    void outputPowerDensityFile( );


protected:
    double slope, initDiam, finalDiam, apparentLength;

private:
    double x,y,z,vx,vy,vz;
};

#endif
