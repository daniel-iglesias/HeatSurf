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
#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file particle.h

  \brief Represents each of the finite charged particles.

  Very basic class with stored the x-y position in the plane section and the energy of the particle.

  \author Daniel Iglesias <daniel.iglesias@ciemat.es>

 */
//////////////////////////////////////////// Doxygen file documentation (end)

/**
 @author Daniel Iglesias <daniel.iglesias@ciemat.es>
 */
class Particle
{
public:
    Particle();

    Particle ( double&, double&,
               double&, double&,
               double&, double&,
//              double&, double&,
               double&/*, double&*/ );

    ~Particle();

    double& getX() { /*std::cout << "X = " << x << std::endl;*/ return x;}
    double& getY() { /*std::cout << "Y = " << y << std::endl;*/ return y;}
    double& getZ() { /*std::cout << "Y = " << y << std::endl;*/ return z;}

    double& getXdiv() { return xdiv;}
    double& getYdiv() { return ydiv;}
    double& getZdiv() { return zdiv;}

    double& getEnergy() {return energy;}


private:
    double x, xdiv, y, ydiv, z, zdiv, /*phase, time, */energy/*, loss*/;

};

#endif
