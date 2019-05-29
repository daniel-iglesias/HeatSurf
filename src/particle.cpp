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
#include "particle.h"

Particle::Particle()
{
}

Particle::Particle ( double& x_in, double& xdiv_in,
                     double& y_in, double& ydiv_in,
                     double& z_in, double& zdiv_in,
//                   double& phase_in, double& time_in,
                     double& energy_in//, double& loss_in
                   )
        : x ( x_in ), xdiv ( xdiv_in )
        , y ( y_in ), ydiv ( ydiv_in )
        , z ( z_in ), zdiv ( zdiv_in )
//  , phase(phase_in), time(time_in)
        , energy ( energy_in ) //, loss(loss_in)
{
}


Particle::~Particle()
{
}


