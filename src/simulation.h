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
#ifndef SIMULATION_H
#define SIMULATION_H

#include <fstream>
#include <vector>
#include <map>

#include "cone.h"
#include "cone2.h"
#include "cone3.h"
#include "cylinder.h"
#include "newCylinder.h"
#include "ring.h"
#include "ogive.h"
#include "plate.h"
//#include "revolute.h"
#include "twoplates.h"
#include "particle.h"
#include "window.h"
//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file simulation.h

  \brief Procedures for the program workflow.

  Reads files, creates all of the objects, lauches the computation and
  creates the output.

  \author Daniel Iglesias <daniel.iglesias@ciemat.es>

 */
//////////////////////////////////////////// Doxygen file documentation (end)

/**
 @author Daniel Iglesias <daniel.iglesias@ciemat.es>
 */
class Simulation
{
public:
    Simulation();

    ~Simulation();

    void read ( char* );

    void compute();

    void output();


private:
    void readGeometry ( std::ifstream & );

    void readParticles ( std::ifstream & );

    void createParticles ( std::ifstream & );

    void createDivertedParticles ( std::ifstream & );

    void readBeamParameters ( std::ifstream & );

    void readDivertedParameters ( std::ifstream & );

private:
    std::vector<Geometry*> geometries;
    std::vector< Particle* > particles; //particles[section][i]
    std::vector<double> theoricParameters;
    Window theWindow;
    double particle_shifting_x,
    particle_shifting_y,
    beam_opening_x,
    beam_opening_y,
    beam_steering_x,
    beam_steering_y,
    backscattering_section_distance;
    bool output_backscattering;
};

#endif
