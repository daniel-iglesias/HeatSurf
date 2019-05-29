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

#include "cone2.h"

// #include <iomanip.h>

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>

Cone2::Cone2()
{
}

Cone2::Cone2 ( std::string type_in,
             double z0_in,
             double length_in,
             double initDiam_in,
             double finalDiam_in,
             int par3_in )
        : Cone ( type_in, z0_in, length_in, initDiam_in, par3_in )
        , finalDiam (finalDiam_in)
{
    slope = atan2 ( (initDiam-finalDiam)/2., length );
    gridWidth = initDiam/2.;
}


Cone2::~Cone2()
{
}

// Trying same as cone really
void Cone2::residue ( lmx::Vector<double>& res, lmx::Vector<double>& conf )
{
    double t = conf.readElement ( 0 );
    res.writeElement (
        ( pow ( x+vx*t,2 ) +pow ( y+vy*t,2 ) ) *pow ( cos ( slope ),2 )
        -pow ( t- ( length+z0 ),2 ) *pow ( sin ( slope ),2 )
        , 0
    );
}

// Trying same as cone really
void Cone2::jacobian ( lmx::Matrix<double>& jac, lmx::Vector<double>& conf )
{
    double t = conf.readElement ( 0 );
    jac.writeElement (
        ( 2*vx* ( x+vx*t ) +2*vy* ( y+vy*t ) ) *pow ( cos ( slope ),2 )
        -2* ( t- ( length+z0 ) ) *pow ( sin ( slope ),2 )
        , 0
        , 0
    );
}

void Cone2::computeIntersection ( Particle* particle )
{
    x = particle->getX();
    y = particle->getY();
    z = particle->getZ();
    vx = particle->getXdiv();
    vy = particle->getYdiv();
    vz = particle->getZdiv();
    lmx::Vector<double> initialGuess ( 1 ); // zero
    lmx::NLSolver<Cone2> theSolver;
    theSolver.setInfo ( 0 );
    theSolver.setInitialConfiguration ( initialGuess );
    theSolver.setSystem ( *this );
    theSolver.setResidue ( &Cone2::residue );
    theSolver.setJacobian ( &Cone2::jacobian );
//   theSolver.setConvergence( &Cone::myConvergence );
//   theSolver.setMaxIterations( 100 );
    theSolver.solve ( 100 );
//   cout << theSolver.getSolution().readElement(0) << endl;
    paramTrajectories.push_back ( theSolver.getSolution().readElement ( 0 ) );

}


