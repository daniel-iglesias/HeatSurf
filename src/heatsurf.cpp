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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <unistd.h>

#include <iostream>
#include <cstdlib>
#include <LMX/lmx_base_stopwatch.h>

#include "simulation.h"

using namespace std;

int main ( int argc, char *argv[] )
{
    if ( argc != 2 )
    {
        cout<< "::: ERROR :::\n Need one paramater with name of input file"
            << std::endl;
    }
    else
    {
        char temp[100];
        cout<< "EXECUTING 'beam_dump "<< argv[1]<< "'\n";
        cout<< "Working directory: "<< getcwd ( temp,100 ) <<"\n";
        Simulation theSimulation;
        {
            lmx::ExactStopwatch sw;
            theSimulation.read ( argv[1] );
        }
        {
            lmx::ExactStopwatch sw;
            theSimulation.compute();
        }
        theSimulation.output();
    }

    cout << "Done!" << endl;

    return EXIT_SUCCESS;
}
