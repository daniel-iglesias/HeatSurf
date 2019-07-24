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
#include "simulation.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
extern "C"
{
#include <gsl/gsl_rng.h>
#include "gsl/gsl_randist.h"
}

Simulation::Simulation()
        : particle_shifting_x ( 0. )
        , particle_shifting_y ( 0. )
        , beam_opening_x ( 1. )
        , beam_opening_y ( 1. )
        , beam_steering_x ( 0. )
        , beam_steering_y ( 0. )
        , backscattering_section_distance ( 0. )
{
}


Simulation::~Simulation()
{
}

void Simulation::read ( char * file_name )
{
    std::ifstream input ( file_name );
    std::string keyword;

    while ( input >> keyword )
    {
        if ( keyword == "GEOMETRY" ) readGeometry ( input );
        else if ( keyword == "PARTICLES" ) readParticles ( input );
//     else if(keyword == "THEORETICAL") readBeamParameters(input);
        else if ( keyword == "PARTICLEGENERATOR" ) createParticles ( input );
        else if ( keyword == "SHIFTPARTICLES" )
            input >> particle_shifting_x >> particle_shifting_y;
        else if ( keyword == "STEERBEAM" )
        {
            input >> beam_steering_x >> beam_steering_y;
            double distance_from_quadrupole = 2500.;
            particle_shifting_x = distance_from_quadrupole * beam_steering_x*1E-3;
            particle_shifting_y = distance_from_quadrupole * beam_steering_y*1E-3;
        }
        else if ( keyword == "OPENBEAM" ) input >> beam_opening_x >> beam_opening_y;
        else if ( keyword == "BACKSCATTERING" )
            input >> backscattering_section_distance;
    }
}

void Simulation::readGeometry ( std::ifstream & input )
{
    std::string keyword;

    input >> keyword;
    if ( keyword == "CONE" )
    {
        double z0, length, initDiam;
        int sectors, section_distance;

        input >> z0 >> length >> initDiam >> sectors >> section_distance; // [mm]
        cout << "NEW CONE: "
		<< z0 << ", "
		<< length << ", "
		<< initDiam << ", "
		<< sectors << ", "
		<< section_distance << endl; // [mm]

        geometries.push_back ( new Cone ( "CONE", z0, length, initDiam, sectors ) );
        geometries.back()->setSections ( section_distance );
    }
    if ( keyword == "CONE2" )
    {
        double z0, length, initDiam, finalDiam;
        int sectors, section_distance;

        input >> z0 >> length >> initDiam >> finalDiam >> sectors >> section_distance; // [mm]
        cout << "NEW CONE2: "
		<< z0 << ", "
		<< length << ", "
		<< initDiam << ", "
		<< finalDiam << ", "
		<< sectors << ", "
		<< section_distance << endl; // [mm]

        geometries.push_back ( new Cone2 ( "CONE2", z0, length, initDiam, finalDiam, sectors ) );
        geometries.back()->setSections ( section_distance );
    }

    //********************************************************************************************8
    //Dan input here

        if ( keyword == "CONE3" )
    {
        double z0, length, initDiam, finalDiam,initialAngle, finalAngle;
        int sectors, section_distance;

        input >> z0 >> length >> initDiam >> finalDiam >> initialAngle >> finalAngle >>  sectors >> section_distance; // [mm]
        cout << "NEW CONE3: "
		<< z0 << ", "
		<< length << ", "
		<< initDiam << ", "
		<< finalDiam << ", "
		<< initialAngle << ", "
		<< finalAngle << ", "
		<< sectors << ", "
		<< section_distance << endl; // [mm]

        geometries.push_back ( new Cone3 ( "CONE3", z0, length, initDiam, finalDiam, initialAngle, finalAngle,  sectors ) );
        geometries.back()->setSections ( section_distance );
    }

    else if(keyword == "RING"){
        double z0, length, internalDiam, externalDiam;
        int sectors, section_distance;

        input >> z0 >> externalDiam >> internalDiam
              >> sectors >> section_distance; // [mm]

        geometries.push_back
            ( new Ring("RING", z0, externalDiam, internalDiam, sectors ) );
        geometries.back()->setSections(section_distance);
      }


    else if ( keyword == "CYLINDER" )
    {
        double z0, length, initDiam;
        int sectors, section_distance;

        input >> z0 >> length >> initDiam >> sectors >> section_distance; // [mm]
        cout << z0 << ", "
             << length << ", "
             << initDiam << ", "
             << sectors << ", "
             << section_distance <<  endl; // [mm]

        geometries.push_back
        ( new Cylinder ( "CYLINDER", z0, length, initDiam, sectors ) );
        geometries.back()->setSections ( section_distance );
    }

    //***********************************************************************************
    //Dan input here
    else if ( keyword == "NEWCYLINDER" )
    {
        double z0, length, initDiam, initialAngle, finalAngle;
        int sectors, section_distance;

        input >> z0 >> length >> initDiam >> initialAngle >> finalAngle >> sectors >> section_distance; // [mm]
        cout << z0 << ", "
             << length << ", "
             << initDiam << ", "
             << initialAngle << ", "
             << finalAngle << ", "
             << sectors << ", "
             << section_distance <<  endl; // [mm]

        geometries.push_back
        ( new newCylinder ( "NEWCYLINDER", z0, length, initDiam, initialAngle, finalAngle, sectors ) );
        geometries.back()->setSections ( section_distance );
    }

    else if ( keyword == "OGIVE" )
    {
        double z0, length, initDiam, ogive_radius;
        int sectors, section_distance;

        input >> z0 >> length >> initDiam >> sectors >> section_distance
        >> ogive_radius; // [mm]

        geometries.push_back
        ( new Ogive ( "OGIVE", z0, length, initDiam, sectors, ogive_radius ) );
        geometries.back()->setSections ( section_distance );
    }

//  else if(keyword == "REVOLUTE"){
//    double length, initDiam;
//    int sectors;
//
//    input >> length >> initDiam >> sectors; // units [mm]
//
//    geometries.push_back( new Revolute("REVOLUTE", length, initDiam, sectors ) );
//  }


   else if(keyword == "PLATE")
    {
        double z0, length, initDiam, angle;
        int sectors, section_distance;

        input >> z0 >> length >> initDiam >> angle >> sectors >> section_distance; // [mm]
        cout << z0 << ", "
             << length << ", "
             << initDiam << ", "
             << angle << ", "
             << sectors << ", "
             << section_distance <<  endl; // [mm]

        geometries.push_back
        ( new Plate ( "PLATE", z0, length, initDiam, angle, sectors ) );
        geometries.back()->setSections ( section_distance );
    }


   else if(keyword == "2PLATES"){
     double length, initDiam, width, orientation;
     int sectors;

     input >> length >> initDiam >> sectors >> width >> orientation; // units [mm]

     geometries.push_back( new TwoPlates("2PLATES", length, initDiam, sectors, width, orientation ) );
   }
}


void Simulation::readParticles ( std::ifstream & input )
{
    std::string preffix, suffix;
    std::ostringstream file_name;

    // preffix first_num last_num delta_num suffix...
    input >> preffix >> suffix;

//   for( int i=first_num; i<=last_num; i+=delta_num ){

    file_name << preffix << suffix;
    std::cout << file_name.str() << std::endl;

    std::ifstream input_p ( file_name.str().c_str() );
    file_name.str ( "" );

    std::ofstream output_p ( "particles_out.txt" );

    char a;
    do {input_p.get ( a ); output_p.put ( a );}
    while ( a!='\n' ); //thrash first line

    std::string keyword;
    double x, xdiv, y, ydiv, z, zdiv, time, phase, energy, loss;
    while ( input_p >> keyword )
    {
        x = atof ( keyword.c_str() );
        input_p >> xdiv >> y >> ydiv >> z >> zdiv >> time >> phase >> energy >> loss; //units: [mm]
        x *= beam_opening_x;
        y *= beam_opening_y;
        x += particle_shifting_x;
        y += particle_shifting_y;
//        z /*+*/= z;
        xdiv += beam_steering_x;
        ydiv += beam_steering_y;
        xdiv *= beam_opening_x * 1E-3;
        ydiv *= beam_opening_y * 1E-3;
        zdiv *= 1E-3;
        particles.push_back ( new Particle ( x, xdiv,
                                             y, ydiv,
                                             z, zdiv,
//                                          time, phase,
                                             energy/*, loss*/ )
                            );
	output_p << std::setiosflags( std::ios::scientific )
		 << std::setiosflags( std::ios::showpos );
      	output_p << std::setprecision(6);

	output_p << x << '\t' << xdiv*1E3 << '\t' << y << '\t' << ydiv*1E3 << '\t' << z << '\t' << zdiv*1E3 << '\t' << time << '\t' << phase << '\t' << energy << '\t' << loss << "\n";
//       particles[i].push_back( new Particle( x, y, z, energy) );
//     }
    }
    std::cout << "Particle size: " << particles.size() <<   std::endl;
}


void Simulation::createParticles ( std::ifstream & input )
{
    std::ofstream output_p ( "particles.txt" );
    double number_of_particles;

    input >> number_of_particles;
    std::cout << "number_of_particles: " << number_of_particles << endl;

    this->readBeamParameters ( input );
    // read energy, sigma_x, sigma_y, div_x, div_y, corr_x_xdiv, corr_y_ydiv

    double z_i = 0.;
    double s_x, s_y, sdiv_x, sdiv_y, corr_div_x, corr_div_y ;
    double max_x, max_y;
    double x, xdiv, y, ydiv, z ( 0 ), zdiv ( 0 ), time ( 0 ), phase ( 0 ), energy ( 9 ), loss ( 0 );
    int i, j=0;
    const gsl_rng_type * T;
    gsl_rng * r;

    // first line of file:
    output_p << "x(mm), xdiv(mrad), y(mm), ydiv(mrad) "
	     << "z(mm), zdiv(mrad), time(s), phase, energy(eV), loss" << endl;

    /* create a generator chosen by the
    environment variable GSL_RNG_TYPE */
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc ( T );

    energy = theoricParameters[0];
    s_x = theoricParameters[1];
    s_y = theoricParameters[2];
    sdiv_x = theoricParameters[3];
    sdiv_y = theoricParameters[4];
    corr_div_x = theoricParameters[5];
    corr_div_y = theoricParameters[6];
    cout << "corrdiv_x = " << corr_div_x << endl;
    max_x = 6*s_x;
    max_y = 6*s_y;
    for ( i=0; i<number_of_particles; ++i )
    {
        gsl_ran_bivariate_gaussian ( r, s_x, sdiv_x, corr_div_x, &x,  &xdiv );
        gsl_ran_bivariate_gaussian ( r, s_y, sdiv_y, corr_div_y, &y,  &ydiv );
        // trying only with circular beam...
//         r = (min + (max-min) * (double)rand()/RAND_MAX);
//         x = (1./2.*3.1416*std::pow(s_x,2))*std::exp(-std::pow(r,2)/(2.*std::pow(s_x,2)));
//         r = (min + (max-min) * (double)rand()/RAND_MAX);
//         y = (1./2.*3.1416*std::pow(s_y,2))*std::exp(-std::pow(r,2)/(2.*std::pow(s_y,2)));

//     if( x < max_x || y < max_y ){
        output_p << x << " " << xdiv << " "
        << y << " " << ydiv << " "
        << z << " " << zdiv << " "
        << time << " " << phase << " "
        << energy << " " << loss << endl; //units: [mm]
        x += particle_shifting_x;
        y += particle_shifting_y;
        xdiv *= 1E-3;
        ydiv *= 1E-3;
        zdiv *= 1E-3;
        particles.push_back ( new Particle ( x, xdiv,
                                             y, ydiv,
                                             z, zdiv,
//                                          time, phase,
                                             energy/*, loss*/ )
                            );
//     }
//     else // particle too far, drop it!
//       --i;
    }
    gsl_rng_free ( r );
    std::cout << "Particle size: " << particles.size() << std::endl;
}

void Simulation::createDivertedParticles ( std::ifstream & input )
{
    std::ofstream output_p ( "diverted-particles.txt" );
    double number_of_particles;

    input >> number_of_particles;
    std::cout << "number_of_particles: " << number_of_particles << endl;

    this->readBeamParameters ( input );
    // read energy, sigma_x, sigma_y, div_x, div_y, corr_x_xdiv, corr_y_ydiv

    double z_i = 0.;
    double s_x, s_y, sdiv_x, sdiv_y, corr_div_x, corr_div_y ;
    double max_x, max_y;
    double x, xdiv, y, ydiv, z ( 0 ), zdiv ( 0 ), time ( 0 ), phase ( 0 ), energy ( 9 ), loss ( 0 );
    int i, j=0;
    const gsl_rng_type * T;
    gsl_rng * r;

    // first line of file:
    output_p << "x(mm), xdiv(mrad), y(mm), ydiv(mrad) "
	     << "z(mm), zdiv(mrad), time(s), phase, energy(eV), loss" << endl;

    /* create a generator chosen by the
    environment variable GSL_RNG_TYPE */
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc ( T );

    energy = theoricParameters[0];
    s_x = theoricParameters[1];
    s_y = theoricParameters[2];
    sdiv_x = theoricParameters[3];
    sdiv_y = theoricParameters[4];
    corr_div_x = theoricParameters[5];
    corr_div_y = theoricParameters[6];
    cout << "corrdiv_x = " << corr_div_x << endl;
    max_x = 6*s_x;
    max_y = 6*s_y;
    for ( i=0; i<number_of_particles; ++i )
    {
        gsl_ran_bivariate_gaussian ( r, s_x, sdiv_x, corr_div_x, &x,  &xdiv );
        gsl_ran_bivariate_gaussian ( r, s_y, sdiv_y, corr_div_y, &y,  &ydiv );
        // trying only with circular beam...
//         r = (min + (max-min) * (double)rand()/RAND_MAX);
//         x = (1./2.*3.1416*std::pow(s_x,2))*std::exp(-std::pow(r,2)/(2.*std::pow(s_x,2)));
//         r = (min + (max-min) * (double)rand()/RAND_MAX);
//         y = (1./2.*3.1416*std::pow(s_y,2))*std::exp(-std::pow(r,2)/(2.*std::pow(s_y,2)));

//     if( x < max_x || y < max_y ){
        output_p << x << " " << xdiv << " "
        << y << " " << ydiv << " "
        << z << " " << zdiv << " "
        << time << " " << phase << " "
        << energy << " " << loss << endl; //units: [mm]
        x += particle_shifting_x;
        y += particle_shifting_y;
        xdiv *= 1E-3;
        ydiv *= 1E-3;
        zdiv *= 1E-3;
        particles.push_back ( new Particle ( x, xdiv,
                                             y, ydiv,
                                             z, zdiv,
//                                          time, phase,
                                             energy/*, loss*/ )
                            );
//     }
//     else // particle too far, drop it!
//       --i;
    }
    gsl_rng_free ( r );
    std::cout << "Particle size: " << particles.size() << std::endl;
}

void Simulation::readBeamParameters ( std::ifstream & input )
{
//   std::string preffix, suffix;
//   std::ostringstream file_name;

//   input >> preffix >> suffix;
//   file_name << preffix << suffix;
//   std::cout << file_name.str() << std::endl;

//   std::ifstream input_p(file_name.str().c_str());
//   file_name.str("");

//   char a;
//   do{input_p.get(a);} while(a!='\n'); //thrash first line

    double energy, sx0, sy0, div_x, div_y, div_corr_x, div_corr_y;
    input >> energy >>
    sx0 >> sy0 >>
    div_x >> div_y >>
    div_corr_x >> div_corr_y;
    theoricParameters.push_back ( energy );
    theoricParameters.push_back ( sx0 );
    theoricParameters.push_back ( sy0 );
    theoricParameters.push_back ( div_x );
    theoricParameters.push_back ( div_y );
    theoricParameters.push_back ( div_corr_x );
    theoricParameters.push_back ( div_corr_y );
    std::cout
        << energy << " "
        << sx0 << " "
        << sy0 << " "
        << div_x << " "
        << div_y << " "
        << div_corr_x << " "
        << div_corr_y << std::endl;
}

void Simulation::readDivertedParameters ( std::ifstream & input )
{
    double energy, sx0, sy0, div_x, div_y, div_corr_x, div_corr_y;
    input >> energy >>
    sx0 >> sy0 >>
    div_x >> div_y >>
    div_corr_x >> div_corr_y;
    theoricParameters.push_back ( energy );
    theoricParameters.push_back ( sx0 );
    theoricParameters.push_back ( sy0 );
    theoricParameters.push_back ( div_x );
    theoricParameters.push_back ( div_y );
    theoricParameters.push_back ( div_corr_x );
    theoricParameters.push_back ( div_corr_y );
    std::cout
        << energy << " "
        << sx0 << " "
        << sy0 << " "
        << div_x << " "
        << div_y << " "
        << div_corr_x << " "
        << div_corr_y << std::endl;
}


void Simulation::compute()
{
    std::vector< Geometry* >::iterator it_geom;
    std::vector< Particle* >::iterator it_part;
    for ( it_geom = geometries.begin();
            it_geom!= geometries.end();
            ++it_geom
        )
    {
        ( *it_geom )->computeGeometry();

        for ( it_part = particles.begin();
                it_part!= particles.end();
                ++it_part
            )
        {
            ( *it_geom )->computeIntersection ( *it_part );

            if ( particles.size() > 0 )
            {
                ( *it_geom )->computeNodalPower ( *it_part );
            }

//       else{
//         (*it_geom)->computeNodalPower( theoricParameters );
//       }
        }
        ( *it_geom )->computeGrid3D( );
        ( *it_geom )->computePowerDensity ( particles.size() );
    }

}

void Simulation::output()
{
    std::vector<Geometry*>::iterator it_geom;
    for ( it_geom = geometries.begin();
            it_geom!= geometries.end();
            ++it_geom
        )
    {
        ( *it_geom )->outputPowerFile ( particles.size() );

        ( *it_geom )->outputPowerDensityFile();

        theWindow.addGeometry( *it_geom );

//        ( *it_geom )->drawGeometry();

//        ( *it_geom )->drawScalar( true );

        ( *it_geom )->outputTable();

//        ( *it_geom )->drawGrid();

//        ( *it_geom )->drawGridScalar();

        ( *it_geom )->outputAnsys3D();

        ( *it_geom )->outputVTKfiles();

        if ( backscattering_section_distance > 0. )
            ( *it_geom )->outputBackscattering ( backscattering_section_distance );
    }
    theWindow.drawGeometry();
    theWindow.drawScalar( true );
}

