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

#include "cone3.h"

// #include <iomanip.h>

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>

Cone3::Cone3()
{
}

Cone3::Cone3 ( std::string type_in,
             std::string name_in,
             double z0_in,
             double length_in,
             double initDiam_in,
             double finalDiam_in,
             double initial_Angle_in,
             double final_Angle_in,
             int sectors_in)
        : Geometry ( type_in, name_in, z0_in, sectors_in, length_in )
        , initDiam ( initDiam_in )
        , finalDiam (finalDiam_in)
        , initialAngle ( initial_Angle_in )
        , finalAngle ( final_Angle_in )
{
    slope = atan2 ( (initDiam-finalDiam)/2., length );
    gridWidth = initDiam/2.;
    apparentLength = (initDiam/2.)/tan(slope);
    initialAngle = initialAngle*(3.1416/180);
    finalAngle = finalAngle*(3.1416/180);
    std::cout << initDiam << endl << finalDiam << endl << slope << endl << length <<endl;
}


Cone3::~Cone3()
{
}


void Cone3::setSections ( double distance )
{
    sections = length /distance + 1 ;
    std::cout << "Sections = " << sections << std::endl;
}


void Cone3::computeGeometry()
{
    int i, j;
    double pi = 3.1416;
    double theta;
    double radius, ra, rb;
    double a, b, c;
    double x, y, z, zlocal;
    double sectionDif = length / ( sections - 1 );
    double segmentSize = finalAngle - initialAngle; //only creating a geometry the size between the two finalAngles
//    std::cout << "Segment size = " << (segmentSize)*(180/pi) << " Degrees" << std::endl;
    // nodes map is going to be filled to a size of i*j where
    // i is the number of longitudinal sections and
    // j is the number of transversal (angular) sectors.
    for ( i=0; i<sections; ++i )
    {
        z = i * sectionDif + z0;
        radius = ( initDiam/2. ) - ( z-z0 ) *tan ( slope );
//     radius = (length + z0 - z) * tan(slope);
        for ( j=0; j<sectors; ++j )   //closed chain
        {
            theta = (segmentSize ) * j / sectors + initialAngle; // finalAngle between -pi and pi
            x = radius * cos ( theta );
            y = radius * sin ( theta );
            nodes[i*sectors + j] = new Node ( x, y, z );
        }
    }

    theta = (segmentSize) / sectors; // finalAngle between -pi and pi
    for ( i=0; i<sections-1; ++i )
    {
        zlocal = i * sectionDif;
        ra = ( initDiam/2. ) - zlocal*tan ( slope );
        rb = ( initDiam/2. ) - ( zlocal+sectionDif ) *tan ( slope );
//     ra = (length - zlocal) * tan(slope);
//     rb = (length - (zlocal+sectionDif)) * tan(slope);
        a = ra * theta;
        b = rb * theta;
        c = sqrt ( pow ( ra-rb, 2 ) + pow ( sectionDif, 2 ) );
        for ( j=0; j<sectors-1; ++j )  //last sector is different
        {
            elements.push_back ( new Element ( 0 ) // type quad
                               );
            elements.back()->setNumber ( i*sectors + j );
            elements.back()->addNode ( i*sectors + j );
            elements.back()->addNode ( i*sectors + j+1 );
            elements.back()->addNode ( ( i+1 ) *sectors + j+1 );
            elements.back()->addNode ( ( i+1 ) *sectors + j );
            elements.back()->generateGeometry();
            elements.back()->setArea ( ( a + b ) /2 * c );
//       cout << "a = "<< a << ", b = "<< b << ", Area = " << elements.back()->getArea() << endl;
        }
        // the last sector:
        elements.push_back ( new Element ( 0 ) // type quad
                           );
        elements.back()->setNumber ( i*sectors + j );
        elements.back()->addNode ( i*sectors + j );
        elements.back()->addNode ( i*sectors );
        elements.back()->addNode ( ( i+1 ) *sectors );
        elements.back()->addNode ( ( i+1 ) *sectors + j );
        elements.back()->generateGeometry();
        elements.back()->setArea ( ( a + b ) /2 * c );
//     cout << "a = "<< a << ", b = "<< b << ", Area = " << elements.back()->getArea() << endl;
    }

    gridPoints = vtkPoints::New();
    gridPoints->SetNumberOfPoints ( nodes.size() );

    for ( std::map<int, Node*>::iterator it = nodes.begin();
            it != nodes.end();
            ++it
        )
    {
        gridPoints->InsertPoint ( it->first, it->second->getX(), it->second->getY(), it->second->getZ() );
    }

    grid = vtkUnstructuredGrid::New();
    grid->Allocate ( 1000,1000 );

    for ( std::vector<Element*>::iterator it=elements.begin();
            it!= elements.end();
            ++it
        )
    {
        grid->InsertNextCell ( ( *it )->geometry->GetCellType(), ( *it )->geometry->GetPointIds() );
    }

    grid->SetPoints ( gridPoints );
}


void Cone3::residue ( lmx::Vector<double>& res, lmx::Vector<double>& conf )
{
    double t = conf.readElement ( 0 );
    double x_p, y_p, z_p;

    x_p = x+vx*t - shift_x * shift_mag; // Particles need to be shifted back to the original position
    y_p = y+vy*t - shift_y * shift_mag;
    z_p = z+t-z0;
    std::vector<double> rotated_coordinates = inverseRotatePointAroundGeometryAxis(x_p, y_p, z_p);

    res.writeElement (
        ( pow ( rotated_coordinates[0],2 ) +pow ( rotated_coordinates[1],2 ) - pow ( initDiam/2.,2 ) )
        , 0
        );

//    res.writeElement (
//        (( pow ( x+vx*t - shift_x * shift_mag,2 ) +pow ( y+vy*t - shift_y * shift_mag,2 ) ) *pow ( cos ( slope ),2 )
//        -pow ( t- ( apparentLength+z0 ),2 )*pow ( sin ( slope ),2 ))
//        , 0
//    );
}

void Cone3::jacobian ( lmx::Matrix<double>& jac, lmx::Vector<double>& conf )
{
    double t = conf.readElement ( 0 );
    jac.writeElement (
        ( 2*vx* ( x+vx*t - shift_x * shift_mag ) +2*vy* ( y+vy*t - shift_y * shift_mag ) ) *pow ( cos ( slope ),2 )
        -2* ( t- ( apparentLength+z0 ) ) *pow ( sin ( slope ),2 )
        , 0
        , 0
    );
}

void Cone3::computeIntersection ( Particle* particle )
{
    x = particle->getX();
    y = particle->getY();
    z = particle->getZ();
    vx = particle->getXdiv();
    vy = particle->getYdiv();
    vz = particle->getZdiv();
    lmx::Vector<double> initialGuess ( 1 ); // zero
    lmx::NLSolver<Cone3> theSolver;
    theSolver.setInfo ( 0 );
    theSolver.setInitialConfiguration ( initialGuess );
    theSolver.setSystem ( *this );
    theSolver.setResidue ( &Cone3::residue );
    theSolver.setJacobian ( &Cone3::jacobian );
//   theSolver.setConvergence( &Cone3::myConvergence );
//   theSolver.setMaxIterations( 100 );
    theSolver.solve ( 100 );
//   cout << theSolver.getSolution().readElement(0) << endl;
    paramTrajectories.push_back ( theSolver.getSolution().readElement ( 0 ) );

}

void Cone3::computeNodalPower ( Particle* particle )
{
        double pi = 3.1416;
    if ( paramTrajectories.back() > (z0 /*- shift_z * shift_mag*/) &&
            paramTrajectories.back() <= ( z0 + length /*+ shift_z * shift_mag*/ ) )
    {

        // First we compute the intersection point of the last particle computed
        // in the "computeIntersection" function:
        x = particle->getX() + particle->getXdiv() * paramTrajectories.back();
        y = particle->getY() + particle->getYdiv() * paramTrajectories.back();
        z = /*particle->getZ() +*/ paramTrajectories.back() - z0;
//     cout << "x = " << x << '\t';
//     cout << "y = " << y << '\t';
//     cout << "z = " << z << endl;
        //   cout << "particle->getZdiv() = " << particle->getZdiv() << endl;
//     cout << "paramTrajectories.back() = " << paramTrajectories.back() << endl;

        // We compute the power (= energy * current) of the particle:
        double particle_power = particle->getEnergy() * 0.1255;
        // Now we add to the adjacent points:
        // The intersection shall be between two sections:
        int section_back = floor ( z / ( length/ ( sections-1 ) ) );
//     cout << "section_back = " << section_back << endl;
        // ... defining a proximity factor:
        double section_back_factor
        = 1. - fmod ( z , length/ ( sections-1 ) ) / ( length/ ( sections-1 ) );
//     cout << "section_back_factor = " << section_back_factor << endl;
        // We search the sector where it lies (between two nodes):
        double theta = pi + (atan2 ( x,y )); //adding pi here to make theta range from 0 to 2pi instead of -pi to pi (matching initalAngle and finalAngle)
        if (theta >=initialAngle && theta <= finalAngle){// The particle collides with the segmented cone. If not, do nothing (ignore collision)
            //cout << theta << endl;
            double pi = 3.1416;
            double segmentSize = finalAngle - initialAngle;
            int sector_back = floor( ( (theta) - initialAngle) / ( finalAngle -initialAngle ) /*+ pi*/ * sectors ); // BUG: not clear that this rotation is associated with any parameter. Should be clearer
            //cout << "sector_back = " << sector_back << endl;
            // ... defining another proximity factor:
            double sector_back_factor =
            1. - fmod ( ( theta /* + pi*/) * sectors, ( segmentSize ) ) / ( segmentSize );
    //     cout << "sector_back_factor = " << sector_back_factor << endl;
            // Now we distribute the energy of the particle between
            // the nodes in the back section:
    //     cout << nodes.size() << ", trying: " << section_back*sectors + sector_back << endl;
            nodes[section_back*sectors + sector_back]
            ->addScalar ( particle_power
                          * section_back_factor
                          * sector_back_factor
                        );
            // Doing the same for the next node:
            int sector_front = sector_back+1;
            // check if we got into the last node of the section (last sector)
            if ( sector_front >= sectors ) sector_front = 0; // Turn around singularity
            nodes[section_back*sectors + sector_front]
            ->addScalar ( particle_power
                          * section_back_factor
                          * ( 1.-sector_back_factor )
                        );
            // And the same for the next section (two other nodes):
            int section_front = section_back+1;
            // Repeating the same:
            // sector_back = floor ( ( theta + pi ) * sectors / ( finalAngle ) ); // Commented as it's doing the same calculation as before (???)
            nodes[section_front*sectors + sector_back]
            ->addScalar ( particle_power
                          * ( 1.-section_back_factor )
                          * sector_back_factor
                        );
            // (bis) Doing the same for the next node:
            sector_front = sector_back+1;
            // (bis) check if we got into the last node of the section (last sector)
            if ( sector_front >= sectors ) sector_front = 0; // Turn around singularity
    //     cout << nodes.size() << ", trying: " << section_front*sectors + sector_front << endl;
            nodes[section_front*sectors + sector_front]
            ->addScalar ( particle_power
                          * ( 1.-section_back_factor )
                          * ( 1.-sector_back_factor )
                        );
         }
      }
    }



void Cone3::outputTable()
{
    std::ofstream out ( name + "top_table.txt" );
}


void Cone3::outputPowerFile ( int particles )
{
    std::ofstream outFile ( name + "total_power.dat" );
    double power2D, z, totalPower ( 0. );
    int i,j;
    double sectionDif = length / ( sections - 1 );

    outFile << "# s (m)\t Power (W)\t TotalPower (W)" << "\n";

    for ( i=0; i<sections; ++i )
    {
        z = i * sectionDif;
        power2D = 0.;
//     radius = (length - z) * tan(slope);
        for ( j=0; j<sectors; ++j )   //closed chain
        {
            power2D += nodes[i*sectors + j]->getScalar();
        }
        totalPower += power2D;
        outFile << z/1000. << "\t" << power2D << "\t"<< totalPower << "\n";
    }
    cout << "Total Power * particles = " << totalPower << endl;
    cout << "Total Power = " << totalPower / particles << " MW" << endl;
}

void Cone3::outputPowerDensityFile()
{
//    std::ofstream outFile ( "power.dat" );
//    std::ofstream outFileParts ( "power_particles.dat" );
//    std::ofstream outFileAnsys ( "Ansys_power_1D.dat" );
//
//    double power2D, z, totalPower ( 0. );
//    int i,j;
//    double sectionDif = length / ( sections - 1 );
//
//    outFile << "# s (m)\t Power Density (MW/m^2)\t Power in x=0 plane\t Power in y=0 plane" << "\n";
//    outFileAnsys << "# s (m)\t Power Density (W/m^2)" << "\n";
//
//    for ( i=0; i<sections; ++i )
//    {
//        z = i * sectionDif + z0;
//        power2D = 0.;
////     radius = (length - z) * tan(slope);
//        for ( j=0; j<sectors; ++j )   //closed chain
//        {
//            power2D += nodes[i*sectors + j]->getDensity();
//            outFileParts << "Node: " << i*sectors + j
//            << ", density = " << nodes[i*sectors + j]->getDensity()
//            << endl;
//        }
//        totalPower += power2D / sectors / sections;
//        outFile << z/1000. << "\t" << power2D / sectors << "\t";
//        outFile << nodes[i*sectors+(sectors+1)/4]->getDensity() << "\t";
//        outFile << nodes[i*sectors]->getDensity() << "\t";
//        outFile << nodes[i*sectors+(sectors+1)/2]->getDensity() << "\t";
//        outFile << nodes[i*sectors+(sectors+1)*3/4]->getDensity() << "\t";
//        outFile << "\n";
//        outFileAnsys << z/1000. << "\t" << power2D*1E6 / sectors << "\n";
//    }
//    cout << "Mean Power density = " << totalPower << " MW/m^2" << endl;
//      std::map< double, std::vector<double> >::iterator it_powers = powers.begin();
//      std::vector< Element* >::iterator it_elements = elements.begin();
//      std::vector<int> connectivity;
//      double sectPower,totalPower=0.;
//      std::ofstream outFile("power.dat");
//      std::ofstream outFile3D("power3D.dat");
//      double x,y,z;
//      int numberOfnodes = (*it_elements)->getNumberOfNodes(); // = 4
//
//      outFile3D << setiosflags( ios::fixed );
//      cout<<setprecision(10);
//
//      if (it_powers != powers.end() ){
//        for( it_powers;
//             it_powers != powers.end();
//             ++it_powers
//           )
//        {
//          sectPower = 0.;
//      //     sectArea = 0.;
//          for(int i=0; i<sectors; ++i){
//            x = y = z = 0.;
//            totalPower += it_powers->second.operator[](i);
//            sectPower += it_powers->second.operator[](i) / (*it_elements)->getArea();
//      //       sectArea += (*it_elements)->getArea();
//            connectivity = (*it_elements)->getConnectivity();
//            for( int j=0; j<connectivity.size(); ++j ){
//              if( connectivity[j] < nodes.size() ){
//    //          x += (1./numberOfNodes) * nodes[ connectivity[j] ]->getX();
//                x += 0.25 * nodes[ connectivity[j] ]->getX();
//                y += 0.25 * nodes[ connectivity[j] ]->getY();
//                z += 0.25 * nodes[ connectivity[j] ]->getZ();
//              }
//            }
//            outFile3D
//                << setw(20) << x << "\t"
//                << setw(20) << y << "\t"
//                << setw(20) << z << "\t"
//                << setw(20) << it_powers->second.operator[](i) / (*it_elements)->getArea()
//                << endl;
//            ++it_elements;
//          }
//          sectPower /= sectors;
//          outFile << it_powers->first << "\t" << sectPower << endl;
//        }
//      }
//      else{ // computed by theoric approach
//        std::vector< Element* >::iterator it_elements = elements.begin();
//        std::vector<int> connectivity;
//        std::map<int, Node*>::iterator it_nodes = nodes.begin();
//
//        for( it_elements;
//             it_elements != elements.end();
//             ++it_elements
//           )
//        {
//          connectivity = (*it_elements)->getConnectivity();
//          for( int j=0; j<connectivity.size(); ++j ){
//            totalPower += 0.25 * nodes[ connectivity[j] ]->getScalar() * (*it_elements)->getArea();
//          }
//        }
//
//        for( it_nodes;
//             it_nodes != nodes.end();
//             ++it_nodes
//           )
//        {
//          outFile << it_nodes->second->getX() << "\t"
//              << it_nodes->second->getY() << "\t"
//              << it_nodes->second->getZ() << "\t"
//              << it_nodes->second->getScalar() << std::endl;
//          ++it_nodes; //get only even nodes
//        }
//      }
//      cout << "TOTAL POWER = " << totalPower << endl;
}


