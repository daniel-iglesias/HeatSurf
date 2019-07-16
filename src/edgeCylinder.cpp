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

#include "edgeCylinder.h"

// #include <iomanip.h>

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>

edgeCylinder::edgeCylinder()
{
}

edgeCylinder::edgeCylinder ( std::string type_in,
                     double par0_in,
                     double par1_in,
                     double par2_in,
                     int par3_in )
        : Geometry ( type_in, par0_in, par3_in, par1_in )
        , initDiam ( par2_in )
{
    slope = 0.;
    gridWidth = initDiam/2.;
}


edgeCylinder::~edgeCylinder()
{
}


void edgeCylinder::setSections ( double distance )
{
    sections = length /distance + 1 ;
    std::cout << "Sections = " << sections << std::endl;
}


void edgeCylinder::computeGeometry()
{
    int i, j;
    double pi = 3.1416;
    double theta;
    double radius, ra, rb;
    double a, b, c;
    double x, y, z, zlocal;
    double sectionDif = length / ( sections - 1 );

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
            theta = 0.3915 * j / sectors - pi; // angle between -pi and pi
            x = radius * cos ( theta );
            y = radius * sin ( theta );
            nodes[i*sectors + j] = new Node ( x, y, z );
//       cout << x<< ", "<< y<<", "<< z<<endl;
        }
    }

    theta = 0.3915 / sectors; // angle between -pi and pi

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


void edgeCylinder::residue ( lmx::Vector<double>& res, lmx::Vector<double>& conf )
{
    double t = conf.readElement ( 0 );
    res.writeElement (
        ( pow ( x+vx*t,2 ) +pow ( y+vy*t,2 ) - pow ( initDiam/2.,2 ) )
        , 0
    );
}

void edgeCylinder::jacobian ( lmx::Matrix<double>& jac, lmx::Vector<double>& conf )
{
    double t = conf.readElement ( 0 );
    jac.writeElement (
        ( 2*vx* ( x+vx*t ) +2*vy* ( y+vy*t ) )
        , 0
        , 0
    );
}

void edgeCylinder::computeIntersection ( Particle* particle )
{
    x = particle->getX();
    y = particle->getY();
    z = particle->getZ();
    vx = particle->getXdiv();
    vy = particle->getYdiv();
    vz = particle->getZdiv();
    lmx::Vector<double> initialGuess ( 1 ); // zero
    lmx::NLSolver<edgeCylinder> theSolver;
    theSolver.setInfo ( 0 );
    theSolver.setInitialConfiguration ( initialGuess );
    theSolver.setSystem ( *this );
    theSolver.setResidue ( &edgeCylinder::residue );
    theSolver.setJacobian ( &edgeCylinder::jacobian );
//   theSolver.setConvergence( &edgeCylinder::myConvergence );
//   theSolver.setMaxIterations( 100 );
    theSolver.solve ( 100 );
//   cout << theSolver.getSolution().readElement(0) << endl;
    paramTrajectories.push_back ( theSolver.getSolution().readElement ( 0 ) );

}

void edgeCylinder::computeNodalPower ( Particle* particle )
{
    if ( paramTrajectories.back() > ( z0-particle->getZ() ) &&
            paramTrajectories.back() <= ( z0-particle->getZ() +length ) )
    {
        // First we compute the intersection point of the last particle computed
        // in the "computeIntersection" function:
        x = particle->getX() + particle->getXdiv() * paramTrajectories.back();
        y = particle->getY() + particle->getYdiv() * paramTrajectories.back();
        z = particle->getZ() + paramTrajectories.back() - z0;
//     cout << "particle->getZ() = " << particle->getZ() << endl;
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
        double theta = atan2 ( x,y );
        double pi = 3.1416;
        int sector_back = floor ( ( theta + pi ) * sectors / ( 2*pi ) );
//     cout << "sector_back = " << sector_back << endl;
        // ... defining another proximity factor:
        double sector_back_factor
        = 1. - fmod ( ( theta + pi ) * sectors, ( 2*pi ) ) / ( 2*pi );
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
        sector_back = floor ( ( theta + pi ) * sectors / ( 2*pi ) );
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


void edgeCylinder::outputTable()
{
    std::ofstream out ( "top_table.txt" );
}


void edgeCylinder::outputPowerFile ( int particles )
{
    std::ofstream outFile ( "total_power.dat" );
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

void edgeCylinder::outputPowerDensityFile()
{
    std::ofstream outFile ( "power.dat" );
    std::ofstream outFileParts ( "power_particles.dat" );
    std::ofstream outFileAnsys ( "Ansys_power_1D.dat" );

    double power2D, z, totalPower ( 0. );
    int i,j;
    double sectionDif = length / ( sections - 1 );

    outFile << "# s (m)\t Power Density (MW/m^2)" << "\n";
    outFileAnsys << "# s (m)\t Power Density (W/m^2)" << "\n";

    for ( i=0; i<sections; ++i )
    {
        z = i * sectionDif + z0;
        outFileParts << "Section: " << z << endl;
        power2D = 0.;
//     radius = (length - z) * tan(slope);
        for ( j=0; j<sectors; ++j )   //closed chain
        {
            power2D += nodes[i*sectors + j]->getDensity();
            outFileParts << "Node: " << i*sectors + j
            << ", x = " << nodes[i*sectors + j]->getX()
            << ", y = " << nodes[i*sectors + j]->getY()
            << ", scalar = " << nodes[i*sectors + j]->getScalar()
            << ", density = " << nodes[i*sectors + j]->getDensity()
            << endl;
        }
        totalPower += power2D / sectors / sections;
        outFile << z/1000. << "\t" << power2D / sectors << "\n";
        outFileAnsys << z/1000. << "\t" << power2D*1E6 / sectors << "\n";
    }
    cout << "Mean Power density = " << totalPower << " MW/m^2" << endl;
    /*  std::map< double, std::vector<double> >::iterator it_powers = powers.begin();
      std::vector< Element* >::iterator it_elements = elements.begin();
      std::vector<int> connectivity;
      double sectPower,totalPower=0.;
      std::ofstream outFile("power.dat");
      std::ofstream outFile3D("power3D.dat");
      double x,y,z;
      int numberOfnodes = (*it_elements)->getNumberOfNodes(); // = 4

      outFile3D << setiosflags( ios::fixed );
      cout<<setprecision(10);

      if (it_powers != powers.end() ){
      for( it_powers;
      it_powers != powers.end();
      ++it_powers
           )
      {
      sectPower = 0.;
      //     sectArea = 0.;
      for(int i=0; i<sectors; ++i){
      x = y = z = 0.;
      totalPower += it_powers->second.operator[](i);
      sectPower += it_powers->second.operator[](i) / (*it_elements)->getArea();
      //       sectArea += (*it_elements)->getArea();
      connectivity = (*it_elements)->getConnectivity();
      for( int j=0; j<connectivity.size(); ++j ){
      if( connectivity[j] < nodes.size() ){
    //          x += (1./numberOfNodes) * nodes[ connectivity[j] ]->getX();
      x += 0.25 * nodes[ connectivity[j] ]->getX();
      y += 0.25 * nodes[ connectivity[j] ]->getY();
      z += 0.25 * nodes[ connectivity[j] ]->getZ();
    }
    }
      outFile3D
      << setw(20) << x << "\t"
      << setw(20) << y << "\t"
      << setw(20) << z << "\t"
      << setw(20) << it_powers->second.operator[](i) / (*it_elements)->getArea()
      << endl;
      ++it_elements;
    }
      sectPower /= sectors;
      outFile << it_powers->first << "\t" << sectPower << endl;
    }
    }
      else{ // computed by theoric approach
      std::vector< Element* >::iterator it_elements = elements.begin();
      std::vector<int> connectivity;
      std::map<int, Node*>::iterator it_nodes = nodes.begin();

      for( it_elements;
      it_elements != elements.end();
      ++it_elements
           )
      {
      connectivity = (*it_elements)->getConnectivity();
      for( int j=0; j<connectivity.size(); ++j ){
      totalPower += 0.25 * nodes[ connectivity[j] ]->getScalar() * (*it_elements)->getArea();
    }
    }

      for( it_nodes;
      it_nodes != nodes.end();
      ++it_nodes
           )
      {
      outFile << it_nodes->second->getX() << "\t"
      << it_nodes->second->getY() << "\t"
      << it_nodes->second->getZ() << "\t"
      << it_nodes->second->getScalar() << std::endl;
      ++it_nodes; //get only even nodes
    }
    }
      cout << "TOTAL POWER = " << totalPower << endl;*/
}


