// /***************************************************************************
//  *   Copyright (C) 2007 by Daniel Iglesias   *
//  *   daniel.iglesias@ciemat.es   *
//  *                                                                         *
//  *   This program is free software; you can redistribute it and/or modify  *
//  *   it under the terms of the GNU General Public License as published by  *
//  *   the Free Software Foundation; either version 2 of the License, or     *
//  *   (at your option) any later version.                                   *
//  *                                                                         *
//  *   This program is distributed in the hope that it will be useful,       *
//  *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
//  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
//  *   GNU General Public License for more details.                          *
//  *                                                                         *
//  *   You should have received a copy of the GNU General Public License     *
//  *   along with this program; if not, write to the                         *
//  *   Free Software Foundation, Inc.,                                       *
//  *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
//  ***************************************************************************/
 #include "twoplates.h"
 #include <vtkPoints.h>
 #include <vtkUnstructuredGrid.h>
 #include <vtkCell.h>

 TwoPlates::TwoPlates()
  : Geometry()
 {
 }


 TwoPlates::TwoPlates(std::string type_in,
                      double length,
                      double initDiam,
                      int sectors,
                      double width_in,
                      double orientation_in)
   : Geometry(type_in, 0.0, sectors, length)
     , initSep(initDiam)
     , width(width_in)
     , orientation(orientation_in)
 {
   slope = atan2( initSep/2, length );
   std::cout << "Slope = " << slope << std::endl;
 }


 TwoPlates::~TwoPlates()
 {
 }


 void TwoPlates::computeEnergy(double section, std::vector< Particle* > & particles)
 {
   double height = initSep * (length - section) / ( 2. * length );
   double leftX, rightX, y_max, x, y;
   std::vector< Particle* >::iterator it;

   cout << "Height = " << height << endl;

   for(int i=0; i<sectors; ++i){
     energies[section].push_back(0.);
     energies[section].push_back(0.);
     powers[section].push_back(0.);
     powers[section].push_back(0.);
     leftX =  i     * width / sectors - width/2;
     rightX = (i+1) * width / sectors - width/2;

     for( it = particles.begin();
          it!= particles.end();
          ++it
        )
     {
       x = (*it)->getX();
       y = (*it)->getY();
       if( y > y_max) y_max = y;
       if( y < height && y >= 0 ){
         if( x >= leftX && x < rightX ){
           energies[section][2*i] += (*it)->getEnergy(); // total energy
 //           ++sectParts;
 //         cout << "energy in section " << section << ", sector " << i << endl;
         }
       }
       else if( y > -height && y < 0 ){
         if( x >= leftX && x < rightX ){
           energies[section][2*i+1] += (*it)->getEnergy(); // total energy
 //           ++sectParts;
 //         cout << "energy in section " << section << ", sector " << i << endl;
         }
       }
     }
     powers[section][2*i] += ( energies[section][2*i] * 0.125 / particles.size() ) * 1E6; // [Watts]
     powers[section][2*i+1] += ( energies[section][2*i+1] * 0.125 / particles.size() ) * 1E6; // [Watts]
   }

   cout << "Particle_max_height = " << y_max << endl;
   sections = energies.size() + 1; // we add the first section here

   cout << "Energies = " << energies[section].size() << endl;
   cout << "Sectors = " << sectors << endl;
 //   sectors *= 2; // for Geometry compatibility
 }


 void TwoPlates::setSections(double distance)
 {
   sections = length /distance + 1 ;
   std::cout << "Sections = " << sections << std::endl;
 }


 void TwoPlates::computeGeometry()
 {
 //   sectors /= 2; // for Geometry compatibility

   int i, j;
   double ra, rb;
   double a, c;
   double x, y, z;
   double sectionDif = length / (sections - 1);

   for(i=0; i<sections-1; ++i){ // last section is different
     z = i * sectionDif;
     y = initSep * ( length - z ) / (2. * length);
     for(j=0; j<=sectors; ++j){ //open chain
       x = j * width / sectors - width/2;
       nodes[2*i*(sectors+1) + 2*j]   = new Node( x,  y, z );
       nodes[2*i*(sectors+1) + 2*j+1] = new Node( x, -y, z );
     }
   }

   i = sections-1; // last section
   z = i * sectionDif;
   y = initSep * ( length - z ) / (2. * length);
   for(j=0; j<=sectors; ++j){ //open chain
     x = j * width / sectors - width/2;
     nodes[2*i*(sectors+1) + 2*j]   = new Node( x,  y, z );
     nodes[2*i*(sectors+1) + 2*j+1] = nodes[2*i*(sectors+1) + 2*j]; //instead of creating new nodes, we assign the vertex nodes (y=0)
   }


   for(i=0; i<sections-1; ++i){
     z = i * sectionDif;
     ra = initSep * ( length - z ) / (2. * length);
     rb = initSep * ( length - (z + sectionDif) ) / (2. * length);
     a = width / sectors;
     c = sqrt( pow(ra-rb, 2) + pow(sectionDif, 2) );
     for(j=0; j<sectors; ++j){//last sector is equal
       elements.push_back( new Element( 0 ) // type quad
                         );
       elements.back()->setNumber( 2*i*sectors + 2*j );
       elements.back()->addNode( 2*i*(sectors+1) + 2*j );
       elements.back()->addNode( 2*i*(sectors+1) + 2*(j+1) );
       elements.back()->addNode( 2*(i+1)*(sectors+1) + 2*(j+1) );
       elements.back()->addNode( 2*(i+1)*(sectors+1) + 2*j );
       elements.back()->generateGeometry();
       elements.back()->setArea( a * c );

       elements.push_back( new Element( 0 ) // type quad
                         );
       elements.back()->setNumber( 2*i*sectors + 2*j+1 );
       elements.back()->addNode( 2*i*(sectors+1) + 2*j+1 );
       elements.back()->addNode( 2*i*(sectors+1) + 2*(j+1)+1 );
       elements.back()->addNode( 2*(i+1)*(sectors+1) + 2*(j+1)+1 );
       elements.back()->addNode( 2*(i+1)*(sectors+1) + 2*j+1 );
       elements.back()->generateGeometry();
       elements.back()->setArea( a * c );
 //       cout << "a = "<< a << ", c = "<< c << ", Area = " << elements.back()->getArea() << endl;
 //       cout << "elem num: " << 2*i*sectors + 2*j << endl;
 //       cout << "elem num: " << 2*i*sectors + 2*j+1 << endl;
     }
   }

   gridPoints = vtkPoints::New();
   gridPoints->SetNumberOfPoints( nodes.size() );

   for(std::map<int, Node*>::iterator it = nodes.begin();
       it != nodes.end();
       ++it
      )
   {
     gridPoints->InsertPoint( it->first, it->second->getX(), it->second->getY(), it->second->getZ() );
 //     cout << "x = "<< it->second->getX() << ", y = "<< it->second->getY() << ", z = " << it->second->getZ() << endl;

   }

   grid = vtkUnstructuredGrid::New();
   grid->Allocate(1000,1000);

   for(std::vector<Element*>::iterator it=elements.begin();
       it!= elements.end();
       ++it
      )
   {
     grid->InsertNextCell( (*it)->geometry->GetCellType(), (*it)->geometry->GetPointIds() );
   }

   grid->SetPoints( gridPoints );

   sectors *= 2; // for Geometry compatibility
 }


void TwoPlates::residue ( lmx::Vector<double>& res, lmx::Vector<double>& conf )
{
    double t = conf.readElement ( 0 );
    res.writeElement (
        ( pow ( x+vx*t,2 ) +pow ( y+vy*t,2 ) ) *pow ( cos ( slope ),2 )
        -pow ( t- ( length+z0 ),2 ) *pow ( sin ( slope ),2 )
        , 0
    );
}

void TwoPlates::jacobian ( lmx::Matrix<double>& jac, lmx::Vector<double>& conf )
{
    double t = conf.readElement ( 0 );
    jac.writeElement (
        ( 2*vx* ( x+vx*t ) +2*vy* ( y+vy*t ) ) *pow ( cos ( slope ),2 )
        -2* ( t- ( length+z0 ) ) *pow ( sin ( slope ),2 )
        , 0
        , 0
    );
}

void TwoPlates::computeIntersection ( Particle* particle )
{
    x = particle->getX();
    y = particle->getY();
    z = particle->getZ();
    vx = particle->getXdiv();
    vy = particle->getYdiv();
    vz = particle->getZdiv();
    lmx::Vector<double> initialGuess ( 1 ); // zero
    lmx::NLSolver<TwoPlates> theSolver;
    theSolver.setInfo ( 0 );
    theSolver.setInitialConfiguration ( initialGuess );
    theSolver.setSystem ( *this );
    theSolver.setResidue ( &TwoPlates::residue );
    theSolver.setJacobian ( &TwoPlates::jacobian );
//   theSolver.setConvergence( &TwoPlates::myConvergence );
//   theSolver.setMaxIterations( 100 );
    theSolver.solve ( 100 );
//   cout << theSolver.getSolution().readElement(0) << endl;
    paramTrajectories.push_back ( theSolver.getSolution().readElement ( 0 ) );

}

void TwoPlates::computeNodalPower ( Particle* particle )
{
}


 void TwoPlates::outputTable()
 {
   std::ofstream out1("top_table.txt");
   std::ofstream out2("low_table.txt");
   std::ofstream out3("top_table.gnu");

   int i,j;
   double y=0, s=0;

   out1<< "0" << "\t";
   out2<< "0" << "\t";
   out3<< "S\tX\tPower Density[MW/m^2" << endl;
   for(i=0; i<=sectors; ++i){
     out1<< nodes[i]->getX()*1E-3 << "\t";
     ++i;
     out2<< nodes[i]->getX()*1E-3 << "\t";
   }
   for(i=0; i<nodes.size(); ++i){
     if( y != nodes[i]->getY() ){
       out1 << "\n";
       out2 << "\n";
       out3 << "\n";
       s = sqrt( pow( length - nodes[i]->getZ(), 2) + pow( nodes[i]->getY(), 2) );
       y = nodes[i]->getY();
       out1 << y*1E-3 << "\t";
       out2 << y*1E-3 << "\t";
     }
     out3 << s*1E-3 << "\t";
     out3 << nodes[i]->getX()*1E-3 << "\t";
     out3<< nodes[i]->getScalar() << "\n";

     out1<< nodes[i]->getScalar()*1E6 << "\t";
     ++i;
     out2<< nodes[i]->getScalar()*1E6 << "\t";
   }
 }


 void TwoPlates::outputPowerFile(int particles)
 {
   std::map< double, std::vector<double> >::iterator it_energies = energies.begin();
   std::vector< Element* >::iterator it_elements = elements.begin();
   std::vector<int> connectivity;
   double sectEnergy, totalEnergy=0.;
   std::ofstream outFile("energy.dat");

   for( it_energies;
        it_energies != energies.end();
        ++it_energies
      )
   {
     sectEnergy=0;
     for(int i=0; i<sectors; ++i){
       totalEnergy += it_energies->second.operator[](i);
       sectEnergy += it_energies->second.operator[](i) / (*it_elements)->getArea();
       ++it_elements;
     }
     outFile << it_energies->first << "\t" << sectEnergy << endl;
   }
   cout << "TOTAL ENERGY = " << totalEnergy << endl;
 }

 void TwoPlates::outputPowerDensityFile()
 {
   std::map< double, std::vector<double> >::iterator it_powers = powers.begin();
   std::vector< Element* >::iterator it_elements = elements.begin();
   std::vector<int> connectivity;
   double sectPower,totalPower=0., totalArea=0.;
   std::ofstream outFile("power.dat");

   if (it_powers != powers.end() ){
     for( it_powers;
          it_powers != powers.end();
          ++it_powers
        )
     {
       sectPower = 0.;
   //     sectArea = 0.;
       for(int i=0; i<sectors; ++i){
         totalPower += it_powers->second.operator[](i);
         sectPower += it_powers->second.operator[](i) / (*it_elements)->getArea();
   //       sectArea += (*it_elements)->getArea();
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
         totalArea += 0.25 * (*it_elements)->getArea() * sin(slope); //verification -> Verified!
       }
     }

     outFile << "x [m]" << "\t"
 //           << it_nodes->second->getY()/1000. << "\t"
         << "z [m]" << "\t"
         << "Power [W/m^2] " << std::endl;

     for( it_nodes;
          it_nodes != nodes.end();
          ++it_nodes
        )
     {
       outFile << it_nodes->second->getX()/1000. << "\t"
 //           << it_nodes->second->getY()/1000. << "\t"
           << it_nodes->second->getZ()/1000. << "\t"
           << it_nodes->second->getScalar()*1E6 << std::endl;
       ++it_nodes; //get only even nodes
     }
   }
   cout << "TOTAL PROJECTED AREA = " << totalArea << endl;
   cout << "TOTAL POWER = " << totalPower << endl;
 }


