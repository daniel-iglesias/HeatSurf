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

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderer.h>

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkFloatArray.h>
#include <vtkDataSetMapper.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkCaptionActor2D.h>
#include <vtkMapper.h>

#include <vtkPlane.h>
#include <vtkCutter.h>
#include <vtkPolyDataMapper.h>

#include <vtkUnstructuredGridWriter.h>

#include "geometry.h"

Geometry::Geometry()
{
}

Geometry::Geometry ( std::string type_in,
                     std::string name_in,
                     double z0_in,
                     int sectors_in,
                     double length_in)
        : type ( type_in )
        , name ( name_in )
        , z0 ( z0_in )
        , sectors ( sectors_in )
        , length ( length_in )
        , shift_x(0)
        , shift_y(0)
        , shift_z(0)
        , xVec_Rotation(0.0001)
        , yVec_Rotation(0.0001)
        , zVec_Rotation(0.0001)
        , angleRotation(0.0)
        , u(1), v(0), w(0)
        , angle(0.0)
{
  rotationMatrix[0][0] = 1.0; // Initalising Rotation matrix
  rotationMatrix[0][1] = 0.0;
  rotationMatrix[0][2] = 0.0;
  rotationMatrix[0][3] = 0.0;

  rotationMatrix[1][0] = 0.0;
  rotationMatrix[1][1] = 1.0;
  rotationMatrix[1][2] = 0.0;
  rotationMatrix[1][3] = 0.0;

  rotationMatrix[2][0] = 0.0;
  rotationMatrix[2][1] = 0.0;
  rotationMatrix[2][2] = 1.0;
  rotationMatrix[2][3] = 0.0;

  rotationMatrix[3][0] = 0.0;
  rotationMatrix[3][1] = 0.0;
  rotationMatrix[3][2] = 0.0;
  rotationMatrix[3][3] = 1.0;
}

Geometry::~Geometry()
{
}


void Geometry::setGeometryShift(
       double geometry_shift_x,
       double geometry_shift_y,
       double geometry_shift_z,
       double geometry_shift_mag )
{
shift_x = geometry_shift_x;
shift_y = geometry_shift_y;
shift_z = geometry_shift_z;
shift_mag = geometry_shift_mag;
}

void Geometry::shiftGeometry( )
{
// Take Node x,y and z coordinates
// Then add the shift in the x,y,z plane multiplied
// with the shift magnitude to the relative coordinate

// Shift nodes
    for (auto it= nodes.begin(); it!= nodes.end(); ++it){
        it->second->setX( it->second->getX() + (shift_x * shift_mag) );
        it->second->setY( it->second->getY() + (shift_y * shift_mag) );
        it->second->setZ( it->second->getZ() + (shift_z * shift_mag) );
        }

// Shift grid3Dnodes
    for (auto it= grid3Dnodes.begin(); it!= grid3Dnodes.end(); ++it){
        it->second->setX( it->second->getX() + (shift_x * shift_mag) );
        it->second->setY( it->second->getY() + (shift_y * shift_mag) );
        it->second->setZ( it->second->getZ() + (shift_z * shift_mag) );
        }

// Shift VTK node objects:

// Shifting gridPoints:
double vec[2];
double NewVecX;
double NewVecY;
double NewVecZ;

    for (int i = 0; i<gridPoints->GetNumberOfPoints(); i++){
        gridPoints->GetPoint(i, vec);
        NewVecX = vec[0] + (shift_x * shift_mag);
        NewVecY = vec[1] + (shift_y * shift_mag);
        NewVecZ = vec[2] + (shift_z * shift_mag);
        gridPoints->SetPoint(i, NewVecX, NewVecY, NewVecZ);
        }

// Shifting grid3DPoints:
double vec2[2];
double NewVecX2;
double NewVecY2;
double NewVecZ2;

    for (int i = 0; i<grid3DPoints->GetNumberOfPoints(); i++){
        grid3DPoints->GetPoint(i, vec2);
        NewVecX2 = vec2[0] + (shift_x * shift_mag);
        NewVecY2 = vec2[1] + (shift_y * shift_mag);
        NewVecZ2 = vec2[2] + (shift_z * shift_mag);
        grid3DPoints->SetPoint(i, NewVecX2, NewVecY2, NewVecZ2);
        }
}

void Geometry::setGeometryRotation(
        double rotationVectorX,
        double rotationVectorY,
        double rotationVectorZ,
        double rotationAngle )
{
xVec_Rotation = rotationVectorX;
yVec_Rotation = rotationVectorY;
zVec_Rotation = rotationVectorZ;
angleRotation = rotationAngle;
}

void Geometry::rotateGeometry ( )
{
// Rotating each individual coordinate around a point,
// around a unit vector plane and a set rotation angle
// Will need to rotate nodes, grid3Dnodes, gridPoints, grid3DPoints

// Shifting nodes
    for (auto it= nodes.begin(); it!= nodes.end(); ++it)
        {
//            inputMatrix[0][0] = it->second->getX();
//            inputMatrix[1][0] = it->second->getY();
//            inputMatrix[2][0] = it->second->getZ();
//            inputMatrix[3][0] = 1.0;
//
//            for(int i = 0; i < 4; i++ ){
//                for(int j = 0; j < 1; j++){
//                    outputMatrix[i][j] = 0;
//                        for(int k = 0; k < 4; k++){
//                            outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
//                    }
//                }
//            }
//        it->second->setX( outputMatrix[0][0] ) ;
//        it->second->setY( outputMatrix[1][0] ) ;
//        it->second->setZ( outputMatrix[2][0] ) ;
            std::vector<double> rotated_xyz = rotatePointAroundGeometryAxis(it->second->getX(), it->second->getY(), it->second->getZ());

            it->second->setX( rotated_xyz[0] ) ;
            it->second->setY( rotated_xyz[1] ) ;
            it->second->setZ( rotated_xyz[2] ) ;


//        cout<< outputMatrix[0][0] /*<< ", "<< outputMatrix[1][0] << ", "<< outputMatrix[2][0] */<< endl;
    }

// Shifting grid3Dnodes
    for (auto it= grid3Dnodes.begin(); it!= grid3Dnodes.end(); ++it)
    {
//        inputMatrix[0][0] = it->second->getX();
//        inputMatrix[1][0] = it->second->getY();
//        inputMatrix[2][0] = it->second->getZ();
//        inputMatrix[3][0] = 1.0;
//
//        for(int i = 0; i < 4; i++ ){
//            for(int j = 0; j < 1; j++){
//                outputMatrix[i][j] = 0;
//                    for(int k = 0; k < 4; k++){
//                        outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
//                }
//            }
//        }
//        it->second->setX( outputMatrix[0][0] ) ;
//        it->second->setY( outputMatrix[1][0] ) ;
//        it->second->setZ( outputMatrix[2][0] ) ;
            std::vector<double> rotated_xyz = rotatePointAroundGeometryAxis(it->second->getX(), it->second->getY(), it->second->getZ());

            it->second->setX( rotated_xyz[0] ) ;
            it->second->setY( rotated_xyz[1] ) ;
            it->second->setZ( rotated_xyz[2] ) ;
    }

// Shifting gridPoints
    double gridPointsVec[2];
    for (int a = 0; a<gridPoints->GetNumberOfPoints(); a++)
    {
        gridPoints->GetPoint(a, gridPointsVec);
//        inputMatrix[0][0] = gridPointsVec[0];
//        inputMatrix[1][0] = gridPointsVec[1];
//        inputMatrix[2][0] = gridPointsVec[2];
//        inputMatrix[3][0] = 1.0;
//
//        for(int i = 0; i < 4; i++ ){
//            for(int j = 0; j < 1; j++){
//                outputMatrix[i][j] = 0;
//                    for(int k = 0; k < 4; k++){
//                        outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
//                }
//            }
//        }
//
//        gridPoints->SetPoint(a, outputMatrix[0][0], outputMatrix[1][0], outputMatrix[2][0]);
//        cout<< inputMatrix[0][0] << ", "<< inputMatrix[1][0] << ", "<< inputMatrix[2][0] << endl;
//        cout<< outputMatrix[0][0] << ", "<< outputMatrix[1][0] << ", "<< outputMatrix[2][0] << endl;
//        cout<< rotationMatrix[0][2] << ", "<< rotationMatrix[1][0] << ", "<< rotationMatrix[2][0] << endl;

            std::vector<double> rotated_xyz = rotatePointAroundGeometryAxis(gridPointsVec[0], gridPointsVec[1], gridPointsVec[2]);

            gridPoints->SetPoint(a, rotated_xyz[0], rotated_xyz[1], rotated_xyz[2]);
    }

// Shifting grid3DPoints
    double grid3DPointsVec[2];
    for (int a = 0; a<grid3DPoints->GetNumberOfPoints(); a++)
    {
        grid3DPoints->GetPoint(a, grid3DPointsVec);
//        inputMatrix[0][0] = grid3DPointsVec[0];
//        inputMatrix[1][0] = grid3DPointsVec[1];
//        inputMatrix[2][0] = grid3DPointsVec[2];
//        inputMatrix[3][0] = 1.0;
//
//        for(int i = 0; i < 4; i++ ){
//            for(int j = 0; j < 1; j++){
//                outputMatrix[i][j] = 0;
//                    for(int k = 0; k < 4; k++){
//                        outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
//                }
//            }
//        }
//
//        grid3DPoints->SetPoint(a, outputMatrix[0][0], outputMatrix[1][0], outputMatrix[2][0]);
            std::vector<double> rotated_xyz = rotatePointAroundGeometryAxis(grid3DPointsVec[0], grid3DPointsVec[1], grid3DPointsVec[2]);

            grid3DPoints->SetPoint(a, rotated_xyz[0], rotated_xyz[1], rotated_xyz[2]);

    }
}


void Geometry::computeGrid3D( )
{
    sectors3D = 2*sectors/3;
    sections3D = 2*sections/3;
    int i,j,k;
    double x, y, z;
    int sectNodes = sectors3D+1;
    double sectionDif = length / ( sections3D - 1 );

    // Define a regular grid of nodes:
    // Delta_z = sectionDif
    // Delta_y = Delta_x = 2*(radius/sectors)
    for ( i=0; i<sections3D; ++i )   // for z = i*sectionDif
    {
        z = (i * sectionDif + z0);
        for ( j=0; j<=sectors3D; ++j )   // for y = radius*( 2*j/sectors - 1 )
        {
            y = (gridWidth* ( 2.*j/sectors3D - 1. ));
            for ( k=0; k<=sectors3D; ++k )   // for x = radius*( 2*k/sectors - 1 )
            {
                x = (gridWidth* ( 2.*k/sectors3D - 1. ));
                grid3Dnodes[i*sectNodes*sectNodes + j*sectNodes + k]
                = new Node ( x, y, z );
//        cout<<x<<","<<y<<","<<z<<endl;

            }
        }
    }

    int elem_number;
    int sec_dim = sectors3D+1;
    // Define the 3D block elements (type-1)
    // Each one has two nodes in plane z and z+i*sectionDif
    // other two nodes in y and y
    for ( i=0; i<sections3D-1; ++i )   // for z = i*sectionDif
    {
        for ( j=0; j<sectors3D; ++j )   // for y = radius( 2*j/sectors - 1 )
        {
            for ( k=0; k<sectors3D; ++k )   // for x = radius( 2*k/sectors - 1 )
            {
                elem_number = i*sectors3D*sectors3D + j*sectors3D + k;
                grid3D.push_back ( new Element ( 1 ) // type brick
                                 );
                grid3D.back()->setNumber ( elem_number );
                // testing the vector numbering (can be deleted):
                if ( grid3D[elem_number]->getNumber() != elem_number )
                    cout << "ERROR storing grid3D elements" << endl;
                grid3D.back()->addNode ( i    *sec_dim*sec_dim +  j   *sec_dim + k );
                grid3D.back()->addNode ( i    *sec_dim*sec_dim +  j   *sec_dim + k+1 );
                grid3D.back()->addNode ( i    *sec_dim*sec_dim + ( j+1 ) *sec_dim + k+1 );
                grid3D.back()->addNode ( i    *sec_dim*sec_dim + ( j+1 ) *sec_dim + k );
                grid3D.back()->addNode ( ( i+1 ) *sec_dim*sec_dim +  j   *sec_dim + k );
                grid3D.back()->addNode ( ( i+1 ) *sec_dim*sec_dim +  j   *sec_dim + k+1 );
                grid3D.back()->addNode ( ( i+1 ) *sec_dim*sec_dim + ( j+1 ) *sec_dim + k+1 );
                grid3D.back()->addNode ( ( i+1 ) *sec_dim*sec_dim + ( j+1 ) *sec_dim + k );

                grid3D.back()->generateGeometry();

//         int node_temp;
//         for(int i=0; i<8; ++i){
//           node_temp = grid3D.back()->getConnectivity().operator[](i);
//           cout << "(" << grid3Dnodes[node_temp]->getX()
//                << "," << grid3Dnodes[node_temp]->getY()
//                << "," << grid3Dnodes[node_temp]->getZ() <<")";
//         }
//         cout << endl;

            }
        }
    }

    // We create the vtk points:
    grid3DPoints = vtkPoints::New();
    grid3DPoints->SetNumberOfPoints ( grid3Dnodes.size() );

    for ( std::map<int, Node*>::iterator it = grid3Dnodes.begin();
            it != grid3Dnodes.end();
            ++it
        )
    {
        grid3DPoints->InsertPoint ( it->first,
                                    it->second->getX(),
                                    it->second->getY(),
                                    it->second->getZ()
                                  );
//        cout<<it->second->getZ()<<endl;
    }
    // And, finally, the vtk cells:
    gridCubes = vtkUnstructuredGrid::New();
    gridCubes->Allocate ( 1000,1000 );

    for ( std::vector<Element*>::iterator it=grid3D.begin();
            it!= grid3D.end();
            ++it
        )
    {
        gridCubes->InsertNextCell ( ( *it )->geometry->GetCellType(),
                                    ( *it )->geometry->GetPointIds() );
    }
    gridCubes->SetPoints ( grid3DPoints );
}

void Geometry::computePowerDensity ( int numParticles )
{
    std::vector< Element* >::iterator it_elements = elements.begin();
    std::vector<int> connectivity;
    // First section nodes Scalar value must be doubled to simulate
    // the half of the element that isn't discretized
    int first_section = 0;
    for ( int j=0; j<sectors; ++j )   //closed chain
    {
        nodes[first_section*sectors + j]->addScalar
        ( nodes[first_section*sectors + j]->getScalar() );
    }

    // Last section nodes Scalar value must be doubled to simulate
    // the half of the element that isn't discretized
    int last_section = sections-1;
    for ( int j=0; j<sectors; ++j )   //closed chain
    {
        nodes[last_section*sectors + j]->addScalar
        ( nodes[last_section*sectors + j]->getScalar() );
    }

    for ( it_elements;
            it_elements != elements.end();
            ++it_elements )
    {
        connectivity = ( *it_elements )->getConnectivity();
        for ( int j=0; j<connectivity.size(); ++j )
        {
            ( *it_elements )->addDensity ( nodes[ connectivity[j] ]->getScalar()
                                           / numParticles
                                           / ( *it_elements )->getArea() *1E6
                                           / ( *it_elements )->getNumberOfNodes()
                                         );
            if ( ( *it_elements )->getDensity() > 1E10 )   // Error!
            {
                cout << "ERROR DETECTED, element density = "
                     << ( *it_elements )->getDensity()
                     << "nodes[ connectivity[j] ]->getScalar() = "
                     << nodes[ connectivity[j] ]->getScalar()
                     << "numParticles = " << numParticles
                     << "(*it_elements)->getArea() = " << ( *it_elements )->getArea()
                     << endl;
            }
        }
        for ( int j=0; j<connectivity.size(); ++j )
        {
            nodes[ connectivity[j] ]->addDensity ( ( *it_elements )->getDensity()
//                                          /(*it_elements)->getNumberOfNodes()
                                                 );
            // Nodes of last section have only two elements connected,
            // so we multiply by two (adding its own value):
            // [Resulting power.dat has a peak on vertex point, must check]
//       if ( connectivity[j] >= (sections-1)*sectors )
//         nodes[ connectivity[j] ]
//             ->addDensity( nodes[ connectivity[j] ]->getDensity()*.5 );

            if ( nodes[ connectivity[j] ]->getDensity() > 1E10 )  // Error!
            {
                cout << "ERROR DETECTED, node density = "
                     << nodes[ connectivity[j] ]->getDensity()
                     << "(*it_elements)->getDensity() = "
                     << ( *it_elements )->getDensity()
                     << "(*it_elements)->getNumberOfNodes() = "
                     << ( *it_elements )->getNumberOfNodes()
                     << endl;
            }
        }
    }

    // Now we transform to the 3D regular grid:
    // The idea is to use the Density values of the surface grid nodes and
    // searching the nearer grid3Dnodes. Then the density can be added to them
    // multiplied by a factor of proximity

    // First we assign the density value to the :
    std::map<int, Node*>::iterator it_nodes;
    std::map<int, Node*>::iterator it_nodes_end = nodes.end();
    double x,y,z;
    int i,j,k;
    double sectionDif = length / ( sections3D - 1 );
    int elem_number;

// Rotating the points back to its original position
    for ( it_nodes = nodes.begin(); it_nodes!= it_nodes_end; ++it_nodes )
    {
            std::vector<double> inverse_rotated_xyz = inverseRotatePointAroundGeometryAxis(it_nodes->second->getX(), it_nodes->second->getY(), it_nodes->second->getZ());
            it_nodes->second->setX( inverse_rotated_xyz[0] ) ;
            it_nodes->second->setY( inverse_rotated_xyz[1] ) ;
            it_nodes->second->setZ( inverse_rotated_xyz[2] ) ;
    }


    for ( it_nodes = nodes.begin(); it_nodes!= it_nodes_end; ++it_nodes )
    { // Finding the nodes in the local coordinate system
      // Need to subtract rotation shift to x, y, and z
      //cout<<"x: "<<it_nodes->second->getX()<<" y: "<<it_nodes->second->getY()<<" z: "<<it_nodes->second->getZ()<<endl;

        x = it_nodes->second->getX() - shift_x * shift_mag;
        y = it_nodes->second->getY() - shift_y * shift_mag;
        z = it_nodes->second->getZ() - shift_z * shift_mag;
//        cout<<"x: "<<x<<" y: "<<y<<" z: "<<z<<endl;
        i = floor ( ( z - z0 ) / sectionDif );
//          y = gridWidth*( 2.*j/sectors - 1. );
//          y / gridWidth +1 = ( 2.*j/sectors  );
        j = floor ( ( y/gridWidth + 1. ) * sectors3D/2. );
//          x / gridWidth = ( 2.*k/sectors - 1. );
        k = floor ( ( x/gridWidth + 1. ) * sectors3D/2. );

        elem_number = i*sectors3D*sectors3D + j*sectors3D + k;
//        cout << "x="<< x <<",y="<< y <<", z="<< z << endl;
//        cout << "i="<< i <<",j="<< j <<", k="<< k << endl;
//              <<"; elem_number="<< elem_number
//              <<"; elem_size="<< grid3D.size() << endl;
        if ( i < floor ( length/sectionDif ) )
            grid3D[elem_number]->addDensity3D ( it_nodes->second->getDensity() );
//        cout<<"got here"<<endl;
//        cout<<floor (length/sectionDif)<<endl;
    }

    // Now the mean density is passed to the eight nodes in the vertices of each
    // brick element using the special function node::addDensity3D(double).
    it_elements = grid3D.begin();
    for ( it_elements;
            it_elements != grid3D.end();
            ++it_elements )
    {
        connectivity = ( *it_elements )->getConnectivity();
        for ( int j=0; j<connectivity.size(); ++j )
        {
//       cout << "element density = " << (*it_elements)->getDensity3D() << endl;
            // Only add non-zero values!
            if ( ( *it_elements )->getDensity3D() > 0. )
                grid3Dnodes[ connectivity[j] ]
                ->addDensity3D ( ( *it_elements )->getDensity3D() );
        }
    }
}


void Geometry::outputAnsys3D()
{
    std::ofstream outFile3D ( name + "_Ansys_power_3D.dat" );
    int i,j,k; // z, y, x indexes
    std::vector<double> x_coords;
    std::vector<double>::iterator it_x_coords;
    int sectNodes = sectors3D+1;

    // Output in file the info of the grid:
    outFile3D << "# Grid dimension: "
    << sectors3D + 1 << " "
    << sectors3D + 1 << " "
    << sections3D << endl;
    // first we store the values of the x-coordinates of the regular grid:
    for ( k=0; k<=sectors3D; ++k )   // for x = radius*( 2*k/sectors3D - 1 )
    {
        x_coords.push_back (
            grid3Dnodes[k]->getX() *1E-3 );
    }

    for ( i=0; i<sections3D; ++i )   // for z = i*sectionDif
    {
        outFile3D << grid3Dnodes[i*sectNodes*sectNodes]->getZ() *1E-3;
        for ( it_x_coords = x_coords.begin();
                it_x_coords!= x_coords.end();
                ++it_x_coords
            )
            { outFile3D << "\t" << * ( it_x_coords ); }
        outFile3D << "\n";
        for ( j=0; j<=sectors3D; ++j )
        {
            outFile3D <<
            grid3Dnodes[i*sectNodes*sectNodes + j*sectNodes]->getY() *1E-3;
            for ( k=0; k<=sectors3D; ++k )
            {
                outFile3D << "\t"
                << grid3Dnodes[i*sectNodes*sectNodes + j*sectNodes + k]->getDensity() *1E6;
            }
            outFile3D << "\n";
        }
    }
    std::ofstream outFile2D ( name + "_Ansys_power_2D.dat" );
    double max_data = 0.0;

    outFile2D << "0.0"; // only one "plane" of data
    for ( it_x_coords = x_coords.begin();
            it_x_coords!= x_coords.end();
            ++it_x_coords
        )
        { outFile2D << "\t" << * ( it_x_coords ); }
    outFile2D << "\n";
    for ( i=0; i<sections3D; ++i )   // for z = i*sectionDif
    {
        outFile2D << grid3Dnodes[i*sectNodes*sectNodes]->getZ() *1E-3;
        for ( k=0; k<=sectors3D; ++k )
        {
            for ( j=0; j<=sectors3D; ++j )
            {
                max_data = std::max ( max_data,
                                      grid3Dnodes[i*sectNodes*sectNodes + j*sectNodes + k]
                                      ->getDensity() *1E6
                                    );
            }
            outFile2D << "\t" << max_data;
            max_data = 0.0;
        }
        outFile2D << "\n";
    }

    cout << "\"Ansys_power_3D.dat\" created,\n"
         << "Number of lines, columns and planes: \n"
         << sectors3D + 1 << ", "
         << sectors3D + 1 << ", "
         << sections3D
         << endl;

}

void Geometry::outputVTKfiles()
{
    std::string name1,name2;
    vtkUnstructuredGridWriter *writer1 = vtkUnstructuredGridWriter::New();
    writer1->SetInputData ( grid );
    name1 = name + "_grid.vtk";
    writer1->SetFileName ( name1.c_str() );
    writer1->Write();

    vtkUnstructuredGridWriter *writer2 = vtkUnstructuredGridWriter::New();
    writer2->SetInputData ( gridCubes );
    name2 = name + "_gridCubes.vtk";
    writer2->SetFileName ( name2.c_str() );
    writer2->Write();
}

// Not valid for plates geometry!
void Geometry::outputBackscattering ( double precision )
{
    std::ofstream ofile ( name + "_density.dat" );

    //supposing al the sections have the same lenght:
    //  number of division of each section:
    size_type slices = fabs ( length/ ( sections-1 ) ) /precision;
    std::vector< double > densities ( slices );
    std::vector< double > areas ( slices );

    double z, zi, zf, ri, rf; // axial and radial coordiantes for section nodes
    double d, di,df ( 0. ); // density for initial and final slices, computed adding nodes
    double zii, zff, rii, rff; // axial and radial coordiantes for section nodes
    double dii, dff; // density for initial and final nodes of the slice

    double area, power, total_power;

    int i,j,k;

    ofile << "# z(m)\t mean power density (MW/m^2)\t power (W/m)\t Area (mm^2)"
    << endl;

    for ( i=0; i<sections-1; ++i )
    {
        zi = nodes[i*sectors]->getZ();
        zf = nodes[ ( i+1 ) *sectors]->getZ();
        ri = fabs ( nodes[i*sectors]->getX() );
        rf = fabs ( nodes[ ( i+1 ) *sectors]->getX() );
        di = 0.; //could be df but careful with the first section
        df = 0.;
//     section_area = 0.;
//     section_power = 0.;
        for ( j=0; j<sectors; ++j )   //closed chain
        {
//       cout << i << ", " << j << endl;
            di += nodes[i*sectors + j]->getDensity();
            df += nodes[ ( i+1 ) *sectors + j]->getDensity();
//       section_area += elements[i*sectors + j]->getArea();
        }
        di /= sectors;
        df /= sectors;
//     cout << di << ", " << df << endl;
//     section_power = section_area*(di+df)/2.;
        for ( k=0; k<slices; ++k )
        {
            zii = zi+ ( ( double ) ( k ) /slices ) * ( zf-zi );
            zff = zi+ ( ( ( double ) ( k ) +1. ) /slices ) * ( zf-zi );
            rii = ri+ ( ( double ) ( k ) /slices ) * ( rf-ri );
            rff = ri+ ( ( ( double ) ( k ) +1. ) /slices ) * ( rf-ri );
            dii = di+ ( ( double ) ( k ) /slices ) * ( df-di );
            dff = di+ ( ( ( double ) ( k ) +1. ) /slices ) * ( df-di );

            area = 3.1416* ( rii+rff ) *sqrt ( pow ( rii-rff,2 ) + pow ( zff-zii,2 ) );
            power = area* ( dii+dff ) /2.;
            total_power += power;
            ofile << zii/1000. << "\t" //z
            //<< (zii + zff)/2. << "\t" //z
            << ( dii + dff ) /2. << "\t" //density
            << power << "\t" //power
            << area << endl; //area
        }
    }
    ofile.close();
    cout << "slices = " << slices
         << " = " << length<< " / (" << sections-1 << "*" << precision << ")"
         << endl;
    cout << "Output for Backscattering total power = "
         << total_power << endl;
}

void Geometry::calculateScalar()
{
    //   # Define the value of a scalar field (in this case the scalar)
    scalar = vtkFloatArray::New();

    for ( std::map<int,Node*>::iterator it = nodes.begin();
            it != nodes.end();
            ++it )
    {
        scalar->InsertNextValue ( it->second->getDensity() );
    }
    grid->GetPointData()->SetScalars ( scalar );

}

std::vector<double> Geometry::rotatePointAroundGeometryAxis(double x_in, double y_in, double z_in)
{
    double u, v, w, angle;
    double inputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};
    double outputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};
    u = xVec_Rotation;
    v = yVec_Rotation;
    w = zVec_Rotation;
    angle = angleRotation;
    std::vector<double> rotated_coordinates;

    if (u > 0.001 || v > 0.001 || w > 0.001 ) //making sure this isnt used if there is no rotation
    {
// Setting up Rotation matrix
        double L = (u * u + v * v + w * w);
        angle = angle * M_PI / 180.0; //converting to radian value
        double u2 = u * u;
        double v2 = v * v;
        double w2 = w * w;

        rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
        rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
        rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
        rotationMatrix[0][3] = 0.0;

        rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
        rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
        rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
        rotationMatrix[1][3] = 0.0;

        rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
        rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
        rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
        rotationMatrix[2][3] = 0.0;

        rotationMatrix[3][0] = 0.0;
        rotationMatrix[3][1] = 0.0;
        rotationMatrix[3][2] = 0.0;

        inputMatrix[0][0] = x_in;
        inputMatrix[1][0] = y_in;
        inputMatrix[2][0] = z_in;
        inputMatrix[3][0] = 1.0;

        for(int i = 0; i < 4; i++ ){
            for(int j = 0; j < 1; j++){
                outputMatrix[i][j] = 0;
                for(int k = 0; k < 4; k++){
                    outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
                }
            }
        }
        rotated_coordinates.push_back( outputMatrix[0][0] ); // Rounding done since there is a rounding error for '0'
        rotated_coordinates.push_back( outputMatrix[1][0] );
        rotated_coordinates.push_back( outputMatrix[2][0] );
    }
    else{
        rotated_coordinates.push_back( x_in );
        rotated_coordinates.push_back( y_in );
        rotated_coordinates.push_back( z_in );
    }
    return rotated_coordinates;
}

std::vector<double> Geometry::inverseRotatePointAroundGeometryAxis(double x_in, double y_in, double z_in)
{
    double u, v, w, angle;
    double inputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};
    double outputMatrix[4][1] = {0.0, 0.0, 0.0, 0.0};
    u = xVec_Rotation;
    v = yVec_Rotation;
    w = zVec_Rotation;
    angle =  (-angleRotation); // Rotating back to 0
    std::vector<double> inverse_rotated_coordinates;

    if (u > 0.001 || v > 0.001 || w > 0.001 )
    {
        double L = (u * u + v * v + w * w);
        angle = angle * M_PI / 180.0; //converting to radian value
        double u2 = u * u;
        double v2 = v * v;
        double w2 = w * w;

        rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
        rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
        rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
        rotationMatrix[0][3] = 0.0;

        rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
        rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
        rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
        rotationMatrix[1][3] = 0.0;

        rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
        rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
        rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
        rotationMatrix[2][3] = 0.0;

        rotationMatrix[3][0] = 0.0;
        rotationMatrix[3][1] = 0.0;
        rotationMatrix[3][2] = 0.0;
        rotationMatrix[3][3] = 1.0;

        inputMatrix[0][0] = x_in;
        inputMatrix[1][0] = y_in;
        inputMatrix[2][0] = z_in;
        inputMatrix[3][0] = 1.0;

            for(int i = 0; i < 4; i++ ){
                for(int j = 0; j < 1; j++){
                    outputMatrix[i][j] = 0;
                        for(int k = 0; k < 4; k++){
                            outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
                    }
                }
            }

        inverse_rotated_coordinates.push_back( round (outputMatrix[0][0]*1000.0) / 1000.0) ; // Rounding done since there is a rounding error for '0'
        inverse_rotated_coordinates.push_back( round (outputMatrix[1][0]*1000.0) / 1000.0) ;
        inverse_rotated_coordinates.push_back( round (outputMatrix[2][0]*1000.0) / 1000.0) ;
    }
    else{
        inverse_rotated_coordinates.push_back( x_in );
        inverse_rotated_coordinates.push_back( y_in );
        inverse_rotated_coordinates.push_back( z_in );
    }
    return inverse_rotated_coordinates;
}



