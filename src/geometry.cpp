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
                     double z0_in,
                     int sectors_in,
                     double length_in )
        : type ( type_in )
        , z0 ( z0_in )
        , sectors ( sectors_in )
        , length ( length_in )
{
}

Geometry::~Geometry()
{
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
    // Delta_y = Delta_x = 2*(radius/seectors)
    for ( i=0; i<sections3D; ++i )   // for z = i*sectionDif
    {
        z = i * sectionDif + z0;
        for ( j=0; j<=sectors3D; ++j )   // for y = radius*( 2*j/sectors - 1 )
        {
            y = gridWidth* ( 2.*j/sectors3D - 1. );
            for ( k=0; k<=sectors3D; ++k )   // for x = radius*( 2*k/sectors - 1 )
            {
                x = gridWidth* ( 2.*k/sectors3D - 1. );
                grid3Dnodes[i*sectNodes*sectNodes + j*sectNodes + k]
                = new Node ( x, y, z );
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

    for ( it_nodes = nodes.begin();
            it_nodes!= it_nodes_end;
            ++it_nodes )
    {
        x = it_nodes->second->getX();
        y = it_nodes->second->getY();
        z = it_nodes->second->getZ();
        i = floor ( ( z-z0 ) / sectionDif );
//          y = gridWidth*( 2.*j/sectors - 1. );
//          y / gridWidth +1 = ( 2.*j/sectors  );
        j = floor ( ( y/gridWidth + 1. ) * sectors3D/2. );
//          x / gridWidth = ( 2.*k/sectors - 1. );
        k = floor ( ( x/gridWidth + 1. ) * sectors3D/2. );

        elem_number = i*sectors3D*sectors3D + j*sectors3D + k;
//          cout << "x="<< x <<",y="<< y <<", z="<< z << endl;
//          cout << "i="<< i <<",j="<< j <<", k="<< k
//              <<"; elem_number="<< elem_number
//              <<"; elem_size="<< grid3D.size() << endl;
        if ( i < floor ( length/sectionDif ) )
            grid3D[elem_number]->addDensity3D ( it_nodes->second->getDensity() );
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


void Geometry::drawGeometry()
{
    aDataSetMapper = vtkDataSetMapper::New();
    aDataSetMapper->SetInputData ( grid );

    anActor = vtkActor::New();
    anActor->SetMapper ( aDataSetMapper );
    anActor->AddPosition ( 0, 0, 0 );
//   anActor->GetProperty()->SetDiffuseColor(1.0, 0.3, 0.3);

    ren = vtkRenderer::New();
    ren->SetBackground ( 0.3, 0.3, 0.3 );
    ren->AddActor ( anActor );

    ren->ResetCamera();
    ren->GetActiveCamera()->Azimuth ( 30 );
    ren->GetActiveCamera()->Elevation ( 10 );
    ren->ResetCameraClippingRange();

    renWin = vtkRenderWindow::New();
    renWin->AddRenderer ( ren );
    renWin->SetSize ( 800, 600 );

    iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow ( renWin );

    style = vtkInteractorStyleTrackballCamera::New();
    iren->SetInteractorStyle ( style );

    // Axis...
    axes = vtkAxesActor::New();
    widget = vtkOrientationMarkerWidget::New();

    double a[3];

    axes->SetShaftTypeToCylinder();
    axes->SetXAxisLabelText ( "X" );
    axes->SetYAxisLabelText ( "Y" );
    axes->SetZAxisLabelText ( "Z" );
    axes->SetTotalLength ( 20.0, 20.0, 20.0 );
    axes->SetCylinderRadius ( 1.0 * axes->GetCylinderRadius() );
    axes->SetConeRadius ( 0.7 * axes->GetConeRadius() );
    axes->GetNormalizedTipLength ( a );
    axes->SetNormalizedTipLength ( 2.0*a[0],2.0*a[1],2.0*a[2] );

    vtkTextProperty* tprop = axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty();
    tprop->ItalicOn();
    tprop->ShadowOn();
    tprop->SetFontFamilyToTimes();

    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->ShallowCopy ( tprop );
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->ShallowCopy ( tprop );

// this static function improves the appearance of the text edges
// since they are overlaid on a surface rendering of the cube's faces
    vtkMapper::SetResolveCoincidentTopologyToPolygonOffset();

// set up the widget
//   widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
    widget->SetOrientationMarker ( axes );
    widget->SetInteractor ( iren );
    widget->SetViewport ( 0.0, 0.0, 0.3, 0.3 );
    widget->SetEnabled ( 1 );
    widget->InteractiveOff();
    widget->InteractiveOn();

}

void Geometry::drawScalar( bool )
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

    table = vtkLookupTable::New();
    table->SetTableRange ( scalar->GetRange() );
    table->SetHueRange ( 0.7, 0 );
//   table->SetSaturationRange (0, 1);
    table->SetValueRange ( 1, 1 );
    table->SetNumberOfTableValues ( 100 );
    table->Build();

    aDataSetMapper->SetScalarRange ( scalar->GetRange() );
    aDataSetMapper->SetLookupTable ( table );
//   aDataSetMapper->SelectColorArray( 2 );

    barActor = vtkScalarBarActor::New();
    barActor->SetLookupTable ( table );
    barActor->SetTitle ( plotTitle.c_str() );
    barActor->SetOrientationToVertical();
    barActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    barActor->GetPositionCoordinate()->SetValue ( 0.8,0.1 );
    barActor->SetWidth ( 0.10 );
    barActor->SetHeight ( 0.8 );
    barActor->SetNumberOfLabels ( 9 );

    ren->AddActor2D ( barActor );
//   ren->GetRenderWindow()->Render();

    iren->Initialize();
    iren->Start();
}


void Geometry::drawGrid()
{
    grid3DDataSetMapper = vtkDataSetMapper::New();
    grid3DDataSetMapper->SetInputData ( gridCubes );

    grid3DActor = vtkActor::New();
    grid3DActor->SetMapper ( grid3DDataSetMapper );
    grid3DActor->AddPosition ( 0, 0, 0 );
//   grid3DActor->GetProperty()->SetDiffuseColor(1.0, 0.3, 0.3);

    grid3Dren = vtkRenderer::New();
    grid3Dren->SetBackground ( 0.3, 0.3, 0.3 );
    grid3Dren->AddActor ( grid3DActor );

    grid3Dren->ResetCamera();
    grid3Dren->GetActiveCamera()->Azimuth ( 30 );
    grid3Dren->GetActiveCamera()->Elevation ( 10 );
    grid3Dren->ResetCameraClippingRange();

    grid3DrenWin = vtkRenderWindow::New();
    grid3DrenWin->AddRenderer ( grid3Dren );
    grid3DrenWin->SetSize ( 800, 600 );

    grid3Diren = vtkRenderWindowInteractor::New();
    grid3Diren->SetRenderWindow ( grid3DrenWin );

    grid3Dstyle = vtkInteractorStyleTrackballCamera::New();
    grid3Diren->SetInteractorStyle ( grid3Dstyle );

    // Axis...
    grid3Daxes = vtkAxesActor::New();
    grid3Dwidget = vtkOrientationMarkerWidget::New();

    double a[3];

    grid3Daxes->SetShaftTypeToCylinder();
    grid3Daxes->SetXAxisLabelText ( "X" );
    grid3Daxes->SetYAxisLabelText ( "Y" );
    grid3Daxes->SetZAxisLabelText ( "Z" );
    grid3Daxes->SetTotalLength ( 20.0, 20.0, 20.0 );
    grid3Daxes->SetCylinderRadius ( 1.0 * grid3Daxes->GetCylinderRadius() );
    grid3Daxes->SetConeRadius ( 0.7 * grid3Daxes->GetConeRadius() );
    grid3Daxes->GetNormalizedTipLength ( a );
    grid3Daxes->SetNormalizedTipLength ( 2.0*a[0],2.0*a[1],2.0*a[2] );

    vtkTextProperty* grid3Dtprop = grid3Daxes->GetXAxisCaptionActor2D()->GetCaptionTextProperty();
    grid3Dtprop->ItalicOn();
    grid3Dtprop->ShadowOn();
    grid3Dtprop->SetFontFamilyToTimes();

    grid3Daxes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->ShallowCopy ( grid3Dtprop );
    grid3Daxes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->ShallowCopy ( grid3Dtprop );

// this static function improves the appearance of the text edges
// since they are overlaid on a surface rendering of the cube's faces
    vtkMapper::SetResolveCoincidentTopologyToPolygonOffset();

// set up the widget
//   widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
    grid3Dwidget->SetOrientationMarker ( grid3Daxes );
    grid3Dwidget->SetInteractor ( grid3Diren );
    grid3Dwidget->SetViewport ( 0.0, 0.0, 0.3, 0.3 );
    grid3Dwidget->SetEnabled ( 1 );
    grid3Dwidget->InteractiveOff();
    grid3Dwidget->InteractiveOn();

}

void Geometry::drawGridScalar( bool )
{
    //   # Define the value of a scalar field (in this case the scalar)
    grid3Dscalar = vtkFloatArray::New();

    for ( std::map<int,Node*>::iterator it = grid3Dnodes.begin();
            it != grid3Dnodes.end();
            ++it )
    {
        grid3Dscalar->InsertNextValue ( it->second->getDensity() );
    }
    gridCubes->GetPointData()->SetScalars ( grid3Dscalar );

    grid3Dtable = vtkLookupTable::New();
    grid3Dtable->SetTableRange ( grid3Dscalar->GetRange() );
    grid3Dtable->SetHueRange ( 0.7, 0 );
//   table->SetSaturationRange (0, 1);
    grid3Dtable->SetValueRange ( 1, 1 );
    grid3Dtable->SetNumberOfTableValues ( 100 );
    grid3Dtable->Build();

    grid3DDataSetMapper->SetScalarRange ( grid3Dscalar->GetRange() );
    grid3DDataSetMapper->SetLookupTable ( grid3Dtable );
//   aDataSetMapper->SelectColorArray( 2 );

    grid3DbarActor = vtkScalarBarActor::New();
    grid3DbarActor->SetLookupTable ( grid3Dtable );
    grid3DbarActor->SetTitle ( plotTitle.c_str() );
    grid3DbarActor->SetOrientationToVertical();
    grid3DbarActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    grid3DbarActor->GetPositionCoordinate()->SetValue ( 0.8,0.1 );
    grid3DbarActor->SetWidth ( 0.10 );
    grid3DbarActor->SetHeight ( 0.8 );
    grid3DbarActor->SetNumberOfLabels ( 9 );

    plane = vtkPlane::New();
//  plane->SetOrigin( grid3Dscalar->GetCenter() );
    plane->SetOrigin ( 0.,0.,0. );
    plane->SetNormal ( 1., 0., 0. );
    planeCut = vtkCutter::New();
    planeCut->SetInputConnection ( grid3DDataSetMapper->GetOutputPort() );
    planeCut->SetCutFunction ( plane );
    cutMapper = vtkPolyDataMapper::New();
    cutMapper->SetInputConnection ( planeCut->GetOutputPort() );
    cutMapper->SetScalarRange ( grid3Dscalar->GetRange() );
    cutActor = vtkActor::New();
    cutActor->SetMapper ( cutMapper );


    grid3Dren->AddActor2D ( grid3DbarActor );
    grid3Dren->AddActor ( cutActor );
//   ren->GetRenderWindow()->Render();

    grid3Diren->Initialize();
    grid3Diren->Start();
}


void Geometry::outputAnsys3D()
{
    std::ofstream outFile3D ( "Ansys_power_3D.dat" );
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
    std::ofstream outFile2D ( "Ansys_power_2D.dat" );
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
    vtkUnstructuredGridWriter *writer1 = vtkUnstructuredGridWriter::New();
    writer1->SetInputData ( grid );
    writer1->SetFileName ( "grid.vtk" );
    writer1->Write();

    vtkUnstructuredGridWriter *writer2 = vtkUnstructuredGridWriter::New();
    writer2->SetInputData ( gridCubes );
    writer2->SetFileName ( "gridCubes.vtk" );
    writer2->Write();

}

// Not valid for plates geometry!
void Geometry::outputBackscattering ( double precision )
{
    std::ofstream ofile ( "density.dat" );

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
