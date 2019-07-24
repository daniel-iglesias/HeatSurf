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

#include "window.h"
#include "geometry.h"


Window::Window() :
    plotTitle("Power density")
{
}

Window::~Window(){
}

void Window::addGeometry( Geometry* newGeometry )
{
    geometries.push_back( newGeometry );
}


void Window::drawGeometry()
{
    ren = vtkRenderer::New();
    ren->SetBackground ( 0.3, 0.3, 0.3 );
    for (auto it=geometries.begin(); it!=geometries.end(); ++it){
        dataSetMappers.push_back( vtkDataSetMapper::New() );
        dataSetMappers.back()->SetInputData ( (*it)->getGrid() );

        actors.push_back( vtkActor::New() );
        actors.back()->SetMapper ( dataSetMappers.back() );
        actors.back()->AddPosition ( 0, 0, 0 );
    //   anActor->GetProperty()->SetDiffuseColor(1.0, 0.3, 0.3);
    ren->AddActor ( actors.back() );
    }

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

 //this static function improves the appearance of the text edges
// since they are overlaid on a surface rendering of the cube's faces
    vtkMapper::SetResolveCoincidentTopologyToPolygonOffset();

// set up the widget
   widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
    widget->SetOrientationMarker ( axes );
    widget->SetInteractor ( iren );
    widget->SetViewport ( 0.0, 0.0, 0.3, 0.3 );
    widget->SetEnabled ( 1 );
    widget->InteractiveOff();
    widget->InteractiveOn();

}

void Window::drawScalar( bool )
{
    double minimum, maximum;
    table = vtkLookupTable::New();
    geometries[0]->calculateScalar();
    minimum = geometries[0]->getScalar()->GetRange()[0];
    maximum = geometries[0]->getScalar()->GetRange()[1];
    int i=0;
    for (auto it=geometries.begin(); it!=geometries.end(); ++it){
        (*it)->calculateScalar();

            if ((*it)->getScalar()->GetRange()[0] < minimum){
                minimum = (*it)->getScalar()->GetRange()[0];
                }

            if ((*it)->getScalar()->GetRange()[1] > maximum){
                maximum = (*it)->getScalar()->GetRange()[1];
                }
    }

    table->SetTableRange ( minimum , maximum ) ; // BUG: This should get the maximum range amongst all geometries
    table->SetHueRange ( 0.7, 0 );
//   table->SetSaturationRange (0, 1);
//    table->SetValueRange ( 1, 1 );
    table->SetNumberOfTableValues ( 100 );
    table->Build();

//    int i=0;
    for (auto it=geometries.begin(); it!=geometries.end(); ++it){
        dataSetMappers[i]->SetScalarRange ( minimum , maximum );
        dataSetMappers[i]->SetLookupTable ( table );
        i++;
    }

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


void Window::drawGrid()
{
//    grid3DDataSetMapper = vtkDataSetMapper::New();
//    grid3DDataSetMapper->SetInputData ( gridCubes );
//
//    grid3DActor = vtkActor::New();
//    grid3DActor->SetMapper ( grid3DDataSetMapper );
//    grid3DActor->AddPosition ( 0, 0, 0 );
////   grid3DActor->GetProperty()->SetDiffuseColor(1.0, 0.3, 0.3);
//
//    grid3Dren = vtkRenderer::New();
//    grid3Dren->SetBackground ( 0.3, 0.3, 0.3 );
//    grid3Dren->AddActor ( grid3DActor );
//
//    grid3Dren->ResetCamera();
//    grid3Dren->GetActiveCamera()->Azimuth ( 30 );
//    grid3Dren->GetActiveCamera()->Elevation ( 10 );
//    grid3Dren->ResetCameraClippingRange();
//
//    grid3DrenWin = vtkRenderWindow::New();
//    grid3DrenWin->AddRenderer ( grid3Dren );
//    grid3DrenWin->SetSize ( 800, 600 );
//
//    grid3Diren = vtkRenderWindowInteractor::New();
//    grid3Diren->SetRenderWindow ( grid3DrenWin );
//
//    grid3Dstyle = vtkInteractorStyleTrackballCamera::New();
//    grid3Diren->SetInteractorStyle ( grid3Dstyle );
//
//    // Axis...
//    grid3Daxes = vtkAxesActor::New();
//    grid3Dwidget = vtkOrientationMarkerWidget::New();
//
//    double a[3];
//
//    grid3Daxes->SetShaftTypeToCylinder();
//    grid3Daxes->SetXAxisLabelText ( "X" );
//    grid3Daxes->SetYAxisLabelText ( "Y" );
//    grid3Daxes->SetZAxisLabelText ( "Z" );
//    grid3Daxes->SetTotalLength ( 20.0, 20.0, 20.0 );
//    grid3Daxes->SetCylinderRadius ( 1.0 * grid3Daxes->GetCylinderRadius() );
//    grid3Daxes->SetConeRadius ( 0.7 * grid3Daxes->GetConeRadius() );
//    grid3Daxes->GetNormalizedTipLength ( a );
//    grid3Daxes->SetNormalizedTipLength ( 2.0*a[0],2.0*a[1],2.0*a[2] );
//
//    vtkTextProperty* grid3Dtprop = grid3Daxes->GetXAxisCaptionActor2D()->GetCaptionTextProperty();
//    grid3Dtprop->ItalicOn();
//    grid3Dtprop->ShadowOn();
//    grid3Dtprop->SetFontFamilyToTimes();
//
//    grid3Daxes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->ShallowCopy ( grid3Dtprop );
//    grid3Daxes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->ShallowCopy ( grid3Dtprop );
//
//// this static function improves the appearance of the text edges
//// since they are overlaid on a surface rendering of the cube's faces
//    vtkMapper::SetResolveCoincidentTopologyToPolygonOffset();
//
//// set up the widget
////   widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
//    grid3Dwidget->SetOrientationMarker ( grid3Daxes );
//    grid3Dwidget->SetInteractor ( grid3Diren );
//    grid3Dwidget->SetViewport ( 0.0, 0.0, 0.3, 0.3 );
//    grid3Dwidget->SetEnabled ( 1 );
//    grid3Dwidget->InteractiveOff();
//    grid3Dwidget->InteractiveOn();
//
}
//
void Window::drawGridScalar( bool )
{
//    //   # Define the value of a scalar field (in this case the scalar)
//    grid3Dscalar = vtkFloatArray::New();
//
//    for ( std::map<int,Node*>::iterator it = grid3Dnodes.begin();
//            it != grid3Dnodes.end();
//            ++it )
//    {
//        grid3Dscalar->InsertNextValue ( it->second->getDensity() );
//    }
//    gridCubes->GetPointData()->SetScalars ( grid3Dscalar );
//
//    grid3Dtable = vtkLookupTable::New();
//    grid3Dtable->SetTableRange ( grid3Dscalar->GetRange() );
//    grid3Dtable->SetHueRange ( 0.7, 0 );
////   table->SetSaturationRange (0, 1);
//    grid3Dtable->SetValueRange ( 1, 1 );
//    grid3Dtable->SetNumberOfTableValues ( 100 );
//    grid3Dtable->Build();
//
//    grid3DDataSetMapper->SetScalarRange ( grid3Dscalar->GetRange() );
//    grid3DDataSetMapper->SetLookupTable ( grid3Dtable );
////   aDataSetMapper->SelectColorArray( 2 );
//
//    grid3DbarActor = vtkScalarBarActor::New();
//    grid3DbarActor->SetLookupTable ( grid3Dtable );
//    grid3DbarActor->SetTitle ( plotTitle.c_str() );
//    grid3DbarActor->SetOrientationToVertical();
//    grid3DbarActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
//    grid3DbarActor->GetPositionCoordinate()->SetValue ( 0.8,0.1 );
//    grid3DbarActor->SetWidth ( 0.10 );
//    grid3DbarActor->SetHeight ( 0.8 );
//    grid3DbarActor->SetNumberOfLabels ( 9 );
//
//    plane = vtkPlane::New();
////  plane->SetOrigin( grid3Dscalar->GetCenter() );
//    plane->SetOrigin ( 0.,0.,0. );
//    plane->SetNormal ( 1., 0., 0. );
//    planeCut = vtkCutter::New();
//    planeCut->SetInputConnection ( grid3DDataSetMapper->GetOutputPort() );
//    planeCut->SetCutFunction ( plane );
//    cutMapper = vtkPolyDataMapper::New();
//    cutMapper->SetInputConnection ( planeCut->GetOutputPort() );
//    cutMapper->SetScalarRange ( grid3Dscalar->GetRange() );
//    cutActor = vtkActor::New();
//    cutActor->SetMapper ( cutMapper );
//
//
//    grid3Dren->AddActor2D ( grid3DbarActor );
//    grid3Dren->AddActor ( cutActor );
////   ren->GetRenderWindow()->Render();
//
//    grid3Diren->Initialize();
//    grid3Diren->Start();
}
