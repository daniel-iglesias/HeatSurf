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
#ifndef WINDOW_H
#define WINDOW_H

#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include <LMX/lmx.h>
#include <LMX/lmx_nlsolvers.h>

#include "particle.h"
#include "element.h"
#include "node.h"
#include "geometry.h"


class vtkRenderWindow;
class vtkInteractorStyleTrackballCamera;
class vtkRenderWindowInteractor;
class vtkRenderer;
class vtkDataSetMapper;
class vtkActor;
class vtkLookupTable;
class vtkScalarBarActor;

class vtkAxesActor;
class vtkOrientationMarkerWidget;

class vtkPlane;
class vtkCutter;
class vtkPolyDataMapper;


class Window{

public:

    Window();

    ~Window();

    void addGeometry( Geometry* );

    void drawGeometry( );

    void drawScalar( bool );

    void drawGrid( );

    void drawGridScalar( bool );


protected:

    std::string plotTitle;

    std::vector<Geometry*> geometries;
    std::vector<vtkDataSetMapper*> dataSetMappers;
    std::vector<vtkActor*> actors;

    vtkRenderWindow* renWin;
    vtkInteractorStyleTrackballCamera* style;
    vtkRenderWindowInteractor* iren;
    vtkRenderer* ren;
//    vtkDataSetMapper* aDataSetMapper;
//    vtkActor* anActor;

    vtkRenderWindow* grid3DrenWin;
    vtkInteractorStyleTrackballCamera* grid3Dstyle;
    vtkRenderWindowInteractor* grid3Diren;
    vtkRenderer* grid3Dren;
    vtkDataSetMapper* grid3DDataSetMapper;
    vtkActor* grid3DActor;

    vtkLookupTable* table;
    vtkScalarBarActor* barActor;

    vtkLookupTable* grid3Dtable;
    vtkScalarBarActor* grid3DbarActor;

    vtkAxesActor* axes;
    vtkOrientationMarkerWidget* widget;

    vtkAxesActor* grid3Daxes;
    vtkOrientationMarkerWidget* grid3Dwidget;

    vtkPlane* plane;
    vtkCutter* planeCut;
    vtkPolyDataMapper* cutMapper;
    vtkActor* cutActor;

};


#endif
