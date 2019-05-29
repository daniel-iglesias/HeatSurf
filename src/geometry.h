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
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include <LMX/lmx.h>
#include <LMX/lmx_nlsolvers.h>

#include "particle.h"
#include "element.h"
#include "node.h"

//////////////////////////////////////////// Doxygen file documentation entry:
/*!
  \file geometry.h

  \brief Base class for the different geometries.

  Abstract class, cannot be instantiated.

  \author Daniel Iglesias <daniel.iglesias@ciemat.es>

 */
//////////////////////////////////////////// Doxygen file documentation (end)


class vtkPoints;
class vtkFloatArray;
class vtkUnstructuredGrid;

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

/**
 @author Daniel Iglesias <daniel.iglesias@ciemat.es>
*/
class Geometry
{
public:
    Geometry();

    Geometry ( std::string, double, int, double );

    virtual ~Geometry();

    virtual void setSections ( double ) = 0;

    virtual void computeGeometry( ) = 0;

    void computeGrid3D( );

    virtual void computeIntersection ( Particle* ) = 0;

    virtual void computeNodalPower ( Particle* ) = 0;

    void computePowerDensity ( int );

    void drawGeometry( );

    void drawScalar( bool );

    void drawGrid( );

    void drawGridScalar( bool );

    virtual void outputPowerFile ( int ) = 0;

    virtual void outputPowerDensityFile( ) = 0;

    virtual void outputTable( ) = 0;

    std::string getType() { return type; }

    void outputAnsys3D();

    void outputVTKfiles();

    void outputBackscattering ( double );

protected:
    std::string type;
    int sections, sections3D;
    int sectors, sectors3D;
    double z0, length;
    double gridWidth, gridSpacing;
    std::vector< double > paramTrajectories;
    std::map<int, Node*> nodes;
    std::map<int, Node*> grid3Dnodes;
    std::vector< Element* > elements;
    std::vector< Element* > grid3D;
    std::string plotTitle;

    vtkPoints* gridPoints;
    vtkFloatArray* scalar;
    vtkUnstructuredGrid* grid;

    vtkPoints* grid3DPoints;
    vtkFloatArray* grid3Dscalar;
    vtkUnstructuredGrid* gridCubes;

    vtkRenderWindow* renWin;
    vtkInteractorStyleTrackballCamera* style;
    vtkRenderWindowInteractor* iren;
    vtkRenderer* ren;
    vtkDataSetMapper* aDataSetMapper;
    vtkActor* anActor;

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
