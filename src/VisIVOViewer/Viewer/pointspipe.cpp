/***************************************************************************
 *   Copyright (C) 2008 by Gabriella Caniglia *
 *  gabriella.caniglia@oact.inaf.it *
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
#include <cstdlib>
#include <cstring>

#include "pointspipe.h"

#include "visivoutils.h"
#include "luteditor.h"

#include "extendedglyph3d.h"

#include <sstream>
#include <algorithm>

#include "vtkSphereSource.h"
#include "vtkConeSource.h"
#include "vtkCylinderSource.h"
#include "vtkCubeSource.h"

#include "vtkCamera.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkLookupTable.h"

#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include"vtkGlyph3D.h"
#include "vtkScalarBarActor.h"
#include "vtkOutlineCornerFilter.h"
#include "vtkProperty.h"

#include "vtkGenericRenderWindowInteractor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "vtkAxesActor.h"
#include "vtkImageCast.h"
#include "vtkImageMathematics.h"
#include "vtkTreeCompositer.h"
#include "vtkCompositeRenderManager.h"
#include "vtkMPICommunicator.h"
#include "vtkWindowToImageFilter.h"
#include "vtkProcess.h"
#include "vtkCompositeZPass.h"
#include "vtkSequencePass.h"
#include "vtkSynchronizedRenderWindows.h"
#include "vtkMultiProcessController.h"
#include "vtkCubeAxesActor2D.h"

#include "vtkCompositedSynchronizedRenderers.h"

#include "vtkPNGWriter.h"
#include "vtkAppendPolyData.h"
#include "vtkImageIterator.h"

#include <mpi.h>
#define VTK_CREATE(type, name) vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

//---------------------------------------------------------------------
PointsPipe::PointsPipe ( VisIVOServerOptions options, vtkMPIController* contr, int n)
//---------------------------------------------------------------------
{
    m_visOpt=options;
    constructVTK();
    m_glyphFilter   = ExtendedGlyph3D::New();
    m_glyph         = vtkGlyph3D::New();
    m_pConeActor    = vtkActor::New();
    m_polyData      = vtkPolyData::New();
    m_pConeMapper   = vtkPolyDataMapper::New();
    num = n;
    m_pController   = contr;
}
//---------------------------------
PointsPipe::~PointsPipe()
//---------------------------------
{
    destroyVTK();
    if ( m_glyph!=0)
        m_glyph->Delete() ;
    if ( m_glyphFilter!=0)
        m_glyphFilter->Delete() ;
    if ( m_pConeMapper != 0 )
        m_pConeMapper->Delete();
    if ( m_pConeActor != 0 )
        m_pConeActor->Delete();
    if ( m_polyData!=0)
        m_polyData->Delete() ;
    
    
}
//---------------------------------
void PointsPipe::destroyAll()
//---------------------------------
{
    destroyVTK();
    if ( m_glyph!=0)
        m_glyph->Delete() ;
    if ( m_glyphFilter!=0)
        m_glyphFilter->Delete() ;
    if ( m_pConeMapper != 0 )
        m_pConeMapper->Delete();
    if ( m_pConeActor != 0 )
        m_pConeActor->Delete();
    if ( m_polyData!=0)
        m_polyData->Delete() ;
    
    
}

namespace
{
    
    class PointProcess : public vtkProcess{
        protected:

const double INVALID_CAM = -123456789.31;
            PointProcess(){
            }
        public:
            static PointProcess* New();
            
            void SetArgs(int anArgc, char* anArgv[]){
                this->Argc = anArgc;
                this->Argv = anArgv;
            }

            void SetOptions(VisIVOServerOptions visOpt){
                m_visOpt = visOpt;
            }

            void Execute(){
                std::clog<<"ENTER EXEC" << std::endl;
                int rank = this->Controller->GetLocalProcessId();
                int size = this->Controller->GetNumberOfProcesses();

                std::clog<<"RenderWin" << std::endl;
                renWin = vtkRenderWindow::New();
                renWin->DoubleBufferOn();
                renderer = vtkRenderer::New();
                
                std::clog<<"SyncWind" << std::endl;
                vtkSynchronizedRenderWindows* syncWindows = vtkSynchronizedRenderWindows::New();
                syncWindows->SetRenderWindow(renWin);
                syncWindows->SetParallelController(this->Controller);
                syncWindows->SetIdentifier(1);

                std::clog<<"CompositeRen" << std::endl;
                vtkCompositedSynchronizedRenderers* syncRenderers = vtkCompositedSynchronizedRenderers::New();
                syncRenderers->SetRenderer(renderer);
                syncRenderers->SetParallelController(this->Controller);
                m_lut           = vtkLookupTable::New();

                m_glyphFilter   = ExtendedGlyph3D::New();
                m_glyph         = vtkGlyph3D::New();
                m_pConeActor    = vtkActor::New();
                m_polyData      = vtkPolyData::New();
                m_pConeMapper   = vtkPolyDataMapper::New();
                
                std::ifstream inFile;
                
                vtkFloatArray *radiusArrays = vtkFloatArray::New();
                vtkFloatArray *xAxis = vtkFloatArray::New();
                vtkFloatArray *yAxis = vtkFloatArray::New();
                vtkFloatArray *zAxis = vtkFloatArray::New();
                
                xAxis->SetNumberOfTuples(m_visOpt.nRows);
                yAxis->SetNumberOfTuples(m_visOpt.nRows);
                zAxis->SetNumberOfTuples(m_visOpt.nRows);
                
                xAxis->SetName(m_visOpt.xField.c_str());
                yAxis->SetName(m_visOpt.yField.c_str());
                zAxis->SetName(m_visOpt.zField.c_str());
                
                std::clog<<"AfterAxis" << std::endl;
                int xIndex, yIndex, zIndex;
                int i;
                if (m_visOpt.dataRead && m_visOpt.goodAllocation)
                {
                    std::map<std::string, int>::iterator p;
                    for (p = m_visOpt.columns.begin(); p != m_visOpt.columns.end(); p++)
                    {
                        if (p->first == m_visOpt.xField)
                            xIndex = p->second;
                        if (p->first == m_visOpt.yField)
                            yIndex = p->second;
                        if (p->first == m_visOpt.zField)
                            zIndex = p->second;
                    }
                    
                    // Divide work among MPI processes
                    int startRow = rank * (m_visOpt.nRows / (size));
                    int endRow = (rank+1) * (m_visOpt.nRows / (size));
                    
                    for (i = startRow; i < endRow; i++)
                    {
                        xAxis->SetValue(i, m_visOpt.tableData[xIndex][i]);
                        yAxis->SetValue(i, m_visOpt.tableData[yIndex][i]);
                        zAxis->SetValue(i, m_visOpt.tableData[zIndex][i]);
                    }
                }
                else
                {
                    // Existing code for reading from file
                    inFile.open(m_visOpt.path.c_str(), std::ios::binary);
                    if (!inFile.is_open())
                        return;
                    
                    float tmpAxis[1];
                    
                    // Divide work among MPI processes
                    int startRow = rank * (m_visOpt.nRows / size);
                    int endRow = (rank + 1) * (m_visOpt.nRows / size);
                    if(rank>0){
                    for (i = startRow; i < endRow; i++)
                    {
                        inFile.seekg((m_visOpt.x * m_visOpt.nRows + i) * sizeof(float));
                        inFile.read((char *)(tmpAxis), sizeof(float));
                        if (m_visOpt.needSwap)
                            tmpAxis[0] = floatSwap((char *)(&tmpAxis[0]));
                        
                        xAxis->SetValue(i, tmpAxis[0]);
                        
                        inFile.seekg((m_visOpt.y * m_visOpt.nRows + i) * sizeof(float));
                        inFile.read((char *)(tmpAxis), sizeof(float));
                        if (m_visOpt.needSwap)
                            tmpAxis[0] = floatSwap((char *)(&tmpAxis[0]));
                        
                        yAxis->SetValue(i, tmpAxis[0]);
                        
                        inFile.seekg((m_visOpt.z * m_visOpt.nRows + i) * sizeof(float));
                        inFile.read((char *)(tmpAxis), sizeof(float));
                        if (m_visOpt.needSwap)
                            tmpAxis[0] = floatSwap((char *)(&tmpAxis[0]));
                        
                        zAxis->SetValue(i, tmpAxis[0]);
                    }
                    }
                    inFile.close();
                }
                //else
                
                xAxis->GetRange(m_xRange);  //!minimum and maximum value
                yAxis->GetRange(m_yRange);  //!minimum and maximum value
                zAxis->GetRange(m_zRange);  //!minimum and maximum value
                
                std::clog<<"AfterRange" << std::endl;
                SetXYZ(xAxis,yAxis,zAxis);
                
                m_polyData->SetPoints(m_points);
                
                m_points->Delete();
                
                
                
                // connect m_pRendererderer and m_pRendererder window and configure m_pRendererder window
                renWin->AddRenderer ( renderer );
                
                int nPoints = m_polyData->GetNumberOfPoints();
                //  // std::clog<<nPoints<<std::endl;
                
                vtkCellArray *newVerts = vtkCellArray::New();
                newVerts->EstimateSize (nPoints,1 );
                //  newVerts->InsertNextCell ( nPoints );
                
                for ( int i = 0; i < nPoints; i++ )
                {
                    newVerts->InsertNextCell(1);
                    newVerts->InsertCellPoint ( i );
                }
                m_polyData->SetVerts ( newVerts );
                
                xAxis->Delete();
                yAxis->Delete();
                zAxis->Delete();
                
                float tmp[1];
                
                if(m_visOpt.radiusscalar!="none"&& m_visOpt.scaleGlyphs!="none"&& m_visOpt.nGlyphs!=0 )
                {
                    if(m_visOpt.color=="none")
                    {
                        m_visOpt.color="yes";
                        m_visOpt.nColorTable=22; //default is whyte
                        
                        if(m_visOpt.oneColor=="yellow") m_visOpt.nColorTable=19;
                        if(m_visOpt.oneColor=="red") m_visOpt.nColorTable=24;
                        if(m_visOpt.oneColor=="green") m_visOpt.nColorTable=25;
                        if(m_visOpt.oneColor=="blue") m_visOpt.nColorTable=26;
                        if(m_visOpt.oneColor=="cyan") m_visOpt.nColorTable=20;
                        if(m_visOpt.oneColor=="violet") m_visOpt.nColorTable=21;
                        if(m_visOpt.oneColor=="black") m_visOpt.nColorTable=23;
                        
                        m_visOpt.colorScalar=m_visOpt.radiusscalar;
                        m_visOpt.showLut=false;
                    }
                    radiusArrays->SetNumberOfTuples(m_visOpt.nRows);
                    if(m_visOpt.dataRead && m_visOpt.goodAllocation)
                    {
                        int iRadius;
                        std::map<std::string, int>::iterator p;
                        for(p=m_visOpt.columns.begin();p!=m_visOpt.columns.end();p++)
                            if(p->first==m_visOpt.radiusscalar) iRadius=p->second;
                        for(i=0; i<m_visOpt.nRows;i++)
                            radiusArrays->  SetValue(i,m_visOpt.tableData[iRadius][i]);
                    } else
                    {
                        for (i=0;i<m_visOpt.nRows;i++)
                        {
                            inFile.seekg((m_visOpt.nRadius * m_visOpt.nRows+i )* sizeof(float));
                            inFile.read((char *)(tmp ),  sizeof(float));
                            
                            if(m_visOpt.needSwap)
                                tmp[0]=floatSwap((char *)(&tmp[0]));
                            
                            radiusArrays->  SetValue(i,tmp[0]);
                        }
                    }//else
                    
                    
                    /*    for (i=0;i<m_visOpt.nRows;i++)
                    {
                    inFile.seekg((m_visOpt.nRadius * m_visOpt.nRows+i )* sizeof(float));
                    inFile.read((char *)(tmp ),  sizeof(float));
                    
                    if(m_visOpt.needSwap)
                    tmp[0]=floatSwap((char *)(&tmp[1]));
                    
                    radiusArrays->  SetValue(i,tmp[0]);
                    }*/
                    radiusArrays->SetName(m_visOpt.radiusscalar.c_str());
                    m_polyData->GetPointData()->SetScalars(radiusArrays);   //!adda scalar to vtkpolidata
                    //radiusArrays->Delete();
                    
                }
                
                
                if(m_visOpt.heightscalar!="none"&& m_visOpt.scaleGlyphs!="none"&& m_visOpt.nGlyphs!=0 && m_visOpt.nGlyphs!=1 )
                {
                    vtkFloatArray *heightArrays = vtkFloatArray::New();
                    heightArrays->SetNumberOfTuples(m_visOpt.nRows);
                    
                    if(m_visOpt.dataRead && m_visOpt.goodAllocation)
                    {
                        int iheight;
                        std::map<std::string, int>::iterator p;
                        for(p=m_visOpt.columns.begin();p!=m_visOpt.columns.end();p++)
                            if(p->first==m_visOpt.heightscalar) iheight=p->second;
                        for(i=0; i<m_visOpt.nRows;i++)
                            heightArrays->  SetValue(i,m_visOpt.tableData[iheight][i]);
                    } else
                    {
                        for (i=0;i<m_visOpt.nRows;i++)
                        {
                            inFile.seekg((m_visOpt.nHeight * m_visOpt.nRows+i )* sizeof(float));
                            inFile.read((char *)(tmp ),  sizeof(float));
                            
                            if(m_visOpt.needSwap)
                                tmp[0]=floatSwap((char *)(&tmp[0]));
                            
                            heightArrays->  SetValue(i,tmp[0]);
                        }
                    }//else
                    
                    heightArrays->SetName(m_visOpt.heightscalar.c_str());
                    m_polyData->GetPointData()->SetScalars(heightArrays);
                    heightArrays->Delete();
                }
                
                if(m_visOpt.colorScalar!="none")
                {
                    vtkFloatArray *lutArrays =vtkFloatArray::New();
                    lutArrays->SetNumberOfTuples(m_visOpt.nRows);
                    
                    if(m_visOpt.dataRead && m_visOpt.goodAllocation)
                    {
                        int ilut;
                        std::map<std::string, int>::iterator p;
                        for(p=m_visOpt.columns.begin();p!=m_visOpt.columns.end();p++)
                            if(p->first==m_visOpt.colorScalar) ilut=p->second;
                        for(i=0; i<m_visOpt.nRows;i++)
                            lutArrays->  SetValue(i,m_visOpt.tableData[ilut][i]);
                        
                    } else
                    {
                        for (i=0;i<m_visOpt.nRows;i++)
                        {
                            inFile.seekg((m_visOpt.nColorScalar * m_visOpt.nRows+i )* sizeof(float));
                            inFile.read((char *)(tmp ),  sizeof(float));
                            
                            if(m_visOpt.needSwap)
                                tmp[0]=floatSwap((char *)(&tmp[0]));
                            lutArrays->  SetValue(i,tmp[0]);
                        }
                    } //else
                    
                    lutArrays->SetName(m_visOpt.colorScalar.c_str());
                    m_polyData->GetPointData()->SetScalars(lutArrays);
                    
                    double range [2];
                    lutArrays->GetRange(range);
                    if(range[0]<=0)
                        m_visOpt.uselogscale="none";
                    
                    lutArrays->Delete();
                }
                
                inFile.close();
                    
                setBoundingBox (m_polyData);
                
                m_pConeMapper->SetInputData (m_polyData );
                //m_pConeMapper->SetScalarModeToUseCellFieldData();
                m_pConeMapper->SetScalarRange(0,size-1);
                m_pConeMapper->SetPiece(rank);
                m_pConeMapper->SetNumberOfPieces(size);
                m_pConeMapper->Update();
                m_pConeActor->SetMapper ( m_pConeMapper );
                
                if(m_visOpt.color=="none")
                {
                    vtkProperty *P = vtkProperty::New();
                    P->SetColor(1, 1 ,1);  //default is white
                    if(m_visOpt.oneColor=="yellow") P->SetColor(1, 1 ,0);
                    if(m_visOpt.oneColor=="red") P->SetColor(1, 0 ,0);
                    if(m_visOpt.oneColor=="green") P->SetColor(0, 1 ,0);
                    if(m_visOpt.oneColor=="blu") P->SetColor(0, 0 ,1);
                    if(m_visOpt.oneColor=="cyane") P->SetColor(0, 1 ,1);
                    if(m_visOpt.oneColor=="violet") P->SetColor(1, 0 ,1);
                    if(m_visOpt.oneColor=="black") P->SetColor(0, 0 ,0);
                    // 110 (giallo);  100(rosso), 010 Verde, 001 Blu, 011 Cyane, 101 Violet, 111 bianco ,000 (nero)
                    m_pConeActor->SetProperty(P);
                    P->Delete();
                }
                
                renderer->AddActor ( m_pConeActor );
                
                if (m_visOpt.opacity<0 )
                    m_visOpt.opacity=0;
                
                else if( m_visOpt.opacity>1)
                    m_visOpt.opacity=1;
                
                m_pConeActor->GetProperty()->SetOpacity ( m_visOpt.opacity);
                
                if (m_visOpt.nGlyphs!=0)
                    setGlyphs (  );
                
                if(m_visOpt.scaleGlyphs!="none")
                    setScaling ();
                
                
                if ( m_visOpt.colorScalar!="none" && m_visOpt.color!="none")
                    setLookupTable ( );
                
                if(m_visOpt.scaleGlyphs!="none"&& m_visOpt.nGlyphs!=0 ||( (m_visOpt.heightscalar!="none" ||  (m_visOpt.radiusscalar!="none" && m_visOpt.nGlyphs!=1))  ))
                    m_glyph->ScalingOn();
                else
                    m_glyph->ScalingOff();
                
                
                //renderer->SetBackground ( 0.0,0.0,0.0 );
                if(m_visOpt.backColor=="yellow") renderer->SetBackground (1, 1 ,0);
                if(m_visOpt.backColor=="red")renderer->SetBackground (1, 0 ,0);
                if(m_visOpt.backColor=="green") renderer->SetBackground (0, 1 ,0);
                if(m_visOpt.backColor=="blue") renderer->SetBackground (0, 0 ,1);
                if(m_visOpt.backColor=="cyan") renderer->SetBackground (0, 1 ,1);
                if(m_visOpt.backColor=="violet") renderer->SetBackground (1, 0 ,1);
                if(m_visOpt.backColor=="white") renderer->SetBackground (1, 1 ,1);
                
                if(m_visOpt.imageSize=="small")
                    renWin->SetSize ( 512,365 );
                else if(m_visOpt.imageSize=="large")
                    renWin->SetSize ( 1024,731 );
                else
                    renWin->SetSize ( 792,566 );
                
                renWin->SetWindowName ("VisIVOServer View");
                
                if(m_visOpt.stereo)
                {
                    renWin->StereoRenderOn();
                    if(m_visOpt.stereoMode=="RedBlue")
                    {
                        renWin->SetStereoTypeToRedBlue();;
                    }
                    else if(m_visOpt.stereoMode=="Anaglyph")
                    {
                        renWin->SetStereoTypeToAnaglyph();
                        renWin->SetAnaglyphColorSaturation(m_visOpt.anaglyphsat);
                        m_visOpt.anaglyphmask=trim(m_visOpt.anaglyphmask);
                        std::stringstream anatmp(m_visOpt.anaglyphmask);
                        int rightColorMask=4, leftColorMask=3;
                        anatmp>>rightColorMask;
                        anatmp>>leftColorMask;
                        renWin->SetAnaglyphColorMask(rightColorMask,leftColorMask);
                    }
                    else
                    {
                        renWin->SetStereoTypeToCrystalEyes();
                        
                        if(m_visOpt.stereoImg==0)
                        {
                            renWin->SetStereoTypeToRight();
                        }else
                        {
                            renWin->SetStereoTypeToLeft();
                        }
                    }
                    renWin->StereoUpdate();
                }
                int retVal = 1;

                            if(rank == 0){
                                vtkRenderWindowInteractor* iren = vtkGenericRenderWindowInteractor::New();
                                iren->SetRenderWindow(renWin);
                                iren->Initialize();
                                iren->Start();

                            std::clog<<"Rank 0 before TriggerRMI" << std::endl;
                                this->Controller->TriggerBreakRMIs();
                            std::clog<<"Rank 0 after TriggerRMI" << std::endl;
                                for (int i = 1; i < size; i++)
                                {
                                this->Controller->Send(&retVal, 1, i, 33);
                                }
                                setCamera (NULL);
                                double *bounds;
                                bounds=new double[6];
                                bounds[0]=m_xRange[0];
                                bounds[1]=m_xRange[1];
                                bounds[2]=m_yRange[0];
                                bounds[3]=m_yRange[1];
                                bounds[4]=m_zRange[0];
                                bounds[5]=m_zRange[1];

                                    if(m_visOpt.showAxes) setAxes (m_polyData,bounds );
                                delete [] bounds;
                                iren->ExitEvent();
                            iren->Delete();
                            saveImageAsPng(0);
                            
                            }
                            else{
                            std::clog<<"Rank 1 before after processRMI" << std::endl;
                                this->Controller->ProcessRMIs();

                            std::clog<<"Rank 1 after processRMI" << std::endl;
                                this->Controller->Receive(&retVal, 1, 0, 33);
                            }
                //std::clog<< "Process " << this->Controller->GetLocalProcessId() << " exiting rendering " << std::endl;
                renderer->Delete();
                renWin->Delete();
                syncWindows->Delete();
                syncRenderers->Delete();
            }

void setAxes ( vtkDataSet *data, double *bounds  )
//---------------------------------------------------------------------
{
  vtkCubeAxesActor2D* axesActor=vtkCubeAxesActor2D::New();
  	/* VTK9 migration
  axesActor->SetInput ( data );
		replaced
  axesActor->SetInputData ( data );
	*/
  
  axesActor->SetInputData ( data );
  axesActor->UseRangesOn();
   
  axesActor->SetBounds ( bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5]);
  axesActor->SetRanges (  bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5] );

  axesActor->SetViewProp ( NULL );
  axesActor->SetScaling ( 0 );
  axesActor->SetPickable ( 0 );
  axesActor->SetCamera ( renderer->GetActiveCamera() );
  axesActor->SetCornerOffset ( 0.1 );
  axesActor->SetLabelFormat ( "%6.5g" );
  axesActor->SetInertia ( 100 );
  axesActor->SetFlyModeToOuterEdges();
  axesActor->SetVisibility ( true );

  axesActor->SetXLabel (m_visOpt.xField.c_str() );
  axesActor->SetYLabel (m_visOpt.yField.c_str() );
  axesActor->SetZLabel ( m_visOpt.zField.c_str() );

  axesActor->Modified();

  renderer->AddActor2D ( axesActor );
  
  if(axesActor!=0)
    axesActor->Delete();
  
}
    //---------------------------------------------------------------------
void setCamera (SplotchCamera *splCamera)
//---------------------------------------------------------------------
{

    
    
	m_camera =renderer->GetActiveCamera();
	

//	m_camera->SetClippingRange(m_visOpt.cliprange[0],m_visOpt.cliprange[1]);
	
	if(m_visOpt.setCameraPos)
	  if(m_visOpt.cameraPos[0]!=INVALID_CAM && m_visOpt.cameraPos[1]!=INVALID_CAM && m_visOpt.cameraPos[2]!=INVALID_CAM)
	      m_camera->SetPosition(m_visOpt.cameraPos[0],m_visOpt.cameraPos[1],m_visOpt.cameraPos[2]);
	  else if(m_visOpt.cameraPosPrev[0]!=INVALID_CAM && m_visOpt.cameraPosPrev[1]!=INVALID_CAM &&
	          m_visOpt.cameraPosPrev[2]!=INVALID_CAM)
	      m_camera->SetPosition(m_visOpt.cameraPosPrev[0],m_visOpt.cameraPosPrev[1],m_visOpt.cameraPosPrev[2]);
	     
	  
	if(m_visOpt.setCameraFocalPoint)
	  if(m_visOpt.cameraFocalPoint[0]!=INVALID_CAM && m_visOpt.cameraFocalPoint[1]!=INVALID_CAM && m_visOpt.cameraFocalPoint[2]!=INVALID_CAM)
	    m_camera->SetFocalPoint(m_visOpt.cameraFocalPoint[0],m_visOpt.cameraFocalPoint[1],m_visOpt.cameraFocalPoint[2]);
	  else if(m_visOpt.cameraFocalPointPrev[0]!=INVALID_CAM && m_visOpt.cameraFocalPointPrev[1]!=INVALID_CAM &&
	          m_visOpt.cameraFocalPointPrev[2]!=INVALID_CAM)
	      m_camera->SetFocalPoint(m_visOpt.cameraFocalPointPrev[0],m_visOpt.cameraFocalPointPrev[1],m_visOpt.cameraFocalPointPrev[2]);

    
    
    
	if(m_visOpt.setCameraRoll)
	  if(m_visOpt.cameraRoll!=INVALID_CAM) 
	      m_camera->SetRoll(m_visOpt.cameraRoll);
	  else if(m_visOpt.cameraRollPrev[0]!=INVALID_CAM)
	       m_camera->SetRoll(m_visOpt.cameraRollPrev[0]);
    
    if(m_visOpt.fovIsGiven)
        m_camera->SetViewAngle(m_visOpt.fov);


	
	
	if(m_visOpt.numImage ==0)  // QUI con Splotch è sbagliato
	{
		m_camera->Azimuth(0);
		m_camera->Elevation(0);
		m_camera->Zoom ( 1 );
	}
	else if(m_visOpt.numImage==1)
	{
		m_camera->Azimuth(90);
		m_camera->Elevation(0);
		m_camera->Zoom (1 );
	}
	else if(m_visOpt.numImage==2)
	{
		m_camera->Azimuth(0);
		m_camera->Elevation(89.90);
		m_camera->Zoom (1 );
	}
	else if(m_visOpt.numImage==3)
	{
		m_camera->Azimuth(45);
		m_camera->Elevation(45);
		m_camera->Zoom (1 );
	}
	else if(m_visOpt.numImage==4)
	{
		if(m_visOpt.elevation==90.0)
			m_visOpt.elevation=89.90;
		if(m_visOpt.elevation==-90.0)
			m_visOpt.elevation=-89.90;
		if(m_visOpt.elevation>90)
		{
			std::cerr<<"Invalid elevation "<<m_visOpt.elevation<<std::endl;	
			m_visOpt.elevation=89.90;
		}
		if(m_visOpt.elevation<-90)
		{
			std::cerr<<"Invalid elevation "<<m_visOpt.elevation<<std::endl;	
			m_visOpt.elevation=-89.90;
		}
 
		m_camera->Azimuth(m_visOpt.azimuth);
		m_camera->Elevation(m_visOpt.elevation);
		m_camera->Zoom ( m_visOpt.zoom ); 

	}
//    double d[2];
//    m_camera->GetClippingRange(d);
//    std::clog<<d[0]<<" "<<d[1]<<std::endl;
    
       if(m_visOpt.clipset) m_camera->SetClippingRange(m_visOpt.cliprange[0],m_visOpt.cliprange[1]);
/*    
	if(m_visOpt.setCameraPos) {
        if(m_visOpt.cameraPos[0]!=INVALID_CAM && m_visOpt.cameraPos[1]!=INVALID_CAM && m_visOpt.cameraPos[2]!=INVALID_CAM) {
            m_camera->SetPosition(m_visOpt.cameraPos[0],m_visOpt.cameraPos[1],m_visOpt.cameraPos[2]);
        }
    } else if(m_visOpt.cameraPosPrev[0]!=INVALID_CAM && m_visOpt.cameraPosPrev[1]!=INVALID_CAM &&
	          m_visOpt.cameraPosPrev[2]!=INVALID_CAM)
    {
        m_camera->SetPosition(m_visOpt.cameraPosPrev[0],m_visOpt.cameraPosPrev[1],m_visOpt.cameraPosPrev[2]);
	}
    
	if(m_visOpt.setCameraFocalPoint) {
        if(m_visOpt.cameraFocalPoint[0]!=INVALID_CAM && m_visOpt.cameraFocalPoint[1]!=INVALID_CAM && m_visOpt.cameraFocalPoint[2]!=INVALID_CAM) {
            m_camera->SetFocalPoint(m_visOpt.cameraFocalPoint[0],m_visOpt.cameraFocalPoint[1],m_visOpt.cameraFocalPoint[2]);
        }
    } else if(m_visOpt.cameraFocalPointPrev[0]!=INVALID_CAM && m_visOpt.cameraFocalPointPrev[1]!=INVALID_CAM &&
	          m_visOpt.cameraFocalPointPrev[2]!=INVALID_CAM)
    {
        m_camera->SetFocalPoint(m_visOpt.cameraFocalPointPrev[0],m_visOpt.cameraFocalPointPrev[1],m_visOpt.cameraFocalPointPrev[2]);
    }
    
	if(m_visOpt.setCameraRoll) {
        if(m_visOpt.cameraRoll!=INVALID_CAM) {
            m_camera->SetRoll(m_visOpt.cameraRoll);
        }
    } else if(m_visOpt.cameraRollPrev[0]!=INVALID_CAM) {
        m_camera->SetRoll(m_visOpt.cameraRollPrev[0]);
    }
	  
*/	  
	double *gettmp;
	
	gettmp= m_camera->GetPosition();
	m_visOpt.cameraPosPrev[0]=gettmp[0];
	m_visOpt.cameraPosPrev[1]=gettmp[1];
	m_visOpt.cameraPosPrev[2]=gettmp[2];
    
    
	
	gettmp= m_camera->GetFocalPoint();
	m_visOpt.cameraFocalPointPrev[0]=gettmp[0];
	m_visOpt.cameraFocalPointPrev[1]=gettmp[1];
	m_visOpt.cameraFocalPointPrev[2]=gettmp[2];

	m_visOpt.cameraRollPrev[0]=m_camera->GetRoll();


    
    
	// m_camera->SetPosition(32+1*10,32+5*10,32+70);
//  m_camera->SetViewUp(1,0,0);
// m_camera->SetViewAngle(160);
// vtkIndent indent;
// std::ofstream planeof("cam.txt");
//m_camera->PrintSelf(planeof,indent);
// planeof.close();

		if(splCamera!=NULL)
		{
		vtkIndent indent;
/*		std::ofstream cameraof("camera.txt");
		m_camera->PrintSelf(cameraof,indent);*/
		double a,b,c;
		m_camera->GetPosition(a,b,c);
		splCamera->position[0]=a;
		splCamera->position[1]=b;
		splCamera->position[2]=c;
            
		m_camera->GetFocalPoint(a,b,c);
		splCamera->lookat[0]=a;
		splCamera->lookat[1]=b;
		splCamera->lookat[2]=c;
        
        
		//splCamera->roll=m_camera->GetRoll()+45;
            
          //  std::cout<<"roll: "<<splCamera->roll<<std::endl;

            m_camera->GetViewUp(a,b,c);
            
           // m_camera->SetParallelScale(400);
            
            splCamera->sky[0]=a;
            splCamera->sky[1]=b;
            splCamera->sky[2]=c;

            for(int i=0;i<3;i++)
		  if(splCamera->lookat[i]==splCamera->position[i]) splCamera->position[i]=splCamera->position[i]*1.01;
		splCamera->fov=m_camera->GetViewAngle();
		m_camera->GetClippingRange(a,b);
		splCamera->clip[0]=a;
		splCamera->clip[1]=b;
 
/*		m_camera->SetPosition(80,80,280);
		m_camera->Roll(180);
//		m_camera->SetFocalPoint(-10,-10,-10);
		m_camera->PrintSelf(cameraof,indent);*/
//		cameraof.close();
		}
  
	return;
}

void setBoundingBox ( vtkDataObject *data )
//---------------------------------------------------------------------
{
	vtkOutlineCornerFilter*corner= vtkOutlineCornerFilter::New();
	/* VTK9 migration
	corner->SetInput(data);
		replaced
	corner->SetInputData(data);
		*/
	corner->SetInputData(data);
	corner->ReleaseDataFlagOn();

	vtkPolyDataMapper *outlineMapper = vtkPolyDataMapper::New();
		/* VTK9 migration
	outlineMapper->SetInput ( corner->GetOutput() );
		replaced
	outlineMapper->SetInputData ( corner->GetOutput() );
		*/
	outlineMapper->SetInputConnection(0, corner->GetOutputPort(0));

	vtkProperty *outlineProperty = vtkProperty::New();
	outlineProperty->SetColor ( 1.0,1.0,1.0 ); // Set the color to white
	outlineProperty->SetAmbient ( 1 );
	if(m_visOpt.showBox)
		outlineProperty->SetOpacity ( 0.999 );
	else
        outlineProperty->SetOpacity ( 0.0 );

	outlineProperty->SetRepresentationToWireframe();
	outlineProperty->SetInterpolationToFlat();

	vtkActor*outlineActor = vtkActor::New();
	outlineActor->SetMapper ( outlineMapper );
	outlineActor->SetProperty ( outlineProperty );
	outlineActor->SetPickable ( false );
	outlineActor->SetVisibility ( true );

	renderer->AddActor ( outlineActor );
  
	if (outlineActor!=0)
		outlineActor->Delete();
	if (outlineMapper!=0)
		outlineMapper->Delete();
	if (outlineProperty!=0)
		outlineProperty->Delete();
	if (corner!=0)
		corner->Delete();
}
            //-------------------------------------------------------------------------
bool SetXYZ(vtkFloatArray *xField, vtkFloatArray *yField, vtkFloatArray *zField  )
//-------------------------------------------------------------------------
{
    double scalingFactors[3];
    scalingFactors[0]=scalingFactors[1]=scalingFactors[2]=0;
    
    m_points=vtkPoints::New();
    m_points->SetNumberOfPoints(m_visOpt.nRows);
    
    
    if(xField->GetNumberOfComponents() != yField->GetNumberOfComponents())
    {
        if(zField && (xField->GetNumberOfComponents() != zField->GetNumberOfComponents() \
                      || yField->GetNumberOfComponents() != zField->GetNumberOfComponents()))
        {
            false;
        }
        false; // component mismatch, do nothing
    }
    
    
    if(m_visOpt.scale=="yes")
    {
        
        double size = 0;
        
        
        size = (m_xRange[1] - m_xRange[0] != 0 ? m_xRange[1] - m_xRange[0] : m_xRange[1]);
        scalingFactors[0] = size * 0.1;
        
        
        
        size = (m_yRange[1] - m_yRange[0] != 0 ? m_yRange[1] - m_yRange[0] : m_yRange[1]);
        scalingFactors[1] = size * 0.1;
        
        
        size = (m_zRange[1] - m_zRange[0] != 0 ? m_zRange[1] - m_zRange[0] : m_zRange[1]);
        scalingFactors[2] = size * 0.1;
    }
    
    double scalingFactorsInv[3];
    
    int i = 0;
    for(i = 0; i < 3; i++)
        scalingFactorsInv[i] = ((scalingFactors && scalingFactors[i] != 0) ? 1/scalingFactors[i] : 0);
    
    // Set the points data
    if(m_visOpt.scale=="yes")
    {
        for(i = 0; i < m_visOpt.nRows; i++)
        {
            float inPoint[3];
            float outPoint[3];
            inPoint[0] = outPoint[0] = xField->GetValue(i) * scalingFactorsInv[0];
            inPoint[1] = outPoint[1] = yField->GetValue(i) * scalingFactorsInv[1];
            inPoint[2] = outPoint[2] = zField->GetValue(i) * scalingFactorsInv[2];
            
            m_points->SetPoint(i,outPoint);
        }
    }
    else
        for(i = 0; i < m_visOpt.nRows; i++)
        {
            float outPoint[3];
            
            outPoint[0] = xField->GetValue(i) ;
            outPoint[1] = yField->GetValue(i) ;
            outPoint[2] = zField->GetValue(i) ;
            
            m_points->SetPoint(i,outPoint);
        }
    
    return true;
}

//---------------------------------------------------------------------
void setScaling ()
//---------------------------------------------------------------------
{
    m_glyphFilter->SetUseSecondScalar(true);
    m_glyphFilter->SetUseThirdScalar(true);
    
    m_glyphFilter->SetScaling(1);
    
    if( m_visOpt.heightscalar!="none" && m_visOpt.scaleGlyphs!="none" && m_visOpt.nGlyphs!=0 && m_visOpt.nGlyphs!=1)
        m_glyphFilter->SetInputScalarsSelectionY(m_visOpt.heightscalar.c_str());
    
    if( m_visOpt.radiusscalar!="none" && m_visOpt.scaleGlyphs!="none" && m_visOpt.nGlyphs!=0)
        m_glyphFilter->SetInputScalarsSelectionXZ(m_visOpt.heightscalar.c_str());
    
    
    if( m_visOpt.nGlyphs!=0)
        m_glyphFilter->SetScaleModeToScaleByScalar();
    else
        m_glyphFilter->ScalarVisibilityOff();
    
    
}
//---------------------------------------------------------------------
void setGlyphs ( )
//---------------------------------------------------------------------
{
    int max=1000;
    
    if ( m_visOpt.nRows<max )
    {
        /* VTK9 migration
         m_glyph->SetInput (m_polyData );
         replaced
         m_glyph->SetInputData (m_polyData );
         
         */
        m_glyph->SetInputData (m_polyData );
        
        
        if (m_visOpt.scale=="yes")
            m_glyph->SetScaleFactor ( 0.04 );
        
        else
            m_glyph->SetScaleFactor ( 2.5 );
        
        m_pConeMapper->SetInputConnection( m_glyph->GetOutputPort() );
        
        
        if (m_visOpt.nGlyphs==1)
        {
            m_sphere   = vtkSphereSource::New();
            setResolution ( );
            setRadius ();
            /* VTK9 migration
             m_glyph->SetSource ( m_sphere->GetOutput() );
             replaced
             m_glyph->SetSourceData ( m_sphere->GetOutput() );
             */
            m_glyph->SetSourceData ( m_sphere->GetOutput() );
            m_sphere->Delete();
        }
        
        else if (m_visOpt.nGlyphs==2)
        {
            m_cone   = vtkConeSource::New();
            setResolution ( );
            setRadius ();
            /* VTK9 migration
             m_glyph->SetSource ( m_cone->GetOutput() );
             replaced
             m_glyph->SetSourceData ( m_cone->GetOutput() );
             */
            m_glyph->SetSourceData ( m_cone->GetOutput() );
            m_cone->Delete();
        }
        
        else if (m_visOpt.nGlyphs==3)
        {
            m_cylinder   = vtkCylinderSource::New();
            setResolution ( );
            setRadius ();
            /* VTK9 migration
             m_glyph->SetSource ( m_cylinder->GetOutput() );
             replaced
             m_glyph->SetSourceData ( m_cylinder->GetOutput() );
             */
            m_glyph->SetSourceData ( m_cylinder->GetOutput() );
            m_cylinder->Delete();
        }
        
        else if (m_visOpt.nGlyphs==4)
        {
            m_cube   = vtkCubeSource::New();
            setRadius ();
            /* VTK9 migration
             m_glyph->SetSource ( m_cube->GetOutput() );
             replaced
             m_glyph->SetSourceData ( m_cube->GetOutput() );
             */
            m_glyph->SetSourceData ( m_cube->GetOutput() );
            m_cube->Delete();
        }
        
        
    }
    return ;
}




//---------------------------------------------------------------------
void setLookupTable ()
//---------------------------------------------------------------------
{
    
    double b[2];
    m_polyData->GetPointData()->SetActiveScalars(m_visOpt.colorScalar.c_str());
    
    m_polyData->GetPointData()->GetScalars(m_visOpt.colorScalar.c_str())->GetRange(b);
    
    
    
    m_lut->SetTableRange(m_polyData->GetPointData()->GetScalars()->GetRange());
    m_lut->GetTableRange(b);
    if(m_visOpt.isColorRangeFrom) b[0]=m_visOpt.colorRangeFrom;
    if(m_visOpt.isColorRangeTo) b[1]=m_visOpt.colorRangeTo;
    if(b[1]<=b[0]) b[1]=b[0]+0.0001;
    m_lut->SetTableRange(b[0],b[1]);
    
    if(m_visOpt.uselogscale=="yes")
        m_lut->SetScaleToLog10();
    else
        m_lut->SetScaleToLinear();
    
    m_lut->Build();
    
    SelectLookTable(&m_visOpt, m_lut);
    
    m_pConeMapper->SetLookupTable(m_lut);
    m_pConeMapper->SetScalarVisibility(1);
    m_pConeMapper->UseLookupTableScalarRangeOn();
    m_pConeMapper->Update();
    
    m_pConeActor->SetMapper(m_pConeMapper);
    
    if(m_visOpt.showLut)  colorBar();
    
}
void colorBar ()
//---------------------------------------------------------------------
{
	vtkScalarBarActor *scalarBar=vtkScalarBarActor::New();
	scalarBar->SetTitle (  m_visOpt.colorScalar.c_str() );
	scalarBar->SetLabelFormat ( "%.3g" );
	scalarBar->SetOrientationToHorizontal();
	scalarBar->SetPosition ( 0.1,0 );
	scalarBar->SetPosition2 ( 0.8,0.1 );
	scalarBar->SetLookupTable (m_lut );
	scalarBar->SetVisibility(1);

	renderer->AddActor ( scalarBar );
  
	if (scalarBar!=0)
		scalarBar->Delete();
}


//---------------------------------------------------------------------
void setRadius ()
//---------------------------------------------------------------------
{
    if (m_visOpt.nGlyphs==1)
        m_sphere->SetRadius ( m_visOpt.radius);
    
    
    else if (m_visOpt.nGlyphs==2)
    {
        
        m_cone->SetRadius ( m_visOpt.radius );
        m_cone->SetHeight (m_visOpt.height );
    }
    
    else if (m_visOpt.nGlyphs==3)
    {
        
        m_cylinder->SetRadius (m_visOpt.radius );
        m_cylinder->SetHeight ( m_visOpt.height );
    }
    else if (m_visOpt.nGlyphs==4)
    {
        m_cube->SetXLength ( m_visOpt.radius );
        m_cube->SetYLength ( m_visOpt.height );
        m_cube->SetZLength ( 1 );
        
    }
}

//---------------------------------------------------------------------
void setResolution ()
//---------------------------------------------------------------------
{
    if (m_visOpt.nGlyphs==1)
    {
        m_sphere->SetPhiResolution ( 10 );
        m_sphere->SetThetaResolution ( 20 );
    }
    
    else if (m_visOpt.nGlyphs==2)
        m_cone->SetResolution ( 10 );
    
    else if (m_visOpt.nGlyphs==3)
        m_cylinder->SetResolution ( 10);
    
    
}

std::string saveImageAsPng(int num )
//------------------------------------
{
	if(m_visOpt.numImage==4 && m_visOpt.azimuth==0 && m_visOpt.elevation==0 && m_visOpt.zoom==1 && !m_visOpt.cycle)
		return "";
    

	else
	{

		int magnification=1;
		vtkWindowToImageFilter *w2i=vtkWindowToImageFilter::New();
		std::string path;
		std::string numero;
		std::stringstream ss;
		ss<<num;
		ss>>numero;
		std::string fileName;
		std::clog <<"Before setinput"<<std::endl;
		w2i->SetInput(renWin);
		std::clog <<"After setinput"<<std::endl;
		/*VTK9 migration
		w2i->SetMagnification(magnification);
		replaced with
		w2i->SetScale(magnification);
		*/
		w2i->SetScale(magnification);

		std::clog <<"Before update"<<std::endl;
		w2i->Update();
		std::clog <<"After update"<<std::endl;
		std::clog <<"out "<< m_visOpt.imageName.c_str()<<std::endl;
		if (m_visOpt.imageName!="VisIVOServerImage")
		{
			if (m_visOpt.imageName.find("/")!=std::string::npos)
			{

				path=getDir(m_visOpt.imageName);
				if (path.at(0)!='/')
					path="./"+ path;

				m_visOpt.imageName=getName(m_visOpt.imageName);
			}

			else 
				path="./";

			

		
		}

		else

			path="./";

	

		if (m_visOpt.noDefault=="yes")
		{  
		     size_t found=m_visOpt.imageName.find(".png");

		    if(found== std::string::npos)
			fileName=path+m_visOpt.imageName+".png";
		    else
			fileName=path+m_visOpt.imageName; 

		  
//			fileName=path+m_visOpt.imageName+".png";
		}
		else
		{
		     size_t found=m_visOpt.imageName.find(".png");

		    if(found== std::string::npos)
			fileName=path+m_visOpt.imageName+numero+".png";
		    else
		        fileName=path+m_visOpt.imageName.substr(0,found)+numero+".png";

		}	
		std::clog <<"fielname "<< fileName.c_str()<<std::endl;

		vtkPNGWriter *w=vtkPNGWriter::New();
		/* VTK9 migration
		w->SetInput(w2i->GetOutput());
		replaced
					w->SetInputData(w2i->GetOutput());
		*/
			w->SetInputData(w2i->GetOutput());

		w->SetFileName(fileName.c_str());

		w->Write();  
    
		if ( w2i != 0 )
			w2i->Delete();
  
		if ( w!= 0 )
			w->Delete();
        
        return fileName;


	}

  	//this->destroyVTK();
	return "";

}
void SetNum(int n) {num=n;}
void SetOpt(VisIVOServerOptions opt) {m_visOpt = opt;}


    int Argc;
  char** Argv;
  double m_xRange[2] ,m_yRange[2] , m_zRange[2];
  VisIVOServerOptions m_visOpt;
      vtkPolyDataMapper  *m_pConeMapper;
    vtkActor           *m_pConeActor;
    vtkPolyData       *m_polyData;
    vtkGlyph3D *m_glyph ;
    vtkSphereSource   *m_sphere;
    vtkConeSource   *m_cone;
    vtkCylinderSource   *m_cylinder;
    vtkCubeSource   *m_cube;
    vtkPoints *m_points;
    vtkLookupTable      *m_lut;
    vtkRenderer* renderer;
    vtkRenderWindow* renWin;
    ExtendedGlyph3D *m_glyphFilter;
     vtkCamera          *m_camera;
     int num;
    };
    vtkStandardNewMacro(PointProcess);
}

//-----------------------------------------------------------------------------------
int PointsPipe::createPipe ()
//------------------------------------------------------------------------------------
{
    
    int rank, size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i = 0;
    //open view
    //------------------------------------------------------------------
    vtkMultiProcessController::SetGlobalController(m_pController);
    //--------------------------------------------------------------------
    //  m_pRenderer->Render();
    PointProcess* p = PointProcess::New();
    //settings..
    p->SetOpt(m_visOpt);
    p->SetNum(num);
    m_pController->SetSingleProcessObject(p);
    std::clog << "Before exec";
    m_pController->SingleMethodExecute();
    std::clog << "AFTER exec";
    
    //if(m_visOpt.showAxes) setAxes (m_polyData,bounds );
    //delete [] bounds;
    //if(rank ==0) saveImageAsPng(num);
    p->Delete();
    //contr->Finalize();
    //contr->Delete();
    //if(newVerts!=0)
        //newVerts->Delete();
    //radiusArrays->Delete();
       
    return 0;
}