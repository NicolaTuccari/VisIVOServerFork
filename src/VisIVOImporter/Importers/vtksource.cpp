#include "vtksource.h"

//---------------------------------------------------------------------
VTKSource::VTKSource()
//---------------------------------------------------------------------
{
  m_pointsFileName = "";
  m_pointsBinaryName="";
  m_nRows = 0;
  m_nCols = 0;
  m_fieldNames.clear();
   
}


//---------------------------------------------------------------------
int VTKSource::readHeader()
//---------------------------------------------------------------------
{
  return 1;
}

//---------------------------------------------------------------------
int VTKSource::readData()
//---------------------------------------------------------------------
{
  return 1;
}

