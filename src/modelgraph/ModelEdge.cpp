//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "ModelEdge.h"
using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CModelEdge::CModelEdge(int start, int end, int id)
{
    m_Key[0]=start;
    m_Key[1]=end;
    m_Next=NULL; //m_Kids[0]=m_Kids[1]=NULL;
    m_id=id;
}
