//------------------------------------------------------------------------------
//  Copyright 2007-2019 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#pragma once

#include "polygon.h"

namespace masc{
namespace polygon{

///////////////////////////////////////////////////////////////////////////////
// This convex hull implemetation realizes the idea from
// A. Melkman, "On-line construction of the convex hull of a simple polygon",
// Info. Proc. Letters 25, 11-12 (1987)
///////////////////////////////////////////////////////////////////////////////


//
// compute the convex hull of the given (subset) polygon
//
// e mush be reachable from s
//
void hull2d(ply_vertex * s, ply_vertex * e, list<ply_vertex*>& hull );

}//end namespace polygon
}//end namespace masc
