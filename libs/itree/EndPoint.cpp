//------------------------------------------------------------------------------
//  Copyright 2010-2012 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "EndPoint.h"

ostream & operator << (ostream & out, const EndPoint& endpoint) 
{
    out << endpoint.getValue();
    return out;
}
