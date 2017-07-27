//------------------------------------------------------------------------------
//  Copyright 2010-2012 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "RectKD.h"

//unsigned int Rect::rect_size=0;

///////////////////////////////////////////////////////
// Rect Implementation

void Rect::reportIntersections() 
{

   if (!intersection.empty()) {
        cout << vid <<" : ";
        cout <<intersection.size()<<" intersections\n";
        sort(intersection.begin(),intersection.end());
        for (Intersect::iterator i = intersection.begin(); i != intersection.end(); ++i){
		if( i!=intersection.begin() ) if( *i==*(i-1) ) continue;
                 cout << "(" << *i <<")";
        }
        cout << endl;
   }
}

template < class Pt, int K >
ostream & operator << (ostream & out, const RectKD<Pt,K>& rectKD)
{
    out << "Vertex ID: " << rectKD.getVID();
    return out;
}

