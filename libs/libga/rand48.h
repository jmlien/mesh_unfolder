
#pragma once

///Added by Jyh-Ming Lien

#ifdef _WIN32

////////////////////////////////////////////////////////////////////////////////////////
// Following functions define M_PI and drand48, which are not starndard c library and 
// definitions. In addition, rint used to round off float points to int is also here.
/////////////////////////////////////////////////////////////////////////////////////////

#define M_PI 3.1415926 //reference PI above
namespace masc {
	namespace ga {
		extern "C" {
			//Implementation of these functions are located in util.cpp
			float drand48();
			float erand48(register unsigned short *xsubi);
			long irand48(register unsigned short m);
			long krand48(register unsigned short *xsubi, unsigned short m);
			long lrand48();
			long mrand48();
			static void next();
			void srand48(long seedval);
			unsigned short * seed48(unsigned short seed16v[3]);
			void lcong48(unsigned short param[7]);
			long nrand48(register unsigned short *xsubi);
			long jrand48(register unsigned short *xsubi);

			/**Round to closest integer.
			*The rint() function rounds x to an integer value according
			*to the prevalent rounding mode.  The default rounding mode
			*is to round to the nearest integer.
			*@return The  rint() function returns the integer value as a float-
			*ing-point number.
			*/
			float rint(float x);

		} //end extern "C"
	}
}

#endif //_WIN32