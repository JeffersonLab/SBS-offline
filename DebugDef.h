//////////////////////////////////////////////////////////////////////////
//
// DebugDef.h
//
// debug handling routines for root-based classes.
//
// This file defines a number of predefined debug handling routines: 
// assert, DEBUG_<LEVEL> and DEBUG_LINE_<LEVEL>
// Where <LEVEL> could be WARNING, INFO or MASSINFO
// 
//////////////////////////////////////////////////////////////////////////
//	
// Author: Jin Huang <mailto:jinhuang@jlab.org>   Feb 2008
//	Modify History:
//      Feb 2008    DEBUG_LEVEL_RELATED_PERFORMACE_CHECKER
//
//////////////////////////////////////////////////////////////////////////
/*
                      REQUIREMENTS

 1.In the head file of a root-based class, insert the following lines BELOW 
  any other header files:

  //------------------------------------------------------//
  //
  //	Debug Definitions
  //	place this section below any other head files
  //
  //------------------------------------------------------//
  #ifdef DEBUG_LEVEL
  #   undef DEBUG_LEVEL
  #endif

  //	DEBUG_LEVEL;	
  //	=0	or not define: no debug, full speed 
  //	>=1	enable debug extra warning (suggested setting)
  //	>=2	above + enable debug assert
  //	>=3	above + enable debug extra info
  //	>=4	above + massive info (in a for or while)
  #define DEBUG_LEVEL   1

  #include    "DebugDef.h"
  //------------------------------------------------------//
  
  Where DEBUG_LEVEL will control which debug handling routines to be 
  execused as discussed in next section

 2.in the .cxx file the the same class, include the header file of that class
  BELOW any other header files

 3.(Optional, Hall A Analyzer Only, Suggested) 
  add to the constructor of the class the following line:
  DEBUG_LEVEL_RELATED_PERFORMACE_CHECKER;

 4.(Optional, Hall A Analyzer Only) 
  add to the constructor of the class the following line:
  DEBUG_HALL_A_ANALYZER_DEBUGER_INIT;

                          USAGE
  
  Once above requirement is satisified, following routines should be working:

 1.assert(exp)
  only effect when DEBUG_LEVEL>=2, other wise do nothing
  check whether exp is ture. If not, stop the code and print exp on screen

 2.DEBUG_WARNING(location,...)
  only effect when DEBUG_LEVEL>=1, other wise do nothing
  DEBUG_WARNING works just like Warning(); function of TObject. 
  See its description on
  http://root.cern.ch/root/html/TObject.html#TObject:Warning

 3.DEBUG_INFO(location,...)
  only effect when DEBUG_LEVEL>=3, other wise do nothing
  DEBUG_INFO works just like Info(); function of TObject. 
  See its description on
  http://root.cern.ch/root/html/TObject.html#TObject:Info

 4.DEBUG_MASSINFO(location,...)
  only effect when DEBUG_LEVEL>=4, other wise do nothing
  suggested usage is for step by step TRACE or info in a LOOP 
  DEBUG_MASSINFO works just like Info(); function of TObject. 
  See its description on
  http://root.cern.ch/root/html/TObject.html#TObject:Info

 5.DEBUG_LINE_LEVEL(level,exp) and DEBUG_LINE_<LEVEL>(exp)
  only execute exp when DEBUG_LEVEL>=level. Similarly, you can demonstrate
  the levle by its name using DEBUG_LINE_<LEVEL>(exp)
  Where <LEVEL> could be WARNING, INFO or MASSINFO
  
 6.Of course if you have more than one line to be controled by DEBUG_LEVEL
  #if DEBUG_LEVEL>= <some level> 
	    ...
  #endif

 7.DEBUG_LEVEL_OVERRIDE
  You can overide all DEBUG_LEVEL settings in any class using   
  thisheader file by change the value DEBUG_LEVEL_OVERRIDE to a  
  desiderated DEBUG_LEVEL value. Only a non-negative value is 
  effective.
  You can also define a DEBUG_LEVEL_OVERRIDE in MakeFile.


               GUIDELINE FOR DEBUG_LEVEL

	DEBUG_LEVEL;	
	<1  or not define: no debug, full speed 
	>=1 enable debug extra warning (suggested setting for most release code)
	>=2 above + enable debug assert
	>=3 above + enable debug extra info
	>=4 above + massive info (step by step tracing or in info within a loop)
*/
//


#if !defined(DEBUG_LEVEL_OVERRIDE) 
//------------------------------------------
//  change DEBUG_LEVEL_OVERRIDE value here
// 
//  You can overide all DEBUG_LEVEL settings in any class using   
//  thisheader file by change the value DEBUG_LEVEL_OVERRIDE to a  
//  desiderated DEBUG_LEVEL value. Only a non-negative value is 
//  effective.
//  You can also define a DEBUG_LEVEL_OVERRIDE in MakeFile.
#define DEBUG_LEVEL_OVERRIDE -1
//------------------------------------------
#endif

#if DEBUG_LEVEL_OVERRIDE>=0
#   undef DEBUG_LEVEL
#   define DEBUG_LEVEL DEBUG_LEVEL_OVERRIDE
#endif

#undef  DEBUG_WARNING
#undef  DEBUG_LINE_WARNING
#if DEBUG_LEVEL>=1//start show warning 
#   define  DEBUG_WARNING(location,...)  Warning(location,__VA_ARGS__);
#   define  DEBUG_LINE_WARNING(exp)  exp;
#else
#   define  DEBUG_WARNING(location,...)     ((void)0)
#   define  DEBUG_LINE_WARNING(exp)         ((void)0)
#endif

#undef NDEBUG
#if DEBUG_LEVEL<2
#   define NDEBUG
#endif
#include <cassert>

#undef  DEBUG_INFO
#undef  DEBUG_LINE_INFO
#if DEBUG_LEVEL>=3//start show info 
#   define  DEBUG_INFO(location,...)  Info(location,__VA_ARGS__);
#   define  DEBUG_LINE_INFO(exp)  exp;
#else
#   define  DEBUG_INFO(location,...)    ((void)0)
#   define  DEBUG_LINE_INFO(exp)        ((void)0)
#endif

#undef  DEBUG_MASSINFO
#undef  DEBUG_LINE_MASSINFO
#if DEBUG_LEVEL>=4//start show massive info. 
                  //Suitable for step by step trace or info in a loop 
#   define  DEBUG_MASSINFO(location,...)  Info(location,__VA_ARGS__);
#   define  DEBUG_LINE_MASSINFO(exp)  exp;
#else
#   define  DEBUG_MASSINFO(location,...)    ((void)0)
#   define  DEBUG_LINE_MASSINFO(exp)        ((void)0)
#endif

#undef  DEBUG_LINE_LEVEL
#define DEBUG_LINE_LEVEL(level,exp) if (DEBUG_LEVEL>=level) {exp;}

#define DEBUG_HALL_A_ANALYZER_DEBUGER_INIT \
    if (DEBUG_LEVEL>=4)\
	    SetDebug(2);\
    else if (DEBUG_LEVEL>=3)\
	    SetDebug(1);\
    else\
        SetDebug(0);

#undef DEBUG_LEVEL_RELATED_PERFORMACE_CHECKER
#if DEBUG_LEVEL>2
#   define DEBUG_LEVEL_RELATED_PERFORMACE_CHECKER \
    Warning(Here("Constructor"),\
        "\tThis class is compliled with DEBUG_LEVEL=%d."\
        "\n\tChange DEBUG_LEVEL to <=1 to gain full speed.",\
        DEBUG_LEVEL);
#else
#   define DEBUG_LEVEL_RELATED_PERFORMACE_CHECKER ((void)0)  
#endif

