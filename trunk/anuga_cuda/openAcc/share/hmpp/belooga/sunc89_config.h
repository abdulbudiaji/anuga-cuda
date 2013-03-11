/*
 * Copyright (C) 2008-2013 CAPS entreprise.  All Rights Reserved.
 * 
 * The source code contained or described herein and all documents  related
 * to the source code ("Material") are owned  by  CAPS  entreprise  or  its
 * suppliers or licensors.
 * 
 * Title to the Material remains with CAPS entreprise or its suppliers  and
 * licensors.  The Material contains trade secrets and proprietary and con-
 * fidential information of CAPS entreprise or its suppliers and licensors.
 * 
 * The Material is protected by the French intellectual property code,  in-
 * tellectual property laws and international treaties.  No part of the Ma-
 * terial may be used, copied, reproduced, modified,  published,  uploaded,
 * posted, transmitted, distributed or disclosed in any  way  without  CAPS
 * entreprise's prior express written permission.
 * 
 * No license under any patent, copyright, trade secret or other  intellec-
 * tual property right is granted to or conferred upon you by disclosure or
 * delivery of the Material, either expressly, by implication,  inducement,
 * estoppel or otherwise.
 * 
 * Any license under such intellectual property rights  must  be  expressed
 * and approved by CAPS entreprise in writing.
 */
/*
 * Configuration file for Belooga for both 32bit and 64bit intel architecture
 *
 */
#ifdef _BELOOGA_

#ifdef _BELOOGA_

#pragma BELOOGA_HLINE_MODE	1;

#if __i386
#pragma BELOOGA_TARGET_NAME	"sunc89-i686"
#elif __x86_64  
#pragma BELOOGA_TARGET_NAME	"sunc89-x86_64"
#else
#pragma BELOOGA_TARGET_NAME	"sunc89-unknown"
#endif


#define __extension__		/%* __extension__ *%/
#define __attribute__(ATT)	/%* __attribute__(ATT) *%/
#define __attribute(ATT)	/%* __attribute(ATT) *%/
#define __builtin_va_list	/%* __builtin_va_list *%/ char **
#define __PRETTY_FUNCTION__	/%* __PRETTY_FUNCTION__ *%/ 1

#define __asm__			/%* __asm__ *%/ asm
#define __asm			/%* __asm *%/ asm
#define __asm__volatile		/%* __asm__volatile *%/ asm volatile
#define __volatile		/%* __volatile *%/ volatile
#define __volatile__		/%* __volatile__ *%/ volatile
#define __typeof		/%* __typeof *%/ typeof
#define __typeof__		/%* __typeof__ *%/ typeof
#define __alignof		/%* __alignof *%/ __builtin_alignof
#define __alignof__		/%* __alignof__ *%/ __builtin_alignof

#include <sys/cdefs.h>

#ifndef __attribute__
#  define __attribute__(ATT)	/%* __attribute__(ATT) *%/
#endif

/* ************************************************************************************* */
/* Type Definitions :  		token,		category or variant			*/
#pragma BELOOGA_DEFINE_TYPE_NODE	"__const";	865;
#pragma BELOOGA_DEFINE_TYPE_NODE	"__const__";	865;
#pragma BELOOGA_DEFINE_TYPE_NODE	"__restrict";	868;
#pragma BELOOGA_DEFINE_TYPE_NODE	"__restrict__";	868;
#pragma BELOOGA_DEFINE_TYPE_NODE	"__inline";	904;
#pragma BELOOGA_DEFINE_TYPE_NODE	"__inline__";	904;
#pragma BELOOGA_DEFINE_TYPE_NODE	"__signed";	833;
#pragma BELOOGA_DEFINE_TYPE_NODE	"__signed__";	833;

#pragma BELOOGA_DEFINE_SYMBOL_NODE	"__func__";	810;	921;

#pragma BELOOGA_UNPARSE_REGPLACE "__builtin_va_list char  **" "__builtin_va_list "
#pragma BELOOGA_UNPARSE_REGPLACE "__PRETTY_FUNCTION__ 1" "__PRETTY_FUNCTION__"

#pragma BELOOGA_UNPARSE_REGPLACE "__asm__ asm" "__asm__"
#pragma BELOOGA_UNPARSE_REGPLACE "__asm asm" "__asm"
#pragma BELOOGA_UNPARSE_REGPLACE "__asm__volatile asm volatile" "__asm__volatile"
#pragma BELOOGA_UNPARSE_REGPLACE "__volatile volatile" "__volatile"
#pragma BELOOGA_UNPARSE_REGPLACE "__volatile *volatile" "*__volatile"
#pragma BELOOGA_UNPARSE_REGPLACE "__volatile__ volatile" "__volatile__"
#pragma BELOOGA_UNPARSE_REGPLACE "__alignof __builtin_alignof" "__alignof"
#pragma BELOOGA_UNPARSE_REGPLACE "__alignof__ __builtin_alignof" "__alignof__"

/*    BUILTIN_BOOL_8  		= 1 ,  8bit  boolean */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     1; 818;
						/* TAG_T_BOOL  */

/*    BUILTIN_INT_8   		= 5 ,  8bit character (not numerical) */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     5; 810;
						/* TAG_T_CHAR  */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     5; 810; 833;
 					/* TAG_T_CHAR <- TAG_T_SIGNED  */

/*     BUILTIN_INT_16  		= 6 ,  16bit signed integer (two-complement)  */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     6; 836;
						/* TAG_T_SHORT */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     6; 811;
						/* TAG_T_SHORT_INT  */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     6; 812; 836;
					/* TAG_T_INT <- TAG_T_SHORT */

#pragma BELOOGA_DEFINE_BUILTIN_TYPE     6; 836; 833;
					/* TAG_T_SHORT <- TAG_T_SIGNED*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     6; 833; 836;
					/* TAG_T_SIGNED <- TAG_T_SHORT*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     6; 811;	833;
					/* TAG_T_SHORT_INT <- TAG_T_SIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     6; 833;	811;
					/* TAG_T_SIGNED <- TAG_T_SHORT_INT*/

#pragma BELOOGA_DEFINE_BUILTIN_TYPE     6; 812; 836; 833;
				/* TAG_T_INT <- TAG_T_SHORT <- TAG_T_SIGNED*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     6; 812; 833; 836;
				/* TAG_T_INT <- TAG_T_SIGNED <- TAG_T_SHORT */

/*     BUILTIN_INT_32  		= 7 ,  32bit signed integer (two-complement)  */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 800;
						/* TAG_T_DEFAULT  */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 812;
						/* TAG_T_INT  */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 812; 833;
				        /* TAG_T_INT <- TAG_T_SIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 833;
				                /* TAG_T_SIGNED */

//#if __LONG_MAX__ == __INT_MAX__
#if 0
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 813;
						/* TAG_T_LONG_INT */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 813; 833;
					/* TAG_T_LONG_INT <- TAG_T_SIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 812; 840;
					/* TAG_T_INT <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 812; 840; 833;
				/* TAG_T_INT <- TAG_T_LONG <- TAG_T_SIGNED*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 812; 833; 840;
				/* TAG_T_INT <- TAG_T_SIGNED <- TAG_T_LONG*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 840;
						/* TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 840; 833;
					/* TAG_T_LONG <- TAG_T_SIGNED*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     7; 833; 840;
 					/* TAG_T_SIGNED <- TAG_T_LONG*/
#endif

/*     BUILTIN_INT_64  		= 8 ,  64bit signed integer (two-complement) */

//#if __LONG_MAX__ == __LONG_LONG_MAX__
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 813;
						/* TAG_T_LONG_INT */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 813; 833;
					/* TAG_T_LONG_INT <- TAG_T_SIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 812; 840;
					/* TAG_T_INT <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 812; 840; 833;
				/* TAG_T_INT <- TAG_T_LONG <- TAG_T_SIGNED*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 812; 833; 840;
				/* TAG_T_INT <- TAG_T_SIGNED <- TAG_T_LONG*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 840;
						/* TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 840; 833;
					/* TAG_T_LONG <- TAG_T_SIGNED*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 833; 840;
 					/* TAG_T_SIGNED <- TAG_T_LONG*/
//#endif

#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 814;						
/* TAG_T_LLONG_INT */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 833; 814;					
/* TAG_T_LLONG_INT <- TAG_T_SIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 813; 840;					
/* TAG_T_LONG_INT <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 813; 840; 833;				
/* TAG_T_LONG_INT <- TAG_T_LONG <- TAG_T_SIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 813; 833; 840;				
/* TAG_T_LONG_INT <- TAG_T_SIGNED <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 840; 840;					
/* TAG_T_LONG <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 840; 840; 833;				
/* TAG_T_LONG <- TAG_T_LONG <- TAG_T_SIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 840; 833; 840;				
/* TAG_T_LONG <- TAG_T_SIGNED <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 833; 840; 840;				
/* TAG_T_SIGNED <- TAG_T_LONG <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 812; 840; 840;				
/* TAG_T_INT <- TAG_T_LONG <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 812; 840; 840; 833;				
/* TAG_T_INT <- TAG_T_LONG <- TAG_T_LONG <- TAG_T_SIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 812; 840; 833; 840;				
/* TAG_T_INT <- TAG_T_LONG <- TAG_T_SIGNED <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE     8; 812; 833; 840; 840;				
/* TAG_T_INT <- TAG_T_SIGNED <- TAG_T_LONG <- TAG_T_LONG*/

/*     BUILTIN_UNSIGNED_8   	= 10 ,  8bit unsigned integer */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    10; 810; 834;					
/* TAG_T_CHAR <- TAG_T_UNSIGNED  */

/*     BUILTIN_UNSIGNED_16  	= 11 ,  16bit unsigned integer */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    11; 811; 834;				
	/* TAG_T_SHORT_INT <- TAG_T_UNSIGNED  */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    11; 812; 836; 834;				
/* TAG_T_INT <- TAG_T_SHORT <- TAG_T_UNSIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    11; 812; 834; 836;				
/* TAG_T_INT <- TAG_T_UNSIGNED <- TAG_T_SHORT */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    11; 836; 834;					
/* TAG_T_SHORT <- TAG_T_UNSIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    11; 834; 836;					
/* TAG_T_UNSIGNED <- TAG_T_SHORT */

/*     BUILTIN_UNSIGNED_32  	= 12 ,  32bit unsigned integer */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    12; 834;						
/* TAG_T_UNSIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    12; 812; 834;					
/* TAG_T_INT <- TAG_T_UNSIGNED */

//#if __LONG_MAX__ == __INT_MAX__
#if 0
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    12; 813; 834;					
/* TAG_T_LONG_INT <- TAG_T_UNSIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    12; 812; 840; 834;				
/* TAG_T_INT <- TAG_T_LONG <- TAG_T_UNSIGNED*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    12; 812; 834; 840;				
/* TAG_T_INT <- TAG_T_UNSIGNED <- TAG_T_LONG*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    12; 840; 834;					
/* TAG_T_LONG <- TAG_T_UNSIGNED*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    12; 834; 840; 					
/* TAG_T_UNSIGNED <- TAG_T_LONG*/
#endif

/*     BUILTIN_UNSIGNED_64  	= 13 ,  64bit unsigned integer */
//#if __LONG_MAX__ == __LONG_LONG_MAX__
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 813; 834;					
/* TAG_T_LONG_INT <- TAG_T_UNSIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 812; 840; 834;				
/* TAG_T_INT <- TAG_T_LONG <- TAG_T_UNSIGNED*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 812; 834; 840;				
/* TAG_T_INT <- TAG_T_UNSIGNED <- TAG_T_LONG*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 840; 834;					
/* TAG_T_LONG <- TAG_T_UNSIGNED*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 834; 840; 					
/* TAG_T_UNSIGNED <- TAG_T_LONG*/
//#endif

#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 834; 814;					
/* TAG_T_LLONG_INT <- TAG_T_UNSIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 813; 840; 834;				
/* TAG_T_LONG_INT <- TAG_T_LONG <- TAG_T_UNSIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 813; 834; 840;				
/* TAG_T_LONG_INT <- TAG_T_UNSIGNED <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 840; 840; 834;				
/* TAG_T_LONG <- TAG_T_LONG <- TAG_T_UNSIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 840; 834; 840;				
/* TAG_T_LONG <- TAG_T_UNSIGNED <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 834; 840; 840;				
/* TAG_T_UNSIGNED <- TAG_T_LONG <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 812; 840; 840; 834;			
	/* TAG_T_INT <- TAG_T_LONG <- TAG_T_LONG <- TAG_T_UNSIGNED */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 812; 840; 834; 840;				
/* TAG_T_INT <- TAG_T_LONG <- TAG_T_UNSIGNED <- TAG_T_LONG */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 812; 834; 840; 840;				
/* TAG_T_INT <- TAG_T_UNSIGNED <- TAG_T_LONG <- TAG_T_LONG*/
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 834; 814;					
/* TAG_T_UNSIGNED <- TAG_T_LLONG_INT  */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    13; 814; 834;				
	/* TAG_T_LLONG_INT <- T_UNSIGNED  */

/*    BUILTIN_FLOAT_32    	= 19 ,  32bit IEEE 754 floating point  */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    19; 815;			                        
/* TAG_T_FLOAT */

/*    BUILTIN_FLOAT_64 	        = 20 ,  64bit IEEE 754 floating point  */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    20; 816;			                        
/* TAG_T_DOUBLE */

/*    BUILTIN_COMPLEX_32     	= 24 ,  2*32bit IEEE 754 floating point complex */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    24; 815; 849;			                
/* TAG_T_FLOAT <- TAG_T_COMPLEX */

/*    BUILTIN_COMPLEX_64 	= 25 ,  2*64bit IEEE 754 floating point complex */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    25; 849;				                
/* TAG_T_COMPLEX */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    25; 816; 849;			                
/* TAG_T_DOUBLE <- TAG_T_COMPLEX */

//#if __LDBL_MANT_DIG__ == 64 

#ifdef __i386
/*    BUILTIN_FLOAT_80_12   	= 22 ,  80bit floating point stored in 12 bytes (Intel 32bit) */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    22; 817;			                        
/* TAG_T_LONG_DOUBLE */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    22; 816; 840;			                
/* TAG_T_DOUBLE <- TAG_T_LONG */
/*    BUILTIN_COMPLEX_80_12 	= 27 ,  2*80bit Intel floating point complex (Intel) */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    27; 817; 849;			                
/* TAG_T_LONG_DOUBLE <- TAG_T_COMPLEX */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    27; 816; 840; 849;		                
/* TAG_T_DOUBLE <- TAG_T_LONG <- TAG_T_COMPLEX */
#endif 

//#ifdef __amd64 
/*    BUILTIN_FLOAT_80_16   	= 23 ,  80bit floating point stored in 12 bytes (Intel 64bit) */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    23; 817;			                        
/* TAG_T_LONG_DOUBLE */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    23; 816; 840;			                
/* TAG_T_DOUBLE <- TAG_T_LONG */
/*    BUILTIN_COMPLEX_80_16 	= 28 ,  2*80bit Intel floating point complex (Intel) */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    28; 817; 849;			                
/* TAG_T_LONG_DOUBLE <- TAG_T_COMPLEX */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    28; 816; 840; 849;		                
/* TAG_T_DOUBLE <- TAG_T_LONG <- TAG_T_COMPLEX */
//#endif

//#endif

/*    BUILTIN_IMAGINARY_32     	= 24 ,  2*32bit IEEE 754 floating point complex */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    32; 815; 850;			                /* TAG_T_FLOAT <- TAG_T_IMAGINARY */

/*    BUILTIN_IMAGINARY_64 	= 25 ,  2*64bit IEEE 754 floating point complex */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    33; 850;				                /* TAG_T_IMAGINARY */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    33; 816; 850;			                /* TAG_T_DOUBLE <- TAG_T_IMAGINARY */

//#if __LDBL_MANT_DIG__ == 64 

#ifdef __i386
/*    BUILTIN_IMAGINARY_80_12 	= 27 ,  2*80bit Intel floating point complex (Intel) */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    35; 817; 850;			                /* TAG_T_LONG_DOUBLE <- TAG_T_IMAGINARY */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    35; 816; 840; 850;		                /* TAG_T_DOUBLE <- TAG_T_LONG <- TAG_T_IMAGINARY */
#endif 

//#ifdef __amd64 
/*    BUILTIN_IMAGINARY_80_16 	= 28 ,  2*80bit Intel floating point complex (Intel) */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    36; 817; 850;			                /* TAG_T_LONG_DOUBLE <- TAG_T_IMAGINARY */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    36; 816; 840; 850;		                /* TAG_T_DOUBLE <- TAG_T_LONG <- TAG_T_IMAGINARY */
//#endif

//#endif

//#if __LONG_MAX__ == __INT_MAX__
#if 0
/*    BUILTIN_POINTER_32     	= 29,   32bit pointer */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    29; 802; 921;			                
/* TAG_T_VOID <- TAG_T_POINTER */
#endif

//#if __LONG_MAX__ == __LONG_LONG_MAX__
/*    BUILTIN_POINTER_64     	= 30 ,  64bit pointer */
#pragma BELOOGA_DEFINE_BUILTIN_TYPE    30; 802; 921;			                
/* TAG_T_VOID <- TAG_T_POINTER */
//#endif

#endif


#endif
