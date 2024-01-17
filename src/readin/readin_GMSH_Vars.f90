!============================================================================================================ xX =================
!        _____     _____    _______________    _______________   _______________            .xXXXXXXXx.       X
!       /    /)   /    /)  /    _____     /)  /    _____     /) /    _____     /)        .XXXXXXXXXXXXXXx  .XXXXx
!      /    //   /    //  /    /)___/    //  /    /)___/    // /    /)___/    //       .XXXXXXXXXXXXXXXXXXXXXXXXXx
!     /    //___/    //  /    //   /    //  /    //___/    // /    //___/    //      .XXXXXXXXXXXXXXXXXXXXXXX`
!    /    _____     //  /    //   /    //  /    __________// /    __      __//      .XX``XXXXXXXXXXXXXXXXX`
!   /    /)___/    //  /    //   /    //  /    /)_________) /    /)_|    |__)      XX`   `XXXXX`     .X`
!  /    //   /    //  /    //___/    //  /    //           /    //  |    |_       XX      XXX`      .`
! /____//   /____//  /______________//  /____//           /____//   |_____/)    ,X`      XXX`
! )____)    )____)   )______________)   )____)            )____)    )_____)   ,xX`     .XX`
!                                                                           xxX`      XXx
! Copyright (C) 2017 Claus-Dieter Munz <munz@iag.uni-stuttgart.de>
! This file is part of HOPR, a software for the generation of high-order meshes.
!
! HOPR is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! HOPR is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with HOPR. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "hopr.h"
MODULE MOD_Readin_GMSH_Vars
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC

INTEGER                            :: bOrd
INTEGER,ALLOCATABLE                :: tetMapGMSH(:,:,:),  pyrMapGMSH(:,:,:),  priMapGMSH(:,:,:),  hexMapGMSH(:,:,:)
INTEGER                            :: tetMapCGNSToGMSH(4),pyrMapCGNSToGMSH(5),priMapCGNSToGMSH(6),hexMapCGNSToGMSH(8)
INTEGER                            :: quadMapCGNSToGMSH(4)
INTEGER                            :: GMSH_TYPES(6,131)
INTEGER                            :: nBCs_GMSH
INTEGER,ALLOCATABLE                :: MapEntityToBC(:)          ! MeshFormat: 4.1
INTEGER,ALLOCATABLE                :: MapBC(:),MapBCInd(:)      ! MeshFormat: 2.2

CONTAINS
SUBROUTINE buildTypes()
!===================================================================================================================================
! Build list of GMSH Types
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                   :: tmp(UBOUND(GMSH_TYPES,1)*UBOUND(GMSH_TYPES,2))  ! ?
!===================================================================================================================================
! BaseElement (GMSH), ElemType(our), nNodes, Degree, Complete, Dim 
tmp=&
(/2 ,  2   , 2   ,  1 ,  1, 2,& ! MSH_LIN_2   1  
  3 ,  3   , 3   ,  1 ,  1, 2,& ! MSH_TRI_3   2  
  4 ,  4   , 4   ,  1 ,  1, 2,& ! MSH_QUA_4   3  
  5 ,  104 , 4   ,  1 ,  1, 3,& ! MSH_TET_4   4  
  8 ,  108 , 8   ,  1 ,  1, 3,& ! MSH_HEX_8   5  
  7 ,  106 , 6   ,  1 ,  1, 3,& ! MSH_PRI_6   6  
  6 ,  105 , 5   ,  1 ,  1, 3,& ! MSH_PYR_5   7  
  2 ,  202 , 3   ,  2 ,  1, 2,& ! MSH_LIN_3   8  
  3 ,  6   , 6   ,  2 ,  1, 2,& ! MSH_TRI_6   9  
  4 ,  7   , 9   ,  2 ,  1, 2,& ! MSH_QUA_9   10 
  5 ,  204 , 10  ,  2 ,  1, 3,& ! MSH_TET_10  11 
  8 ,  208 , 27  ,  2 ,  1, 3,& ! MSH_HEX_27  12 
  7 ,  206 , 18  ,  2 ,  1, 3,& ! MSH_PRI_18  13 
  6 ,  205 , 14  ,  2 ,  1, 3,& ! MSH_PYR_14  14 
  1 ,  1   , 1   ,  1 ,  1, 1,& ! MSH_PNT     15 
  4 ,  7   , 8   ,  2 ,  0, 2,& ! MSH_QUA_8   16 
  8 ,  208 , 20  ,  2 ,  0, 3,& ! MSH_HEX_20  17 
  7 ,  206 , 15  ,  2 ,  0, 3,& ! MSH_PRI_15  18 
  6 ,  205 , 13  ,  2 ,  0, 3,& ! MSH_PYR_13  19 
  3 ,  5   , 9   ,  3 ,  0, 2,& ! MSH_TRI_9   20 
  3 ,  5   , 10  ,  3 ,  1, 2,& ! MSH_TRI_10  21 
  3 ,  5   , 12  ,  4 ,  0, 2,& ! MSH_TRI_12  22 
  3 ,  5   , 15  ,  4 ,  1, 2,& ! MSH_TRI_15  23 
  3 ,  5   , 15  ,  5 ,  0, 2,& ! MSH_TRI_15I 24 
  3 ,  5   , 21  ,  5 ,  1, 2,& ! MSH_TRI_21  25 
  2 ,  202 , 4   ,  3 ,  1, 2,& ! MSH_LIN_4   26 
  2 ,  202 , 5   ,  3 ,  1, 2,& ! MSH_LIN_5   27 
  2 ,  202 , 6   ,  3 ,  1, 2,& ! MSH_LIN_6   28 
  4 ,  204 , 20  ,  3 ,  1, 3,& ! MSH_TET_20  29 
  4 ,  204 , 35  ,  4 ,  1, 3,& ! MSH_TET_35  30 
  4 ,  204 , 56  ,  5 ,  1, 3,& ! MSH_TET_56  31 
  4 ,  204 , 34  ,  4 ,  0, 3,& ! MSH_TET_34  32 
  4 ,  204 , 52  ,  5 ,  0, 3,& ! MSH_TET_52  33 
  9 , -1   ,-1   , -1 , -1,-2,& ! MSH_POLYG_  34 
  10, -1   ,-1   , -1 , -1,-3,& ! MSH_POLYH_  35 
  4 ,  7   , 16  ,  3 ,  1, 2,& ! MSH_QUA_16  36 
  4 ,  7   , 25  ,  4 ,  1, 2,& ! MSH_QUA_25  37 
  4 ,  7   , 36  ,  5 ,  1, 2,& ! MSH_QUA_36  38 
  4 ,  7   , 12  ,  3 ,  0, 2,& ! MSH_QUA_12  39 
  4 ,  7   , 16  ,  4 ,  0, 2,& ! MSH_QUA_16I 40 
  4 ,  7   , 20  ,  5 ,  0, 2,& ! MSH_QUA_20  41 
  3 ,  5   , 28  ,  6 ,  1, 2,& ! MSH_TRI_28  42 
  3 ,  5   , 36  ,  7 ,  1, 2,& ! MSH_TRI_36  43 
  3 ,  5   , 45  ,  8 ,  1, 2,& ! MSH_TRI_45  44 
  3 ,  5   , 55  ,  9 ,  1, 2,& ! MSH_TRI_55  45 
  3 ,  5   , 66  ,  10,  1, 2,& ! MSH_TRI_66  46 
  4 ,  7   , 49  ,  6 ,  1, 2,& ! MSH_QUA_49  47 
  4 ,  7   , 64  ,  7 ,  1, 2,& ! MSH_QUA_64  48 
  4 ,  7   , 81  ,  8 ,  1, 2,& ! MSH_QUA_81  49 
  4 ,  7   , 100 ,  9 ,  1, 2,& ! MSH_QUA_100 50 
  4 ,  7   , 121 ,  10,  1, 2,& ! MSH_QUA_121 51 
  3 ,  5   , 18  ,  6 ,  0, 2,& ! MSH_TRI_18  52 
  3 ,  5   , 21  ,  7 ,  0, 2,& ! MSH_TRI_21I 53 
  3 ,  5   , 24  ,  8 ,  0, 2,& ! MSH_TRI_24  54 
  3 ,  5   , 27  ,  9 ,  0, 2,& ! MSH_TRI_27  55 
  3 ,  5   , 30  ,  10,  0, 2,& ! MSH_TRI_30  56 
  4 ,  7   , 24  ,  6 ,  0, 2,& ! MSH_QUA_24  57 
  4 ,  7   , 28  ,  7 ,  0, 2,& ! MSH_QUA_28  58 
  4 ,  7   , 32  ,  8 ,  0, 2,& ! MSH_QUA_32  59 
  4 ,  7   , 36  ,  9 ,  0, 2,& ! MSH_QUA_36I 60 
  4 ,  7   , 40  ,  10,  0, 2,& ! MSH_QUA_40  61 
  2 ,  202 , 7   ,  6 ,  1, 2,& ! MSH_LIN_7   62 
  2 ,  202 , 8   ,  7 ,  1, 2,& ! MSH_LIN_8   63 
  2 ,  202 , 9   ,  8 ,  1, 2,& ! MSH_LIN_9   64 
  2 ,  202 , 10  ,  9 ,  1, 2,& ! MSH_LIN_10  65 
  2 ,  202 , 11  ,  10,  1, 2,& ! MSH_LIN_11  66 
  2 , -1   ,-1   , -1 , -1,-2,& ! MSH_LIN_B   67 
  3 , -1   ,-1   , -1 , -1,-2,& ! MSH_TRI_B   68 
  9 , -1   ,-1   , -1 , -1,-2,& ! MSH_POLYG_B 69 
  2 , -1   ,-1   , -1 , -1,-2,& ! MSH_LIN_C   70 
  5 ,  204 , 84  ,  6 ,  1, 3,& ! MSH_TET_84  71 
  5 ,  204 , 120 ,  7 ,  1, 3,& ! MSH_TET_120 72 
  5 ,  204 , 165 ,  8 ,  1, 3,& ! MSH_TET_165 73 
  5 ,  204 , 220 ,  9 ,  1, 3,& ! MSH_TET_220 74 
  5 ,  204 , 286 ,  10,  1, 3,& ! MSH_TET_286 75 
  1 , -1   ,-1   , -1 , -1,-9,& ! dummy       76
  1 , -1   ,-1   , -1 , -1,-9,& ! dummy       77
  1 , -1   ,-1   , -1 , -1,-9,& ! dummy       78
  5 ,  204 , 74  ,  6 ,  0, 3,& ! MSH_TET_74  79 
  5 ,  204 , 100 ,  7 ,  0, 3,& ! MSH_TET_100 80 
  5 ,  204 , 130 ,  8 ,  0, 3,& ! MSH_TET_130 81 
  5 ,  204 , 164 ,  9 ,  0, 3,& ! MSH_TET_164 82 
  5 ,  204 , 202 ,  10,  0, 3,& ! MSH_TET_202 83 
  2 , -1   ,-1   , -1 , -1,-2,& ! MSH_LIN_1   84 
  3 , -1   ,-1   , -1 , -1,-2,& ! MSH_TRI_1   85 
  4 , -1   ,-1   , -1 , -1,-2,& ! MSH_QUA_1   86 
  5 , -1   ,-1   , -1 , -1,-3,& ! MSH_TET_1   87 
  8 , -1   ,-1   , -1 , -1,-3,& ! MSH_HEX_1   88 
  7 , -1   ,-1   , -1 , -1,-3,& ! MSH_PRI_1   89 
  7 ,  206 , 40  ,  3 ,  1, 3,& ! MSH_PRI_40  90 
  7 ,  206 , 75  ,  4 ,  1, 3,& ! MSH_PRI_75  91 
  8 ,  208 , 64  ,  3 ,  1, 3,& ! MSH_HEX_64  92 
  8 ,  208 , 125 ,  4 ,  1, 3,& ! MSH_HEX_125 93 
  8 ,  208 , 216 ,  5 ,  1, 3,& ! MSH_HEX_216 94 
  8 ,  208 , 343 ,  6 ,  1, 3,& ! MSH_HEX_343 95 
  8 ,  208 , 512 ,  7 ,  1, 3,& ! MSH_HEX_512 96 
  8 ,  208 , 729 ,  8 ,  1, 3,& ! MSH_HEX_729 97 
  8 ,  208 , 1000,  9 ,  1, 3,& ! MSH_HEX_100 98 
  8 ,  208 , 56  ,  3 ,  0, 3,& ! MSH_HEX_56  99 
  8 ,  208 , 98  ,  4 ,  0, 3,& ! MSH_HEX_98  100
  8 ,  208 , 152 ,  5 ,  0, 3,& ! MSH_HEX_152 101
  8 ,  208 , 222 ,  6 ,  0, 3,& ! MSH_HEX_222 102
  8 ,  208 , 296 ,  7 ,  0, 3,& ! MSH_HEX_296 103
  8 ,  208 , 386 ,  8 ,  0, 3,& ! MSH_HEX_386 104
  8 ,  208 , 488 ,  9 ,  0, 3,& ! MSH_HEX_488 105
  7 ,  206 , 126 ,  5 ,  1, 3,& ! MSH_PRI_126 106
  7 ,  206 , 196 ,  6 ,  1, 3,& ! MSH_PRI_196 107
  7 ,  206 , 288 ,  7 ,  1, 3,& ! MSH_PRI_288 108
  7 ,  206 , 405 ,  8 ,  1, 3,& ! MSH_PRI_405 109
  7 ,  206 , 550 ,  9 ,  1, 3,& ! MSH_PRI_550 110
  7 ,  206 , 38  ,  3 ,  0, 3,& ! MSH_PRI_38  111
  7 ,  206 , 66  ,  4 ,  0, 3,& ! MSH_PRI_66  112
  7 ,  206 , 102 ,  5 ,  0, 3,& ! MSH_PRI_102 113
  7 ,  206 , 146 ,  6 ,  0, 3,& ! MSH_PRI_146 114
  7 ,  206 , 198 ,  7 ,  0, 3,& ! MSH_PRI_198 115
  7 ,  206 , 258 ,  8 ,  0, 3,& ! MSH_PRI_258 116
  7 ,  206 , 326 ,  9 ,  0, 3,& ! MSH_PRI_326 117
  6 ,  205 , 30  ,  3 ,  1, 3,& ! MSH_PYR_30  118
  6 ,  205 , 55  ,  4 ,  1, 3,& ! MSH_PYR_55  119
  6 ,  205 , 91  ,  5 ,  1, 3,& ! MSH_PYR_91  120
  6 ,  205 , 140 ,  6 ,  1, 3,& ! MSH_PYR_140 121
  6 ,  205 , 204 ,  7 ,  1, 3,& ! MSH_PYR_204 122
  6 ,  205 , 285 ,  8 ,  1, 3,& ! MSH_PYR_285 123
  6 ,  205 , 385 ,  9 ,  1, 3,& ! MSH_PYR_385 124
  6 ,  205 , 29  ,  3 ,  0, 3,& ! MSH_PYR_29  125
  6 ,  205 , 50  ,  4 ,  0, 3,& ! MSH_PYR_50  126
  6 ,  205 , 77  ,  5 ,  0, 3,& ! MSH_PYR_77  127
  6 ,  205 , 110 ,  6 ,  0, 3,& ! MSH_PYR_110 128
  6 ,  205 , 149 ,  7 ,  0, 3,& ! MSH_PYR_149 129
  6 ,  205 , 194 ,  8 ,  0, 3,& ! MSH_PYR_194 130
  6 ,  205 , 245 ,  9 ,  0, 3/) ! MSH_PYR_245 131
GMSH_TYPES=RESHAPE(tmp,(/UBOUND(GMSH_TYPES,1),UBOUND(GMSH_TYPES,2)/))

END SUBROUTINE buildTypes

SUBROUTINE getGMSHVolumeMapping()
!===================================================================================================================================
! Define mapping from GMSH .msh structure to tensorproduct (i,j,k) structure
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
ALLOCATE(tetMapGMSH(0:bOrd-1,0:bOrd-1,0:bOrd-1),pyrMapGMSH(0:bOrd-1,0:bOrd-1,0:bOrd-1))
ALLOCATE(priMapGMSH(0:bOrd-1,0:bOrd-1,0:bOrd-1),hexMapGMSH(0:bOrd-1,0:bOrd-1,0:bOrd-1))
tetMapCGNSToGMSH = (/1,4,2,3/)
pyrMapCGNSToGMSH = (/1,2,3,4,5/)
priMapCGNSToGMSH = (/1,2,3,4,5,6/)
hexMapCGNSToGMSH = (/4,3,7,8,1,2,6,5/)

quadMapCGNSToGMSH = (/2,3,4,1/)

SELECT CASE(bOrd)
CASE(2)
  tetMapGMSH(0,0,0)= 1
  tetMapGMSH(1,0,0)= 4
  tetMapGMSH(0,1,0)= 2
  tetMapGMSH(0,0,1)= 3
  !
  pyrMapGMSH(0,0,0)= 1
  pyrMapGMSH(1,0,0)= 2
  pyrMapGMSH(0,1,0)= 4
  pyrMapGMSH(1,1,0)= 3
  pyrMapGMSH(0,0,1)= 5

  !priMapGMSH(0,0,0)= 1
  !priMapGMSH(1,0,0)= 2
  !priMapGMSH(0,1,0)= 3
  !priMapGMSH(0,0,1)= 4
  !priMapGMSH(1,0,1)= 5
  !priMapGMSH(0,1,1)= 6
  !
  hexMapGMSH(0,0,0)= 4
  hexMapGMSH(1,0,0)= 3
  hexMapGMSH(0,1,0)= 8
  hexMapGMSH(1,1,0)= 7
  hexMapGMSH(0,0,1)= 1
  hexMapGMSH(1,0,1)= 2
  hexMapGMSH(0,1,1)= 5
  hexMapGMSH(1,1,1)= 6

CASE(3)
  tetMapGMSH(0,0,0)=  1
  tetMapGMSH(1,0,0)=  8
  tetMapGMSH(2,0,0)=  4
  tetMapGMSH(0,1,0)=  5
  tetMapGMSH(1,1,0)= 10
  tetMapGMSH(0,2,0)=  2
  tetMapGMSH(0,0,1)=  7
  tetMapGMSH(1,0,1)=  9
  tetMapGMSH(0,1,1)=  6
  tetMapGMSH(0,0,2)=  3

  !pyrMapGMSH(0,0,0)=  1
  !pyrMapGMSH(2,0,0)=  2
  !pyrMapGMSH(2,2,0)=  3
  !pyrMapGMSH(0,2,0)=  4
  !pyrMapGMSH(0,0,2)=  5
  !pyrMapGMSH(1,0,0)=  6
  !pyrMapGMSH(0,1,0)=  7
  !pyrMapGMSH(0,0,1)=  8
  !pyrMapGMSH(2,1,0)=  9
  !pyrMapGMSH(1,0,1)= 10
  !pyrMapGMSH(1,2,0)= 11
  !pyrMapGMSH(1,1,1)= 12
  !pyrMapGMSH(0,1,1)= 13
  !pyrMapGMSH(1,1,0)= 14

  !priMapGMSH(0,0,0)=  1 
  !priMapGMSH(2,0,0)=  2 
  !priMapGMSH(0,2,0)=  3 
  !priMapGMSH(0,0,2)=  4 
  !priMapGMSH(2,0,2)=  5 
  !priMapGMSH(0,2,2)=  6 
  !priMapGMSH(1,0,0)=  7 
  !priMapGMSH(0,1,0)=  8 
  !priMapGMSH(0,0,1)=  9 
  !priMapGMSH(1,1,0)= 10
  !priMapGMSH(2,0,1)= 11
  !priMapGMSH(0,2,1)= 12
  !priMapGMSH(1,0,2)= 13
  !priMapGMSH(0,1,2)= 14
  !priMapGMSH(1,1,2)= 15
  !priMapGMSH(1,0,1)= 16
  !priMapGMSH(0,1,1)= 17
  !priMapGMSH(1,1,1)= 18

  hexMapGMSH(0,0,0)=  4
  hexMapGMSH(1,0,0)= 14
  hexMapGMSH(2,0,0)=  3
  hexMapGMSH(0,1,0)= 16
  hexMapGMSH(1,1,0)= 25
  hexMapGMSH(2,1,0)= 15
  hexMapGMSH(0,2,0)=  8
  hexMapGMSH(1,2,0)= 20
  hexMapGMSH(2,2,0)=  7
  hexMapGMSH(0,0,1)= 10
  hexMapGMSH(1,0,1)= 21
  hexMapGMSH(2,0,1)= 12
  hexMapGMSH(0,1,1)= 23
  hexMapGMSH(1,1,1)= 27
  hexMapGMSH(2,1,1)= 24
  hexMapGMSH(0,2,1)= 18
  hexMapGMSH(1,2,1)= 26
  hexMapGMSH(2,2,1)= 19
  hexMapGMSH(0,0,2)=  1
  hexMapGMSH(1,0,2)=  9
  hexMapGMSH(2,0,2)=  2
  hexMapGMSH(0,1,2)= 11
  hexMapGMSH(1,1,2)= 22
  hexMapGMSH(2,1,2)= 13
  hexMapGMSH(0,2,2)=  5
  hexMapGMSH(1,2,2)= 17
  hexMapGMSH(2,2,2)=  6
CASE(4)
  tetMapGMSH(0,0,0)=  1 
  tetMapGMSH(1,0,0)= 12 
  tetMapGMSH(2,0,0)= 11 
  tetMapGMSH(3,0,0)=  4 
  tetMapGMSH(0,1,0)=  5 
  tetMapGMSH(1,1,0)= 18 
  tetMapGMSH(2,1,0)= 15 
  tetMapGMSH(0,2,0)=  6 
  tetMapGMSH(1,2,0)= 16 
  tetMapGMSH(0,3,0)=  2
  tetMapGMSH(0,0,1)= 10
  tetMapGMSH(1,0,1)= 19
  tetMapGMSH(2,0,1)= 13
  tetMapGMSH(0,1,1)= 17
  tetMapGMSH(1,1,1)= 20
  tetMapGMSH(0,2,1)=  7
  tetMapGMSH(0,0,2)=  9
  tetMapGMSH(1,0,2)= 14
  tetMapGMSH(0,1,2)=  8
  tetMapGMSH(0,0,3)=  3

  hexMapGMSH(0,0,0)=  4
  hexMapGMSH(1,0,0)= 20
  hexMapGMSH(2,0,0)= 19
  hexMapGMSH(3,0,0)=  3
  hexMapGMSH(0,1,0)= 23
  hexMapGMSH(1,1,0)= 50
  hexMapGMSH(2,1,0)= 49
  hexMapGMSH(3,1,0)= 21
  hexMapGMSH(0,2,0)= 24
  hexMapGMSH(1,2,0)= 51
  hexMapGMSH(2,2,0)= 52
  hexMapGMSH(3,2,0)= 22
  hexMapGMSH(0,3,0)=  8
  hexMapGMSH(1,3,0)= 32
  hexMapGMSH(2,3,0)= 31
  hexMapGMSH(3,3,0)=  7
  hexMapGMSH(0,0,1)= 12
  hexMapGMSH(1,0,1)= 34
  hexMapGMSH(2,0,1)= 35
  hexMapGMSH(3,0,1)= 16
  hexMapGMSH(0,1,1)= 44
  hexMapGMSH(1,1,1)= 60
  hexMapGMSH(2,1,1)= 59
  hexMapGMSH(3,1,1)= 46
  hexMapGMSH(0,2,1)= 43
  hexMapGMSH(1,2,1)= 64
  hexMapGMSH(2,2,1)= 63
  hexMapGMSH(3,2,1)= 47
  hexMapGMSH(0,3,1)= 28
  hexMapGMSH(1,3,1)= 56
  hexMapGMSH(2,3,1)= 55
  hexMapGMSH(3,3,1)= 30
  hexMapGMSH(0,0,2)= 11
  hexMapGMSH(1,0,2)= 33
  hexMapGMSH(2,0,2)= 36
  hexMapGMSH(3,0,2)= 15
  hexMapGMSH(0,1,2)= 41
  hexMapGMSH(1,1,2)= 57
  hexMapGMSH(2,1,2)= 58
  hexMapGMSH(3,1,2)= 45
  hexMapGMSH(0,2,2)= 42
  hexMapGMSH(1,2,2)= 61
  hexMapGMSH(2,2,2)= 62
  hexMapGMSH(3,2,2)= 48
  hexMapGMSH(0,3,2)= 27
  hexMapGMSH(1,3,2)= 53
  hexMapGMSH(2,3,2)= 54
  hexMapGMSH(3,3,2)= 29
  hexMapGMSH(0,0,3)=  1
  hexMapGMSH(1,0,3)=  9
  hexMapGMSH(2,0,3)= 10
  hexMapGMSH(3,0,3)=  2
  hexMapGMSH(0,1,3)= 13
  hexMapGMSH(1,1,3)= 37
  hexMapGMSH(2,1,3)= 38
  hexMapGMSH(3,1,3)= 17
  hexMapGMSH(0,2,3)= 14
  hexMapGMSH(1,2,3)= 40
  hexMapGMSH(2,2,3)= 39
  hexMapGMSH(3,2,3)= 18
  hexMapGMSH(0,3,3)=  5
  hexMapGMSH(1,3,3)= 25
  hexMapGMSH(2,3,3)= 26
  hexMapGMSH(3,3,3)=  6
CASE(5)
  hexMapGMSH(0,0,0) = 4
  hexMapGMSH(1,0,0) = 26
  hexMapGMSH(2,0,0) = 25
  hexMapGMSH(3,0,0) = 24
  hexMapGMSH(4,0,0) = 3
  hexMapGMSH(0,1,0) = 30
  hexMapGMSH(1,1,0) = 82
  hexMapGMSH(2,1,0) = 85
  hexMapGMSH(3,1,0) = 81
  hexMapGMSH(4,1,0) = 27
  hexMapGMSH(0,2,0) = 31
  hexMapGMSH(1,2,0) = 86
  hexMapGMSH(2,2,0) = 89
  hexMapGMSH(3,2,0) = 88
  hexMapGMSH(4,2,0) = 28
  hexMapGMSH(0,3,0) = 32
  hexMapGMSH(1,3,0) = 83
  hexMapGMSH(2,3,0) = 87
  hexMapGMSH(3,3,0) = 84
  hexMapGMSH(4,3,0) = 29
  hexMapGMSH(0,4,0) = 8
  hexMapGMSH(1,4,0) = 44
  hexMapGMSH(2,4,0) = 43
  hexMapGMSH(3,4,0) = 42
  hexMapGMSH(4,4,0) = 7
  hexMapGMSH(0,0,1) = 14
  hexMapGMSH(1,0,1) = 46
  hexMapGMSH(2,0,1) = 50
  hexMapGMSH(3,0,1) = 47
  hexMapGMSH(4,0,1) = 20
  hexMapGMSH(0,1,1) = 66
  hexMapGMSH(1,1,1) = 102
  hexMapGMSH(2,1,1) = 112
  hexMapGMSH(3,1,1) = 101
  hexMapGMSH(4,1,1) = 73
  hexMapGMSH(0,2,1) = 69
  hexMapGMSH(1,2,1) = 114
  hexMapGMSH(2,2,1) = 123
  hexMapGMSH(3,2,1) = 113
  hexMapGMSH(4,2,1) = 77
  hexMapGMSH(0,3,1) = 65
  hexMapGMSH(1,3,1) = 106
  hexMapGMSH(2,3,1) = 118
  hexMapGMSH(3,3,1) = 105
  hexMapGMSH(4,3,1) = 74
  hexMapGMSH(0,4,1) = 38
  hexMapGMSH(1,4,1) = 93
  hexMapGMSH(2,4,1) = 96
  hexMapGMSH(3,4,1) = 92
  hexMapGMSH(4,4,1) = 41
  hexMapGMSH(0,0,2) = 13
  hexMapGMSH(1,0,2) = 49
  hexMapGMSH(2,0,2) = 53
  hexMapGMSH(3,0,2) = 51
  hexMapGMSH(4,0,2) = 19
  hexMapGMSH(0,1,2) = 70
  hexMapGMSH(1,1,2) = 108
  hexMapGMSH(2,1,2) = 119
  hexMapGMSH(3,1,2) = 110
  hexMapGMSH(4,1,2) = 76
  hexMapGMSH(0,2,2) = 71
  hexMapGMSH(1,2,2) = 121
  hexMapGMSH(2,2,2) = 125
  hexMapGMSH(3,2,2) = 122
  hexMapGMSH(4,2,2) = 80
  hexMapGMSH(0,3,2) = 68
  hexMapGMSH(1,3,2) = 116
  hexMapGMSH(2,3,2) = 124
  hexMapGMSH(3,3,2) = 117
  hexMapGMSH(4,3,2) = 78
  hexMapGMSH(0,4,2) = 37
  hexMapGMSH(1,4,2) = 97
  hexMapGMSH(2,4,2) = 98
  hexMapGMSH(3,4,2) = 95
  hexMapGMSH(4,4,2) = 40
  hexMapGMSH(0,0,3) = 12
  hexMapGMSH(1,0,3) = 45
  hexMapGMSH(2,0,3) = 52
  hexMapGMSH(3,0,3) = 48
  hexMapGMSH(4,0,3) = 18
  hexMapGMSH(0,1,3) = 63
  hexMapGMSH(1,1,3) = 99
  hexMapGMSH(2,1,3) = 107
  hexMapGMSH(3,1,3) = 100
  hexMapGMSH(4,1,3) = 72
  hexMapGMSH(0,2,3) = 67
  hexMapGMSH(1,2,3) = 109
  hexMapGMSH(2,2,3) = 120
  hexMapGMSH(3,2,3) = 111
  hexMapGMSH(4,2,3) = 79
  hexMapGMSH(0,3,3) = 64
  hexMapGMSH(1,3,3) = 103
  hexMapGMSH(2,3,3) = 115
  hexMapGMSH(3,3,3) = 104
  hexMapGMSH(4,3,3) = 75
  hexMapGMSH(0,4,3) = 36
  hexMapGMSH(1,4,3) = 90
  hexMapGMSH(2,4,3) = 94
  hexMapGMSH(3,4,3) = 91
  hexMapGMSH(4,4,3) = 39
  hexMapGMSH(0,0,4) = 1
  hexMapGMSH(1,0,4) = 9
  hexMapGMSH(2,0,4) = 10
  hexMapGMSH(3,0,4) = 11
  hexMapGMSH(4,0,4) = 2
  hexMapGMSH(0,1,4) = 15
  hexMapGMSH(1,1,4) = 54
  hexMapGMSH(2,1,4) = 58
  hexMapGMSH(3,1,4) = 55
  hexMapGMSH(4,1,4) = 21
  hexMapGMSH(0,2,4) = 16
  hexMapGMSH(1,2,4) = 61
  hexMapGMSH(2,2,4) = 62
  hexMapGMSH(3,2,4) = 59
  hexMapGMSH(4,2,4) = 22
  hexMapGMSH(0,3,4) = 17
  hexMapGMSH(1,3,4) = 57
  hexMapGMSH(2,3,4) = 60
  hexMapGMSH(3,3,4) = 56
  hexMapGMSH(4,3,4) = 23
  hexMapGMSH(0,4,4) = 5
  hexMapGMSH(1,4,4) = 33
  hexMapGMSH(2,4,4) = 34
  hexMapGMSH(3,4,4) = 35
  hexMapGMSH(4,4,4) = 6
CASE(6)
  hexMapGMSH(0,0,0) = 4
  hexMapGMSH(1,0,0) = 32
  hexMapGMSH(2,0,0) = 31
  hexMapGMSH(3,0,0) = 30
  hexMapGMSH(4,0,0) = 29
  hexMapGMSH(5,0,0) = 3
  hexMapGMSH(0,1,0) = 37
  hexMapGMSH(1,1,0) = 122
  hexMapGMSH(2,1,0) = 126
  hexMapGMSH(3,1,0) = 125
  hexMapGMSH(4,1,0) = 121
  hexMapGMSH(5,1,0) = 33
  hexMapGMSH(0,2,0) = 38
  hexMapGMSH(1,2,0) = 127
  hexMapGMSH(2,2,0) = 134
  hexMapGMSH(3,2,0) = 133
  hexMapGMSH(4,2,0) = 132
  hexMapGMSH(5,2,0) = 34
  hexMapGMSH(0,3,0) = 39
  hexMapGMSH(1,3,0) = 128
  hexMapGMSH(2,3,0) = 135
  hexMapGMSH(3,3,0) = 136
  hexMapGMSH(4,3,0) = 131
  hexMapGMSH(5,3,0) = 35
  hexMapGMSH(0,4,0) = 40
  hexMapGMSH(1,4,0) = 123
  hexMapGMSH(2,4,0) = 129
  hexMapGMSH(3,4,0) = 130
  hexMapGMSH(4,4,0) = 124
  hexMapGMSH(5,4,0) = 36
  hexMapGMSH(0,5,0) = 8
  hexMapGMSH(1,5,0) = 56
  hexMapGMSH(2,5,0) = 55
  hexMapGMSH(3,5,0) = 54
  hexMapGMSH(4,5,0) = 53
  hexMapGMSH(5,5,0) = 7
  hexMapGMSH(0,0,1) = 16
  hexMapGMSH(1,0,1) = 58
  hexMapGMSH(2,0,1) = 63
  hexMapGMSH(3,0,1) = 64
  hexMapGMSH(4,0,1) = 59
  hexMapGMSH(5,0,1) = 24
  hexMapGMSH(0,1,1) = 92
  hexMapGMSH(1,1,1) = 156
  hexMapGMSH(2,1,1) = 172
  hexMapGMSH(3,1,1) = 171
  hexMapGMSH(4,1,1) = 155
  hexMapGMSH(5,1,1) = 106
  hexMapGMSH(0,2,1) = 98
  hexMapGMSH(1,2,1) = 175
  hexMapGMSH(2,2,1) = 202
  hexMapGMSH(3,2,1) = 201
  hexMapGMSH(4,2,1) = 173
  hexMapGMSH(5,2,1) = 111
  hexMapGMSH(0,3,1) = 97
  hexMapGMSH(1,3,1) = 176
  hexMapGMSH(2,3,1) = 203
  hexMapGMSH(3,3,1) = 204
  hexMapGMSH(4,3,1) = 174
  hexMapGMSH(5,3,1) = 112
  hexMapGMSH(0,4,1) = 91
  hexMapGMSH(1,4,1) = 160
  hexMapGMSH(2,4,1) = 184
  hexMapGMSH(3,4,1) = 183
  hexMapGMSH(4,4,1) = 159
  hexMapGMSH(5,4,1) = 107
  hexMapGMSH(0,5,1) = 48
  hexMapGMSH(1,5,1) = 140
  hexMapGMSH(2,5,1) = 146
  hexMapGMSH(3,5,1) = 145
  hexMapGMSH(4,5,1) = 139
  hexMapGMSH(5,5,1) = 52
  hexMapGMSH(0,0,2) = 15
  hexMapGMSH(1,0,2) = 62
  hexMapGMSH(2,0,2) = 70
  hexMapGMSH(3,0,2) = 71
  hexMapGMSH(4,0,2) = 65
  hexMapGMSH(5,0,2) = 23
  hexMapGMSH(0,1,2) = 99
  hexMapGMSH(1,1,2) = 164
  hexMapGMSH(2,1,2) = 186
  hexMapGMSH(3,1,2) = 187
  hexMapGMSH(4,1,2) = 168
  hexMapGMSH(5,1,2) = 110
  hexMapGMSH(0,2,2) = 104
  hexMapGMSH(1,2,2) = 196
  hexMapGMSH(2,2,2) = 212
  hexMapGMSH(3,2,2) = 211
  hexMapGMSH(4,2,2) = 198
  hexMapGMSH(5,2,2) = 118
  hexMapGMSH(0,3,2) = 103
  hexMapGMSH(1,3,2) = 195
  hexMapGMSH(2,3,2) = 216
  hexMapGMSH(3,3,2) = 215
  hexMapGMSH(4,3,2) = 199
  hexMapGMSH(5,3,2) = 119
  hexMapGMSH(0,4,2) = 96
  hexMapGMSH(1,4,2) = 180
  hexMapGMSH(2,4,2) = 208
  hexMapGMSH(3,4,2) = 207
  hexMapGMSH(4,4,2) = 182
  hexMapGMSH(5,4,2) = 113
  hexMapGMSH(0,5,2) = 47
  hexMapGMSH(1,5,2) = 147
  hexMapGMSH(2,5,2) = 152
  hexMapGMSH(3,5,2) = 151
  hexMapGMSH(4,5,2) = 144
  hexMapGMSH(5,5,2) = 51
  hexMapGMSH(0,0,3) = 14
  hexMapGMSH(1,0,3) = 61
  hexMapGMSH(2,0,3) = 69
  hexMapGMSH(3,0,3) = 72
  hexMapGMSH(4,0,3) = 66
  hexMapGMSH(5,0,3) = 22
  hexMapGMSH(0,1,3) = 100
  hexMapGMSH(1,1,3) = 163
  hexMapGMSH(2,1,3) = 185
  hexMapGMSH(3,1,3) = 188
  hexMapGMSH(4,1,3) = 167
  hexMapGMSH(5,1,3) = 109
  hexMapGMSH(0,2,3) = 101
  hexMapGMSH(1,2,3) = 193
  hexMapGMSH(2,2,3) = 209
  hexMapGMSH(3,2,3) = 210
  hexMapGMSH(4,2,3) = 197
  hexMapGMSH(5,2,3) = 117
  hexMapGMSH(0,3,3) = 102
  hexMapGMSH(1,3,3) = 194
  hexMapGMSH(2,3,3) = 213
  hexMapGMSH(3,3,3) = 214
  hexMapGMSH(4,3,3) = 200
  hexMapGMSH(5,3,3) = 120
  hexMapGMSH(0,4,3) = 95
  hexMapGMSH(1,4,3) = 179
  hexMapGMSH(2,4,3) = 205
  hexMapGMSH(3,4,3) = 206
  hexMapGMSH(4,4,3) = 181
  hexMapGMSH(5,4,3) = 114
  hexMapGMSH(0,5,3) = 46
  hexMapGMSH(1,5,3) = 148
  hexMapGMSH(2,5,3) = 149
  hexMapGMSH(3,5,3) = 150
  hexMapGMSH(4,5,3) = 143
  hexMapGMSH(5,5,3) = 50
  hexMapGMSH(0,0,4) = 13
  hexMapGMSH(1,0,4) = 57
  hexMapGMSH(2,0,4) = 68
  hexMapGMSH(3,0,4) = 67
  hexMapGMSH(4,0,4) = 60
  hexMapGMSH(5,0,4) = 21
  hexMapGMSH(0,1,4) = 89
  hexMapGMSH(1,1,4) = 153
  hexMapGMSH(2,1,4) = 161
  hexMapGMSH(3,1,4) = 162
  hexMapGMSH(4,1,4) = 154
  hexMapGMSH(5,1,4) = 105
  hexMapGMSH(0,2,4) = 93
  hexMapGMSH(1,2,4) = 165
  hexMapGMSH(2,2,4) = 189
  hexMapGMSH(3,2,4) = 190
  hexMapGMSH(4,2,4) = 169
  hexMapGMSH(5,2,4) = 116
  hexMapGMSH(0,3,4) = 94
  hexMapGMSH(1,3,4) = 166
  hexMapGMSH(2,3,4) = 192
  hexMapGMSH(3,3,4) = 191
  hexMapGMSH(4,3,4) = 170
  hexMapGMSH(5,3,4) = 115
  hexMapGMSH(0,4,4) = 90
  hexMapGMSH(1,4,4) = 157
  hexMapGMSH(2,4,4) = 177
  hexMapGMSH(3,4,4) = 178
  hexMapGMSH(4,4,4) = 158
  hexMapGMSH(5,4,4) = 108
  hexMapGMSH(0,5,4) = 45
  hexMapGMSH(1,5,4) = 137
  hexMapGMSH(2,5,4) = 141
  hexMapGMSH(3,5,4) = 142
  hexMapGMSH(4,5,4) = 138
  hexMapGMSH(5,5,4) = 49
  hexMapGMSH(0,0,5) = 1
  hexMapGMSH(1,0,5) = 9
  hexMapGMSH(2,0,5) = 10
  hexMapGMSH(3,0,5) = 11
  hexMapGMSH(4,0,5) = 12
  hexMapGMSH(5,0,5) = 2
  hexMapGMSH(0,1,5) = 17
  hexMapGMSH(1,1,5) = 73
  hexMapGMSH(2,1,5) = 77
  hexMapGMSH(3,1,5) = 78
  hexMapGMSH(4,1,5) = 74
  hexMapGMSH(5,1,5) = 25
  hexMapGMSH(0,2,5) = 18
  hexMapGMSH(1,2,5) = 84
  hexMapGMSH(2,2,5) = 85
  hexMapGMSH(3,2,5) = 86
  hexMapGMSH(4,2,5) = 79
  hexMapGMSH(5,2,5) = 26
  hexMapGMSH(0,3,5) = 19
  hexMapGMSH(1,3,5) = 83
  hexMapGMSH(2,3,5) = 88
  hexMapGMSH(3,3,5) = 87
  hexMapGMSH(4,3,5) = 80
  hexMapGMSH(5,3,5) = 27
  hexMapGMSH(0,4,5) = 20
  hexMapGMSH(1,4,5) = 76
  hexMapGMSH(2,4,5) = 82
  hexMapGMSH(3,4,5) = 81
  hexMapGMSH(4,4,5) = 75
  hexMapGMSH(5,4,5) = 28
  hexMapGMSH(0,5,5) = 5
  hexMapGMSH(1,5,5) = 41
  hexMapGMSH(2,5,5) = 42
  hexMapGMSH(3,5,5) = 43
  hexMapGMSH(4,5,5) = 44
  hexMapGMSH(5,5,5) = 6
CASE(8)
  hexMapGMSH(0,0,0) = 4
  hexMapGMSH(1,0,0) = 44
  hexMapGMSH(2,0,0) = 43
  hexMapGMSH(3,0,0) = 42
  hexMapGMSH(4,0,0) = 41
  hexMapGMSH(5,0,0) = 40
  hexMapGMSH(6,0,0) = 39
  hexMapGMSH(7,0,0) = 3
  hexMapGMSH(0,1,0) = 51
  hexMapGMSH(1,1,0) = 226
  hexMapGMSH(2,1,0) = 232
  hexMapGMSH(3,1,0) = 231
  hexMapGMSH(4,1,0) = 230
  hexMapGMSH(5,1,0) = 229
  hexMapGMSH(6,1,0) = 225
  hexMapGMSH(7,1,0) = 45
  hexMapGMSH(0,2,0) = 52
  hexMapGMSH(1,2,0) = 233
  hexMapGMSH(2,2,0) = 246
  hexMapGMSH(3,2,0) = 250
  hexMapGMSH(4,2,0) = 249
  hexMapGMSH(5,2,0) = 245
  hexMapGMSH(6,2,0) = 244
  hexMapGMSH(7,2,0) = 46
  hexMapGMSH(0,3,0) = 53
  hexMapGMSH(1,3,0) = 234
  hexMapGMSH(2,3,0) = 251
  hexMapGMSH(3,3,0) = 258
  hexMapGMSH(4,3,0) = 257
  hexMapGMSH(5,3,0) = 256
  hexMapGMSH(6,3,0) = 243
  hexMapGMSH(7,3,0) = 47
  hexMapGMSH(0,4,0) = 54
  hexMapGMSH(1,4,0) = 235
  hexMapGMSH(2,4,0) = 252
  hexMapGMSH(3,4,0) = 259
  hexMapGMSH(4,4,0) = 260
  hexMapGMSH(5,4,0) = 255
  hexMapGMSH(6,4,0) = 242
  hexMapGMSH(7,4,0) = 48
  hexMapGMSH(0,5,0) = 55
  hexMapGMSH(1,5,0) = 236
  hexMapGMSH(2,5,0) = 247
  hexMapGMSH(3,5,0) = 253
  hexMapGMSH(4,5,0) = 254
  hexMapGMSH(5,5,0) = 248
  hexMapGMSH(6,5,0) = 241
  hexMapGMSH(7,5,0) = 49
  hexMapGMSH(0,6,0) = 56
  hexMapGMSH(1,6,0) = 227
  hexMapGMSH(2,6,0) = 237
  hexMapGMSH(3,6,0) = 238
  hexMapGMSH(4,6,0) = 239
  hexMapGMSH(5,6,0) = 240
  hexMapGMSH(6,6,0) = 228
  hexMapGMSH(7,6,0) = 50
  hexMapGMSH(0,7,0) = 8
  hexMapGMSH(1,7,0) = 80
  hexMapGMSH(2,7,0) = 79
  hexMapGMSH(3,7,0) = 78
  hexMapGMSH(4,7,0) = 77
  hexMapGMSH(5,7,0) = 76
  hexMapGMSH(6,7,0) = 75
  hexMapGMSH(7,7,0) = 7
  hexMapGMSH(0,0,1) = 20
  hexMapGMSH(1,0,1) = 82
  hexMapGMSH(2,0,1) = 89
  hexMapGMSH(3,0,1) = 90
  hexMapGMSH(4,0,1) = 91
  hexMapGMSH(5,0,1) = 92
  hexMapGMSH(6,0,1) = 83
  hexMapGMSH(7,0,1) = 32
  hexMapGMSH(0,1,1) = 156
  hexMapGMSH(1,1,1) = 300
  hexMapGMSH(2,1,1) = 328
  hexMapGMSH(3,1,1) = 327
  hexMapGMSH(4,1,1) = 326
  hexMapGMSH(5,1,1) = 325
  hexMapGMSH(6,1,1) = 299
  hexMapGMSH(7,1,1) = 190
  hexMapGMSH(0,2,1) = 168
  hexMapGMSH(1,2,1) = 333
  hexMapGMSH(2,2,1) = 418
  hexMapGMSH(3,2,1) = 422
  hexMapGMSH(4,2,1) = 421
  hexMapGMSH(5,2,1) = 417
  hexMapGMSH(6,2,1) = 329
  hexMapGMSH(7,2,1) = 197
  hexMapGMSH(0,3,1) = 167
  hexMapGMSH(1,3,1) = 334
  hexMapGMSH(2,3,1) = 423
  hexMapGMSH(3,3,1) = 430
  hexMapGMSH(4,3,1) = 429
  hexMapGMSH(5,3,1) = 428
  hexMapGMSH(6,3,1) = 330
  hexMapGMSH(7,3,1) = 198
  hexMapGMSH(0,4,1) = 166
  hexMapGMSH(1,4,1) = 335
  hexMapGMSH(2,4,1) = 424
  hexMapGMSH(3,4,1) = 431
  hexMapGMSH(4,4,1) = 432
  hexMapGMSH(5,4,1) = 427
  hexMapGMSH(6,4,1) = 331
  hexMapGMSH(7,4,1) = 199
  hexMapGMSH(0,5,1) = 165
  hexMapGMSH(1,5,1) = 336
  hexMapGMSH(2,5,1) = 419
  hexMapGMSH(3,5,1) = 425
  hexMapGMSH(4,5,1) = 426
  hexMapGMSH(5,5,1) = 420
  hexMapGMSH(6,5,1) = 332
  hexMapGMSH(7,5,1) = 200
  hexMapGMSH(0,6,1) = 155
  hexMapGMSH(1,6,1) = 304
  hexMapGMSH(2,6,1) = 352
  hexMapGMSH(3,6,1) = 351
  hexMapGMSH(4,6,1) = 350
  hexMapGMSH(5,6,1) = 349
  hexMapGMSH(6,6,1) = 303
  hexMapGMSH(7,6,1) = 191
  hexMapGMSH(0,7,1) = 68
  hexMapGMSH(1,7,1) = 264
  hexMapGMSH(2,7,1) = 276
  hexMapGMSH(3,7,1) = 275
  hexMapGMSH(4,7,1) = 274
  hexMapGMSH(5,7,1) = 273
  hexMapGMSH(6,7,1) = 263
  hexMapGMSH(7,7,1) = 74
  hexMapGMSH(0,0,2) = 19
  hexMapGMSH(1,0,2) = 88
  hexMapGMSH(2,0,2) = 102
  hexMapGMSH(3,0,2) = 107
  hexMapGMSH(4,0,2) = 108
  hexMapGMSH(5,0,2) = 103
  hexMapGMSH(6,0,2) = 93
  hexMapGMSH(7,0,2) = 31
  hexMapGMSH(0,1,2) = 169
  hexMapGMSH(1,1,2) = 312
  hexMapGMSH(2,1,2) = 354
  hexMapGMSH(3,1,2) = 359
  hexMapGMSH(4,1,2) = 360
  hexMapGMSH(5,1,2) = 355
  hexMapGMSH(6,1,2) = 320
  hexMapGMSH(7,1,2) = 196
  hexMapGMSH(0,2,2) = 176
  hexMapGMSH(1,2,2) = 388
  hexMapGMSH(2,2,2) = 452
  hexMapGMSH(3,2,2) = 468
  hexMapGMSH(4,2,2) = 467
  hexMapGMSH(5,2,2) = 451
  hexMapGMSH(6,2,2) = 402
  hexMapGMSH(7,2,2) = 210
  hexMapGMSH(0,3,2) = 182
  hexMapGMSH(1,3,2) = 394
  hexMapGMSH(2,3,2) = 471
  hexMapGMSH(3,3,2) = 498
  hexMapGMSH(4,3,2) = 497
  hexMapGMSH(5,3,2) = 469
  hexMapGMSH(6,3,2) = 407
  hexMapGMSH(7,3,2) = 215
  hexMapGMSH(0,4,2) = 181
  hexMapGMSH(1,4,2) = 393
  hexMapGMSH(2,4,2) = 472
  hexMapGMSH(3,4,2) = 499
  hexMapGMSH(4,4,2) = 500
  hexMapGMSH(5,4,2) = 470
  hexMapGMSH(6,4,2) = 408
  hexMapGMSH(7,4,2) = 216
  hexMapGMSH(0,5,2) = 175
  hexMapGMSH(1,5,2) = 387
  hexMapGMSH(2,5,2) = 456
  hexMapGMSH(3,5,2) = 480
  hexMapGMSH(4,5,2) = 479
  hexMapGMSH(5,5,2) = 455
  hexMapGMSH(6,5,2) = 403
  hexMapGMSH(7,5,2) = 211
  hexMapGMSH(0,6,2) = 164
  hexMapGMSH(1,6,2) = 344
  hexMapGMSH(2,6,2) = 436
  hexMapGMSH(3,6,2) = 442
  hexMapGMSH(4,6,2) = 441
  hexMapGMSH(5,6,2) = 435
  hexMapGMSH(6,6,2) = 348
  hexMapGMSH(7,6,2) = 201
  hexMapGMSH(0,7,2) = 67
  hexMapGMSH(1,7,2) = 277
  hexMapGMSH(2,7,2) = 284
  hexMapGMSH(3,7,2) = 290
  hexMapGMSH(4,7,2) = 289
  hexMapGMSH(5,7,2) = 283
  hexMapGMSH(6,7,2) = 272
  hexMapGMSH(7,7,2) = 73
  hexMapGMSH(0,0,3) = 18
  hexMapGMSH(1,0,3) = 87
  hexMapGMSH(2,0,3) = 106
  hexMapGMSH(3,0,3) = 114
  hexMapGMSH(4,0,3) = 115
  hexMapGMSH(5,0,3) = 109
  hexMapGMSH(6,0,3) = 94
  hexMapGMSH(7,0,3) = 30
  hexMapGMSH(0,1,3) = 170
  hexMapGMSH(1,1,3) = 311
  hexMapGMSH(2,1,3) = 358
  hexMapGMSH(3,1,3) = 366
  hexMapGMSH(4,1,3) = 367
  hexMapGMSH(5,1,3) = 361
  hexMapGMSH(6,1,3) = 319
  hexMapGMSH(7,1,3) = 195
  hexMapGMSH(0,2,3) = 183
  hexMapGMSH(1,2,3) = 395
  hexMapGMSH(2,2,3) = 460
  hexMapGMSH(3,2,3) = 482
  hexMapGMSH(4,2,3) = 483
  hexMapGMSH(5,2,3) = 464
  hexMapGMSH(6,2,3) = 406
  hexMapGMSH(7,2,3) = 214
  hexMapGMSH(0,3,3) = 188
  hexMapGMSH(1,3,3) = 400
  hexMapGMSH(2,3,3) = 492
  hexMapGMSH(3,3,3) = 508
  hexMapGMSH(4,3,3) = 507
  hexMapGMSH(5,3,3) = 494
  hexMapGMSH(6,3,3) = 414
  hexMapGMSH(7,3,3) = 222
  hexMapGMSH(0,4,3) = 187
  hexMapGMSH(1,4,3) = 399
  hexMapGMSH(2,4,3) = 491
  hexMapGMSH(3,4,3) = 512
  hexMapGMSH(4,4,3) = 511
  hexMapGMSH(5,4,3) = 495
  hexMapGMSH(6,4,3) = 415
  hexMapGMSH(7,4,3) = 223
  hexMapGMSH(0,5,3) = 180
  hexMapGMSH(1,5,3) = 392
  hexMapGMSH(2,5,3) = 476
  hexMapGMSH(3,5,3) = 504
  hexMapGMSH(4,5,3) = 503
  hexMapGMSH(5,5,3) = 478
  hexMapGMSH(6,5,3) = 409
  hexMapGMSH(7,5,3) = 217
  hexMapGMSH(0,6,3) = 163
  hexMapGMSH(1,6,3) = 343
  hexMapGMSH(2,6,3) = 443
  hexMapGMSH(3,6,3) = 448
  hexMapGMSH(4,6,3) = 447
  hexMapGMSH(5,6,3) = 440
  hexMapGMSH(6,6,3) = 347
  hexMapGMSH(7,6,3) = 202
  hexMapGMSH(0,7,3) = 66
  hexMapGMSH(1,7,3) = 278
  hexMapGMSH(2,7,3) = 291
  hexMapGMSH(3,7,3) = 296
  hexMapGMSH(4,7,3) = 295
  hexMapGMSH(5,7,3) = 288
  hexMapGMSH(6,7,3) = 271
  hexMapGMSH(7,7,3) = 72
  hexMapGMSH(0,0,4) = 17
  hexMapGMSH(1,0,4) = 86
  hexMapGMSH(2,0,4) = 105
  hexMapGMSH(3,0,4) = 113
  hexMapGMSH(4,0,4) = 116
  hexMapGMSH(5,0,4) = 110
  hexMapGMSH(6,0,4) = 95
  hexMapGMSH(7,0,4) = 29
  hexMapGMSH(0,1,4) = 171
  hexMapGMSH(1,1,4) = 310
  hexMapGMSH(2,1,4) = 357
  hexMapGMSH(3,1,4) = 365
  hexMapGMSH(4,1,4) = 368
  hexMapGMSH(5,1,4) = 362
  hexMapGMSH(6,1,4) = 318
  hexMapGMSH(7,1,4) = 194
  hexMapGMSH(0,2,4) = 184
  hexMapGMSH(1,2,4) = 396
  hexMapGMSH(2,2,4) = 459
  hexMapGMSH(3,2,4) = 481
  hexMapGMSH(4,2,4) = 484
  hexMapGMSH(5,2,4) = 463
  hexMapGMSH(6,2,4) = 405
  hexMapGMSH(7,2,4) = 213
  hexMapGMSH(0,3,4) = 185
  hexMapGMSH(1,3,4) = 397
  hexMapGMSH(2,3,4) = 489
  hexMapGMSH(3,3,4) = 505
  hexMapGMSH(4,3,4) = 506
  hexMapGMSH(5,3,4) = 493
  hexMapGMSH(6,3,4) = 413
  hexMapGMSH(7,3,4) = 221
  hexMapGMSH(0,4,4) = 186
  hexMapGMSH(1,4,4) = 398
  hexMapGMSH(2,4,4) = 490
  hexMapGMSH(3,4,4) = 509
  hexMapGMSH(4,4,4) = 510
  hexMapGMSH(5,4,4) = 496
  hexMapGMSH(6,4,4) = 416
  hexMapGMSH(7,4,4) = 224
  hexMapGMSH(0,5,4) = 179
  hexMapGMSH(1,5,4) = 391
  hexMapGMSH(2,5,4) = 475
  hexMapGMSH(3,5,4) = 501
  hexMapGMSH(4,5,4) = 502
  hexMapGMSH(5,5,4) = 477
  hexMapGMSH(6,5,4) = 410
  hexMapGMSH(7,5,4) = 218
  hexMapGMSH(0,6,4) = 162
  hexMapGMSH(1,6,4) = 342
  hexMapGMSH(2,6,4) = 444
  hexMapGMSH(3,6,4) = 445
  hexMapGMSH(4,6,4) = 446
  hexMapGMSH(5,6,4) = 439
  hexMapGMSH(6,6,4) = 346
  hexMapGMSH(7,6,4) = 203
  hexMapGMSH(0,7,4) = 65
  hexMapGMSH(1,7,4) = 279
  hexMapGMSH(2,7,4) = 292
  hexMapGMSH(3,7,4) = 293
  hexMapGMSH(4,7,4) = 294
  hexMapGMSH(5,7,4) = 287
  hexMapGMSH(6,7,4) = 270
  hexMapGMSH(7,7,4) = 71
  hexMapGMSH(0,0,5) = 16
  hexMapGMSH(1,0,5) = 85
  hexMapGMSH(2,0,5) = 101
  hexMapGMSH(3,0,5) = 112
  hexMapGMSH(4,0,5) = 111
  hexMapGMSH(5,0,5) = 104
  hexMapGMSH(6,0,5) = 96
  hexMapGMSH(7,0,5) = 28
  hexMapGMSH(0,1,5) = 172
  hexMapGMSH(1,1,5) = 309
  hexMapGMSH(2,1,5) = 353
  hexMapGMSH(3,1,5) = 364
  hexMapGMSH(4,1,5) = 363
  hexMapGMSH(5,1,5) = 356
  hexMapGMSH(6,1,5) = 317
  hexMapGMSH(7,1,5) = 193
  hexMapGMSH(0,2,5) = 173
  hexMapGMSH(1,2,5) = 385
  hexMapGMSH(2,2,5) = 449
  hexMapGMSH(3,2,5) = 457
  hexMapGMSH(4,2,5) = 458
  hexMapGMSH(5,2,5) = 450
  hexMapGMSH(6,2,5) = 401
  hexMapGMSH(7,2,5) = 209
  hexMapGMSH(0,3,5) = 177
  hexMapGMSH(1,3,5) = 389
  hexMapGMSH(2,3,5) = 461
  hexMapGMSH(3,3,5) = 485
  hexMapGMSH(4,3,5) = 486
  hexMapGMSH(5,3,5) = 465
  hexMapGMSH(6,3,5) = 412
  hexMapGMSH(7,3,5) = 220
  hexMapGMSH(0,4,5) = 178
  hexMapGMSH(1,4,5) = 390
  hexMapGMSH(2,4,5) = 462
  hexMapGMSH(3,4,5) = 488
  hexMapGMSH(4,4,5) = 487
  hexMapGMSH(5,4,5) = 466
  hexMapGMSH(6,4,5) = 411
  hexMapGMSH(7,4,5) = 219
  hexMapGMSH(0,5,5) = 174
  hexMapGMSH(1,5,5) = 386
  hexMapGMSH(2,5,5) = 453
  hexMapGMSH(3,5,5) = 473
  hexMapGMSH(4,5,5) = 474
  hexMapGMSH(5,5,5) = 454
  hexMapGMSH(6,5,5) = 404
  hexMapGMSH(7,5,5) = 212
  hexMapGMSH(0,6,5) = 161
  hexMapGMSH(1,6,5) = 341
  hexMapGMSH(2,6,5) = 433
  hexMapGMSH(3,6,5) = 437
  hexMapGMSH(4,6,5) = 438
  hexMapGMSH(5,6,5) = 434
  hexMapGMSH(6,6,5) = 345
  hexMapGMSH(7,6,5) = 204
  hexMapGMSH(0,7,5) = 64
  hexMapGMSH(1,7,5) = 280
  hexMapGMSH(2,7,5) = 281
  hexMapGMSH(3,7,5) = 285
  hexMapGMSH(4,7,5) = 286
  hexMapGMSH(5,7,5) = 282
  hexMapGMSH(6,7,5) = 269
  hexMapGMSH(7,7,5) = 70
  hexMapGMSH(0,0,6) = 15
  hexMapGMSH(1,0,6) = 81
  hexMapGMSH(2,0,6) = 100
  hexMapGMSH(3,0,6) = 99
  hexMapGMSH(4,0,6) = 98
  hexMapGMSH(5,0,6) = 97
  hexMapGMSH(6,0,6) = 84
  hexMapGMSH(7,0,6) = 27
  hexMapGMSH(0,1,6) = 153
  hexMapGMSH(1,1,6) = 297
  hexMapGMSH(2,1,6) = 305
  hexMapGMSH(3,1,6) = 306
  hexMapGMSH(4,1,6) = 307
  hexMapGMSH(5,1,6) = 308
  hexMapGMSH(6,1,6) = 298
  hexMapGMSH(7,1,6) = 189
  hexMapGMSH(0,2,6) = 157
  hexMapGMSH(1,2,6) = 313
  hexMapGMSH(2,2,6) = 369
  hexMapGMSH(3,2,6) = 373
  hexMapGMSH(4,2,6) = 374
  hexMapGMSH(5,2,6) = 370
  hexMapGMSH(6,2,6) = 321
  hexMapGMSH(7,2,6) = 208
  hexMapGMSH(0,3,6) = 158
  hexMapGMSH(1,3,6) = 314
  hexMapGMSH(2,3,6) = 380
  hexMapGMSH(3,3,6) = 381
  hexMapGMSH(4,3,6) = 382
  hexMapGMSH(5,3,6) = 375
  hexMapGMSH(6,3,6) = 322
  hexMapGMSH(7,3,6) = 207
  hexMapGMSH(0,4,6) = 159
  hexMapGMSH(1,4,6) = 315
  hexMapGMSH(2,4,6) = 379
  hexMapGMSH(3,4,6) = 384
  hexMapGMSH(4,4,6) = 383
  hexMapGMSH(5,4,6) = 376
  hexMapGMSH(6,4,6) = 323
  hexMapGMSH(7,4,6) = 206
  hexMapGMSH(0,5,6) = 160
  hexMapGMSH(1,5,6) = 316
  hexMapGMSH(2,5,6) = 372
  hexMapGMSH(3,5,6) = 378
  hexMapGMSH(4,5,6) = 377
  hexMapGMSH(5,5,6) = 371
  hexMapGMSH(6,5,6) = 324
  hexMapGMSH(7,5,6) = 205
  hexMapGMSH(0,6,6) = 154
  hexMapGMSH(1,6,6) = 301
  hexMapGMSH(2,6,6) = 337
  hexMapGMSH(3,6,6) = 338
  hexMapGMSH(4,6,6) = 339
  hexMapGMSH(5,6,6) = 340
  hexMapGMSH(6,6,6) = 302
  hexMapGMSH(7,6,6) = 192
  hexMapGMSH(0,7,6) = 63
  hexMapGMSH(1,7,6) = 261
  hexMapGMSH(2,7,6) = 265
  hexMapGMSH(3,7,6) = 266
  hexMapGMSH(4,7,6) = 267
  hexMapGMSH(5,7,6) = 268
  hexMapGMSH(6,7,6) = 262
  hexMapGMSH(7,7,6) = 69
  hexMapGMSH(0,0,7) = 1
  hexMapGMSH(1,0,7) = 9
  hexMapGMSH(2,0,7) = 10
  hexMapGMSH(3,0,7) = 11
  hexMapGMSH(4,0,7) = 12
  hexMapGMSH(5,0,7) = 13
  hexMapGMSH(6,0,7) = 14
  hexMapGMSH(7,0,7) = 2
  hexMapGMSH(0,1,7) = 21
  hexMapGMSH(1,1,7) = 117
  hexMapGMSH(2,1,7) = 121
  hexMapGMSH(3,1,7) = 122
  hexMapGMSH(4,1,7) = 123
  hexMapGMSH(5,1,7) = 124
  hexMapGMSH(6,1,7) = 118
  hexMapGMSH(7,1,7) = 33
  hexMapGMSH(0,2,7) = 22
  hexMapGMSH(1,2,7) = 136
  hexMapGMSH(2,2,7) = 137
  hexMapGMSH(3,2,7) = 141
  hexMapGMSH(4,2,7) = 142
  hexMapGMSH(5,2,7) = 138
  hexMapGMSH(6,2,7) = 125
  hexMapGMSH(7,2,7) = 34
  hexMapGMSH(0,3,7) = 23
  hexMapGMSH(1,3,7) = 135
  hexMapGMSH(2,3,7) = 148
  hexMapGMSH(3,3,7) = 149
  hexMapGMSH(4,3,7) = 150
  hexMapGMSH(5,3,7) = 143
  hexMapGMSH(6,3,7) = 126
  hexMapGMSH(7,3,7) = 35
  hexMapGMSH(0,4,7) = 24
  hexMapGMSH(1,4,7) = 134
  hexMapGMSH(2,4,7) = 147
  hexMapGMSH(3,4,7) = 152
  hexMapGMSH(4,4,7) = 151
  hexMapGMSH(5,4,7) = 144
  hexMapGMSH(6,4,7) = 127
  hexMapGMSH(7,4,7) = 36
  hexMapGMSH(0,5,7) = 25
  hexMapGMSH(1,5,7) = 133
  hexMapGMSH(2,5,7) = 140
  hexMapGMSH(3,5,7) = 146
  hexMapGMSH(4,5,7) = 145
  hexMapGMSH(5,5,7) = 139
  hexMapGMSH(6,5,7) = 128
  hexMapGMSH(7,5,7) = 37
  hexMapGMSH(0,6,7) = 26
  hexMapGMSH(1,6,7) = 120
  hexMapGMSH(2,6,7) = 132
  hexMapGMSH(3,6,7) = 131
  hexMapGMSH(4,6,7) = 130
  hexMapGMSH(5,6,7) = 129
  hexMapGMSH(6,6,7) = 119
  hexMapGMSH(7,6,7) = 38
  hexMapGMSH(0,7,7) = 5
  hexMapGMSH(1,7,7) = 57
  hexMapGMSH(2,7,7) = 58
  hexMapGMSH(3,7,7) = 59
  hexMapGMSH(4,7,7) = 60
  hexMapGMSH(5,7,7) = 61
  hexMapGMSH(6,7,7) = 62
  hexMapGMSH(7,7,7) = 6

CASE DEFAULT
  CALL abort(__STAMP__,&
             'Elements of specified or higher order are not implemented yet. Order: ',bOrd)
END SELECT
END SUBROUTINE getGMSHVolumeMapping


END MODULE MOD_Readin_GMSH_Vars
