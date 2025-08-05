!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!Oct 2024: Alireza Moghaddasi created for Noah-MP.4.0.1 GPS displacement
#include "LIS_misc.h"
module noahmp401_dagpsdisp_Mod
!BOP
!
! !MODULE: noahmp401_dagpsdisp_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:
!
! !USES:        
  use ESMF
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: noahmp401_dagpsdisp_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: noahmp401_dasm_gpsdisp_struc
!EOP

 type, public :: dasm_gpsdisp_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type dasm_gpsdisp_dec

  type(dasm_gpsdisp_dec), allocatable :: noahmp401_dasm_gpsdisp_struc(:)

contains
!BOP
! 
! !ROUTINE: noahmp401_dagpsdisp_init
! \label{noahmp401_dagpsdisp_init}
! 
! !INTERFACE:
  subroutine noahmp401_dagpsdisp_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP

    implicit none
    integer                :: k
    integer                :: n
    integer                :: status
    integer                :: ngrid

    if(.not.allocated(noahmp401_dasm_gpsdisp_struc)) then
       allocate(noahmp401_dasm_gpsdisp_struc(LIS_rc%nnest))
    endif

  end subroutine noahmp401_dagpsdisp_init
end module noahmp401_dagpsdisp_Mod
