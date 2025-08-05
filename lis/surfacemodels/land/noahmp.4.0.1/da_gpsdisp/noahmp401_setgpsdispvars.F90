!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_setgpsdispvars
!  \label{noahmp401_setgpsdispvars}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp401_setgpsdispvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp401_lsmMod
  use module_sf_noahlsm_36
  use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE,SMCWLT_TABLE
  use noahmp401_gpsdisp_DAlogMod, only : noahmp401_gpsdisp_DAlog

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP
  real, parameter        :: MIN_GWS_THRESHOLD = 0.00
  real, parameter        :: MAX_GWS_THRESHOLD = 7000.0
  real, parameter        :: MAX_WA = 7000.0
  real, parameter        :: ZSOIL = 2 !mm
  real, parameter        :: ROUS = 0.2 ! specific yield
  real                   :: MIN_THRESHOLD 
  real                   :: MAX_THRESHOLD
  real                   :: sm_threshold
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  type(ESMF_Field)       :: gwField
  type(ESMF_Field)       :: sweField

  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real, pointer          :: gws(:)
  real, pointer          :: swe(:)
  integer                :: t, gid
  integer                :: status
  real                   :: delta1
  real                   :: delta
  logical                :: update_flag(LIS_rc%ngrid(n))
  integer                :: SOILTYP           ! soil type index [-]
  real                   :: dsneqv,dsnowh,swe_old
  
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in noahmp401_setgpsdispvars")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 2 failed in noahmp401_setgpsdispvars")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 3 failed in noahmp401_setgpsdispvars")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 4 failed in noahmp401_setgpsdispvars")
  call ESMF_StateGet(LSM_State,"Groundwater Storage",gwField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp401_setgpsdispvars")
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: SWE failed in noahmp401_setgpsdispvars")


  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noahmp401_setgpsdispvars")
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in noahmp401_setgpsdispvars")
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in noahmp401_setgpsdispvars")
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in noahmp401_setgpsdispvars")
  call ESMF_FieldGet(gwField,localDE=0,farrayPtr=gws,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Groundwater Storage failed in noahmp401_setgpsdispvars")
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: SWE failed in noahmp401_setgpsdispvars")


  update_flag = .true. 

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
  
     SOILTYP = NOAHMP401_struc(n)%noahmp401(t)%soiltype        
     MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP)
     MIN_THRESHOLD = 0.02  ! Use same as GRACE DA
     sm_threshold = SMCMAX_TABLE(SOILTYP) - 0.02
     
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row) 
     
     delta1 = soilm1(t)-NOAHMP401_struc(n)%noahmp401(t)%smc(1)

     if(NOAHMP401_struc(n)%noahmp401(t)%sh2o(1)+delta1.gt.MIN_THRESHOLD .and.&
          NOAHMP401_struc(n)%noahmp401(t)%sh2o(1)+delta1.lt.&
          sm_threshold) then 
        update_flag(gid) = update_flag(gid).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
     endif
     
     delta1 = soilm2(t)-NOAHMP401_struc(n)%noahmp401(t)%smc(2)

     if(NOAHMP401_struc(n)%noahmp401(t)%sh2o(2)+delta1.gt.MIN_THRESHOLD .and.&
          NOAHMP401_struc(n)%noahmp401(t)%sh2o(2)+delta1.lt.&
          sm_threshold) then
        update_flag(gid) = update_flag(gid).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
     endif

     delta1 = soilm3(t)-NOAHMP401_struc(n)%noahmp401(t)%smc(3)

     if(NOAHMP401_struc(n)%noahmp401(t)%sh2o(3)+delta1.gt.MIN_THRESHOLD .and.&
          NOAHMP401_struc(n)%noahmp401(t)%sh2o(3)+delta1.lt.&
          sm_threshold) then
        update_flag(gid) = update_flag(gid).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
     endif

     delta1 = soilm4(t)-NOAHMP401_struc(n)%noahmp401(t)%smc(4)

     if(NOAHMP401_struc(n)%noahmp401(t)%sh2o(4)+delta1.gt.MIN_THRESHOLD .and.&
          NOAHMP401_struc(n)%noahmp401(t)%sh2o(4)+delta1.lt.&
          sm_threshold) then
        update_flag(gid) = update_flag(gid).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
     endif
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     if(update_flag(gid)) then
        delta = soilm1(t) - NOAHMP401_struc(n)%noahmp401(t)%smc(1)
              NOAHMP401_struc(n)%noahmp401(t)%smc(1) = soilm1(t)
        NOAHMP401_struc(n)%noahmp401(t)%sh2o(1) = &
             NOAHMP401_struc(n)%noahmp401(t)%sh2o(1) + delta

        delta = soilm2(t) - NOAHMP401_struc(n)%noahmp401(t)%smc(2)
                 NOAHMP401_struc(n)%noahmp401(t)%smc(2) = soilm2(t)
        NOAHMP401_struc(n)%noahmp401(t)%sh2o(2) = &
             NOAHMP401_struc(n)%noahmp401(t)%sh2o(2) + delta

        delta = soilm3(t) - NOAHMP401_struc(n)%noahmp401(t)%smc(3)
                 NOAHMP401_struc(n)%noahmp401(t)%smc(3) = soilm3(t)
        NOAHMP401_struc(n)%noahmp401(t)%sh2o(3) = &
             NOAHMP401_struc(n)%noahmp401(t)%sh2o(3) + delta

        delta = soilm4(t) - NOAHMP401_struc(n)%noahmp401(t)%smc(4)
                 NOAHMP401_struc(n)%noahmp401(t)%smc(4) = soilm4(t)
        NOAHMP401_struc(n)%noahmp401(t)%sh2o(4) = &
             NOAHMP401_struc(n)%noahmp401(t)%sh2o(4) + delta
        
     endif

     NOAHMP401_struc(n)%noahmp401(t)%wa = gws(t)

  enddo
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     swe_old = noahmp401_struc(n)%noahmp401(t)%sneqv

     dsneqv = swe(t)*1000.0 - noahmp401_struc(n)%noahmp401(t)%sneqv !in mm
     dsnowh = 0.0  ! GPS DA does not update snow depth, only SWE

     call noahmp401_snow_update(n, t, dsneqv, dsnowh)
     
     if(noahmp401_struc(n)%noahmp401(t)%sneqv.lt.0) then 
        write(LIS_logunit,*) '[WARN] GPS DA update resulted in negative snow values, rejecting update'
        write(LIS_logunit,*) '[WARN] dsneqv:', dsneqv
        write(LIS_logunit,*) '[WARN] target swe(t):', swe(t)
        write(LIS_logunit,*) '[WARN] current sneqv:', &
             noahmp401_struc(n)%noahmp401(t)%sneqv
        ! Reset to original values instead of stopping
        noahmp401_struc(n)%noahmp401(t)%sneqv = swe_old
     endif

  enddo

end subroutine noahmp401_setgpsdispvars

