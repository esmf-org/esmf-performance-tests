!==============================================================================
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module NUOPC_COMP1

  !-----------------------------------------------------------------------------
  ! NUOPC_COMP1: Generic Component 1 (basic atmosphere)
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    modelSS     => SetServices

  implicit none

  private

  public SetServices

  real(ESMF_KIND_R8),parameter :: zeroValue    =      0.0_ESMF_KIND_R8
  real(ESMF_KIND_R8),parameter :: missingValue = 999999.0_ESMF_KIND_R8
  real(ESMF_KIND_R8),parameter :: islandLoc(4) = (/ 289.5_ESMF_KIND_R8, &
                                                    295.5_ESMF_KIND_R8, &
                                                      9.5_ESMF_KIND_R8, &
                                                     25.5_ESMF_KIND_R8 /)
  character(len=ESMF_MAXSTR) :: testType = "empty"

  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    logical                    :: attrPresent   
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! derive from NUOPC_Model
    call NUOPC_CompDerive(model, modelSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    ! specialize model
    call NUOPC_CompSpecialize(model, specLabel=label_Advertise, &
      specRoutine=Advertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    call NUOPC_CompSpecialize(model, specLabel=label_RealizeProvided, &
      specRoutine=Realize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
      specRoutine=Advance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    call NUOPC_CompSpecialize(model, specLabel=label_Finalize, &
      specRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

        
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Advertise(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)        :: importState, exportState
    logical                    :: attrPresent   
    
    rc = ESMF_SUCCESS

    ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    ! Read in test information
    call NUOPC_CompAttributeGet(model, name="testType", &
      isPresent=attrPresent, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    if (attrPresent) then
      call NUOPC_CompAttributeGet(model, name="testType", &
        value=testType, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return
    else
      testType = "empty"
    endif

        
    write(*,*) "Comp1 testType=",trim(testType)

    
    ! Make sure testtype is valid
    if ((testType /= "empty") .and. &
        (testType /= "geom") .and. &
        (testType /= "regrid")) then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
               msg="Invalid testType: "//trim(testType), &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
          return
    endif

    ! Advertise Fields, if we need to for the testType
    if (testType == "regrid") then
    
       ! importable field: sea_surface_temperature
       call NUOPC_Advertise(importState, &
            StandardName="sea_surface_temperature", name="sst", rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       
       ! exportable field: air_pressure_at_sea_level
       call NUOPC_Advertise(exportState, &
            StandardName="air_pressure_at_sea_level", name="pmsl", rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       
       ! exportable field: surface_net_downward_shortwave_flux
       call NUOPC_Advertise(exportState, &
            StandardName="surface_net_downward_shortwave_flux", name="rsns", rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Realize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR)         :: name
    integer                        :: localPet
    type(ESMF_State)               :: importState, exportState
    logical                        :: attrPresent
    type(ESMF_Field)               :: field
    type(ESMF_Grid)                :: defaultGrid
    type(ESMF_Mesh)                :: modelMesh
    type(ESMF_Geom)                :: modelGeom
    integer                        :: tlb(2), tub(2)
    real(ESMF_KIND_R8), pointer    :: lon_fptr(:)
    real(ESMF_KIND_R8), pointer    :: lat_fptr(:)
    integer(ESMF_KIND_I4), pointer :: msk_fptr(:,:)
    integer                        :: i,j
    character(ESMF_MAXSTR)         :: fieldCount
    character(len=ESMF_MAXSTR) :: geomValue = "default"
    character(len=ESMF_MAXSTR) :: geomFileValue = "(none)"

    ! Init to success
    rc = ESMF_SUCCESS

    ! query the Component for its name
    call ESMF_GridCompGet(model, name=name, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    ! query for importState and exportState
    call NUOPC_ModelGet(model, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    
    ! If the test type calls for it build a geom
    if ((testType == "geom") .or. &
        (testType == "regrid")) then

       write(*,*) "Comp1 build a geom"

       
       ! Get geom information
       call NUOPC_CompAttributeGet(model, name="Geom", &
            isPresent=attrPresent, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       if (attrPresent) then
          call NUOPC_CompAttributeGet(model, name="Geom", &
               value=geomValue, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return
       else
          geomValue = "default"
       endif
       
       
       ! Get geom file
       call NUOPC_CompAttributeGet(model, name="GeomFile", &
            isPresent=attrPresent, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       if (attrPresent) then
          call NUOPC_CompAttributeGet(model, name="GeomFile", &
               value=geomFileValue, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return
       else
          geomFileValue = "(none)"
       endif

       ! Output geom       
       if (localPet.eq.0) then
          write (*,"(A,A,A)") trim(name),": ", &
               "Geom = "//trim(geomValue)
          write (*,"(A,A,A)") trim(name),": ", &
               "GeomFile = "//trim(geomFileValue)
       endif
       
       ! Create Geom
       select case (geomValue)
       case ('default') ! By default, just build a small Grid       
          
          ! create a Grid object for Fields
          defaultGrid = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/24, 25/), &
               minCornerCoord=(/245._ESMF_KIND_R8, -5._ESMF_KIND_R8/), &
               maxCornerCoord=(/350._ESMF_KIND_R8, 55._ESMF_KIND_R8/), &
               coordSys=ESMF_COORDSYS_SPH_DEG, &
               staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
               rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return
          
          ! get grid coordinates
          call ESMF_GridGetCoord(defaultGrid, coordDim=1, &
               staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=lon_fptr, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return
          call ESMF_GridGetCoord(defaultGrid, coordDim=2, &
               staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=lat_fptr, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return
          
          ! add mask and island
          call ESMF_GridAddItem(defaultGrid, itemflag=ESMF_GRIDITEM_MASK, &
               staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return
          call ESMF_GridGetItem(defaultGrid, itemflag=ESMF_GRIDITEM_MASK, &
               staggerLoc=ESMF_STAGGERLOC_CENTER, &
               totalLBound=tlb, totalUBound=tub, farrayPtr=msk_fptr, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return
          msk_fptr = 0
          do j=tlb(2), tub(2)
             do i=tlb(1), tub(1)
                if ((lon_fptr(i).ge.islandLoc(1)) .AND. &
                     (lon_fptr(i).le.islandLoc(2)) .AND. &
                     (lat_fptr(j).ge.islandLoc(3)) .AND. &
                     (lat_fptr(j).le.islandLoc(4)) ) then
                   msk_fptr(i,j) = 1
                endif
             enddo
          enddo
          modelGeom = ESMF_GeomCreate(defaultGrid, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return
          
       case ('ESMF_FILEFORMAT_ESMFMESH')
          
          modelMesh = ESMF_MeshCreate(geomFileValue, &
               fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return  ! bail out
          modelGeom = ESMF_GeomCreate(modelMesh, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, &
               file=__FILE__)) &
               return
          
       case default
          
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
               msg="Invalid Geom value: "//trim(geomValue), &
               line=__LINE__, &
               file=__FILE__, &
               rcToReturn=rc)
          return
          
       endselect       
    endif

    ! If test type appropriate then add Fields
    if (testType == "regrid") then

       write(*,*) "Comp1 Adding Fields"

       
       ! importable field: sea_surface_temperature
       field = ESMF_FieldCreate(name="sst", geom=modelGeom, &
            typekind=ESMF_TYPEKIND_R8, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       call NUOPC_Realize(importState, field=field, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       ! fill import field sst with missingValue
       call ESMF_FieldFill(field, dataFillScheme="const", &
            const1=missingValue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       
       ! exportable field: air_pressure_at_sea_level
       field = ESMF_FieldCreate(name="pmsl", geom=modelGeom, &
            typekind=ESMF_TYPEKIND_R8, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       call NUOPC_Realize(exportState, field=field, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       
       ! fill export field pmsl with 95000
       call ESMF_FieldFill(field, dataFillScheme="const", &
            const1=95000.0_ESMF_KIND_R8, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       
       ! exportable field: surface_net_downward_shortwave_flux
       field = ESMF_FieldCreate(name="rsns", geom=modelGeom, &
            typekind=ESMF_TYPEKIND_R8, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       call NUOPC_Realize(exportState, field=field, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
       ! fill export field rsns with 200
       call ESMF_FieldFill(field, dataFillScheme="const", &
            const1=200.0_ESMF_KIND_R8, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
    endif
    
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Advance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    character(len=160)          :: msgString

    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

    ! Because of the way that the internal Clock was set by default,
    ! its timeStep is equal to the parent timeStep. As a consequence the
    ! currTime + timeStep is equal to the stopTime of the internal Clock
    ! for this call of the Advance() routine.

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing COMP1 from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    call ESMF_ClockPrint(clock, options="stopTime", &
      preString="---------------------> to: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Finalize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    ! local variables
    type(ESMF_VM)                           :: vm
    character(ESMF_MAXSTR)                  :: name
    character(len=8)                        :: value
    logical                                 :: zeroValues
    logical                                 :: missingValues
    integer                                 :: localPet
    type(ESMF_State)                        :: importState
    integer                                 :: itemCount, i
    integer                                 :: fieldCount
    character(len=80), allocatable          :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable  :: itemTypeList(:)
    type(ESMF_Field)                        :: field
    integer                                 :: rank
    real(ESMF_KIND_R8), pointer             :: farrayPtr1d(:)
    real(ESMF_KIND_R8), pointer             :: farrayPtr2d(:,:)
    integer                                 :: lclZero(1)
    integer                                 :: gblZero(1)
    integer                                 :: lclMissing(1)
    integer                                 :: gblMissing(1)

    rc = ESMF_SUCCESS

    ! query the Component for its vm
    call ESMF_GridCompGet(model, vm=vm, name=name, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    ! get diagnostic settings
    call ESMF_AttributeGet(model, name="zeroValues", &
      value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    value = ESMF_UtilStringLowerCase(value, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    select case (value)
      case ('true','t','yes')
        zeroValues=.true.
      case default
        zeroValues=.false.
    endselect

    call ESMF_AttributeGet(model, name="missingValues", &
      value=value, defaultValue="false", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    value = ESMF_UtilStringLowerCase(value, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    select case (value)
      case ('true','t','yes')
        missingValues=.true.
      case default
        missingValues=.false.
    endselect

#if 0    
    ! query the Component for its importState
    call NUOPC_ModelGet(model, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    ! fill import state with zeros
    call ESMF_StateGet(importState, nestedFlag=.true., &
      itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    allocate(itemNameList(itemCount), itemTypeList(itemCount))
    call ESMF_StateGet(importState, nestedFlag=.true., &
      itemNameList=itemNameList, itemTypeList=itemTypeList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    do i=1, itemCount
      if (itemTypeList(i)==ESMF_STATEITEM_FIELD) then
        call ESMF_StateGet(importState, field=field, itemName=itemNameList(i), &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return
        call ESMF_FieldGet(field, rank=rank, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return
        if (rank .eq. 1) then
          call ESMF_FieldGet(field, farrayPtr=farrayPtr1d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
          lclZero(1) = COUNT(farrayPtr1d(:).eq.zeroValue)
          lclMissing(1) = COUNT(farrayPtr1d(:).eq.missingValue)
        else
          call ESMF_FieldGet(field, farrayPtr=farrayPtr2d, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return
          lclZero(1) = COUNT(farrayPtr2d(:,:).eq.zeroValue)
          lclMissing(1) = COUNT(farrayPtr2d(:,:).eq.missingValue)
        endif
        call ESMF_VMReduce(vm, lclZero, gblZero, &
          reduceflag=ESMF_REDUCE_SUM, count=1, rootPet=0, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return
        call ESMF_VMReduce(vm, lclMissing, gblMissing, &
          reduceflag=ESMF_REDUCE_SUM, count=1, rootPet=0, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return
        if (localPet.eq.0) then
          select case (itemNameList(i))
            case ('sst')
              write (*,"(A,A,A,A,I0)") trim(name),": ", &
                trim(itemNameList(i)), &
                " zero values = ", gblZero(1)
              write (*,"(A,A,A,A,I0)") trim(name),": ", &
                trim(itemNameList(i)), &
                " missing values = ", gblMissing(1)
              if ((gblZero(1).gt.0).neqv.zeroValues) then
                write (*,"(A,A,A,L1,A)") "ERROR: ",trim(name), &
                  " zero values must be ", zeroValues, &
                  ", see [runconfig] file "
                call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
                  msg="Zero values does not match config settings", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
                return
              elseif ((gblMissing(1).gt.0).neqv.missingValues) then
                write (*,"(A,A,A,L1,A)") "ERROR: ",trim(name), &
                  " missing values must be ", missingValues, &
                  ", see [runconfig] file "
                call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
                  msg="Missing values does not match config settings", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
                return
              endif
            case default
              write (*,"(A)") "Field is unknown "//trim(itemNameList(i))
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
                msg="Field is unknown "//trim(itemNameList(i)), &
                line=__LINE__, &
                file=__FILE__, &
                rcToReturn=rc)
              return
          endselect
        endif
      endif
    enddo
    deallocate(itemNameList, itemTypeList)

#endif
    
    call ESMF_LogWrite("COMP1: Field check passed.", ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    
  end subroutine

  !-----------------------------------------------------------------------------

end module
