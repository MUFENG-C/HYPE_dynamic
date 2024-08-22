!> \file statedata.f90
!> Contains module state_datamodule.

!> \brief Load and save model states.
!>
!> Procedures for loading and saving initial states from file. Also 
!> processing them for submodel.
MODULE STATE_DATAMODULE
  !Copyright 2014-2021 SMHI
  !
  !This file is part of HYPE.
  !
  !HYPE is free software: you can redistribute it and/or modify it under
  !the terms of the Lesser GNU General Public License as published by
  !the Free Software Foundation, either version 3 of the License, or (at
  !your option) any later version. HYPE is distributed in the hope that
  !it will be useful, but WITHOUT ANY WARRANTY; without even the implied
  !warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
  !the Lesser GNU General Public License for more details. You should
  !have received a copy of the Lesser GNU General Public License along
  !with HYPE. If not, see <http://www.gnu.org/licenses/>.
  !
  !--------------------------------------------------------------------
  USE STATETYPE_MODULE
  USE LibDate
  !Subroutines also uses modules convert, modelmodule, modvar, readwrite_routines and worldvar.
  
  IMPLICIT NONE
  PRIVATE
! Private procedures
! -------------------
! divide_large_array  
! write_state_check
! read_and_perform_state_check
! check_basic_states
! check_substance_states
  PUBLIC :: initiate_state_for_submodel,&
            load_saved_state,&
            finalize_outstate, &
            read_state_check,&
            reset_soil_state,&
            save_soil_state_file

CONTAINS
  

  !>Initiate state variables for submodel simulation
  !----------------------------------------------------------------------
  SUBROUTINE initiate_state_for_submodel(dir,basemodel,   &
                frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 

    USE MODELMODULE, ONLY : initiate_model,   &
                            initiate_model_state
    USE MODVAR, ONLY : nsub,                     &
                       conduct, &
                       statesize
    USE WORLDVAR, ONLY : instate, &
                         resetstate, &
                         BASEMODELTYPE

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir      !<file directory
    TYPE(BASEMODELTYPE), INTENT(IN) :: basemodel  !<Information on basemodel
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate   !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states

    !Local variables
    TYPE(snowicestatetype) :: frozenstate2   !Temporary snow and ice states
    TYPE(soilstatetype)    :: soilstate2   !Temporary soil states
    TYPE(aquiferstatetype) :: aquiferstate2  !<Aquifer states
    TYPE(riverstatetype)   :: riverstate2  !Temporary river states
    TYPE(lakestatetype)    :: lakestate2   !Temporary lake states
    TYPE(miscstatetype)    :: miscstate2   !Temporary misc states

    !>\b Algoritm
    CALL deallocate_model_states(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
    !>If statefiles exist: read and store states temporary
    IF(instate%fromfile)THEN
      CALL allocate_model_states(basemodel%nsub,statesize,conduct, &
                               frozenstate2,soilstate2,aquiferstate2,riverstate2,lakestate2,miscstate2) 
      CALL load_saved_state(dir,basemodel%nsub,frozenstate2,soilstate2,aquiferstate2,riverstate2,lakestate2,miscstate2)
    ENDIF
    !>Reallocate state variables to submodel size
    CALL allocate_model_states(nsub,statesize,conduct, &
                               frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
    !>If statefiles exist: Initiate state variables from those and deallocate temporary storage
    IF(instate%fromfile)THEN
      CALL initiate_modelstate_submodel(nsub,statesize,conduct,basemodel%iindex,   &
                                        frozenstate,frozenstate2,soilstate,soilstate2, &
                                        aquiferstate,aquiferstate2,riverstate,riverstate2, &
                                        lakestate,lakestate2,miscstate,miscstate2)
      CALL deallocate_model_states(frozenstate2,soilstate2,aquiferstate2,riverstate2,lakestate2,miscstate2)
      IF(resetstate) CALL save_soil_state_file(dir,nsub,soilstate,miscstate)
    !>Else: Initiate state variables with default values
    ELSE
      CALL initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
    ENDIF
    !>Initiate other model variables and parameters
    CALL initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)  

  END SUBROUTINE initiate_state_for_submodel

  !>Load starting state from file and initiate state variables
  !------------------------------------------------------------------
  SUBROUTINE load_saved_state(dir,ns,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    USE MODVAR, ONLY : conduct, &
                       seconds_per_timestep, &
                       statesize, &
                       STATECONFIGURATIONTYPE
    USE WORLDVAR, ONLY : bdate, &
                         fileunit_temp, &   
                         instate, &
                         maxcharpath, &
                         seconds_per_day
    USE READWRITE_ROUTINES, ONLY : decompress_file_in_place, &
                                   read_array_from_file
 
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<file directory
    INTEGER, INTENT(IN)  :: ns            !<number of subbasins
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
    
    !Local variables
    INTEGER ffunit               
    INTEGER ios
    INTEGER dim,idim
    INTEGER ipiece,npiece
    INTEGER readsubst
    TYPE(STATECONFIGURATIONTYPE) :: statefileconfig
    INTEGER, ALLOCATABLE :: sectionlimits(:,:)
    REAL, ALLOCATABLE :: array(:)
    CHARACTER(LEN=28) filename   
    CHARACTER(LEN=16) bdatestr 
     
    !>\b Algorithm \n
    !>Set the file name 
    IF(seconds_per_timestep==seconds_per_day)THEN
      CALL format_date(bdate,'yyyymmdd',bdatestr)
    ELSE
      CALL format_date(bdate,'yyyymmddHHMM',bdatestr)
    ENDIF
    filename = 'state_save'//TRIM(ADJUSTL(bdatestr))//'.txt'                
    
    !>Decompress the state file if it is compressed
    IF(instate%compress) &
       CALL decompress_file_in_place(dir,filename)
    
    !>Open the file and check it
    ffunit = fileunit_temp
    OPEN(FILE=TRIM(dir)//TRIM(filename),UNIT=ffunit,STATUS='old',FORM='formatted',IOSTAT=ios,ACTION='read')
    IF(ios/=0) THEN
      WRITE(6,*) 'ERROR: Statefile ', TRIM(dir)//TRIM(filename), ' not found'
      STOP 1
    ENDIF
    CALL read_and_perform_state_check(ffunit,readsubst,statefileconfig)

    !>Get the frozen state variables from file
    CALL get_frozenstate_variables_arraysize(ns,readsubst,statesize,statefileconfig,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_frozenstate_variables_from_array(ns,readsubst, &
                  statesize,conduct,frozenstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Get the soil state variables from file
    CALL get_soilstate_variables_arraysize(ns,readsubst,statesize,statefileconfig,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_soilstate_variables_from_array(ns,readsubst,statesize,conduct,soilstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Get the aquifer state variables from file
    CALL get_aquiferstate_variables_arraysize(readsubst,statesize,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_aquiferstate_variables_from_array(readsubst,statesize,aquiferstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Get the river state variables from file
    CALL get_riverstate_variables_arraysize(ns,readsubst,statesize,  &
              statefileconfig,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_riverstate_variables_from_array(ns,readsubst,statesize,conduct, &
                  statefileconfig,riverstate,sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Get the lake state variables from file
    CALL get_lakestate_variables_arraysize(ns,readsubst,statesize,statefileconfig,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_lakestate_variables_from_array(ns,readsubst,statesize, &
                  statefileconfig,lakestate,sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Get the miscellaneuos state variables from file
    CALL get_miscstate_variables_arraysize(ns,statesize,statefileconfig,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_miscstate_variables_from_array(ns,readsubst,statesize, &
                  statefileconfig,conduct,miscstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Close the file, delete the uncomressed txt-file and end the subroutine
    IF(instate%compress)THEN
      CLOSE(ffunit,STATUS='DELETE')
    ELSE
      CLOSE(ffunit)
    ENDIF
    WRITE(6,*) 'File read: ', TRIM(dir)//TRIM(filename)

  END SUBROUTINE load_saved_state

  !>Saves state values for later use as starting state
  !---------------------------------------------------
  SUBROUTINE finalize_outstate(dir,ns,stateoutdate,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 

    USE MODVAR, ONLY : seconds_per_timestep,    &
                       conduct,  &
                       statesize
    USE WORLDVAR, ONLY : fileunit_temp, &
                         maxcharpath, &
                         outstate, &
                         seconds_per_day
    USE READWRITE_ROUTINES, ONLY : compress_and_delete_file, &
                                   write_array_to_file

    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir            !<file directory
    INTEGER, INTENT(IN)  :: ns            !<number of subbasins
    TYPE(DateType), INTENT(IN) :: stateoutdate     !<date for writing state           
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
    
    !Local variables
    INTEGER ipiece,npiece
    INTEGER ffunit
    INTEGER dim,idim
    INTEGER,ALLOCATABLE :: sectionlimits(:,:)
    REAL,ALLOCATABLE :: array(:)
    CHARACTER(LEN=16) stateoutdatestr
    CHARACTER(LEN=28) filename  
    
    !Local parameters
    !INTEGER, PARAMETER :: seconds_per_day  = 86400 

    !>\b Algorithm \n
    !>Set the file name 
    ffunit = fileunit_temp
    IF(seconds_per_timestep==seconds_per_day)THEN
      CALL format_date(stateoutdate,'yyyymmdd',stateoutdatestr)
    ELSE
      CALL format_date(stateoutdate,'yyyymmddHHMM',stateoutdatestr)
    ENDIF
    filename = 'state_save'//TRIM(ADJUSTL(stateoutdatestr))//'.txt'

    !>Open the file and write file checks to it
    OPEN(FILE=TRIM(dir)//TRIM(filename),UNIT=ffunit,STATUS='unknown',FORM='formatted',ACTION='write',ERR=900)
    CALL write_state_check(ffunit)  !uses nsub,conduct and statesize mm

    !>Write the frozen state variables to file
    CALL get_frozenstate_variables_arraysize(ns,statesize%substance,statesize,conduct,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_frozenstate_variables_to_array(ns,statesize,conduct,frozenstate,  &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Write the soil state variables to file
    CALL get_soilstate_variables_arraysize(ns,statesize%substance,statesize,conduct,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_soilstate_variables_to_array(ns,statesize,conduct,soilstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Write the aquifer state variables to file
    CALL get_aquiferstate_variables_arraysize(statesize%substance,statesize,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_aquiferstate_variables_to_array(statesize,aquiferstate,  &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Write the river state variables to file
    CALL get_riverstate_variables_arraysize(ns,statesize%substance,statesize,conduct,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_riverstate_variables_to_array(ns,statesize,conduct,riverstate,  &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Write the lake state variables to file
    CALL get_lakestate_variables_arraysize(ns,statesize%substance,statesize,conduct,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_lakestate_variables_to_array(ns,statesize,conduct,lakestate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Write the miscellaneuos state variables to file
    CALL get_miscstate_variables_arraysize(ns,statesize,conduct,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_miscstate_variables_to_array(ns,statesize,  &
                  conduct,miscstate,sectionlimits(1,ipiece), &
                  sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !>Close the file
    CLOSE(ffunit)
    
    !>If asked for: compress the state file
    IF(outstate%compress) &
      CALL compress_and_delete_file(dir,filename)
    RETURN
    
    !Error handling
900 CONTINUE 
    WRITE(6,*) 'ERROR: Open file '//TRIM(dir)//TRIM(filename)//' for writing.'
    STOP 1
    
  END SUBROUTINE finalize_outstate

  !>Read soil starting state from file and reinitiate some soil state variables
  !------------------------------------------------------------------
  SUBROUTINE reset_soil_state(dir,ns,soilstate,miscstate)

    USE MODVAR, ONLY : conduct, &
                       statesize
    USE WORLDVAR, ONLY : filename_resetstate, &
                         fileunit_temp
    USE READWRITE_ROUTINES, ONLY : read_array_from_file2
 
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<file directory
    INTEGER, INTENT(IN)  :: ns            !<number of subbasins (submodel)
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
    
    !Local variables
    INTEGER dim,idim
    INTEGER ipiece,npiece
    INTEGER i,ios,reclvar
    INTEGER, ALLOCATABLE :: sectionlimits(:,:)
    REAL, ALLOCATABLE :: array(:)
  INTEGER,PARAMETER :: maxchunksize = 100
     
    !Open file to read
    INQUIRE(IOLENGTH=reclvar) (0., i=1,maxchunksize)
    OPEN(FILE=TRIM(dir)//filename_resetstate,UNIT=fileunit_temp,STATUS='old',FORM='unformatted',IOSTAT=ios,ACTION='read',RECL=reclvar)
    IF(ios/=0)THEN
      WRITE(6,*) 'ERROR: Open file for reading: ',TRIM(dir)//filename_resetstate
      STOP 1
    ENDIF

    !Read soil states
    CALL get_reset_soilstate_variables_arraysize(ns,statesize,conduct,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file2(fileunit_temp,maxchunksize,idim,array)
        CALL set_reset_soilstate_variables_from_array(ns,statesize,conduct,soilstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !Read misc states 
    CALL get_reset_miscstate_variables_arraysize(ns,statesize,conduct,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file2(fileunit_temp,maxchunksize,idim,array)
        CALL set_reset_miscstate_variables_from_array(ns,statesize,conduct,miscstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !Finish the subroutine
    CLOSE(fileunit_temp)

  END SUBROUTINE reset_soil_state

  !>Saves selected state values for later use to reset states during simulation
  !---------------------------------------------------
  SUBROUTINE save_soil_state_file(dir,ns,soilstate,miscstate) 

    USE MODVAR, ONLY : conduct,  &
                       statesize
    USE WORLDVAR, ONLY : fileunit_temp, &
                         filename_resetstate
    USE READWRITE_ROUTINES, ONLY : write_array_to_file2

    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<file directory
    INTEGER, INTENT(IN)  :: ns            !<number of subbasins
    TYPE(soilstatetype),INTENT(IN) :: soilstate   !<Soil states
    TYPE(miscstatetype),INTENT(IN) :: miscstate   !<Misc states
    
    !Local variables
    INTEGER ipiece,npiece
    INTEGER dim,idim
    INTEGER i,ios,reclvar
    INTEGER,ALLOCATABLE :: sectionlimits(:,:)
    REAL,ALLOCATABLE :: array(:)
    INTEGER,PARAMETER :: maxchunksize = 100
    
    INQUIRE(IOLENGTH=reclvar) (0., i=1,maxchunksize)
    OPEN(FILE=TRIM(dir)//filename_resetstate,UNIT=fileunit_temp,STATUS='unknown',FORM='unformatted',IOSTAT=ios,ACTION='write',RECL=reclvar)
    IF(ios/=0)THEN
      WRITE(6,*) 'ERROR: Open file for writing: ',TRIM(dir)//filename_resetstate
      STOP 1
    ENDIF

    !Write soilstatetype soil states to be reset
    CALL get_reset_soilstate_variables_arraysize(ns,statesize,conduct,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_reset_soilstate_variables_to_array(ns,statesize,conduct,soilstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file2(fileunit_temp,maxchunksize,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    !Write miscstatetype soil states to be reset
    CALL get_reset_miscstate_variables_arraysize(ns,statesize,conduct,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_reset_miscstate_variables_to_array(ns,statesize,  &
                  conduct,miscstate,sectionlimits(1,ipiece), &
                  sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file2(fileunit_temp,maxchunksize,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CLOSE(fileunit_temp)

  END SUBROUTINE save_soil_state_file

  !--------------------------------------------------------------------
  !> Calculates appropriate size sections of large array
  !--------------------------------------------------------------------
  SUBROUTINE divide_large_array(dim,npiece,array)

  !Argument declarations
  INTEGER,INTENT(IN) :: dim
  INTEGER,INTENT(OUT) :: npiece
  INTEGER,ALLOCATABLE,INTENT(OUT) :: array(:,:)
  
  !Local varaibles
  INTEGER i
  INTEGER,PARAMETER :: maxchunksize = 12500000  !corresponds to real array of 50MB
  
  npiece = dim/maxchunksize
  IF(dim-maxchunksize*npiece>0) npiece = npiece + 1
  ALLOCATE(array(2,npiece))
  DO i = 1,npiece
    array(1,i) = 1 + (i-1)*maxchunksize
    array(2,i) = i*maxchunksize
  ENDDO
  array(2,npiece) = dim
  
  END SUBROUTINE divide_large_array
  
  !--------------------------------------------------------------------
  !> Saves values for later use as check if starting state is appropriate
  !--------------------------------------------------------------------
  SUBROUTINE write_state_check(ffunitloc) 

    USE MODVAR, ONLY : conduct, &
                       i_in,i_sp,i_t1,i_t2,i_oc,i_ss,i_ae,i_dsi,i_asi, &
                       modeloption,p_lakeriverice, &
                       nsub, &
                       statesize
    USE CONVERT, ONLY : logical_convert_to_integer
         
    !Argument declarations
    INTEGER, INTENT(IN) :: ffunitloc   !<File unit
    
    !Local variables
    INTEGER log2intvar1,log2intvar2,log2intvar3,log2intvar4,log2intvar5
    INTEGER log2intvar6,log2intvar8,log2intvar9,log2intvar10,log2intvar11
    INTEGER log2intvar12,log2intvar13,log2intvar14,log2intvar15,log2intvar16

    !Transform logical variables to integer
    log2intvar1 = logical_convert_to_integer(conduct%simN)
    log2intvar2 = logical_convert_to_integer(conduct%simP)
    log2intvar3 = logical_convert_to_integer(conduct%simC)
    log2intvar11 = logical_convert_to_integer(conduct%simS)
    log2intvar16 = logical_convert_to_integer(conduct%simSi)
    log2intvar4 = logical_convert_to_integer(conduct%riverwetland)
    log2intvar5 = logical_convert_to_integer(conduct%irrigation)
    log2intvar6 = logical_convert_to_integer(conduct%glacier)
    log2intvar8 = logical_convert_to_integer(conduct%arupdating)
    log2intvar9 = logical_convert_to_integer(conduct%growthdegreeday)
    log2intvar10 = logical_convert_to_integer(conduct%icelens)
    log2intvar12 = logical_convert_to_integer(conduct%watertransfer)
    log2intvar13 = logical_convert_to_integer(conduct%floodplain)
    log2intvar14 = logical_convert_to_integer(conduct%wetland)
    log2intvar15 = logical_convert_to_integer(conduct%snowheat)
    
    !Checkwrite for substances, number of subbasins and slc-classes
    WRITE(ffunitloc,'(10I3,I7,21I5)') statesize%substance, i_in,i_sp,i_t1, &
            i_t2,i_oc,i_ss,i_ae,i_dsi,i_asi,nsub,statesize%slcclass, &
            statesize%soillayer,statesize%riverqueue,statesize%timestep, &
            log2intvar1,log2intvar2,log2intvar3,log2intvar11,log2intvar16, &
            log2intvar4,log2intvar5,log2intvar6,modeloption(p_lakeriverice), &
            log2intvar8,log2intvar9,log2intvar10,log2intvar12,log2intvar13, &
            statesize%aquifer,log2intvar14,log2intvar15

  END SUBROUTINE write_state_check

  !> Check if starting state is appropriate
  !------------------------------------------------------------
  SUBROUTINE read_and_perform_state_check(ffunitloc,nsubst,config)

    USE MODVAR, ONLY : STATECONFIGURATIONTYPE, &
                       STATEDIMENSIONTYPE, &
                       conduct,statesize

    !Argument declarations
    INTEGER, INTENT(IN)  :: ffunitloc !<File unit
    INTEGER, INTENT(OUT) :: nsubst    !<number of substances in file
    TYPE(STATECONFIGURATIONTYPE), INTENT(OUT) :: config !<state file configuration
    
    !Local variables
    INTEGER :: nsub_file  !Number of subbasin of states in state file
    INTEGER :: lrice_file  !Modeloption lakeriverice in state file
    INTEGER :: sorder_file(9) !Order of substance indices in state file
    INTEGER :: status   !Status of deviation between models
    LOGICAL :: compatiblefile   !Status of statefile
    TYPE(STATEDIMENSIONTYPE) :: dim !state file state dimensions

    !>\b Algorithm \n
    !>Get model set-up for statefile and current model
    CALL read_state_check(ffunitloc,nsub_file,sorder_file,lrice_file,dim,config)
    
    !>Set output variables
    nsubst = dim%substance
   
    !>Compare model set-up with statefile set-up
    status = 0
    compatiblefile = check_basic_states(nsub_file,lrice_file,dim,config)
    IF(config%arupdating.NEQV.conduct%arupdating)THEN
      status = 1  !no vital difference
    ENDIF
    IF(statesize%substance==0)THEN
      IF(dim%substance>0) status = 1  !simulation of no substances from state with substances is ok
    ELSE
      compatiblefile = compatiblefile .AND. check_substance_states(sorder_file,dim,config)
    ENDIF
    IF(.NOT.compatiblefile) status = 2
    
    !>Write result of check to log-file
    IF(status==2)THEN
      WRITE(6,*) ' '
      WRITE(6,*) 'ERROR: State file not compatible with model set-up and simulation'
      STOP 1
    ELSEIF(status==1)THEN
      WRITE(6,*) ' '
      WRITE(6,*) 'INFO: State file separate from model set-up and simulation regarding '
      WRITE(6,*) 'substances or AR-updating. Simulation continue regardless of these differences.'
    ENDIF
    RETURN
   
  END SUBROUTINE read_and_perform_state_check
  
  !> Check if starting state is appropriate
  !------------------------------------------------------------
  SUBROUTINE read_state_check(ffunitloc,nsub,substanceorder,lrice,dim,config)

    USE MODVAR, ONLY : STATECONFIGURATIONTYPE, &
                       STATEDIMENSIONTYPE
    USE CONVERT, ONLY : integer_convert_to_logical

    !Argument declarations
    INTEGER, INTENT(IN)  :: ffunitloc !<File unit
    INTEGER, INTENT(OUT) :: nsub  !<Number of subbasins in saved state
    INTEGER, INTENT(OUT) :: substanceorder(9)  !<Order of substances in saved state file
    INTEGER :: lrice  !<Modeloption lakeriverice in saved state
    TYPE(STATEDIMENSIONTYPE), INTENT(OUT) :: dim !<state file state dimensions
    TYPE(STATECONFIGURATIONTYPE), INTENT(OUT) :: config !<state file configuration
    
    !Local variables
    INTEGER :: statecheck_file(32)  !Checkrow from saved state file

    !>\b Algorithm \n
    !>Read state file information
    READ(ffunitloc, *,ERR=200,END=201) statecheck_file
    
    !>Set output variables from state file information
    nsub = statecheck_file(11)   !Not in dim, because can be different in allocated states and state file
    substanceorder = statecheck_file(2:10) !Index of substances
    lrice = statecheck_file(24)  !Modeloption lakeriverice in state file
    
    dim%substance = statecheck_file(1)
    dim%slcclass = statecheck_file(12)
    dim%aquifer = statecheck_file(30)      !Number of aquifers in model domain
    dim%river = 2        !Number of rivers per subbasin (=2)
    dim%lake = 2         !Number of lakes per subbasin (=2)
    dim%soillayer = statecheck_file(13)
    dim%riverqueue = statecheck_file(14)
    dim%timestep = statecheck_file(15)
    
    config%simT1 = statecheck_file(4)>0
    config%simT2 = statecheck_file(5)>0
    config%simN = integer_convert_to_logical(statecheck_file(16))
    config%simP = integer_convert_to_logical(statecheck_file(17))
    config%simC = integer_convert_to_logical(statecheck_file(18))
    config%simS = integer_convert_to_logical(statecheck_file(19))
    config%simSi = integer_convert_to_logical(statecheck_file(20))
    config%arupdating = integer_convert_to_logical(statecheck_file(25))
    config%floodplain = integer_convert_to_logical(statecheck_file(29))
    config%glacier = integer_convert_to_logical(statecheck_file(23))
    config%growthdegreeday = integer_convert_to_logical(statecheck_file(26))
    config%icelens = integer_convert_to_logical(statecheck_file(27))
    config%irrigation = integer_convert_to_logical(statecheck_file(22))
    config%lakeriverice = (statecheck_file(24)==1).OR.(statecheck_file(24)==2)
    config%riverwetland = integer_convert_to_logical(statecheck_file(21))
    config%watertransfer = integer_convert_to_logical(statecheck_file(28))
    config%wetland = integer_convert_to_logical(statecheck_file(31))
    config%snowheat = integer_convert_to_logical(statecheck_file(32))
    config%qbank = config%simN .OR. config%simP .OR. config%simS .OR. config%simT1
    
    RETURN
    
200 WRITE(6,*) 'ERROR: reading first line of state-file'
    STOP 1
    RETURN
201 WRITE(6,*) 'ERROR: reading first line of state-file, end of file'
    STOP 1
    RETURN
    
  END SUBROUTINE read_state_check

  !> Check if the basic states in the state file are appropriate for this model
  !------------------------------------------------------------
  LOGICAL FUNCTION check_basic_states(n,lrice,dim,config)

    USE MODVAR, ONLY : STATECONFIGURATIONTYPE, &
                       STATEDIMENSIONTYPE, &
                       conduct, &
                       statesize, &
                       modeloption,p_lakeriverice
    USE WORLDVAR, ONLY : basemodel

    !Argument declarations
    INTEGER, INTENT(IN) :: n  !Number of subbasins in saved state
    INTEGER, INTENT(IN) :: lrice  !Modeloption lakeriverice in saved state
    TYPE(STATEDIMENSIONTYPE), INTENT(IN) :: dim !<state file state dimensions
    TYPE(STATECONFIGURATIONTYPE), INTENT(IN) :: config !<state file configuration
    
    !>\b Algorithm \n
    check_basic_states = .TRUE.
    
    !>Time step
    IF(dim%timestep/=statesize%timestep) check_basic_states = .FALSE.
    !>Basic dimensions
    IF(n/=basemodel%nsub) check_basic_states = .FALSE.
    IF(dim%slcclass/=statesize%slcclass) check_basic_states = .FALSE.
    IF(dim%soillayer/=statesize%soillayer) check_basic_states = .FALSE.
    IF(dim%riverqueue/=statesize%riverqueue) check_basic_states = .FALSE.
    !>Model options
    IF(config%riverwetland.NEQV.conduct%riverwetland) check_basic_states = .FALSE.
    IF(config%irrigation.NEQV.conduct%irrigation) check_basic_states = .FALSE.
    IF(config%glacier.NEQV.conduct%glacier) check_basic_states = .FALSE.
    IF(lrice/=modeloption(p_lakeriverice)) check_basic_states = .FALSE.
    IF(dim%aquifer/=statesize%aquifer) check_basic_states = .FALSE.
    IF(config%icelens.NEQV.conduct%icelens) check_basic_states = .FALSE.
    IF(config%watertransfer.NEQV.conduct%watertransfer) check_basic_states = .FALSE.
    IF(config%floodplain.NEQV.conduct%floodplain) check_basic_states = .FALSE.
    IF(config%wetland.NEQV.conduct%wetland) check_basic_states = .FALSE.
    IF(config%snowheat.NEQV.conduct%snowheat) check_basic_states = .FALSE.
   
  END FUNCTION check_basic_states

  !> Check if substance related states in the state file are appropriate for this model
  !------------------------------------------------------------
  LOGICAL FUNCTION check_substance_states(substanceorder,dim,config)

    USE MODVAR, ONLY : STATECONFIGURATIONTYPE, &
                       STATEDIMENSIONTYPE, &
                       conduct, &
                       statesize, &
                       i_in,i_sp,i_oc,i_ss,i_ae,i_dsi,i_asi,i_t1,i_t2

    !Argument declarations
    INTEGER, INTENT(IN) :: substanceorder(9)  !<state file substance order
    TYPE(STATEDIMENSIONTYPE), INTENT(IN) :: dim !<state file state dimensions
    TYPE(STATECONFIGURATIONTYPE), INTENT(IN) :: config !<state file configuration
    
    !>\b Algorithm \n
    check_substance_states = .TRUE.
    
    !>Basic dimensions
    IF(dim%substance/=statesize%substance) check_substance_states = .FALSE.
    IF(config%simN.NEQV.conduct%simN) check_substance_states = .FALSE.
    IF(config%simP.NEQV.conduct%simP) check_substance_states = .FALSE.
    IF(config%simC.NEQV.conduct%simC) check_substance_states = .FALSE.
    IF(config%simS.NEQV.conduct%simS) check_substance_states = .FALSE.
    IF(config%simSi.NEQV.conduct%simSi) check_substance_states = .FALSE.
    !>Order of substances
    IF(substanceorder(1)/=i_in) check_substance_states = .FALSE.
    IF(substanceorder(2)/=i_sp) check_substance_states = .FALSE.
    IF(substanceorder(3)/=i_t1) check_substance_states = .FALSE.
    IF(substanceorder(4)/=i_t2) check_substance_states = .FALSE.
    IF(substanceorder(5)/=i_oc) check_substance_states = .FALSE.
    IF(substanceorder(6)/=i_ss) check_substance_states = .FALSE.
    IF(substanceorder(7)/=i_ae) check_substance_states = .FALSE.
    IF(substanceorder(8)/=i_dsi) check_substance_states = .FALSE.
    IF(substanceorder(9)/=i_asi) check_substance_states = .FALSE.    
    !>Model options
    IF(config%growthdegreeday.NEQV.conduct%growthdegreeday) check_substance_states = .FALSE.
   
  END FUNCTION check_substance_states

END MODULE STATE_DATAMODULE
