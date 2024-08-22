SUBROUTINE reload_geoclass(fname,n,geoclassdata,numslc)

    USE WORLDVAR, ONLY : fileunit_temp
    USE MODVAR, ONLY : classtype, &
                       classattributedata,&
                       soilthick, & !OUT
                       soildepth, & !OUT
                       set_coded_classes, &
                       slc_iwet,slc_owet, &
                       maxsoillayers
    USE READWRITE_ROUTINES, ONLY : read_geoclass

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: fname    !<File path/name
    INTEGER, INTENT(IN) :: n                 !<Number of subbasin
    TYPE(CLASSTYPE),ALLOCATABLE, INTENT(INOUT) :: geoclassdata(:)  !<class characteristics
    INTEGER, INTENT(OUT) :: numslc           !<Number of new class in GeoClass.txt

    !Local parameters
    INTEGER,PARAMETER :: maxgccol = 3     !Number of data columns in GeoClass: combination	PREVIOUS_CLASS(for soil) NEW_CLASS(for luse) special_code

    !Local variables
    INTEGER j,ss,r                  !class index (slc)
    REAL, ALLOCATABLE :: tempsoilthick, tempsoildepth               !help variable
    REAL, ALLOCATABLE :: slcdata(:,:)   !data from Geoclass-file

    !>\b Algorithm \n
    !>Read the file with the change of class before and after land use change  (GeoClass.txt)
    CALL reread_geoclass(fileunit_temp,TRIM(fname),maxgccol,slcdata,numslc)

    IF(ALLOCATED(geoclassdata)) DEALLOCATE(geoclassdata)
    IF(.NOT.ALLOCATED(geoclassdata)) ALLOCATE(geoclassdata(numslc))
    r=size(soildepth,dim=2)
    IF(.NOT.ALLOCATED(tempsoildepth)) ALLOCATE(tempsoildepth(maxsoillayers,r))
    IF(.NOT.ALLOCATED(tempsoildepth)) ALLOCATE(tempsoildepth(maxsoillayers,r))
    tempsoildepth = soildepth
    tempsoilthick = soilthick
    IF(ALLOCATED(soildepth)) DEALLOCATE(soildepth)
    IF(.NOT.ALLOCATED(soildepth)) ALLOCATE(soildepth(maxsoillayers,nclass))
    IF(ALLOCATED(soilthick)) DEALLOCATE(soilthick)
    IF(.NOT.ALLOCATED(soilthick)) ALLOCATE(soilthick(maxsoillayers,nclass))
    DO j = 1, numslc
        geoclassdata(j)%luse  = classattributedata(NINT(slcdata(2,j)))%luse
        geoclassdata(j)%crop  = classattributedata(NINT(slcdata(2,j)))%crop
        geoclassdata(j)%crop2 = classattributedata(NINT(slcdata(2,j)))%crop2
        geoclassdata(j)%rotation = classattributedata(NINT(slcdata(2,j)))%rotation
        geoclassdata(j)%vegtype = classattributedata(NINT(slcdata(2,j)))%vegtype
        geoclassdata(j)%soil  = classattributedata(NINT(slcdata(1,j)))%soil 
        geoclassdata(j)%tiledepth = classattributedata(NINT(slcdata(1,j)))%tiledepth 
        geoclassdata(j)%streamdepth = classattributedata(NINT(slcdata(1,j)))%streamdepth
        DO ss = 1,maxsoillayers
            geoclassdata(j)%soildepth(ss) = classattributedata(NINT(slcdata(1,j)))%soildepth(ss)
        ENDDO
        soildepth(:,j) = tempsoildepth(:,NINT(slcdata(1,j)))
        soilthick(:,j) = tempsoilthick(:,NINT(slcdata(1,j)))
    ENDDO

    !>Find coded classes
    CALL set_coded_classes(numslc,slcdata(3,1:numslc))

    !Check no tiles in wetland
    IF(slc_iwet>0)THEN
      IF(geoclassdata(slc_iwet)%tiledepth>0.)THEN
        WRITE(6,*) 'ERROR: Tiles not allowed together with wetland. SLC: ',slc_iwet
        STOP 1
      ENDIF
    ENDIF
    IF(slc_owet>0)THEN
      IF(geoclassdata(slc_owet)%tiledepth>0.)THEN
        WRITE(6,*) 'ERROR: Tiles not allowed together with wetland. SLC: ',slc_owet
        STOP 1
      ENDIF
    ENDIF

    ! Set state of new class.
    CALL inherit_state_to_new_class(n,numslc,slcdata(1,1:numslc),frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    !Clean up subroutine
    IF(ALLOCATED(slcdata)) DEALLOCATE(slcdata)
    WRITE(6,*) 'Class information loaded (GeoClass)'

END SUBROUTINE reload_geoclass

! Inherit current state to new class. According to new class's previous class index
SUBROUTINE inherit_state_to_new_class(n,numslc,slcdata,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    USE MODVAR, ONLY : conduct, &
                       statesize,&
                       maxsoillayers

    !Argument declarations
    INTEGER, INTENT(IN) :: n                 !<Number of subbasin
    INTEGER, INTENT(IN) :: numslc            !<Number of class
    REAL,INTENT(IN) :: slcdata(:,:)          !data from Geoclass-file
    TYPE(snowicestatetype),TARGET,INTENT(INOUT) :: frozenstate  !<Snow and ice states
    TYPE(soilstatetype),TARGET,INTENT(INOUT)    :: soilstate    !<Soil states
    TYPE(aquiferstatetype),TARGET,INTENT(INOUT) :: aquiferstate !<Aquifer states
    TYPE(riverstatetype),TARGET,INTENT(INOUT)   :: riverstate   !<River states
    TYPE(lakestatetype),TARGET,INTENT(INOUT)    :: lakestate    !<Lake states
    TYPE(miscstatetype),TARGET,INTENT(INOUT)    :: miscstate    !<Misc states

    !Local variables
    TYPE(snowicestatetype), POINTER :: pfrozen  !Temporary snow and ice states
    TYPE(soilstatetype), POINTER    :: psoil     !Temporary soil states
    TYPE(aquiferstatetype), POINTER :: paquifer !Temporary Aquifer states
    TYPE(riverstatetype), POINTER   :: priver    !Temporary river states
    TYPE(lakestatetype), POINTER    :: plake     !Temporary lake states
    TYPE(miscstatetype) , POINTER   :: pmisc    !Temporary misc states
    TYPE(snowicestatetype) :: frozenstate2   !Temporary snow and ice states
    TYPE(soilstatetype)    :: soilstate2     !Temporary soil states
    TYPE(aquiferstatetype) :: aquiferstate2  !Temporary Aquifer states
    TYPE(riverstatetype)   :: riverstate2    !Temporary river states
    TYPE(lakestatetype)    :: lakestate2     !Temporary lake states
    TYPE(miscstatetype)    :: miscstate2     !Temporary misc states
    INTEGER :: i,j,ss

    !>temp store current state value
    CALL allocate_model_states(n,statesize,conduct, &
                               frozenstate2,soilstate2,aquiferstate2,riverstate2,lakestate2,miscstate2)
    ! use pointer
    pfrozen  => frozenstate
    psoil    => soilstate
    paquifer => aquiferstate
    priver   => riverstate
    pmisc    => miscstate
    ! Set temp state
    frozenstate2  = pfrozen
    soilstate2    = psoil
    aquiferstate2 = paquifer
    riverstate2   = priver
    miscstate2    = pmisc


    !RENEW STATESIZE & CONDUCT NEED CHANGE
    !renew statesize
    statesize%slcclass = numslc
    CALL deallocate_model_states(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
    CALL allocate_model_states(n,statesize,conduct, &
    frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    !> Inherit
    DO i = 1,n
        DO j = 1,numslc 
            !>frozenstate
            frozenstate%snow(j,i) = frozenstate2%snow(slcdata(j),i)
            frozenstate%snowage(j,i) = frozenstate2%snowage(slcdata(j),i)
            frozenstate%snowdepth(j,i) = frozenstate2%snowdepth(slcdata(j),i)
            frozenstate%snowcov(j,i) = frozenstate2%snowcov(slcdata(j),i)
            frozenstate%snowmax(j,i) = frozenstate2%snowmax(slcdata(j),i)
            IF(statesize%substance>0)THEN
              frozenstate%csnow(:,j,i) = frozenstate2%csnow(:,slcdata(j),)
            ENDIF
            IF(conduct%snowheat)THEN
              frozenstate%snowheat(j,i) = frozenstate2%snowheat(slcdata(j),i)
            ENDIF
            frozenstate%snowliq(j,i) = frozenstate2%snowliq(slcdata(j),i)

            !>soilstate
            soilstate%water(:,j,i) = soilstate2%water(:,slcdata(j),i)
            soilstate%temp(:,j,i) = soilstate2%temp(:,slcdata(j),i)
            soilstate%deeptemp(j,i) = soilstate2%deeptemp(slcdata(j),i)
            IF(statesize%substance>0)THEN
                soilstate%conc(:,:,j,i) = soilstate2%conc(:,slcdata(j),i)
                IF(ALLOCATED(soilstate%humusN))THEN
                    soilstate%humusN(:,j,i) = soilstate2%humusN(:,slcdata(j),i)
                    soilstate%fastN(:,j,i) = soilstate2%fastN(:,slcdata(j),i)
                ENDIF
                IF(ALLOCATED(soilstate%humusP))THEN
                    soilstate%partP(:,j,i) = soilstate2%partP(:,slcdata(j),i)
                    soilstate%fastP(:,j,i) = soilstate2%fastP(:,slcdata(j),i)
                    soilstate%humusP(:,j,i) = soilstate2%humusP(:,slcdata(j),i)
                    soilstate%PPrelpool(j,i) = soilstate2%PPrelpool(slcdata(j),i)
                ENDIF
                IF(ALLOCATED(soilstate%humusC))THEN
                    soilstate%fastC(:,j,i) = soilstate2%fastC(:,slcdata(j),i)
                    soilstate%humusC(:,j,i) = soilstate2%humusC(:,slcdata(j),i)
                    soilstate%oldgrw(j,i) = soilstate2%oldgrw(slcdata(j),i)
                ENDIF
                IF(ALLOCATED(soilstate%Srelpool)) soilstate%Srelpool(j,i) = soilstate2%Srelpool(slcdata(j),i)
                IF(ALLOCATED(soilstate%partT1))   soilstate%partT1(:,j,i) = soilstate2%partT1(:,slcdata(j),i)
                IF(ALLOCATED(soilstate%surface))  soilstate%surface(:,j,i) = soilstate2%surface(:,slcdata(j),i)
                IF(ALLOCATED(soilstate%icelens))  soilstate%icelens(j,i) = soilstate2%icelens(slcdata(j),i)
            ENDIF
            !>aquiferstate no
            !>riverstate no
            !>lakestate no
            !>miscstate
            IF(conduct%growthdegreeday)THEN
                miscstate%gdd(:,j,i) = miscstate2%gdd(:,slcdata(j),i)
                miscstate%gsbegin(:,j,i) = miscstate2%gsbegin(:,slcdata(j),i)
            ENDIF
            IF(conduct%irrigation)THEN
                miscstate%nextirrigation(j,i) = miscstate2%nextirrigation(slcdata(j),i)
                IF(ns>0)THEN
                  miscstate%cnextirrigation(:,j,i) = miscstate2%cnextirrigation(:,slcdata(j),i)
                ENDIF
            ENDIF
            IF(conduct%simT1)THEN
                miscstate%partT1sf(j,i) = miscstate2%partT1sf(slcdata(j),i)
            ENDIF
            
        ENDDO
    ENDDO

    CALL deallocate_model_states(frozenstate2,soilstate2,aquiferstate2,riverstate2,lakestate2,miscstate2)

END SUBROUTINE inherit_state_to_new_class

SUBROUTINE reread_geoclass(funit,infile,n,x,dmax)

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit           !<Unit for file
    CHARACTER (LEN=*), INTENT(IN) :: infile !<Name of file to be read
    INTEGER, INTENT(IN)  :: n               !<Maximum number of identifications
    REAL,ALLOCATABLE,INTENT(OUT) :: x(:,:)  !<Identification numbers
    INTEGER, INTENT(OUT) :: dmax            !<Number of combinations read
    
    !Local variables 
    INTEGER d,nslc,status,io
    REAL    y(n)                 !Data on row
    CHARACTER(LEN=80) line

    !>\b Algoritm \n
    !>Count number of classes (actually nslc include all rows)
    WRITE(6,*) 'File opened: ', TRIM(infile)
    CALL count_data_rows(funit,infile,0,nslc,status)
    IF(status/=0) STOP 1
    
    !>Allocate variable for holding class information
    IF(.NOT.ALLOCATED(x)) ALLOCATE(x(n,nslc))
    
    !>Read information of GeoClass.txt
    OPEN(UNIT = funit,FILE = infile, STATUS = 'old',ACTION='read',ERR=202) 
    !Skip comment lines in the beginning of the file
    DO 
      READ(funit,'(a)',END=200) line    
      IF(line(1:1)=='!')THEN 
      ELSE
        EXIT
      ENDIF
    ENDDO

    !Read class information (requires slc in order in file)
    d = 0
    DO
      y = 0
      READ(line,*,END=200,ERR=201) d,y(1:n)
      x(:,d)=y(:)
      DO
        READ(funit,'(a)',END=200,ERR=201,IOSTAT=io) line
          IF(line(1:1)=='!')THEN
            CYCLE
          ELSE
            EXIT
          ENDIF
      ENDDO
    ENDDO

200 CONTINUE
    CLOSE(funit)
    dmax = d
    IF(d==0)THEN
      WRITE(6,*) 'Error: No classes found in ',TRIM(infile)
      STOP 1
    ENDIF
    RETURN
    
201 CONTINUE
    WRITE(6,*) 'Error reading ',TRIM(infile)
    WRITE(6,*) 'Error io=',io
    CLOSE(funit)
    STOP 1
    RETURN
202 CONTINUE
    WRITE(6,*) 'Error opening ',TRIM(infile)
    STOP 1
    RETURN

END SUBROUTINE reread_geoclass

SUBROUTINE reallocate_nclass_based_variables()
    USE ATMOSPHERIC_PROCESSES, ONLY : calculate_class_wind_transformation_factor
    USE SOIL_PROCESSES, ONLY : initiate_soil_water
    USE MODVAR, ONLY : conductregest, &
                       genpar, &
                       basin, &
                       classbasin, &
                       nsub, &
                       nclass, &
                       numsubstances, &
                       timesteps_per_day, &
                       noutvarclass,   &
                       max_classoutvar 
    USE HYPE_INDATA, ONLY : set_regest_parameter
    USE HYPEVARIABLES, ONLY : m_pcelevth, &
                              m_pcelevadd,  &
                              m_pcelevmax,  &
                              m_pcelevstd,  &
                              n_pcet,n_pcea,n_pcem,  &
                              basinpreccorr, &  !OUT
    USE WORLDVAR, ONLY: outvarclassinfo, &
                        outvarclassinfotemp
    
    !Local variables
    REAL preccorr_hight
    REAL classheight  !masl
    REAL pcelevth,pcelevadd,pcelevmax   !loop-values
    INTEGER i,j

     IF(ALLOCATED(irrtype)) DEALLOCATE(irrtype)
     IF(.NOT.ALLOCATED(irrtype)) ALLOCATE(irrtype(nclass))

     IF(ALLOCATED(fieldneed)) THEN
        IF(ALLOCATED(irrigation)) DEALLOCATE(irrigation)
        IF(.NOT.ALLOCATED(irrigation)) ALLOCATE(irrigation(nclass,nsub))
        IF(ALLOCATED(cirrigation)) DEALLOCATE(cirrigation)
        IF(.NOT.ALLOCATED(cirrigation)) ALLOCATE(cirrigation(numsubstances,nclass,nsub))
        IF(ALLOCATED(fieldneed)) DEALLOCATE(fieldneed)
        IF(.NOT.ALLOCATED(fieldneed)) ALLOCATE(fieldneed(nclass,nsub))
     ENDIF

     IF(writeload)THEN
        IF(ALLOCATED(accload%slcclass)) DEALLOCATE(accload%slcclass)
        IF(.NOT.ALLOCATED(accload%slcclass)) ALLOCATE(accload%slcclass(nclass,max_classoutvar,numsubstances,nsub))
     ENDIF

     IF(ALLOCATED(windi)) CALL calculate_class_wind_transformation_factor(windtrans) 

     !>Set precipitation height correction
    IF(ALLOCATED(basinpreccorr)) DEALLOCATE(basinpreccorr)
    IF(.NOT.ALLOCATED(basinpreccorr)) ALLOCATE(basinpreccorr(nsub,nclass))
    pcelevth  = genpar(m_pcelevth)
    pcelevadd = genpar(m_pcelevadd)
    pcelevmax = genpar(m_pcelevmax)
    DO i = 1,nsub
      DO j= 1,nclass
        !Replace parameter values with regional parameter estimates
        IF(conductregest)THEN
          CALL set_regest_parameter(i,n_pcet,pcelevth)
          CALL set_regest_parameter(i,n_pcea,pcelevadd)
          CALL set_regest_parameter(i,n_pcem,pcelevmax)
        ENDIF
        classheight = basin(i)%elev+classbasin(i,j)%deltah  !masl

        !Calculate height correction
        IF(classheight>pcelevth)THEN
          preccorr_hight = (classheight-pcelevth) / 100. *  pcelevadd
          preccorr_hight = preccorr_hight + basin(i)%selev / 100. * genpar(m_pcelevstd)
          IF(preccorr_hight > pcelevmax) preccorr_hight = pcelevmax !maxmimum for prec hight correction
          basinpreccorr(i,j) = 1. + preccorr_hight
        ELSE
          basinpreccorr(i,j) = 1.
        ENDIF
      ENDDO
    ENDDO
    
    CALL initiate_soil_water()

    !Initiate soil temperature parameters and variables
    avertemp =(/5.,10.,20.,30./)*timesteps_per_day    !Number of timestep over which meantemp is calculated
    IF(ALLOCATED(soilmem)) DEALLOCATE(soilmem)
    IF(.NOT.ALLOCATED(soilmem)) ALLOCATE(soilmem(maxsoillayers,nclass))
    soilmem = 0.
    DO j= 1,nclass
      DO k = 1,maxsoillayers
        IF(k>1)THEN
          soilmem(k,j) = timesteps_per_day*landpar(m_surfmem,classdata(j)%luse)*EXP(landpar(m_depthrel,classdata(j)%luse)*(soildepth(k-1,j)+(soilthick(k,j) / 2.)))
        ELSE  
          soilmem(k,j) = timesteps_per_day*landpar(m_surfmem,classdata(j)%luse)*EXP(landpar(m_depthrel,classdata(j)%luse)*(soilthick(k,j) / 2.))
        ENDIF
      ENDDO
    ENDDO

    CALL calculate_landarea(nsub,landarea)

    !Initialize internal landclass wetland (iwet).
    IF(slc_iwet>0)THEN
        IF(ALLOCATED(iwetnoninflow)) DEALLOCATE(iwetnoninflow(nsub))
        IF(.NOT.ALLOCATED(iwetnoninflow)) ALLOCATE(iwetnoninflow(nsub))
        !iwetcatch is fraction of basinarea
        landfraction = 1.
        landfraction = landfraction - classbasin(:,slc_iwet)%part
        IF(slc_ilake>0)  landfraction = landfraction - classbasin(:,slc_ilake)%part
        IF(slc_lriver>0) landfraction = landfraction - classbasin(:,slc_lriver)%part
        IF(slc_olake>0)  landfraction = landfraction - classbasin(:,slc_olake)%part !floodplain not included in flow to local river
        IF(slc_mriver>0) landfraction = landfraction - classbasin(:,slc_mriver)%part !floodplain not included in flow to local river
        IF(slc_owet>0) landfraction = landfraction - classbasin(:,slc_owet)%part
        iwetnoninflow = 1.- basin%iwetcatch/landfraction  !fraction of runoff
        WHERE(iwetnoninflow<0.) iwetnoninflow = 0.    !Check in test.
      ENDIF

    !No initiation of floodwater to zero here, have inital value from state-file.
    !Calculate floodplain reference level (deepest bottom of floodplain), if not set in file
    IF(ALLOCATED(floodindex))THEN
        DO i = 1,nsub
          IF(floodindex(i)>0)THEN
            IF(flooding(floodindex(i))%hrefl==missing_value)THEN
              IF(slc_olake>0) flooding(floodindex(i))%hrefl = basin(i)%elev + classbasin(i,slc_olake)%deltah - basin(i)%lakedepth(2) + flooding(floodindex(i))%floll - flooding(floodindex(i))%flolp
            ENDIF
            IF(flooding(floodindex(i))%hrefr==missing_value)THEN
              IF(slc_mriver>0)THEN
                help = 0.   !river depth
                IF(deadriver(2,i)>0.) help = (basin(i)%area*classbasin(i,slc_mriver)%part*(1.-flooding(floodindex(i))%fpfmr))/deadriver(2,i)
                flooding(floodindex(i))%hrefr = basin(i)%elev + classbasin(i,slc_mriver)%deltah - help + flooding(floodindex(i))%flmrr - flooding(floodindex(i))%flmrp
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      IF(modeloption(p_deepgroundwater)==1)THEN
        CALL initiate_regional_groundwater_flow(nsub,numsubstances,slc_ilake,slc_olake,slc_lriver,slc_mriver,slc_iwet,slc_owet)
      ENDIF
      
      IF(conductload)THEN
        IF(ALLOCATED(Latmdep)) DEALLOCATE(Latmdep)
        IF(.NOT.ALLOCATED(Latmdep))  ALLOCATE(Latmdep(nclass,2,numsubstances))   !Substances in order IN, ON, SP, PP
        IF(ALLOCATED(Lcultiv)) DEALLOCATE(Lcultiv)
        IF(.NOT.ALLOCATED(Lcultiv))  ALLOCATE(Lcultiv(nclass,2,numsubstances))   !1=fertiliser, 2=plantdecay
        IF(ALLOCATED(Lirrsoil)) DEALLOCATE(Lirrsoil)
        IF(.NOT.ALLOCATED(Lirrsoil)) ALLOCATE(Lirrsoil(nclass,numsubstances))    !irrigation on soil
        IF(ALLOCATED(Lrurala)) DEALLOCATE(Lrurala)
        IF(.NOT.ALLOCATED(Lrurala))  ALLOCATE(Lrurala(nclass,numsubstances))     !rural a
        IF(ALLOCATED(Lstream)) DEALLOCATE(Lstream)
        IF(.NOT.ALLOCATED(Lstream))  ALLOCATE(Lstream(nclass,numsubstances))     !runoff from soil to stream
     ENDIF 
    IF(simulatesubstances)THEN
       IF(ALLOCATED(Lgrwsoil)) DEALLOCATE(Lgrwsoil)
       IF(.NOT.ALLOCATED(Lgrwsoil)) ALLOCATE(Lgrwsoil(nclass,numsubstances))       !regional groundwaterflow to soil
       IF(ALLOCATED(Lgrwclass)) DEALLOCATE(Lgrwclass)
       IF(.NOT.ALLOCATED(Lgrwclass)) ALLOCATE(Lgrwclass(nclass,numsubstances,nsub))       !regional groundwater outflow from soil
     ENDIF

     !Initate calculations for water balance output
    !IF(conductwb)THEN
    
    IF(ALLOCATED(outvarclassdata)) DEALLOCATE(outvarclassdata)
    IF(.NOT.ALLOCATED(outvarclassdata)) ALLOCATE(outvarclassdata(nclass,nsub,noutvarclass))
    IF(ALLOCATED(outvarclassfraction)) DEALLOCATE(outvarclassfraction)
    IF(.NOT.ALLOCATED(outvarclassfraction)) ALLOCATE(outvarclassfraction(nclass,nsub,noutvarclass))
    IF(conductload) THEN
      IF(ALLOCATED(outvar_classload)) DEALLOCATE(outvar_classload )
      IF(.NOT.ALLOCATED(outvar_classload)) ALLOCATE(outvar_classload(nclass,max_classoutvar,numsubstances,nsub))
    ENDIF
    IF(ALLOCATED(outvarclassinfo)) ALLOCATE(outvarclassinfo)
    IF(.NOT.ALLOCATED(outvarclassinfo)) ALLOCATE(outvarclassinfo(noutvarclass))
    DO i = 1,noutvarclass
      outvarclassinfo(i)=outvarclassinfotemp(I)
    ENDDO
    DEALLOCATE(outvarclassinfotemp)

END SUBROUTINE reallocate_nclass_based_variables