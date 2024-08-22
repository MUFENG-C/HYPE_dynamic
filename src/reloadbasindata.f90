
SUBROUTINE reload_basindata(dir,n,status,year) 

    USE WORLDVAR, ONLY : fileunit_temp
    USE MODVAR, ONLY : classdata,         & !OUT 
                       nluse,nsoil,       & !OUT
                       nclass,            & !OUT
                       slc_iwet,slc_owet, &
                       conductdamconstruction, &   !OUT                
                       modeloption, &  !OUT
                       p_wetland, &
                       conduct, & !OUT
                       conductwarning, &
                       numrotations    !OUT

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir      !<File directory (model files)
    INTEGER, INTENT(IN)  :: year             !<current year for change geo 
    INTEGER, INTENT(IN) :: n                 !<Number of subbasins
    INTEGER, INTENT(OUT) :: status           !<Error status
    
    !Local parameters
    CHARACTER(LEN=12) :: curyear            !Year for file name
    CHARACTER(LEN=17) :: lfile              !Name of GeoClass file 
    CHARACTER(LEN=18) :: lfile2             !Name of column heading ClassData file to be read instead of GeoClass
    CHARACTER(LEN=16) :: infile             !Name of GeoData file 
    INTEGER,PARAMETER :: rcols = 50         !Number of data columns
    INTEGER,PARAMETER :: maxgccol = 13      !Maximum number of data columns in GeoClass
    INTEGER :: ncols                        !Total number of columns (rcols+number of s-lu-combinations (slc,dhslc,cr2) 
    
    !Local variables
    INTEGER,ALLOCATABLE :: lakedataid(:)    !lakedataid from GeoData

    !> \b Algorithm \n
    status = 0
    ! comfirm reload which geo file
    write(curyear,"(I4)") year
    lfile='GeoClass_'//trim(curyear)//'.txt'   !Name of GeoClass file 
    lfile2 = 'ClassData_'//trim(curyear)//'.txt'    !Name of column heading ClassData file to be read instead of GeoClass
    infile='GeoData_'//trim(curyear)//'.txt'   !Name of GeoData file 

    !>Read the class characteristics; either from ClassData or GeoClass
    CALL load_classdata(fileunit_temp,TRIM(dir)//lfile2,status)
    IF(status/=0)THEN
      status = 0
      CALL load_geoclass(TRIM(dir)//lfile,classdata,nclass,nluse,nsoil,numrotations)
    ENDIF
    
    IF(slc_iwet>0.OR.slc_owet>0)THEN
      conduct%wetland = .TRUE.
      IF(modeloption(p_wetland)/=2 .AND. conductwarning) WRITE(6,*) 'WARNING: Wetlands will be simulated although modeloption wetlandmodel 2 not set'
    ENDIF

    !  below read subbasin,only classbasin needs reload
    !>Count number of columns in GeoData
    CALL count_data_cols(fileunit_temp,TRIM(dir)//infile,0,ncols,status)
    IF(status/=0)RETURN

    !>Read the slc information of the basins (GeoData.txt) 
    CALL read_only_slc(fileunit_temp,TRIM(dir)//infile,n,ncols)

    !>Read additional description of important lakes (LakeData.txt) (possibly overwrite GeoData information)
    conductdamconstruction = .FALSE.   !Default, may be set in load_lakedata or load_damdata
    CALL load_lakedata(fileunit_temp,dir,n,lakedataid,status)
    IF(status/=0)RETURN
    IF(ALLOCATED(lakedataid)) DEALLOCATE(lakedataid)

    !>Read additional description of important dams (DamData.txt)
    CALL load_damdata(fileunit_temp,dir,n,status)
    IF(status/=0)RETURN

    CALL load_deposition(fileunit_temp,dir,n,status) 
    IF(status/=0)RETURN
    
END SUBROUTINE reload_basindata 
 
 
! pre-read GeoData.txt, get infomation except slc
subroutine pre_load_basindata(dir,fdir,n,status)
  
    USE WORLDVAR, ONLY : fileunit_temp
    USE MODVAR, ONLY : basin,             & !OUT

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir      !<File directory (model files)
    CHARACTER(LEN=*), INTENT(IN) :: fdir     !<File directory (forcing data)
    INTEGER, INTENT(OUT) :: n                !<Number of subbasins
    INTEGER, INTENT(OUT) :: status           !<Error status
    
    !Local parameters
    CHARACTER(LEN=11),PARAMETER :: infile='GeoData.txt'   !Name of GeoData file for the first simulate period
    INTEGER,PARAMETER :: rcols = 50         !Number of data columns
    INTEGER :: ncols                        !Total number of columns (rcols+number of s-lu-combinations (slc,dhslc,cr2) 

    CALL read_and_calc_basindata_withoutslc(fileunit_temp,TRIM(dir)//infile,n,lakedataid,ncols)


end subroutine pre_load_basindata

!without slc
SUBROUTINE read_and_calc_basindata_withoutslc(funit,infile,n,lakedataid,maxcol) 

    USE WORLDVAR, ONLY : i_str, &
                         i_intg, &
                         i_real, &
                         subweightcrit  !OUT
    USE MODVAR, ONLY : timesteps_per_day, &
                       conductbasinlocsoil, & !OUT
                       conductwarning, &
                       simulate, &  !OUT (%ilakecatch)
                       deposition, &
                       i_in,i_on,i_sp,i_pp,i_t1,i_t2,i_ss,i_ae,  &
                       i_open,i_forest,i_water, &
                       max_pstype, &
                       basin,        &    !OUT
                       pathsubid,    &    !OUT
                       load,         &    !OUT
                       wetland,      &    !OUT
                       petmodel,     &    !OUT
                       nregiondivisions, &
                       nregions, &        !OUT  !Number of parameter regions
                       lakesectiondata    !OUT
    USE MODEL_TEST_ROUTINES, ONLY : propagate_external_msg,e_geodata,e_error,e_info,e_warning

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                  !<Unit for file
    CHARACTER (LEN=*), INTENT(IN) :: infile        !<Name of characteristics file to be read
    INTEGER, INTENT(IN)  :: n                      !<Number of subbasins
    INTEGER, INTENT(OUT) :: lakedataid(n)          !<lakedataid from GeoData
    INTEGER, INTENT(IN)  :: maxcol                 !<Maximum number of data columns
!    INTEGER, INTENT(OUT) :: mcols                  !<Actual number of columns
    
    !Local parameters
    INTEGER, PARAMETER :: letters = 11  !Max number of letters (read) in heading, keep to 10 anyway

    !Local variables
    INTEGER i
    INTEGER idlks
    INTEGER xi(n,maxcol)               !Integer data read from file
    INTEGER code(maxcol)               !Code for column variable
    INTEGER rindex(maxcol)             !Index for column real variables
    INTEGER iindex(maxcol)             !Index for column integer variables
    INTEGER status                     !Error status
    INTEGER wsfdir                     !direction for Winstral coefficient
    INTEGER mcols                      !Actual number of columns
    LOGICAL deploadnfound,drydep_nfound,wetdep_nfound !Help variable for deallocating deposition-variables if not used
    LOGICAL locsoil_found              !Help variable for loc_soil column in GeoData
    LOGICAL lks_nfound,lks_dfound,lks_ifound,lks_afound !Help variables for deallocating lakesectiondata if not used
    REAL    onetspday                  !Reciprocal of timesteps per day
    REAL    xr(n,maxcol)               !Real data read from file
    REAL    localsource(n,8)           !Read rural source info (subbasins,column)
    CHARACTER(LEN=letters) str(maxcol)      !Content string
    CHARACTER(LEN=150) :: propagate_str                  !Propagate string
    CHARACTER(LEN=:),ALLOCATABLE :: propagate_str_extra  !Propagate string extra

    !Start of subroutine
    status = 0
    localsource = 0.
    deploadnfound = .FALSE.
    drydep_nfound = .FALSE.
    wetdep_nfound = .FALSE.
    locsoil_found = .FALSE.
    lks_nfound = .FALSE.
    lks_dfound = .FALSE.
    lks_ifound = .FALSE.
    lks_afound = .FALSE.
    onetspday = 1./REAL(timesteps_per_day)
    
    OPEN(UNIT = funit,FILE = infile, STATUS = 'old', ACTION='read')     !Open GeoData-file
    WRITE(6,*) 'File opened: ', TRIM(infile)

    !Reads the column headings from file
    CALL read_column_headings(funit,maxcol,letters,str,mcols,status)
    IF(status.NE.0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
      propagate_str = 'reading file: '//TRIM(infile)
      CALL propagate_external_msg(e_geodata,e_error,propagate_str)
      RETURN
    ENDIF

    IF(.NOT.ALLOCATED(propagate_str_extra)) ALLOCATE(CHARACTER(LEN=(13*maxcol)) :: propagate_str_extra)
    propagate_str_extra = 'columns ignored: -'

    !Code variables for easy finding of variable type
    code=i_str    !string, ignore
    DO i = 1,mcols
      IF(str(i)(1:letters)=='area       ') code(i) = i_real
      IF(str(i)(1:letters)=='subid      ') code(i) = i_intg
      IF(str(i)(1:letters)=='xcoord     ') code(i) = i_real
      IF(str(i)(1:letters)=='ycoord     ') code(i) = i_real
      IF(str(i)(1:letters)=='longitude  ') code(i) = i_real
      IF(str(i)(1:letters)=='latitude   ') code(i) = i_real
      IF(str(i)(1:letters)=='elev_mean  ') code(i) = i_real
      IF(str(i)(1:letters)=='elev_std   ') code(i) = i_real
      IF(str(i)(1:letters)=='slope_mean ') code(i) = i_real
      IF(str(i)(1:letters)=='slope_std  ') code(i) = i_real
      IF(str(i)(1:letters)=='ps1_vol    '.OR.str(i)(1:letters)=='ps2_vol    '.OR.str(i)(1:letters)=='ps3_vol    ')THEN
        IF(conductwarning) WRITE(6,*) 'WARNING: point sources no longer read from GeoData'
      ENDIF
      IF(str(i)(1:letters)=='loc_tp     ') code(i) = i_real
      IF(str(i)(1:letters)=='loc_sp     ') code(i) = i_real
      IF(str(i)(1:letters)=='loc_tn     ') code(i) = i_real
      IF(str(i)(1:letters)=='loc_in     ') code(i) = i_real
      IF(str(i)(1:letters)=='loc_t1     ') code(i) = i_real
      IF(str(i)(1:letters)=='loc_t2     ') code(i) = i_real
      IF(str(i)(1:letters)=='loc_ts     ') code(i) = i_real
      IF(str(i)(1:letters)=='loc_ss     ') code(i) = i_real
      IF(str(i)(1:letters)=='loc_vol    ') code(i) = i_real
      IF(str(i)(1:10)=='locsoil   ')THEN
        code(i) = i_real
        locsoil_found = .TRUE.
      ENDIF
      IF(str(i)(1:letters)=='region     ') code(i) = i_intg
      IF(str(i)(1:letters)=='petmodel   ') code(i) = i_intg
      IF(str(i)(1:letters)=='lakeregion ') code(i) = i_intg
      IF(str(i)(1:letters)=='parreg_4   ') code(i) = i_intg
      IF(str(i)(1:letters)=='parreg_5   ') code(i) = i_intg
      IF(str(i)(1:letters)=='ilregion   ') code(i) = i_intg
      IF(str(i)(1:letters)=='olregion   ') code(i) = i_intg
      IF(str(i)(1:letters)=='lake_depth ') code(i) = i_real
      IF(str(i)(1:letters)=='lake_whigh '.OR.str(i)(1:letters)=='lake_qavg  '.OR.   &
         str(i)(1:letters)=='lake_qamp  '.OR.str(i)(1:letters)=='lake_rate  '.OR.   &
         str(i)(1:letters)=='lake_qhigh '.OR.str(i)(1:letters)=='lake_exp   '.OR.   &
         str(i)(1:letters)=='mq         '.OR.str(i)(1:letters)=='regvol     '.OR.   &
         str(i)(1:letters)=='lake_qpha  '.OR.str(i)(1:letters)=='lake_wref  ')THEN
        IF(conductwarning) WRITE(6,*) 'WARNING: specific lake parameter no longer read from GeoData',str
      ENDIF
      IF(str(i)(1:letters)=='icatch     ') code(i) = i_real
      IF(str(i)(1:letters)=='iwetcatch  ') code(i) = i_real
      IF(str(i)(1:letters)=='mainfl     ') WRITE(6,*) 'WARNING: mainfl in GeoData is an column name no longer used, remove' 
      IF(str(i)(1:letters)=='branch     ') WRITE(6,*) 'WARNING: branch in GeoData is an column name no longer used, remove'
      IF(str(i)(1:letters)=='grwflow1   ') WRITE(6,*) 'WARNING: grwflow1 in GeoData is an column name no longer used, remove'
      IF(str(i)(1:letters)=='branchdown ') WRITE(6,*) 'WARNING: branchdown in GeoData is an column name no longer used, remove'
      IF(str(i)(1:letters)=='maindown   ') code(i) = i_intg
      IF(str(i)(1:letters)=='grwdown    ') code(i) = i_intg
      IF(str(i)(1:letters)=='grwolake   ') code(i) = i_real
      IF(str(i)(1:letters)=='parreg     ') code(i) = i_intg
      IF(str(i)(1:letters)=='wqparreg   ') code(i) = i_intg
      IF(str(i)(1:letters)=='rivlen     ') code(i) = i_real
      IF(str(i)(1:letters)=='loc_rivlen ') code(i) = i_real
      IF(str(i)(1:letters)=='wetdep_n   ')THEN
        code(i) = i_real
        wetdep_nfound = .TRUE.  ! For deallocating deposition%inwetconc if not used
      ENDIF
      IF(str(i)(1:letters)=='drydep_n1  ')THEN
        code(i) = i_real
        drydep_nfound = .TRUE.  ! For deallocating deposition%indryload if not used
      ENDIF
      IF(str(i)(1:letters)=='drydep_n2  ')THEN
        code(i) = i_real
        drydep_nfound = .TRUE.  ! For deallocating deposition%indryload if not used
      ENDIF
      IF(str(i)(1:letters)=='drydep_n3  ')THEN
        code(i) = i_real
        drydep_nfound = .TRUE.  ! For deallocating deposition%indryload if not used
      ENDIF
      IF(str(i)(1:letters)=='deploadn1  ') THEN
        code(i) = i_real
        deploadnfound = .TRUE.  ! For deallocating deposition% if total IN load on water is not used
      ENDIF
      IF(str(i)(1:letters)=='deploadn2  ') code(i) = i_real
      IF(str(i)(1:letters)=='deploadn3  ') code(i) = i_real
      IF(str(i)(1:letters)=='deploadn4  ') code(i) = i_real
      IF(str(i)(1:letters)=='deploadn5  ') code(i) = i_real
      IF(str(i)(1:letters)=='deploadn6  ') code(i) = i_real
      IF(str(i)(1:letters)=='deploadn7  ') code(i) = i_real
      IF(str(i)(1:letters)=='deploadn8  ') code(i) = i_real
      IF(str(i)(1:letters)=='deploadn9  ') code(i) = i_real
      IF(str(i)(1:letters)=='deploadn10 ') code(i) = i_real
      IF(str(i)(1:letters)=='deploadn11 ') code(i) = i_real
      IF(str(i)(1:letters)=='deploadn12 ') code(i) = i_real
      IF(str(i)(1:letters)=='lrwet_area ') code(i) = i_real
      IF(str(i)(1:letters)=='mrwet_area ') code(i) = i_real
      IF(str(i)(1:letters)=='lrwet_dep  ') code(i) = i_real
      IF(str(i)(1:letters)=='mrwet_dep  ') code(i) = i_real
      IF(str(i)(1:letters)=='lrwet_part ') code(i) = i_real
      IF(str(i)(1:letters)=='mrwet_part ') code(i) = i_real
      IF(str(i)(1:letters)=='buffer     ') code(i) = i_real
      IF(str(i)(1:letters)=='eroindex   ') code(i) = i_real
      IF(str(i)(1:letters)=='close_w    ') code(i) = i_real
      IF(str(i)(1:letters)=='weathcorr  ') code(i) = i_real
      IF(str(i)(1:letters)=='suspchannel') code(i) = i_real
      IF(str(i)(1:letters)=='lakedataid ') code(i) = i_intg   !id i LakeData
      IF(str(i)(1:letters)=='pobsid     ') WRITE(6,*) 'WARNING: pobsid etc in GeoData in no longer used, remove confusing columns'
      IF(str(i)(1:letters)=='tobsid     ') WRITE(6,*) 'WARNING: tobsid etc in GeoData in no longer used, remove confusing columns'
      IF(str(i)(1:letters)=='cloud_jan  ') code(i) = i_real
      IF(str(i)(1:letters)=='cloud_feb  ') code(i) = i_real
      IF(str(i)(1:letters)=='cloud_mar  ') code(i) = i_real
      IF(str(i)(1:letters)=='cloud_apr  ') code(i) = i_real
      IF(str(i)(1:letters)=='cloud_may  ') code(i) = i_real
      IF(str(i)(1:letters)=='cloud_jun  ') code(i) = i_real
      IF(str(i)(1:letters)=='cloud_jul  ') code(i) = i_real
      IF(str(i)(1:letters)=='cloud_aug  ') code(i) = i_real
      IF(str(i)(1:letters)=='cloud_sep  ') code(i) = i_real
      IF(str(i)(1:letters)=='cloud_oct  ') code(i) = i_real
      IF(str(i)(1:letters)=='cloud_nov  ') code(i) = i_real
      IF(str(i)(1:letters)=='cloud_dec  ') code(i) = i_real
      IF(str(i)(1:letters)=='clay       ') code(i) = i_real   !CB soil clay fraction
      IF(str(i)(1:letters)=='sand       ') code(i) = i_real   !CB soil sand fraction
      IF(str(i)(1:letters)=='silt       ') code(i) = i_real   !CB soil silt fraction
      IF(str(i)(1:letters)=='weight_sub ') code(i) = i_real
      IF(str(i)(1:letters)=='lks_num    ') THEN
        code(i) = i_intg  !number of lake type 1 sections in this subbasin
        lks_nfound = .TRUE.
      ENDIF
      IF(str(i)(1:7)=='lks_dp_')THEN
        code(i) = i_real  !lake section depth above lakedepth at outlet
        lks_dfound = .TRUE.
      ENDIF
      IF(str(i)(1:7)=='lks_fi_')THEN
        code(i) = i_real  !lake section fraction of icatch
        lks_ifound = .TRUE.
      ENDIF
      IF(str(i)(1:7)=='lks_fa_')THEN
        code(i) = i_real  !lake section fraction of lakearea
        lks_afound = .TRUE.
      ENDIF
      IF(code(i) == i_str) THEN
        propagate_str_extra = propagate_str_extra//','//TRIM(str(i)(1:letters))
      ENDIF
    ENDDO

    CALL propagate_external_msg(e_geodata,e_info,propagate_str_extra,code)
    IF(ALLOCATED(propagate_str_extra)) DEALLOCATE(propagate_str_extra)

    IF(.NOT.deploadnfound) DEALLOCATE(deposition%inloadwater)
    IF(.NOT.wetdep_nfound) DEALLOCATE(deposition%inwetconc)
    IF(.NOT.drydep_nfound)THEN
      DEALLOCATE(deposition%indryload)
    ELSE
      deposition%indryload = 0.
    ENDIF
    IF(.NOT.lks_nfound .OR. .NOT.lks_ifound .OR. .NOT.lks_afound .OR. .NOT.lks_dfound)THEN
      IF(ALLOCATED(lakesectiondata)) DEALLOCATE(lakesectiondata)
    ENDIF
    
    !Read all data
    CALL read_basindata5(funit,infile,maxcol,n,mcols,code,rindex,iindex,xi,xr)
    CLOSE(UNIT=funit)

    DO i = 1,mcols
      IF(str(i)(1:letters)=='area       ')   basin(1:n)%area      = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='subid      ')   basin(1:n)%subid     = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='xcoord     ')   basin(1:n)%xcoord    = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='ycoord     ')   basin(1:n)%ycoord    = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='longitude  ')   basin(1:n)%longitude = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='latitude   ')   basin(1:n)%latitude  = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='elev_mean  ')   basin(1:n)%elev      = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='elev_std   ')   basin(1:n)%selev     = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='slope_mean ')   basin(1:n)%slope     = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='slope_std  ')   basin(1:n)%sslope    = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='loc_tp     ')   localsource(1:n,1) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='loc_sp     ')   localsource(1:n,2) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='loc_tn     ')   localsource(1:n,3) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='loc_in     ')   localsource(1:n,4) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='loc_t1     ')   localsource(1:n,5) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='loc_t2     ')   localsource(1:n,6) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='loc_ts     ')   localsource(1:n,7) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='loc_ss     ')   localsource(1:n,8) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='loc_vol    ')  load(1:n)%volloc    = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='locsoil    ')  load(1:n)%locsoil   = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='region     ')  basin(1:n)%region   = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='petmodel   ')  petmodel(1:n) = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='lakeregion ')  basin(1:n)%parregion(6) = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='ilregion   ')  basin(1:n)%parregion(4) = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='olregion   ')  basin(1:n)%parregion(5) = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='parreg_4   ')  basin(1:n)%parregion(4) = xi(1:n,iindex(i)) !Alternative
      IF(str(i)(1:letters)=='parreg_5   ')  basin(1:n)%parregion(5) = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='lake_depth ')  basin(1:n)%lakedepth(2) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='icatch     ')  basin(1:n)%ilakecatch = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='iwetcatch  ')  basin(1:n)%iwetcatch  = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_jan  ')  basin(1:n)%cloudiness(1) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_feb  ')  basin(1:n)%cloudiness(2) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_mar  ')  basin(1:n)%cloudiness(3) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_apr  ')  basin(1:n)%cloudiness(4) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_may  ')  basin(1:n)%cloudiness(5) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_jun  ')  basin(1:n)%cloudiness(6) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_jul  ')  basin(1:n)%cloudiness(7) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_aug  ')  basin(1:n)%cloudiness(8) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_sep  ')  basin(1:n)%cloudiness(9) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_oct  ')  basin(1:n)%cloudiness(10) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_nov  ')  basin(1:n)%cloudiness(11) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='cloud_dec  ')  basin(1:n)%cloudiness(12) = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='maindown   ')  pathsubid(1:n)%main   = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='grwdown    ')  pathsubid(1:n)%grw1   = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='grwolake   ')  pathsubid(1:n)%grwtolake = xr(1:n,rindex(i))
      IF(str(i)(1:4)=='lks_' .AND. ALLOCATED(lakesectiondata))THEN
        IF(str(i)(1:letters)=='lks_num    ') basin(1:n)%lakesection = xi(1:n,iindex(i))
        IF(str(i)(1:7)=='lks_dp_')THEN
          idlks=0
          IF(ICHAR(str(i)(8:8))>=49 .AND. ICHAR(str(i)(8:8))<=57)THEN
            idlks = ICHAR(str(i)(8:8)) - 48
          ENDIF
          IF(ICHAR(str(i)(9:9))>=48 .AND. ICHAR(str(i)(9:9))<=57)THEN
            IF(idlks>0) idlks = idlks * 10
            idlks = idlks + ICHAR(str(i)(8:8)) - 48
          ENDIF
          lakesectiondata(1:n,idlks)%depth = xr(1:n,rindex(i))
        ENDIF
        IF(str(i)(1:7)=='lks_fi_')THEN
          idlks=0
          IF(ICHAR(str(i)(8:8))>=49 .AND. ICHAR(str(i)(8:8))<=57)THEN
            idlks = ICHAR(str(i)(8:8)) - 48
          ENDIF
          IF(ICHAR(str(i)(9:9))>=48 .AND. ICHAR(str(i)(9:9))<=57)THEN
            IF(idlks>0) idlks = idlks * 10
            idlks = idlks + ICHAR(str(i)(8:8)) - 48
          ENDIF
          lakesectiondata(1:n,idlks)%ficatch = xr(1:n,rindex(i))
        ENDIF
        IF(str(i)(1:7)=='lks_fa_')THEN
          idlks=0
          IF(ICHAR(str(i)(8:8))>=49 .AND. ICHAR(str(i)(8:8))<=57)THEN
            idlks = ICHAR(str(i)(8:8)) - 48
          ENDIF
          IF(ICHAR(str(i)(9:9))>=48 .AND. ICHAR(str(i)(9:9))<=57)THEN
            IF(idlks>0) idlks = idlks * 10
            idlks = idlks + ICHAR(str(i)(8:8)) - 48
          ENDIF
          lakesectiondata(1:n,idlks)%farea = xr(1:n,rindex(i))
        ENDIF
      ENDIF
      IF(str(i)(1:letters)=='parreg     ')  basin(1:n)%parregion(1) = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='wqparreg   ')  basin(1:n)%parregion(3) = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='loc_rivlen ')  basin(1:n)%rivlen(1)    = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='rivlen     ')  basin(1:n)%rivlen(2)    = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='wetdep_n   ')  deposition%inwetconc(1:n) = xr(1:n,rindex(i))*1.E-3    !ug/L -> mg/L
      IF(str(i)(1:letters)=='drydep_n1  ')  deposition%indryload(1:n,i_open) = xr(1:n,rindex(i))*onetspday  !kg/km2/d->kg/km2/ts
      IF(str(i)(1:letters)=='drydep_n2  ')  deposition%indryload(1:n,i_forest) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='drydep_n3  ')  deposition%indryload(1:n,i_water) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn1  ')  deposition%inloadwater(1:n,1) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn2  ')  deposition%inloadwater(1:n,2) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn3  ')  deposition%inloadwater(1:n,3) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn4  ')  deposition%inloadwater(1:n,4) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn5  ')  deposition%inloadwater(1:n,5) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn6  ')  deposition%inloadwater(1:n,6) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn7  ')  deposition%inloadwater(1:n,7) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn8  ')  deposition%inloadwater(1:n,8) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn9  ')  deposition%inloadwater(1:n,9) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn10 ')  deposition%inloadwater(1:n,10) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn11 ')  deposition%inloadwater(1:n,11) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='deploadn12 ')  deposition%inloadwater(1:n,12) = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:letters)=='lrwet_area ')  wetland(1:n,1)%area     = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='mrwet_area ')  wetland(1:n,2)%area     = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='lrwet_dep  ')  wetland(1:n,1)%depth    = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='mrwet_dep  ')  wetland(1:n,2)%depth    = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='lrwet_part ')  wetland(1:n,1)%part     = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='mrwet_part ')  wetland(1:n,2)%part     = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='buffer     ')  basin(1:n)%buffer       = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='eroindex   ')  basin(1:n)%eroindex     = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='close_w    ')  basin(1:n)%closewater   = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='weathcorr  ')  basin(1:n)%weatheringfactor = xr(1:n,rindex(i))
      IF(str(i)(1:letters)=='clay       ')  basin(1:n)%clay         = xr(1:n,rindex(i)) !subbasin soil clay fraction
      IF(str(i)(1:letters)=='sand       ')  basin(1:n)%sand         = xr(1:n,rindex(i)) !subbasin soil sand fraction
      IF(str(i)(1:letters)=='silt       ')  basin(1:n)%silt         = xr(1:n,rindex(i)) !subbasin soil silt fraction
      IF(str(i)(1:letters)=='suspchannel')  basin(1:n)%channelfactor = xr(1:n,rindex(i)) !river resuspension/erosion factor
      IF(str(i)(1:letters)=='lakedataid ')  lakedataid(1:n)         = xi(1:n,iindex(i))
      IF(str(i)(1:letters)=='weight_sub ')THEN
        IF(.NOT.ALLOCATED(subweightcrit)) ALLOCATE(subweightcrit(n))
        subweightcrit(1:n) = xr(1:n,rindex(i))
      ENDIF
    ENDDO

    !Check groundwater flow direction      
    IF(SUM(pathsubid(1:n)%grw1)==0) pathsubid(1:n)%grw1 = pathsubid(1:n)%main
    !Calculate number of parameter regions
    DO i=1,nregiondivisions
      nregions(i) = MAXVAL(basin(1:n)%parregion(i))
    ENDDO
    !Set which subbasins will use parameter for ilakecatch
    IF(.NOT.ALLOCATED(simulate%ilakecatpar)) ALLOCATE(simulate%ilakecatpar(n))
    simulate%ilakecatpar(1:n) = (basin(1:n)%ilakecatch<0.)
    !Set 
    DO i = 1,n
      !Calculate rural household loads (loc load in GeoData) or zeroing the loads.
      load(i)%locconc = 0.  !default/initialise
      IF(i_in>0)THEN
        load(i)%locconc(i_in) = localsource(i,3)*localsource(i,4)
        load(i)%locconc(i_on) = localsource(i,3)*(1. - localsource(i,4))
      ENDIF
      IF(i_sp>0)THEN
        load(i)%locconc(i_sp) = localsource(i,1)*localsource(i,2)
        load(i)%locconc(i_pp) = localsource(i,1)*(1. - localsource(i,2))
      ENDIF
      IF(i_t1>0) load(i)%locconc(i_t1) = localsource(i,5)
      IF(i_t2>0) load(i)%locconc(i_t2) = localsource(i,6)
      IF(i_ss>0)THEN
        load(i)%locconc(i_ss) = localsource(i,7)*localsource(i,8)
        load(i)%locconc(i_ae) = localsource(i,7)*(1. - localsource(i,8))
      ENDIF
    ENDDO
    conductbasinlocsoil = locsoil_found
    !Check subbasin criteria weight allocated and set
    IF(.NOT.ALLOCATED(subweightcrit))THEN
      ALLOCATE(subweightcrit(n))
      subweightcrit(1:n) = 1.   !default
    ENDIF
    WRITE(6,*) 'Geographical information(without SLC) loaded'

  END SUBROUTINE read_and_calc_basindata_withoutslc
  


  SUBROUTINE read_only_slc(funit,infile,n,maxcol) 

    USE WORLDVAR, ONLY : i_str, &
                         i_intg, &
                         i_real
    USE MODVAR, ONLY : conductwarning, &
                       classbasin,   &    !OUT
                       nclass
    USE MODEL_TEST_ROUTINES, ONLY : propagate_external_msg,e_geodata,e_error,e_info,e_warning

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                  !<Unit for file
    CHARACTER (LEN=*), INTENT(IN) :: infile        !<Name of characteristics file to be read
    INTEGER, INTENT(IN)  :: n                      !<Number of subbasins
    INTEGER, INTENT(IN)  :: maxcol                 !<Maximum number of data columns

    
    !Local parameters
    INTEGER, PARAMETER :: letters = 11  !Max number of letters (read) in heading, keep to 10 anyway

    !Local variables
    INTEGER i
    INTEGER idslc
    INTEGER xi(n,maxcol)               !Integer data read from file
    INTEGER code(maxcol)               !Code for column variable
    INTEGER rindex(maxcol)             !Index for column real variables
    INTEGER iindex(maxcol)             !Index for column integer variables
    INTEGER status                     !Error status
    INTEGER wsfdir                     !direction for Winstral coefficient
    INTEGER mcols                      !Actual number of columns
    REAL    xr(n,maxcol)               !Real data read from file
    CHARACTER(LEN=letters) str(maxcol)      !Content string
    CHARACTER(LEN=150) :: propagate_str                  !Propagate string
    CHARACTER(LEN=:),ALLOCATABLE :: propagate_str_extra  !Propagate string extra

    !Start of subroutine
    status = 0

    ! clean previous data
    IF(ALLOCATED(classbasin)) DEALLOCATE(classbasin)
    IF(.NOT.ALLOCATED(classbasin)) ALLOCATE(classbasin(n,nclass))

    OPEN(UNIT = funit,FILE = infile, STATUS = 'old', ACTION='read')     !Open GeoData-file
    WRITE(6,*) 'File opened: ', TRIM(infile)

    !Reads the column headings from file
    CALL read_column_headings(funit,maxcol,letters,str,mcols,status)
    IF(status.NE.0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
      propagate_str = 'reading file: '//TRIM(infile)
      CALL propagate_external_msg(e_geodata,e_error,propagate_str)
      RETURN
    ENDIF

    IF(.NOT.ALLOCATED(propagate_str_extra)) ALLOCATE(CHARACTER(LEN=(13*maxcol)) :: propagate_str_extra)
    propagate_str_extra = 'columns ignored: -'

    !Code variables for easy finding of variable type
    code=i_str    !string, ignore
    DO i = 1,mcols
      IF(str(i)(1:4)=='slc_')        code(i) = i_real
      IF(str(i)(1:6)=='dhslc_')      code(i) = i_real
      IF(str(i)(1:4)=='scr_')        code(i) = i_real 
      IF(str(i)(1:3)=='ws_')         code(i) = i_real
    ENDDO

    CALL propagate_external_msg(e_geodata,e_info,propagate_str_extra,code)
    IF(ALLOCATED(propagate_str_extra)) DEALLOCATE(propagate_str_extra)


    !Read all data
    CALL read_basindata5(funit,infile,maxcol,n,mcols,code,rindex,iindex,xi,xr)
    CLOSE(UNIT=funit)

    DO i = 1,mcols
      IF(str(i)(1:4)=='slc_') THEN
        idslc = 0
        IF(ICHAR(str(i)(5:5))>=49 .AND. ICHAR(str(i)(5:5))<=57)THEN
          idslc = ICHAR(str(i)(5:5)) - 48
        ENDIF
        IF(ICHAR(str(i)(6:6))>=48 .AND. ICHAR(str(i)(6:6))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(6:6)) - 48
        ENDIF
        IF(ICHAR(str(i)(7:7))>=48 .AND. ICHAR(str(i)(7:7))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(7:7)) - 48
        ENDIF
        IF(ICHAR(str(i)(8:8))>=48 .AND. ICHAR(str(i)(8:8))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(8:8)) - 48
        ENDIF
        IF(idslc<=nclass)THEN
          classbasin(1:n,idslc)%part=xr(1:n,rindex(i))
        ELSE
          IF(conductwarning) WRITE(6,*) 'WARNING: slc larger than nclass skipped: ', TRIM(str(i))
          propagate_str = 'slc larger than nclass skipped: '//TRIM(str(i))
          CALL propagate_external_msg(e_geodata,e_warning,propagate_str)
        ENDIF          
      ENDIF
      IF(str(i)(1:6)=='dhslc_') THEN
        idslc = 0
        IF(ICHAR(str(i)(7:7))>=49 .AND. ICHAR(str(i)(7:7))<=57)THEN
          idslc = ICHAR(str(i)(7:7)) - 48
        ENDIF
        IF(ICHAR(str(i)(8:8))>=48 .AND. ICHAR(str(i)(8:8))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(8:8)) - 48
        ENDIF
        IF(ICHAR(str(i)(9:9))>=48 .AND. ICHAR(str(i)(9:9))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(9:9)) - 48
        ENDIF
        IF(ICHAR(str(i)(10:10))>=48 .AND. ICHAR(str(i)(10:10))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(10:10)) - 48
        ENDIF
        IF(idslc<=nclass)THEN
          classbasin(1:n,idslc)%deltah = xr(1:n,rindex(i))
        ELSE
          IF(conductwarning) WRITE(6,*) 'WARNING: dhslc larger than nclass skipped: ', TRIM(str(i))
          propagate_str = 'dhslc larger than nclass skipped: '//TRIM(str(i))
          CALL propagate_external_msg(e_geodata,e_warning,propagate_str)
        ENDIF          
      ENDIF
      IF(str(i)(1:4)=='scr_') THEN    !area part with secondary crop
        idslc = 0
        IF(ICHAR(str(i)(5:5))>=49 .AND. ICHAR(str(i)(5:5))<=57)THEN
          idslc = ICHAR(str(i)(5:5)) - 48
        ENDIF
        IF(ICHAR(str(i)(6:6))>=48 .AND. ICHAR(str(i)(6:6))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(6:6)) - 48
        ENDIF
        IF(ICHAR(str(i)(7:7))>=48 .AND. ICHAR(str(i)(7:7))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(7:7)) - 48
        ENDIF
        IF(ICHAR(str(i)(8:8))>=48 .AND. ICHAR(str(i)(8:8))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(8:8)) - 48
        ENDIF
        IF(idslc<=nclass)THEN
          classbasin(1:n,idslc)%part2cr = xr(1:n,rindex(i))
        ELSE
          IF(conductwarning) WRITE(6,*) 'WARNING: scr larger than nclass skipped: ', TRIM(str(i))
          propagate_str = 'scr larger than nclass skipped: '//TRIM(str(i))
          CALL propagate_external_msg(e_geodata,e_warning,propagate_str)
        ENDIF          
      ENDIF
      IF(str(i)(1:3)=='ws_') THEN   !wsf
        idslc = 0
        IF(ICHAR(str(i)(4:4))>=49 .AND. ICHAR(str(i)(4:4))<=57)THEN 
          idslc = ICHAR(str(i)(4:4)) - 48
          wsfdir = ICHAR(str(i)(6:6)) - 48 
        ENDIF
        IF(ICHAR(str(i)(5:5))>=48 .AND. ICHAR(str(i)(5:5))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(5:5)) - 48
          wsfdir = ICHAR(str(i)(7:7)) - 48 
          IF(ICHAR(str(i)(6:6))>=48 .AND. ICHAR(str(i)(6:6))<=57)THEN
            IF(idslc>0) idslc = idslc * 10
            idslc = idslc + ICHAR(str(i)(6:6)) - 48
            wsfdir = ICHAR(str(i)(8:8)) - 48 
            IF(ICHAR(str(i)(7:7))>=48 .AND. ICHAR(str(i)(7:7))<=57)THEN
              IF(idslc>0) idslc = idslc * 10
              idslc = idslc + ICHAR(str(i)(7:7)) - 48
              wsfdir = ICHAR(str(i)(9:9)) - 48 
            ENDIF        
          ENDIF        
        ENDIF
        IF(idslc<=nclass)THEN
          classbasin(1:n,idslc)%wsf(wsfdir) = xr(1:n,rindex(i))
        ELSE
          IF(conductwarning) WRITE(6,*) 'WARNING: ws larger than nclass skipped: ', TRIM(str(i))
          propagate_str = 'ws larger than nclass skipped: '//TRIM(str(i))
          CALL propagate_external_msg(e_geodata,e_warning,propagate_str)
        ENDIF          
      ENDIF
    ENDDO
    WRITE(6,*) 'SLC information loaded'

  END SUBROUTINE read_only_slc