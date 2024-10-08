!> \file atm_proc.f90
!> Contains module atmospheric_processes.

!>Subroutines for calculating current atmospheric forcing
MODULE ATMOSPHERIC_PROCESSES

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

  !Uses modvar, hypevariables, hype_indata

  IMPLICIT NONE
  PRIVATE 
  !----------------------------------------------
  ! Private procedures 
  !----------------------------------------------
  !apply_classelevation_temperature_correction
  !get_snowfall_fraction
  !saturationpressure_function
  !netlongwaveradiation_function
  !calculate_vapour_pressures
  !calculate_shortwave_radiation
  !calculate_net_radiation
  !airpressure_elevationfunction
  !latentheat_tempfunction
  !psychrometric_constant
  !----------------------------------------------
  PUBLIC :: set_atmospheric_parameters_corrections, &
            calculate_class_atmospheric_forcing,  &
            calculate_subbasin_temperature, &
            calculate_subbasin_precipitation, &
            currently_snowing, &
            calculate_rain_snow_from_precipitation, & 
            calculate_extraterrestrial_radiation, &
            deltasaturationpressure_function, &
            set_class_precipitation_concentration_and_load, &
            calculate_class_wind_transformation_factor, &
            calculate_daylength, &
            calculate_winddirspeed, &
            calculate_snowfall_distribution, &
            get_current_cloudiness, &
            calculate_potential_evaporation
  CONTAINS

  !>Calculate basin temperature and precipitation corrections
  !>
  !>\b Consequences Module hypevariables variable basintcalt, basintempadd, 
  !>basinpreccorr, basinpcurain and basinpcusnow are allocated and set.
  !>
  !>\b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !----------------------------------------------------------
  SUBROUTINE set_atmospheric_parameters_corrections()
    
    USE MODVAR, ONLY : conductregest, &
                       genpar, &
                       regpar, &
                       regiondivision, &
                       tobselevation, &
                       basin, &
                       classbasin, &
                       nsub, &
                       nclass
    USE HYPE_INDATA, ONLY : set_regest_parameter
    USE HYPEVARIABLES, ONLY : m_tcalt,    &
                              m_pcelevth, &
                              m_pcelevadd,  &
                              m_pcelevmax,  &
                              m_pcelevstd,  &
                              m_tcobselev,   &
                              m_tcadd, m_tcelevadd, &
                              m_pcurain, &
                              m_pcusnow, &
                              n_tcalt,n_tcea,n_tcorr,  &
                              n_pcet,n_pcea,n_pcem,  &
                              n_pcur,n_pcus,  &
                              basintcalt,   &   !OUT
                              basintempadd, &   !OUT
                              basinpreccorr, &  !OUT
                              basinpcurain,  &  !OUT
                              basinpcusnow      !OUT

    !Local variables
    INTEGER i,j,isb
    REAL temppar1,temppar2
    REAL classheight  !masl
    REAL pcelevth,pcelevadd,pcelevmax   !loop-values
    REAL preccorr_hight
     
    !>\b Algoritm \n

    !>Set temperature altitude correction
    IF(.NOT.ALLOCATED(basintcalt)) ALLOCATE(basintcalt(nsub))
    basintcalt = genpar(m_tcalt)
    IF(conductregest)THEN
      DO i = 1,nsub
        CALL set_regest_parameter(i,n_tcalt,basintcalt(i))
      ENDDO
    ENDIF

    !>Set temperature correction
    IF(.NOT.ALLOCATED(basintempadd)) ALLOCATE(basintempadd(nsub))
    basintempadd = 0.
    IF(genpar(m_tcobselev)/=0.AND. ALLOCATED(tobselevation)) basintempadd = basintempadd - genpar(m_tcobselev)*(basin%elev-tobselevation)*0.01    !adjust for Tobs elevation
    temppar1 = genpar(m_tcelevadd)
    DO i = 1,nsub
      IF(basin(i)%parregion(regiondivision(m_tcadd))>0)THEN
        temppar2 = regpar(m_tcadd,basin(i)%parregion(regiondivision(m_tcadd)))
      ELSE
        temppar2  = 0.
      ENDIF
      
      !Replace parameter values with regional parameter estimates
      IF(conductregest)THEN
        CALL set_regest_parameter(i,n_tcea,temppar1)
        CALL set_regest_parameter(i,n_tcorr,temppar2)
      ENDIF
      
      !Adjust subbasin air temperature for observation elevation, subbasin elevation and regional correction.
      basintempadd(i) = basintempadd(i) + temppar2 - temppar1*basin(i)%elev*0.01       !Corrected subbasin temperature
    ENDDO

    !>Set precipitation height correction
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

    !>Set subbasin rain and snow parameters
    IF(.NOT.ALLOCATED(basinpcurain))THEN
      ALLOCATE(basinpcurain(nsub))
      ALLOCATE(basinpcusnow(nsub))
    ENDIF
    basinpcurain = genpar(m_pcurain)
    basinpcusnow = genpar(m_pcusnow)
    !>Replace parameter values with regional parameter estimates
    IF(conductregest)THEN
      DO isb=1,nsub
        CALL set_regest_parameter(isb,n_pcur,basinpcurain(isb))
        CALL set_regest_parameter(isb,n_pcus,basinpcusnow(isb))
      ENDDO
    ENDIF

  END SUBROUTINE set_atmospheric_parameters_corrections

  !>Calculate class temperature and precipitation
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !----------------------------------------------------------
  SUBROUTINE calculate_class_atmospheric_forcing(i,j,radext,cloud,  &
             temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,icpevap2,netrad,wind,swpot)
    
    USE MODVAR, ONLY : tempi,     &
                       preci,     &
                       tmini,     &
                       tmaxi,     &
                       windi,     &
                       humidi,    &
                       shortwavei,  &
                       genpar,    &
                       landpar,   &
                       basin,     &
                       classbasin,  &
                       classdata, &
                       missing_value
    USE HYPEVARIABLES, ONLY : m_pcluse,   &
                              m_alb,      &
                              m_mwind,    &
                              basinpreccorr, &
                              windtrans,  &
                              calcSWRAD,  &
                              calcVAPOUR, &
                              calcWIND
    
    !Argument declarations
    INTEGER,INTENT(IN) :: i       !<current index of subbasin
    INTEGER,INTENT(IN) :: j       !<current index of class
    REAL, INTENT(IN)  :: radext   !<subbasin extraterrestrial solar radiation (MJ/m2/day)
    REAL, INTENT(IN)  :: cloud    !<cloudiness for subbasin
    REAL, INTENT(OUT) :: temp     !<current class temperature (C)
    REAL, INTENT(OUT) :: prec     !<current class precipitation (mm/timestep)
    REAL, INTENT(OUT) :: tmin     !<current daily min temperature (C)
    REAL, INTENT(OUT) :: tmax     !<current daily max temperature (C)
    REAL, INTENT(OUT) :: swrad    !<daily mean shortwave radiation (MJ/m2/day)
    REAL, INTENT(OUT) :: rhmin    !<daily min relative humidity [fraction 0-1]
    REAL, INTENT(OUT) :: actvap   !<actual vapour pressure [kPa]
    REAL, INTENT(OUT) :: satvap   !<saturated vapour pressure [kPa]
    REAL, INTENT(OUT) :: icpevap  !<interception losses (evaporation) [mm]
    REAL, INTENT(OUT) :: icpevap2 !<interception losses (evaporation and condensation) (can be negative) [mm]
    REAL, INTENT(OUT) :: netrad   !<net downward radiation [MJ/m2/day]
    REAL, INTENT(OUT) :: wind     !<current wind speed (m/s)
    REAL, INTENT(OUT) :: swpot    !<clearsky radiation (m/s)
   
    !Local variables
    REAL rhmax,rhmean   !for future use: daily max/mean relative humidity
    REAL albedo         !current albedo (-)  
    REAL relswrad       !daily relative shortwave radiation (MJ/m2/day)
     
    !>\b Algoritm \n
    !>Set default output values for optional forcing; missing
    tmin = missing_value
    tmax = missing_value
    swrad = missing_value
    rhmin = missing_value
    actvap = missing_value
    satvap = missing_value
    icpevap = 0.
    netrad = missing_value
    wind = missing_value
    swpot = missing_value
      
    !Initiation of other variables
    rhmax = missing_value
    rhmean = missing_value
    relswrad = missing_value

    !>Calculate class temperature adjusted for elevation
    temp = apply_classelevation_temperature_correction(i,j,tempi(i))
    IF(ALLOCATED(tmini)) tmin = apply_classelevation_temperature_correction(i,j,tmini(i)) ! Adds SLC-hight dependent correction for minimum temperature
    IF(ALLOCATED(tmaxi)) tmax = apply_classelevation_temperature_correction(i,j,tmaxi(i)) ! Adds SLC-hight dependent correction for maximum temperature

    !>Calculate class precipitation adjusted for elevation and landuse bias, and calculate interception loss
    prec = preci(i) * basinpreccorr(i,j)
    IF(landpar(m_pcluse,classdata(j)%luse)>0.) icpevap = prec * landpar(m_pcluse,classdata(j)%luse)
    icpevap2 = prec * landpar(m_pcluse,classdata(j)%luse)
    prec = prec * (1. - landpar(m_pcluse,classdata(j)%luse))
      
    !>Calculate class specific shortwave radiation, if needed
    IF(calcSWRAD)THEN
      IF(ALLOCATED(shortwavei)) swrad = shortwavei(i)
      CALL calculate_shortwave_radiation(tmin,tmax,(basin(i)%elev + classbasin(i,j)%deltah),radext,cloud,swrad,relswrad)
      IF(relswrad.GT.0)THEN
        swpot = swrad/relswrad
      ELSE
        swpot = 0.
      ENDIF
    ENDIF

    !>Calculate class specific vapor pressures and net radiation, if needed (and possibly tmin and tmax)
    IF(calcVAPOUR)THEN
      IF(ALLOCATED(humidi)) rhmean = humidi(i)
      CALL calculate_vapour_pressures(temp,tmin,tmax,rhmean,rhmin,rhmax,swrad,radext,actvap,satvap)
      albedo = landpar(m_alb,classdata(j)%luse)
      IF(albedo.LE.0.) albedo = 0.23 ! FAO standard albedo
      CALL calculate_net_radiation(temp,tmin,tmax,albedo,actvap,swrad,relswrad,netrad)
    ENDIF
      
    !>Set class specific windspeed
    IF(calcWIND)THEN
      wind = genpar(m_mwind)
      IF(ALLOCATED(windi)) wind = windtrans(j) * windi(i)
    ENDIF

  END SUBROUTINE calculate_class_atmospheric_forcing

  !>Apply corrections to get subbasin average temperature
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !----------------------------------------------------------
  REAL FUNCTION apply_classelevation_temperature_correction(i,j,temp)
    
    USE HYPEVARIABLES, ONLY : basintcalt
    USE MODVAR, ONLY : classbasin
       
    !Argument declaration
    INTEGER, INTENT(IN) :: i    !<current subbasin
    INTEGER, INTENT(IN) :: j    !<current class
    REAL, INTENT(IN)    :: temp !<temperature to apply correction to
    
    !>\b Algoritm \n
    apply_classelevation_temperature_correction = temp
    
    !>Add SLC-hight dependent correction for temperature
    IF(basintcalt(i)/=0.)THEN
      apply_classelevation_temperature_correction = temp - basintcalt(i)*classbasin(i,j)%deltah*0.01
    ENDIF
    
  END FUNCTION apply_classelevation_temperature_correction

  !>Apply corrections to get subbasin average temperature
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !----------------------------------------------------------
  SUBROUTINE calculate_subbasin_temperature(n,month,temparr)
    
    USE HYPEVARIABLES, ONLY : basintempadd, &
                              m_mlapse
    USE MODVAR, ONLY : monthpar,  &
                       basin
       
    !Argument declaration
    INTEGER, INTENT(IN) :: n           !<number of subbasins
    INTEGER, INTENT(IN) :: month       !<current month
    REAL, INTENT(INOUT) :: temparr(n)  !<subbasin average temperature
    
    !>\b Algoritm \n
    
    !>Add temperature correction to current air temperature
    temparr = temparr + basintempadd - monthpar(m_mlapse,month)*basin%elev*0.01
    
  END SUBROUTINE calculate_subbasin_temperature

  !>Apply corrections to get subbasin average precipitation
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !-----------------------------------------------------------
  SUBROUTINE calculate_subbasin_precipitation(n,temppobs,precarr,pcorricep)
    
    USE MODVAR, ONLY : genpar,regpar, &
                       regiondivision, &
                       basin, &
                       classbasin, &
                       classdata, &
                       nclass
    USE HYPEVARIABLES, ONLY : m_preccorr,   &
                              m_pcaddg, & 
                              basinpcurain,basinpcusnow

    !Argument declaration
    INTEGER, INTENT(IN) :: n                 !<number of subbasins
    REAL, INTENT(IN)    :: temppobs(n)       !<temperature at precipitation input level
    REAL, INTENT(INOUT) :: precarr(n)        !<subbasin average precipitation [mm/timestep]
    REAL, INTENT(OUT)   :: pcorricep(n)      !<interception evaporation due to negative preccorr-parameter [mm/timestep]
    
    !Local variables
    INTEGER i,j
    REAL preccorr, temp
    REAL snowfall, snowfrac, preccorr_undercatch
 
    pcorricep = 0.    !Initiate default value, no interception evaporation
    !>\b Algoritm \n
    DO i = 1,n
      IF(basinpcurain(i)/=0. .OR. basinpcusnow(i)/=0.)THEN
        !>Calculate average snowfall fraction
      !DG 20141113 undercatch correction of precipitation, dependent on snowfall fraction, applied before general and regional corrections:
        !IF(ALLOCATED(snowfraci)) sffrac = snowfraci(i) !set snowfall fraction from input data, if existing
        ! loop over classes to derive average snowfall fraction dependent on class landuse and areal fractions in the subbasin.
        !   In the subbasin, the class snowfall fraction will also depend on the classtemperature, however, for the undercatch
        !   correction, we derive the snowfall fraction for each class, as if they were located at the elevation of 
        !   the precipitation input data. Ideally, we should also know what type of landuse the observation represent. But for the moment,
        !   we assume that the subbasin landuse distribution is relevant also for the input data. It could potentially give much more error
        !   to derive the snowfall fraction using the elevation corrected class temperature - so the current solution should be more conservative.
        !CP20200218 changed to use class temperature for snowfall fraction calculation
        snowfrac = 0.
        DO j=1,nclass
          IF(classbasin(i,j)%part.GT.0.)THEN
            temp = apply_classelevation_temperature_correction(i,j,temppobs(i))
            ! calculate snowfall fraction with existing subroutine, using unit precipitation
            snowfall = get_snowfall_fraction(i,classdata(j)%luse,1.,temp)
            ! summing areally weighted snowfall fraction
            snowfrac = snowfrac + snowfall * classbasin(i,j)%part
          ENDIF
        ENDDO
      
        !>Undercatch correction factor give first correction of precarr(i)
        preccorr_undercatch = basinpcurain(i) * (1. - snowfrac) + basinpcusnow(i) * snowfrac
        precarr(i) = precarr(i)*(1.+preccorr_undercatch)
      ENDIF
      
      !>Calculate adjustment of subbasin precipitation for regional correction and general bias.
      IF(basin(i)%parregion(regiondivision(m_preccorr))>0)THEN
        preccorr = 1. + regpar(m_preccorr,basin(i)%parregion(regiondivision(m_preccorr)))   !Regional correction of precipitation
      ELSE
        preccorr  = 1.
      ENDIF
      !>Calculate interception as negative preccorr
      IF(preccorr<1.) pcorricep(i) = precarr(i)*(1.+genpar(m_pcaddg)) * (1.-preccorr)
      !>Set new precipitation
      precarr(i) = precarr(i)*(1.+genpar(m_pcaddg)) * preccorr
    ENDDO 

  END SUBROUTINE calculate_subbasin_precipitation

  !>Subroutine for calculation snowfall fraction of precipitaion
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !--------------------------------------------------------------------------
  REAL FUNCTION get_rainfall_fraction(i,iluse,prec,temp)

    USE MODVAR, ONLY : genpar,      &
                       landpar,     &
                       modeloption, &
                       p_snowfall,  &
                       snowfraci
    USE HYPEVARIABLES, ONLY : m_ttpd,   &
                              m_ttpi,   &
                              m_ttmp

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<index of subbasin
    INTEGER, INTENT(IN) :: iluse    !<index of landuse
    REAL, INTENT(IN)    :: prec     !<precipitation (mm/timestep)
    REAL, INTENT(IN)    :: temp     !<air temperature (C)
   
    !Local variables
    REAL tt       !threshold temperature for rain/snow precip (C)
    REAL dtt      !temperature interval for rain/snow precip (C)
    REAL arain    !fraction of precipitation as snow (-)

    !> \b Algorithm \n
    !Set local parameter values
    dtt = genpar(m_ttpi)              !temperature interval rain/snow precip
    tt  = landpar(m_ttmp,iluse) + genpar(m_ttpd) !threshold temperature for rain/snow precip (ttmp for melt/evap + diff)
   
    !>Select model for rain/snow separation and calculate rain fraction
    SELECT CASE(modeloption(p_snowfall))
      CASE(1)       !Fraction of precipitation as snowfall is given in input data
        arain = 1. - snowfraci(i)
      CASE DEFAULT  !Original model based on threshold temperatures
        arain = 0.  !Pure snow
        IF(prec>0.)THEN
          IF(temp >= tt + dtt) THEN 
            arain = 1. !Pure rain
          ELSEIF(temp > tt - dtt .AND. temp < tt + dtt) THEN
            arain = (temp - (tt - dtt)) / (2 * dtt) !snow and rain mix
          ENDIF
        ENDIF
      END SELECT
      get_rainfall_fraction = arain
      
  END FUNCTION get_rainfall_fraction

  !>Subroutine for calculation snowfall fraction of precipitaion
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !--------------------------------------------------------------------------
  REAL FUNCTION get_snowfall_fraction(i,iluse,prec,temp)

    USE MODVAR, ONLY : genpar,      &
                       landpar,     &
                       modeloption, &
                       p_snowfall,  &
                       snowfraci
    USE HYPEVARIABLES, ONLY : m_ttpd,   &
                              m_ttpi,   &
                              m_ttmp

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<index of subbasin
    INTEGER, INTENT(IN) :: iluse    !<index of landuse
    REAL, INTENT(IN)    :: prec     !<precipitation (mm/timestep)
    REAL, INTENT(IN)    :: temp     !<air temperature (C)
   
    !Local variables
    REAL tt       !threshold temperature for rain/snow precip (C)
    REAL dtt      !temperature interval for rain/snow precip (C)
    REAL arain    !fraction of precipitation as rain (-)
    REAL asnow    !fraction of precipitation as snow (-)

    !> \b Algorithm \n
    !Set local parameter values
    dtt = genpar(m_ttpi)              !temperature interval rain/snow precip
    tt  = landpar(m_ttmp,iluse) + genpar(m_ttpd) !threshold temperature for rain/snow precip (ttmp for melt/evap + diff)
   
    !>Select model for rain/snow separation and calculate rain fraction
    SELECT CASE(modeloption(p_snowfall))
      CASE(1)       !Fraction of precipitation as snowfall is given in input data
        asnow = snowfraci(i)
      CASE DEFAULT  !Original model based on threshold temperatures
        arain = 0.  !Pure snow
        IF(prec>0.)THEN
          IF(temp >= tt + dtt) THEN 
            arain = 1. !Pure rain
          ELSEIF(temp > tt - dtt .AND. temp < tt + dtt) THEN
            arain = (temp - (tt - dtt)) / (2 * dtt) !snow and rain mix
          ENDIF
        ENDIF
        asnow = 1. - arain
        !asnow = 0.  !Pure rain �r inte identiskt
        !IF(prec>0.)THEN
        !  IF(temp < tt - dtt) THEN 
        !    asnow = 1. !Pure snow
        !  ELSEIF(temp > tt - dtt .AND. temp < tt + dtt) THEN
        !    asnow = (tt + dtt - temp) / (2 * dtt) !snow and rain mix
        !  ENDIF
        !ENDIF
      END SELECT
      get_snowfall_fraction = asnow
      
  END FUNCTION get_snowfall_fraction

  !>Subroutine for checking if it is currently snowing in the subbasin. A parameter 
  !>defines the limit for how large part of the subbasin must have snow for it to be sain to be snowing.
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !--------------------------------------------------------------------------
  LOGICAL FUNCTION currently_snowing(i,prec,tempin,lim)

    USE MODVAR, ONLY : classbasin, &
                       classdata, &
                       modeloption, &
                       p_snowfall,  &
                       nclass, &
                       snowfraci

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<index of subbasin
    REAL, INTENT(IN)    :: prec     !<precipitation (mm/timestep)
    REAL, INTENT(IN)    :: tempin   !<subbasin average temperature (C)
    REAL, INTENT(IN)    :: lim      !<subbasin area fraction defining it as snowing
   
    !Local variables
    INTEGER j
    REAL temp      !temperature of class (C)
    REAL snowfall    !fraction of precipitation as snow (-)

    !> \b Algorithm \n
    currently_snowing = .FALSE.   !Default not snowing
    
    !>If there is no precipitation there is no snowfall
    IF(prec==0.) RETURN
    
    !>Check if forcingdata says it is snowing
    IF(modeloption(p_snowfall)==1)THEN
      IF(snowfraci(i)>0.)THEN
        currently_snowing = .TRUE.
      ENDIF
      RETURN
    ENDIF
    
    !>Check how much of subbasin has snowfall
    snowfall = 0.
    DO j= 1,nclass
      IF(classbasin(i,j)%part>0.)THEN
        temp = apply_classelevation_temperature_correction(i,j,tempin)
        snowfall = snowfall + get_snowfall_fraction(i,classdata(j)%luse,1.,temp)*classbasin(i,j)%part
      ENDIF
    ENDDO
    IF(snowfall>=lim)THEN
      currently_snowing = .TRUE.
    ENDIF

  END FUNCTION currently_snowing

  !>Subroutine for calculation amount of rainfall and snowfall from precipitaion
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !--------------------------------------------------------------------------
  SUBROUTINE calculate_rain_snow_from_precipitation(i,iluse,prec,temp,snowfall,rainfall)

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<index of subbasin
    INTEGER, INTENT(IN) :: iluse    !<index of landuse
    REAL, INTENT(IN)    :: prec     !<precipitation (mm/timestep)
    REAL, INTENT(IN)    :: temp     !<air temperature (C)
    REAL, INTENT(OUT)   :: snowfall !<Precipitation as snow (mm/timestep)
    REAL, INTENT(OUT)   :: rainfall !<Precipitation as rain (mm/timestep)
   
    !Local variables
    REAL asnow    !fraction of precipitation as snow (-)

    !> \b Algorithm \n
    asnow = get_snowfall_fraction(i,iluse,prec,temp)
    
    !>Calculate rainfall and snowfall
    rainfall = prec * (1.-asnow)
    snowfall = prec * asnow

  END SUBROUTINE calculate_rain_snow_from_precipitation

  
  !>Get monthly cloudiness climatology for subbasins (Geodata input)
  !>----------------------------------------------------------------
  SUBROUTINE get_current_cloudiness(n,cmonth,cloud)

    USE MODVAR, ONLY : basin
    
    !Argument declarations
    INTEGER, INTENT(IN) :: n      !<number of subbasins
    INTEGER, INTENT(IN) :: cmonth  !<current month
    REAL, INTENT(OUT) :: cloud(n) !<cloudiness (fraction)
  
    cloud(1:n) = basin(1:n)%cloudiness(cmonth)
  
  END SUBROUTINE get_current_cloudiness

  !>Calculates extraterrestrial solar radiation as a function of
  !>latitude and day of year for all subbasins in the model domain
  !---------------------------------------------------------------
  SUBROUTINE calculate_extraterrestrial_radiation(n,jday,radext)

    USE MODVAR, ONLY : basin, &
                       pi,    & !=3.1415927
                       solar    !=0.0820 MJ/m2/min, solar constant
     
    !Argument declarations
    INTEGER,INTENT(IN) :: n         !<number of subbasins   
    INTEGER,INTENT(IN) :: jday      !<day of year (1-366)
    REAL,INTENT(OUT)   :: radext(n) !<extraterrestrial solar radiation [MJ/m2/day]
     
    !Local variables
    INTEGER i   !loop index
    REAL latrad !latitude in radians
    REAL dr     !inverse relative distance Earth-Sun
    REAL d      !Solar declination
    REAL omega  !Sunset hour angle
     
    !> \b Algorithm \n
    !>Calculate variables needed for extraterrestrial radiation calculation; distance to sun and declination.  
    dr = 1 + 0.033 * cos(2. * pi * jday / 365.)    !Inverse relative distance Earth-Sun (FAO, ekv 23)
    d = 0.409 * sin(2. * pi * jday / 365. - 1.39)  !Solar declination (FAO, ekv 24)

    !>For every subbasin:
    DO i = 1,n     
      latrad = basin(i)%latitude * pi / 180.0     !Transform basin latitude from degrees to radians 

      !>\li Calculate sunset hour angle, with special care for high latitudes
      omega = -tan(latrad) * tan(d) ! first check value of omega=cos(omega)
      IF(omega.GT.1)THEN
        omega = 0.              !Polar night, cos(omega)>1, set omega = 0
      ELSEIF(omega.LT.-1)THEN
        omega = pi              !Midnight sun, cos(omega)<1, set omega = pi
      ELSE
        omega = acos(-tan(latrad) * tan(d))    !Sunset hour angle,(FAO, ekv 25)
      ENDIF
       
      !>\li Calculate extraterrestrial radiation
      IF(omega.GT.0.)THEN
        radext(i) = 1440. / pi * solar * dr * (omega * sin(latrad)*sin(d)+cos(latrad)*cos(d)*sin(omega))
      ELSE
        radext(i) = 0.
      ENDIF
    ENDDO
    
  END SUBROUTINE calculate_extraterrestrial_radiation
   
  !>Calculates solar radiation at land surface.
  !>Calculates relative shortwave radiation and (if missing) daily mean shortwave radiation 
  !---------------------------------------------------------------
  SUBROUTINE calculate_shortwave_radiation(tmin,tmax,elev,radext,cloud,swrad,relswrad)

    USE MODVAR, ONLY: genpar,missing_value
    USE HYPEVARIABLES, ONLY: m_krs

    !Argument declarations
    REAL, INTENT(IN)  :: tmin           !<daily minimum air temperature [C]
    REAL, INTENT(IN)  :: tmax           !<daily maximum air temperature [C]
    REAL, INTENT(IN)  :: elev           !<elevation [m.a.s.l]
    REAL, INTENT(IN)  :: radext         !<extraterrestrial solar radiation [MJ/m2/day]
    REAL, INTENT(IN)  :: cloud          !<cloudiness [fraction]
    REAL, INTENT(INOUT) :: swrad        !<daily mean shortwave radiation [MJ/m2/day]
    REAL, INTENT(OUT) :: relswrad       !<relative shortwave radiation (actual/clearsky) [-]
  
    !Local parameters
    REAL, PARAMETER :: turbmin = 0.25   ! �ngstr�m formula, minimum turbidity
    REAL, PARAMETER :: turbmax = 0.75   ! �ngstr�m formula, maximum turbidity (at sea level)
    
    !Local variables
    REAL turbidity                      ! Turbidity (surface_radiation/top_of_atmosphere_radiation)
    REAL turbidity_clearsky             ! �ngstr�ms clear sky turbidity
    REAL turbidity_angstrom             ! �ngstr�ms cloudiness turbidity
    REAL turbidity_hargreaves           ! Hargreaves turbidity krs * (Tmax-Tmin)^0.5
 
    !Initialize turbidity as missing
    turbidity = missing_value
     
    !Clear sky turbidity, �ngstr�m formula, with elevation correction from FAO ekv 37:
    turbidity_clearsky = min(1.,turbmax + elev * 2.E-5)
   
    !Estimate turbidity from the input radiation data - if available - limited by �ngstr�m turbidity
    IF(swrad.NE.missing_value .AND. radext.GT.0.)THEN
      turbidity = max(turbmin,min(turbidity_clearsky, swrad / radext))
    ENDIF     
   
    !Estimate turbidity from cloudiness - if available - by �ngstr�m turbidity function
    IF(cloud.NE.missing_value)THEN
      turbidity_angstrom = ((turbidity_clearsky-turbmin)*(1.- cloud) + turbmin)
      IF(turbidity.EQ.missing_value)turbidity = turbidity_angstrom
    ENDIF
   
    !Estimate Hargreaves turbidity from tmin and tmax, if available
    IF(tmin.NE.missing_value .AND. tmax.GT.tmin)THEN
      !Hargreaves turbidity function, limited by the �ngstr�ms function
      turbidity_hargreaves = genpar(m_krs) * (tmax - tmin)**0.5
      turbidity_hargreaves = max(turbmin,min(turbidity_clearsky,turbidity_hargreaves))
      !Use Hargreaves turbidity if swrad was missing
      IF(turbidity.EQ.missing_value)turbidity = turbidity_hargreaves
    ENDIF
    
    !Assume clear sky turbidity if all of swrad, cloudiness and tmin/tmax were missing
    IF(turbidity.EQ.missing_value) turbidity = turbidity_clearsky

    !Calculate downward?? shortwave radiation, if missing
    IF(swrad.EQ.missing_value) swrad = radext * turbidity
    
    !Calculate relative shortwave radiation (shortwave/clearsky_shortwave)
    IF(turbidity_clearsky.GT.0 .AND. turbidity.GT.0)THEN
      relswrad = turbidity / turbidity_clearsky
    ELSE
      relswrad = 0.
    ENDIF
  END SUBROUTINE calculate_shortwave_radiation
  
  !>Calculates daily mean actual and saturated vapour pressure
  !!depending on data availability, following recommended FAO
  !!procedures.
  !---------------------------------------------------------------
  SUBROUTINE calculate_vapour_pressures(tmean,tmin,tmax,rhmean,rhmin,rhmax,swrad,radext,actvap,satvap)

    USE MODVAR, ONLY: genpar,missing_value
    USE HYPEVARIABLES, ONLY: m_krs

    !Argument declarations
    REAL, INTENT(IN)    :: tmean     !<daily mean temperature [C]? temperature of timestep 
    REAL, INTENT(INOUT) :: tmin      !<daily min temperature [C]
    REAL, INTENT(INOUT) :: tmax      !<daily max temperature [C]
    REAL, INTENT(IN)    :: rhmean    !<daily mean relative humidity [fraction 0-1]
    REAL, INTENT(IN)    :: rhmin     !<daily min relative humidity [fraction 0-1]
    REAL, INTENT(IN)    :: rhmax     !<daily max relative humidity [fraction 0-1]
    REAL, INTENT(IN)    :: swrad     !<daily mean shortwave radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: radext    !<extraterrestrial solar radiation [MJ/m2/day]
    REAL, INTENT(OUT)   :: actvap    !<actual vapour pressure [kPa]
    REAL, INTENT(OUT)   :: satvap    !<saturated vapour pressure [kPa]
     
    !Local variables
    REAL turbidity, trange
     
    !Initialize output variables as missing
    actvap = missing_value ; satvap = missing_value
     
    !Saturated vapor pressure from Tmin and Tmax (if avalable) or from Tmean
    IF(tmin.NE.missing_value .AND. tmax.NE.missing_value)THEN
      satvap = 0.5 * (saturationpressure_function(tmax)+saturationpressure_function(tmin))
    ELSE
      IF(tmean.NE.missing_value) satvap = saturationpressure_function(tmean)
    ENDIF
     
    !Actual vapour pressure, following FAO recommended procedure and function/data priority
    IF(tmin.NE.missing_value .AND. tmax.NE.missing_value .AND. rhmin.NE.missing_value .AND. rhmax.NE.missing_value)THEN
      !FAO, ekv 17, using rhmin, rhmax, tmin, and tmax
      actvap = 0.5 * (saturationpressure_function(tmax) * rhmin + & 
                      saturationpressure_function(tmin) * rhmax)
    ELSEIF(rhmax.NE.missing_value .AND. tmin.NE.missing_value)THEN
      !FAO, ekv 18, using rhmax and tmin
      actvap = saturationpressure_function(tmin) * rhmax
    ELSEIF(rhmean.NE.missing_value .AND. satvap.NE.missing_value)THEN
      !FAO ekv 19, using rhmean and saturation pressure from tmin/max or tmean
      actvap = rhmean * satvap 
    ELSEIF(tmin.NE.missing_value)THEN
      !Final solution, taking actual vapour pressure = saturation pressure at tmin
      actvap = saturationpressure_function(tmin)
    ELSEIF(swrad>=0. .AND. radext>0.)THEN
      !As a finalfinal solution, Tmin and Tmax is inferred from 
      ! the Hargreaves turbidity function (if swrad and radext is available), 
      ! and then actvap is estimated from the estimated Tmin
     
      !Turbidity
      turbidity = swrad / radext
      IF(genpar(m_krs)>0.)THEN
        !Diurnal temperature range inferred from Hargreaves turbidity function
        trange = (turbidity/genpar(m_krs))**2
        !Assume min and max is evenly distributed around the mean
        tmin = tmean - 0.5 * trange
        tmax = tmean + 0.5 * trange
      ELSE
        tmin= tmean   !CP added 181108 to handle divide by zero, better formula to use? Add to test!
      ENDIF
      !Actual vapour pressure from tmin
      actvap = saturationpressure_function(tmin)
    ENDIF
     
    !Finally, make sure actvap <= satvap
    IF(actvap.NE.missing_value .AND. satvap.NE.missing_value) actvap = min(actvap,satvap)
       
  END SUBROUTINE calculate_vapour_pressures
  
  !>Calculates net radiation at land surface 
  !---------------------------------------------------------------
  SUBROUTINE calculate_net_radiation(tmean,tmin,tmax,albedo,actvap,swrad,relswrad,netrad)
    USE MODVAR, ONLY: missing_value
    
    !Argument declarations
    REAL, INTENT(IN)    :: tmean      !<daily mean temperature [C]? temperature of timestep 
    REAL, INTENT(IN)    :: tmin       !<daily min temperature [C]
    REAL, INTENT(IN)    :: tmax       !<daily max temperature [C]
    REAL, INTENT(IN)    :: albedo     !<albedo [fraction, 0-1]
    REAL, INTENT(IN)    :: actvap     !<actual vapour pressure [kPa]
    REAL, INTENT(IN)    :: swrad      !<daily mean shortwave radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: relswrad   !<relative shortwave radiation (actual/clearsky) [-]
    REAL, INTENT(INOUT) :: netrad     !<net downward radiation [MJ/m2/day]
      
    !Local variables
    REAL radnetshort, radnetlong
   
    !Estimate net radiation if missing, following FAO recommended procedure
    IF(netrad.EQ.missing_value)THEN
      !Net shortwave radiation
      radnetshort = swrad * (1.-albedo)

      !Net downward longwave radiation using tmin/max, actual vapour pressure, and relative shortwave, if available
      IF(tmin.NE.missing_value .AND. tmax.NE.missing_value .AND. relswrad.NE.missing_value .AND. actvap.NE.missing_value)THEN
        radnetlong = netlongwaveradiation_function(tmax,tmin,actvap,relswrad)
      ELSEIF(tmean.NE.missing_value .AND. relswrad.NE.missing_value .AND. actvap.NE.missing_value)THEN
        !Use Tmean if Tmin and Tmax is missing
        radnetlong = netlongwaveradiation_function(tmean,tmean,actvap,relswrad)
      ELSE
        radnetlong = 0.
      ENDIF

      !net radiation
      netrad = radnetshort - radnetlong
    ENDIF
    
  END SUBROUTINE calculate_net_radiation
 
  !>Saturation pressure (kPa) as a function of temperature in deg C, following FAO
  !-------------------------------------------------------------------------------
  REAL FUNCTION saturationpressure_function(temp)
    
    !Argument declarations
    REAL,INTENT(IN) :: temp
    
    saturationpressure_function = 0.6108 * EXP((17.27 * temp)/(temp+237.3))
    
  END FUNCTION saturationpressure_function

  !>Slope of the saturation pressure temperature function (kPa/C)
  !-------------------------------------------------------------------------------
  REAL FUNCTION deltasaturationpressure_function(temp)
    
    !Argument declarations
    REAL,INTENT(IN) :: temp 
    
    deltasaturationpressure_function = 4098 * 0.6108 * EXP((17.27 * temp)/(temp+237.3)) / (temp + 237.3)**2
    
  END FUNCTION deltasaturationpressure_function
  
  !>Net (upward) longwave radiation, FAO, ekv 39
  !-------------------------------------------------------------------------------
  REAL FUNCTION netlongwaveradiation_function(tmax,tmin,actvap,relshortwave)
    
    !Argument declarations
    REAL,INTENT(IN) :: tmax,tmin,actvap,relshortwave
     
    !Local parameters
    REAL, PARAMETER :: sigma = 4.903E-9 !Stefan-Boltzmann, MJ K^-4 m^-2 day^-1
    REAL, PARAMETER :: vpar1 = 0.34     !vapour pressure constants
    REAL, PARAMETER :: vpar2 = 0.14
    REAL, PARAMETER :: rpar1 = 1.35     !cloudiness constants
    REAL, PARAMETER :: rpar2 = 0.35
    
    !Net upward longwave radiation
    netlongwaveradiation_function = sigma * 0.5 * ((tmax+273.15)**4 + (tmin+273.15)**4) * &
       (vpar1-vpar2*(actvap)**0.5) * (rpar1 * MIN(1.,relshortwave)- rpar2)
    
  END FUNCTION netlongwaveradiation_function

  !>Calculate class concentration of precipitation and nutrient load
  !!
  !>\b Reference ModelDescription Chapter Processes above ground (Atmospheric deposition of nitrogen and phosphorus)
  !---------------------------------------------------------------------
  SUBROUTINE set_class_precipitation_concentration_and_load(i,iluse,iveg,ns,area, &
                                      precorg,temp,prec,cprec,precload, &
                                      snowfall,rainfall,watertemp,icecover)

    USE MODVAR, ONLY : genpar, &
                       get_class_wet_deposition, &
                       i_t1,i_t2,i_in,i_sp, &
                       simulate, &
                       xobsi, &
                       xobsindex
    USE HYPEVARIABLES, ONLY : o_cprecT1, &
                              o_cprecIN, &
                              o_cprecSP, &
                              m_atmload, &
                              m_wetsp

    !Argument declarations
    INTEGER, INTENT(IN)   :: i            !<index of current subbasin
    INTEGER, INTENT(IN)   :: iluse        !<index of current landuse
    INTEGER, INTENT(IN)   :: iveg         !<index of current vegtype
    INTEGER,INTENT(IN)    :: ns           !<number of substances, array dimension
    REAL,   INTENT(IN)    :: area         !<class area (km2)
    REAL,   INTENT(IN)    :: precorg      !<precipitation from Pobs.txt (mm/timestep)
    REAL,   INTENT(IN)    :: temp         !<temperature of class (degree Celsius)
    REAL,   INTENT(IN)    :: prec         !<precipitation of class (mm/timestep)
    REAL,   INTENT(OUT)   :: cprec(ns)    !<concentration of precipitation (prec)
    REAL,   INTENT(OUT)   :: precload(ns) !<nutrient load of precipitation (kg/timestep)
    REAL,OPTIONAL,INTENT(IN) :: snowfall   !<snow fall (mm)
    REAL,OPTIONAL,INTENT(IN) :: rainfall   !<rain fall (mm)
    REAL,OPTIONAL,INTENT(IN) :: watertemp  !<temperature of water (C)
    REAL,OPTIONAL,INTENT(IN) :: icecover   !<ice cover (-)

    !Local variables
    INTEGER isubst
    LOGICAL xobsset(ns)
    REAL wetconc(ns)  !wet deposition concentration from AtmdepData.txt

    cprec = 0
    precload = 0.
    IF(ns==0) RETURN
    xobsset = .FALSE.
    
    !> \b Algorithm \n
    !>Set precipitation concentration from observed time series
    IF(simulate%substance(i_t1) .AND. xobsindex(o_cprecT1,i)>0)THEN
      cprec(i_t1) = xobsi(xobsindex(o_cprecT1,i))
      xobsset(i_t1)=.TRUE.
    ENDIF
    IF(simulate%substance(i_in))THEN
      IF(xobsindex(o_cprecIN,i)>0)THEN
        cprec(i_in) = xobsi(xobsindex(o_cprecIN,i))*1.E-3
        xobsset(i_in)=.TRUE.
      ENDIF
    ENDIF
    IF(simulate%substance(i_sp))THEN
      IF(xobsindex(o_cprecSP,i)>0)THEN
        cprec(i_sp) = xobsi(xobsindex(o_cprecSP,i))*1.E-3
        xobsset(i_sp)=.TRUE.
      ENDIF
    ENDIF
    IF(simulate%substance(i_t2))THEN
      xobsset(i_t2)=.TRUE.        !T2 is always set this way, never by AtmdepData.txt
    ENDIF
    
    !>Set wet deposition from input data (NEW)
    CALL get_class_wet_deposition(i,ns,iluse,iveg,wetconc)
    DO isubst = 1,ns
      IF(.NOT.xobsset(isubst))THEN
        cprec(isubst) = wetconc(isubst)
      ENDIF
    ENDDO
    
    !>Adjust concentration due to precipitation correction to keep wet deposition load
    IF(prec>0)THEN
      IF(genpar(m_atmload)==1)THEN
        IF(precorg/=prec)THEN
          cprec = cprec*(precorg/prec)
        ENDIF
      ENDIF    

      !>Temperature not adjusted, but set after adjustment. For lake/river depend on precipitation and ice cover
      IF(simulate%substance(i_t2))THEN
        cprec(i_t2)= MAX(0.,temp)     !Temp.conc. for landclasses [DG/JS Temp.model May-2013]
        IF(PRESENT(snowfall).AND.PRESENT(rainfall).AND.PRESENT(watertemp).AND.PRESENT(icecover))THEN  !lake/river
          CALL set_T2_concentration_in_precipitation_on_water(prec,temp,snowfall,rainfall,watertemp,icecover,cprec(i_t2))
        ENDIF
      ENDIF
      
      IF(ns>0) precload(:) = cprec(:)*prec*area
    ENDIF

    where (cprec < 0) 
      cprec=0
    ENDWHERE

  END SUBROUTINE set_class_precipitation_concentration_and_load

  !>Calculate temperature(T2) "concentration" in lake/river precipitation
  !>due to ice presens. Changes the default T2 concentration, set for class.
  !----------------------------------------------------------
  SUBROUTINE set_T2_concentration_in_precipitation_on_water(prec,temp,snowfall,rainfall,watertemp,icecover,cprect2)
  
    USE MODVAR, ONLY : cwater,cice,Lfreezing
   
    !Argument declaration
    REAL, INTENT(IN)      :: prec       !<precipitation
    REAL, INTENT(IN)      :: temp       !<air temperature
    REAL, INTENT(IN)      :: snowfall   !<snow fall
    REAL, INTENT(IN)      :: rainfall   !<rain fall
    REAL, INTENT(IN)      :: watertemp  !<temperature of water
    REAL, INTENT(IN)      :: icecover   !<ice cover
    REAL, INTENT(INOUT)   :: cprect2    !<T2 concentration of precipitation

    !This is now much more straight forward, using the fractional ice cover:
    ! Rainfall has always cprec = airtemp (but not lower than freezing point)
    ! Snowfall on the ice-free fraction has cprec = latentheat of freezing + sensible heat content
    ! Snowfall in the ice-covered fraction has cprec = laketemp
    IF(prec>0.)THEN
      !Rainfall temperature = max(0,air temp)
      cprect2 = rainfall * MAX(0.0,temp)
      
      !Snowfalltemp on ice   = watertemp (temporary), negative latent heat is added later when snow is melting in the ice routine
      !Snowfalltemp on water = airtemp + negative latent heat, taking into account diff spec.heat of ice and water
      cprect2 = cprect2 + snowfall * (watertemp * icecover + (MIN(temp,0.0)*(cice/cwater) - Lfreezing/cwater)*(1.-icecover))

      !Weighting by total precipitation
      cprect2 = cprect2/prec  
    ENDIF
    
  END SUBROUTINE set_T2_concentration_in_precipitation_on_water
  
  !>Calculate transformation factor for wind speed to different height than observations
  !---------------------------------------------------------------------
  SUBROUTINE calculate_class_wind_transformation_factor(windtrans)
  
  USE HYPEVARIABLES, ONLY : m_zwind,    &
                            m_zwish,    &
                            m_zpdh,     &
                            m_roughness
  USE MODVAR, ONLY : nclass,  &
                     genpar

  !Argument declarations
  REAL,ALLOCATABLE, INTENT(INOUT) :: windtrans(:)
  
  !Local parameters
  INTEGER j
  REAL lnz0,d0,zwind,zwish

  !Allocate wind transformation factor 
  IF(ALLOCATED(windtrans)) DEALLOCATE(windtrans)
  IF(.NOT.(ALLOCATED(windtrans))) ALLOCATE(windtrans(nclass))

  !Default: no transformation of wind 
  IF(genpar(m_zwind)==genpar(m_zwish) .OR. &
     (genpar(m_zwind)==0. .OR. genpar(m_zwish)==0.) .AND. &
     genpar(m_roughness)==0. .AND. genpar(m_zpdh)==0.)THEN
    windtrans = 1.
  ELSE
    !Wind is transformed to wanted height
    zwind = genpar(m_zwind)
    zwish = genpar(m_zwish)
    DO j = 1,nclass
      lnz0 = LOG(genpar(m_roughness))   !Could be land use parameters
      d0 = genpar(m_zpdh)
      windtrans(j) = (LOG(zwind-d0)-lnz0)/(LOG(zwish-d0)-lnz0)
    ENDDO
  ENDIF
  
  END SUBROUTINE calculate_class_wind_transformation_factor
  
  !>Calculate daylength based on latitude and julian day
  !---------------------------------------------------------------------
  SUBROUTINE calculate_daylength(jday,lat,length)
  
  USE MODVAR, ONLY : pi

  !Argument declarations
  INTEGER,INTENT(IN) :: jday     !<Current julian day number
  REAL, INTENT(IN)   :: lat      !<latitude
  REAL, INTENT(OUT)  :: length   !<day length (hours)
  
  !Local parameters
  REAL dec,a1

    dec = -23.45 * COS(pi*(REAL(jday) + 10.173)/182.61)
    a1 = MIN(1.,MAX(-1.,(SIN(lat * pi / 180.) * SIN(dec * pi/180)) / (COS(lat * pi / 180.) * COS(dec * pi/180.))))
    length = (1440. - 120./(15.*pi/180.) * ACOS(a1)) /60.
  
  END SUBROUTINE calculate_daylength
  
  !>Calculate wind direction (integer) and mean velocity from u(x) and v(y) components consistent with wind shelter classes
  !-------------------------------------------------------------------------------
  SUBROUTINE calculate_winddirspeed(numdir,uwind,vwind,winddirint,winddirdeg,windspeed) 
  
    USE MODVAR, ONLY : pi
    
    !Argument declarations
    INTEGER,INTENT(IN)  :: numdir     !<Number of wind direction classes
    REAL,INTENT(IN)     :: uwind      !<U-wind, west to east, m/s
    REAL,INTENT(IN)     :: vwind      !<V-wind, south to north, m/s
    INTEGER,INTENT(OUT) :: winddirint !<Wind direction, integer class, 1 to numdir (1=northerly)
    REAL,INTENT(OUT)    :: winddirdeg !<Wind direction, degree from north, 0 to 360
    REAL,INTENT(OUT)    :: windspeed  !<Wind speed in wind direction, m/s 

    !Local parameter
    REAL, PARAMETER :: radians_to_degree = 180./pi
    
    !>Wind speed in wind direction
    windspeed = SQRT(uwind ** 2 + vwind ** 2)
    
    IF(windspeed.EQ.0.)THEN
      !>Zero windspeed - nothing to do
      winddirdeg = 360.
      winddirint = 1
      RETURN
    ELSE
      !>Wind direction in degrees, clockwise from North (meteorological definition)
      IF(uwind.EQ.0.)THEN
        IF(vwind.GT.0.)THEN
          winddirdeg = 180. !southerly
        ELSE
          winddirdeg = 360.   !northerly
        ENDIF
      ELSEIF(vwind.EQ.0.)THEN
        IF(uwind.GT.0.)THEN
          winddirdeg = 270. !westerly
        ELSE
          winddirdeg = 90.  !easterly
        ENDIF
      ELSE
        IF(uwind.LT.0. .AND. vwind.LT.0.)THEN
          !NE winds
          winddirdeg = 90. - radians_to_degree * ATAN(vwind/uwind)
        ELSEIF(uwind.LT.0. .AND. vwind.GT.0.)THEN
          !SE winds
          winddirdeg = 90. + radians_to_degree * ATAN(-vwind/uwind)
        ELSEIF(uwind.GT.0. .AND. vwind.GT.0.)THEN
          !SW winds
          winddirdeg = 270. - radians_to_degree * ATAN(vwind/uwind)
        ELSE
          !NW winds
          winddirdeg = 270. + radians_to_degree * ATAN(-vwind/uwind)
        ENDIF
      ENDIF
      
      !>Wind direction class, 1 to numdir, where class 1 is centered to the north
      IF(numdir .EQ. 8)THEN
        IF((winddirdeg .GT. 337.5) .OR. (winddirdeg .LE. 22.5))THEN
          winddirint = 1
        ELSEIF((winddirdeg .GT. 22.5) .AND. (winddirdeg .LE. 67.5))THEN
          winddirint = 2
        ELSEIF((winddirdeg .GT. 67.5) .AND. (winddirdeg .LE. 112.5))THEN
          winddirint = 3
        ELSEIF((winddirdeg .GT. 112.5) .AND. (winddirdeg .LE. 157.5))THEN
          winddirint = 4
        ELSEIF((winddirdeg .GT. 157.5) .AND. (winddirdeg .LE. 202.5))THEN
          winddirint = 5
        ELSEIF((winddirdeg .GT. 202.5) .AND. (winddirdeg .LE. 247.5))THEN
          winddirint = 6
        ELSEIF((winddirdeg .GT. 247.5) .AND. (winddirdeg .LE. 292.5))THEN
          winddirint = 7
        ELSE
          winddirint = 8
        ENDIF
      ELSEIF(numdir .EQ. 4)THEN
        IF((winddirdeg .GT. 315.) .OR. (winddirdeg .LE. 45.))THEN
          winddirint = 1
        ELSEIF((winddirdeg .GT. 45.) .AND. (winddirdeg .LE. 135.))THEN
          winddirint = 2
        ELSEIF((winddirdeg .GT. 135.) .AND. (winddirdeg .LE. 225.))THEN
          winddirint = 3
        ELSE
          winddirint = 4
        ENDIF
      ELSE
        WRITE(6,*) 'Only 4 or 8 directions for wind is allowed, check parameter winddirs in par.txt'
        winddirint = 1
        RETURN
      ENDIF
    ENDIF
  END SUBROUTINE calculate_winddirspeed
  
  !>Calculate distribution of falling snow using Winstral coefficients (WSF)
  !>  
  !>  Original linear function (see HUVA reports Gustafsson et al, 2015):
  !>     relative snowfall = 1. + WSFscale * (WSF + WSFbias)
  !>
  !>  Alternative log-linear option (recommended), see manuscript Clemenzi et al 2020
  !>     relative snowfall = 10 ** (WSFscale * (WSF + WSFbias))
  !>
  !>  In both cases WSFscale can be scaled by landuse WSFluse.
  !>
  !>  The relative snowfall is normalized within each subbasin so that the subbasin
  !>  mean snowfall is preserved.
  !>
  !>  A placeholder for windspeed influence is also in the code, but still inactive.
  !------------------------------------------------------------------------
  SUBROUTINE calculate_snowfall_distribution(i, iluse, windspeed, wsf, sfdist)

    USE modvar, ONLY : nclass, &
                       classbasin, &
                       genpar, &
                       landpar, &
                       modeloption, &
                       p_snowfalldist
    USE hypevariables, only: m_wsfbias, m_wsfscale, m_wsfluse, m_sfdmax
    
    !Argument declarations
    INTEGER,INTENT(IN) :: i              !<subbasin number
    INTEGER,INTENT(IN) :: iluse(nclass)  !<land use
    REAL,INTENT(IN)    :: windspeed      !<windspeed
    REAL,INTENT(IN)    :: wsf(nclass)    !<Winstral factors for subbasin i, all classes
    REAL,INTENT(INOUT) :: sfdist(nclass) !<Snowfall distribution factor for class (area weighted for mass conservation)

    !Local variables
    REAL sum_dist       !Summation of factor(j)*slcArea
    REAL sum_area       !Summation of slcArea (nonforest and glacier)
    REAL a(nclass)      !class area
    REAL wsfwind        !windspeed scaling variable
    INTEGER j
    !variables for parameters
    REAL wsfbias,wsfscale,sfdmax,wsfluse

    !general parameters
    wsfscale = genpar(m_wsfscale)
    wsfbias  = genpar(m_wsfbias)
    sfdmax   = genpar(m_sfdmax)
    
    !skip this routine if scale factor (wsfscale) or max factor (sfrmax) equals 0. OR windspeed eq. 0 (=> no redistribution)
    IF(wsfscale.GT.0. .AND. sfdmax.GT.0. .AND. windspeed.GT.0.)THEN
    
      !class areas from subbasin
      a = classbasin(i,:)%part
    
      !initialize snowfall redistribution factors and summation area 
      sfdist(:) = 0; sum_dist = 0
      sum_area = 0.
      
      !Select snowfall distribution model (linear or log-linear)
      SELECT CASE(modeloption(p_snowfalldist))
      CASE(1) ! linear
        !Loop over slc-classes and calculate unscaled redistribution factor: sfdist = 1 + WSFscale * WSFluse * (WSF + WSFbias)
        DO j=1,nclass
          IF(a(j)>0)THEN
            !landuse scaling parameter
            wsfluse = landpar(m_wsfluse,iluse(j))
            !windspeed scaling parameter?
            wsfwind = 1.
            !sfred equation
            sfdist(j) = 1. + wsfwind * wsfscale * wsfluse * (wsf(j) + wsfbias)
            !Upper and lower limit for sfred 
            IF(sfdist(j) .LT. (1. - sfdmax))THEN
              sfdist(j) = 1. - sfdmax
            ELSEIF(sfdist(j) .GT. (1. + sfdmax))THEN
              sfdist(j) = 1. + sfdmax
            ENDIF
            IF(sfdist(j) .LT. 0.)sfdist(j) = 0. !may happen if sfrmax > 1
            !summation for weighting
            sum_dist = sum_dist + sfdist(j)*a(j)
            sum_area = sum_area + a(j)
          ENDIF
        ENDDO
      CASE(2) !log-linear
        !Loop over slc-classes and calculate unscaled redistribution factor: sfdist = 10 ** (WSFscale * WSFluse * (WSF + WSFbias))
        DO j=1,nclass
          IF(a(j)>0)THEN
            !landuse scaling parameter
            wsfluse = landpar(m_wsfluse,iluse(j))
            !windspeed scaling parameter?
            wsfwind = 1.
            !sfred equation
            sfdist(j) = 10**(wsfwind * wsfscale * wsfluse * (wsf(j) + wsfbias))
            !sfred truncation - upper and lower limit for sfred 
            IF(sfdist(j).GT.sfdmax)sfdist(j) = sfdmax
            IF(sfdist(j).LT.0.)sfdist(j)=0.
            !summation for weighting
            sum_dist = sum_dist + sfdist(j)*a(j)
            sum_area = sum_area + a(j)
          ENDIF
        ENDDO
      CASE DEFAULT
      END SELECT
      ! weight to make the areally weighted distribution factors within subbasin == 1
      DO j=1,nclass 
        IF (a(j)>0) THEN
          sfdist(j) = sfdist(j) * sum_area / sum_dist 
        ELSE
          sfdist(j) = 1.
        ENDIF      
      ENDDO
    ELSE
      sfdist=1.
    ENDIF
    RETURN
  END SUBROUTINE calculate_snowfall_distribution
  
  !>Calculate air pressure (kPa) as a function of elevation, FAO(7)
  !-------------------------------------------------------------------------------
  REAL FUNCTION airpressure_elevationfunction(elev)
     
    !Argument decalaration
    REAL, INTENT(IN) :: elev   !<elevation
     
    airpressure_elevationfunction = 101.3 * ((293. - 0.0065 * elev)/293. ) ** 5.26
     
  END FUNCTION airpressure_elevationfunction
  
  !>Calculate latent heat of vaporization (MJ kg-1) as a function of temperature
  !-------------------------------------------------------------------------------
  REAL FUNCTION latentheat_tempfunction(temp)

    !Argument decalaration
    REAL, INTENT(IN) :: temp !<temperature (C)
     
    latentheat_tempfunction = 2.501 - 0.002361 * temp  ![MJ/kg]
     
  END FUNCTION latentheat_tempfunction
  
  !> Calculate psychrometric constant (kPa C^-1) as a function of
  !! air pressure and latent heat of vaporization, FAO
  !-------------------------------------------------------------------------------
  REAL FUNCTION psychrometric_constant(airpressure,lambda)
  
    !Argument declarations
    REAL, INTENT(IN) :: airpressure !<air pressure [kPa]
    REAL, INTENT(IN) :: lambda      !<latent heat of vaporization [MJ/kg]
    
    !Parameter declaration
    REAL, PARAMETER :: cp = 0.001013  !specific heat of moist air at constant pressure (MJ kg^-1 C^-1)
 
    psychrometric_constant = cp * airpressure / (0.622 * lambda)
 
  END FUNCTION psychrometric_constant

  !>Calculates potential evaporation or uses a value supplied as input
  !
  !> \b Reference ModelDescription Processes above ground (Evaporation)
  !--------------------------------------------------------------
  SUBROUTINE calculate_potential_evaporation(i,j,temp,epot,radext,swrad,netrad,actvap,satvap,wind,epotsnow)
  
    USE MODVAR, ONLY : basin,classdata, &
                       landpar, &
                       genpar, &
                       get_current_petmodel, &
                       pi, &
                       xobsi,xobsindex, &
                       classbasin, &
                       current_time, &
                       timesteps_per_day
    USE HYPEVARIABLES, ONLY : o_reepot, &
                              m_ttmp,m_cevp, &
                              basincevpcorr, &
                              basincevpam, &
                              basincevpph, &
                              m_kc, &
                              m_krs,m_jhtadd,m_jhtscale,m_alfapt,m_fepotsnow

    !Argument declarations
    INTEGER, INTENT(IN) :: i      !<index of current subbasin
    INTEGER, INTENT(IN) :: j      !<index of current class 
    REAL, INTENT(IN)    :: temp   !<air temperature
    REAL, INTENT(OUT)   :: epot   !<potential evapotranspiration [mm/timestep]
    REAL, INTENT(IN)    :: radext !<extraterrestrial solar radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: swrad  !<downward shortwave radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: netrad !<net downward radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: actvap !<actual vapor pressure [kPa]
    REAL, INTENT(IN)    :: satvap !<saturated vapour pressure [kPa]
    REAL, INTENT(IN)    :: wind   !<wind speed [m/s]
    REAL, INTENT(OUT)   :: epotsnow !<potential evapotranspiration for snow mm/timestep
    
    !Local variables
    REAL tt       !threshold temperature for melting (C)
    REAL ce       !coefficient for evaporation (mm/C/timestep)
    REAL dsatvap  !Slope of saturation pressure curve [kPa/C]
    REAL gamma    !psychrometric constant
    REAL lambda   !latent heat of evaporation [MJ/kg]
    REAL airpressure  !atmospheric pressure [kPa]
    REAL kc       !crop coefficient used for the new PET functions
    REAL elev     !elevation
    REAL turbidity ! atmospheric turbidity
    REAL fepotsnow ! fraction of potential evaporation used for snow
    INTEGER current_petmodel
    
    !>\b Algorithm \n
    !>Set local parameters and corrections
    tt = landpar(m_ttmp,classdata(j)%luse)       !Threshold temperature for snow melt and evaporation
    ce = landpar(m_cevp,classdata(j)%luse)       !Coefficient for potential evaporation
    fepotsnow = landpar(m_fepotsnow,classdata(j)%luse)       !Coefficient for potential evaporation for snow
    ce = ce * (1 + basincevpam(i)*SIN(2.*pi*(current_time%dayno-1+REAL(current_time%tsofday)/REAL(timesteps_per_day)-basincevpph(i))/365.))
    
    !>Calculate additional input variables for the alternative PET functions
    current_petmodel = get_current_petmodel(i)
    IF(current_petmodel.GT.1)THEN
      dsatvap = deltasaturationpressure_function(temp)  !Slope of saturated vapour pressure curve, using mean temperature
      lambda = latentheat_tempfunction(temp)  !Latent heat of vaporization
      elev = basin(i)%elev+classbasin(i,j)%deltah
      airpressure = airpressure_elevationfunction(elev)  !Air pressure, assuming normal pressure at sea level
      gamma = psychrometric_constant(airpressure,lambda) !Psychrometric constant
      kc = landpar(m_kc(current_petmodel),classdata(j)%luse)  !PET-model specific Landuse scaling parameter, "crop coefficient"
      IF(kc==0) kc = landpar(m_kc(1),classdata(j)%luse)  !Default Landuse scaling parameter, "crop coefficient"
      IF(radext>0.)THEN
        turbidity = swrad / radext
      ELSE
        turbidity = 1.  !any value, PET(3) will be zero due to radext=0
      ENDIF
    ENDIF
      
    !>Calculate potential evaporation with the selected petmodel
    epot = 0.
    SELECT CASE(current_petmodel)
      CASE(0) !HYPE original model (with Xobs replacement, if available)
        IF(xobsindex(o_reepot,i)>0)THEN       
          epot = xobsi(xobsindex(o_reepot,i))
        ELSEIF(temp>tt)THEN
          epot = ce*(temp-tt)
        ELSE
          epot = 0.  
        ENDIF
      CASE(1) !HYPE original model (without Xobs replacement)
        IF(temp>tt)THEN
          epot = ce*(temp-tt)
        ELSE
          epot = 0.  
        ENDIF
      CASE(2) !Modified Jensen-Haise/McGuinness following Oudin et al (2005)
        !parameters suggested by Oudin et al, jhtadd = 5, jhtscale = 100
        epot = kc * MAX(0.,radext / (lambda) * (temp + genpar(m_jhtadd)) / genpar(m_jhtscale))
      CASE(3) !Hargreaves-Samani (known to overpredict in humid areas)
        ! The function is modified by DG to limit the "turbidity-factor" with the �ngstr�m formula:
        ! 
        !   The original Hargreaves function is:
        !     epot = 0.0023 * radext / (lambda*rho) * (Tmax-Tmin)^0.5 * (temp + 17.8)
        !   and the Hargreaves turbidity for estimating swrad = krs * (Tmax-Tmin)^0.5
        !
        !   Thus, by replacing (Tmax-Tmin)^2 with turbidity/krs, we get a reasonable limitation of the Hargreaves (tmax-Tmin) impact
        !  (furthermore, if Tmax-min was missing, we actually use the clearsky turbidity at this point)
        !
        ! also note that rho = 1 and excluded in equations below...
        epot = max(0.,kc * 0.0023 * radext /(lambda) * turbidity / genpar(m_krs) * (temp+17.8))
      CASE(4) ! Priestly Taylor (known to underpredict in arid and semi-arid areas)
        epot = max(0.,kc * genpar(m_alfapt) * dsatvap * netrad / (lambda * (dsatvap+gamma)))
      CASE(5) ! FAO Penman Monteith reference crop evapotranspiration
        epot = max(0., kc * ((0.408 * dsatvap * netrad + gamma*900./(temp+273.)*wind*(satvap-actvap))/(dsatvap+gamma*(1.+0.34*wind))))
      END SELECT
      !>Calculate potential evaporation for snow evaporation (sublimation)
      epotsnow = fepotsnow * epot
      
      !>Adjust potential evaporations with regional parameter cevpcorr
      epot = epot * basincevpcorr(i)
      epotsnow = epotsnow * basincevpcorr(i)

  END SUBROUTINE calculate_potential_evaporation

END MODULE ATMOSPHERIC_PROCESSES
