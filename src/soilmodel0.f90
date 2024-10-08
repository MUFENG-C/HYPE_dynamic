!> \file soilmodel0.f90
!> Contains module soilmodel_default.

!>HYPE default soil model
MODULE SOILMODEL_DEFAULT

  !Copyright 2012-2021 SMHI
  !
  !This file is part of HYPE.

  !HYPE is free software: you can redistribute it and/or modify it under
  !the terms of the Lesser GNU General Public License as published by
  !the Free Software Foundation, either version 3 of the License, or (at
  !your option) any later version. HYPE is distributed in the hope that
  !it will be useful, but WITHOUT ANY WARRANTY; without even the implied
  !warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
  !the Lesser GNU General Public License for more details. You should
  !have received a copy of the Lesser GNU General Public License along
  !with HYPE. If not, see <http://www.gnu.org/licenses/>.

  !-----------------------------------------------------------------------------------------
  !Procedures in this module
  !-----------------------------------------------------------------------------------------
  ! soilmodel_0
  !-----------------------------------------------------------------------------------------

  USE STATETYPE_MODULE
  USE MODVAR,  ONLY: missing_value,pi,   &
                     current_time, &
                     genpar,landpar,soilpar,  &
                     load,basin,path, &
                     numsubstances,nsub,maxsoillayers, &
                     classdata,soilthick,soildepth, &
                     i_t1,i_t2,i_in,i_on,i_sp,i_pp,i_oc,&
                     modeloption,p_deepgroundwater,p_infiltration, &
                     p_soilleakage,p_snowheat, &
                     slc_iwet,slc_owet,&
                     conduct,simulate, &
                     simulatesubstances
  USE HYPEVARIABLES, ONLY : wpmm,fcmm,epmm,pwmm, &
                            m_srrcs, &
                            m_wetsp,m_drypp, m_ponatm,m_cfrost,m_sfrost, &
                            m_sswcorr,m_immdep,m_iwdfrac,m_wdpar, &
                            m_ripz,m_rips,m_ripe, & 
                            soilmem,m_deepmem, &
                            epotdist
  USE GENERAL_WATER_CONCENTRATION, ONLY : remove_water, &
                                          error_remove_water
  USE ATMOSPHERIC_PROCESSES, ONLY : calculate_potential_evaporation,          &
                                    calculate_rain_snow_from_precipitation
  USE SOIL_PROCESSES, ONLY : calculate_snow, &
                             calculate_actual_soil_evapotranspiration, &
                             add_macropore_flow,    &
                             calculate_tile_drainage,  &
                             calculate_soil_runoff,    &
                             calculate_infiltration_flow_diversion, &
                             add_infiltration,   &
                             percolation, &
                             calculate_groundwater_table,  &
                             calculate_soiltemp,   &
                             calculate_frostdepth,  &
                             calculate_unfrozen_soil_water, &
                             calculate_liquid_water_fraction, &
                             calculate_soil_moisture_deficit
  USE NPC_SOIL_PROCESSES, ONLY : add_dry_deposition_to_landclass,  &
                                 calculate_plant,        &
                                 soil_substance_processes, &
                                 particle_processes_for_runoff, &
                                 local_diffuse_source,   &
                                 class_riparian_zone_processes
  USE TRACER_PROCESSES, ONLY : soil_tracer_processes, &
                               set_soil_T2_from_soiltemp_model
  USE IRRIGATION_MODULE, ONLY : apply_irrigation, &
                                calculate_irrigation_water_demand
  USE REGIONAL_GROUNDWATER_MODULE, ONLY : add_regional_groundwater_flow_to_soil

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: soilmodel_0
  !Private parameters, global in this module
  CHARACTER(LEN=27) :: errstring(1)  !error message for location of remove_water call
  PARAMETER (errstring = (/'surface runoff, soillayer 1'/))

CONTAINS

  !----------------------------------------------------------------
  !>\brief Default soilmodel for land classes
  !!Calculate snow and soil processes for a land class. 
  !----------------------------------------------------------------
  SUBROUTINE soilmodel_0(i,j,isoil,iluse,subid,pdayno,classarea,prec,cprec,temp, & 
       daylength,mintemp,maxtemp,swrad,  &
       radext,netrad,actvap,satvap,wind,rrcscorr,  &
       frozenstate,soilstate,miscstate,surfaceflow,csrunoff,crunoffd,   &
       cropuptakein,nitrif,denitrif,epot,gwat,frostdepth,    &
       smdef,evap,cevap,crunoff1,crunoff2,crunoff3,  &
       pwneed,irrappl,irrsources,snowfall,rainfall,cropsources,ruralaload,rgrwload,atmdepload,infiltrationflows,evapflows,  &
       runofflows,verticalflows,cverticalflows,horizontalflows,horizontalflows2,evapsnow,cruralflow,snowtemp,snowsurftemp)  

    INTEGER, INTENT(IN) :: i        !<index for current subbasin
    INTEGER, INTENT(IN) :: j        !<index for current class 
    INTEGER, INTENT(IN) :: isoil    !<index of soil type
    INTEGER, INTENT(IN) :: iluse    !<index of landuse
    INTEGER, INTENT(IN) :: subid    !<subbasin id
    INTEGER, INTENT(IN) :: pdayno    !<pseudo dayno for use in soil model subroutines
    REAL, INTENT(IN) :: classarea   !<class area [km2]
    REAL, INTENT(IN) :: prec        !<precipitation (mm/timestep)
    REAL, INTENT(IN) :: cprec(numsubstances)        !<concentration of precipitation
    REAL, INTENT(IN) :: temp        !<temperature
    REAL, INTENT(IN) :: daylength   !<day length (hours)
    REAL, INTENT(IN) :: mintemp     !<current daily min temperature (C)
    REAL, INTENT(IN) :: maxtemp     !<current daily max temperature (C)
    REAL, INTENT(IN) :: swrad       !<downward shortwave radiation [MJ/m2/day]
    REAL, INTENT(IN) :: radext      !<extraterrestrial solar radiation [MJ/m2/day]
    REAL, INTENT(IN) :: netrad      !<net downward radiation [MJ/m2/day]
    REAL, INTENT(IN) :: actvap      !<actual vapor pressure [kPa]
    REAL, INTENT(IN) :: satvap      !<saturated vapour pressure [kPa]
    REAL, INTENT(IN) :: wind        !<wind speed [m/s]
    REAL, INTENT(IN) :: rrcscorr    !<correction of recession coefficients
    TYPE(snowicestatetype),INTENT(INOUT)  :: frozenstate   !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states
    REAL, INTENT(OUT) :: surfaceflow(2)  !<saturated overflow and surface excess infilt
    REAL, INTENT(OUT) :: csrunoff(numsubstances)   !<concentration surface flow
    REAL, INTENT(OUT) :: crunoffd (numsubstances)  !<concentration tile runoff
    REAL, INTENT(OUT) :: cropuptakein  !<crop uptake of IN      
    REAL, INTENT(OUT) :: nitrif     !<nitrification
    REAL, INTENT(OUT) :: denitrif(maxsoillayers)   !<denitrification
    REAL, INTENT(OUT) :: epot       !<potential evaporation (mm/timestep)
    REAL, INTENT(OUT) :: gwat       !<groundwater table (m)
    REAL, INTENT(OUT) :: frostdepth   !<soil frost depth 
    REAL, INTENT(OUT) :: smdef        !<soil moisture deficit (mm)
    REAL, INTENT(OUT) :: evap       !<evapotranspiration (mm) weighted sum of evap(snowfree)+evapsnow
    REAL, INTENT(OUT) :: cevap(numsubstances)   !<concentration of evapotranspiration
    REAL, INTENT(OUT) :: crunoff1(numsubstances)   !<concentration of runoff from soil layer 1 (mg/L)
    REAL, INTENT(OUT) :: crunoff2(numsubstances)   !<concentration of runoff from soil layer 2 (mg/L)
    REAL, INTENT(OUT) :: crunoff3(numsubstances)   !<concentration of runoff from soil layer 3 (mg/L)
    REAL, INTENT(OUT) :: pwneed                    !<irrigation water demand for this classe (m3)
    REAL, INTENT(INOUT) :: irrappl   !<applied irrigation (mm), for summation basin output
    REAL, INTENT(INOUT) :: irrsources(numsubstances)     !<Load from irrigation to soil (kg/timestep)
    REAL, INTENT(OUT) :: snowfall     !<Precipitation as rain (mm)
    REAL, INTENT(OUT) :: rainfall     !<Precipitation as snow (mm)
    REAL, INTENT(INOUT) :: cropsources(2,numsubstances)  !<Load from fertiliser and resudues (kg/timestep)
    REAL, INTENT(OUT) :: ruralaload(numsubstances)   !<Load from rural households (kg/timestep)
    REAL, INTENT(INOUT) :: rgrwload(numsubstances)     !<Load from regional groundwater flow to soil (kg/timestep)
    REAL, INTENT(INOUT) :: atmdepload(numsubstances)   !<Load of atmospheric dry deposition (kg/timestep)
    REAL, INTENT(OUT) :: infiltrationflows(7)  !<several infiltration flows [mm]
    REAL, INTENT(OUT) :: evapflows(4)  !<evaporation from soillayers, snow and (glacier) [mm]
    REAL, INTENT(OUT) :: runofflows(7) !<different runoff flows:1-3=soil runoff sl 1-3,4-6=tile runoff sl 1-3,7=saturated surface runoff
    REAL, INTENT(OUT) :: verticalflows(6) !<vertical flows:1-2=percolation,3-4=upwelling due to rural,5-6=upwelling due to reg. grw flows
    REAL, INTENT(OUT) :: cverticalflows(2,numsubstances) !<concentration of vertical flows:1-2=percolation
    REAL, INTENT(OUT) :: horizontalflows(3)  !<horizontal flows:1-3=recieved rural load flow
    REAL, INTENT(INOUT) :: horizontalflows2(3,nsub) !<horizontal flows:1-3=division of regional groundwater flows to grwdown
    REAL, INTENT(OUT) :: evapsnow !<actual evaporation from snow
    REAL, INTENT(OUT) :: cruralflow(numsubstances) !<concentration of rural flow
    REAL, INTENT(OUT) :: snowtemp !<snowpack temperature
    REAL, INTENT(OUT) :: snowsurftemp !<snow surface temperature
    
    !Local variables
    INTEGER status  !error status of subroutine call
    LOGICAL wetland
    REAL sc       !coefficient for runoff recession surface runoff(no unit)
    REAL plantuptake(2,2)          !uptake of plant (what they want), kg NP/km2/timestep

    !Variables for class values
    REAL ginfilt,infilt   !gross infiltration (rain+snowrunoff), actual infiltration (after removed surfaceflow and macroporeflow)
    REAL cginfilt(numsubstances),cinfilt(numsubstances)   !concentration of infiltration
    REAL totalsurfaceflow
    REAL satoverflow    !surface flow from saturated soil
    REAL excessinfilt   !infiltration excess surface runoff 
    REAL macroflow
    REAL cmacroflow(numsubstances), cexcessinfilt(numsubstances)          !concentration in infiltration excess runoff and macropore flow
    REAL melt, cmelt(numsubstances)      !snow melt,concentration in snow melt water (same in snowrunoff)
    REAL cevapsnow(numsubstances) !concentration of snow evaporation
    REAL soilrunoff(maxsoillayers)    !soil runoff
    REAL csoilrunoff(numsubstances,maxsoillayers)    !concentration of soil runoff
    REAL cweights(maxsoillayers) ! weigths to calculate T2 conc. in tile drainage 
    REAL effsnowcov       !effective snow cover with respect to evaporation evap = evap(snowfree)*(1-effsnowcov)+evapsnow
    REAL epotsnow         !potential evaporation from snow
    REAL T1release !release of T1 from surface sources
    REAL runoffd      !tile drainage runoff
    REAL soilwater(maxsoillayers) !soilwater temporary variable
    REAL snowheat     !snow heat content, temporary variable
    REAL liqfrac(maxsoillayers) !liquid fraction of soil water
    REAL frozenvol(maxsoillayers) !<Frozen water volume (mm)
    REAL snowrunoff !<Runoff from snow layer, may be smaller than snowmelt!

    !Output, default values
    infiltrationflows = 0.
    runofflows = 0.
    evapflows = 0.
    verticalflows=0.
    cverticalflows=0.
    horizontalflows=0.
    nitrif = 0.; denitrif=0.
    cropuptakein = 0.
    pwneed = 0.
    ruralaload = 0.
    cruralflow = 0.
    
    !Check for special class wetland
    wetland = .FALSE.
    IF(j==slc_iwet.OR.j==slc_owet) wetland = .TRUE.

    !Locally defined variables for indata to be used in this subroutine
    liqfrac = 1.
    frozenvol = 0.
    T1release = 0.

!    IF(modeloption(p_soilleakage)==0)THEN
    IF(modeloption(p_soilleakage)/=1)THEN
      !Atmospheric deposition, add to soil or snow
      CALL add_dry_deposition_to_landclass(i,j,iluse,classdata(j)%vegtype, &
           classarea,atmdepload,frozenstate,soilstate)
    
      !Calculate plant growth (for uptake of nutrients)
      CALL calculate_plant(i,j,temp,daylength,plantuptake,miscstate)
    ENDIF

    !Irrigate the soil     
    IF(.NOT.wetland.AND.conduct%irrigation) CALL apply_irrigation(i,j,soilstate,irrappl,irrsources)

    !Potential evapotranspiration (before snow calculations, to calculate snow evaporation)
    CALL calculate_potential_evaporation(i,j,temp,epot,radext,swrad,netrad,actvap,satvap,wind,epotsnow)
    
    !Update snow pack; snowfall
    CALL calculate_rain_snow_from_precipitation(i,iluse,prec,temp,snowfall,rainfall) !form of precipitation
    
    !Update snow pack; melting, evaporation (sublimation)
    snowheat = 0.
    IF(modeloption(p_snowheat)>=1) snowheat = frozenstate%snowheat(j,i)
    CALL calculate_snow(i,j,basin(i)%subid,iluse,snowfall,cprec,frozenstate%snow(j,i),   &
              frozenstate%csnow(:,j,i),temp,melt,cmelt,swrad, &
              frozenstate%snowage(j,i),frozenstate%snowmax(j,i),frozenstate%snowdepth(j,i),  &
              frozenstate%snowcov(j,i),snowheat,epotsnow,evapsnow,cevapsnow,effsnowcov,snowtemp,snowsurftemp,frozenstate%snowliq(j,i),snowrunoff)
    evapflows(3) = evapsnow
    IF(modeloption(p_snowheat)>=1) frozenstate%snowheat(j,i) = snowheat
    
    !Gross infiltration
    ginfilt  = rainfall + snowrunoff 
    IF(simulatesubstances)THEN
      IF(ginfilt>0.)THEN
        cginfilt = (cmelt*snowrunoff + cprec*rainfall) / ginfilt
      ELSE 
        cginfilt = 0.
      ENDIF
    ENDIF

    !Calculate soil temperature and frost depth
    CALL calculate_soiltemp(maxsoillayers,temp,frozenstate%snowdepth(j,i),genpar(m_deepmem),soilmem(:,j),soilstate%deeptemp(j,i),soilstate%temp(:,j,i))
    IF(soilthick(3,j)>0)THEN
      CALL calculate_frostdepth(fcmm(:,j),wpmm(:,j),landpar(m_cfrost,iluse),soilpar(m_sfrost,isoil),   &
            soilstate%water(1,j,i)+soilstate%water(2,j,i)+soilstate%water(3,j,i),frostdepth,soilstate%temp(1:2,j,i),soilthick(:,j))
    ELSEIF(soilthick(2,j)>0)THEN
      CALL calculate_frostdepth(fcmm(:,j),wpmm(:,j),landpar(m_cfrost,iluse),soilpar(m_sfrost,isoil),   &
            soilstate%water(1,j,i)+soilstate%water(2,j,i),frostdepth,soilstate%temp(1:2,j,i),soilthick(:,j))
    ELSE
      CALL calculate_frostdepth(fcmm(:,j),wpmm(:,j),landpar(m_cfrost,iluse),soilpar(m_sfrost,isoil),   &
            soilstate%water(1,j,i),frostdepth,soilstate%temp(1:1,j,i),soilthick(:,j))
    ENDIF
    
    !Initial calculation of liquid water fraction and frozen water volume
    CALL calculate_unfrozen_soil_water(i,j,isoil,temp,wpmm(:,j),fcmm(:,j),epmm(:,j),soilstate,frozenvol,liqfrac)

    IF(modeloption(p_infiltration)==0 .OR. modeloption(p_infiltration)==1)THEN
      !Calculate and add infiltration to soil, including calculation of surface flow and macropore flow due to limited infiltration capacity 
      CALL calculate_infiltration_flow_diversion(i,j,isoil,wpmm(:,j),fcmm(:,j),epmm(:,j),ginfilt,cginfilt,temp,mintemp,maxtemp,  &
         infilt,cinfilt,excessinfilt,cexcessinfilt,macroflow,cmacroflow,frozenstate,soilstate)
      IF(wetland)THEN   !If wetland runoff from surface is calculted in routing part
        infilt = infilt + excessinfilt
        excessinfilt = 0.
      ENDIF
      CALL add_infiltration(i,j,iluse,infilt,cinfilt,soilstate)
      IF(ginfilt>0.)THEN
        infiltrationflows(1) = snowrunoff/ginfilt
      ENDIF
      infiltrationflows(2) = infilt
      infiltrationflows(3) = excessinfilt
      CALL calculate_liquid_water_fraction(i,j,soilstate,frozenvol,liqfrac) !update liqfrac
      
      !Percolation down through the soil layers, including N,P,OC-reduction
      CALL percolation(i,j,isoil,subid,wpmm(:,j),fcmm(:,j),pwmm(:,j),soilthick(:,j),liqfrac,verticalflows(1:2),cverticalflows,soilstate)
      CALL calculate_liquid_water_fraction(i,j,soilstate,frozenvol,liqfrac) !update liqfrac
    ENDIF

    !Calculate and remove evapotranspiration from the soil upper two layers
    soilwater = soilstate%water(:,j,i) !set soilwater variable to be used for evapotranspiration, soil runoff and tile runoff
    CALL calculate_actual_soil_evapotranspiration(i,j,2,soilwater,temp,epot, &
           wpmm(:,j),fcmm(:,j),epotdist(:,j),MIN(1.,MAX(1.-effsnowcov,0.)), &
           liqfrac,soilstate,evap,evapflows(1:2),cevap)
    CALL calculate_liquid_water_fraction(i,j,soilstate,frozenvol,liqfrac) !update liqfrac

    !Total evaporation, and Weighted average concentrations and potential evapotranspiration
    evap = evap + evapsnow
    IF(numsubstances.GT.0 .AND. evap.GT.0.) cevap(:) = (cevap(:)*(evap-evapsnow) + cevapsnow(:)*evapsnow)/evap
    epot = epot * (1.-effsnowcov) + epotsnow * effsnowcov
    
    !Calculate and remove soil runoff
    CALL calculate_soil_runoff(i,j,isoil,subid,soilwater,epmm(:,j), &
              classdata(j)%streamdepth,liqfrac,soilstate,soilrunoff,csoilrunoff)
    CALL calculate_liquid_water_fraction(i,j,soilstate,frozenvol,liqfrac) !update liqfrac
    runofflows(1:3) = soilrunoff(1:3)
    crunoff1 = csoilrunoff(:,1)
    crunoff2 = csoilrunoff(:,2)
    crunoff3 = csoilrunoff(:,3)

    !Calculate and remove runoff by tile or drainage pipe
    CALL calculate_tile_drainage(i,j,isoil,subid,soilwater,epmm(:,j),soildepth(:,j), &
              soilthick(:,j),classdata(j)%tiledepth,rrcscorr,liqfrac,soilstate,runoffd, &
              crunoffd,cweights)
    runofflows(4:6) = runoffd*cweights(1:3)
    CALL calculate_liquid_water_fraction(i,j,soilstate,frozenvol,liqfrac) !update liqfrac

    !Regional lateral groundwater flow from other subbasins
    !IF(modeloption(p_deepgroundwater)==1) CALL add_regional_groundwater_flow_to_soil(i,j,classarea,pwmm(:,j),soilstate,rgrwload,verticalflows(5:6),horizontalflows2)
    IF(modeloption(p_deepgroundwater)==1 .AND. .NOT.wetland) CALL add_regional_groundwater_flow_to_soil(i,j,classarea,pwmm(:,j),soilstate,rgrwload,verticalflows(5:6),horizontalflows2)
    
    !Load from local diffuse NP-source to the lowest soil layer
    IF(.NOT.wetland) CALL local_diffuse_source(i,j,pwmm(:,j),classarea,soilstate,ruralaload,verticalflows(3:4),horizontalflows(1:3),cruralflow)
    CALL calculate_liquid_water_fraction(i,j,soilstate,frozenvol,liqfrac) !update liqfrac
    
    IF(modeloption(p_infiltration)==2 .OR. modeloption(p_infiltration)==3)THEN
      !Calculate and add infiltration to soil, including calculation of surface flow and macropore flow due to limited infiltration capacity 
      CALL calculate_infiltration_flow_diversion(i,j,isoil,wpmm(:,j),fcmm(:,j),epmm(:,j),ginfilt,cginfilt,temp,mintemp,maxtemp,  &
         infilt,cinfilt,excessinfilt,cexcessinfilt,macroflow,cmacroflow,frozenstate,soilstate)
      IF(wetland)THEN   !If wetland runoff from surface is calculted in routing part
        infilt = infilt + excessinfilt
        excessinfilt = 0.
      ENDIF
      CALL add_infiltration(i,j,iluse,infilt,cinfilt,soilstate)
      IF(ginfilt>0.)THEN
        infiltrationflows(1) = snowrunoff/ginfilt
      ENDIF
      infiltrationflows(2) = infilt
      infiltrationflows(3) = excessinfilt
      CALL calculate_liquid_water_fraction(i,j,soilstate,frozenvol,liqfrac) !update liqfrac
    ENDIF
    
    IF(modeloption(p_infiltration)==2 .OR. modeloption(p_infiltration)==3)THEN
      !Percolation down through the soil layers, including N,P,OC-reduction
      CALL percolation(i,j,isoil,subid,wpmm(:,j),fcmm(:,j),pwmm(:,j),soilthick(:,j),liqfrac,verticalflows(1:2),cverticalflows,soilstate)
      CALL calculate_liquid_water_fraction(i,j,soilstate,frozenvol,liqfrac) !update liqfrac
    ENDIF

    !Surface runoff from saturated overland flow of uppermost soil layer
    sc=landpar(m_srrcs,iluse)*rrcscorr      !Runoff coefficient for surface runoff
    IF(sc>1.) sc = 1.
    satoverflow = MAX(sc * (soilstate%water(1,j,i)-pwmm(1,j)),0.)
    IF(wetland) satoverflow=0.
    csrunoff = 0.
    IF(satoverflow > 0.) THEN
      CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),satoverflow,soilstate%conc(:,1,j,i),status)
      IF(status.NE.0) CALL error_remove_water(errstring(1),subid,i,j)
      CALL calculate_liquid_water_fraction(i,j,soilstate,frozenvol,liqfrac) !update liqfrac
    ENDIF
    runofflows(7) = satoverflow

    !Total surfaceflow (saturated overland flow and excess infiltration)
    totalsurfaceflow = satoverflow + excessinfilt
    surfaceflow(1) = satoverflow
    surfaceflow(2) = excessinfilt
    IF(totalsurfaceflow > 0. .AND. simulatesubstances) THEN
      csrunoff(:) = (soilstate%conc(:,1,j,i) * satoverflow + cexcessinfilt(:) * excessinfilt) / totalsurfaceflow     !used for satoverflow and excessinfilt
    ENDIF

!    IF(modeloption(p_soilleakage)==0)THEN
    IF(modeloption(p_soilleakage)/=1)THEN
    !Erosion of particles with fastflow (surface flow and macropore flow) including delay in temporary storage.
    IF(conduct%simP.OR.conduct%simS) CALL particle_processes_for_runoff(i,j,isoil,iluse,pdayno,rainfall,totalsurfaceflow,   &
              macroflow,runoffd,runofflows(1)+runofflows(2)+runofflows(3)+runoffd+totalsurfaceflow,  &
              csrunoff,cmacroflow,crunoffd,crunoff1,    &
              crunoff2,crunoff3,frozenstate%snow(j,i),frostdepth,soilstate)  
    ENDIF

    !Add macropore water to soil layer with groundwater level (except the PP)
    CALL add_macropore_flow(i,j,macroflow,cmacroflow,epmm(:,j),pwmm(:,j), &
            soildepth(:,j),soilthick(:,j),infiltrationflows(4:6),soilstate)
    CALL calculate_liquid_water_fraction(i,j,soilstate,frozenvol,liqfrac) !update liqfrac

    IF(modeloption(p_infiltration)==0 .OR. modeloption(p_infiltration)==1)THEN
      !Second percolation down through the soil layers, including N,P,OC-reduction (limited to same maxperc)
      CALL percolation(i,j,isoil,subid,wpmm(:,j),fcmm(:,j),pwmm(:,j),soilthick(:,j),liqfrac,verticalflows(1:2),cverticalflows,soilstate)
      CALL calculate_liquid_water_fraction(i,j,soilstate,frozenvol,liqfrac) !update liqfrac
    ENDIF

    !Groundwater level and soil moisture deficit for output variable
    CALL calculate_groundwater_table(j,soilstate%water(:,j,i),  &
            epmm(:,j),soildepth(:,j),soilthick(:,j),gwat) 
    CALL calculate_soil_moisture_deficit(soilstate%water(:,j,i),wpmm(:,j),  &
            fcmm(:,j),soilthick(:,j),smdef) 

!    IF(modeloption(p_soilleakage)==0)THEN
    IF(modeloption(p_soilleakage)/=1)THEN
    !Soil transformation processes for substances             
    IF(numsubstances>0) CALL soil_substance_processes(i,j,iluse,isoil, &
         conduct,classarea,wpmm(:,j),fcmm(:,j),pwmm(:,j),plantuptake,  &
         soilthick(:,j),nitrif,denitrif,cropuptakein,cropsources,   &
         soilstate)
    CALL soil_tracer_processes(i,j,current_time,soilthick(:,j),soilstate,miscstate,ginfilt,T1release)
    
    !Add released T1 to surface runoff or top soil
    IF(simulate%substance(i_t1))THEN
      IF(T1release>0.)THEN
        IF(totalsurfaceflow>0.)THEN
          csrunoff(i_t1) = csrunoff(i_t1) + (T1release * (totalsurfaceflow/(infilt+totalsurfaceflow))) / totalsurfaceflow
        ENDIF
        soilstate%partT1(1,j,i) = soilstate%partT1(1,j,i) + T1release * (infilt/(infilt+totalsurfaceflow))
      ENDIF
    ENDIF
    ENDIF

    !Calculate irrigation water demand (for next timestep)
    IF(.NOT.wetland.AND.conduct%irrigation) CALL calculate_irrigation_water_demand(i,j,  &
         classarea,genpar(m_sswcorr),genpar(m_immdep),genpar(m_iwdfrac),  &
         genpar(m_wdpar),soilstate%water(:,j,i),wpmm(:,j),fcmm(:,j),  &
         epmm(:,j),epot,pwneed)

!    IF(modeloption(p_soilleakage)==0)THEN
    IF(modeloption(p_soilleakage)/=1)THEN
      !Riparian zone for OC
      IF(conduct%simC) CALL class_riparian_zone_processes(i,j,numsubstances,runofflows(1)+runofflows(2)+runofflows(3), &
           crunoff1(:),crunoff2(:),crunoff3(:),soilstate%oldgrw(j,i),landpar(m_ripz,iluse),        &
           soilstate%temp(:,j,i),genpar(m_ripe),genpar(m_rips),gwat,miscstate%temp10(i),miscstate%temp20(i),     &
           soilstate%water(1,j,i)+soilstate%water(2,j,i)+soilstate%water(3,j,i),SUM(wpmm(:,j)),SUM(pwmm(:,j)),soildepth(maxsoillayers,j))
         
      !Runoff T2 temperature dependent on soiltemp calculation [DG/JS Temp.model, May-2013]
      IF(simulate%substance(i_t2)) &
        CALL set_soil_T2_from_soiltemp_model(i,j,satoverflow,excessinfilt,cexcessinfilt,cweights,& 
                                             crunoff1,crunoff2,crunoff3,csrunoff,crunoffd,soilstate)
    ENDIF

  END SUBROUTINE soilmodel_0


END MODULE SOILMODEL_DEFAULT
