!> \file t_proc.f90
!> Contains module tracer_processes.

!>General substance (T1) and water temperature (T2) processes in HYPE
!Todo: Move T2 processes to this module.
MODULE TRACER_PROCESSES

  !Copyright 2014,2016-2020 SMHI
  !
  !This file is part of HYPE.
  !HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  !HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
  !You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

  !Used modules
  USE GENERAL_FUNCTIONS, ONLY : exponential_decay,sedimentation
  USE GENERAL_WATER_CONCENTRATION, ONLY : production_pool,retention_pool,new_concentration
  USE STATETYPE_MODULE, ONLY : snowicestatetype,soilstatetype,riverstatetype,lakestatetype,miscstatetype
  !Subroutines uses modvar, hypevariables
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------
  ! Private procedures 
  !----------------------------------------------
  !decay_of_tracer
  !decay_of_tracer_sorbedphase
  !sorption_of_tracer
  !sedimentation_resuspension_of_tracer
  !T1_crop_sources
  !T1_incorporation
  !release_from_source
  !----------------------------------------------
  PUBLIC :: add_tracer_point_source_to_river, &
            add_tracer_point_source_to_lake, &
            soil_tracer_processes, &
            tracer_processes_in_river, &
            tracer_processes_in_lake, &
            set_soil_T2_from_soiltemp_model
  
  CONTAINS

  !>Add T1-T2 load from point sources to river inflow
  !>
  !>\b Reference ModelDescription Chapter Water management (Point sources -
  !>Tracer T1 point sources)
  !-----------------------------------------------------------------
  SUBROUTINE add_tracer_point_source_to_river(i,itype,qin,cin)

    USE MODVAR, ONLY : simulate, &
                       i_t1, i_t2, &
                       npsused,psinfo,psload, &
                       tload, &
                       tloadexist, &
                       numsubstances, &
                       find_next_pointsource

    !Argument declarations
    INTEGER, INTENT(IN) :: i                        !<index of subbasin
    INTEGER, INTENT(IN) :: itype                    !<type of river
    REAL, INTENT(INOUT) :: qin                      !<flow into river (m3/s)
    REAL, INTENT(INOUT) :: cin(numsubstances)       !<concentration of flow into river (-)
    
    !Local variables
    INTEGER j   !row in tload
    INTEGER dim
    INTEGER ips,lastps
    REAL qadd
    REAL cadd(numsubstances)

    IF(.NOT.ALLOCATED(psinfo))THEN

      IF(ALLOCATED(tloadexist))THEN
        IF(tloadexist(i))THEN
          !Find load
          dim = SIZE(tload)
          IF(dim==0) RETURN
          DO j=1,dim
            IF(tload(j)%subindex==i)THEN
              IF(tload(j)%sw_code==itype*2-1)EXIT
            ENDIF
          ENDDO
          IF(j>dim) RETURN
    
          !Set source to be added 
          qadd = tload(j)%psvol   !m3/s
          IF(qadd>0)THEN
            cadd = 0.
            IF(simulate%substance(i_t1)) cadd(i_t1) = tload(j)%psconc
            IF(simulate%substance(i_t2)) cadd(i_t2) = tload(j)%pstemp

            !Add source to river inflow    
            IF(qin>0)THEN
              cin = (qin * cin + qadd * cadd)/(qin + qadd)
              qin = qin + qadd
            ELSE
              qin = qadd
              cin = cadd
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    
    ELSE
      !Use time dependent pointsources
      !Calculate source to be added to local river
      qadd = 0.;cadd = 0.
      lastps = 0
      DO
        CALL find_next_pointsource(i,lastps,npsused,psinfo,ips)
        IF(ips==0) EXIT
        lastps = ips
        IF(psinfo(ips)%sw_code/=1) CYCLE !add t1load to local river here only (main river added with add_point_sources_to_main_river)
        IF(psload(ips)%flow<=0.) CYCLE
        qadd = qadd + psload(ips)%flow   !m3/s
        cadd = cadd + psload(ips)%load   !g/s
      ENDDO
      IF(qadd>0)THEN
        !Add source to river      
        IF(qin>0)THEN
          cin = (qin * cin + cadd)/(qin + qadd)
          qin = qin + qadd
        ELSE
          qin = qadd
          cin = cadd / qadd   !mg/L
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE add_tracer_point_source_to_river

  !>Add T1-T2 load from point sources to lake inflow
  !>
  !>\b Reference ModelDescription Chapter Water management (Point sources -
  !>Tracer T1 point sources)
  !-----------------------------------------------------------------
  SUBROUTINE add_tracer_point_source_to_lake(i,itype,qin,cin)

    USE MODVAR, ONLY : i_t1, &
                       i_t2, &
                       npsused,psinfo,psload, &
                       tload, &
                       tloadexist, &
                       simulate, &
                       numsubstances, &
                       find_next_pointsource

    !Argument declarations
    INTEGER, INTENT(IN) :: i                        !<index of subbasin
    INTEGER, INTENT(IN) :: itype                    !<type of lake
    REAL, INTENT(INOUT) :: qin                      !<flow into lake (m3/s)
    REAL, INTENT(INOUT) :: cin(numsubstances)       !<concentration of flow into lake (-)
    
    !Local variables
    INTEGER j   !row in tload
    INTEGER dim
    INTEGER ips,lastps
    REAL qadd
    REAL cadd(numsubstances)

    IF(.NOT.ALLOCATED(psinfo))THEN

      !Periodically constant point source
      IF(ALLOCATED(tloadexist))THEN
        IF(tloadexist(i))THEN
          !Find load
          dim = SIZE(tload)
          IF(dim==0) RETURN
          DO j=1,dim
            IF(tload(j)%subindex==i)THEN
              IF(tload(j)%sw_code==itype*2)EXIT
            ENDIF
          ENDDO
          IF(j>dim) RETURN
    
          !Set source to be added
          qadd = tload(j)%psvol   !m3/s
          IF(qadd>0)THEN
            cadd = 0.
            IF(simulate%substance(i_t1)) cadd(i_t1) = tload(j)%psconc
            IF(simulate%substance(i_t2)) cadd(i_t2) = tload(j)%pstemp

            !Add source to lake inflow     
            IF(qin>0)THEN
              cin = (qin * cin + qadd * cadd)/(qin + qadd)
              qin = qin + qadd
            ELSE
              qin = qadd
              cin = cadd
            ENDIF
          ENDIF
        ENDIF
      ENDIF
        
    ELSE
      !Use time dependent pointsources
      !Calculate source to be added to lake
      qadd = 0.;cadd = 0.
      lastps = 0
      DO
        CALL find_next_pointsource(i,lastps,npsused,psinfo,ips)
        IF(ips==0) EXIT
        lastps = ips
        IF(psinfo(ips)%sw_code/=itype*2) CYCLE !add load to internal or outlet lake (skip river code)
        IF(psload(ips)%flow<=0.) CYCLE
        qadd = qadd + psload(ips)%flow   !m3/s
        cadd = cadd + psload(ips)%load   !g/s
      ENDDO
      IF(qadd>0)THEN
        !Add source to lake inflow      
        IF(qin>0)THEN
          cin = (qin * cin + cadd)/(qin + qadd)
          qin = qin + qadd
        ELSE
          qin = qadd
          cin = cadd / qadd   !mg/L
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE add_tracer_point_source_to_lake

  !>Calculate tracer processes in soil
  !
  !>\b Reference ModelDescription Chapter Tracer simulation (General 
  !>tracer (T1) - Processes)
  !------------------------------------------------------------------------------
  SUBROUTINE soil_tracer_processes(i,j,ct,thickness,soilstate,miscstate,infilt,T1release)    

    USE MODVAR, ONLY : maxsoillayers, &
                       simulate, &
                       i_t1, &
                       genpar, &
                       cropdata, &
                       TIMEINFORMATIONTYPE
    
    USE HYPEVARIABLES, ONLY : m_t1rel
    
    INTEGER, INTENT(IN) :: i                        !<index of subbasin 
    INTEGER, INTENT(IN) :: j                        !<index of class
    TYPE(TIMEINFORMATIONTYPE),INTENT(IN) :: ct      !<current time
    REAL, INTENT(IN)    :: thickness(maxsoillayers) !<thickness of soil layers
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate !<Soil states
    TYPE(miscstatetype),INTENT(INOUT)  :: miscstate !<Miscellaneous states
    REAL, INTENT(IN)    :: infilt                   !<infiltration (snowmelt and precip)
    REAL, INTENT(OUT)   :: T1release                !<release of T1 from surface sources
   
    IF(simulate%substance(i_t1))THEN
      !Add source to fields through cropdata
      IF(ALLOCATED(cropdata)) CALL T1_crop_sources(i,j,ct%date%year,miscstate%partT1sf(j,i))
      !Decay of T1 in different compartments
      CALL decay_of_tracer(1,soilstate%water(1,j,i),soilstate%conc(:,1,j,i)) 
      CALL decay_of_tracer_sorbedphase(1,soilstate%partT1(1,j,i))
      CALL decay_of_tracer_sorbedphase(1,miscstate%partT1sf(j,i))      
      !Release of surface T1 from rain/snowmelt
      IF(infilt>0. .AND. miscstate%partT1sf(j,i)>0.) THEN
        CALL release_from_source(miscstate%partT1sf(j,i),infilt,T1release,genpar(m_t1rel))
      ENDIF
      !Adsorption/desorption of T1
      CALL sorption_of_tracer(soilstate%water(1,j,i),soilstate%conc(:,1,j,i),soilstate%partT1(1,j,i),thickness(1)) 
      IF(thickness(2)>0) THEN
        CALL decay_of_tracer(1,soilstate%water(2,j,i),soilstate%conc(:,2,j,i))
        CALL decay_of_tracer_sorbedphase(1,soilstate%partT1(2,j,i))
        CALL sorption_of_tracer(soilstate%water(2,j,i),soilstate%conc(:,2,j,i),soilstate%partT1(2,j,i),thickness(2)) 
      ENDIF
      IF(thickness(3)>0) THEN
        CALL decay_of_tracer(1,soilstate%water(3,j,i),soilstate%conc(:,3,j,i))
        CALL decay_of_tracer_sorbedphase(1,soilstate%partT1(3,j,i))
        CALL sorption_of_tracer(soilstate%water(3,j,i),soilstate%conc(:,3,j,i),soilstate%partT1(3,j,i),thickness(3)) 
      ENDIF
      !Incorporation of T1 (ploughing down surface T1)
      IF(ALLOCATED(cropdata))THEN
        IF(thickness(2)==0)THEN
          CALL T1_incorporation(i,j,ct%dayno,1,soilstate,miscstate)
        ELSE
          CALL T1_incorporation(i,j,ct%dayno,2,soilstate,miscstate)
        ENDIF    
      ENDIF
    ENDIF    

  END SUBROUTINE soil_tracer_processes

  !>River processes of tracer T1
  !>
  !>\b Reference ModelDescription Chapter Tracer simulation (General 
  !>tracer (T1) - Processes)
  !-----------------------------------------------------------------
  SUBROUTINE tracer_processes_in_river(i,itype,riverarea,depth,transq,qbank,riverstate)

    USE MODVAR, ONLY : simulate, &
                       i_t1, &
                       genpar
    USE HYPEVARIABLES, ONLY : m_t1sedexp
    
    !Argument variables
    INTEGER, INTENT(IN) :: i               !<index of subbasin
    INTEGER, INTENT(IN) :: itype           !<index of river type (local = 1, main = 2)
    REAL, INTENT(IN)    :: riverarea       !<river area
    REAL, INTENT(IN)    :: depth           !<river depth (m)   
    REAL, INTENT(IN)    :: transq          !<flow out of translation box chain (m3/s)
    REAL, INTENT(IN)    :: qbank           !<bank full river flow
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    
    !Calculate T1 river processes
    IF(simulate%substance(i_t1))THEN
      CALL decay_of_tracer(1,riverstate%water(itype,i),riverstate%conc(:,itype,i))
      CALL sedimentation_resuspension_of_tracer(i,itype,riverarea,genpar(m_t1sedexp),riverstate,transq,qbank,depth)
      CALL decay_of_tracer_sorbedphase(1,riverstate%T1sed(itype,i))
    ENDIF

  END SUBROUTINE tracer_processes_in_river

  !>Lake processes of tracer T1
  !>
  !>\b Reference ModelDescription Chapter Tracer simulation (General 
  !>tracer (T1) - Processes)
  !-----------------------------------------------------------------
  SUBROUTINE tracer_processes_in_lake(i,itype,lakearea,pooladd_t1,lakestate)

    USE MODVAR, ONLY : simulate, &
                       i_t1, &
                       genpar
    USE HYPEVARIABLES, ONLY : m_t1sed
    
    !Argument declarations
    INTEGER,INTENT(IN) :: i             !<index of subbasin
    INTEGER,INTENT(IN) :: itype         !<index of lake type (ilake = 1, olake = 2)
    REAL, INTENT(IN)   :: lakearea      !<lake area (m2)
    REAL, INTENT(OUT)   :: pooladd_t1   !<T1 added to lake sediment
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate      !<Lake state

    
    !Calculate T1 lake processes
    IF(simulate%substance(i_t1))THEN
      CALL decay_of_tracer(1,lakestate%water(itype,i),lakestate%conc(:,itype,i))
      CALL sedimentation(lakestate%water(itype,i)*1.E-3,lakestate%conc(i_t1,itype,i),genpar(m_t1sed))
      IF(lakestate%water(itype,i)*1.E-3>genpar(m_t1sed))THEN
        pooladd_t1 = genpar(m_t1sed)*lakearea*lakestate%conc(i_t1,itype,i)*1.E-3   !m*m2*conc(mg/L)/1000->kg <-> U (conc(uU/L)
      ELSE
        pooladd_t1 = lakestate%water(itype,i)*1.E-3*lakearea*lakestate%conc(i_t1,itype,i)*1.E-3
      ENDIF
    ENDIF

  END SUBROUTINE tracer_processes_in_lake

  !>Decay processes of tracer T1
  !>
  !>\b Reference ModelDescription Chapter Tracer simulation (General 
  !>tracer (T1) - Processes)
  !-----------------------------------------------------------------
  SUBROUTINE decay_of_tracer(timesteps,vol,conc)

    USE MODVAR, ONLY : simulate, &
                       i_t1, &
                       genpar, &
                       timesteps_per_day, &
                       numsubstances
    USE HYPEVARIABLES, ONLY : m_t1expdec

    !Argument declarations
    INTEGER, INTENT(IN) :: timesteps            !<number of timesteps
    REAL, INTENT(IN)    :: vol                  !<volume of water 
    REAL, INTENT(INOUT) :: conc(numsubstances)  !<concentration of water (-)
    
    !Local variables
    REAL ts,left

    IF(simulate%substance(i_t1))THEN
      IF(vol>0. .AND. conc(i_t1)>0.)THEN
        IF(genpar(m_t1expdec)>0)THEN
          !Exponential decay with parameter t1expdec [days]
          ts = REAL(timesteps) / REAL(timesteps_per_day)  !time in days
          left = exponential_decay(ts,genpar(m_t1expdec))
          conc(i_t1) = conc(i_t1) * left
        ELSE
          !Other decay equation
        ENDIF
        
      ENDIF
    ENDIF
    
  END SUBROUTINE decay_of_tracer
  
  !>Decay of tracer in sorbed phase
  !>
  !>\b Reference ModelDescription Chapter Tracer simulation (General 
  !>tracer (T1) - Processes)
  !--------------------------------------------------------------------
  SUBROUTINE decay_of_tracer_sorbedphase(timesteps,pool)

    USE MODVAR, ONLY : genpar,  &
                       timesteps_per_day, &
                       simulate, &
                       i_t1
    USE HYPEVARIABLES, ONLY : m_t1expdec

    !Argument declarations
    INTEGER, INTENT(IN) :: timesteps            !<number of timesteps
    REAL, INTENT(INOUT) :: pool                 !<amount of T1
    
    !Local variables
    REAL ts,left

    IF(simulate%substance(i_t1))THEN
      IF(pool>0.)THEN
        IF(genpar(m_t1expdec)>0)THEN
          !Exponential decay with parameter t1expdec [days]
          ts = REAL(timesteps) / REAL(timesteps_per_day)  !time in days
          left = exponential_decay(ts,genpar(m_t1expdec))
          pool = pool * left
        ELSE
          !Other decay equation
        ENDIF
        
      ENDIF
    ENDIF
    
  END SUBROUTINE decay_of_tracer_sorbedphase  

  !>Adsorption of traced to soil particles
  !>
  !>\b Reference ModelDescription Chapter Tracer simulation (General 
  !>tracer (T1) - Processes)
  !--------------------------------------------------------------
  SUBROUTINE sorption_of_tracer(vol,conc,sorbed,layerthick)

    USE MODVAR, ONLY : simulate, &
                       i_t1, &
                       genpar, &
                       numsubstances
    USE HYPEVARIABLES, ONLY : m_T1freuc, bulkdensity
    

    !Argument declarations
    REAL, INTENT(IN)    :: vol                  !<volume of water (mm)
    REAL, INTENT(INOUT) :: conc(numsubstances)  !<concentration in water (-)
    REAL, INTENT(INOUT) :: sorbed               !<T1 sorbed phase
    REAL, INTENT(IN)    :: layerthick           !<thickness of soil layer (m)
    
    !Local variables
    REAL totalT1, sorbed_conc0, equi_solved,coeff, equi_sorbed, adsdes
    
    IF(simulate%substance(i_t1))THEN
      IF((vol>0. .AND. conc(i_t1)>0.) .OR. sorbed>0.)THEN
        totalT1 = sorbed + conc(i_t1) * vol  !Total amount of T1/km2 
        sorbed_conc0 = (sorbed / bulkdensity) / layerthick        !T1 /kg soil  
        coeff = genpar(m_T1freuc) * bulkdensity * layerthick
        equi_solved = totalT1 / (vol + coeff)
        equi_sorbed = genpar(m_T1freuc) * equi_solved
        
        !Calculate new pool and concentration, depends on the equilibrium concentration
        adsdes = (equi_sorbed - sorbed_conc0)* bulkdensity * layerthick
        sorbed = sorbed + adsdes   !new Pool T1
        conc(i_t1) = conc(i_t1) - (adsdes / vol)  !new liquid conc
      ENDIF
    ENDIF
    
  END SUBROUTINE sorption_of_tracer

  !>Sedimentation and resuspension of tracer in river flow
  !>
  !>\b Reference ModelDescription Chapter Tracer simulation (General 
  !>tracer (T1) - Processes)
  !---------------------------------------------------------------------------------------------
  SUBROUTINE sedimentation_resuspension_of_tracer(i,watertype,area,sedexppar,riverstate,riverq,qbank,depth)

    USE HYPEVARIABLES, ONLY : m_addsusp,m_hygeomk,m_hygeomm,m_qbank, &
                              m_sedch,m_suspconT1,m_suspexpT1,m_vpeak
    USE MODVAR, ONLY : basin,genpar,i_t1,modeloption,p_sedresusp,numsubstances

    !Argument declaration
    INTEGER, INTENT(IN)           :: i          !<index of current subbasin
    INTEGER, INTENT(IN)           :: watertype  !<river type (1=local, 2=main)
    REAL, INTENT(IN)              :: area       !<river surface area (m2)
    REAL, INTENT(IN)              :: sedexppar  !<sedimentation/resuspension parameter
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    REAL, INTENT(IN)              :: riverq     !<river discharge
    REAL, INTENT(IN)              :: qbank      !<Q bank full
    REAL, INTENT(IN)              :: depth      !<river depth (m) 

    !Local variables
    REAL, DIMENSION(1) :: T1pool           !pools in water
    REAL, DIMENSION(1) :: sedT1, resuspT1  !changes (1/d)
    REAL, DIMENSION(1) :: tempsed          !temporary variable for T1 in river sediment
    REAL sedresp, help, qbankcorr
    REAL sedch,riverv,vpeak,cmax(numsubstances)   !for model 2

    !Initial check
    IF(area <= 0.) RETURN
    sedch = genpar(m_sedch) !default value, general parameter
    
    !Calculate pool of water and sediment
    T1pool = (riverstate%water(watertype,i) * riverstate%conc(i_t1,watertype,i))* 1.0E-3 
    tempsed(1) = riverstate%T1sed(watertype,i)  

    !>Select current model for river resuspension/sedimentation
    SELECT CASE(modeloption(p_sedresusp))
      
    !>For model 0:
    CASE(0,1)  !original HYPE method, based om current flow related to bankful flow
      IF(sedexppar == 0) RETURN
      IF(qbank <= 0.) RETURN

      !>Calculate sedimentation and resuspension factor (per day)
      IF(modeloption(p_sedresusp)==0)THEN
        qbankcorr = 0.7 * qbank
      ELSE
        qbankcorr = genpar(m_qbank) * qbank
      ENDIF

      !Calculate sedimentation and resuspension
      resuspT1 = 0.
      sedT1 = 0.
      help = 0
      IF(qbankcorr-riverq/=0) help = help + ((qbankcorr-riverq)/qbankcorr)**sedexppar 
      IF(riverq>0) help = help - (riverq/qbankcorr)**sedexppar
      sedresp = MAX(-1., MIN(1.,help))
      IF(sedresp > 0) THEN !sedimentation
        sedT1 = sedresp * (riverstate%conc(i_t1,watertype,i) * MIN(riverstate%water(watertype,i),area * depth)) / 1.0E3  
      ELSE                 !resuspension
        resuspT1 = - sedresp * tempsed 
      ENDIF

      !Move sedimentation/resuspension to the water/sediment pool
      CALL retention_pool(1,T1pool,sedT1)              !sedT1 may change
      CALL production_pool(1,tempsed,sedT1)
      CALL retention_pool(1,tempsed,resuspT1)          !resuspT1 may change
      CALL production_pool(1,T1pool,resuspT1)

    !>For model 2: Simplified Bagnold Equation
    CASE(2)
      IF(basin(i)%channelfactor>=0.) sedch = basin(i)%channelfactor

      !>Calculate river peak velocity across channel
      riverv = genpar(m_hygeomk) * (riverq) ** genpar(m_hygeomm) ! Calculate river velocity
      vpeak = genpar(m_vpeak) * riverv ! Calculate peak channel velocity (m/s) = vpeak * average velocity (across the channel)
          
      !>Calculate amount of tracer to transfer
      cmax(i_t1) = genpar(m_suspconT1) * (vpeak ** genpar(m_suspexpT1)) * 1.0E6 ! Calculate maximum amount of sediment that can be transported and convert from kg/L to mg/L
      IF(riverstate%conc(i_t1,watertype,i) > cmax(i_t1)) THEN
        sedresp = (riverstate%conc(i_t1,watertype,i) - cmax(i_t1)) * MIN(riverstate%water(watertype,i),area * depth) / 1.0E3 ! amount of sediment (kg) to settle when concentration > max transport concentration
      ELSEIF(riverstate%conc(i_t1,watertype,i) < cmax(i_t1)) THEN
        sedresp = -1. * (cmax(i_t1) - riverstate%conc(i_t1,watertype,i)) * MIN(riverstate%water(watertype,i),area * depth) * sedch / 1.0E3 ! amount of sediment (kg) to resuspend when concentration < max transport concentration; multiply by -1. so it will be removed from pool
      ENDIF
      
      !>Transfer SS between sediment and water pools
      IF(sedresp>0.)THEN  !net sedimentation
        sedT1 = sedresp 
        CALL retention_pool(1,T1pool,sedT1)           !sedT1 may change
        CALL production_pool(1,tempsed,sedT1)
      ELSEIF(sedresp<0.)THEN  !net resuspension
        IF(-sedresp > tempsed(1)) THEN! #CB: resuspend full sediment pool and then use addsed (0.0 - 1.0) to control how much erosion in excess of temporary pool can occur
          resuspT1 = tempsed + ((-sedresp - tempsed) * genpar(m_addsusp))
          tempsed = 0.
        ELSE
          resuspT1 = -sedresp
          CALL retention_pool(1,tempsed,resuspT1)         !resuspT1 may NOT change, because already checked
        ENDIF
        CALL production_pool(1,T1pool,resuspT1)
      ENDIF
      
    END SELECT

    !>Update state variables
    riverstate%T1sed(watertype,i) = tempsed(1)
    CALL new_concentration(T1pool(1),riverstate%water(watertype,i)*1.0E-3,riverstate%conc(i_t1,watertype,i))

  END SUBROUTINE sedimentation_resuspension_of_tracer

  !>Add T1 from manure application to above-soil pool of sorbed T1
  !>
  !>\b Reference ModelDescription Chapter Tracer simulation (General 
  !>tracer (T1) - Sources)
  !---------------------------------------------------------------------------------
  SUBROUTINE T1_crop_sources(i,j,yearno,partT1sf)

    USE MODVAR, ONLY : find_croppart, &
                       in_season_period, &
                       cropdata
 
    !Argument declarations
    INTEGER, INTENT(IN) :: i                        !<subbasin
    INTEGER, INTENT(IN) :: j                        !<class
    INTEGER, INTENT(IN) :: yearno                   !<year
    REAL, INTENT(INOUT) :: partT1sf                 !pool of T1 on soil surface
    
    !Local variables
    INTEGER k,kcrop
    INTEGER applicperiod
    REAL common_t1add                                 
    REAL part                           !part of class area that is the current crop
    
    common_t1add   = 0.
    
    DO kcrop = 1, 2                             !main and secondary crop
      CALL find_croppart(i,j,kcrop,k,part)
      IF(k==0) CYCLE                            !no crop given, no uptake
      IF(cropdata(k)%T1amount==0) CYCLE
      IF(cropdata(k)%T1year==0 .OR. cropdata(k)%T1year == yearno) THEN
        IF(part>0)THEN
          !Calculate common fertiliser and manure use
          applicperiod = MAX(1,cropdata(k)%T1numberofdays)
          IF(in_season_period(cropdata(k)%T1day,applicperiod)) THEN
            common_t1add = common_t1add + part * cropdata(k)%T1amount / applicperiod
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    !Add T1 to surface pool
    IF(common_t1add>0.)THEN
      partT1sf = partT1sf + common_t1add
    ENDIF

  END SUBROUTINE T1_crop_sources

  !>Move sorbed T1 from above-soil pool to soil pool (tilling)
  !>
  !>\b Reference ModelDescription Chapter Tracer simulation (General 
  !>tracer (T1) - Sources)
  !-----------------------------------------------------------------
  SUBROUTINE T1_incorporation(i,j,jday,soillayers,soilstate,miscstate)
  
      USE MODVAR, ONLY : find_croppart,   &
                         cropdata
    
    !Argument declarations
    INTEGER, INTENT(IN) :: i                        !<subbasin
    INTEGER, INTENT(IN) :: j                        !<class
    INTEGER, INTENT(IN) :: jday                     !<day number of the year
    INTEGER, INTENT(IN) :: soillayers               !<number of soil layers T1 will be added to
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<soil states
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<miscellaneous states
    
    !Local variables
    INTEGER k,kcrop
    REAL part                           !part of class area that is the current crop
    
    DO kcrop = 1, 2                             !main and secondary crop
      CALL find_croppart(i,j,kcrop,k,part)
      IF(k==0) CYCLE
      IF(part>0)THEN
        IF(jday == cropdata(k)%T1daydown)THEN
          IF(soillayers == 1)THEN             
            soilstate%partT1(1,j,i) = soilstate%partT1(1,j,i) + miscstate%partT1sf(j,i) * cropdata(k)%T1down1
            miscstate%partT1sf(j,i) = miscstate%partT1sf(j,i) * (1. - cropdata(k)%T1down1)
          ELSE
            soilstate%partT1(1,j,i) = soilstate%partT1(1,j,i) + miscstate%partT1sf(j,i) * cropdata(k)%T1down1               
            soilstate%partT1(2,j,i) = soilstate%partT1(2,j,i) + miscstate%partT1sf(j,i) * cropdata(k)%T1down2
            miscstate%partT1sf(j,i) = miscstate%partT1sf(j,i) * (1. - cropdata(k)%T1down1 - cropdata(k)%T1down2)           
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  
  END SUBROUTINE T1_incorporation
  
  !>Release sorbed T1 to dissolved T1 on ground due to rain/snowmelt or surface runoff
  !>
  !>\b Reference ModelDescription Chapter Tracer simulation (General 
  !>tracer (T1) - Processes)
  !---------------------------------------------------------------------------------------
  SUBROUTINE release_from_source(sourceT1,qmm,releaseT1,relpar)
      
    !Argument declarations
    REAL, INTENT(INOUT) :: sourceT1                 !<pool of tracer
    REAL, INTENT(IN)    :: qmm                      !<water in mm in precipitation or surface runoff
    REAL, INTENT(OUT)   :: releaseT1                !<release of tracer from source
    REAL, INTENT(IN)    :: relpar                   !<release parameter
      
    releaseT1 = sourceT1 * (1. - exp(-relpar * qmm))
    sourceT1 = sourceT1 - releaseT1    
  
  END SUBROUTINE release_from_source
  
  !> T2 soil and runoff temperatures are determined by calculated soil temperature
  !---------------------------------------------------------------------------------------
  SUBROUTINE set_soil_T2_from_soiltemp_model(i,j,surfflow1,surfflow2,csurfflow2, &
                                               cweights,crunoff1,crunoff2,crunoff3, & 
                                               csrunoff,crunoffd,soilstate)
      
    USE MODVAR, ONLY : i_t2,numsubstances,maxsoillayers

    !Argument declarations
    INTEGER, INTENT(IN) :: i                      !<current subbasin
    INTEGER, INTENT(IN) :: j                      !<current class
    REAL, INTENT(IN)    :: surfflow1              !<runoff saturaed overland flow
    REAL, INTENT(IN)    :: surfflow2              !<runoff excess infiltration
    REAL, INTENT(IN)    :: csurfflow2(numsubstances) !<concentration of runoff excess infiltration
    REAL, INTENT(IN)    :: cweights(maxsoillayers)   !<soillayer with tile drainage
    REAL, INTENT(INOUT) :: crunoff1(numsubstances)  !<concentration of runoff
    REAL, INTENT(INOUT) :: crunoff2(numsubstances)  !<concentration of runoff
    REAL, INTENT(INOUT) :: crunoff3(numsubstances)  !<concentration of runoff
    REAL, INTENT(INOUT) :: csrunoff(numsubstances)  !<concentration of runoff
    REAL, INTENT(INOUT) :: crunoffd(numsubstances)  !<concentration of runoff
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<soil states

    !Local variables
    REAL trunofftemp(maxsoillayers)
    REAL totalsurfaceflow
    
    !Runoff temperature concentrations dependent on soiltemp calculation [DG/JS Temp.model, May-2013]
    
      !Calculate soil layer temperature (T2)
      trunofftemp(1) = AMAX1(0.,soilstate%temp(1,j,i))
      trunofftemp(2) = AMAX1(0.,soilstate%temp(2,j,i))
      trunofftemp(3) = AMAX1(0.,soilstate%temp(3,j,i))
      
      !Set T2 conc. in soil layers
      soilstate%conc(i_t2,1,j,i) = trunofftemp(1)
      soilstate%conc(i_t2,2,j,i) = trunofftemp(2)
      soilstate%conc(i_t2,3,j,i) = trunofftemp(3)
  
      !Set T2 conc. in runoff from soil layers
      crunoff1(i_t2) = trunofftemp(1)
      crunoff2(i_t2) = trunofftemp(2)
      crunoff3(i_t2) = trunofftemp(3)
      
      !Set T2 conc. in surface runoff, mixture of excess infiltration and saturated topsoil
      totalsurfaceflow = surfflow1 + surfflow2
      IF(totalsurfaceflow>0)THEN
        csrunoff(i_t2) = (trunofftemp(1) * surfflow1 + surfflow2 * csurfflow2(i_t2)) / totalsurfaceflow 
      ELSE
        csrunoff(i_t2) = 0.0
      ENDIF    
      
      !Set T2 conc. in tile drainage (weigthed average over tile drainage depth)
      crunoffd(i_t2) = cweights(1) * trunofftemp(1) + cweights(2) * trunofftemp(2) + &
                       cweights(3) * trunofftemp(3)
    
  END SUBROUTINE set_soil_T2_from_soiltemp_model

END MODULE TRACER_PROCESSES
