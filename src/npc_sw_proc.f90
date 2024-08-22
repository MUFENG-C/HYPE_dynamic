!> \file npc_sw_proc.f90
!> Contains module npc_surfacewater_processes.

!>Nitrogen, phosphorus and organic carbon processes in surface water in HYPE
MODULE NPC_SURFACEWATER_PROCESSES

  !Copyright 2012-2021 SMHI
  !
  !This file is part of HYPE.
  !HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  !HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
  !You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

  !Used modules
  USE STATETYPE_MODULE, ONLY : riverstatetype,lakestatetype
  USE GENERAL_WATER_CONCENTRATION
  USE GENERAL_FUNCTIONS
  !Subroutines also uses modvar and hypevariables
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------
  ! Private procedures 
  !----------------------------------------------
  ! denitrification_water 
  ! production_mineralisation
  ! lake_sedimentation
  ! river_sedimentation_resuspension
  ! calculate_lake_tpmean
  ! calculate_river_tpmean
  ! internal_lake_load 
  ! oc_production_mineralisation
  ! calculate_riverwetland_np 
  ! macrophyte_uptake
  !----------------------------------------------
  PUBLIC :: initiate_river_substance_state, &
       initiate_lake_substance_state, &
       add_deposition_to_river_as_load, &
       add_deposition_to_lake_as_load, &
       substance_processes_in_river, &
       substance_processes_in_lake, &
       sediment_pool_density, &
       add_diffuse_source_to_local_river, &
       add_point_sources_to_main_river, &
       calculate_river_wetland, &
       wetland_substance_processes
CONTAINS

  !>\brief Initiation river variables for nutrients and organic
  !>carbon simulations. Concentration (mg/L)
  !>
  !>\b Consequences Module hypevariables variable Qmax, Q2max,
  !>iQmax, iQ2max may be allocated and set.
  !-----------------------------------------------------------------
  SUBROUTINE initiate_river_substance_state(config,riverstate)

    USE HYPEVARIABLES, ONLY : Qmax,       &   !OUT
                              Q2max,      &   !OUT
                              iQmax,      &   !OUT
                              iQ2max,     &   !OUT
                              m_tpmean,   &
                              m_ldTPmean
    USE MODVAR, ONLY : STATECONFIGURATIONTYPE, &
                       nsub, &
                       basin, &
                       regpar, &
                       regiondivision, &
                       lakedatapar, &
                       lakedataparindex

    !Argument declarations
    TYPE(STATECONFIGURATIONTYPE), INTENT(IN) :: config !<state file configuration
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states

    !Local variables
    INTEGER isb

    !>\b Algorithm \n
    !>Allocate and initialize river sediment variables
    IF(config%qbank)THEN !For calculation of bankful flow
      IF(.NOT.ALLOCATED(Qmax))   ALLOCATE(Qmax(2,nsub))
      IF(.NOT.ALLOCATED(Q2max))  ALLOCATE(Q2max(2,nsub))
      IF(.NOT.ALLOCATED(iQmax))  ALLOCATE(iQmax(2,nsub))
      IF(.NOT.ALLOCATED(iQ2max)) ALLOCATE(iQ2max(2,nsub))
      riverstate%Q365 = 0.0001
      riverstate%Qdayacc = 0.000
      Qmax = 0.0001; Q2max = 0.0001
      iQmax = 365; iQ2max = 364
    ENDIF

    !>Set TPmean-variable if phosphorus is not calculated by HYPE
    IF((config%simN.OR.config%simC.OR.config%simS.OR.config%simSi).AND.(.NOT.config%simP))THEN
      DO isb = 1,nsub
        IF(basin(isb)%parregion(regiondivision(m_tpmean))>0) riverstate%TPmean(:,isb) = regpar(m_tpmean,basin(isb)%parregion(regiondivision(m_tpmean)))
        IF(ALLOCATED(lakedatapar)) riverstate%TPmean(1,isb) = lakedatapar(lakedataparindex(isb,1),m_ldtpmean)
        IF(ALLOCATED(lakedatapar)) riverstate%TPmean(2,isb) = lakedatapar(lakedataparindex(isb,2),m_ldtpmean)
      ENDDO
    ENDIF

  END SUBROUTINE initiate_river_substance_state

  !>Initiation lake for nutrients and organic carbon. Concentration in 
  !(mg/L)
  !!
  !>\b Consequences Module hypevariables variable slowlakeini
  !> may be allocated and set.
  !-----------------------------------------------------------------
  SUBROUTINE initiate_lake_substance_state(config,lakestate)

    USE HYPEVARIABLES, ONLY : m_tpmean,     &
                              m_tnmean,     &
                              m_tocmean,    &
                              m_ldtpmean,   &
                              m_ldtnmean,   &
                              m_ldtocmean,  &
                              m_iniT2,m_iniSi
    USE MODVAR, ONLY : STATECONFIGURATIONTYPE, &
                       basin, &
                       nsub, &
                       genpar, &
                       regpar, &
                       regiondivision, &
                       lakedatapar, &
                       lakedataparindex, &
                       i_in,i_on, &
                       i_pp,i_oc,i_t2,i_dsi

    !Argument declarations
    TYPE(STATECONFIGURATIONTYPE), INTENT(IN) :: config !<state file configuration
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    
    !Local variables
    INTEGER isb     !loop-variables

    !Initialize lake concentration (mg/L)
    IF(config%simN.OR.config%simP)THEN
      DO isb = 1,nsub
        IF(config%simP)THEN
          IF(ALLOCATED(lakedatapar))THEN
            lakestate%conc(i_pp,:,isb) = lakedatapar(lakedataparindex(isb,:),m_ldtpmean)
          ELSEIF(basin(isb)%parregion(regiondivision(m_tpmean))>0)THEN
            lakestate%conc(i_pp,:,isb) = regpar(m_tpmean,basin(isb)%parregion(regiondivision(m_tpmean)))
          ENDIF  
        ENDIF
        IF(config%simN)THEN
          IF(ALLOCATED(lakedatapar))THEN
            lakestate%conc(i_on,:,isb) = lakedatapar(lakedataparindex(isb,:),m_ldtnmean)*0.5
            lakestate%conc(i_in,:,isb) = lakedatapar(lakedataparindex(isb,:),m_ldtnmean)*0.5
          ELSEIF(basin(isb)%parregion(regiondivision(m_tnmean))>0)THEN
            lakestate%conc(i_on,:,isb) = regpar(m_tnmean,basin(isb)%parregion(regiondivision(m_tnmean)))*0.5
            lakestate%conc(i_in,:,isb) = regpar(m_tnmean,basin(isb)%parregion(regiondivision(m_tnmean)))*0.5
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    IF(config%simC)THEN
      DO isb = 1,nsub
        IF(ALLOCATED(lakedatapar))THEN 
          lakestate%conc(i_oc,:,isb) = lakedatapar(lakedataparindex(isb,:),m_ldtocmean)
        ELSEIF(basin(isb)%parregion(regiondivision(m_tocmean))>0)THEN
          lakestate%conc(i_oc,:,isb) = regpar(m_tocmean,basin(isb)%parregion(regiondivision(m_tocmean)))
        ENDIF
      ENDDO
    ENDIF

    IF(config%simSi)THEN
      DO isb = 1,nsub
        lakestate%conc(i_dsi,:,isb) = genpar(m_iniSi) !0.5mg/L
      ENDDO
    ENDIF

    IF(config%simT2)THEN
      DO isb = 1,nsub
        lakestate%conc(i_t2,:,isb) = genpar(m_iniT2)
      ENDDO
      lakestate%lowertemp(:,:) = genpar(m_iniT2)
      lakestate%uppertemp(:,:) = genpar(m_iniT2)
    ENDIF

    !Set TPmean-variable if phosphorus is not calculated by HYPE
    IF((config%simN.OR.config%simC.OR.config%simS.OR.config%simSi).AND..NOT.config%simP)THEN
      DO isb = 1,nsub
        IF(ALLOCATED(lakedatapar))THEN
          lakestate%TPmean(1,isb) = lakedatapar(lakedataparindex(isb,1),m_ldtpmean)
          lakestate%TPmean(2,isb) = lakedatapar(lakedataparindex(isb,2),m_ldtpmean)
        ELSEIF(basin(isb)%parregion(regiondivision(m_tpmean))>0)THEN
          lakestate%TPmean(:,isb) = regpar(m_tpmean,basin(isb)%parregion(regiondivision(m_tpmean)))
        ENDIF  
      ENDDO
    ENDIF

  END SUBROUTINE initiate_lake_substance_state

  !>\brief Calculate atmospheric deposition of N and P and add it
  !>to lakewater as load
  !>
  !>\b Reference ModelDescription Chapter Processes above ground (Atmospheric deposition of nitrogen and phosphorus)
  !----------------------------------------------------------------------------
  SUBROUTINE add_deposition_to_lake_as_load(i,iluse,pooltype,veg,areaij,sourcedry,lakestate)

    USE MODVAR, ONLY : numsubstances, &
                       set_class_deposition

    !Argument declarations
    INTEGER, INTENT(IN) :: i                          !<index of subbasin
    INTEGER, INTENT(IN) :: iluse                      !<index of land use
    INTEGER, INTENT(IN) :: pooltype                   !<laketype: 1=ilake, 2=olake
    INTEGER, INTENT(IN) :: veg                        !<vegetation index
    REAL, INTENT(IN)    :: areaij                     !<classarea (km2)
    REAL, INTENT(OUT)   :: sourcedry(numsubstances)   !<dry deposition (kg/timestep)    
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate    !<Lake state

    !Local variables
    REAL :: sourcedd(numsubstances)

    sourcedry = 0.

    !Prepare deposition
    sourcedd = 0.
    CALL set_class_deposition(i,numsubstances,iluse,veg,sourcedd) !kg/km2/timestep

    !Add deposition to lake
    IF(lakestate%water(pooltype,i)>0)THEN
      IF(SUM(sourcedd)>0.) CALL add_source_to_water(lakestate%water(pooltype,i),numsubstances,lakestate%conc(:,pooltype,i),sourcedd)
    ENDIF

    !Calculate atmospheric deposition loads (kg/timestep)
    sourcedry = sourcedd*areaij

  END SUBROUTINE add_deposition_to_lake_as_load
  
  !>\brief Calculate atmospheric deposition of N and P and add it
  !>to river water components as load
  !>
  !\b Reference ModelDescription Chapter Processes above ground (Atmospheric 
  !>deposition of nitrogen and phosphorus)
  !----------------------------------------------------------------------------
  SUBROUTINE add_deposition_to_river_as_load(i,iluse,pooltype,veg,areaij,sourcedry,riverstate)

    USE MODVAR, ONLY : numsubstances, &
                       set_class_deposition
    USE HYPEVARIABLES, ONLY : ttpart,ttstep

    !Argument declarations
    INTEGER, INTENT(IN) :: i                          !<index of subbasin
    INTEGER, INTENT(IN) :: iluse                      !<index of land use
    INTEGER, INTENT(IN) :: pooltype                   !<rivertype: 1=lriver, 2=mriver
    INTEGER, INTENT(IN) :: veg                        !<vegetation index
    REAL, INTENT(IN)    :: areaij                     !<classarea (km2)
    REAL, INTENT(OUT)   :: sourcedry(numsubstances)   !< dry deposition (kg/timestep)    
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River state

    !Local variables
    INTEGER l
    REAL sourcedd(numsubstances)
    REAL totvol
    REAL deposition(numsubstances)
    REAL fractionw,fractionpart,notused   !fractions of volume in river compartments
    REAL,ALLOCATABLE :: fractionqueue(:)  !fractions of volume in river compartments

    !>\b Algorithm \n
    sourcedry = 0.
    
    !>Prepare deposition
    CALL set_class_deposition(i,numsubstances,iluse,veg,sourcedd)   !kg/km2/timestep (for NP)
    
    sourcedry = sourcedd*areaij   !kg/timestep, OUTPUT
    sourcedd = sourcedry*1000.    !g=mg/L*m3, for adding
    
    IF(SUM(sourcedd)>0.)THEN
      !>Calculate fractions to be added to river water compartments
      ALLOCATE(fractionqueue(ttstep(pooltype,i)))
      CALL calculate_water_fractions(ttstep(pooltype,i),ttpart(pooltype,i),riverstate%water(pooltype,i),0.,riverstate%qqueue(1:ttstep(pooltype,i)+1,pooltype,i),totvol,fractionw,notused,fractionqueue,fractionpart)
      IF(totvol<=0)THEN
        IF(ALLOCATED(fractionqueue)) DEALLOCATE(fractionqueue)
        sourcedry = 0.
        RETURN
      ENDIF

      !>Add deposition of Inorg-N and PartP/SoluteP to river watercourse compartments
      IF(fractionw>0)THEN
        deposition = fractionw*sourcedd
        CALL add_source_to_water(riverstate%water(pooltype,i),numsubstances,riverstate%conc(:,pooltype,i),deposition)
      ENDIF
      DO l = 1,ttstep(pooltype,i)
        IF(fractionqueue(l)>0)THEN
          deposition = fractionqueue(l)*sourcedd
          CALL add_source_to_water(riverstate%qqueue(l,pooltype,i),numsubstances,riverstate%cqueue(:,l,pooltype,i),deposition)
        ENDIF
      ENDDO
      IF(fractionpart>0)THEN
        l = ttstep(pooltype,i) + 1
        deposition = fractionpart*sourcedd
        CALL add_source_to_water(riverstate%qqueue(l,pooltype,i),numsubstances,riverstate%cqueue(:,l,pooltype,i),deposition)
      ENDIF
    ENDIF
    IF(ALLOCATED(fractionqueue)) DEALLOCATE(fractionqueue)

  END SUBROUTINE add_deposition_to_river_as_load  
  
  !>\brief Calculate nutrient, organic carbon, silica and sediment processes in river 
  !> This include denitrification, mineralisation, primary production,
  !!sedimentation, resuspension, exchange with sediment
  !>
  !>\b Reference ModelDescription Chapter Nitrogen and phosphorus processes in rivers and lakes (Common things in lakes and river,
  !> Denitrification, Primary production and mineralization, and Sedimentation/Resuspension) and  Chapter Organic carbon (River and Lakes - Primary production and mineralization)
  !---------------------------------------------------------------------------
  SUBROUTINE substance_processes_in_river(i,itype,area,depth,transq,Qbank,maxSSconc,sedresSS,riverstate)   

    USE HYPEVARIABLES, ONLY : m_denitwr, m_denitwrl,m_incorr, &
                              m_hsatINwater,m_ldwprodn,m_ldwprodp,m_ldwprodc, &
                              m_wprodsi,m_sitmpexp,m_hsatTP,m_sedexp, &
                              m_limsedpp,m_plimsi, &
                              m_muptnriver,m_muptpriver,m_muptdepriver

    USE MODVAR, ONLY : basin, &
                       conduct, &
                       lakedatapar, &
                       lakedataparindex, &
                       regiondivision, &
                       genpar,regpar, &
                       simulate, &
                       i_in,i_sp,i_oc

    !Argument declarations
    INTEGER, INTENT(IN) :: i         !<index of current subbasin
    INTEGER, INTENT(IN) :: itype     !<river type (local or main)
    REAL, INTENT(IN)    :: area      !<river area (m2)
    REAL, INTENT(IN)    :: depth     !<river depth (m)   
    REAL, INTENT(IN)    :: transq    !<flow out of translation box chain (m3/s)
    REAL, INTENT(IN)    :: qbank     !<bank full river flow
    REAL, INTENT(OUT)   :: maxSSconc !<river SS maximum transport concentration (mg/L)
    REAL, INTENT(OUT)   :: sedresSS  !<river sedimentation/resuspension of SS (kg/timestep)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states

    !Local parameters
    INTEGER, PARAMETER :: systemtype = 2    !river system
    
    !Local variables
    REAL :: incorr    !correction factor
    REAL :: denpar    !model parameter denitrification rate (kg/m2/day)
    REAL :: hsatINpar !model parameter half saturation IN (mg/L)
    REAL :: prodNpar  !model parameter production ON 
    REAL :: prodPpar  !model parameter production PP
    REAL :: prodCpar  !model parameter production OC
    REAL :: prodSipar  !model parameter production Algal Si
    REAL :: tmpexpSipar   !model parameter temperature dependence production Algal Si
    REAL :: hsatTPpar !model parameter half saturation TP (mg/L)
    REAL :: sedexppar !sedimentation/resuspension parameter (mg/L)
    REAL :: limpppar  !limitation of sedimentation parameter (mg/L)
    REAL :: plimsipar   !P limitation of Si production parameter (mg/L)
    REAL :: muptNpar   !macrophyte N uptake
    REAL :: muptPpar   !macrophyte P uptake
    REAL :: muptdeppar !depth above which macrophyte uptake can occur

    !Default output 
    maxSSconc = 0.;sedresSS = 0.
    
    !Short notation for parameters dependent on regional corrections
    IF(basin(i)%parregion(regiondivision(m_incorr))>0)THEN
      incorr = (1. + regpar(m_incorr,basin(i)%parregion(regiondivision(m_incorr))))     !Correction of inorganic nitrogen
    ELSE
      incorr = 1.
    ENDIF
    IF(itype==1)THEN
      denpar = (2.-incorr)*genpar(m_denitwrl)   !denitrification opposite correction, aka 1-incorrpar=2-(1+incorrpar)
    ELSE
      denpar = (2.-incorr)*genpar(m_denitwr)
    ENDIF
    hsatINpar = genpar(m_hsatINwater)
    prodNpar = lakedatapar(lakedataparindex(i,itype),m_ldwprodn)
    prodPpar = lakedatapar(lakedataparindex(i,itype),m_ldwprodp)
    prodCpar = lakedatapar(lakedataparindex(i,itype),m_ldwprodc)
    prodSipar = genpar(m_wprodsi)
    tmpexpSipar = genpar(m_sitmpexp)
    hsatTPpar = genpar(m_hsatTP)
    sedexppar = genpar(m_sedexp)
    limpppar = genpar(m_limsedpp)
    plimsipar = genpar(m_plimsi)
    muptNpar = genpar(m_muptnriver)
    muptPpar = genpar(m_muptpriver)
    muptdeppar = genpar(m_muptdepriver)

    !Calculate the nutrient processes
    IF(area>0)THEN
      IF(i_sp>0) CALL calculate_river_tpmean(i,itype,riverstate)
      
      IF(simulate%substance(i_in)) CALL denitrification_water(i,itype,systemtype,area,denpar,hsatINpar,RIVERSTATE=riverstate)  !denitrification in rivers
      CALL production_mineralisation(i,itype,systemtype,area,prodNpar,prodPpar,prodSipar,tmpexpSipar,hsatTPpar,limpppar,plimsipar,RIVERSTATE=riverstate,DEPTH=depth) !mineraliation and primary production of NP in rivers  
      IF(conduct%simC) CALL oc_production_mineralisation(systemtype,area,prodCpar,hsatTPpar,limpppar,       &
              riverstate%water(itype,i),riverstate%conc(i_oc,itype,i),  &
              riverstate%temp(itype,i),riverstate%TPmean(itype,i),      &
              riverstate%temp10(itype,i),riverstate%temp20(itype,i),depth) 
      CALL macrophyte_uptake(i,itype,systemtype,area,muptNpar,muptPpar,hsatTPpar,muptdeppar,limpppar,RIVERSTATE=riverstate)
      IF(conduct%simP.OR.conduct%simS) CALL river_sedimentation_resuspension(i,itype,area,sedexppar,transq,Qbank,depth,maxSSconc,sedresSS,riverstate) !sedimentation and resuspension of PP and SS in rivers  
       !    No sediment SRP exchange. The concentration smoothed with deadvolume
    ENDIF

  END SUBROUTINE substance_processes_in_river


  !>\brief Calculate nutrient, organic carbon, silica and sediment processes in lake: 
  !!denitrification, mineralisation, primary production, macrophyte uptake,
  !!sedimentation, internal load
  !>
  !>\b Reference ModelDescription Chapter Nitrogen and phosphorus processes in rivers and lakes (Common things in lakes and river,
  !> Denitrification, Primary production and mineralization, Sedimentation/Resuspension and Internal load) and Chapter Organic carbon (River and Lakes) 
  !------------------------------------------------------------------
  SUBROUTINE substance_processes_in_lake(i,itype,area,lakestate,resuspdown, &
                                pooladd_ts,pooladd_tn,pooladd_tp,pooladd_oc, &
                                pooladd_si,poolnet_ts)

    USE HYPEVARIABLES, ONLY : m_lddenitwl,m_incorr,m_oncorr, &
                              m_hsatINwater,m_ldwprodn,m_ldwprodp,m_ldwprodc, &
                              m_ldwprodsi,m_sitmpexp,m_hsatTP, &
                              m_limsedpp,m_plimsi,m_limsedon,m_limsedss, &
                              m_ldsedon,m_ldsedpp,m_ldsedoc,m_ldsedss,m_ldsedsi,m_sedae, &
                              m_ldmuptn,m_ldmuptp,m_muptdep

    USE MODVAR, ONLY : basin, &
                       conduct, &
                       lakedatapar, &
                       lakedataparindex, &
                       numsubstances, &
                       regiondivision, &
                       genpar,regpar, &
                       i_on,i_pp,i_oc,i_ss,i_ae,i_asi

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<index of subbasin
    INTEGER, INTENT(IN) :: itype      !<lake type (ilake or olake)
    REAL, INTENT(IN)    :: area       !<lake area (m2)
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    REAL, INTENT(OUT)   :: resuspdown(numsubstances)  !<resuspended substance to be transported downstream (kg?/ts)
    REAL, INTENT(OUT)   :: pooladd_ts !<sediment added to lake sediment pool (SS+AE) (kg)
    REAL, INTENT(OUT)   :: pooladd_tn !<nitrogen added to lake sediment (ON) (kg)
    REAL, INTENT(OUT)   :: pooladd_tp !<phosphorus added to lake sediment (PP) (kg)
    REAL, INTENT(OUT)   :: pooladd_oc !<org. carbon added to lake sediment (OC) (kg)
    REAL, INTENT(OUT)   :: pooladd_si !<silica added to lake sediment (ASi) (kg)
    REAL, INTENT(OUT)   :: poolnet_ts !<net change in sediment in lake (SS+AE) (kg)
    
    !Local parameters
    INTEGER, PARAMETER :: systemtype = 1    !lake

    !Local variables
    REAL incorr,oncorr
    REAL denpar     !model parameter denitrification rate (kg/m2/day)
    REAL hsatINpar  !model parameter half saturation IN (mg/L)
    REAL prodNpar   !model parameter production ON 
    REAL prodPpar   !model parameter production PP 
    REAL prodCpar   !model parameter production OC
    REAL prodSipar  !model parameter production Algal Si
    REAL tmpexpSipar   !model parameter temperature dependence production Algal Si
    REAL sedonpar   !ON sedimentation rate  (lakes)
    REAL sedpppar   !PP sedimentation rate  (lakes)
    REAL sedocpar   !OC sedimentation rate  (lakes)
    REAL sedsspar   !SS sedimentation rate  (lakes)
    REAL sedsipar   !Algal Si sedimentation rate  (lakes)    
    REAL hsatTPpar  !model parameter half saturation TP (mg/L)
    REAL limpppar   !limitation of sedimentation PP parameter (mg/L)
    REAL plimsipar  !P limitation of Si production parameter (mg/L)
    REAL ss_frac    !Fraction of siltation depth that is ss (before sedimentation and flushing)
    REAL pooladd_ss, pooladd_ae !amount of substance added to lake sediment
    REAL resusptot(numsubstances)    ! #CB: amount of substance removed from lake sediment

    !>\b Algorithm \n
    pooladd_ts = 0. !default output
    pooladd_tn = 0.
    pooladd_tp = 0.
    pooladd_oc = 0.
    pooladd_si = 0.
    poolnet_ts = 0.
    
    !>Short notation for parameters dependent on regional corrections, and some more
    IF(basin(i)%parregion(regiondivision(m_incorr))>0)THEN
      incorr = (1. + regpar(m_incorr,basin(i)%parregion(regiondivision(m_incorr))))     !Correction of inorganic nitrogen
    ELSE
      incorr = 1.
    ENDIF
    IF(basin(i)%parregion(regiondivision(m_oncorr))>0)THEN
      oncorr = (1. + regpar(m_oncorr,basin(i)%parregion(regiondivision(m_oncorr))))     !Correction of organic nitrogen
    ELSE
      oncorr = 1.
    ENDIF
    denpar = (2.-incorr)*lakedatapar(lakedataparindex(i,itype),m_lddenitwl)   !denitrification opposite correction, aka 1-incorrpar=2-(1+incorrpar)
    hsatINpar = genpar(m_hsatINwater)
    prodNpar = lakedatapar(lakedataparindex(i,itype),m_ldwprodn)
    prodPpar = lakedatapar(lakedataparindex(i,itype),m_ldwprodp)
    prodCpar = lakedatapar(lakedataparindex(i,itype),m_ldwprodc)
    prodSipar = lakedatapar(lakedataparindex(i,itype),m_ldwprodsi)
    tmpexpSipar = genpar(m_sitmpexp)
    hsatTPpar = genpar(m_hsatTP)
    limpppar = genpar(m_limsedpp)
    plimsipar = genpar(m_plimsi)
    sedonpar = (2.-oncorr)*lakedatapar(lakedataparindex(i,itype),m_ldsedon)   !sedimentation opposite correction, aka 1-oncorrpar=2-(1+oncorrpar)
    sedpppar = lakedatapar(lakedataparindex(i,itype),m_ldsedpp)
    sedocpar = lakedatapar(lakedataparindex(i,itype),m_ldsedoc)
    sedsspar = lakedatapar(lakedataparindex(i,itype),m_ldsedss)
    sedsipar = lakedatapar(lakedataparindex(i,itype),m_ldsedsi)
              
    !>Calculate the nutrient degradation and uptake processes
    IF(conduct%simP) CALL calculate_lake_tpmean(i,itype,lakestate)
    IF(conduct%simN)THEN
      CALL denitrification_water(i,itype,systemtype,area,denpar,hsatINpar,LAKESTATE=lakestate) !denitrification in lakes
    ENDIF
    CALL production_mineralisation(i,itype,systemtype,area,prodNpar,prodPpar,prodSipar,tmpexpSipar,hsatTPpar,limpppar,plimsipar,LAKESTATE=lakestate)  !primary production and mineralisation in lakes
    IF(conduct%simC) CALL oc_production_mineralisation(systemtype,area,prodCpar,hsatTPpar,limpppar,   &
            lakestate%water(itype,i),lakestate%conc(i_oc,itype,i),  &
            lakestate%temp(itype,i),lakestate%TPmean(itype,i),              &
            lakestate%temp10(itype,i),lakestate%temp20(itype,i))
    CALL macrophyte_uptake(i,itype,systemtype,area,lakedatapar(lakedataparindex(i,itype),m_ldmuptn), &
              lakedatapar(lakedataparindex(i,itype),m_ldmuptp),hsatTPpar,genpar(m_muptdep),limpppar,LAKESTATE=lakestate)
    
    !>Calculate the sedimentation and resuspension processes
    IF(conduct%siltation) CALL sediment_fraction(i,itype,lakestate,ss_frac)
    IF(conduct%simN) CALL lake_sedimentation(i,i_on,itype,sedonpar,genpar(m_limsedon),lakestate,pooladd_tn)
    IF(conduct%simP) CALL lake_sedimentation(i,i_pp,itype,sedpppar,limpppar,lakestate,pooladd_tp)
    IF(conduct%simC) CALL lake_sedimentation(i,i_oc,itype,sedocpar,0.,lakestate,pooladd_oc)
    IF(conduct%simS)THEN
      CALL lake_sedimentation(i,i_ss,itype,sedsspar,genpar(m_limsedss),lakestate,pooladd_ss)
      CALL lake_sedimentation(i,i_ae,itype,genpar(m_sedae),0.,lakestate,pooladd_ae)
      pooladd_ts = (pooladd_ss + pooladd_ae) * area ! #CB: sediment added to lake production pool (SS+AE) (kg) = sediment added to lake production pool (SS+AE) (kg/m2) * lake area (m2)
    ENDIF
    IF(conduct%simSi) CALL lake_sedimentation(i,i_asi,itype,sedsipar,0.,lakestate,pooladd_si)
    CALL lake_siltation_and_flushing(i,itype,ss_frac,lakestate,resuspdown,resusptot)
    resuspdown = resuspdown * area  !resuspended substance to be transported downstream (kg?/ts)
    IF(conduct%simS) poolnet_ts = pooladd_ts - ((resusptot(i_ss) + resusptot(i_ae)) * area) ! #CB: net change in sediment in lake production pool (SS+AE) (kg) = sediment added (kg) - (sediment removed(kg/m2) * lake area (m2))
    IF(conduct%simN) pooladd_tn = pooladd_tn * area !kg
    IF(conduct%simP) pooladd_tp = pooladd_tp * area
    IF(conduct%simC) pooladd_oc = pooladd_oc * area
    IF(conduct%simSi) pooladd_si = pooladd_si * area
    
    IF(conduct%simP) CALL internal_lake_load(i,itype,systemtype,area,lakestate)  !internal load of phosphorus

  END SUBROUTINE substance_processes_in_lake

  !>\brief Calculates the denitrification in river and lakes
  !!Lake processes take place in whole lake volume
  !>
  !>\b Reference ModelDescription Chapter Nitrogen and phosphorus processes in rivers and lakes (Denitrification)
  !-----------------------------------------------------------------------
  SUBROUTINE denitrification_water(i,watertype,systemtype,area,denpar,halfsatINwater,riverstate,lakestate)

    USE MODVAR, ONLY : i_in
    USE HYPEVARIABLES, ONLY : maxdenitriwater   

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<index of current subbasin
    INTEGER, INTENT(IN) :: watertype  !<Lake or river type (1=local, 2=main/outlet)
    INTEGER, INTENT(IN) :: systemtype !<aquatic system type (1=lake, 2=river)
    REAL, INTENT(IN)    :: area       !<lake surface area/river bottom area (m2)
    REAL, INTENT(IN)    :: denpar     !<model parameter denitrification rate (kg/m2/day)
    REAL, INTENT(IN)    :: halfsatINwater  !<model parameter half saturation IN (mg/L)
    TYPE(riverstatetype),INTENT(INOUT),OPTIONAL :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT),OPTIONAL  :: lakestate   !<Lake states

    !Local variables
    REAL, DIMENSION(1) :: denitri_water, inorganicNpool
    REAL tmpfcn, concfcn, watertemp, waterconc,vol
    
    !Local parameters
    INTEGER, PARAMETER :: pooldim = 1

    !Initial pools and values
    IF(systemtype==1) THEN   !lakes
      vol = lakestate%water(watertype,i) * area * 1.0E-6    !0.001 m3
      waterconc = lakestate%conc(i_in,watertype,i)          !mg/L
      inorganicNpool = vol * waterconc                          !kg
      watertemp = lakestate%temp(watertype,i)
    ELSE                     !rivers
      vol = riverstate%water(watertype,i) * 1.0E-3
      waterconc = riverstate%conc(i_in,watertype,i)
      inorganicNpool = vol * waterconc     !kg      
      watertemp = riverstate%temp(watertype,i)
    ENDIF

    !Temperature and concentration dependence factor
    tmpfcn  = tempfactor(watertemp)
    concfcn = halfsatconcfactor(waterconc,halfsatINwater)

    !Denitrification    
    denitri_water = denpar * area * concfcn * tmpfcn   !kg  
    denitri_water = MIN(maxdenitriwater*inorganicNpool, denitri_water)    !max 50% kan be denitrified
    CALL retention_pool(pooldim, inorganicNPool, denitri_water)
    IF(systemtype==1) THEN   !lakes
      CALL new_concentration(inorganicNpool(1),vol,lakestate%conc(i_in,watertype,i))
    ELSE                     !rivers
      IF(riverstate%water(watertype,i) > 0.) THEN
        CALL new_concentration(inorganicNpool(1),vol,riverstate%conc(i_in,watertype,i))
      ENDIF
    ENDIF

  END SUBROUTINE denitrification_water

  !>\brief Calculates transformation between IN/ON and SRP/PP in water.
  !!Also simulated algae production (=ON-production) for sediment simulation.
  !!Simulating the combined processes of primary production and
  !!mineralisation. Lake process in whole lake volume. 
  !>
  !>\b Reference ModelDescription Chapter Nitrogen and phosphorus processes 
  !>in rivers and lakes (Primary production and mineralization)  
  !-----------------------------------------------------------------------
  SUBROUTINE production_mineralisation(i,watertype,systemtype,area, &
                                       prodNpar,prodPpar,prodSipar,tmpexpSipar,halfsatTPwater,  &
                                       limpppar,plimsipar,riverstate,lakestate,depth)

    USE MODVAR, ONLY : conduct, &
                       i_in,i_on,i_sp,i_pp,i_ae,i_asi,i_dsi
    USE HYPEVARIABLES, ONLY : maxprodwater,   &
                              maxdegradwater, &
                              NPratio

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<index of current subbasin
    INTEGER, INTENT(IN) :: watertype  !<Lake or river type (1=local, 2=main/outlet)
    INTEGER, INTENT(IN) :: systemtype !<aquatic system type (1=lake, 2=river)
    REAL, INTENT(IN)    :: area       !<lake surface area/ river bottom area (m2)
    REAL, INTENT(IN)    :: prodNpar   !<model parameter production rate ON in water
    REAL, INTENT(IN)    :: prodPpar   !<model parameter production rate PP in water
    REAL, INTENT(IN)    :: prodSipar  !<model parameter production Algal Si
    REAL, INTENT(IN)    :: tmpexpSipar   !<model parameter temperature dependence production Algal Si
    REAL, INTENT(IN)    :: halfsatTPwater !<model parameter half saturation TP (mg/L)
    REAL, INTENT(IN)    :: limpppar   !<limitation of sedimentation parameter (mg/L)
    REAL, INTENT(IN)    :: plimsipar   !<P limitation of Si production parameter (mg/L)
    TYPE(riverstatetype),INTENT(INOUT),OPTIONAL :: riverstate !<River states
    TYPE(lakestatetype),INTENT(INOUT),OPTIONAL  :: lakestate  !<Lake states
    REAL, INTENT(IN), OPTIONAL :: depth      !<river depth (m) 

    !Local variables
    REAL, DIMENSION(1) :: ONpool, INpool, SRPpool,PPpool,AEpool,DSipool,ASipool
    REAL, DIMENSION(1) :: minprodN, minprodP,minprodAE,minprodSi
    REAL watertemp,waterTPmean,temp10,temp20
    REAL tmpfcn, tmpfcn1, tmpfcn2, TPfcn, actpconc, TPfcnSi
    REAL vol
    REAL waterdepth     !m

    !Local parameters
    INTEGER, PARAMETER :: pooldim = 1

    IF(.NOT.(conduct%simN.OR.conduct%simP.OR.conduct%simS.OR.conduct%simSi)) RETURN

    !Pools of nutrients in the water, water temperature and fraction of depth of water volume that is active
    IF (systemtype==1)THEN   !lakes
      vol = lakestate%water(watertype,i) * area / 1.0E6
      IF(conduct%simN)THEN 
        INpool = vol * lakestate%conc(i_in,watertype,i)   !kg
        ONpool = vol * lakestate%conc(i_on,watertype,i)   !kg
      ENDIF
      IF(conduct%simP)THEN
        SRPpool = vol * lakestate%conc(i_sp,watertype,i) !kg
        PPpool  = vol * lakestate%conc(i_pp,watertype,i)  !kg
      ENDIF
      IF(conduct%simS) AEpool = vol * lakestate%conc(i_ae,watertype,i) !kg N
      IF(conduct%simSi)THEN
        DSipool = vol * lakestate%conc(i_dsi,watertype,i) !kg
        ASipool = vol * lakestate%conc(i_asi,watertype,i) !kg
      ENDIF      
   ELSE                     !rivers
      vol = riverstate%water(watertype,i) / 1.0E3
      IF(conduct%simN)THEN
        INpool = vol * riverstate%conc(i_in,watertype,i) !kg
        ONpool = vol * riverstate%conc(i_on,watertype,i) !kg
      ENDIF
      IF(conduct%simP)THEN
        SRPpool = vol * riverstate%conc(i_sp,watertype,i)    !kg    
        PPpool  = vol * riverstate%conc(i_pp,watertype,i)    !kg
      ENDIF
      IF(conduct%simS) AEpool = vol * riverstate%conc(i_ae,watertype,i)    !kg N 
      IF(conduct%simSi)THEN
        DSipool = vol * riverstate%conc(i_dsi,watertype,i)    !kg    
        ASipool = vol * riverstate%conc(i_asi,watertype,i)    !kg
      ENDIF
    ENDIF

    !Set help variables
    IF (systemtype==1) THEN   !lakes
      watertemp = lakestate%temp(watertype,i)
      waterdepth = lakestate%water(watertype,i) * 1.E-3
      waterTPmean = lakestate%TPmean(watertype,i)
      temp10 = lakestate%temp10(watertype,i)
      temp20 = lakestate%temp20(watertype,i)
    ELSE                     !rivers
      watertemp = riverstate%temp(watertype,i)  
      waterdepth = depth
      waterTPmean = riverstate%TPmean(watertype,i)
      temp10 = riverstate%temp10(watertype,i)
      temp20 = riverstate%temp20(watertype,i)
    ENDIF

    IF(watertemp > 0.)THEN
      !Total phosphorus concentration dependent factor; one general and one for silica
      TPfcn = halfsatconcfactor(MAX(waterTPmean-limpppar,0.),halfsatTPwater)
      IF(i_dsi>0)THEN
        IF(i_sp>0)THEN
          IF (systemtype==1) THEN 
            actpconc = MAX(0.,lakestate%conc(i_sp,watertype,i) + lakestate%conc(i_pp,watertype,i) - plimsipar)
          ELSE
            actpconc = MAX(0.,riverstate%conc(i_sp,watertype,i) + riverstate%conc(i_pp,watertype,i) - plimsipar)
          ENDIF
        ELSE
          actpconc = MAX(waterTPmean-limpppar,0.)
        ENDIF
        TPfcnSi = halfsatconcfactor(actpconc,halfsatTPwater)
      ENDIF

      !Temperature dependent factor
      tmpfcn1 = watertemp / 20.    
      tmpfcn2 = (temp10 - temp20) / 5.
      tmpfcn = tmpfcn1*tmpfcn2

      !Production/mineralisation
      IF(conduct%simN.OR.conduct%simS)THEN
        minprodN = prodNpar * TPfcn * tmpfcn * waterdepth * area  !kg
        IF(conduct%simN)THEN
          IF(minprodN(1) > 0.) THEN  !production        
            minprodN = MIN(maxprodwater * INpool, minprodN)
          ELSE                       !mineralisation
            minprodN = MAX(-maxdegradwater * ONpool, minprodN)
          ENDIF
          CALL retention_pool(pooldim,INpool,minprodN)   !minprodN may be negative
          CALL production_pool(pooldim,ONpool,minprodN)
        ENDIF
        IF(conduct%simS)THEN
          minprodAE = MAX(-AEpool, minprodN)  !AE can be produced from zero, thus no max depending on pool.
          CALL production_pool(pooldim,AEpool,minprodAE) !AEpool is a part of ON pool
        ENDIF
      ENDIF
      IF(conduct%simP)THEN
        minprodP = prodPpar * NPratio * TPfcn * tmpfcn * waterdepth * area  !kg  
        IF(minprodP(1) > 0.) THEN  !production        
          minprodP = MIN(maxprodwater * SRPpool,minprodP)
        ELSE                       !mineralisation
          minprodP = MAX(-maxdegradwater * PPpool,minprodP)
        ENDIF
        CALL retention_pool(pooldim,SRPpool,minprodP)    !minprodP may be negative
        CALL production_pool(pooldim,PPpool,minprodP)
      ENDIF
      IF(conduct%SimSi)THEN
        minprodSi = prodSipar * TPfcnSi * tmpfcn1**tmpexpSipar *tmpfcn2 * waterdepth * area  !kg
        IF(minprodSi(1) > 0.) THEN  !production        
          minprodSi = MIN(maxprodwater * DSipool, minprodSi)
        ELSE                       !mineralisation
          minprodSi = MAX(-maxdegradwater * ASipool, minprodSi)
        ENDIF
        CALL retention_pool(pooldim,DSipool,minprodSi)   !minprodSi may be negative
        CALL production_pool(pooldim,ASipool,minprodSi)
      ENDIF

      !New concentration due to changes in pools
      IF(systemtype==1) THEN            !lakes
        IF(conduct%simN) CALL new_concentration(INpool(1),vol,lakestate%conc(i_in,watertype,i))
        IF(conduct%simN) CALL new_concentration(ONpool(1),vol,lakestate%conc(i_on,watertype,i))
        IF(conduct%simP) CALL new_concentration(SRPpool(1),vol,lakestate%conc(i_sp,watertype,i))
        IF(conduct%simP) CALL new_concentration(PPpool(1),vol,lakestate%conc(i_pp,watertype,i))
        IF(conduct%simS) CALL new_concentration(AEpool(1),vol,lakestate%conc(i_ae,watertype,i))
        IF(conduct%simSi) CALL new_concentration(DSipool(1),vol,lakestate%conc(i_dsi,watertype,i))
        IF(conduct%simSi) CALL new_concentration(ASipool(1),vol,lakestate%conc(i_asi,watertype,i))       
      ELSE                              !rivers
        IF(riverstate%water(watertype,i) > 0.) THEN
          IF(conduct%simN) CALL new_concentration(INpool(1),vol,riverstate%conc(i_in,watertype,i))
          IF(conduct%simN) CALL new_concentration(ONpool(1),vol,riverstate%conc(i_on,watertype,i))
          IF(conduct%simP) CALL new_concentration(SRPpool(1),vol,riverstate%conc(i_sp,watertype,i))
          IF(conduct%simP) CALL new_concentration(PPpool(1),vol,riverstate%conc(i_pp,watertype,i))
          IF(conduct%simS) CALL new_concentration(AEpool(1),vol,riverstate%conc(i_ae,watertype,i))          
          IF(conduct%simSi) CALL new_concentration(DSipool(1),vol,riverstate%conc(i_dsi,watertype,i))
          IF(conduct%simSi) CALL new_concentration(ASipool(1),vol,riverstate%conc(i_asi,watertype,i))
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE production_mineralisation

  !>\brief Calculates macrophyte uptake of IN/SP in surface water.
  !>
  !>\b Reference ModelDescription Chapter Nitrogen and phosphorus processes 
  !>in rivers and lakes (Macrophyte uptake)  
  !-----------------------------------------------------------------------
  SUBROUTINE macrophyte_uptake(i,watertype,systemtype,area, &
                               muptNpar,muptPpar,halfsatTPwater,  &
                               proddeppar,limpppar,lakestate,riverstate)

    USE MODVAR, ONLY : conduct, &
                       i_in,i_sp

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<index of current subbasin
    INTEGER, INTENT(IN) :: watertype  !<Lake or river type (1=local, 2=main/outlet)
    INTEGER, INTENT(IN) :: systemtype !<aquatic system type (1=lake, 2=river)
    REAL, INTENT(IN)    :: area       !<lake surface area/ river bottom area (m2)
    REAL, INTENT(IN)    :: muptNpar   !<model parameter production rate ON in water
    REAL, INTENT(IN)    :: muptPpar   !<model parameter production rate PP in water
    REAL, INTENT(IN)    :: halfsatTPwater !<model parameter half saturation TP (mg/L)
    REAL, INTENT(IN)    :: proddeppar   !<depth above which macrophyte grows (m)
    REAL, INTENT(IN)    :: limpppar    !<(mg/L)
    TYPE(lakestatetype),INTENT(INOUT),OPTIONAL  :: lakestate  !<Lake states
    TYPE(riverstatetype),INTENT(INOUT),OPTIONAL :: riverstate !<River states

    !Local variables
    REAL, DIMENSION(1) :: INpool, SRPpool
    REAL, DIMENSION(1) :: macrouptake_in, macrouptake_sp
    REAL watertemp,waterTPmean,temp20
    REAL tmpfcn, tmpfcn1, tmpfcn2, TPfcn
    REAL vol
    REAL fracarea

    !Local parameters
    INTEGER, PARAMETER :: pooldim = 1

    IF(.NOT.(conduct%simN.OR.conduct%simP)) RETURN
    IF((muptNpar==0. .AND. muptPpar==0.) .OR. proddeppar==0.)RETURN
    
    !Pools of nutrients in the water, water temperature and fraction of depth of water volume that is active
    IF(systemtype==1)THEN   !lakes
      vol = lakestate%water(watertype,i) * area / 1.0E6   !km*m2
      IF(conduct%simN)THEN 
        INpool = vol * lakestate%conc(i_in,watertype,i)   !kg
      ENDIF
      IF(conduct%simP)THEN
        SRPpool = vol * lakestate%conc(i_sp,watertype,i) !kg
      ENDIF
    ELSEIF(systemtype==2)THEN   !rivers
      vol = riverstate%water(watertype,i) / 1.0E3   !km*m2
      IF(conduct%simN)THEN
        INpool = vol * riverstate%conc(i_in,watertype,i) !kg
      ENDIF
      IF(conduct%simP)THEN
        SRPpool = vol * riverstate%conc(i_sp,watertype,i)    !kg    
      ENDIF
    ENDIF

    !Set help variables
    IF(systemtype==1) THEN   !lakes
      watertemp = lakestate%temp(watertype,i)
      waterTPmean = lakestate%TPmean(watertype,i)
      temp20 = lakestate%temp20(watertype,i)
    ELSE                     !rivers
      watertemp = riverstate%temp(watertype,i)  
      waterTPmean = riverstate%TPmean(watertype,i)
      temp20 = riverstate%temp20(watertype,i)
    ENDIF
    !Temperature dependent factor
    tmpfcn1 = (MAX(0.,watertemp) / 20.)**0.3  !genpar(m_wltmpexp) (maybe parameter?)
    tmpfcn2 = (watertemp - temp20) / 5.
    tmpfcn = MAX(0.,tmpfcn1*tmpfcn2)

    IF(watertemp > 0.)THEN
      !Total phosphorus concentration dependent factor
      TPfcn = halfsatconcfactor(MAX(waterTPmean-limpppar,0.),halfsatTPwater)
     
      !Estimated fraction of lake area above production depth
      fracarea = MIN(1.0, proddeppar * (0.001/((vol/area)*2)))  !m*/(km*m2/m2)*10^-3 -> unitless
      
      !Macrophyte uptake
      IF(conduct%simN)THEN
        macrouptake_in = muptNpar * tmpfcn * fracarea * area * TPfcn
        macrouptake_in = MIN(0.5 * INpool, macrouptake_in)
        CALL retention_pool(pooldim,INpool,macrouptake_in)
      ENDIF
      IF(conduct%simP)THEN
        macrouptake_sp = muptPpar * tmpfcn * fracarea * area * TPfcn
        macrouptake_sp = MIN(0.5 * SRPpool, macrouptake_sp)
        CALL retention_pool(pooldim,SRPpool,macrouptake_sp)   
      ENDIF

      !New concentration due to changes in pools
      IF(systemtype==1) THEN            !lakes
        IF(conduct%simN) CALL new_concentration(INpool(1),vol,lakestate%conc(i_in,watertype,i))
        IF(conduct%simP) CALL new_concentration(SRPpool(1),vol,lakestate%conc(i_sp,watertype,i))
      ELSE                              !rivers
        IF(riverstate%water(watertype,i) > 0.) THEN
          IF(conduct%simN) CALL new_concentration(INpool(1),vol,riverstate%conc(i_in,watertype,i))
          IF(conduct%simP) CALL new_concentration(SRPpool(1),vol,riverstate%conc(i_sp,watertype,i))
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE macrophyte_uptake

  !>\brief Calculate sedimentation of substance in lake
  !>
  !>\b Reference ModelDescription Chapter Nitrogen and phosphorus processes in 
  !>rivers and lakes (Sedimentation/Resuspension) and Organic carbon (River and lakes - Sedimentation)
  !--------------------------------------------------------------------------
  SUBROUTINE lake_sedimentation(i,substance,watertype,sedrate,limsedpar,lakestate,pooladd)

    USE MODVAR, ONLY : conduct, i_ss

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<current index of subbasin
    INTEGER, INTENT(IN) :: substance  !<current index of substance (PP,ON,SS,AE)
    INTEGER, INTENT(IN) :: watertype  !<Lake type (1=local, 2=outlet)
    REAL, INTENT(IN)    :: sedrate    !<sedimentation rate  (lakes) (m/ts)
    REAL, INTENT(IN)    :: limsedpar  !<concentration limit for sedimentation (mg/L)
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    REAL, INTENT(OUT), OPTIONAL :: pooladd  ! #CB: amount of substance added to pool (kg/m2/ts)
    
    !Local variables
    REAL, DIMENSION(1) :: pool     !substance pool in water per area (g/m2) ((kg))
    REAL, DIMENSION(1) :: delta    !sedimentation (g/m2/ts)
    REAL, DIMENSION(1), PARAMETER :: onetimestep = 1.    !age increase during a timestep (timesteps)
    REAL depth                     !lake average depth (m)

    !>\b Algoritm
    !>Calculate water volume, amount of substance in water (pool) and the sedimentation
    depth = lakestate%water(watertype,i)  * 1.0E-3        !m
    pool = depth * lakestate%conc(substance,watertype,i) !mg/l*m
    delta = sedrate * MAX(lakestate%conc(substance,watertype,i)-limsedpar,0.) !m/ts*mg/l

    !>Remove sedimentation from the water pool
    CALL retention_pool(1,pool,delta)
    
    !>Amount of substance added to sediment pool
    IF(PRESENT(pooladd)) pooladd = delta(1) / 1.0E3 ! convert from g/m2/ts to kg/m2/ts

    !>Add sedimentation to the sediment pool
    IF(conduct%siltation)THEN
      CALL production_pool(1,lakestate%sedpool(substance,watertype,i),delta / 1.0E3)
      IF(substance==i_ss) CALL production_pool(1,lakestate%sedage(watertype,i),onetimestep) !Update sediment age only once
    ENDIF

    !>Calculate the new concentration in the water due to the change in the water pool
    CALL new_concentration(pool(1),depth,lakestate%conc(substance,watertype,i))

  END SUBROUTINE lake_sedimentation

  !>\brief Calculate density of sediment in reservoir sediment pool
  !>
  !> Calculate the density of substance in sediment (SS and AE) with different options:
  !> Siltation Option 1: No Compaction + General Density - Use density (kg/m3) of substance (SS or AE) from general parameter in par.txt
  !> Siltation Option 2: No Compaction + Density from Soil Fraction & Reservoir Operation - Calculate the density (kg/m3) of substance (SS or AE) from GeoData.txt soil fractions and reservoir operation mode
  !> Siltation Option 3: Compaction + Density from Soil Fraction & Reservoir Operation - Calculate the compacted density (kg/m3) of substance (SS or AE) from GeoData.txt soil fractions and reservoir operation mode
  !>
  !>\b Reference ModelDescription Chapter Sediment
  !--------------------------------------------------------------------------
  SUBROUTINE sediment_pool_density(watertype,res_mode,siltation_option,age,frac_clay,frac_silt,frac_sand,density)
  
  USE MODVAR, ONLY : genpar, timesteps_per_day
  USE HYPEVARIABLES, ONLY : m_lseddens
  
    !Argument declarations
    INTEGER, INTENT(IN) :: watertype    !<Lake type (1=local, 2=outlet)
    INTEGER, INTENT(IN) :: res_mode     !<reservoir operation mode (1=sediment always submerged or nearly submerged, 2=normally moderate to considerable reservoir drawdown, 3=reservoir normally empty, 4=riverbed sediments)
    INTEGER, INTENT(IN) :: siltation_option     !<siltation density option (1=constant, 2=depentent on soil, 3=compactation depending on soil and age
    REAL, INTENT(IN)    :: age          !<number of timesteps sediment has been compacting
    REAL, INTENT(IN)    :: frac_clay    !<fraction of clay in incoming sediment
    REAL, INTENT(IN)    :: frac_silt    !<fraction of silt in incoming sediment
    REAL, INTENT(IN)    :: frac_sand    !<fraction of sand in incoming sediment
    REAL, INTENT(OUT)   :: density      !<new sediment production pool density (kg/m3)
    
    !Local variables
    INTEGER mode  !reservoir operation mode
    REAL agedays  !age in days
    REAL d_clay !density clay, initial (kg/m3)
    REAL d_silt !density silt, initial (kg/m3)
    REAL d_sand !density sand, initial (kg/m3)
    REAL k_clay !compaction coefficient, clay
    REAL k_silt !compaction coefficient, silt
    
    ! Siltation Option 1: No Compaction + General Density - Use density (kg/m3) of substance (SS or AE) from general parameter in par.txt
    IF(siltation_option==1) THEN
      density = genpar(m_lseddens)
    ELSE
      
    ! Siltation Option 2 and 3: Density from Soil Fraction & Reservoir Operation - Calculate the density (kg/m3) of substance (SS or AE) from GeoData.txt soil fractions and reservoir operation mode
      
      !If lake is an olake (watertype==2), then use reservoir operation mode from LakeData; If lake is an ilake, then use riverbed sediments operation mode
      IF(watertype==2 .AND. (res_mode>0 .AND. res_mode<5)) THEN
        mode = res_mode
      ELSE
        mode = 4
      ENDIF
        
      !Get densities and compaction coefficients for clay/silt/sand based on reservoir operation mode
      IF(mode==1) THEN
        d_clay = 416.
        d_silt = 1120.
        d_sand = 1550.
        k_clay = 256.
        k_silt = 91.
      ELSEIF(mode==2) THEN
        d_clay = 561.
        d_silt = 1140.
        d_sand = 1550.
        k_clay = 135.
        k_silt = 29.
      ELSEIF(mode==3) THEN
        d_clay = 641.
        d_silt = 1150.
        d_sand = 1550.
        k_clay = 0.     !No compactation
        k_silt = 0.
      ELSEIF(mode==4) THEN
        d_clay = 961.
        d_silt = 1170.
        d_sand = 1550.
        k_clay = 0.     !No compatation
        k_silt = 0.
      ENDIF
    
      !Calculate initial density of all deposits
      density = ((frac_clay * d_clay) + (frac_silt * d_silt) + (frac_sand * d_sand))
      agedays = age/timesteps_per_day
      
      !Calculate average sediment bulk density of all deposits after compaction (kg/m3) during years of operation
      IF(agedays>365. .AND. siltation_option==3) THEN
        density = density + (0.4343 * ((frac_clay * k_clay) + (frac_silt * k_silt)) * ((((agedays / 365.) / ((agedays / 365.) - 1.)) * LOG(agedays / 365.)) - 1.)) ! Convert age from days to years by dividing  by 365
      ENDIF
    ENDIF
    
  END SUBROUTINE sediment_pool_density

  !>\brief Calculate siltation of lake with sespended sediments (and algae) and flushing
  !>
  !>\b Reference ModelDescription Chapter Sediment
  !--------------------------------------------------------------------------
  SUBROUTINE lake_siltation_and_flushing(i,watertype,ss_frac,lakestate,resuspdown,resusptot)

    USE MODVAR, ONLY : basin,conduct,conductwarning,current_time, &
                       lake,lakeindex,numsubstances, &
                       elake,lakebasin,lakebasinindex, &
                       modeloption, p_siltation, i_ss, i_ae, &
                       timesteps_per_day
    USE HYPEVARIABLES, ONLY : thresholddepth

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<current index of subbasin
    INTEGER, INTENT(IN) :: watertype  !<Lake type (1=local, 2=outlet)
    REAL, INTENT(IN)    :: ss_frac    !<ss fraction of sediment pool depth
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    REAL, INTENT(OUT)   :: resuspdown(numsubstances)  !<resuspended/flushed substance to be transported downstream (kg/m2?/ts)
    REAL, INTENT(OUT)   :: resusptot(numsubstances) ! #CB: amount of sediment removed from lake sediment pool kg/m2
    
    !Local variables
    INTEGER current_elake
    INTEGER flushmode,resopmode    !current lake data
    INTEGER destpar                !current lake data
    REAL density                   !sediment density
    REAL seddepth(numsubstances)   !sediment depth (m)
    REAL flush_rate(numsubstances) !flush rate (kg/m2/day)
    REAL noflushperiod             !current lake data, days before flushing may start: (frequency-1) + flushperiod
    REAL ratepar,startpar,stoppar   !current lake data
    LOGICAL reservoir_flushstatus   !status of reservoir, flushing possible
    LOGICAL lake_in_lakedata        !status of lake information
    LOGICAL elake_in_lakedata       !status of lake information, multi-basin lake
    LOGICAL last_lakebasin

      !>\b Algoritm
      !>Initiate resuspension and sediment flush rates to 0 - set to 0 for non-sediment substances
      resuspdown = 0. ! Set substance flushing to zero
      resusptot = 0.  ! Set substance levaing sediment to zero
      IF(.NOT.conduct%siltation) RETURN
      IF(.NOT.conduct%simS) RETURN
      flush_rate = 0.
        
      !>Check if lake information is available and if flushing is possible for this lake/reservoir (set help variables)
      lake_in_lakedata = .FALSE.
      elake_in_lakedata = .FALSE.
      last_lakebasin = .FALSE.
      reservoir_flushstatus = .FALSE.
      IF(watertype == 2)THEN !No Lakedata regulation of ilakes
        IF(ALLOCATED(lakeindex))THEN
          IF(lakeindex(i)>0)THEN    !No flushing for lakes not in LakeData
            lake_in_lakedata = .TRUE.
          ENDIF
        ENDIF
        IF(ALLOCATED(lakebasinindex))THEN
          IF(lakebasinindex(i)>0)THEN
            elake_in_lakedata = .TRUE.
            current_elake = lakebasin(lakebasinindex(i))%ilk
            last_lakebasin = lakebasin(lakebasinindex(i))%last
          ENDIF
        ENDIF
      ENDIF
      IF(watertype == 2)THEN !No flushing of ilakes
        IF(lake_in_lakedata)THEN !Flushing of olake in LakeData?
          IF(lake(lakeindex(i))%resmode .NE. 4)THEN    !No flushing for reservoir that are always filled with water
            IF(lake(lakeindex(i))%sedmgmt%mode .NE. 0) THEN  !Zero is no flushing mode
              reservoir_flushstatus = .TRUE.
            ENDIF
          ENDIF
        ELSEIF(elake_in_lakedata.AND.last_lakebasin)THEN !Flushing of last lakebasin of multi-basin lake?
          IF(elake(current_elake)%resmode .NE. 4)THEN    !No flushing for reservoir that are always filled with water
            IF(elake(current_elake)%sedmgmt%mode .NE. 0) THEN  !Zero is no flushing mode
              reservoir_flushstatus = .TRUE.
            ENDIF
          ENDIF
        ENDIF
      ENDIF
        
        
      !>\li Sediment is accumulating
      !>Calculate the density of substance (SS or AE) and update lake depth
      IF(lakestate%sedflush(watertype,i) == 0) THEN
        IF(lake_in_lakedata)THEN
          CALL sediment_pool_density(watertype,lake(lakeindex(i))%resmode,modeloption(p_siltation),lakestate%sedage(watertype,i),basin(i)%clay,basin(i)%silt,basin(i)%sand,density)
        ELSEIF(elake_in_lakedata)THEN
          CALL sediment_pool_density(watertype,elake(current_elake)%resmode,modeloption(p_siltation),lakestate%sedage(watertype,i),basin(i)%clay,basin(i)%silt,basin(i)%sand,density)
        ELSE  !if operation mode not defined use riverbed sediments 
          CALL sediment_pool_density(watertype,4,modeloption(p_siltation),lakestate%sedage(watertype,i),basin(i)%clay,basin(i)%silt,basin(i)%sand,density)
        ENDIF
        lakestate%seddens(watertype,i) = density  !update state for no flush time steps
              
        ! Calculate new depth of substance
        seddepth = lakestate%sedpool(:,watertype,i) / density   !
    
        ! Update lake depth - set to zero if sediment has completely filled lake
        thresholddepth(watertype,i) = basin(i)%lakedepth(watertype) - seddepth(i_ss) - seddepth(i_ae)
        IF(thresholddepth(watertype,i)<0.)THEN
          !WRITE(6,*) 'WARNING: Lake has filled up by siltation. Setting lake depth'
          !WRITE(6,*) 'WARNING: to zero for subbasin ',basin(i)%subid
          !Take the removed sediment and add to resuspension (as Conrad did).
          resuspdown(i_ss) = - thresholddepth(watertype,i)*density
          resuspdown(i_ae) = resuspdown(i_ss) * seddepth(i_ae)/(seddepth(i_ss)+seddepth(i_ae))
          resuspdown(i_ss) = resuspdown(i_ss) - resuspdown(i_ae)
          thresholddepth(watertype,i) = 0.
          lakestate%sedpool(i_ss,watertype,i) = lakestate%sedpool(i_ss,watertype,i) - resuspdown(i_ss)
          lakestate%sedpool(i_ae,watertype,i) = lakestate%sedpool(i_ae,watertype,i) - resuspdown(i_ae)
        ENDIF
      ENDIF
            
      !>\li Sediment is flushing
      !>Calculate the flushing of substance (SS or AE) and update lake depth
      IF(lakestate%sedflush(watertype,i) == 1) THEN
        !Get lake data parameter values
        IF(lake_in_lakedata)THEN
          flushmode = lake(lakeindex(i))%sedmgmt%mode
          resopmode = lake(lakeindex(i))%resmode
          ratepar = lake(lakeindex(i))%sedmgmt%rate
          startpar = lake(lakeindex(i))%sedmgmt%start
          stoppar = lake(lakeindex(i))%sedmgmt%sstop
          destpar = lake(lakeindex(i))%sedmgmt%dest
        ELSEIF(elake_in_lakedata)THEN
          flushmode = elake(current_elake)%sedmgmt%mode
          resopmode = elake(current_elake)%resmode
          ratepar = elake(current_elake)%sedmgmt%rate
          startpar = elake(current_elake)%sedmgmt%start
          stoppar = elake(current_elake)%sedmgmt%sstop
          destpar = elake(current_elake)%sedmgmt%dest
        ENDIF

        ! Calculate SS and AE flush rate (kg/m2/ts) = SS density * SS flush rate
        flush_rate(i_ss) = lakestate%seddens(watertype,i) * ratepar * basin(i)%lakedepth(watertype) * ss_frac / timesteps_per_day
        flush_rate(i_ae) = lakestate%seddens(watertype,i) * ratepar * basin(i)%lakedepth(watertype) * (1. - ss_frac) / timesteps_per_day 
                  
        ! Remove SS and AE from production pool (kg/m2) - set to 0 if sediment has been completely flushed
        CALL retention_pool(numsubstances,lakestate%sedpool(:,watertype,i),flush_rate)  !flush_rate may change
                
        ! Add sediment flush rate to resuspension downstream if lake is an olake, reservoir is not completely full of sediment, and flushing occurs (reservoir operation mode isn't river sediments, flushing mode is not set to 0 (Don't Flush), and flushing destination is set to 1 (flow downstream))
        IF(watertype == 2 .AND. (resopmode .NE. 4) .AND. (flushmode .NE. 0) .AND. (destpar == 1)) THEN
          resuspdown = resuspdown + flush_rate ! kg/m2/timestep
        ENDIF
        ! Export amount of substance removed from the pool (kg/m2)
        resusptot = flush_rate
                    
        ! Calculate new depth of substance - Assume constant sediment density while flushing
        seddepth = lakestate%sedpool(:,watertype,i) / lakestate%seddens(watertype,i)

        ! Update lake depth - set to initial lake depth if sediment has been completely flushed
        thresholddepth(watertype,i) = basin(i)%lakedepth(watertype) - seddepth(i_ss) - seddepth(i_ae)
        IF(thresholddepth(watertype,i)>basin(i)%lakedepth(watertype))THEN
          IF(conductwarning) WRITE(6,*) 'WARNING: Lake has emptied too much during flushing. Setting depth to initial.'
          thresholddepth(watertype,i) = basin(i)%lakedepth(watertype)
        ENDIF
        IF(thresholddepth(watertype,i)<0.)THEN
          WRITE(6,*) 'Lake has filled up during flushing!?. Setting depth to zero.'
          !Take the removed sediment and add to resuspension as Conrad did (above).
          resuspdown(i_ss) = - thresholddepth(watertype,i)*lakestate%seddens(watertype,i)
          resuspdown(i_ae) = resuspdown(i_ss) * seddepth(i_ae)/(seddepth(i_ss)+seddepth(i_ae))
          resuspdown(i_ss) = resuspdown(i_ss) - resuspdown(i_ae)
          thresholddepth(watertype,i) = 0.
          lakestate%sedpool(i_ss,watertype,i) = lakestate%sedpool(i_ss,watertype,i) - resuspdown(i_ss)
          lakestate%sedpool(i_ae,watertype,i) = lakestate%sedpool(i_ae,watertype,i) - resuspdown(i_ae)
        ENDIF

      ENDIF
        
      !>Determine whether or not to flush sediments next timestep and reset sediment age
      IF(reservoir_flushstatus) THEN
        IF(lake_in_lakedata)THEN
          ! Flush Mode 1: Flush according to reservoir capacity
          IF(lake(lakeindex(i))%sedmgmt%mode == 1) THEN
            ! Sediments have accumulated past threshold, so start flushing
            IF((lakestate%sedflush(watertype,i) == 0) .AND. (thresholddepth(watertype,i) <= (lake(lakeindex(i))%sedmgmt%start * basin(i)%lakedepth(watertype)))) THEN
              lakestate%sedflush(watertype,i) = 1 ! Change state to indicate flushing
            ! Sediments have been flushed past threshold, so stop flushing
            ELSEIF((lakestate%sedflush(watertype,i) == 1) .AND. (((lakestate%sedpool(i_ss,2,i) == 0) .AND. (lakestate%sedpool(i_ae,2,i) == 0)) .OR. (thresholddepth(watertype,i) >= (lake(lakeindex(i))%sedmgmt%sstop * basin(i)%lakedepth(watertype))))) THEN ! SS and AE pool are both 0 (If flushing to initial depth) or flush past threshold
              lakestate%sedflush(watertype,i) = 0 ! Change state to indicate not flushing
              lakestate%sedage(watertype,i) = 1 ! Reset sediment pool age
            ENDIF

          ! Flush Mode 2: Flush according to day of year
          ELSEIF(lake(lakeindex(i))%sedmgmt%mode == 2) THEN
            noflushperiod = lake(lakeindex(i))%sedmgmt%sstop - lake(lakeindex(i))%sedmgmt%start
            IF(noflushperiod < 0.) noflushperiod = noflushperiod + 365.
            noflushperiod = noflushperiod + (lake(lakeindex(i))%sedmgmt%frequency-1)*365.
            ! Time of flushing has come, so start flushing
            IF((lakestate%sedflush(watertype,i) == 0) .AND. (current_time%dayno==lake(lakeindex(i))%sedmgmt%start) .AND. (lakestate%sedage(watertype,i)/timesteps_per_day > noflushperiod)) THEN
              lakestate%sedflush(watertype,i) = 1 ! Change state to indicate flushing
            ! Sediments have been flushed past threshold, so stop flushing
            ELSEIF((lakestate%sedflush(watertype,i) == 1) .AND. (current_time%dayno==lake(lakeindex(i))%sedmgmt%sstop)) THEN
              lakestate%sedflush(watertype,i) = 0 ! Change state to indicate not flushing
              lakestate%sedage(watertype,i) = 1 ! Reset sediment pool age
            ENDIF
          ENDIF
            
        ! Only last lakebasin may be flushed in a multibasin lake
        ELSEIF(elake_in_lakedata)THEN
          ! Flush Mode 1: Flush according to reservoir capacity
          IF(elake(current_elake)%sedmgmt%mode == 1) THEN
            ! Sediments have accumulated past threshold, so start flushing
            IF((lakestate%sedflush(watertype,i) == 0) .AND. (thresholddepth(watertype,i) <= (elake(current_elake)%sedmgmt%start * basin(i)%lakedepth(watertype)))) THEN
              lakestate%sedflush(watertype,i) = 1 ! Change state to indicate flushing
            ! Sediments have been flushed past threshold, so stop flushing
            ELSEIF((lakestate%sedflush(watertype,i) == 1) .AND. (((lakestate%sedpool(i_ss,2,i) == 0) .AND. (lakestate%sedpool(i_ae,2,i) == 0)) .OR. (thresholddepth(watertype,i) >= (elake(current_elake)%sedmgmt%sstop * basin(i)%lakedepth(watertype))))) THEN ! SS and AE pool are both 0 (If flushing to initial depth) or flush past threshold
              lakestate%sedflush(watertype,i) = 0 ! Change state to indicate not flushing
              lakestate%sedage(watertype,i) = 1 ! Reset sediemnt pool age
            ENDIF

          ! Flush Mode 2: Flush according to day of year
          ELSEIF(elake(current_elake)%sedmgmt%mode == 2) THEN
            noflushperiod = elake(current_elake)%sedmgmt%sstop - elake(current_elake)%sedmgmt%start
            IF(noflushperiod < 0.) noflushperiod = noflushperiod + 365.
            noflushperiod = noflushperiod + (elake(current_elake)%sedmgmt%frequency-1)*365.
            ! Time of flushing has come, so start flushing
            IF((lakestate%sedflush(watertype,i) == 0) .AND. (current_time%dayno==elake(current_elake)%sedmgmt%start) .AND. (lakestate%sedage(watertype,i)/timesteps_per_day > noflushperiod)) THEN
              lakestate%sedflush(watertype,i) = 1 ! Change state to indicate flushing
            ! Sediments have been flushed past threshold, so stop flushing
            ELSEIF((lakestate%sedflush(watertype,i) == 1) .AND. (current_time%dayno==elake(current_elake)%sedmgmt%sstop)) THEN
              lakestate%sedflush(watertype,i) = 0 ! Change state to indicate not flushing
              lakestate%sedage(watertype,i) = 1 ! Reset sediment pool age
            ENDIF
          ENDIF
        ENDIF !lake_in_lakedata/elake_in_lakedata
      ENDIF !reservoir_flushstatus
    
  END SUBROUTINE lake_siltation_and_flushing
  
  !>\brief Calculate fraction of siltation depth that is due to SS
  !>
  !>\b Reference ModelDescription Chapter Sediment
  !--------------------------------------------------------------------------
  SUBROUTINE sediment_fraction(i,watertype,lakestate,ss_frac)

    USE MODVAR, ONLY : conduct, i_ss, i_ae

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<current index of subbasin
    INTEGER, INTENT(IN) :: watertype  !<Lake type (1=local, 2=outlet)
    TYPE(lakestatetype),INTENT(IN) :: lakestate  !<Lake state
    REAL, INTENT(OUT)   :: ss_frac  !<fraction of resuspended substance that is SS (-)
    
      !>\b Algoritm
      ss_frac = 0.    !Default, set SS substance fraction of flushing to zero
      IF(.NOT.conduct%siltation) RETURN
      IF(.NOT.conduct%simS) RETURN
      IF(watertype==1) RETURN !No flushing for ilake

      !Sediment pool density is calculated for SS/AE/TS only. So far no other sediment pools has a density defined. 
      !Since density is the same, the depth relation is the same as the mass relation.
      IF(lakestate%sedpool(i_ss,watertype,i) + lakestate%sedpool(i_ae,watertype,i) >0.) &
      ss_frac = lakestate%sedpool(i_ss,watertype,i) / (lakestate%sedpool(i_ss,watertype,i) + lakestate%sedpool(i_ae,watertype,i)) ! Save fraction of sediment depth that is SS
              
  END SUBROUTINE sediment_fraction

  !>\brief Calculate sedimentation and resuspension of PP and SS in rivers.
  !>
  !>\b Reference ModelDescription Chapter Nitrogen and phosphorus processes in rivers and lakes (Sedimentation/Resuspension)
  !-----------------------------------------------------------------
  SUBROUTINE river_sedimentation_resuspension(i,watertype,area,sedexppar,riverq,qbank,depth,maxSSconc,sedresSS,riverstate)

    USE HYPEVARIABLES, ONLY : m_qbank,m_vpeak,m_addsusp,m_suspconSS, &
                              m_suspexpSS,m_suspconPP,m_suspexpPP, &
                              m_sedch,m_hygeomm,m_hygeomk
    USE MODVAR, ONLY : basin,simulate,i_pp,i_ss, &
                       genpar,modeloption,p_sedresusp, &
                       numsubstances,simulate

    !Argument declaration
    INTEGER, INTENT(IN) :: i          !<index of current subbasin
    INTEGER, INTENT(IN) :: watertype  !<river type (1=local, 2=main)
    REAL, INTENT(IN)    :: area       !<river surface area (m2)
    REAL, INTENT(IN)    :: sedexppar  !<sedimentation/resuspension parameter
    REAL, INTENT(IN)    :: riverq     !<river discharge (m3/s)
    REAL, INTENT(IN)    :: qbank      !<bank full flow (m3/s)
    REAL, INTENT(IN)    :: depth      !<river depth (m) 
    REAL, INTENT(OUT)   :: maxSSconc  !<river SS maximum transport concentration (mg/L)
    REAL, INTENT(OUT)   :: sedresSS   !<river sedimentation/resuspension of SS (kg/timestep)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states

    !Local variables
    REAL, DIMENSION(1) :: PPpool,SSpool    !pools in water (kg)
    REAL, DIMENSION(1) :: transport        !change in pools (kg/d)
    REAL, DIMENSION(1) :: tempsed,tempsed2 !temparary variables for river sediment pools (kg)
    REAL sedresp, sedrespPP,help, qbankcorr
    REAL riverv, vpeak, cmax(numsubstances)    !variables for SS maximum concentration and PP
    REAL sedch  !modification factor of resuspension (erodability/channel cover effects)

    !>\b Algorithm
    !Default output 
    maxSSconc = 0.;sedresSS = 0.
    !>Initial check if calculation is to be made
    IF(area<=0.) RETURN 
    sedch = genpar(m_sedch) !default value, general parameter

    !Sedimentation and resuspension of substance T1 is taken care of in other routine. It has only one model.
    !Could make that one smaller and us it to get transport. Then add/remove here with the other.
    !Warning:I cannot do this routine for all substances as a vector (as it is now) because T1 has its own subroutin to call. Maybe with T1=0?
    
    !>Assess pool of PP and SS in water and sediment
    IF(simulate%substance(i_pp))THEN
      tempsed(1) = riverstate%Psed(watertype,i)
      PPpool = (riverstate%water(watertype,i) * riverstate%conc(i_pp,watertype,i))* 1.0E-3 !kg
    ENDIF
    IF(simulate%substance(i_ss))THEN
      tempsed2(1) = riverstate%Ssed(watertype,i)
      SSpool = (riverstate%water(watertype,i) * riverstate%conc(i_ss,watertype,i))* 1.0E-3 !kg
    ENDIF

    !>Select current model for river resuspension/sedimentation
    SELECT CASE(modeloption(p_sedresusp))
      
    !>For model 0:
    CASE(0,1)  !original HYPE method, based om current flow related to bankful flow
      IF(sedexppar == 0) RETURN
      IF(qbank == 0.) RETURN

      !>Calculate sedimentation and resuspension factor (per day)
      IF(modeloption(p_sedresusp)==0)THEN
        qbankcorr = 0.7 * qbank
      ELSE
        qbankcorr = genpar(m_qbank) * qbank
      ENDIF
      help = 0.
      IF(qbankcorr-riverq>0.) help = help + ((qbankcorr-riverq)/qbankcorr)**sedexppar   !sedimentation at low flow
      IF(riverq>0) help = help - (riverq/qbankcorr)**sedexppar  !resuspension at all flows
      sedresp = MAX(-1., MIN(1.,help))
      
      !>Transfer PP and SS between sediment and water pools
      IF(sedresp>0.)THEN  !net sedimentation
        IF(simulate%substance(i_pp))THEN
          transport = sedresp * (riverstate%conc(i_pp,watertype,i) * MIN(riverstate%water(watertype,i),area * depth)) / 1.0E3  
          CALL retention_pool(1,PPpool,transport)           !transfer may change
          CALL production_pool(1,tempsed,transport)
        ENDIF
        IF(simulate%substance(i_ss))THEN
          transport = sedresp * (riverstate%conc(i_ss,watertype,i) * MIN(riverstate%water(watertype,i),area * depth)) / 1.0E3  
          CALL retention_pool(1,SSpool,transport)           !transfer may change
          CALL production_pool(1,tempsed2,transport)
          sedresSS = transport(1)   !output
        ENDIF
      ELSE                !net resuspension
        IF(simulate%substance(i_pp))THEN
          transport = - sedresp * tempsed 
          CALL retention_pool(1,tempsed,transport)          !transfer may change
          CALL production_pool(1,PPpool,transport)
        ENDIF
        IF(simulate%substance(i_ss))THEN
          transport = - sedresp * tempsed2
          CALL retention_pool(1,tempsed2,transport)         !transfer may change
          CALL production_pool(1,SSpool,transport)
          sedresSS = - transport(1)   !output
        ENDIF
      ENDIF
      
    !>For model 2: Simplified Bagnold Equation
    CASE(2)
      sedrespPP = 0.
      IF(basin(i)%channelfactor>=0.) sedch = basin(i)%channelfactor

      !>Calculate river peak velocity across channel
      riverv = genpar(m_hygeomk) * (riverq) ** genpar(m_hygeomm) ! Calculate river velocity
      vpeak = genpar(m_vpeak) * riverv ! Calculate peak channel velocity (m/s) = vpeak * average velocity (across the channel)
          
      !>Calculate amount of sediment (kg) to transfer, and/or amount of PP
      IF(simulate%substance(i_ss))THEN
        cmax(i_ss) = genpar(m_suspconSS) * (vpeak ** genpar(m_suspexpSS)) * 1.0E6 ! Calculate maximum amount of sediment that can be transported and convert from kg/L to mg/L
        maxSSconc = cmax(i_ss)  !output
        IF(riverstate%conc(i_ss,watertype,i) > cmax(i_ss)) THEN
          sedresp = (riverstate%conc(i_ss,watertype,i) - cmax(i_ss)) * MIN(riverstate%water(watertype,i),area * depth) / 1.0E3 ! amount of sediment (kg) to settle when concentration > max transport concentration
        ELSEIF(riverstate%conc(i_ss,watertype,i) < cmax(i_ss)) THEN
          sedresp = -1. * (cmax(i_ss) - riverstate%conc(i_ss,watertype,i)) * MIN(riverstate%water(watertype,i),area * depth) * sedch / 1.0E3 ! amount of sediment (kg) to resuspend when concentration < max transport concentration; multiply by -1. so it will be removed from pool
        ENDIF
      ENDIF
      IF(simulate%substance(i_pp))THEN
        cmax(i_pp) = genpar(m_suspconPP) * (vpeak ** genpar(m_suspexpPP)) * 1.0E6 ! Calculate maximum amount of sediment that can be transported and convert from kg/L to mg/L
        IF(riverstate%conc(i_pp,watertype,i) > cmax(i_pp)) THEN
          sedrespPP = (riverstate%conc(i_pp,watertype,i) - cmax(i_pp)) * MIN(riverstate%water(watertype,i),area * depth) / 1.0E3 ! amount of sediment (kg) to settle when concentration > max transport concentration
        ELSEIF(riverstate%conc(i_pp,watertype,i) < cmax(i_pp)) THEN
          sedrespPP = -1. * (cmax(i_pp) - riverstate%conc(i_pp,watertype,i)) * MIN(riverstate%water(watertype,i),area * depth) * sedch / 1.0E3 ! amount of sediment (kg) to resuspend when concentration < max transport concentration; multiply by -1. so it will be removed from pool
        ENDIF
      ENDIF
      
      !>Transfer SS between sediment and water pools
      IF(simulate%substance(i_ss))THEN
        IF(sedresp>0.)THEN  !net sedimentation
          transport = sedresp 
          CALL retention_pool(1,SSpool,transport)           !transport may change
          CALL production_pool(1,tempsed2,transport)
          sedresSS = transport(1)    !output
        ELSEIF(sedresp<0.)THEN  !net resuspension
          IF(-sedresp > tempsed2(1)) THEN! #CB: resuspend full sediment pool and then use addsed (0.0 - 1.0) to control how much erosion in excess of temporary pool can occur
            transport = tempsed2 + ((-sedresp - tempsed2) * genpar(m_addsusp))
            tempsed2 = 0.
          ELSE
            transport = -sedresp
            CALL retention_pool(1,tempsed2,transport)         !transfer may NOT change, because already checked
          ENDIF
          CALL production_pool(1,SSpool,transport)
          sedresSS = - transport(1)   !output
        ENDIF
      ENDIF
      
      !>Transfer PP between sediment and water pools
      IF(simulate%substance(i_pp))THEN
        IF(sedrespPP>0.)THEN  !net sedimentation
          transport = sedrespPP 
          CALL retention_pool(1,PPpool,transport)           !transport may change
          CALL production_pool(1,tempsed,transport)
        ELSEIF(sedrespPP<0.)THEN  !net resuspension
          IF(-sedrespPP > tempsed(1)) THEN  !resuspend full temporary pool and then use addsed (0.0 - 1.0) to control how much erosion in excess of temporary pool can occur
            transport = tempsed + ((-sedrespPP - tempsed) * genpar(m_addsusp))
            tempsed = 0.
          ELSE
            transport = -sedrespPP
            CALL retention_pool(1,tempsed,transport)         !transfer may NOT change, because already checked
          ENDIF
          CALL production_pool(1,PPpool,transport)
        ENDIF
      ENDIF

    END SELECT

    !>Update state variables
    IF(simulate%substance(i_pp))THEN
      riverstate%Psed(watertype,i) = tempsed(1)
      CALL new_concentration(PPpool(1),riverstate%water(watertype,i)*1.0E-3,riverstate%conc(i_pp,watertype,i))
    ENDIF
    IF(simulate%substance(i_ss))THEN
      riverstate%Ssed(watertype,i) = tempsed2(1)
      CALL new_concentration(SSpool(1),riverstate%water(watertype,i)*1.0E-3,riverstate%conc(i_ss,watertype,i))
    ENDIF

  END SUBROUTINE river_sedimentation_resuspension

  !>\brief Calculates straight 365-day running average mean of TP
  !>concentration in lake
  !>
  !>\b Reference ModelDescription Chapter Nitrogen and phosphorus processes in rivers and lakes (Primary production 
  !> and mineralization) and Organic carbon (River and Lakes - Primary production and mineralization)
  !-----------------------------------------------------------------------
  SUBROUTINE calculate_lake_tpmean(i,watertype,lakestate)

    USE MODVAR, ONLY : i_sp,i_pp, &
                       timesteps_per_day

    !Argument declarations
    INTEGER, INTENT(IN) :: i         !<index of current subbasin
    INTEGER, INTENT(IN) :: watertype !<Lake type (1=local, 2=outlet)
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    
    IF(lakestate%water(watertype,i)>0.)THEN
      lakestate%TPmean(watertype,i) = lakestate%TPmean(watertype,i) + (lakestate%conc(i_sp,watertype,i)+lakestate%conc(i_pp,watertype,i) - lakestate%TPmean(watertype,i))/(365.* timesteps_per_day)
    ENDIF

  END SUBROUTINE calculate_lake_tpmean

  !>\brief Calculates straight 365-day running average mean of TP
  !>concentration in river
  !>
  !>\b Reference ModelDescription Chapter  Nitrogen and phosphorus processes in rivers and lakes (Primary production 
  !> and mineralization) and Organic carbon (River and Lakes - Primary production and mineralization)
  !-------------------------------------------------------------------
  SUBROUTINE calculate_river_tpmean(i,watertype,riverstate)

    USE MODVAR, ONLY : i_sp,i_pp, &
                       timesteps_per_day

    !Argument declarations
    INTEGER, INTENT(IN) :: i         !<index of current subbasin
    INTEGER, INTENT(IN) :: watertype !<River type (1=local, 2=main)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states

    riverstate%TPmean(watertype,i) = riverstate%TPmean(watertype,i) + (riverstate%conc(i_sp,watertype,i) &
         + riverstate%conc(i_pp,watertype,i) - riverstate%TPmean(watertype,i))/(365.* timesteps_per_day)

  END SUBROUTINE calculate_river_tpmean

  !>\brief Calculates and add internal load of phosphorus for lakes
  !>
  !>\b Reference ModelDescription Chapter Nitrogen and phosphorus processes in rivers and lakes (Internal load)
  !-------------------------------------------------------------
  SUBROUTINE internal_lake_load(i,watertype,systemtype,area,lakestate)

    USE MODVAR, ONLY : lakeindex,   &
                       lakedatapar, &
                       lakedataparindex,  &
                       i_sp,i_pp,   &
                       numsubstances
    USE HYPEVARIABLES, ONLY : m_ldprodpp,  &
                              m_ldprodsp       

    !Argument declarations
    INTEGER, INTENT(IN) :: i           !<index of current subbasin
    INTEGER, INTENT(IN) :: watertype   !<Lake or river type (1=local, 2=main/outlet)
    INTEGER, INTENT(IN) :: systemtype  !<aquatic system type (1=lake, 2=river)
    REAL, INTENT(IN)    :: area        !<lake surface area/ river bottom area (m2)
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state

    !Local variables
    REAL prodPP, prodSP
    REAL tmpfcn, TPfcn
    REAL vol
    REAL pppar,sppar
    REAL :: sourceP(numsubstances)
    
    !Local parameters
    INTEGER, PARAMETER :: pooldim = 1

    !>\b Algorithm \n
    !>Check if internal phosphorus load is to be calculated
    IF(systemtype==2) RETURN   !river
    IF(watertype==1) RETURN    !local
    IF(.NOT.ALLOCATED(lakeindex)) RETURN  !no special lakes
    pppar = lakedatapar(lakedataparindex(i,watertype),m_ldprodpp)
    sppar = lakedatapar(lakedataparindex(i,watertype),m_ldprodsp)
    IF(pppar==0 .AND. sppar==0) RETURN
    sourceP = 0.

    !>Calculate pool of P, and concentration and temperature dependent factors
    TPfcn = 0.1 !mg/L
    tmpfcn = 0.86**(ABS(lakestate%temp(watertype,i)-15.))   !laketemp=T20 for olake

    !> Calculate internal load of phosphorus
    prodPP = pppar * TPfcn * tmpfcn * area / 1000.  !kg/d
    prodSP = sppar * TPfcn * tmpfcn * area / 1000.  !kg/d
    sourceP(i_pp) = prodPP
    sourceP(i_sp) = prodSP

    !>Add internal load of phosphorus to lake water
    vol = lakestate%water(watertype,i) * area / 1.0E6
    CALL add_source_to_water(vol,numsubstances,lakestate%conc(:,watertype,i),sourceP)

  END SUBROUTINE internal_lake_load

  !>\brief Calculates transformation between OC/DIC in water 
  !!Simulating the combined processes of primary production and
  !!mineralisation.
  !>
  !>\b Reference ModelDescription Organic carbon (River and lakes - Primary production 
  !> and mineralization)
  !----------------------------------------------------------------
  SUBROUTINE oc_production_mineralisation(systemtype,area,prodpar,halfsatTPwater, &
                             limpppar,water,conc,watertemp,waterTPmean,temp10,temp20,depth)

    USE HYPEVARIABLES, ONLY : maxdegradwater, &
                              NCratio

    !Argument declarations
    INTEGER, INTENT(IN)        :: systemtype  !<aquatic system type (1=lake, 2=river)
    REAL, INTENT(IN)           :: area        !<lake surface area/ river bottom area (m2)
    REAL, INTENT(IN)           :: prodpar     !<model parameter production rate OC in water
    REAL, INTENT(IN)           :: halfsatTPwater !<model parameter half saturation TP (mg/L)
    REAL, INTENT(IN)           :: limpppar    !<limitation of sedimentation parameter (mg/L)
    REAL, INTENT(IN)           :: water       !<river or lake water (mm or m3)
    REAL, INTENT(INOUT)        :: conc        !<OC concentration of river or lake
    REAL, INTENT(IN)           :: watertemp   !<water temperature
    REAL, INTENT(IN)           :: waterTPmean !<water TP mean
    REAL, INTENT(IN)           :: temp10      !<10-day water temperature
    REAL, INTENT(IN)           :: temp20      !<20-day water temperature
    REAL, INTENT(IN), OPTIONAL :: depth       !<river depth (m) 
    
    !Local variables
    REAL, DIMENSION(1) :: OCpool, minprodN, minprodC,minC,prodC
    REAL tmpfcn, tmpfcn1, tmpfcn2, TPfcn
    REAL vol
    REAL waterdepth !(m)
    
    !Local parameter
    INTEGER, PARAMETER :: pooldim = 1

    !>\b Algorithm \n
    !>Calculate pools of organic carbon in the water, water temperature 
    !>and fraction of depth of water volume that is active
    IF(systemtype==1) THEN   !lakes
      OCpool = (water * area * conc) /1.0E6  !kg
      waterdepth = water/1000.
    ELSE                     !rivers
      OCpool = (water * conc)/ 1.0E3 !kg
      waterdepth=depth
    ENDIF

    !>Calculate dependency factors (Tot-P and temperature)
    TPfcn = halfsatconcfactor(waterTPmean-limpppar,halfsatTPwater)
    IF(watertemp >= 0.) THEN
      tmpfcn1 = watertemp / 20.    
    ELSE 
      tmpfcn1 = 0.
    ENDIF
    tmpfcn2 = (temp10 - temp20) / 5.
    tmpfcn = tmpfcn1*tmpfcn2

    !>Calculate production/mineralisation of organic carbon
    minprodN = 0.
    IF(watertemp > 0. ) THEN 
      minprodN = prodpar * TPfcn * tmpfcn * waterdepth * area  !kg  
      IF(minprodN(1) > 0.) THEN  !production        
        minprodC = minprodN * NCratio
      ELSE                       !mineralisation
        minprodC = MAX(-maxdegradwater * OCpool, minprodN * NCratio)
      ENDIF
    ENDIF
    minC = -minprodC
    prodC = minprodC
    IF(minprodC(1)>0.) CALL production_pool(pooldim,OCpool,prodC)
    IF(minprodC(1)<0.) CALL retention_pool(pooldim,OCpool,minC)

    !>Set new concentration due to changes in pools
    IF(systemtype==1) THEN            !lakes
      vol = water * area / 1.0E6
      CALL new_concentration(OCpool(1),vol,conc)
    ELSE                                 !rivers
      IF(water > 0.) THEN
        vol = water / 1.0E3
        CALL new_concentration(OCpool(1),vol,conc)
      ENDIF
    ENDIF

  END SUBROUTINE oc_production_mineralisation

  !>Add load from local diffuse sources to local river inflow
  !>
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land 
  !!routines (Nutrient sources - Rural household diffuse source)
  !-----------------------------------------------------------------
  SUBROUTINE add_diffuse_source_to_local_river(i,qin,cin,source,addedflow)

    USE MODVAR, ONLY : load,                &
                       genpar,              &
                       numsubstances,       &
                       seconds_per_timestep, &
                       conductbasinlocsoil
    USE HYPEVARIABLES, ONLY : m_locsoil

    !Argument declarations
    INTEGER, INTENT(IN) :: i                      !<index of subbasin
    REAL, INTENT(INOUT) :: qin                    !<flow in local river (m3/s)
    REAL, INTENT(INOUT) :: cin(numsubstances)     !<concentration of flow into local river (mg/L)
    REAL, INTENT(OUT)   :: source(numsubstances)  !<local source added to local river (kg/timestep)
    REAL, INTENT(OUT)   :: addedflow              !<added flow (m3/timestep)
    
    !Local variables
    REAL qhelp
    REAL qadd
    REAL cadd(numsubstances)

    !Initiation
    source = 0.
    cadd = 0.

    !> \b Algorithm \n
    !>Calculate diffuse source from rural households to local river
    IF(conductbasinlocsoil)THEN
      addedflow = (1. - load(i)%locsoil) * load(i)%volloc   !m3/ts
    ELSE
      addedflow = (1. - genpar(m_locsoil)) * load(i)%volloc   !m3/ts
    ENDIF
    qadd = addedflow / seconds_per_timestep   !m3/s !fel vid korta tidsteg!
    IF(qadd>0)THEN
      qhelp = qadd * seconds_per_timestep * 1.E-3   !1000m3/timestep
      cadd = load(i)%locconc
      source = cadd * qhelp    !Diffuse load, ruralB, (kg/timestep for NPS, U/ts for T1)

      !>Add diffuse source to inflow to local river flow
      IF(qin>0)THEN
        cin = (qin * cin + qadd * cadd)/(qin + qadd)
        qin = qin + qadd
      ELSE
        qin = qadd
        cin = cadd
      ENDIF
    ENDIF

  END SUBROUTINE add_diffuse_source_to_local_river

  !>Add load from point sources to main river inflow
  !>
  !> \b Reference ModelDescription Chapter Water management (Point sources)
  !-----------------------------------------------------------------
  SUBROUTINE add_point_sources_to_main_river(isb,qin,cin,source,addedflow)

    USE MODVAR, ONLY : max_pstype,          &
                       load,                &
                       npsused,psinfo,psload, &
                       numsubstances,       &
                       seconds_per_timestep, &
                       find_next_pointsource

    !Argument declarations
    INTEGER, INTENT(IN) :: isb                      !<index of subbasin
    REAL, INTENT(INOUT) :: qin                      !<flow into main river (m3/s)
    REAL, INTENT(INOUT) :: cin(numsubstances)       !<concentration of flow into main river (mg/L)
    REAL, INTENT(OUT)   :: source(numsubstances,max_pstype)  !<point sources added to main river (kg/timestep)
    REAL, INTENT(OUT)   :: addedflow                !<added flow (m3/timestep)
    
    !Local variables
    INTEGER k,j
    INTEGER ips,lastps
    REAL divvolps
    REAL qadd
    REAL cadd(numsubstances)

    !Initiation
    source = 0.
    addedflow = 0.
    cadd = 0.
    qadd = 0.

    IF(.NOT.ALLOCATED(psinfo))THEN
      IF(.NOT.ALLOCATED(load(isb)%psvol)) RETURN
    
      !Calculate permanent source to be added to river
      DO k = 1,max_pstype
        qadd = qadd + load(isb)%psvol(k)   !m3/s
      ENDDO
      addedflow = qadd * seconds_per_timestep
      IF(qadd>0)THEN
        divvolps = 1000./qadd/seconds_per_timestep                    !kg/ts,m3/s->mg/L
        DO k = 1,max_pstype
          DO j = 1,numsubstances
            cadd(j) = cadd(j) + load(isb)%psload(k,j)
            source(j,k) = load(isb)%psload(k,j)        !Point source k, substance j
          ENDDO
        ENDDO
        IF(numsubstances>0) cadd(:) = cadd(:) * divvolps    !mg/L

        !Add source to river      
        IF(qin>0)THEN
          cin = (qin * cin + qadd * cadd)/(qin + qadd)
          qin = qin + qadd
        ELSE
          qin = qadd
          cin = cadd
        ENDIF
      ENDIF
    ELSE
      !Use time dependent pointsources
      !Calculate source to be added to river
      lastps = 0
      DO
        CALL find_next_pointsource(isb,lastps,npsused,psinfo,ips)
        IF(ips==0) EXIT
        lastps = ips
        k = psinfo(ips)%pstype
        IF(psinfo(ips)%sw_code/=3) CYCLE
        IF(psload(ips)%flow<=0.) CYCLE
        qadd = qadd + psload(ips)%flow   !m3/s
        cadd = cadd + psload(ips)%load   !g/s
        IF(numsubstances>0) source(:,psinfo(ips)%pstype) = source(:,psinfo(ips)%pstype) + psload(ips)%load(:)*seconds_per_timestep*1.E-3 !Point source kg/timestep
      ENDDO
      addedflow = qadd * seconds_per_timestep
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

  END SUBROUTINE add_point_sources_to_main_river

  !>Calculate effect of river wetland constructed for nutrient removal
  !>
  !>\b Reference ModelDescription Chapter Water management (Constructed wetlands)  
  !-----------------------------------------------------------------
  SUBROUTINE calculate_river_wetland(i,itype,n,temp5,temp30,qin,cin,cwetland)

    USE MODVAR, ONLY : wetland,     &
                       seconds_per_timestep

    !Argument declarations
    INTEGER, INTENT(IN) :: i           !<index of subbasin
    INTEGER, INTENT(IN) :: itype       !<index of river type (local or main)
    INTEGER, INTENT(IN) :: n           !<number of substances
    REAL, INTENT(IN)    :: temp5       !<temperature (5-day-mean) (degree Celsius)
    REAL, INTENT(IN)    :: temp30      !<temperature (30-day-mean) (degree Celsius)
    REAL, INTENT(IN)    :: qin         !<flow into/out of river wetland (m3/s)
    REAL, INTENT(INOUT) :: cin(n)      !<concentration of flow into/out of river wetland (mg/L)
    REAL, INTENT(INOUT) :: cwetland(n) !<concentration of river wetland (mg/L)
    
    !Local variables
    REAL wetlandvol     !m3 (constant)
    REAL wetlandinflow  !m3/timestep

    !Start of calculations
    IF(wetland(i,itype)%area==0) RETURN   !no wetland

    wetlandvol = wetland(i,itype)%area * wetland(i,itype)%depth   !m3
    wetlandinflow = qin * wetland(i,itype)%part * seconds_per_timestep     !m3/timestep
    CALL calculate_riverwetland_np(n,wetlandinflow,cin,wetland(i,itype)%area,wetlandvol,cwetland,temp5,temp30)
    IF(qin>0) cin = cin * (1. - wetland(i,itype)%part) + cwetland * wetland(i,itype)%part           !New concentration

  END SUBROUTINE calculate_river_wetland

  !>\brief Calculate nutrient processes in river wetland. 
  !!Retention is limited to 99.9% of the pool.
  !>
  !>\b Reference ModelDescription Chapter Water management (Constructed wetlands)  
  !------------------------------------------------------------------------
  SUBROUTINE calculate_riverwetland_np(n,qin,cin,area,vol,cvol,temp5,temp30)

    USE MODVAR, ONLY : i_in,i_sp,i_pp,  &
                       conduct

    !Argument declarations
    INTEGER, INTENT(IN) :: n       !<number of substances
    REAL, INTENT(IN)    :: qin     !<flow into wetland (m3/d)
    REAL, INTENT(IN)    :: cin(n)  !<concentration of river flow (mg/l) (before and after wetland processes
    REAL, INTENT(IN)    :: area    !<area of wetland (m2)
    REAL, INTENT(IN)    :: vol     !<volume of wetland (m3)
    REAL, INTENT(INOUT) :: cvol(n) !<concentration of wetland volume (mg/l) (before and after wetland processes
    REAL, INTENT(IN)    :: temp5   !<temperature (5-day-mean) (degree Celsius)
    REAL, INTENT(IN)    :: temp30  !<temperature (30-day-mean) (degree Celsius)
    
    !Local variables
    REAL wetlandnutrient(n), wetlandconc(n)
    REAL retention(n)
    REAL retention_tp, production_tp
    REAL wetland_tp,srpfrac
    
    !Local parameters
    REAL, PARAMETER :: teta = 1.2
    REAL, PARAMETER :: tkoeff = 20.   !temperature coefficient (degree Celsius)
    REAL, PARAMETER :: inpar = 2.3    !model parameter for inorganic nitrogen retention (mm/d/degree Celsius)
    REAL, PARAMETER :: sedpar = 0.09  !model parameter for phosphorus sedimentation (m/d)
    REAL, PARAMETER :: uptpar = 0.1   !model parameter for phosphorus uptake (m/d)

    !Calculate the nutrient processes
    wetlandnutrient = vol*cvol+qin*cin         !g
    wetlandconc = wetlandnutrient /(vol+qin)   !mg/l
    retention = 0.
    IF(conduct%simN)THEN
      IF(temp5>0) retention(i_in) = inpar * wetlandconc(i_in) * area * temp5 * 1.E-3         !g/d denitrification
      IF(retention(i_in)<0) retention(i_in) = 0.
      IF(retention(i_in)>0.999*wetlandnutrient(i_in)) retention(i_in) = 0.999 * wetlandnutrient(i_in)
    ENDIF
    IF(conduct%simP)THEN
      retention_tp = sedpar * (wetlandconc(i_pp) + wetlandconc(i_sp)) * area                 !g/d sedimentation
      IF(retention_tp<0) retention_tp = 0.
      production_tp = uptpar * (cin(i_pp) + cin(i_sp)) * (teta ** (temp30 - tkoeff)) * area  !g/d uptake
      IF(production_tp<0) production_tp = 0.
      wetland_tp = wetlandnutrient(i_pp) + wetlandnutrient(i_sp)    !g
      IF(retention_tp - production_tp < 0.999 * wetland_tp)THEN
        srpfrac = wetlandnutrient(i_sp) / wetland_tp
        retention(i_sp) = srpfrac * (retention_tp - production_tp)
        retention(i_pp) = (1.-srpfrac) * (retention_tp - production_tp)
      ELSE
        retention_tp = 0.999 * wetland_tp
        IF(wetland_tp>0)THEN
          srpfrac = wetlandnutrient(i_sp)/wetland_tp
        ELSE
          srpfrac = 0.
        ENDIF
        retention(i_sp) = srpfrac * retention_tp
        retention(i_pp) = (1.-srpfrac) * retention_tp
      ENDIF
    ENDIF
    cvol = (wetlandnutrient - retention)/(vol+qin)    !New concentration of wetland volume

  END SUBROUTINE calculate_riverwetland_np

  !>\brief Calculate processes for substances in wetland. 
  !!Retention is limited to 99.9% (sed) or 50% (uptake) of the pool.
  !>
  !>\b Reference ModelDescription Chapter Water management (Constructed wetlands)  
  !------------------------------------------------------------------------
  SUBROUTINE wetland_substance_processes(n,area,vol,cvol,temp5,temp30,fastN,fastP,humusN,humusP,partP)

    USE MODVAR, ONLY : i_in,i_on,i_sp,i_pp,i_ss,  &
                       conduct,timesteps_per_day, &
                       genpar
    USE HYPEVARIABLES, ONLY : m_wlsed,m_wlproddep,m_wlmphuptin,m_wlmphuptsp,m_wlfastfrac,m_wlpartfrac,m_wltmpexp

    !Argument declarations
    INTEGER, INTENT(IN) :: n       !<number of substances
    REAL, INTENT(IN)    :: area    !<area of wetland (m2)
    REAL, INTENT(IN)    :: vol     !<volume of wetland (m3)
    REAL, INTENT(INOUT) :: cvol(n) !<concentration of wetland volume (mg/l) (before and after wetland processes)
    REAL, INTENT(IN)    :: temp5   !<temperature (5-day-mean) (degree Celsius)
    REAL, INTENT(IN)    :: temp30  !<temperature (30-day-mean) (degree Celsius)
    REAL, INTENT(INOUT) :: fastN   !<immobile fast turnover nitrogen in upper soil layer(kg/km2)
    REAL, INTENT(INOUT) :: fastP   !<immobile fast turnover phosphorus in upper soil layer (kg/km2)
    REAL, INTENT(INOUT) :: humusN  !<immobile slow turnover nitrogen in upper soil layer(kg/km2)
    REAL, INTENT(INOUT) :: humusP  !<immobile slow turnover phosphorus in upper soil layer (kg/km2)
    REAL, INTENT(INOUT) :: partP   !<immobile particulate phosphorus in upper soil layer (kg/km2)
    
    !Local variables
    REAL sedvel
    REAL wetlandnutrient(n), wetlandconc(n)
    REAL retention(n)
    REAL sedimentation_pp,sedimentation_ss,sedimentation_on, macrouptake_in, macrouptake_sp
    REAL tmpfcn1, tmpfcn2, tmpfcn,fracarea,waterTPmean,TPfcn
    
    !Local parameters
    REAL, PARAMETER :: halfsatTPwater = 0.05
    
    IF(vol<=0.)RETURN
    
    !Initialisation of wetland variables
    wetlandnutrient = vol*cvol         !g
    wetlandconc = cvol   !mg/l
    retention = 0.
    
    !Fractional area of wetland with macrophyte uptake
    fracarea = MIN(1.0, genpar(m_wlproddep) * (1/((vol/area)*2)))
    !Sedimentation velocity
    sedvel = genpar(m_wlsed)/REAL(timesteps_per_day)   !1/ts
    !Temperature dependent factor
    tmpfcn1 = (MAX(0.,temp5) / 20.)**genpar(m_wltmpexp)
    tmpfcn2 = (temp5 - temp30) / 5.
    tmpfcn = MAX(0.,tmpfcn1*tmpfcn2)
    !Total phosphorus concentration dependent factor
    IF(conduct%simP)THEN
      waterTPmean = wetlandconc(i_pp) + wetlandconc(i_sp)
      TPfcn = halfsatconcfactor(waterTPmean,halfsatTPwater)
    ELSE
      TPfcn = 0.5
    ENDIF
    
    IF(conduct%simN)THEN
      !Denitrification in soil water (including wetland water volume) calculated in the soil routines.
      
      !Sedimentation
      sedimentation_on = sedvel * (wetlandconc(i_on)) * area                 !g/ts sedimentation
      IF(sedimentation_on > 0.999 * wetlandnutrient(i_on)) sedimentation_on = 0.999 * wetlandnutrient(i_on)
      retention(i_on) = retention(i_on) + sedimentation_on
      fastN = fastN + 1.E3*sedimentation_on/area
     
      !Uptake of IN by macrophytes
      macrouptake_in = genpar(m_wlmphuptin) * tmpfcn * fracarea * area * TPfcn
      IF(macrouptake_in > 0.5 * wetlandnutrient(i_in)) macrouptake_in = 0.5 * wetlandnutrient(i_in)   
      retention(i_in) = retention(i_in) + macrouptake_in
      fastN = fastN + genpar(m_wlfastfrac)*(1.E3*macrouptake_in/area)
      humusN = humusN + (1-genpar(m_wlfastfrac))*(1.E3*macrouptake_in/area)
    ENDIF
    
    IF(conduct%simP)THEN
      !Sedimentation
      sedimentation_pp = sedvel * (wetlandconc(i_pp)) * area                 !g/ts sedimentation
      IF(sedimentation_pp > 0.999 * wetlandnutrient(i_pp)) sedimentation_pp = 0.999 * wetlandnutrient(i_pp)
      retention(i_pp) = retention(i_pp) + sedimentation_pp
      fastP = fastP + (1.-genpar(m_wlpartfrac))*1.E3*sedimentation_pp/area
      partP = partP + genpar(m_wlpartfrac)*1.E3*sedimentation_pp/area
      
      !Uptake of SP by macrophytes
      macrouptake_sp = genpar(m_wlmphuptsp) * tmpfcn * fracarea * area * TPfcn
      IF(macrouptake_sp > 0.5 * wetlandnutrient(i_sp)) macrouptake_sp = 0.5 * wetlandnutrient(i_sp)   
      retention(i_sp) = retention(i_sp) + macrouptake_sp
      fastP = fastP + genpar(m_wlfastfrac)*(1.E3*macrouptake_sp/area)
      humusP = humusP + (1-genpar(m_wlfastfrac))*(1.E3*macrouptake_sp/area)
    ENDIF
    
    IF(conduct%simS)THEN
      !Sedimentation
      sedimentation_ss = sedvel * (wetlandconc(i_ss)) * area                 !g/ts sedimentation
      IF(sedimentation_ss > 0.999 * wetlandnutrient(i_ss)) sedimentation_ss = 0.999 * wetlandnutrient(i_ss)
      retention(i_ss) = retention(i_ss) + sedimentation_ss
    ENDIF
    
    cvol = (wetlandnutrient - retention)/vol    !New concentration of wetland water volume due to sedimentation

  END SUBROUTINE wetland_substance_processes


END MODULE NPC_SURFACEWATER_PROCESSES
