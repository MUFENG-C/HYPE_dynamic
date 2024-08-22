!> \file main.f90
!> Contains main program of HYSS-HYPE

!> Main program for HYSS - Hydrological Simulation System
!> 
!> The system features hydrological simulation of soil and surface water system, 
!> calibration of parameters and criteria calculations, updating of flow and state 
!> to observed values, ensemble simulation, and more.
PROGRAM HYSS

!Copyright 2011-2021 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

!-----------------------------------------------------------------------------------------

!Used modules
!      USE IFQWIN                               !for nobutton
      USE MODELMODULE, ONLY : model, &
                              model_version_information, &
                              define_output_variables, &
                              set_modelconfig_from_parameters, &
                              initiate_model, &
                              initiate_model_state, &
                              define_model_parameters, &
                              load_modeldefined_input, &
                              reallocate_nclass_based_variables
      USE STATETYPE_MODULE
      USE WORLDVAR, ONLY :  writematlab,       &
                            writeload,         &
                            output,            &
                            noutput,           &
                            nacrit,            &
                            nsubCrit,          &
                            ndt,               &
                            fileunit_temp,     &
                            fileunit_tests,    &
                            filename_best,     &
                            filename_MC,       &
                            maxcharpath,       &
                            simsequence,       &
                            parseq,            &
                            instate, &
                            infodir,           &
                            modeldir,          &
                            resdir,            &
                            forcingdir,        &
                            logdir,            &
                            bdate,             &
                            sdate,             &
                            outstartdate,      &
                            dtskip,            &
                            doens,doopt,       &
                            numoptimpar,       &   
                            deallocate_worldvar,      &
                            deallocate_MCvariables,   &
                            optim,             &
                            bestMCoptcrit,     &
                            bestMCperformance, &
                            bestMCparameters,  &
                            maxperf,           &
                            maxsubass,         &
                            simsubmodel,       &
                            basemodel,     &
                            psdates,    &
                            optimStartTime, &      
                            optimFuncCall,  &
                            lineSearchCallCount, &
                            doassimilation, &
                            resetstate, &
                            noutreg, &
                            da_allocate_accumulation,  &
                            allocate_accumulation,     &
                            reallocate_outvar_information, &
                            usestop84, &
                            indatacheckonoff, &
                            indatachecklevel, &
                            outstate, &
                            status_write_outstate, &
                            steplen, &
                            readadddate
      USE MODVAR, ONLY : nsub,ncrop,        &
                         nclass,numsubstances, &
                         naquifers,   &
                         dim_update, &
                         maxsoillayers,   &
                         max_classoutvar,   &
                         max_basinoutvar,   &
                         max_noutvar, &
                         conduct,statesize,  &
                         conductxoms,conductregest,   &
                         simulate, &
                         allocate_outvar,&
                         deallocate_modvar,  &
                         current_time, &
                         noutvar,noutvarclass, &
                         nrivertypes,nlaketypes, &
                         outvar   !for safety not lost
      USE COMPOUT, ONLY : compute_mapoutput,          &
                          compute_outloads,         &
                          prepare_to_compute_crit,  &
                          calculate_criteria
      USE TIMEROUTINES, ONLY : calculate_time_for_model
      USE READWRITE_ROUTINES
      USE LIBDATE, ONLY : DateType, format_date, AddDates, OPERATOR(.EQ.)
      USE DATAMODULE
      USE OPTIMIZATION
      USE STATE_DATAMODULE, ONLY : initiate_state_for_submodel, &
                                   load_saved_state,  &
                                   save_soil_state_file, &
                                   reset_soil_state, &
                                   finalize_outstate
      USE MODEL_TEST_ROUTINES, ONLY : stop_simulation_and_finalize_tests,setup_for_hype_tests,finalize_hype_tests, &
                                     run_hype_tests,run_hype_observation_tests
#ifdef _ASSIMILATION
      !use the Data Assimilation modules (pre-compiler flag ASSIMILATION to be set in Visual Studio project settings or in the makefile)
      USE ASSIMILATION_INTERFACE
      USE ASSIMILATION_ROUTINES
      USE ASSIMILATION_VARIABLES
#endif

      IMPLICIT NONE

!Parameter declarations
      INTEGER, PARAMETER :: maxreadstates = 100   !Max number of dates for reading soil states

!Variable declarations
      TYPE(DateType) d            !Current time
      INTEGER idt                 !Current timestep
      INTEGER ivar                !Current output variable
      INTEGER iens                !Current ensemble being simulated
      INTEGER iout                !Current output
!      INTEGER nobutton            !No exit window
      LOGICAL pwrite              !Flag for periodend, time to write to file
      LOGICAL conductinibin       !Flag for initializing assimilation from bin-files
      LOGICAL alive               !change geo or not
      CHARACTER(LEN=maxcharpath+25) filename  !hyss filename
      CHARACTER(LEN=maxcharpath+25) filename_year  !geodata/geoclass change filename
      CHARACTER(LEN=12) :: current_year
      CHARACTER(LEN=16) :: strdate  !Date for printing
      CHARACTER(LEN=8)  :: logdate  !Date for log-file name
      CHARACTER(LEN=10) :: logtime  !Time for log-file name
      CHARACTER(LEN=3)  :: logseq   !Seqnr for log-file name
      INTEGER :: datim(8) 
      INTEGER :: oldyear            !year of last time step
     
      REAL, ALLOCATABLE :: par(:)
      REAL optcrit, condcrit, condthres
      REAL, ALLOCATABLE :: basincrit(:,:,:)   !R2, CC, RE, RSDE, QC, QR, STDC, STDR, MAE, RMSE, Bias, STDbias, KGE, KGEpartSTD, KGEpartMM, NRMSE per subbasin och kriterie
      REAL, ALLOCATABLE :: simperformance(:,:)   !rr2,sr2,wr2,rmae,sbias,rrve,wrve,rra,sra,meanRA,tau,medianr2,medianra,meanrs,meancc,mediankg,meanabsre
      REAL, ALLOCATABLE :: ensemble_parameters(:,:)   !<Parameter values to be used for ensemble simulation
      INTEGER, ALLOCATABLE :: subincrit(:)      !Subbasins to be included in criteria calculations
      INTEGER :: status       !Subroutine return status
      INTEGER npar            !Number of parameters to be calibrated (couted in file)
      INTEGER :: nmapperiod   !Number of periods for map print out
      INTEGER :: numreadstates                            !Number of dates for reading soil state
      TYPE(DateType) :: readstatedate(maxreadstates)      !Dates for reading soil state
      
!Variables for updating of Q, W and concentrations
      LOGICAL :: update_allstations(dim_update), update_nostations(dim_update)
      
!Model state variables and other saved variables declaration
      TYPE(SNOWICESTATETYPE) :: frozenstate
      TYPE(SOILSTATETYPE)    :: soilstate      
      TYPE(AQUIFERSTATETYPE) :: aquiferstate      
      TYPE(RIVERSTATETYPE)   :: riverstate      
      TYPE(LAKESTATETYPE)    :: lakestate
      TYPE(MISCSTATETYPE)    :: miscstate

#ifdef _ASSIMILATION
      !Declaration of some data assimilation variables
      INTEGER assim_ens_size  !number of ensemble members
      INTEGER iassim     !loop-variable for ensembles members
      TYPE(STATEINFOTYPE),ALLOCATABLE :: stateinfo(:)
#endif

!Program start
!>\b Algorithm \n

      CALL DATE_AND_TIME (logdate, logtime,values=datim)

!Current model domain
      CALL model_version_information(0)   !Print model version information on screen
      CALL get_hyss_arguments(infodir,simsequence,parseq)   
      WRITE(logseq,'(I3.3)') simsequence
      WRITE(filename,'(a)') TRIM(infodir)//'hyss_'//logseq(1:3)//'_'//logdate(3:8)//'_'//logtime(1:8)//'.log'
      OPEN(UNIT=6,FILE=TRIM(filename),STATUS = 'replace',ACTION='write',ERR=900)
      WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
               ' Job start date: ',datim(1),'-',datim(2),'-',datim(3),    &
               '  time: ',datim(5),':',datim(6),':',datim(7)
      WRITE(6,*) '---------------------------------------------------'
!      nobutton = SETEXITQQ(qwin$exitnopersist)    !nobutton version
      CALL model_version_information(6)   !Print model version information in logg-file
      CALL define_output_variables()
      CALL define_model_parameters()
!>Get information for this simulation
      CALL load_coded_info(infodir,status,bdate,sdate,outstartdate,dtskip,ndt,numsubstances,  &
                           maxreadstates,numreadstates,readstatedate,modeldir,resdir,forcingdir,logdir, &
                           conductinibin,update_allstations,update_nostations,subincrit)
!Test setup
      WRITE(filename,'(a)') TRIM(logdir)//'tests_'//logseq(1:3)//'_'//logdate(3:8)//'_'//logtime(1:8)//'.log'
      CALL setup_for_hype_tests(fname=filename,onoff=indatacheckonoff,level=indatachecklevel)
      
      IF(status/=0) CALL stop_simulation_and_finalize_tests(1, 'loading coded info')
      nmapperiod = 0 !initialization needed if no mapoutput

!>Read input data (fot the first landuse period)
      CALL load_basindata(modeldir,forcingdir,basemodel%nsub,status)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, 'loading basin data')
      CALL load_cropdata(modeldir,ncrop,status)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, 'loading crop data')
      CALL load_pointsourcedata(modeldir,'PointSourceData.txt',basemodel%nsub,status) 
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, 'loading pointsource data')
      CALL load_permanent_soilleakage(forcingdir,basemodel%nsub,bdate,status) 
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, 'loading soil leakage')
      CALL load_branchdata(modeldir,status)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, 'loading branch data')
      CALL load_aquiferdata(modeldir,basemodel%nsub,naquifers,status) 
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, 'loading aquifer data')
      CALL load_glacierdata(modeldir,basemodel%nsub,status) 
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, 'loading glacier data')
      CALL initiate_model_parameters(basemodel%nsub,status)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, 'initiate model params')

!Test setup
!      WRITE(filename,'(a)') TRIM(logdir)//'tests_'//logseq(1:3)//'_'//logdate(3:8)//'_'//logtime(1:8)//'.log'
!      CALL setup_for_hype_tests(fname=filename,onoff=indatacheckonoff,level=indatachecklevel)
!Test obeservation data
      CALL run_hype_observation_tests(status)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(2, 'checked observations')

      CALL load_submodel_information(infodir,simsubmodel,basemodel,nsub,status)   !allocation and initialisation of basemodel
      IF(status/=0) CALL stop_simulation_and_finalize_tests(1, 'loading submodel information')
      CALL load_output_regions(modeldir,noutreg,status)    !Prepare for outregions, read Outregions.txt
      IF(status/=0) CALL stop_simulation_and_finalize_tests(1, 'loading output regions')
      
      CALL load_observations(forcingdir,basemodel%nsub,bdate,sdate,ndt,status)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, 'loading observations')
      
      CALL load_otest(modeldir,'otest.txt',status)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, "loading otest")

      CALL load_parameters(modeldir,basemodel,status,seqflag=parseq)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, "loading parameters")
      CALL finish_lakedata_table(basemodel%nsub) !Set parameter values for missing lakes with general values
      CALL finish_atmdep_file_deposition(modeldir)
      
!Read base model defined input data and set base model configuration
      CALL set_model_base_configuration(basemodel,modeldir,forcingdir,statesize,status)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, "setting model base conf")
      
      IF(simsubmodel)THEN
        CALL reform_inputdata_for_submodel(basemodel,nsub)
      ENDIF
      CALL calculate_path(nsub)

!Read model defined input data and set model configuration
      CALL load_modeldefined_input(modeldir,forcingdir,basemodel%nsub,nsub,basemodel%iindex,bdate,sdate,conductxoms,conductregest,status)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, "loading model defined input")
      CALL set_model_configuration(conduct,simulate)
      CALL set_modelconfig_from_parameters()
      
!Test input data and validate options
      CALL run_hype_tests(status)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(2, "input testing")

      CALL prepare_for_update(modeldir,update_allstations,update_nostations,basemodel%nsub,nsub,status)
      IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, "loading update input")
      CALL finalize_hype_tests()

!>Initialisations for memory allocation (states and output)
      CALL allocate_model_states(basemodel%nsub,statesize,conduct, &
              frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

      CALL set_outvar(noutput,noutvar,noutvarclass)
      IF(.NOT.doassimilation)THEN
        CALL allocate_outvar(noutvar,noutvarclass)
        CALL reallocate_outvar_information(noutvar,noutvarclass)
        CALL set_outvar_test_information(noutvar)
        CALL allocate_accumulation(nsub,nclass,numsubstances,max_classoutvar,max_basinoutvar)
      ENDIF

!Allocate local variables
      nsubCrit = nsub
      ALLOCATE(basincrit(nsubCrit,maxsubass,nacrit))
      ALLOCATE(simperformance(maxperf,nacrit))

!Preparations for subbasin output
      CALL prepare_subbasin_output(subincrit,status)
      IF(status/=0) STOP 1
      
!For data assimilation simulation, skip optimization and simulation code
      IF(.NOT.doassimilation)THEN

        CALL DATE_AND_TIME (values=datim)
        WRITE(6,*)
        WRITE(6,*) '---------------------------------------------------'
        WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
                   ' Initialisations finished, calculations starts: ',datim(1),'-',datim(2),'-',datim(3),    &
                   '  time: ',datim(5),':',datim(6),':',datim(7)
        WRITE(6,*) '---------------------------------------------------'
        FLUSH 6


!>Optimization
        IF(doopt)THEN
          CALL load_optpar(infodir)       !Reads optpar.txt and set all numerical optimization variables accordingly

          IF(optim%task_MC)THEN
            WRITE(6,*) 'MonteCarlo simulation with',optim%nruns_MC,'runs'
            CALL MonteCarlo_simulation(resdir,optim%task_writeall,frozenstate,  &
                 soilstate,aquiferstate,riverstate,lakestate,miscstate,npar)
            CALL set_optim_modpar(npar,npar,bestMCparameters(1,:))    !save the best parameters to modpar
            CALL save_respar(resdir,npar,nsub)                       !and to file
          ENDIF

          IF(optim%task_boundps)THEN
            WRITE(6,*) 'Reduce bounds of parameter space for continued MonteCarlo simulation'
            CALL bounded_MonteCarlo_simulation(optim%task_MC,frozenstate,soilstate, &
                 aquiferstate,riverstate,lakestate,miscstate,npar)
            CALL set_optim_modpar(npar,npar,bestMCparameters(1,:))    !save the best parameters to modpar
            CALL save_respar(resdir,npar,nsub)                       !and to file
          ENDIF

          lineSearchCallCount = 0
          optimFuncCall = 0
          CALL CPU_TIME(optimStartTime)

          IF(optim%task_Scanning) CALL param_scanning(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

          IF(optim%task_stageMC)THEN
            CALL stage_MonteCarlo(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
            CALL set_optim_modpar(numoptimpar,numoptimpar,bestMCparameters(1,:))    !save the best to modpar
            CALL save_respar(resdir,numoptimpar,nsub)                               !and file
          ENDIF
        
          IF(optim%task_BrentNew .OR. optim%task_stpstDesc .OR. optim%task_DFP .OR. optim%task_BFGS)THEN
            ALLOCATE(par(numoptimpar))
            CALL linesearch_methods_calibration(numoptimpar,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,par)
            CALL set_optim_modpar(numoptimpar,numoptimpar,par)    !save the best to modpar
            CALL save_respar(resdir,numoptimpar,nsub)                            !and file
          ENDIF
          IF(optim%task_DEMC)THEN
            WRITE(6,*) 'DEMC Differential-Evolution Markov Chain, with',optim%DEMC_npop, &
            'populations, and ',optim%DEMC_ngen,'generations (',optim%DEMC_ngen-1,  &
            'evolution steps)',', in total',optim%DEMC_npop*(optim%DEMC_ngen),'simulations.'
          
            CALL DEMC_simulation(resdir,optim%task_writeall,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,npar)
            CALL set_optim_modpar(npar,npar,bestMCparameters(1,:))         !save the MEDIAN parameters to modpar (1 row in bestMC)
            CALL save_respar(resdir,npar,nsub)        
          ENDIF
          IF(optim%task_runens)THEN       !write best ensemble resultat to file
            CALL save_ensemble_simulations(resdir,numoptimpar,maxperf,nacrit, &
                 optim%nruns_best,bestMCoptcrit,bestMCperformance,bestMCparameters)
          ENDIF
        ENDIF

        IF(doens)THEN
          CALL load_optpar(infodir)       !Reads optpar.txt and set all numerical optimization variables accordingly
          CALL count_ensemble_simulations(optim%task_ensbest,optim%task_ensall,optim%nruns_simloop) !optim%nruns_simloop set
          WRITE(6,*) 'Parameter ensemble simulation with',optim%nruns_simloop,'simulations.'
          IF(optim%task_ensall.OR.optim%task_ensbest)THEN
            CALL set_optpar_index()
            IF(.NOT.ALLOCATED(ensemble_parameters)) ALLOCATE(ensemble_parameters(optim%nruns_simloop,numoptimpar))
            IF(optim%task_ensall) CALL load_ensemble_simulations(modeldir,filename_MC,numoptimpar,maxperf,nacrit,&
                                          optim%nruns_simloop,ensemble_parameters)
            IF(optim%task_ensbest) CALL load_ensemble_simulations(modeldir,filename_best,numoptimpar,maxperf,nacrit,&
                                          optim%nruns_simloop,ensemble_parameters)
          !ELSE !Use par_nnn.txt files, read in ensemble loop
          ENDIF
        ENDIF
        
        IF(.NOT.optim%task_writesim)THEN
!Ensemble loop, simulate all ensemble members
          ensemble_loop:    &
       &  DO iens = 1, optim%nruns_simloop

!>Simulation start
            !>Initial model calculations; initial states, parameters, output
            ! If output for class is not asked, do not change; otherwise should be reset with changing slc info according to year.
            CALL prepare_outputfiles(resdir,nsub,naquifers,iens,optim%task_runens,optim%task_ensall)
            IF(optim%task_runens)THEN
              IF(doopt)THEN
                CALL set_optim_modpar(numoptimpar,numoptimpar,bestMCparameters(iens,:))   !set model parameters
              ELSEIF(doens .AND. (optim%task_ensall.OR.optim%task_ensbest))THEN
                CALL set_optim_modpar(numoptimpar,numoptimpar,ensemble_parameters(iens,:))   !set model parameters
              ELSEIF(doens)THEN
                CALL load_parameters(modeldir,basemodel,status,iens=iens)
                IF(status/=0) STOP 1
              ENDIF
              CALL finish_lakedata_table(nsub)
              CALL finish_atmdep_file_deposition(modeldir)
            ENDIF
            CALL initiate_output_routines()     !All output accumulation variables zeroed
            CALL set_modelconfig_from_parameters()
            IF(simsubmodel)THEN
              CALL initiate_state_for_submodel(forcingdir,basemodel,  &
                  frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
            ELSE
              IF(instate%fromfile) THEN
                CALL load_saved_state(forcingdir,nsub,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
                WRITE(6,*)
                WRITE(6,*) 'Loading saved state.'
                IF(resetstate) CALL save_soil_state_file(forcingdir,nsub,soilstate,miscstate)
              ELSE
                CALL initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
              ENDIF
              CALL initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 
            ENDIF
            oldyear = 0

!>Time Loop:
            time_loop:    &
       &    DO idt = 1,ndt
              CALL get_current_forcing(idt,nsub,d)
              ! refresh slc
              ! read recent year
              CALL calculate_time_for_model(idt,d)  !sets current_time (to d)
              IF(current_time%date%year>oldyear .AND. oldyear/=0) THEN
                write(current_year,"(I4)") current_time%date%year
                filename_year=TRIM(modeldir)//'GeoClass_'//trim(current_year)//'.txt'
                INQUIRE(file=filename_year,exist=alive)
                IF (alive) THEN
                   CALL reload_basindata(modeldir,basemodel%nsub,current_time%date%year,status) 
                   IF(status.NE.0) CALL stop_simulation_and_finalize_tests(1, 'reloading basin data')              
                   WRITE(6,*) 'Geodata and GeoClass change to year:',current_time%date%year
                   CALL inherit_state_to_new_class(modeldir,current_time%date%year,nsub,nclass,frozenstate,soilstate,miscstate)
                   WRITE(6,*) 'State inherit success'
                   CALL reallocate_nclass_based_variables()
                   WRITE(6,*) 'reallocate variables success'
                ENDIF
              END IF
              !>\li Get current input data    
              CALL log_progress(oldyear,current_time%date%year)          
              IF(ALLOCATED(psdates)) CALL get_current_pointsources(modeldir,'PointSourceData.txt',nsub,current_time%date,status) 
              IF(status.NE.0) STOP 1
              CALL get_current_soilleakage(forcingdir,nsub,current_time%date,status)
              IF(status.NE.0) STOP 1
              IF(resetstate)THEN    !Reset soil state
                DO ivar = 1,numreadstates
                  IF(current_time%date.EQ.readstatedate(ivar)) THEN
                    CALL reset_soil_state(forcingdir,nsub,soilstate,miscstate)
                    CALL format_date(readstatedate(ivar),'yyyy-mm-dd HH:MM',strdate)
                    WRITE(6,*) 'Resetting soil states '//strdate
                  ENDIF
                ENDDO
              ENDIF
              DO ivar=1,10
                IF(readadddate(ivar)==d)THEN
                  CALL reload_deposition(fileunit_temp,modeldir,basemodel,nsub,status,d)
                  IF(status/=0) STOP 1
                ENDIF
              ENDDO

              CALL initiate_outvar(idt)
              !>\li Calculate flows and update states
              CALL model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
        
              IF(status_write_outstate(current_time%date,outstate)) &    !Write state output
                CALL finalize_outstate(resdir,nsub,AddDates(current_time%date,steplen),  &
                frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 

              CALL revise_outvar(idt)    !calculate regional outvar and upstream and classes
              CALL prepare_to_compute_crit(current_time,idt,ndt)

              !>\li Calculate and write time dependent output
              DO iout = 1,noutput
                IF(output(iout)%fileformat==6) CALL write_subbasinfiles_class(iout,idt,ndt,current_time)
                IF(output(iout)%fileformat==5) CALL write_timefiles_class(iout,idt,ndt,current_time)
                IF(output(iout)%fileformat==4) CALL write_regionfiles(iout,idt,ndt,current_time)
                IF(output(iout)%fileformat==1) CALL write_subbasinfiles(iout,idt,ndt,current_time)
                IF(output(iout)%fileformat==3) CALL write_timefiles(iout,idt,ndt,current_time)
                IF(output(iout)%fileformat==2) CALL compute_mapoutput(current_time,iout,idt,ndt,writematlab,nmapperiod)  !Save data for map output
              ENDDO
              IF(writeload)THEN
                CALL compute_outloads(current_time%date,pwrite,idt)       !Write yearly load for all subbasins
                IF(pwrite) CALL save_loadfiles(resdir,current_time%date%year)
              ENDIF

            ENDDO time_loop

!>Compute and write criteria
            IF(nacrit/=0) THEN
              CALL calculate_criteria(optcrit,basincrit,simperformance,condcrit,condthres)
              CALL write_simulation_assessment(resdir,iens,nacrit,optcrit,     &
                  simperformance,optim%task_runens,condcrit,condthres)
              CALL write_subbasin_assessment(resdir,nsubCrit,nacrit,basincrit,iens,optim%task_runens)
            ENDIF

!Save and close files or prepare them for next ensemble member simulation
            CALL close_outputfiles(nsub,naquifers,iens)
            IF(iens==optim%nruns_simloop)THEN
              CALL close_observations(forcingdir)
            ELSE
              CALL reset_observations(forcingdir,status)
              IF(status.NE.0) STOP 1
            ENDIF

!>Write results to files
            CALL save_mapfiles(resdir,nmapperiod,iens,optim%task_runens,optim%task_ensall)

          ENDDO ensemble_loop   !simulate ensembles
        ENDIF !.NOT.optim%task_writesim
      ENDIF !.NOT.doassimilate

  
#ifdef _ASSIMILATION      
!===========================================================================
!Data assimilation version of simulation; time and ensemble loops
!===========================================================================
      IF(doAssimilation)THEN
        WRITE(6,*)
        WRITE(6,*) '---------------------------------------------------'
        WRITE(6,*) 'HYPE Data Assimilation Simulation'
        WRITE(6,*) '---------------------------------------------------'
        WRITE(6,*)
        WRITE(6,*) ' [Data Assimilation] Initialization started...'
        WRITE(6,*)

!HYSS: Initializations of hyss variables
!Note. The counter, "iens", is used in the standard code above for ensemble simulations and is set to 1 here as if this was a normal determinstic simulation.
!In the assimilation section below, we use a dedicated counter called "iassim" for looping over the assimilation ensemble members.
        iens=1      
        
!HYSS: Initial model calculations; initial states, parameters
        CALL set_modelconfig_from_parameters()
        IF(simsubmodel)THEN
          CALL initiate_state_for_submodel(forcingdir,basemodel,  &
              frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
        ELSE
          IF(instate%fromfile) THEN
            CALL load_saved_state(forcingdir,nsub,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
            WRITE(6,*)
            WRITE(6,*) 'Loading saved state.'
          ELSE
            CALL initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
          ENDIF
          CALL initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 
        ENDIF
        oldyear = 0
        
!ASSIMILATION: Initialization of Data Assimilation variables
        CALL set_stateinfo(nsub,statesize,stateinfo, &
              frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
        CALL assim_initialize(forcingdir,modeldir,resdir,conductinibin,myAssimData,assim_ens_size,stateinfo,noutvar,status)
        IF(status/=0) STOP 1
        
!Prepare files for output
!The standard output files are the mean (or median) of the ensemble data.
!Additional output in sequence named files are optional; output files with
!ensemble statistics (min,max,quantiles) or all ensemble members individually. 
        CALL reallocate_outvar_information(noutvar,noutvarclass)
        CALL da_allocate_accumulation(nsub,myAssimData%info%nstatout)
        DO iens = 2,myAssimData%info%nstatout+1           !this is statistics and ensemble members
          CALL prepare_outputfiles(resdir,nsub,naquifers,iens,.TRUE.,.FALSE.,ensstat=myAssimData%info%nstatout)
        ENDDO
        iens = 1    !reset
        CALL prepare_outputfiles(resdir,nsub,naquifers,iens,.FALSE.,.FALSE.)  !for standard output meanORmedian results
        CALL initiate_output_routines()     !All output accumulation variables zeroed
        
!ASSIMILATION: Initialization finished log message
        CALL DATE_AND_TIME (values=datim)
        WRITE(6,*)
        WRITE(6,*) '---------------------------------------------------'
        WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
                 ' [Data Assimilation] Initialisations finished, calculations starts: ',datim(1),'-',datim(2),'-',datim(3),    &
                 '  time: ',datim(5),':',datim(6),':',datim(7)
        WRITE(6,*) '---------------------------------------------------'

!ASSIMILATION: Time Loop:
        assim_time_loop: DO idt = 1,ndt

!HYSS: Get current input data
          CALL get_current_forcing(idt,nsub,d)
          CALL calculate_time_for_model(idt,d)
          CALL log_progress(oldyear,current_time%date%year)
          WRITE(6,*)' [Data Assimilation] Start of new timestep at day no:', current_time%dayno  !dayno set in calculate_time_for_model
          IF(ALLOCATED(psdates)) CALL get_current_pointsources(modeldir,'PointSourceData.txt',nsub,current_time%date,status) 
          IF(status.NE.0) STOP 1
          CALL get_current_soilleakage(forcingdir,nsub,current_time%date,status)
          IF(status.NE.0) STOP 1

!ASSIMILATION: Generate input ensembles for this timestep
          CALL generate_forcing_ensemble(myAssimData)

!ASSIMILATION: Generate observation ensemble if analysis is enabled for this time step
          !IF(idt>dtskip)THEN
          !  CALL generate_observation_ensemble(myAssimData)
          !ENDIF

!ASSIMILATION: Loop over assimilation ensemble members for this timestep
          assim_ensemble_loop: DO iassim = 1,assim_ens_size

!ASSIMILATION: Write states, forcing, parameters etc, for the current (iassim) ensemble member to the HYPE model data
            CALL ensemble_to_model2(iassim,myAssimData%info%nA,myAssimData%info%nF, & 
                                    myAssimData%X,myAssimData%A,myAssimData%F,stateinfo) !,myAssimData%info%nX
                
!HYSS: Initiate the outvars for the current (iassim) ensemble member
            CALL initiate_outvar(idt)

!HYSS: Run the model for the current (iassim) ensemble member
            CALL model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
  	            
!HYSS: revise some outvars (regional outvars for instance) - DG20200212: it has to be inside ensemble loop to prepare outvars for predicted obs-ensemble
            CALL revise_outvar(idt)
            
!ASSIMILATION: Write model 'forecast' back to the ensemble data
            CALL model_to_ensemble2(iassim,myAssimData%X,myAssimData%A,myAssimData%info%nA,stateinfo,.FALSE.)
            
!ASSIMILATION: Generate observation ensemble if analysis is enabled for this time step
!ASSIMILATION: Write the predicted observations to the ensemble if analysis is enabled (idt>dtskip)
            IF(idt>dtskip)THEN
              IF(iassim==1) CALL generate_observation_ensemble(myAssimData)
              CALL modelobservations_to_ensemble(iassim,myAssimData)
            ENDIF
          
          ENDDO assim_ensemble_loop

!ASSIMILATION, Ensemble Kalman Filter Analysis (later generalize to any DA filter method) - ensemble state matrix will be updated if ensemble size>1 enkf analysis enabled (idt>dtskip) 
          IF(idt>dtskip.AND.myAssimData%info%nE.GT.1)THEN
            CALL enkf_analysis_main(myAssimData)
          ENDIF

!ASSIMILATION, Update ensemble statistics
          CALL updateEnsembleStatistics(myAssimData)
            
!ASSIMILATION: Write ensemble mean (or median) to model variables, to be included in the standard output files.
          CALL meanORmedian_to_model2(myAssimData%info%nA,myAssimData%X,myAssimData%A,myAssimData%info%meanout,stateinfo)
            
!ASSIMILATION: Reset selected assimilation ensembles to ensemble mean (all variables not included in the "control vector" will be reset to ensemble mean after each timestep) 
          IF(myAssimData%info%collapseNonControlled) &
            CALL meanORmedian_to_ensemble(myAssimData%info%nX,myAssimData%X,myAssimData%info%meanout)
            
!HYSS: Write state output files
          IF(status_write_outstate(current_time%date,outstate)) &
            CALL finalize_outstate(resdir,nsub,AddDates(current_time%date,steplen),  &
            frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 
!ASSIMILATION: Write state ensemble output files for later restart
          CALL finalize_ensemble_states_outstate(resdir,nsub,current_time%date,myAssimData)

!HYSS: Preparatory calculations for criteria (update various sums, square sums and counters with current timestep data) 
!DG20200212: revise_outvar has to be inside ensemble-loop          CALL revise_outvar(idt)
          CALL prepare_to_compute_crit(current_time,idt,ndt)
	
!HYSS: Calculate and write time dependent standard output
          DO iout = 1,noutput
            IF(output(iout)%fileformat==4) CALL write_regionfiles(iout,idt,ndt,current_time)
            IF(output(iout)%fileformat==1) CALL write_subbasinfiles(iout,idt,ndt,current_time)
            IF(output(iout)%fileformat==3) CALL write_timefiles(iout,idt,ndt,current_time)
            IF(output(iout)%fileformat==5) CALL write_timefiles_class(iout,idt,ndt,current_time)
            IF(output(iout)%fileformat==6) CALL write_subbasinfiles_class(iout,idt,ndt,current_time)
            IF(output(iout)%fileformat==2) CALL compute_mapoutput(current_time,iout,idt,ndt,writematlab,nmapperiod)
          ENDDO
          IF(writeload)THEN
            CALL compute_outloads(current_time%date,pwrite,idt)       !Write yearly load total for all subbasins
            IF(pwrite) CALL save_loadfiles(resdir,current_time%date%year)
          ENDIF
          
!ASSIMILATION: Output of "statistical" and ensemble member simulations result
          IF(idt>dtskip)THEN
            DO iens = 2,myAssimData%info%nstatout+1
              CALL statistics_to_modeloutput(myAssimData%info%nA,myAssimData%A,iens) !,myAssimData%info%nE
              DO iout = 1,noutput
                IF(output(iout)%fileformat==1) CALL write_subbasinfiles_in_parallel(iout,idt,ndt,iens,current_time)
                IF(output(iout)%fileformat==3) CALL write_timefiles_in_parallel(iout,idt,ndt,iens,current_time)
                IF(output(iout)%fileformat==5) CALL write_timefiles_class_in_parallel(iout,idt,ndt,iens,current_time)
                IF(output(iout)%fileformat==4) CALL write_regionfiles_in_parallel(iout,idt,ndt,iens,current_time)
                IF(output(iout)%fileformat==6) CALL write_subbasinfiles_class_in_parallel(iout,idt,ndt,iens,current_time)
              ENDDO
            ENDDO
            iens = 1
            CALL statistics_to_modeloutput(myAssimData%info%nA,myAssimData%A,iens,myAssimData%info%meanout) !,myAssimData%info%nE
          ENDIF

!ASSIMILATION: End of assimilation time loop
        ENDDO assim_time_loop

!HYSS: Compute and write criteria (nb! iens should still be set to 1) 
        IF(nacrit/=0) THEN
          CALL calculate_criteria(optcrit,basincrit,simperformance,condcrit,condthres)
          CALL write_simulation_assessment(resdir,iens,nacrit,optcrit,     &
                simperformance,.FALSE.,condcrit,condthres)
          CALL write_subbasin_assessment(resdir,nsubCrit,nacrit,basincrit,iens,.FALSE.)
        ENDIF

!HYSS: Save and close files
        DO iens = 1,myAssimData%info%nstatout+1
          CALL close_outputfiles(nsub,naquifers,iens,ensstat=myAssimData%info%nstatout)
        ENDDO
        CALL close_observations(forcingdir)
        iens = 1    !reset
        
!HYSS: Write results to map-files
        CALL save_mapfiles(resdir,nmapperiod,iens,.FALSE.,.FALSE.)

!ASSIMILATION: (some deallocation), closing of files, (cleaning up, and final log messages to be added?)
        CALL finalize_assimilation(myAssimData)
        
      ENDIF !doassimilation
!=====================================================================
!End of Data Assimilation section
!=====================================================================
#endif
  
  
!>Deallocate variables
      CALL deallocate_worldvar()
      CALL deallocate_modvar(nsub)
      CALL deallocate_model_states(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
      IF(doopt.AND.optim%task_MC) CALL deallocate_MCvariables()
      IF(ALLOCATED(ensemble_parameters)) DEALLOCATE(ensemble_parameters)
      IF(ALLOCATED(par)) DEALLOCATE(par)
      IF(ALLOCATED(simperformance)) DEALLOCATE(simperformance)
      IF(ALLOCATED(basincrit)) DEALLOCATE(basincrit)
      IF(ALLOCATED(subincrit)) DEALLOCATE(subincrit)

!Write stop time on log-file
      CALL DATE_AND_TIME (values=datim)
      WRITE(6,*)
      WRITE(6,*) '---------------------------------------------------'
      WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
               ' Job finished date: ',datim(1),'-',datim(2),'-',datim(3),    &
               '  time: ',datim(5),':',datim(6),':',datim(7)

      CLOSE(6)
      IF(usestop84) STOP 84
      !CALL stop_simulation_and_finalize_tests(0, "successfully") !CP tests already finalized above
      WRITE(0,*) '-> Simulation will halt: successfully'
      STOP 0
      
      !Error handling
900   CONTINUE
      WRITE(6,*) 'ERROR: Open file: '//TRIM(filename)//' for writing.'
      STOP 1
      
      END PROGRAM
