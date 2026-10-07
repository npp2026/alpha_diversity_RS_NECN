run_pipeline_until <- function(target_stage = "validation", config_path = Sys.getenv("BREAK_CONFIG", "config/v5_6.yml"), project_dir = getwd()) {
  stages <- c("primary","sensitivity","uncertainty","spatial","tables","validation")
  if(!target_stage %in% stages)stop("Unknown pipeline stage: ",target_stage,call.=FALSE)
  cfg <- read_config(config_path, project_dir); assert_packages(c("terra","yaml","digest","jsonlite")); ctx <- create_run_context(cfg); write_resolved_config(cfg,ctx)
  inventory<-rho_tab<-rho_by_response<-primary<-boundary<-rho_sens<-multiplicity<-pixel<-trajectory<-spatial<-NULL
  stage_reached <- "setup"
  tryCatch({
    inventory <- build_input_inventory(cfg); write_csv(inventory$manifest,file.path(ctx$root,"config","input_manifest.csv")); log_msg("Input fingerprint: ",inventory$fingerprint,.file=ctx$log)
    if(!isTRUE(cfg$rho_estimation$enabled))stop("v5.6 currently requires rho_estimation.enabled=true for empirical primary calibration.",call.=FALSE)
    rho_tab <- estimate_rho_all(cfg,inventory,ctx); rho_by_response <- rho_list_from_table(rho_tab)
    m <- as.integer(cfg$analysis$primary$min_segment); primary <- run_scenario(cfg,inventory,rho_by_response,m,"PRIMARY",file.path(ctx$root,"primary"),ctx)
    stage_reached <- "primary"; write_run_manifest(cfg,inventory,rho_tab,primary,ctx,"running",stage=stage_reached)
    if(target_stage=="primary"){write_run_manifest(cfg,inventory,rho_tab,primary,ctx,"completed_partial",stage=stage_reached);return(list(cfg=cfg,ctx=ctx,inventory=inventory,rho=rho_tab,primary=primary))}

    boundary <- run_boundary_sensitivity(cfg,inventory,rho_by_response,ctx,primary)
    rho_sens <- run_rho_sensitivity(cfg,inventory,rho_by_response,ctx,primary)
    multiplicity <- build_multiplicity_summary(cfg,primary,ctx); stage_reached <- "sensitivity"
    if(target_stage=="sensitivity"){write_run_manifest(cfg,inventory,rho_tab,primary,ctx,"completed_partial",stage=stage_reached);return(list(cfg=cfg,ctx=ctx,inventory=inventory,rho=rho_tab,primary=primary,boundary=boundary,rho_sensitivity=rho_sens,multiplicity=multiplicity))}

    pixel <- run_pixel_bootstrap(cfg,inventory,primary,rho_by_response,ctx)
    trajectory <- run_regional_trajectory_analysis(cfg,inventory,ctx); stage_reached <- "uncertainty"
    if(target_stage=="uncertainty"){write_run_manifest(cfg,inventory,rho_tab,primary,ctx,"completed_partial",stage=stage_reached);return(list(cfg=cfg,ctx=ctx,inventory=inventory,rho=rho_tab,primary=primary,boundary=boundary,rho_sensitivity=rho_sens,multiplicity=multiplicity,pixel=pixel,trajectory=trajectory))}

    spatial <- run_spatial_validation(cfg,primary,boundary,ctx); stage_reached <- "spatial"
    if(target_stage=="spatial"){write_run_manifest(cfg,inventory,rho_tab,primary,ctx,"completed_partial",stage=stage_reached);return(list(cfg=cfg,ctx=ctx,inventory=inventory,rho=rho_tab,primary=primary,boundary=boundary,rho_sensitivity=rho_sens,multiplicity=multiplicity,pixel=pixel,trajectory=trajectory,spatial=spatial))}

    build_publication_tables(cfg,primary,boundary,multiplicity,rho_tab,pixel,trajectory,spatial,ctx); stage_reached <- "tables"
    if(target_stage=="tables"){write_run_manifest(cfg,inventory,rho_tab,primary,ctx,"completed_partial",stage=stage_reached);return(list(cfg=cfg,ctx=ctx,primary=primary))}

    if(isTRUE(cfg$validation$run_unit_tests)){source(file.path(project_dir,"validation","unit_tests.R"));run_unit_tests_v56(project_dir)}
    if(isTRUE(cfg$validation$run_smoke_tests)){source(file.path(project_dir,"validation","smoke_test.R"));run_smoke_test_v56(project_dir)}
    simfiles<-character()
    if(isTRUE(cfg$validation$run_simulation_null)){source(file.path(project_dir,"validation","simulation_null.R"));z<-run_null_simulation_v56(project_dir,B_sim=as.integer(cfg$validation$simulation_B));f<-file.path(ctx$root,"validation","null_simulation.csv");write_csv(z,f);simfiles<-c(simfiles,f)}
    if(isTRUE(cfg$validation$run_simulation_power)){source(file.path(project_dir,"validation","simulation_power.R"));z<-run_power_simulation_v56(project_dir,B_sim=as.integer(cfg$validation$simulation_B));f<-file.path(ctx$root,"validation","power_simulation.csv");write_csv(z,f);simfiles<-c(simfiles,f)}
    if(isTRUE(cfg$validation$run_bootstrap_coverage)){source(file.path(project_dir,"validation","bootstrap_coverage.R"));z<-run_bootstrap_coverage_v56(project_dir,B_outer=as.integer(cfg$validation$bootstrap_coverage_outer),B_boot=as.integer(cfg$validation$bootstrap_coverage_inner));f<-file.path(ctx$root,"validation","bootstrap_coverage.csv");write_csv(z,f);simfiles<-c(simfiles,f)}
    build_validation_report(cfg,inventory,rho_tab,primary,boundary,multiplicity,trajectory,spatial,ctx,simfiles)
    file.copy(file.path(project_dir,"validation","risk_register.csv"),file.path(ctx$root,"validation","risk_register.csv"),overwrite=TRUE)
    stage_reached <- "validation"; write_run_manifest(cfg,inventory,rho_tab,primary,ctx,"completed",stage=stage_reached)
    list(cfg=cfg,ctx=ctx,inventory=inventory,rho=rho_tab,primary=primary,boundary=boundary,rho_sensitivity=rho_sens,multiplicity=multiplicity,pixel=pixel,trajectory=trajectory,spatial=spatial)
  }, error=function(e){
    write_run_manifest(cfg,inventory,rho_tab,primary,ctx,"failed",error_message=conditionMessage(e),stage=stage_reached)
    stop(e)
  })
}
