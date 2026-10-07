write_run_manifest <- function(cfg, inventory = NULL, rho_tab = NULL, primary = NULL, ctx, status = "completed", error_message = NULL, stage = NULL) {
  pkgs <- c("yaml","terra","digest","jsonlite")
  versions <- stats::setNames(lapply(pkgs,function(p)if(requireNamespace(p,quietly=TRUE))as.character(utils::packageVersion(p)) else NA_character_),pkgs)
  x <- list(version=cfg$version,run_id=ctx$run_id,status=status,stage=stage,
            started_or_updated=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),config_hash=ctx$config_hash,code_hash=ctx$code_hash,
            input_fingerprint=if(!is.null(inventory))inventory$fingerprint else NULL,
            primary_scenario=if(!is.null(primary))primary$id else NULL,
            rho=if(!is.null(rho_tab)&&nrow(rho_tab))stats::setNames(as.list(rho_tab$rho_for_calibration),rho_tab$response)else NULL,
            supf=list(B=cfg$supf$mc_B,B_by_response=cfg$supf$mc_B_by_response %||% NULL,p_rule=cfg$supf$p_rule,mc_chunk=cfg$supf$mc_chunk,mc_rng_scheme="chunk_key_v1"),multiple_testing=cfg$multiple_testing,rng=cfg$rng,parallel=cfg$runtime$parallel %||% list(enabled=FALSE),
            package_versions=versions,R_version=R.version.string,platform=R.version$platform,error_message=error_message)
  write_json(x,file.path(ctx$root,"config","run_manifest.json"))
  if(!is.null(inventory))write_csv(inventory$manifest,file.path(ctx$root,"config","input_manifest.csv"));invisible(x)
}
