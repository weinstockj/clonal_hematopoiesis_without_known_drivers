library(targets)

# Set target options:
tar_option_set(
    packages = c(
        "tibble",
        "reticulate",
        "dplyr",
        "glue",
        "logger",
        "ggplot2",
        "torch",
        "BSgenome.Hsapiens.UCSC.hg38"
    ), # packages that your targets need to run
    format = "rds" # default storage format
      # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
# options(clustermq.scheduler = "slurm")
# options(clustermq.template = "clustermq.tmpl")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
list(
  tar_target(
    name = control_variants,
    command = read_in_control_variants(),
    format = "parquet"
  ),
  tar_target(
    name = somatic_variants,
    command = read_in_somatic_variants(),
    format = "parquet"
  ),
  tar_target(
    name = common_dups,
    command = read_common_dups()
  ),
  tar_target(
    name = hg19_diff,
    command = read_in_hg19_diff()
  ),
  tar_target(
    name = umap,
    command = read_in_umap_mappability()
  ),
  tar_target(
    name = exclude_bad_mutations,
    command = exclude_beds(
        somatic_variants,
        list(
            "hg19diff" = hg19_diff,
            "dups" = common_dups
        )
    ),
    format = "parquet"
  ),
  tar_target(
    name = filter_to_high_mappability,
    command = overlap_beds(
        exclude_bad_mutations,
        list(
            "umap" = umap
        )
    ),
    format = "parquet"
  ),
  tar_target(
    name = somatic_variants_filtered,
    command = prune_somatic_variants_by_sample_list_and_vaf(filter_to_high_mappability),
    format = "parquet"
  ),
  # tar_target(
  #   name = sample_control_variants,
  #   command = sample_variants(control_variants),
  #   format = "parquet"
  # ),
  # tar_target(
  #   name = sample_somatic_variants,
  #   command = sample_variants(somatic_variants),
  #   format = "parquet"
  # ),
  tar_target(
    name = training_sample,
    command = dplyr::bind_rows(control_variants, somatic_variants_filtered),
    format = "parquet"
  ),
  tar_target(
    name = write_to_vcf_,
    command = write_to_vcf(training_sample),
    format = "file"
  ),
  tar_target(
    name = execute_vep_,
    command = execute_vep(write_to_vcf_),
    format = "file"
  ),
  tar_target(
    name = training_annotated,
    command = read_in_training_annotated(execute_vep_, training_sample),
    format = "parquet"
  ),
  tar_target(name = cd34, command = read_in_cd34()),
  tar_target(
    name = annotate_chromatin,
    command = join_chromatin(training_annotated, cd34)
  ),
  tar_target(
    name = annotate_chromatin_sequence,
    command = attach_sequence_context(annotate_chromatin),
    format = "parquet"
  ),
  tar_target(
    name = imputed_training,
    command = impute_and_transform(annotate_chromatin_sequence),
    format = "parquet"
  ),
  tar_target(
    name = train_model,
    command = create_recipe_and_bake(imputed_training),
    format = "qs",
    resources = tar_resources(
        qs = tar_resources_qs(preset = "balanced")
    )
  ),
  tar_target(
    name = training_with_predictions,
    command = predict_training_future(train_model, imputed_training, somatic_variants_filtered),
    format = "parquet"
  ),
  tar_target(
    name = vip_plots,
    command = plot_vip(train_model$final),
    format = "file"
  ),
  tar_target(
    name = pred_grid_impact_af,
    command = create_pred_grid_impact_af(train_model)
  ),
  tar_target(
    name = pred_grid_mut_af,
    command = create_pred_grid_mut_af(train_model)
  ),
  tar_target(
    name = pred_grid_cadd_af,
    command = create_pred_grid_cadd_af(train_model)
  ),
  tar_target(
    name = pdp_plots_impact_af,
    command = plot_pdp(pred_grid_impact_af, "AF", IMPACT, AF),
    format = "file"
  ),
  tar_target(
    name = pdp_plots_mut_af,
    command = plot_pdp(pred_grid_mut_af, "MUT", MUT_TYPE, AF),
    format = "file"
  ),
  tar_target(
    name = pdp_plots_cadd_af,
    command = plot_pdp(pred_grid_cadd_af, "CADD", CADD_PHRED, AF),
    format = "file"
  ),
  tar_target(
    name = plot_auc_,
    command = plot_auc(train_model),
    format = "file"
  ),
  tar_target(
    name = sample_meta,
    command = read_in_sample_meta(),
    format = "parquet"
  ),
  tar_target(
    name = encore_input,
    command = write_encore_output(sample_meta, training_with_predictions),
    format = "file"
  ),
  tar_target(
    name = torch_input,
    command = create_data_for_torch(train_model, imputed_training, somatic_variants_filtered, sample_meta),
    format = "qs",
    resources = tar_resources(
        qs = tar_resources_qs(preset = "balanced")
    )
  ),
  tar_target(
    name = fitted_model,
    command = torch_model(
        torch_input$x[, torch_predictors()], 
        torch_input$y,
        torch_input$sample_map,
        lr = 0.08,
        iters = 200,
        15,
        2000
    )
  ),
  tar_target(
    name = plot_fitted_model_beta_,
    command = plot_fitted_model_beta(fitted_model),
    format = "file"
  ),
  tar_target(
    name = mCA_manifest,
    command = read_in_mCA_manifest(),
    format = "parquet"
  ),
  tar_target(
    name = write_encore_input_semi_supervised_,
    command = write_semi_supervised_encore_output(sample_meta, fitted_model, torch_input, mCA_manifest),
    format = "file"
  )
)
