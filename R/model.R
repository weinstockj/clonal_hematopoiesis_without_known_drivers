create_recipe_and_bake = function(training, n_sample = 1e5) {

    set.seed(1)

    training = dplyr::sample_n(training, size = n_sample)
    train_test_split = rsample::initial_split(training)

    other = rsample::training(train_test_split) # other = train + validation

    test = rsample::testing(train_test_split)

    val_set = rsample::vfold_cv(other)

    logger::log_info("now building recipe")

    formula = as.formula(
                glue::glue(
                "FILTER ~ IMPACT + DISTANCE + SOMATIC + HIGH_INF_POS + CADD_PHRED + MUT_TYPE + ",
                # "FILTER ~ CADD_PHRED ",
                # "FILTER ~ MPACT + CADD_PHRED ",
                paste(glue::glue("seq_context_{1:36}"), collapse = " + "),
                " + VAF + AF + state" 
                # " + VAF " 
                )
            )

    rec = recipes::recipe(
            formula,
            data = other
    ) %>%
        recipes::step_dummy(recipes::all_nominal_predictors(), one_hot = TRUE)

    # rec = recipes::prep(rec, other, retain = TRUE)
    # other_baked = recipes::bake(rec, new_data = NULL)
    # test_baked = recipes::bake(rec, new_data = test)

    # browser()

    # model = parsnip::mlp(
    #             hidden_units = 5L,
    #             penalty = tune::tune(),
    #             # dropout = tune::tune(),
    #             # learn_rate = tune::tune(),
    #             epochs = 10,
    #             # mixture = 0.5,
    #             activation = "relu"
    #     ) %>%
    #     # parsnip::set_engine("brulee") %>%
    #     parsnip::set_engine("keras") %>%
    #     parsnip::set_mode("classification")

    # grid = dials::grid_regular(
    #     dials::penalty(),
    #     # dials::dropout(),
    #     # dials::learn_rate(),
    #     levels = 2
    # )

    model <- parsnip::boost_tree(
                    trees = 1000, 
                    tree_depth = tune::tune(), min_n = tune::tune(), 
                    loss_reduction = tune::tune(),                     ## first three: model complexity
                    sample_size = tune::tune(), mtry = tune::tune(),         ## randomness
                    learn_rate = tune::tune(),                         ## step size
                ) %>% 
                parsnip::set_engine("xgboost") %>% 
                parsnip::set_mode("classification")

    grid <- dials::grid_latin_hypercube(
                dials::tree_depth(),
                dials::min_n(),
                dials::loss_reduction(),
                sample_size = dials::sample_prop(),
                dials::finalize(dials::mtry(), training),
                dials::learn_rate(),
                size = 10
        )

    logger::log_info("now building workflow")

    model_workflow = workflows::workflow() %>%
        workflows::add_model(model) %>%
        workflows::add_recipe(rec) 

    logger::log_info("now training")


    num_workers = 20L
    cl <- parallel::makePSOCKcluster(num_workers)
    doParallel::registerDoParallel(cl)

    workflow_result = model_workflow %>%
        # parsnip::fit(training)
        tune::tune_grid(
            resamples = val_set,
            grid = grid
        )

    parallel::stopCluster(cl)
    foreach::registerDoSEQ()

    logger::log_info("done training workflow")

    # train_baked = recipes::bake(trained_rec, new_data = other)
    # test_baked = recipes::bake(trained_rec, new_data = test)

    best = workflow_result %>%
        tune::select_best(metric = "roc_auc")

    logger::log_info("now finalizing workflow")

    final = model_workflow %>%
        tune::finalize_workflow(best) %>%
        parsnip::fit(other)

    logger::log_info("done finalizing workflow")

    trained_rec = recipes::prep(rec, data = other)

    return(
        list(
            # "train" = train_baked,
            # "test" = test_baked,
            "other" = other,
            "test" = test,
            "recipe" = trained_rec,
            "workflow" = model_workflow,
            "result" = workflow_result,
            "final" = final
        )
    )
}

figure_dir = function() {
    "output/figures"
}

plot_vip = function(model) {

    plot = workflows::extract_fit_parsnip(model) %>%
        vip::vip(.) +
        cowplot::theme_cowplot(font_size = 13)

    fname = file.path(figure_dir(), "vip_xgboost.pdf")

    ggplot2::ggsave(
        filename = fname,
        plot = plot,
        width = 6,
        height = 4,
        units = "in"
    )

    return(fname)
}

plot_auc = function(train) {
    
    fnames = c(
        "roc_somatic" = file.path(figure_dir(), "somatic_auc.pdf"),
        "pr_somatic" = file.path(figure_dir(), "somatic_pr.pdf"),
        "confusion" = file.path(figure_dir(), "confusion_mat.pdf")
    )

    preds_probs = workflows:::predict.workflow(
        train$final,
        new_data = train$test,
        type = "prob"
    )

    preds_class = workflows:::predict.workflow(
        train$final,
        new_data = train$test
    )

    preds_probs$FILTER = as.factor(train$test$FILTER) 
    preds_class$FILTER = as.factor(train$test$FILTER) 


    plot = preds_probs %>%
        yardstick::roc_curve(FILTER, .pred_artifact, .pred_germline, .pred_somatic) %>%
        autoplot(.) +
        cowplot::theme_cowplot(font_size = 11)

    ggsave(fnames["roc_somatic"], plot, width = 6, height = 4, units = "in")

    plot = preds_probs %>%
        yardstick::pr_curve(FILTER, .pred_artifact, .pred_germline, .pred_somatic) %>%
        autoplot(.) +
        cowplot::theme_cowplot(font_size = 11)

    ggsave(fnames["pr_somatic"], plot, width = 6, height = 4, units = "in")

    plot = preds_class %>%
        yardstick::conf_mat(FILTER, .pred_class) %>%
        autoplot(type = "heatmap") +
        cowplot::theme_cowplot(font_size = 11)

    ggsave(fnames["confusion"], plot, width = 6, height = 4, units = "in")

    return(fnames)
}

create_pred_grid_impact_af = function(train) {

    base_df = train$other %>%
        dplyr::slice(1)

    grid = tidyr::expand_grid(
        VAF = seq(.03, .7, length = 20),
        AF = c(1e-6, 1e-4, 1e-3, 1e-2, .3),
        state = c("Quies"),
        # state = c("Quies", "TxWk"),
        IMPACT = c("missing", "MODIFIER")
    )

    grids = purrr::pmap_dfr(
        list(grid$VAF, grid$AF, grid$state, grid$IMPACT),
        function(a, b, c, d) {
            row = base_df
            row$raw_VAF = a
            row$VAF = (a - .37) / .23
            row$IMPACT = d
            row$AF = b
            row$state = c
            return(row)
        }
    )

    preds = workflows:::predict.workflow(train$final, new_data = grids, type = "prob")

    return(dplyr::bind_cols(grids, preds))
}

create_pred_grid_cadd_af = function(train) {

    base_df = train$other %>%
        dplyr::slice(1)

    grid = tidyr::expand_grid(
        VAF = seq(.03, .7, length = 20),
        AF = c(1e-6, 1e-4, 1e-3, 1e-2, .3),
        state = "Quies",
        IMPACT = c("MODIFIER"),
        CADD_PHRED = c(0, 10, 20),
        MUT_TYPE = c("T.C")
    )

    grids = purrr::pmap_dfr(
        list(grid$VAF, grid$AF, grid$state, grid$IMPACT, grid$MUT_TYPE, grid$CADD_PHRED),
        function(a, b, c, d, e, f) {
            row = base_df
            row$raw_VAF = a
            row$VAF = (a - .37) / .23
            row$IMPACT = d
            row$MUT_TYPE = e
            row$AF = b
            row$CADD_PHRED = f
            row$state = c
            return(row)
        }
    )

    preds = workflows:::predict.workflow(train$final, new_data = grids, type = "prob")

    return(dplyr::bind_cols(grids, preds))

}

create_pred_grid_mut_af = function(train) {

    base_df = train$other %>%
        dplyr::slice(1)

    grid = tidyr::expand_grid(
        VAF = seq(.03, .7, length = 20),
        AF = c(1e-6, 1e-4, 1e-3, 1e-2, .3),
        state = "Quies",
        IMPACT = c("MODIFIER"),
        MUT_TYPE = c("C.T", "G.A")
    )

    grids = purrr::pmap_dfr(
        list(grid$VAF, grid$AF, grid$state, grid$IMPACT, grid$MUT_TYPE),
        function(a, b, c, d, e) {
            row = base_df
            row$raw_VAF = a
            row$VAF = (a - .37) / .23
            row$IMPACT = d
            row$MUT_TYPE = e
            row$AF = b
            row$state = c
            return(row)
        }
    )

    preds = workflows:::predict.workflow(train$final, new_data = grids, type = "prob")

    return(dplyr::bind_cols(grids, preds))

}

plot_grid = function(grid, ...) {

    grid %>%
        tidyr::pivot_longer(names_to = "class", values_to = "probability", starts_with(".pred")) %>%
        dplyr::mutate(
            class = stringr::str_remove(class, ".pred_")
        ) %>%
        dplyr::select(-VAF) %>%
        dplyr::rename(VAF = raw_VAF) %>%
        ggplot(data = ., aes(x = VAF, y = probability, color = class)) +
            geom_point() + 
            geom_smooth() + 
            scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
            facet_wrap(vars(...), nrow = 2) +
            cowplot::theme_cowplot(font_size = 12) +
            cowplot::panel_border()

}

pdp_fnames = function(tag) {

    fnames = c(
        "VAF" = file.path(figure_dir(), "pdp_vaf_plot.pdf"),
        "AF" = file.path(figure_dir(), "pdp_af_plot.pdf"),
        "CADD" = file.path(figure_dir(), "pdp_cadd_af_plot.pdf"),
        "MUT" = file.path(figure_dir(), "pdp_mut_plot.pdf")
    )

    return(fnames[tag])
}

plot_pdp = function(grid, tag, ...) {
    
    plot = plot_grid(grid, ...)
    ggsave(pdp_fnames(tag), plot, width = 10, height = 6, units = "in")

    return(pdp_fnames(tag))
}


predict_training = function(model, training, variants) {

    meta = variants %>%
        dplyr::select(ID, NWD_ID) %>%
        dtplyr::lazy_dt(key_by = "ID")

    log_info("now generating predictions across all training data")
    
    preds_probs = workflows:::predict.workflow(
        model$final,
        new_data = training %>% 
            # dplyr::mutate(SOMATIC = as.integer(SOMATIC)) %>% 
            dplyr::filter(FILTER == "somatic"),
        type = "prob"
    ) %>%
        dplyr::bind_cols(training %>% dplyr::select(ID, CHROM, POS, REF, ALT, VAF, AF, state, CADD_PHRED)) %>%
        dtplyr::lazy_dt(key_by = "ID")

    log_info("done generating predictions across all training data")

    preds_probs = preds_probs %>%
        dplyr::inner_join(meta, by = "ID")

    return(preds_probs)
}

predict_training_future = function(model, training, variants, n_workers = 45L) {

    meta = variants %>%
        dplyr::select(ID, NWD_ID) %>%
        dtplyr::lazy_dt(key_by = "ID")

    log_info("now generating predictions across all training data")

    slurm_template = "/net/topmed2/working/jweinstk/count_singletons/annotation_singletons/classifier/batchtools.slurm.tmpl"

    # future::plan(
    #         future.batchtools::batchtools_slurm, 
    #         template = slurm_template,
    #         resources = list(
    #             partition = "topmed", 
    #             array_throttle = n_workers, 
    #             ncpus = 2, memory = 5000, 
    #             walltime = 1200 * 60, 
    #             omp.threads = 1, 
    #             blas.threads = 1
    #         )
    # )

    future::plan(
        future::multiprocess,
        workers = n_workers
    )

    new_data = training %>% 
        dplyr::mutate(SOMATIC = as.integer(SOMATIC)) %>% 
        dplyr::filter(FILTER == "somatic") %>%
        dplyr::mutate(
            row_idx = 1:n(),
            predict_group = row_idx %/% 1e5
        )

    new_data_grouped = new_data %>%
        dplyr::group_by(predict_group) %>%
        tidyr::nest()

    log_info(glue::glue("new data grouped has {nrow(new_data_grouped)} rows"))

    preds_probs = furrr::future_imap_dfr(
            new_data_grouped$data,
            ~{
                log_info(glue::glue("now working on chunk {.y}"))
                workflows:::predict.workflow(model$final, .x, type = "prob") %>%
                dplyr::mutate(
                    row_idx = .x$row_idx
                )
            },
            .progress = TRUE
        ) %>%
        dplyr::inner_join(
            training %>% 
                dplyr::filter(FILTER == "somatic") %>%
                dplyr::select(ID, CHROM, POS, REF, ALT, VAF, AF, state, CADD_PHRED) %>%
                dplyr::mutate(row_idx = 1:n()),
            by = "row_idx"
        ) %>%
        dtplyr::lazy_dt(key_by = "ID")

    future::plan(future::sequential)
    
    # preds_probs = workflows:::predict.workflow(
    #     model$final,
    #     new_data = training %>% 
    #         dplyr::mutate(SOMATIC = as.integer(SOMATIC)) %>% 
    #         dplyr::filter(FILTER == "somatic"),
    #     type = "prob"
    # ) %>%
    #     dplyr::bind_cols(training %>% dplyr::select(ID, CHROM, POS, REF, ALT, VAF, AF, state, CADD_PHRED)) %>%
    #     dtplyr::lazy_dt(key_by = "ID")

    log_info("done generating predictions across all training data")

    preds_probs = preds_probs %>%
        dplyr::inner_join(meta, by = "ID")

    return(preds_probs %>% tibble::as_tibble(.))
}


create_data_for_torch = function(trained_model, training, variants, samples) {

    variant_meta = variants %>%
        dplyr::select(ID, NWD_ID) %>%
        dplyr::inner_join(
            samples %>%
                dplyr::select(NWD_ID, age),
            by = "NWD_ID"
        ) %>%
        dtplyr::lazy_dt(key_by = "ID")

    new_data = training %>% 
        dplyr::mutate(SOMATIC = as.integer(SOMATIC)) %>% 
        dplyr::filter(FILTER == "somatic") %>%
        dtplyr::lazy_dt(key_by = "ID")

    log_info("sorting sample meta and variant meta data")
    # sort variants
    new_data = dplyr::semi_join(new_data, variant_meta, by = "ID") %>%
        dplyr::arrange(ID)
    variant_meta = dplyr::semi_join(variant_meta, new_data, by = "ID") %>%
        dplyr::arrange(ID)

    stopifnot(all(new_data$ID == variant_meta$ID))

    x_mat = recipes::bake(trained_model$recipe, new_data) %>%
        dplyr::select(-FILTER) %>%
        dplyr::mutate(
            AF = log10(AF),
            AF = replace_Inf(AF, log10(1e-7)), # a recode monomorphic values.. 
            AF = (AF - 5) / .8,
            DISTANCE = log10(DISTANCE)
        ) %>%
        as.matrix

    stopifnot(is.numeric(x_mat))

    y = variant_meta %>% 
            dplyr::distinct(NWD_ID, age) %>%
            dplyr::pull(age)

    log_info("creating sample index map")

    sample_map_grouped = variant_meta %>%
        dplyr::mutate(idx = 1:dplyr::n()) %>%
        dplyr::select(NWD_ID, idx) %>%
        tibble::as_tibble(.) %>%
        dplyr::group_by(NWD_ID) %>%
        tidyr::nest() %>%
        dplyr::ungroup(.)

    sample_map = sample_map_grouped %>%
        dplyr::pull(data) %>%
        purrr::map(purrr::pluck, 1) %>%
        purrr::set_names(sample_map_grouped$NWD_ID)

    return(
        list(
            "x" = x_mat,
            "y" = y,
            "sample_map" = sample_map
        )
    )
}

torch_predictors = function() {
    bases = c("A", "G", "C", "T")

    muts = tidyr::expand_grid(ref = bases, alt = bases) %>%
        dplyr::filter(ref != alt) %>%
        dplyr::mutate(mut = glue::glue("{ref}.{alt}"))

    return(c(
        "SOMATIC",
        "VAF",
        "CADD_PHRED",
        "IMPACT_LOW",
        "state_Quies",
        "state_EnhG2",
        "state_ReprPCWk",
        "state_TssA",
        "state_TxWk",
        "state_Tx",
        "state_Het",
        "state_missing",
        glue::glue("MUT_TYPE_{muts$mut}")
        # glue::glue("seq_context_{1:36}")
    ))
}

torch_model = function(x, y, sample_map, lr = 0.1, iters = 180, num_threads = 20L, size = 2000) {

    tictoc::tic()

    torch_set_num_threads(num_threads)

    stopifnot(is.matrix(x))
    stopifnot(is.vector(y))
    stopifnot(ncol(x) == length(torch_predictors()))
    stopifnot(length(sample_map) == length(y))

    logger::log_info("here 1")

    ids = sample(1:length(y), size = size, replace = FALSE)
    all_nwds = names(sample_map)
    nwds = all_nwds[ids]
    sample_map = purrr::map(ids, ~sample_map[[.x]]) %>%
                    purrr::set_names(nwds)
        
    # rows = unlist(sample_map)
    logger::log_info("now selecting {length(unlist(sample_map))} mutations")

    # n x p
    x = torch_tensor(x, dtype = torch_float64())
    # k 
    y = torch_tensor((y[ids] - 60) / 11, dtype = torch_float64())

    logger::log_info("done converting to torch_tensors")

    # p x 1
    beta = torch_zeros(ncol(x), dtype = torch_float64(), requires_grad = TRUE)
    # scalar
    mu = torch_tensor(log(0.005), requires_grad = TRUE)
    # scalar
    log_sigma = torch_tensor(log(0.1), requires_grad = TRUE)
    beta_prior = distr_normal(0, 1.0)
    log_sigma_prior = distr_normal(log(0.8), 2.0)
    mu_prior = distr_normal(log(0.005), 1.0)
    intercept = torch_tensor(-0.1, requires_grad = TRUE)
    inner_intercept = torch_tensor(-1.4, requires_grad = TRUE)
    pred = torch_zeros(length(y))
    loss_tracker = torch_tensor(1000)

    likelihood = torch_tensor(0.0, requires_grad = TRUE)
    # optimizer <- optim_sgd(c(beta, mu, log_sigma), lr = lr, momentum = 0.0)
    optimizer <- optim_adam(c(beta, mu, log_sigma, intercept, inner_intercept), lr)
    # optimizer <- optim_lbfgs(c(beta, mu, log_sigma, intercept), lr, max_iter = 5)


    for(j in 1:iters) {

        logger::log_info(glue::glue("Iteration: {j}"))

        optimizer$zero_grad()
        # active = torch_sigmoid(torch_matmul(x, beta))

        # pred$zero_()
        likelihood = 0
        for(i in 1:length(ids)) {
            # pred[i] = exp(mu) * torch_sum(active[sample_map[[i]]]) 
            # pred[i] = torch_exp(mu) * torch_sum(
            #     torch_sigmoid(
            #         torch_matmul(x[sample_map[[i]], ], beta)
            #     )
            # ) 
            # expected = torch_exp(mu) * torch_sum(torch_sigmoid(
            #                                         torch_matmul(x[sample_map[[i]], ], beta)
            #                                     )
            #                             ) 
            # expected = torch_exp(mu) + torch_sum(torch_matmul(x[sample_map[[i]], ], beta)) 
            expected = intercept + torch_exp(mu) * torch_log2(torch_sum(
                torch_sigmoid(
                    inner_intercept + torch_matmul(x[sample_map[[i]], ], beta)
                    # torch_matmul(x[sample_map[[i]], ], beta)
                ))
            ) 
            pred[i] = expected

            if(i == 1) {

                logger::log_info(glue::glue("prediction = {as.numeric(expected)}, actual = {as.numeric(y[i])}"))
            }
            likelihood = likelihood + distr_normal(
                # torch_exp(mu) + torch_sum(torch_matmul(x[sample_map[[i]], ], beta)), 
                # torch_exp(mu) * torch_sum(torch_sigmoid(torch_matmul(x[sample_map[[i]], ], beta))), 
                expected, 
                torch_exp(log_sigma))$log_prob(y[i])
        }
        # browser()

        # print(pred[1])

        # y_dist = distr_normal(pred, torch_exp(log_sigma))

        # likelihood = torch_sum(y_dist$log_prob(y))
        prior = torch_sum(beta_prior$log_prob(beta)) + 
            log_sigma_prior$log_prob(log_sigma) + 
            mu_prior$log_prob(mu)

        loss = -1.0 * (likelihood + prior)
        # loss = -1.0 * (likelihood / length(y) + prior)
        # loss = -1.0 * (likelihood / length(y))

        loss_tracker = torch_cat(list(loss_tracker, loss))
        # logger::log_info(glue::glue("loss = {loss}, likelihood = {likelihood}, prior = {prior}"))
        if(j %% 3 == 0) {
            logger::log_info("loss, likelihood, prior")
            print(loss)
            print(likelihood)
            print(prior)

            print("beta[1]")
            print(beta[2])
            print("mu")
            print(torch_exp(mu))
        }

        # loss$backward(retain_graph = TRUE)
        loss$backward()
        # logger::log_info("backward()")
        optimizer$step()
        # logger::log_info("step()")
    }
    # n x p
    torch_set_num_threads(1L)
    tictoc::toc()

    return(
        list(
            "beta" = as.numeric(beta) %>% setNames(torch_predictors()),
            "sigma" = as.numeric(torch_exp(log_sigma)),
            "intercept" = as.numeric(intercept),
            "inner_intercept" = as.numeric(inner_intercept),
            "x" = as.matrix(x),
            "mu" = as.numeric(mu),
            "fitted" = as.numeric(pred),
            "y" = as.numeric(y),
            "sample_map" = sample_map,
            "nwds" = nwds,
            "indicator" = as.numeric(torch_sigmoid(inner_intercept + torch_matmul(x, beta))),
            "loss" = as.numeric(loss_tracker)
        )
    )
}

plot_fitted_model_beta = function(fitted_model) {

    fname = file.path(figure_dir(), "semisupervised_model_betas.png")

    coefs = tibble::enframe(fitted_model$beta)
    
    plot = ggplot(data = coefs, aes(y = forcats::fct_reorder(name, value, .desc = TRUE), x = value)) +
        geom_col() +
        cowplot::theme_cowplot(font_size = 12) +
        labs(x = "beta") + 
        theme(
            axis.title.y = element_blank()
        )

    ggsave(fname, plot, width = 6, height = 4, units = "in")
}
