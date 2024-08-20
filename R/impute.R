replace_NA_lgl = function(x) {
    dplyr::coalesce(x, FALSE)
}

replace_NA_0 = function(x) {
    dplyr::coalesce(x, 0)
}

replace_NA_missing = function(x) {
    dplyr::coalesce(x, "missing")
}

replace_NA_non_zero = function(x, val) {
    dplyr::coalesce(x, val)
}

replace_Inf = function(x, val) {
    ifelse(is.infinite(x), val, x)
}

impute_training = function(df) {

    mean_distance_somatic = 2696

    df %>%
        dplyr::mutate(
            dplyr::across(where(is.logical), replace_NA_lgl),
            dplyr::across(where(is.character), replace_NA_missing),
            dplyr::across(DISTANCE, ~replace_NA_non_zero(replace_Inf(.x, 1e6), mean_distance_somatic)),
            dplyr::across(CADD_PHRED, ~replace_NA_0(.x)) 
        )
}

replace_labels = function(df) {
    df %>%
        dplyr::mutate(
            FILTER = dplyr::case_when(
                FILTER == "somatic" ~ "somatic",
                FILTER == "PASS" ~ "germline",
                TRUE ~ "artifact" 
            )
        )
}


transform_values = function(df) {
    df %>%
        dplyr::mutate(
            VAF = (VAF - .37) / .23, # zscore
            MAF = pmin(AF, 1 - AF), # now MAF
            MAF = log10(MAF),
            MAF = replace_Inf(MAF, log10(1e-7)), # a recode monomorphic values.. 
            MAF = (MAF - 5) / .8
        ) %>%
        dplyr::mutate(
            SOMATIC = as.integer(SOMATIC),
            MUT_TYPE = glue::glue("{REF}-{ALT}")
        )
}

exclude_indels = function(df) {
    df %>%
        dplyr::filter(nchar(ALT) == 1 & nchar(REF) == 1)
}

impute_and_transform = function(df) {
    df %>%
        exclude_indels %>%
        impute_training %>%
        transform_values %>%
        replace_labels
}
