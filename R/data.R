# PYTHON = "/net/fantasia/home/jweinstk/anaconda3/bin/python3.7"
# reticulate::use_python(PYTHON)
# reticulate::use_condaenv("cyvcf2")
# reticulate::source_python("python/parse.py")

get_path = function(label) {
    paths = c(
        "bcf" = "/net/topmed2/working/gt-release/exchange-area/freeze.10b/sites/freeze.10b.chr21.pass_and_fail.sites.bcf",
        "somatic" =  "/net/topmed2/working/jweinstk/count_singletons/singletons_spark_2021_05_25.parquet",
        "training_vcf" = file.path("output", "vcf", "training_vcf_2022_08_29.vcf"),
        "ancestry" = "/net/topmed2/working/gt-release/exchange-area/freeze.8.autosome.local.ancestry/freeze.8.global.ancestry.txt"
    )

    path = paths[label]
    return(path)
}

read_in_control_variants = function() {
    bcf = get_path("bcf")

    log_info(glue("now parsing {bcf}"))

    control_variants = walk_bcf(bcf) %>%
        as_tibble %>%
        mutate(VAF = 1 - ABE)

    log_info(glue("done parsing {bcf}"))

    return(control_variants)
}

read_in_somatic_variants = function(path = get_path("somatic")) {

    cols = c(
        "NWD_ID",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "VAF"
    )

    log_info(glue("now reading in somatic variants in {path}"))

    somatic = arrow::read_parquet(path, col_select = all_of(cols)) %>%
        dplyr::mutate(AF = 1 / (73000 * 2)) %>%
        dplyr::mutate(FILTER = "somatic") %>%
        dplyr::mutate(
            CHROM = glue("chr{CHROM}"),
            CHROM = ifelse(CHROM == "chr23", "chrX", CHROM),
            ID = glue("{CHROM}-{POS}-{REF}-{ALT}")
        )

    log_info(glue("done reading in {path}"))

    return(tibble::as_tibble(somatic))
}

prune_somatic_variants_by_sample_list_and_vaf = function(somatic_variants) {
    somatic_variants %>%
        exclude_high_vaf %>% # 12166206 -> 7498042
        fitler_to_samples # 7498042 -> 4404800
}

exclude_high_vaf = function(somatic_variants, threshold = .26) {

    log_info(glue("before filtering variants for VAF, we have {nrow(somatic_variants)} mutations"))

    somatic_variants = somatic_variants %>%
        dplyr::filter(VAF < threshold)

    log_info(glue("after filtering variants for VAF, we have {nrow(somatic_variants)} mutations"))

    return(somatic_variants)
}

fitler_to_samples = function(somatic_variants) {

    log_info(glue("before filtering samples, we have {nrow(somatic_variants)} mutations"))

    NWD_IDs = vroom::vroom("/net/topmed2/working/jweinstk/count_singletons/annotation_singletons/encore_input/chromhmm_mutation_counts_filtered_studies_ad_gte_2_lte_6_2022_04_04.tsv") %>%
        dplyr::pull(NWD_ID)

    somatic_variants = somatic_variants %>%
        dplyr::filter(NWD_ID %in% NWD_IDs)
        
    log_info(glue("after filtering samples, we have {nrow(somatic_variants)} mutations"))

    return(somatic_variants)
}

resource_dir = function() {
    "/net/topmed2/working/jweinstk/count_singletons/new_drivers/resources"
}

read_in_hg19_diff = function() {

    bed = vroom::vroom(
            file.path(resource_dir(), "hg19diff.assembly.bed.gz"),
            delim = "\t",
            col_names = c("bin", "CHROM", "START", "END", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"),
            skip = 1
            ) %>%
        dplyr::filter(score < 1000) %>% # score 1000 contigs are not so bad
        dplyr::select(
                seqnames = CHROM,
                start = START,
                end = END,
                score
                ) %>%
        plyranges::as_granges(.)

        return(bed)
}

read_in_umap_mappability = function() {
    umap = data.table::fread(
            file.path(resource_dir(), "k100.umap.bed.gz"),
            sep = "\t",
            col.names = c("CHROM", "START", "END", "type", "mappability", "strand"),
            skip = 1
            ) %>%
        tibble::as_tibble(.) %>%
        dplyr::select(
                seqnames = CHROM,
                start = START,
                end = END
                ) %>%
        plyranges::as_granges(.)

        return(umap)
}

read_common_dups = function(dir = "/net/topmed2/working/jweinstk/count_singletons/new_drivers/mappability_analysis/output/high_maf_dups") {

    files = list.files(dir, full.names = TRUE)

    dups = purrr::map_dfr(
            files,
            ~vroom::vroom(.x, col_names = c("seqnames", "start", "end", "ID"))
            ) %>%
    dplyr::mutate(strand = "+") %>%
    plyranges::as_granges(.)

    dups
}

exclude_beds = function(variants, list_of_beds) {
    # hg19diff 22028631 -> 15093246
    # dups 15093246 -> 12301208

    stopifnot(is.list(list_of_beds))

    logger::log_info(glue::glue("Identified {length(list_of_beds)} beds to exclude"))

    final = variants %>%
        rename_to_genomic_ranges %>%
        plyranges::as_granges(.)

    for(i in 1:length(list_of_beds)) {
        label = names(list_of_beds)[i]
            logger::log_info(glue::glue("Prior to filtering by {label}, there are {length(final)} variants"))
            final = plyranges::filter_by_non_overlaps(final, list_of_beds[[i]])
            logger::log_info(glue::glue("After filtering by {label}, there are {length(final)} variants"))
    }

    return(final %>% GenomicRanges::as.data.frame(.) %>% rename_from_genomic_ranges)
}

overlap_beds = function(variants, list_of_beds) {
    # 12301208 -> 12166206

    stopifnot(is.list(list_of_beds))

    logger::log_info(glue::glue("Identified {length(list_of_beds)} beds to exclude"))

    final = variants %>%
        rename_to_genomic_ranges %>%
        plyranges::as_granges(.)

    for(i in 1:length(list_of_beds)) {
        label = names(list_of_beds)[i]
            logger::log_info(glue::glue("Prior to filtering by {label}, there are {length(final)} variants"))
            final = plyranges::filter_by_overlaps(final, list_of_beds[[i]])
            logger::log_info(glue::glue("After filtering by {label}, there are {length(final)} variants"))
    }

    return(final %>% GenomicRanges::as.data.frame(.) %>% rename_from_genomic_ranges)
}

sample_variants = function(variants, sample_size = 20000) {

    set.seed(1)

    variants %>%
        dplyr::sample_n(size = sample_size, replace = FALSE)
}


write_to_vcf = function(df, path = get_path("training_vcf")) {

    res = df %>%
        dplyr::select(
            `#CHROM` = CHROM,
            POS,
            ID,
            REF,
            ALT
        ) %>%
        readr::write_tsv(path)

    return(path)
}

execute_vep = function(in_path) {

    log_info("now running VEP")
    cmd = "bash annotate_vep.sh"
    system(cmd)
    log_info("done running VEP")

    date = "2022_08_29"
    out_path = file.path("output", "txt", glue("training_vcf_annotated_{date}.txt"))

    return(out_path)
}

filter_pick = function(x) {
    x %>%
        dplyr::filter(!is.na(PICK))
}

make_variant_meta = function(df) {

    df %>%
        dplyr::select(ID, CHROM, POS, REF, ALT, FILTER, VAF, AF) %>%
        dtplyr::lazy_dt(key_by = ID)

}

read_in_training_annotated = function(in_path, variants) {

    log_info(glue("now reading in {in_path}"))

    df = vroom::vroom(
        in_path, 
        comment = "##", 
        na = c("", "NA", "-"),
        col_select = c(
            `#Uploaded_variation`, 
            Location,
            Gene,
            Feature,
            Feature_type,
            Consequence,
            Existing_variation,
            IMPACT,
            DISTANCE,
            FLAGS,
            PICK,
            SYMBOL,
            CLIN_SIG,
            SOMATIC,
            PHENO,
            MOTIF_NAME,
            TRANSCRIPTION_FACTORS,
            HIGH_INF_POS,
            CADD_PHRED
        )
    ) %>%
        dplyr::rename(ID = `#Uploaded_variation`) %>%
        dtplyr::lazy_dt(key_by = ID)

    variants = variants %>%
        make_variant_meta 

    print(df)
    log_info(glue("done reading in {in_path}; read {nrow(as.data.frame(df))} rows"))

    df = df %>%
        filter_pick

    log_info(glue("done filtering based on FILTER_PICK; now {nrow(as.data.frame(df))} rows"))

    df = df %>%
        dplyr::right_join(variants, by = "ID")

    log_info(glue("done joining to labels; now {nrow(as.data.frame(df))} rows"))

    return(tibble::as_tibble(df))
}

rename_to_genomic_ranges = function(df) {
    if(any(c("END", "end") %in% names(df))) {
        result = df %>%
            dplyr::rename(
                seqnames = CHROM,
                start = POS
            ) 
    } else {
        result = df %>%
            dplyr::rename(
                seqnames = CHROM,
                start = POS
            )  %>%
            dplyr::mutate(
                end = start
            )
    }

    return(result)
}

rename_from_genomic_ranges = function(df) {
    df %>%
        dplyr::rename(
            POS = start,
            CHROM = seqnames
        )
}

read_in_cd34 = function() {
    chromHMM = vroom::vroom(
                        "/net/topmed2/working/jweinstk/count_singletons/new_drivers/resources/BSS00233_18_CALLS_segments.bed.gz",
                        col_names = c("CHROM", "START", "END", "state", "v5", "v6", "v7", "v8", "v9")
                ) %>%
                dplyr::select(CHROM, START, END, state) %>%
                dplyr::mutate(dplyr::across(c(START, END), as.integer)) %>%
                dplyr::rename(
                    seqnames = CHROM,
                    start = START,
                    end = END
                ) %>%
                plyranges::as_granges(.)

    return(chromHMM)
}

join_chromatin = function(training, cd34) {

    log_info("now merging in CD34 chromatin annotations")
    dfm = training %>%
        rename_to_genomic_ranges %>%
        plyranges::as_granges(.) %>%
        plyranges::join_overlap_left(cd34) %>%
        GenomicRanges::as.data.frame(.) %>%
        tibble::as_tibble(.) %>%
        rename_from_genomic_ranges

    log_info("done merging in CD34 chromatin annotations")
    return(dfm)
}

read_in_mCAs = function() {

    df = vroom::vroom("/net/topmed2/working/jweinstk/count_singletons/mCA_calls/TOPMed_mCA_calls.txt") %>%
        dplyr::rename(NWD_ID = sample_id)

    return(df)
}

read_in_LOY = function() {

    df = vroom::vroom("/net/topmed2/working/jweinstk/count_singletons/mCA_calls/LoY_phenotypes_v4.txt") %>%
        dplyr::rename(NWD_ID = NWDID)

    return(df)
}

read_in_mCA_manifest = function() {

    mCA_calls = read_in_mCAs()

    mCA_samples = readLines("/net/topmed2/working/jweinstk/count_singletons/mCA_calls/samples.inFinal.Analyses.txt")

    LOY_calls = read_in_LOY()

    LOY_samples = LOY_calls %>% dplyr::pull(NWD_ID)

    stopifnot(length(setdiff(LOY_samples, mCA_samples)) == 0)

    loss_of_X = mCA_calls %>%
        dplyr::filter(chrom == "chrX")

    exclude_loss_of_X = mCA_calls %>%
        dplyr::filter(chrom != "chrX")

    loss_of_Y = LOY_calls %>%
        dplyr::filter(loy == 1)

    df = tibble::tibble(
                NWD_ID = mCA_samples
            ) %>%
            dplyr::mutate(
                autosomal_mCA = NWD_ID %in% exclude_loss_of_X$NWD_ID,
                LOX = NWD_ID %in% loss_of_X$NWD_ID,
                any_mCA_LOX = LOX | autosomal_mCA,
                LOY = NWD_ID %in% loss_of_Y$NWD_ID,
                any_mCA_LOX_LOY = any_mCA_LOX | LOY
            )

    return(df)
}

attach_sequence_context = function(training) {

    num_workers = 8L
    future::plan(future::multicore, workers = num_workers)
    N_FLANK = 4L

    log_info("now extracting flanking sequence")

    training = training %>%
        dplyr::group_by(CHROM) %>%
        tidyr::nest() %>%
        dplyr::mutate(
            sequence = furrr::future_map2(CHROM, data, ~{
                BSgenome::getSeq(
                    BSgenome.Hsapiens.UCSC.hg38, 
                    names = .x,
                    start = .y$POS - N_FLANK,
                    end   = .y$POS + N_FLANK,
                    as.character = TRUE
                )
            })
        ) %>%
        tidyr::unnest(cols = c(data, sequence)) %>%
        dplyr::ungroup(.)

    future::plan(future::sequential)

    log_info("done extracting flanking sequence; now converting to one hot")

    Rcpp::sourceCpp("src/encode.cpp")
    encodings = one_hot(training$sequence) %>%
        tibble::as_tibble(.) %>%
        setNames(glue::glue("seq_context_{1:((N_FLANK * 2 + 1) * 4L)}"))

    training = dplyr::bind_cols(training, encodings)

    log_info("done converting to one hot")

    return(training)
}

join_labels = function(training, labels) {
    training %>%
       dplyr::inner_join(labels %>% dplyr::select(ID, FILTER, VAF, AF), by = "ID") 
}


encore_input = function() {
    "../encore_input"
}

read_in_sample_meta = function() {
    vroom::vroom(file.path(encore_input(), "chromhmm_mutation_counts_filtered_studies_ad_gte_2_lte_6_2022_04_04.tsv"))
}

create_scatter = function(df, x, y) {
    ggplot(data = df, aes({{x}}, {{y}})) +
        geom_point(alpha = .3) + 
        geom_smooth() +
        scale_x_continuous(trans = "log2") +
        scale_y_continuous(trans = "log2") +
        labs(x = "Number of estimated mutations", y = "Numberof of observed mutations") +
        cowplot::theme_cowplot(font_size = 12)
}

create_scatter_color = function(df, x, y) {
    ggplot(data = df, aes({{x}}, {{y}}, color = age)) +
        geom_point(alpha = .3) + 
        # geom_smooth() +
        scale_colour_viridis_c(alpha = .3) +
        scale_x_continuous(trans = "log2") +
        scale_y_continuous(trans = "log2") +
        labs(x = "Number of estimated mutations", y = "Numberof of observed mutations") +
        cowplot::theme_cowplot(font_size = 12)
}


write_encore_output = function(sample_meta, predictions) {

    fnames = c(
        "tsv" = file.path(encore_input(), "classifier_output_2022_09_11.tsv"),
        "plot" = file.path(figure_dir(), "scatter_classifier.png")
    )

    estimated_counts = predictions %>%
        dplyr::group_by(NWD_ID) %>%
        dplyr::summarize(
            observed_counts = dplyr::n(),
            classifier_counts = sum(.pred_somatic),
            residual = observed_counts - classifier_counts,
            residual_per_variant = residual / dplyr::n()
        ) %>%
        dplyr::ungroup(.)

    sample_meta = sample_meta %>%
        dplyr::inner_join(estimated_counts, by = "NWD_ID")

    sample_meta %>%
        readr::write_tsv(fnames["tsv"])

    ggsave(
        fnames["plot"],
        create_scatter_color(sample_meta, classifier_counts, observed_counts),
        width = 6, height = 4, units = "in"
    )

    return(fnames)

}
write_semi_supervised_encore_output = function(sample_meta, fitted_model, torch_input, mCA, global_ancestry) {

    date = "2025_3_13"
    prefix = "semi_supervised"

    fnames = c(
        "tsv_all" = file.path(encore_input(), glue::glue("{prefix}_output_{date}.tsv")),
        "tsv_female" = file.path(encore_input(), glue::glue("{prefix}_female_output_{date}.tsv")),
        "tsv_male" = file.path(encore_input(), glue::glue("{prefix}_male_output_{date}.tsv")),
        "tsv_exclude_LOX" = file.path(encore_input(), glue::glue("{prefix}_exclude_LOX_output_{date}.tsv")),
        "tsv_exclude_LOY" = file.path(encore_input(), glue::glue("{prefix}_exclude_LOY_output_{date}.tsv")),
        "tsv_contains_LOX" = file.path(encore_input(), glue::glue("{prefix}_contains_LOX_output_{date}.tsv")),
        "tsv_contains_LOY" = file.path(encore_input(), glue::glue("{prefix}_contains_LOY_output_{date}.tsv")),
        "tsv_exclude_autosomal_mCA" = file.path(encore_input(), glue::glue("{prefix}_exclude_autosomal_mCA_output_{date}.tsv")),
        "tsv_contains_autosomal_mCA" = file.path(encore_input(), glue::glue("{prefix}_contains_autosomal_mCA_output_{date}.tsv")),
        # "tsv_exclude_any_mCA_LOX" = file.path(encore_input(),glue::glue("{prefix}_exclude_any_mCA_LOX_output_{date}.tsv")),
        # "tsv_contains_any_mCA" = file.path(encore_input(),glue::glue("{prefix}_contains_any_mCA_output_{date}.tsv")),
        "tsv_contains_any_mCA_LOX_LOY" = file.path(encore_input(),glue::glue("{prefix}_contains_any_mCA_LOX_LOY_output_{date}.tsv")),
        "tsv_exclude_any_mCA_LOX_LOY" = file.path(encore_input(),glue::glue("{prefix}_exclude_any_mCA_LOX_LOY_output_{date}.tsv")),
        "tsv_EUR" = file.path(encore_input(), glue::glue("{prefix}_only_EUR_{date}.tsv")),
        "tsv_AFR" = file.path(encore_input(), glue::glue("{prefix}_only_AFR_{date}.tsv")),
        "tsv_EAS" = file.path(encore_input(), glue::glue("{prefix}_only_EAS_{date}.tsv")),
        "tsv_AMR" = file.path(encore_input(), glue::glue("{prefix}_only_AMR_{date}.tsv")),
        "plot" = file.path(figure_dir(), glue::glue("scatter_{prefix}_{date}.png"))
    )

    estimated_counts = purrr::imap_dfr(
       torch_input$sample_map,
       ~tibble::tibble(
           estimated_counts = sum(fitted_model$indicator[.x]),
           observed_counts = length(.x),
           NWD_ID = .y
       )
    ) %>%
        dplyr::mutate(
            residual = observed_counts - estimated_counts
        )

    sample_meta = sample_meta %>%
        dplyr::inner_join(estimated_counts, by = "NWD_ID") %>%
        dplyr::left_join(mCA, by = "NWD_ID") 

    sample_meta %>%
        readr::write_tsv(fnames["tsv_all"])

    sample_meta %>%
        dplyr::filter(INFERRED_SEX == 2) %>%
        readr::write_tsv(fnames["tsv_female"])

    sample_meta %>%
        dplyr::filter(INFERRED_SEX == 1) %>%
        readr::write_tsv(fnames["tsv_male"])

    sample_meta %>%
        dplyr::filter(!is.na(any_mCA_LOX_LOY)) %>%
        dplyr::filter(LOX) %>%
        readr::write_tsv(fnames["tsv_contains_LOX"])

    sample_meta %>%
        dplyr::filter(!is.na(any_mCA_LOX_LOY)) %>%
        dplyr::filter(LOY) %>%
        readr::write_tsv(fnames["tsv_contains_LOY"])

    sample_meta %>%
        dplyr::filter(!is.na(any_mCA_LOX_LOY)) %>%
        dplyr::filter(!LOX) %>%
        readr::write_tsv(fnames["tsv_exclude_LOX"])

    sample_meta %>%
        dplyr::filter(!is.na(any_mCA_LOX_LOY)) %>%
        dplyr::filter(!LOY) %>%
        readr::write_tsv(fnames["tsv_exclude_LOY"])

    sample_meta %>%
        dplyr::filter(!is.na(any_mCA_LOX_LOY)) %>%
        dplyr::filter(autosomal_mCA) %>%
        readr::write_tsv(fnames["tsv_contains_autosomal_mCA"])

    sample_meta %>%
        dplyr::filter(!is.na(any_mCA_LOX_LOY)) %>%
        dplyr::filter(!autosomal_mCA) %>%
        readr::write_tsv(fnames["tsv_exclude_autosomal_mCA"])

    sample_meta %>%
        dplyr::filter(!is.na(any_mCA_LOX_LOY)) %>%
        dplyr::filter(any_mCA_LOX_LOY) %>%
        readr::write_tsv(fnames["tsv_contains_any_mCA_LOX_LOY"])

    sample_meta %>%
        dplyr::filter(!is.na(any_mCA_LOX_LOY)) %>%
        dplyr::filter(!any_mCA_LOX_LOY) %>%
        readr::write_tsv(fnames["tsv_exclude_any_mCA_LOX_LOY"])

    logger::log_info("Now here 1")

    sample_meta %>%
        dplyr::inner_join(global_ancestry, by = "NWD_ID") %>%
        dplyr::filter(ancestry == "europe") %>%
        readr::write_tsv(fnames["tsv_EUR"])

    sample_meta %>%
        dplyr::inner_join(global_ancestry, by = "NWD_ID") %>%
        dplyr::filter(ancestry == "sub-saharan africa") %>%
        readr::write_tsv(fnames["tsv_AFR"])

    sample_meta %>%
        dplyr::inner_join(global_ancestry, by = "NWD_ID") %>%
        dplyr::filter(ancestry == "east asia") %>%
        readr::write_tsv(fnames["tsv_EAS"])

    sample_meta %>%
        dplyr::inner_join(global_ancestry, by = "NWD_ID") %>%
        dplyr::filter(ancestry == "native america") %>%
        readr::write_tsv(fnames["tsv_AMR"])

    ggsave(
        fnames["plot"],
        create_scatter_color(
                sample_meta, 
                estimated_counts, 
                observed_counts
        ),
        width = 6, height = 4, units = "in"
    )

    return(fnames)
}

read_in_ancestry = function() {
    
    path = get_path("ancestry")

    global_ancestry = data.table::fread(path) %>% 
        as.data.frame %>%
        setNames(
            c('NWD_ID',
              'sub-saharan africa',
              'central and south asia', 
              'east asia', 
              'europe', 
              'native america', 
              'oceania', 
              'middle east', 
              'pi', 
              'stage'
            )
        ) %>%
    as_tibble %>%
    dplyr::select(-pi, -stage)

    global_ancestry_labels = global_ancestry %>% 
        tidyr::gather(ancestry, proportion, -NWD_ID) %>%
        dplyr::group_by(NWD_ID) %>%
        dplyr::filter(proportion == max(proportion)) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(
            NWD_ID,
            ancestry
        )

    return(global_ancestry_labels)
}
