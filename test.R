func <- function(data, select_col1, select_col2) {
    select_col1 <- rlang::sym(select_col1)
    select_col2 <- rlang::sym(select_col2)
    dplyr::select(
        .data = data,
        sample_id = {{ select_col1 }},
        project = {{ select_col2 }}
    )
}
