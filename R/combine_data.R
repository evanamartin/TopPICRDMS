#' Combine data sets
#'
#' Combines a group of data sets output from TopPIC. This function also adds
#' several variables that will be used throughout the top down workflow and
#' removes variables that are not used.
#'
#' @param x A list of data frames. All data frames must have the same number of
#'   columns and the names of the columns must match exactly.
#'
#' @return A data table with the following variables added:
#'   \tabular{ll}{
#'   \code{ProjID} \tab The project ID for each individual. \cr
#'   \code{LC_Column} \tab The LC column the sample was run on. \cr
#'   \code{mz} \tab The m/z value calculated from the \code{Precursor mass}
#'   and \code{Charge} variables. \cr
#'   \code{Gene} \tab The gene, extracted from the \code{Protein description}
#'   variable, a protein belongs to. \cr
#'   \code{isDecoy} \tab A logical value denoting whether or not a protein
#'   is scrambled. \cr
#'   \code{RTmin} \tab The retention time in minutes. \cr
#'   }
#'
#' @export
#'
combine_data <- function (x = list()) {

  # Create a character vector of variables that are not needed at any point in
  # the top down workflow. They will be removed with dplyr::any_of. The function
  # will not throw an error if any of the variables listed are not present in
  # the data set.
  useless <- c("Data file name", "Prsm ID", "Spectrum ID", "Fragmentation",
               "Retention time", "#peaks", "Proteoform ID", "Feature score",
               "MIScore", "#variable PTMs", "#matched peaks",
               "#matched fragment ions", "Q-value (spectral FDR)",
               "Proteoform FDR", "First residue", "Last residue")

  x <- rbindlist(x) %>%
    # Remove data.table class.
    # tibble() %>%
    # Add useful variables.
    dplyr::filter(!grepl("Contaminant", `Protein accession`)) %>%
    dplyr::filter(grepl("GN=", `Protein description`)) %>%
    dplyr::mutate(ProjID = stringr::str_extract(Dataset, "\\d{8}")) %>%
    dplyr::mutate(LC_Column = stringr::str_sub(Dataset, start = -8)) %>%
    dplyr::mutate(
      mz = (`Precursor mass` + Charge * 1.007276466621) / Charge
    ) %>%
    dplyr::mutate(Gene = sub(".*GN=(\\S+).*","\\1",`Protein description`)) %>%
    dplyr::mutate(isDecoy = grepl("^XXX", `Protein accession`)) %>%
    dplyr::mutate(RTmin = `Retention time` / 60) %>%
    # Remove not so useful variables. In other words, remove variables that will
    # never be used or thought of again.
    dplyr::select(-dplyr::any_of(useless))

  return (x)

}
