#' @title Create a 4way plot
#' @order 2
#' @usage NULL
#' @description Create a 4way plot to compare the logFC values from
#'  two contrasts of differential gene expression.
#' @param DGEdata The object to plot from:
#' \itemize{
#' \item \code{limma:} A \link[limma:MArrayLM-class]{MArrayLM} object from
#'  \link[limma]{eBayes} or \link[limma]{treat}
#' \item \code{edgeR:} A list of \link[edgeR:DGELRT-class]{DGELRT} objects
#'  from \link[edgeR]{glmQLFTest}, \link[edgeR]{glmTreat},
#'  or \link[edgeR]{glmLRT}
#' \item \code{DESeq2:} a \link[DESeq2:DESeqDataSet-class]{DESeqDataSet} from
#'  \link[DESeq2]{DESeq} or a list of
#'  \link[DESeq2:DESeqResults-class]{DESeqResults} from \link[DESeq2]{results}
#'  \item \code{Other packages:} A list of data.frames,
#'   see details section for more infromation
#' }
#' @param x Character specifying name of DGE results within object
#'  for the x-axis
#' @param y Character specifying name of DGE results within object
#'  for the y-axis
#' @param sep Character specifying the separator between conditions
#'  for the contrast
#' @param ID Column name for gene ID
#' @param symbol Column name for gene symbol description
#' @param logFC Column name for logFC values
#' @param FDR Column name for FDR values
#' @param FDRcutoff Numeric for the FDR cut-off for DEGs, default is 0.05
#' @param logFCcutoff Numeric for the absolute Log2FC cut-off for DEGs,
#'  default is 1
#' @param label Character vector specifying genes to label
#'  (FALSE for none, TRUE for all blue)
#' @param colorVector Character vector of colors in the following order:
#' "not significant", "significant in x", "significant in y",
#'  "significant in both"
#' @param lineColor Color of lines
#' @param textSize Numeric specifying size of text with gene
#'  overlap category totals
#' @param textNudge Numeric specifying nudge of text with gene
#'  overlap category totals
#' @param ... Support for additional arguments used internally by
#'  \code{gg4way.MArrayLM}, \code{gg4way.list},
#'  and \code{gg4way.DESeqDataSet}
#' @details
#' When a list of data.frames is provided to the \code{DGEdata} argument,
#' they should have the following column names and data:
#'  \tabular{ll}{
#'   \code{ID} \tab Character vector with the feature ID (i.e. EnsemblID) \cr
#'   \tab \cr
#'   \code{symbol} \tab Optional character vector with gene symbol for labels \cr
#'   \tab \cr
#'   \code{LogFC} \tab Numeric with the logFC \cr
#'   \tab \cr
#'   \code{FDR} \tab Numeric with the FDR \cr
#'   }
#'
#' The correlation coefficient is useful for comparing across multiple plots.
#'  It's important to consider whether there are any common factors when
#'  comparing values, since that can result in a larger value.
#' @return A \link[ggplot2]{ggplot}
#' @export
#'
gg4way <- function(DGEdata,
                   ...){
    UseMethod("gg4way")
}

#' @rdname gg4way
#' @order 1
#' @name gg4way
#' @importFrom rlang warn sym
#' @importFrom glue glue_collapse
#' @importFrom dplyr filter pull
#' @importFrom janitor tabyl
#' @importFrom purrr set_names
#' @examples
#' data("airwayFit")
#' airwayFit |>
#'     gg4way(x = "N61311 vs N052611",
#'            y = "N061011 vs N052611")
#' @export
#'
gg4way.default <- function(DGEdata,
                           x = NULL,
                           y = NULL,
                           ID = "ID",
                           symbol = "symbol",
                           logFC = "logFC",
                           FDR = "adj.P.Val",
                           sep = " vs ",
                           FDRcutoff = 0.05,
                           logFCcutoff = 1,
                           label = FALSE,
                           colorVector = c("grey80", "firebrick",
                                           "forestgreen", "mediumblue"),
                           lineColor = "grey60",
                           textSize = 4,
                           textNudge = 0.25,
                           ...) {
    stopifnot(!is.null(x) | !is.null((y)))
    if (is.null(symbol)) {
        symbol <- "ID"
    }
    stopifnot(c(x,y) %in% names(DGEdata))

    missingFeatures <- c(DGEdata[[x]]$ID, DGEdata[[y]]$ID) %>%
        janitor::tabyl() %>%
        dplyr::filter(n == 1) %>%
        dplyr::pull(1)

    if (!identical(missingFeatures, character(0))) {
        rlang::warn(paste(x, "and", y, "don't have some genes IDs in common.",
                          length(missingFeatures), "IDs will be filtered out:",
                          glue::glue_collapse({missingFeatures},
                                              sep = ", ",
                                              width = 37)),
                    use_cli_format = TRUE)
    }

    DGEtibble <- .prepareData(DGEdata = DGEdata,
                              x = x,
                              y = y,
                              ID = ID,
                              symbol = symbol,
                              logFC = logFC,
                              FDR = FDR,
                              logFCcutoff = logFCcutoff,
                              FDRcutoff = FDRcutoff)

    corRes <- .testCor(DGEtibble = DGEtibble)

    totalTibble <- .totalCounts(DGEtibble = DGEtibble,
                                x = x,
                                y = y,
                                logFCcutoff = logFCcutoff)

    colorKey <- colorVector |>
        purrr::set_names(c("Not Significant",
                           paste("Significant in", x),
                           paste("Significant in", y),
                           "Significant in Both"))

    textKey <- .prepareAnnotations(totalTibble = totalTibble,
                                   colorKey = colorKey,
                                   textNudge = textNudge)

    p1 <- .plot4way(DGEtibble = DGEtibble,
                    x = x,
                    y = y,
                    sep = sep,
                    logFCcutoff = logFCcutoff,
                    lineColor = lineColor,
                    colorKey = colorKey,
                    corRes = corRes,
                    textKey = textKey,
                    hjust = hjust,
                    vjust = vjust,
                    textSize = textSize,
                    label = label)

    class(p1) <- c("gg4way", class(p1))

    p1
}

#' @rdname gg4way
#' @usage NULL
#' @importFrom limma topTable
#' @importFrom purrr set_names map
#' @importFrom stringr str_replace
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom magrittr %>%
#' @export
#'
gg4way.MArrayLM <- function(DGEdata,
                            ID = ID,
                            symbol = symbol,
                            logFC = logFC,
                            FDR = FDR,
                            ...){
    if (missing(ID)) {ID <- "ID"}
    if (missing(symbol)) {symbol <- "symbol"}
    if (missing(logFC)) {logFC <- "logFC"}
    if (missing(FDR)){ FDR <- "adj.P.Val"}

    ## magrittr pipe used due to need for unnamed placeholder
    DGEdata$contrasts %>%
        colnames() %>%
        purrr::set_names() %>%
        purrr::map(~ DGEdata %>%
                       {if("treat.lfc" %in% names(DGEdata)){
                           limma::topTreat(fit = .,
                                           n = Inf,
                                           coef = .x)
                       }else{
                           limma::topTable(fit = .,
                                           n = Inf,
                                           coef = .x)
                       }} %>%
                       tibble::rownames_to_column() %>%
                       tibble::as_tibble()) %>%
        purrr::set_names(names(.) %>%
                             stringr::str_replace("[-]", "vs")) %>%
        gg4way.default(DGEdata = .,
                       ID = ID,
                       symbol = symbol,
                       logFC = logFC,
                       FDR = FDR,
                       ...)
}

#' @rdname gg4way
#' @usage NULL
#' @importFrom edgeR topTags
#' @importFrom purrr set_names map
#' @importFrom stringr str_replace
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @export
#'
gg4way.list <- function(DGEdata,
                        ID = ID,
                        symbol = symbol,
                        logFC = logFC,
                        FDR = FDR,
                        ...){
    if (all(vapply(DGEdata, is, logical(1), "DGELRT"))) {
        if (missing(ID)) {ID <- "ID"}
        if (missing(symbol)) {symbol <- "symbol"}
        if (missing(logFC)) {logFC <- "logFC"}
        if (missing(FDR)) {FDR <- "FDR"}

    ## magrittr pipe used due to need for unnamed placeholder
        DGEdata %>%
            names() %>%
            purrr::set_names() %>%
            purrr::map(~ DGEdata[[.x]] %>%
                           edgeR::topTags(n = Inf) %>%
                           .$table %>%
                           tibble::as_tibble()) %>%
            purrr::set_names(names(.) %>%
                                 stringr::str_replace("[-]", "vs")) %>%
            gg4way.default(DGEdata = .,
                           ID = ID,
                           symbol = symbol,
                           logFC = logFC,
                           FDR = FDR,
                           ...)

    }else if (all(vapply(DGEdata, is, logical(1), "DESeqResults"))) {
        if (missing(ID)) {ID <- "ID"}
        if (missing(symbol)) {symbol <- "ID"}
        if (missing(logFC)) {logFC <- "log2FoldChange"}
        if (missing(FDR)) {FDR <- "padj"}

        DGEdata %>%
            purrr::map(~ .x %>%
                           as.data.frame() %>%
                           tibble::rownames_to_column("ID") %>%
                           tibble::as_tibble()) %>%
            gg4way.default(DGEdata = .,
                           ID = ID,
                           symbol = symbol,
                           logFC = logFC,
                           FDR = FDR,
                           ...)

    }else if (all(vapply(DGEdata, is.data.frame, logical(1)))) {
        if (missing(ID)) {ID <- "ID"}
        if (missing(symbol)) {symbol <- "symbol"}
        if (missing(logFC)) {logFC <- "logFC"}
        if (missing(FDR)) {FDR <- "adj.P.Val"}

        DGEdata %>%
            gg4way.default(DGEdata = .,
                           ID = ID,
                           symbol = symbol,
                           logFC = logFC,
                           FDR = FDR,
                           ...)
    }else{
        stop("Objects in list or their combination are not supported")
    }
}

#' @rdname gg4way
#' @usage NULL
#' @importFrom DESeq2 results
#' @importFrom purrr set_names map
#' @importFrom stringr str_subset str_replace_all str_remove str_trim
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom magrittr %>%
#' @export
#'
gg4way.DESeqDataSet <- function(DGEdata,
                                ID = ID,
                                symbol = symbol,
                                logFC = logFC,
                                FDR = FDR,
                                ...){
    if (missing(ID)) {ID <- "ID"}
    if (missing(symbol)) {symbol <- "ID"}
    if (missing(logFC)) {logFC <- "log2FoldChange"}
    if (missing(FDR)) {FDR <- "padj"}

    ## magrittr pipe used due to need for unnamed placeholder
    DGEdata %>%
        DESeq2::resultsNames() %>%
        stringr::str_subset("Intercept", negate = TRUE) %>%
        purrr::set_names() %>%
        purrr::map(~ DGEdata %>%
                       DESeq2::results(name = .x)) %>%
        purrr::set_names(names(.) %>%
                             stringr::str_replace_all("_", " ") %>%
                             stringr::str_remove(stringr::word(., 1)) %>%
                             stringr::str_trim()) %>%
        purrr::map(~ .x %>%
                       as.data.frame() %>%
                       tibble::rownames_to_column("ID") %>%
                       tibble::as_tibble()) %>%
        gg4way.default(DGEdata = .,
                       ID = ID,
                       symbol = symbol,
                       logFC = logFC,
                       FDR = FDR,
                       ...)
}

#' @rdname gg4way
#' @usage NULL
#' @importFrom purrr map
#' @importFrom tibble rownames_to_column
#' @export
#'
gg4way.SimpleList <- function(DGEdata = .,
                              ID = ID,
                              symbol = symbol,
                              logFC = logFC,
                              FDR = FDR,
                              ...){
    if (missing(ID)) {ID <- "ID"}
    if (missing(symbol)) {symbol <- "symbol"}
    if (missing(logFC)) {logFC <- "LogFC"}
    if (missing(FDR)) {FDR <- "FDR"}

    DGEdata |>
        purrr::map(~ .x |>
                       as.data.frame() |>
                       tibble::rownames_to_column("ID")) |>
        gg4way.default(DGEdata = _,
                       ID = ID,
                       symbol = symbol,
                       logFC = logFC,
                       FDR = FDR,
                       ...)
}
