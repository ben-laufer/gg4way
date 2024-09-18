#' @title Missing names check
#' @description Check for missing names in the DGEdata object
#' @keywords internal
#' @inheritParams gg4way
#' @return A character
#' @importFrom magrittr extract
#' @importFrom purrr map imap list_c
#' @importFrom dplyr setdiff
#' @importFrom glue glue glue_collapse
#' @importFrom rlang abort
#'
.checkNames <- function(DGEdata = DGEdata,
                        x = x,
                        y = y,
                        ID = ID,
                        symbol = symbol,
                        logFC = logFC,
                        FDR = FDR){
    missingList <- DGEdata |>
        magrittr::extract(c(x, y)) |>
        purrr::map(\(contrastData)
                   dplyr::setdiff(c(ID, symbol, logFC, FDR),
                                  names(contrastData)))

    errorMessage <- missingList |>
        purrr::imap(\(missingNames, listName)
                    if (length(missingNames) > 0) {
                        glue::glue("{listName} is missing column name(s): {missing}",
                                   missing = glue::glue_collapse(missingNames,
                                                                 sep = ", "))
                    }) |>
        purrr::list_c()

    if( length(errorMessage) > 0) {
        rlang::abort(errorMessage)
    }
}

#' @title Missing features warning
#' @description Warn about features not shared between x and y
#' @keywords internal
#' @inheritParams gg4way
#' @return A character
#' @importFrom dplyr pull filter
#' @importFrom janitor tabyl
#' @importFrom rlang warn
#' @importFrom glue glue_collapse
#'
.checkFeatures <- function(DGEdata = DGEdata,
                           x = x,
                           y = y,
                           ID = ID){
    missingFeatures <- c(DGEdata[[x]] |>
                             dplyr::pull(!!ID),
                         DGEdata[[y]] %>%
                             dplyr::pull(!!ID)) |>
        janitor::tabyl() |>
        dplyr::filter(n == 1) |>
        dplyr::pull(1)

    if (!identical(missingFeatures, character(0))) {
        rlang::warn(paste(x, "and", y, "don't have some genes IDs in common.",
                          length(missingFeatures), "IDs will be filtered out:",
                          glue::glue_collapse({missingFeatures},
                                              sep = ", ",
                                              width = 37)))
    }
}

#' @title Prepare data
#' @description Prepare data for a 4way plot
#' @keywords internal
#' @inheritParams gg4way
#' @return A \link[tibble:tbl_df-class]{tibble}
#' @importFrom magrittr extract
#' @importFrom dplyr bind_rows mutate select case_when any_of arrange
#' @importFrom tidyr pivot_wider drop_na
#' @importFrom rlang sym
#'
.prepareData <- function(DGEdata = DGEdata,
                         x = x,
                         y = y,
                         ID = ID,
                         symbol = symbol,
                         logFC = logFC,
                         FDR = FDR,
                         logFCcutoff = logFCcutoff,
                         FDRcutoff = FDRcutoff){

    DGEdata |>
        magrittr::extract(c(x, y)) |>
        dplyr::bind_rows(.id = "contrast") |>
        dplyr::mutate(ID = !!rlang::sym(ID),
                      symbol = !!rlang::sym(symbol),
                      LogFC = !!rlang::sym(logFC),
                      FDR = !!rlang::sym(FDR)) |>
        dplyr::mutate(Direction = dplyr::case_when(LogFC > 0 ~ "Up",
                                                   LogFC < 0 ~ "Down")) |>
        dplyr::select("contrast", "ID", "LogFC", "FDR", "Direction",
                      dplyr::any_of(c("symbol"))) |>
        dplyr::mutate(FDRpass = dplyr::case_when(abs(LogFC) > logFCcutoff &
                                                     FDR < FDRcutoff ~ TRUE,
                                                 .default = FALSE)) |>
        tidyr::pivot_wider(names_from = contrast,
                           names_glue = "{contrast} {.value}",
                           values_from = c(LogFC, FDR, FDRpass, Direction)) |>
        tidyr::drop_na(!dplyr::any_of(c("ID", "symbol"))) |>
        dplyr::mutate(Significant = dplyr::case_when(
            !!rlang::sym(paste(x, "FDRpass")) == FALSE &
                !!rlang::sym(paste(y, "FDRpass")) == FALSE ~
                "Not Significant",
            !!rlang::sym(paste(x, "FDRpass")) == TRUE &
                !!rlang::sym(paste(y, "FDRpass")) == FALSE ~
                paste("Significant in", x),
            !!rlang::sym(paste(x, "FDRpass")) == FALSE &
                !!rlang::sym(paste(y, "FDRpass")) == TRUE ~
                paste("Significant in", y),
            !!rlang::sym(paste(x, "FDRpass")) == TRUE &
                !!rlang::sym(paste(y, "FDRpass")) == TRUE ~
                "Significant in Both") |>
                factor(levels = c("Not Significant",
                                  "Significant in Both",
                                  paste("Significant in", x),
                                  paste("Significant in", y)))) |>
        dplyr::arrange(Significant)
}

#' @title Correlation test
#' @description Test the correlation between DGE contrasts
#' @keywords internal
#' @param DGEtibble A \link[tibble:tbl_df-class]{tibble} of DGE results
#' @return A numeric of the Pearson correlation
#' @importFrom dplyr select contains filter_all all_vars
#' @importFrom stats cor
#'
.testCor <- function(DGEtibble = DGEtibble){
    corMatrix <- DGEtibble |>
        dplyr::select(dplyr::contains("LogFC")) |>
        dplyr::filter_all(dplyr::all_vars(is.finite(.))) |>
        as.matrix()

    cor(x = corMatrix[, 1], y = corMatrix[, 2])
}

#' @title Summarize counts
#' @description Create a summary table counts for DGE contrast overlaps for
#'  shared (quadrants) and non-shared (lines) DEGs
#' @keywords internal
#' @inheritParams gg4way
#' @inheritParams .testCor
#' @return A \link[tibble:tbl_df-class]{tibble}
#' @importFrom dplyr filter mutate case_when bind_rows select
#' @importFrom rlang sym
#' @importFrom janitor tabyl
#' @importFrom tidyr crossing
#' @importFrom glue glue_data
#'
.totalCounts <- function(DGEtibble = DGEtibble,
                         x = x,
                         y = y,
                         logFCcutoff = logFCcutoff){
    quadCounts <- DGEtibble |>
        dplyr::filter(!!rlang::sym(paste(x, "FDRpass")) == TRUE &
                          !!rlang::sym(paste(y, "FDRpass")) == TRUE) |>
        dplyr::mutate(countGroup = dplyr::case_when(
            !!rlang::sym(paste(x, "LogFC")) > logFCcutoff &
                !!rlang::sym(paste(y, "LogFC")) > logFCcutoff ~ "topRight",
            !!rlang::sym(paste(x, "LogFC")) < -logFCcutoff &
                !!rlang::sym(paste(y, "LogFC")) > logFCcutoff ~ "topLeft",
            !!rlang::sym(paste(x, "LogFC")) > logFCcutoff &
                !!rlang::sym(paste(y, "LogFC")) < -logFCcutoff ~ "bottomRight",
            !!rlang::sym(paste(x, "LogFC")) < -logFCcutoff &
                !!rlang::sym(paste(y, "LogFC")) < -logFCcutoff ~ "bottomLeft") |>
                factor(c("topRight", "topLeft",
                         "bottomRight", "bottomLeft"))) |>
        janitor::tabyl(countGroup)

    stopifnot(!any(is.na(quadCounts$countGroup)))

    lineCounts <- DGEtibble |>
        dplyr::filter(!Significant %in% c("Not Significant",
                                          "Significant in Both")) |>
        dplyr::mutate(countGroup = dplyr::case_when(
            Significant == paste("Significant in", x) ~
                paste0("sigX ", "x", !!rlang::sym(paste(x, "Direction")),
                       " y", !!rlang::sym(paste(y, "Direction"))),
            Significant == paste("Significant in", y) ~
                paste0("sigY ", "x", !!rlang::sym(paste(x, "Direction")),
                       " y", !!rlang::sym(paste(y, "Direction")))) |>
                factor(levels = tidyr::crossing(Sig = c("X", "Y"),
                                                Direction1 = c("Up", "Down"),
                                                Direction2 = Direction1) |>
                           glue::glue_data("sig{Sig} \\
                                           x{Direction1} y{Direction2}"))) |>
        janitor::tabyl(countGroup)

    if(all(is.na(lineCounts$percent))){
        lineCounts <- lineCounts |>
            dplyr::mutate(percent = 0)
    }

    stopifnot(!any(is.na(lineCounts$percent)))

    quadCounts |>
        dplyr::bind_rows(lineCounts) |>
        dplyr::select(countGroup, n)
}

#' @title Prepare annotations
#' @description Prepare text annotations of sums for plotting
#' @keywords internal
#' @param totalTibble A \link[tibble:tbl_df-class]{tibble} of summarized counts
#' @inheritParams gg4way
#' @return A \link[tibble:tbl_df-class]{tibble}
#' @importFrom dplyr left_join mutate recode
#' @importFrom tibble tribble
#' @importFrom stats setNames
#'
.prepareAnnotations <- function(totalTibble = totalTibble,
                                colorKey = colorKey,
                                textNudge = textNudge){
    totalTibble |>
        dplyr::left_join(
            tibble::tribble(
                ~countGroup, ~x, ~y, ~color, ~hjust, ~vjust,
                "sigX xUp yUp", Inf, textNudge, colorKey[2], 1, 0,
                "sigX xUp yDown", Inf, -textNudge, colorKey[2], 1, 1,
                "sigX xDown yUp", -Inf, textNudge, colorKey[2], 0, 0,
                "sigX xDown yDown", -Inf, -textNudge, colorKey[2], 0, 1,
                "sigY xUp yUp", textNudge, Inf, colorKey[3], 0, 1,
                "sigY xUp yDown", textNudge, -Inf, colorKey[3], 0, 0,
                "sigY xDown yUp", -textNudge, Inf, colorKey[3], 1, 1,
                "sigY xDown yDown", -textNudge, -Inf, colorKey[3], 1, 0,
                "topRight", Inf, Inf, colorKey[4], 1, 1,
                "topLeft", -Inf, Inf, colorKey[4], 0, 1,
                "bottomRight", Inf, -Inf, colorKey[4], 1, 0,
                "bottomLeft", -Inf, -Inf, colorKey[4], 0, 0),
            by = "countGroup") |>
        dplyr::mutate(Significant = dplyr::recode(color,
                                                  !!!setNames(names(colorKey),
                                                              colorKey)))

}

#' @title Tidy axis labels
#' @description Process axis labels from contrast names
#' @keywords internal
#' @inheritParams gg4way
#' @return A \link[base]{call}
#' @importFrom stringr str_split str_pad
#' @importFrom rlang expr
#' @importFrom purrr flatten_chr
#'
.tidyLabel <- function(label = NULL,
                       sep = " vs ",
                       labelType = c("x", "y")) {
    stopifnot(!is.null(label) | is.character(label))

    labelParts <- label |>
        stringr::str_split(sep) |>
        purrr::flatten_chr()

    tidyPart1 <- labelParts[1]

    tidyPart2 <- labelParts[2]

    spacer1 <- stringr::str_pad("", nchar(tidyPart2), side = "right")

    spacer2 <- stringr::str_pad("", nchar(tidyPart1), side = "right")

    part1 <- rlang::expr(!!spacer2 ~ "Higher in" ~ !!tidyPart2 ~
                             "" %<->% "" ~
                             "Higher in" ~ !!tidyPart1 ~ !!spacer1)

    if (labelType == "x") {
        label <- rlang::expr(atop(!!part1, paste(!!label, " LogFC")))
    } else if (labelType == "y") {
        label <- rlang::expr(atop(paste(!!label, " LogFC"), !!part1))
    }

    return(label)
}

#' @title gg4way plot
#' @description Creates a 4way plot
#' @keywords internal
#' @inheritParams gg4way
#' @return A \link[ggplot2]{ggplot}
#' @import ggplot2
#' @importFrom stringr str_split str_pad
#' @importFrom rlang expr
#' @importFrom purrr flatten_chr
#' @importFrom scales alpha
#' @importFrom ggrepel geom_label_repel
#'
.plot4way <- function(DGEtibble = DGEtibble,
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
                      label = label){
    p1 <- DGEtibble |>
        dplyr::mutate(alpha = dplyr::case_match(Significant,
                                                "Not Significant" ~ 0.5,
                                                .default = 0.6)) |>
        ggplot2::ggplot(ggplot2::aes(x = DGEtibble |>
                                         purrr::pluck(paste(x, "LogFC")),
                                     y = DGEtibble |>
                                         purrr::pluck(paste(y, "LogFC")),
                                     color = Significant)) +
        ggplot2::geom_point(ggplot2::aes(alpha = alpha), size = 1) +
        ggplot2::geom_vline(xintercept = c(-logFCcutoff, logFCcutoff),
                            linetype = "dashed",
                            color = lineColor) +
        ggplot2::geom_hline(yintercept = c(-logFCcutoff, logFCcutoff),
                            linetype = "dashed",
                            color = lineColor) +
        ggplot2::geom_hline(yintercept = 0, color = lineColor) +
        ggplot2::geom_vline(xintercept = 0, color = lineColor) +
        ggplot2::geom_abline(intercept = 0, slope = 1, color = lineColor) +
        ggplot2::xlab(x |> .tidyLabel(labelType = "x", sep = sep)) +
        ggplot2::ylab(y |> .tidyLabel(labelType = "y", sep = sep)) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::scale_color_manual(
            values = colorKey,
            guide = ggplot2::guide_legend(
                title = glue::glue("r = {tidyCor}",
                                   tidyCor = round(corRes, digits = 2)),
                nrow = 2,
                override.aes = list(size = 2)),
            drop = FALSE) +
        ## Use geom_label to prevent text from getting clipped
        ggplot2::geom_label(data = textKey,
                            ggplot2::aes(x = x,
                                         y = y,
                                         color = Significant,
                                         label = n,
                                         hjust = hjust,
                                         vjust = vjust),
                            size = textSize,
                            alpha = 0.1,
                            fill = NA,
                            show.legend = FALSE,
                            label.padding = ggplot2::unit(0.1, "lines"),
                            label.size = NA) +
        ggplot2::scale_alpha(range = c(0.5, 0.6),
                             guide = "none")

    if (label[1] != FALSE) {
        if (label[1] == TRUE) {
            label <- DGEtibble |>
                dplyr::filter(Significant == "Significant in Both") |>
                dplyr::pull(symbol)
        }

        p1 <- p1 +
            ggrepel::geom_label_repel(
                data = DGEtibble |>
                    dplyr::mutate(symbol = dplyr::case_match(symbol,
                                                             label ~ symbol,
                                                             .default = "")),
                ggplot2::aes(label = symbol),
                arrow = grid::arrow(length = ggplot2::unit(0.01, "npc")),
                min.segment.length = ggplot2::unit(0, "npc"),
                fill = scales::alpha("white", 0.4),
                label.padding = ggplot2::unit(0.01, "lines"),
                label.size = NA,
                show.legend = FALSE,
                max.overlaps = Inf)
    }
    return(p1)
}

utils::globalVariables(c(
    ".", "contrast", "LogFC", "FDR", "FDRpass", "Direction",
    "countGroup", "Significant", "Direction1", "n",
    "color", "hjust", "vjust", "symbol", "atop",
    "N052611", "N061011", "N61311", "airway"))
