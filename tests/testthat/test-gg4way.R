test_that("example data works", {
    data("airwayFit")
    expect_true(is(airwayFit , "MArrayLM"))
})

test_that("correlation works", {
    data("airwayFit")
    airwayFit |>
        gg4way(x = "N61311 vs N052611",
               y = "N061011 vs N052611") |>
        getCor() |>
        expect_equal(0.42,
                     tolerance = 3e-2)
})

test_that("shared table works", {
    data("airwayFit")
    airwayFit |>
        gg4way(x = "N61311 vs N052611",
               y = "N061011 vs N052611") |>
        getShared() |>
        expect_s3_class("tbl_df")
})

test_that("totals tabyl works", {
    data("airwayFit")
    airwayFit |>
        gg4way(x = "N61311 vs N052611",
               y = "N061011 vs N052611") |>
        getTotals() |>
        expect_s3_class("tabyl")
})

test_that("axis labels work", {
    "N61311 vs N052611" |>
        .tidyLabel(labelType = "x") |>
        expect_type("language")

    "N061011 vs N052611" |>
        .tidyLabel(labelType = "y") |>
        expect_type("language")
})

test_that("plot works", {
    data("airwayFit")
    airwayFit |>
        gg4way(x = "N61311 vs N052611",
               y = "N061011 vs N052611") %T>%
        expect_s3_class("ggplot") |>
        expect_doppelganger("plot", fig = _)
})
