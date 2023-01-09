library(shiny)
library(rhandsontable)
library(data.table)
library(dplyr)

xy2well <- function(x, y) {
    sprintf("%s%02d", LETTERS[y], x)
}
# plate data
plate <- list()
# 96well plate
colwise <- 1:96
names(colwise) <- paste0(rep(LETTERS[1:8], times=12), rep(sprintf("%02d", 1:12), each=8))
rowwise <- 1:96
names(rowwise) <- paste0(rep(LETTERS[1:8], each=12), rep(sprintf("%02d", 1:12), times=8))
plate$w96 <-
    expand.grid(x = 1:12, y = 1:8) %>%
    mutate(well = xy2well(x, y), X = "", selected = TRUE) %>%
    mutate(colwise = colwise[well], rowwise = rowwise[well]) %>%
    as.data.table()

# 384well plate
colwise <- 1:384
names(colwise) <- paste0(rep(LETTERS[1:16], times=24), rep(sprintf("%02d", 1:24), each=16))
rowwise <- 1:384
names(rowwise) <- paste0(rep(LETTERS[1:16], each=24), rep(sprintf("%02d", 1:24), times=16))
plate$w384 <-
    expand.grid(x = 1:24, y = 1:16) %>%
    mutate(well = xy2well(x, y), X = "", selected = TRUE) %>%
    mutate(colwise = colwise[well], rowwise = rowwise[well]) %>%
    as.data.table()

channel_select <- function(df, plate_type, ch_row, ch_col) {
    if (plate_type == "96") {
        tmp_fun <- function(.x, .y) expand.grid(x = seq(.x, min(12, .x + ch_col - 1)), y = seq(.y, min(8, .y + ch_row - 1)))
        index <- 
            df %>%
            filter(selected_) %>% 
            {purrr::map2(.$x, .$y, tmp_fun)} %>% 
            bind_rows() %>% 
            {(.$x - 1) * 8 + .$y}
        point_index <- function(x, y) {(x - 1) * 8 + y}
    } else {
        tmp_fun <- function(.x, .y) expand.grid(x = seq(.x, min(24, .x + 2*(ch_col - 1)), 2), y = seq(.y, min(16, .y + 2*(ch_row - 1)), 2))
        index <- 
            df %>% 
            filter(selected_) %>% 
            {purrr::map2(.$x, .$y, tmp_fun)} %>% 
            bind_rows() %>% 
            {(.$x - 1) * 16 + .$y}
        point_index <- function(x, y) {(x - 1) * 16 + y}
    }
    df %>% 
        mutate(.t = point_index(x, y)) %>% 
        mutate(selected_ = (.t %in% index)) %>% 
        select(-.t)
} 

ui <- fluidPage(
    titlePanel("96 or 384 well plate"),
    fluidRow(
        column(
            width = 8,
                wellPanel(
            fluidRow(
                    column(
                        width = 4,
                        radioButtons("plateType", "plate type", c("96", "384")),
                        radioButtons("wellOrder", "well order", c("colwise", "rowwise")),
                    ),
                    column(
                        width = 4,
                        numericInput("channelRow", "channel (row)", value = 1L, min = 1L, max = 16L, step = 1L),
                        numericInput("channelCol", "channel (col)", value = 1L, min = 1L, max = 24L, step = 1L)
                    ),
                    column(
                        width = 4,
                        selectizeInput("colNames", "column names", choices = c("X"),
                                       selected = "X",
                                       multiple = TRUE, options = list(create = TRUE))
                    )
                )
            ),
            plotOutput("plot", click = "click", brush = "brush"),
            tags$script(HTML("
                $('#plot').mousedown(function(e) {
                    var parentOffset = $(this).offset();
                    var relX = e.pageX - parentOffset.left;
                    var relY = e.pageY - parentOffset.top;
                    Shiny.setInputValue('x1', relX);
                    Shiny.setInputValue('y1', relY);
                }).mouseup(function(e) {
                    var parentOffset = $(this).offset();
                    var relX = e.pageX - parentOffset.left;
                    var relY = e.pageY - parentOffset.top;
                    Shiny.setInputValue('x2', relX);
                    Shiny.setInputValue('y2', relY);
                    Shiny.setInputValue('action', Math.random());
                });
            ")),
        ),
        column(
            width = 4,
            rHandsontableOutput("table", height = 650)
        )
    )
)

server <- function(input, output, session) {
    data <- reactiveValues(values = plate)
    observeEvent(input$plateType, {
        updateSelectizeInput(inputId = "colNames", choices = "X", selected = "X")
        data$values <- plate
    })
    observeEvent(input$colNames, {
        req(input$colNames)
        is_preserved_name <- input$colNames %in% c("x", "y", "rowwise", "colwise", "selected", "selected_")
        is_invalid_name <- !stringr::str_detect(input$colNames, "^[A-Za-z]\\w*$")
        if (any(is_preserved_name) || any(is_invalid_name)) {
            showNotification("Invalid column name!", type = "warning")
            .t <- input$colNames[!(is_preserved_name | is_invalid_name)]
            updateSelectizeInput(inputId = "colNames", choices = .t, selected = .t)
            return(NULL)
        }
        if (input$plateType == "96") {
            d <- data$values$w96
        } else {
            d <- data$values$w384
        }
        new_colnames <- setdiff(input$colNames, colnames(d))
        for (nc in new_colnames) {
            d[[nc]] <- ""
        }
        if (input$plateType == "96") {
            data$values$w96 <- d
        } else {
            data$values$w384 <- d
        }
        print(data$values$w96)
                                 
    })

    # table output
    output$table <- renderRHandsontable({

        data$values %>%
            {if (input$plateType == "96") purrr::pluck(., "w96") else purrr::pluck(., "w384")} %>%
            {if (input$wellOrder == "colwise") dplyr::arrange(., colwise) else dplyr::arrange(., rowwise)} %>%
            filter(selected) %>%
            dplyr::select(well, !!input$colNames) %>%
            rhandsontable()
    })

    # input from hansontable
    observeEvent(eventExpr = input$table, {
        tmp <- hot_to_r(input$table)
        if (input$plateType == "96") {
            tmp <-
                data$values$w96 %>%
                left_join(tmp, by = "well", suffix = c("", ".new"))
        } else {
            tmp <-
                data$values$w384 %>%
                left_join(tmp, by = "well", suffix = c("", ".new"))
        }
        tmp_fun <- function(.x, .xnew) {
            x <- if_else(is.na(.xnew) | .xnew == "", .x, .xnew)
            x
        }
        for (cn in input$colNames) {
            cn2 <- paste(cn, "new", sep = ".")
            print(head(tmp))
            tmp[[cn]] <- tmp_fun(tmp[[cn]], tmp[[cn2]])
            tmp <- tmp %>% select(-(!!cn2))
        }
        if (input$plateType == "96") {
            data$values$w96 <- tmp
        } else {
            data$values$w384 <- tmp
        }
    })


    output$plot <- renderPlot({
        session$resetBrush("brush")
        par(mar = c(0, 0, 0, 0))
        if (input$plateType == "96") {
            data <- data$values$w96
            plot(data$x, data$y, axes=FALSE, xlab=NA, ylab=NA, cex=3, asp=1, xlim = c(-1, 14), ylim = c(10, -1), lwd = data$selected + 1)
            rect(-1, -1, 13, 9)
            text(x = 0, y = 1:8, LETTERS[1:8])
            text(x = 1:12, y = 0, as.character(1:12))
        } else {
            data <- data$values$w384
            plot(data$x, data$y, axes=FALSE, xlab=NA, ylab=NA, cex=2, asp=1, xlim = c(-1, 26), ylim = c(18, -1), lwd = data$selected + 1)
            rect(-1, -1, 25, 17)
            text(x = 0, y = 1:16, LETTERS[1:16])
            text(x = 1:24, y = 0, as.character(1:24))
        }
    })

    observeEvent(input$action, {
        req(input$click)
        if ((input$x1 - input$x2)^2 + (input$y1 - input$y2)^2 < 5) {
            cat("click\n")
            if (input$click$x < 0.5 || input$click$y < 0.5) {
                if (input$plateType == "96") {
                    data$values$w96 <-
                        data$values$w96 %>%
                        mutate(selected = FALSE)
                } else {
                    data$values$w384 <-
                        data$values$w384 %>%
                        mutate(selected = FALSE)
                }
            } else {
                updated_data <- if (input$plateType == "96") data$values$w96 else data$values$w384
                updated_data <-
                    updated_data %>%
                    nearPoints(input$click, xvar = "x", yvar = "y", threshold = 10, allRows = TRUE) %>%
                    {if (all(.$selected)) mutate(., selected = FALSE) else .} %>%
                    #mutate(., selected = selected | selected_) %>%
                    channel_select(plate_type = input$plateType, ch_col = input$channelCol, ch_row = input$channelRow) %>% 
                    mutate(., selected = xor(selected, selected_)) %>%
                    select(-selected_)
                if (input$plateType == "96") {
                    data$values$w96 <- updated_data
                } else {
                    data$values$w384 <- updated_data
                }
            }
        } else {
            cat("brush\n")
            # update data
            updated_data <- if (input$plateType == "96") data$values$w96 else data$values$w384
            updated_data <-
                updated_data %>%
                brushedPoints(input$brush, xvar = "x", yvar = "y", allRows = TRUE) %>%
                {if (all(.$selected)) mutate(., selected = FALSE) else .} %>%
                #mutate(selected = selected | selected_) %>%
                channel_select(plate_type = input$plateType, ch_col = input$channelCol, ch_row = input$channelRow) %>% 
                mutate(., selected = xor(selected, selected_)) %>%
                select(-selected_)
            if (input$plateType == "96") {
                data$values$w96 <- updated_data
            } else {
                data$values$w384 <- updated_data
            }

        }
    })
}

shinyApp(ui, server)
