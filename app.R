library(shiny)
library(data.table)
library(limma)
library(caTools)
library(rtracklayer)

# Allow uploads up to 2 GB (ibed and BigWig files can be large)
options(shiny.maxRequestSize = 2 * 1024^3)  # 2 GB
options(shiny.maxRequestSize = 500 * 1024^2) # 500 MB

# Default color palette
plotcols <- c("blue", "red", "green", "black", "cyan", "purple", "orange", "gray")

# ─── Helper: gene parking (avoid overlapping gene labels) ────────────────────
parking <- function(left, right) {
  y <- rep(-1, length(right))
  lengths <- right - left
  for (i in order(lengths, decreasing = TRUE)) {
    otherleft  <- left;  otherleft[i]  <- NA
    otherright <- right; otherright[i] <- NA
    placed <- FALSE; y[i] <- 0
    while (!placed) {
      placed <- sum(
        (right[i] > otherleft[y == y[i]]) & (left[i] < otherright[y == y[i]]),
        na.rm = TRUE
      ) == 0
      if (!placed) y[i] <- y[i] + 1
    }
  }
  return(y)
}

# ─── Gene track plot ──────────────────────────────────────────────────────────
plot_genes <- function(genome_df, x.min, x.max) {
  y_plot <- parking(genome_df$Start, genome_df$End)
  par(mar = c(0, 2, 0, 2) + 0.1)
  plot(c(x.min, x.max), c(1, -max(y_plot) - 0.5),
       col = "white", ylab = "", xlab = "",
       fg = "white", col.axis = "white",
       xaxs = "i", yaxs = "i")
  arrowHeads <- pretty(x.min:x.max, n = 50)
  for (i in seq_len(nrow(genome_df))) {
    x <- c(
      genome_df$Start[i],
      arrowHeads[arrowHeads > genome_df$Start[i] & arrowHeads < genome_df$End[i]],
      genome_df$End[i]
    )
    if (length(x) < 2) next
    if (genome_df$Strand[i] == "-") {
      arrows(x[2:length(x)], -y_plot[i], x[1:(length(x)-1)], col = "blue", length = 0.08)
    } else {
      arrows(x[1:(length(x)-1)], -y_plot[i], x[2:length(x)], col = "blue", length = 0.08)
    }
    text(genome_df$Start[i], -y_plot[i] + 0.4, adj = 0, labels = genome_df$Name[i])
  }
}

# ─── Epigenomic track plot ─────────────────────────────────────────────────────
plot_tracks <- function(epigenome, chr, x.min, x.max, track_names) {
  tmps     <- list()
  plotlim  <- 0
  for (i in seq_along(track_names)) {
    tmp <- epigenome[[track_names[i]]]$track
    tmp <- tmp[as.character(seqnames(tmp)) == chr &
                 start(tmp) >= x.min & end(tmp) <= x.max, ]
    tmps[[i]] <- tmp
    if (length(tmp) > 0) plotlim <- max(plotlim, score(tmp), na.rm = TRUE)
  }
  for (i in seq_along(tmps)) {
    par(mar = c(0, 2, 0, 2) + 0.1)
    plot.new()
    plot.window(xlim = c(x.min, x.max), ylim = c(0, max(plotlim, 1)),
                xaxs = "i", yaxs = "i")
    if (length(tmps[[i]]) > 0) {
      segments(start(tmps[[i]]), 0, end(tmps[[i]]), score(tmps[[i]]),
               col = epigenome[[track_names[i]]]$color)
    }
    title(main = track_names[i],
          col.main = epigenome[[track_names[i]]]$color, adj = 0, line = -2)
  }
}

# ─── Group epigenomic tracks by level ─────────────────────────────────────────
get_track_levels <- function(epigenome) {
  levels <- NULL
  for (i in seq_along(epigenome)) {
    if (epigenome[[i]]$level > 0) levels <- c(levels, epigenome[[i]]$level)
  }
  if (is.null(levels)) return(list())
  levels <- unique(sort(levels))
  plottracklevels <- list()
  for (l in seq_along(levels)) {
    for (i in seq_along(epigenome)) {
      if (epigenome[[i]]$level == levels[l]) {
        if (length(plottracklevels) < l) {
          plottracklevels[[l]] <- names(epigenome)[i]
        } else {
          plottracklevels[[l]] <- c(plottracklevels[[l]], names(epigenome)[i])
        }
      }
    }
  }
  return(plottracklevels)
}

# ─── Core CHiC plot function ──────────────────────────────────────────────────
plot_chic <- function(ibed, gene, min.plot, max.plot, win.plot,
                      settings, ints, genome_df, epigenome) {

  table <- ibed[ibed$Bait_name == gene, ]
  if (nrow(table) == 0) stop(paste("No interactions found for bait gene:", gene))
  table <- table[order(table$start_OE), ]
  ref   <- table$start_Bait[1]

  # Quantile-normalise read columns across all conditions
  usecols <- unlist(lapply(settings, `[[`, "cols"))
  n <- normalizeBetweenArrays(as.matrix(table[, usecols]))
  for (i in seq_along(usecols)) table[, usecols[i]] <- n[, i]

  # Define plot window
  if (!is.na(min.plot) && !is.na(max.plot)) {
    table <- table[table$start_OE >= min.plot & table$end_OE <= max.plot, ]
  } else if (!is.na(win.plot) && !is.null(win.plot)) {
    table <- table[abs(table$start_OE - ref) <= as.integer(win.plot), ]
  } else {
    stop("Provide either a window size or start/end coordinates.")
  }
  if (nrow(table) == 0) stop("No interactions in the specified region.")

  # Count sub-panels
  n.plots <- 1  # CHiC profile
  plottracklevels <- get_track_levels(epigenome)
  n.plots <- n.plots + length(plottracklevels)
  if (!is.null(genome_df)) n.plots <- n.plots + 1
  has_ints <- length(ints) > 0 && any(sapply(ints, `[[`, "plot") == "1")
  if (has_ints) n.plots <- n.plots + 1

  if (n.plots > 1) {
    layout(matrix(seq_len(n.plots), ncol = 1, nrow = n.plots),
           heights = c(20, rep(10, n.plots - 1)))
  }

  # ── Panel 1: CHiC interaction profile ───────────────────────────────────────
  x.min <- min(table$start_OE)
  x.max <- max(table$end_OE)

  avg_list <- lapply(settings, function(s) {
    tab <- as.matrix(table[, s[["cols"]], drop = FALSE])
    caTools::runmean(rowMeans(tab), k = 10, endrule = "mean", align = "center")
  })
  y.max <- max(unlist(avg_list), na.rm = TRUE)

  main.text <- paste0(gene, "  –  ", table$chr_Bait[1], ": ",
                      format(x.min, big.mark = ","), "–", format(x.max, big.mark = ","))
  par(mar = c(4, 4, 3, 2) + 0.1)
  plot(table$start_OE, avg_list[[1]], type = "l",
       col  = settings[[1]]$color, lwd = 2,
       ylab = "Running mean QN signal", xlab = "Genomic coordinate",
       xaxs = "i", yaxs = "i",
       ylim = c(0, y.max * 1.1), main = main.text)

  if (length(settings) > 1) {
    for (i in 2:length(settings)) {
      lines(table$start_OE, avg_list[[i]], col = settings[[i]]$color, lwd = 2)
    }
  }
  legend("topleft", inset = 0.02,
         legend = sapply(settings, `[[`, "name"),
         col    = sapply(settings, `[[`, "color"),
         lty = 1, lwd = 2, title = "Condition", bg = "white")
  abline(v = ref, col = "black", lty = 2)

  # Overlay called interactions as points
  if (has_ints) {
    keep <- Filter(function(x) x$plot == "1", ints)
    for (i in seq_along(keep)) {
      inttab <- keep[[i]]$ibed
      inttab <- inttab[inttab$Bait_name == gene &
                         inttab$start_OE >= x.min &
                         inttab$end_OE   <= x.max, ]
      if (nrow(inttab) > 0) {
        points(inttab$start_OE, inttab$N, col = keep[[i]]$color, pch = 19, cex = 0.8)
      }
    }
    # Interaction legend panel
    par(mar = c(0, 2, 0, 2) + 0.1)
    plot.new()
    legend("center", inset = 0.05,
           legend = sapply(keep, `[[`, "name"),
           col    = sapply(keep, `[[`, "color"),
           pch = 19, title = "Called interactions", bg = "white")
  }

  # ── Gene track ───────────────────────────────────────────────────────────────
  if (!is.null(genome_df)) {
    g <- genome_df[genome_df$Chr == table$chr_Bait[1] &
                     genome_df$Start >= x.min &
                     genome_df$End   <= x.max, ]
    if (nrow(g) > 0) {
      plot_genes(g, x.min, x.max)
    } else {
      par(mar = c(0, 2, 0, 2) + 0.1)
      plot.new()
      text(0.5, 0.5, "No genes in region", col = "gray50")
    }
  }

  # ── Epigenomic tracks ────────────────────────────────────────────────────────
  if (length(plottracklevels) > 0) {
    for (l in seq_along(plottracklevels)) {
      plot_tracks(epigenome, table$chr_Bait[1], x.min, x.max, plottracklevels[[l]])
    }
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-size: 13px; }
    h4   { margin-top: 12px; margin-bottom: 4px; color: #2c3e50; font-size: 14px; }
    h5   { margin-top: 8px;  margin-bottom: 4px; color: #555;    font-size: 12px; }
    .sidebar-panel { background: #f8f9fa; padding: 10px; border-radius: 6px; }
    .shiny-input-container { margin-bottom: 4px; }
    hr { margin: 10px 0; }
    .btn-block { width: 100%; margin-top: 6px; }
  "))),

  titlePanel("ChiCMaxima – Promoter Capture HiC Browser"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      class = "sidebar-panel",

      # ── Data loading ──────────────────────────────────────────────────────
      h4("Data Files"),
      fileInput("ibed_file", "ibed interaction file (.ibed)",
                accept = c(".ibed", "text/plain"), width = "100%"),

      conditionalPanel("output.ibed_loaded == true",
        h5("Assign columns to conditions (0 = exclude):"),
        uiOutput("conditions_ui"),
        actionButton("set_conditions_btn", "Apply Conditions",
                     class = "btn-primary btn-sm btn-block")
      ),

      conditionalPanel("output.conditions_set == true",
        h5("Condition names & colours:"),
        uiOutput("plot_settings_ui")
      ),

      hr(),

      # ── Navigation ────────────────────────────────────────────────────────
      h4("Navigation"),
      textInput("bait",     "Bait gene name",          value = "Myc"),
      fluidRow(
        column(6, numericInput("min_plot", "Start (bp)", value = NA)),
        column(6, numericInput("max_plot", "End (bp)",   value = NA))
      ),
      numericInput("win_plot", "Window (bp) — used when start/end are NA",
                   value = 1000000, step = 100000),
      actionButton("plot_btn", "Plot", class = "btn-success btn-block"),

      hr(),

      # ── Additional tracks ─────────────────────────────────────────────────
      h4("Additional Tracks"),
      fileInput("genes_file", "Gene annotation (.txt)",
                accept = ".txt", width = "100%"),
      fileInput("epi_files",  "Epigenomic tracks (.bw / .bedGraph)",
                accept = c(".bw", ".bedGraph"), multiple = TRUE, width = "100%"),

      conditionalPanel("output.epi_loaded == true",
        h5("Track display level (0 = hidden, shared level = same row):"),
        uiOutput("track_levels_ui"),
        actionButton("apply_tracks_btn", "Apply Track Levels",
                     class = "btn-sm btn-block")
      ),

      fileInput("int_files", "Called interactions (.ibed)",
                accept = ".ibed", multiple = TRUE, width = "100%"),

      conditionalPanel("output.ints_loaded == true",
        h5("Interaction overlays (show / name / colour):"),
        uiOutput("int_settings_ui"),
        actionButton("apply_ints_btn", "Apply Interaction Settings",
                     class = "btn-sm btn-block")
      ),

      hr(),
      downloadButton("save_fig", "Save Figure (PDF)", class = "btn-block")
    ),

    mainPanel(
      width = 9,
      plotOutput("chic_plot", height = "680px"),
      br(),
      wellPanel(verbatimTextOutput("status_msg"))
    )
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# Server
# ══════════════════════════════════════════════════════════════════════════════

server <- function(input, output, session) {

  rv <- reactiveValues(
    ibed      = NULL,
    N         = NULL,
    settings  = NULL,
    ints      = list(),
    genome_df = NULL,
    epigenome = list(),
    status    = "Ready. Load an .ibed file to begin."
  )

  # ── Load ibed ──────────────────────────────────────────────────────────────
  observeEvent(input$ibed_file, {
    req(input$ibed_file)
    tryCatch({
      df <- as.data.frame(
        fread(input$ibed_file$datapath, header = TRUE, fill = TRUE),
        stringsAsFactors = FALSE
      )
      cols <- names(df)
      if (length(cols) < 11) {
        rv$status <- "Error: ibed file must have ≥11 columns (10 standard + read counts)."
        return()
      }
      rv$ibed     <- df
      rv$N        <- cols[-(1:10)]
      rv$settings <- NULL
      rv$status   <- sprintf("Loaded ibed: %d interactions | %d sample columns: %s",
                             nrow(df), length(rv$N), paste(rv$N, collapse = ", "))
    }, error = function(e) rv$status <- paste("Error loading ibed:", e$message))
  })

  output$ibed_loaded <- reactive({ !is.null(rv$ibed) })
  outputOptions(output, "ibed_loaded", suspendWhenHidden = FALSE)

  output$conditions_set <- reactive({ !is.null(rv$settings) })
  outputOptions(output, "conditions_set", suspendWhenHidden = FALSE)

  # ── Dynamic condition-assignment UI ───────────────────────────────────────
  output$conditions_ui <- renderUI({
    req(rv$N)
    tagList(lapply(seq_along(rv$N), function(i) {
      fluidRow(
        column(8, tags$small(rv$N[i])),
        column(4, numericInput(paste0("cond_", i), label = NULL,
                               value = 0, min = 0, max = 20, step = 1, width = "65px"))
      )
    }))
  })

  observeEvent(input$set_conditions_btn, {
    req(rv$N)
    assign <- setNames(
      sapply(seq_along(rv$N), function(i) {
        v <- input[[paste0("cond_", i)]]; if (is.null(v)) 0L else as.integer(v)
      }),
      rv$N
    )
    keep   <- assign[!is.na(assign) & assign > 0]
    if (length(keep) == 0) {
      rv$status <- "No conditions set. Assign at least one column a condition number > 0."
      return()
    }
    levels <- unique(sort(keep))
    rv$settings <- lapply(seq_along(levels), function(l) {
      hits <- which(assign == levels[l])
      list(
        name  = sub("\\..*", "", names(hits)[1]),
        N     = names(hits),
        cols  = hits + 10L,
        color = plotcols[l %% length(plotcols) + 1]
      )
    })
    rv$status <- sprintf("%d condition(s) configured.", length(rv$settings))
  })

  # ── Condition name/colour UI ──────────────────────────────────────────────
  output$plot_settings_ui <- renderUI({
    req(rv$settings)
    tagList(
      lapply(seq_along(rv$settings), function(i) {
        fluidRow(
          column(6, textInput(paste0("cname_", i), label = NULL,
                              value = rv$settings[[i]]$name,
                              placeholder = "Name")),
          column(6, textInput(paste0("ccol_",  i), label = NULL,
                              value = rv$settings[[i]]$color,
                              placeholder = "Colour"))
        )
      }),
      actionButton("apply_plot_settings_btn", "Apply Names & Colours",
                   class = "btn-sm btn-block")
    )
  })

  observeEvent(input$apply_plot_settings_btn, {
    req(rv$settings)
    for (i in seq_along(rv$settings)) {
      rv$settings[[i]]$name  <- input[[paste0("cname_", i)]]
      rv$settings[[i]]$color <- input[[paste0("ccol_",  i)]]
    }
    rv$status <- "Condition names and colours updated."
  })

  # ── Load gene annotation ──────────────────────────────────────────────────
  observeEvent(input$genes_file, {
    req(input$genes_file)
    tryCatch({
      rv$genome_df <- read.table(input$genes_file$datapath,
                                 header = TRUE, stringsAsFactors = FALSE)
      rv$status <- paste("Gene annotation loaded:", nrow(rv$genome_df), "genes.")
    }, error = function(e) rv$status <- paste("Error loading gene file:", e$message))
  })

  # ── Load epigenomic tracks ────────────────────────────────────────────────
  observeEvent(input$epi_files, {
    req(input$epi_files)
    tryCatch({
      epi <- list()
      for (i in seq_len(nrow(input$epi_files))) {
        tname <- sub("\\.[^.]+$", "", input$epi_files$name[i])
        epi[[tname]] <- list(
          level = i,
          color = plotcols[i %% length(plotcols) + 1],
          track = import(input$epi_files$datapath[i])
        )
      }
      rv$epigenome <- epi
      rv$status    <- paste("Loaded", length(epi), "epigenomic track(s).")
    }, error = function(e) rv$status <- paste("Error loading tracks:", e$message))
  })

  output$epi_loaded <- reactive({ length(rv$epigenome) > 0 })
  outputOptions(output, "epi_loaded", suspendWhenHidden = FALSE)

  output$track_levels_ui <- renderUI({
    req(length(rv$epigenome) > 0)
    tagList(lapply(seq_along(rv$epigenome), function(i) {
      fluidRow(
        column(8, tags$small(names(rv$epigenome)[i])),
        column(4, numericInput(paste0("trklev_", i), label = NULL,
                               value = rv$epigenome[[i]]$level,
                               min = 0, max = 20, step = 1, width = "65px"))
      )
    }))
  })

  observeEvent(input$apply_tracks_btn, {
    req(length(rv$epigenome) > 0)
    for (i in seq_along(rv$epigenome)) {
      lev <- input[[paste0("trklev_", i)]]
      if (!is.null(lev)) rv$epigenome[[i]]$level <- as.integer(lev)
    }
    rv$status <- "Track levels updated."
  })

  # ── Load interaction overlays ─────────────────────────────────────────────
  observeEvent(input$int_files, {
    req(input$int_files)
    tryCatch({
      ints <- lapply(seq_len(nrow(input$int_files)), function(i) {
        list(
          name  = sub("\\.[^.]+$", "", input$int_files$name[i]),
          color = plotcols[i %% length(plotcols) + 1],
          plot  = "1",
          ibed  = read.table(input$int_files$datapath[i],
                             header = TRUE, stringsAsFactors = FALSE)
        )
      })
      rv$ints  <- ints
      rv$status <- paste("Loaded", length(ints), "interaction set(s).")
    }, error = function(e) rv$status <- paste("Error loading interactions:", e$message))
  })

  output$ints_loaded <- reactive({ length(rv$ints) > 0 })
  outputOptions(output, "ints_loaded", suspendWhenHidden = FALSE)

  output$int_settings_ui <- renderUI({
    req(length(rv$ints) > 0)
    tagList(lapply(seq_along(rv$ints), function(i) {
      fluidRow(
        column(2, checkboxInput(paste0("intplot_", i), label = NULL, value = TRUE)),
        column(5, textInput(paste0("intname_", i), label = NULL,
                            value = rv$ints[[i]]$name)),
        column(5, textInput(paste0("intcol_",  i), label = NULL,
                            value = rv$ints[[i]]$color))
      )
    }))
  })

  observeEvent(input$apply_ints_btn, {
    req(length(rv$ints) > 0)
    for (i in seq_along(rv$ints)) {
      rv$ints[[i]]$name  <- input[[paste0("intname_", i)]]
      rv$ints[[i]]$color <- input[[paste0("intcol_",  i)]]
      rv$ints[[i]]$plot  <- if (isTRUE(input[[paste0("intplot_", i)]])) "1" else "0"
    }
    rv$status <- "Interaction settings updated."
  })

  # ── Main plot ─────────────────────────────────────────────────────────────
  plot_params <- eventReactive(input$plot_btn, {
    req(rv$ibed, rv$settings)
    list(
      ibed      = rv$ibed,
      gene      = input$bait,
      min.plot  = input$min_plot,
      max.plot  = input$max_plot,
      win.plot  = input$win_plot,
      settings  = rv$settings,
      ints      = rv$ints,
      genome_df = rv$genome_df,
      epigenome = rv$epigenome
    )
  })

  output$chic_plot <- renderPlot({
    pd <- plot_params()
    tryCatch(
      plot_chic(pd$ibed, pd$gene, pd$min.plot, pd$max.plot, pd$win.plot,
                pd$settings, pd$ints, pd$genome_df, pd$epigenome),
      error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Plot error:\n", e$message),
             cex = 1.1, col = "firebrick", adj = 0.5)
      }
    )
  })

  output$status_msg <- renderText({ rv$status })

  # ── PDF download ──────────────────────────────────────────────────────────
  output$save_fig <- downloadHandler(
    filename = function() {
      paste0("ChiCMaxima_", input$bait, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pd <- plot_params()
      pdf(file, width = 14, height = 10)
      tryCatch(
        plot_chic(pd$ibed, pd$gene, pd$min.plot, pd$max.plot, pd$win.plot,
                  pd$settings, pd$ints, pd$genome_df, pd$epigenome),
        error = function(e) {
          plot.new()
          text(0.5, 0.5, paste("Plot error:", e$message), cex = 1.2, col = "firebrick")
        }
      )
      dev.off()
    }
  )
}

shinyApp(ui, server)
