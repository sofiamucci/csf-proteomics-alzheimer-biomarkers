# Interactive viewer for the CSF proteome reanalysis (PXD016278).
# Runs in the browser with Shinylive (webR); data from scripts/10_app_data.R.

library(shiny)
library(bslib)
library(ggplot2)

d <- readRDS("app_data.rds")

GRP_COLS  <- c(control = "#2a78d6", AD = "#eb6834")
TIER_COLS <- c("barrier-robust" = "#eb6834", "main model only" = "#2a78d6",
               "not significant" = "grey75")
MAIN      <- "All (main model)"
COHORTS   <- c(MAIN, "Sweden", "Magdeburg", "Berlin")

prot <- d$proteins
prot$label <- make.unique(prot$gene)
choices <- setNames(prot$protein_group, prot$label)[order(prot$label)]

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank()))

uniprot_url <- function(pg) {
  paste0("https://www.uniprot.org/uniprotkb/", sub("[;-].*", "", pg), "/entry")
}

# ---- UI -----------------------------------------------------------------------
ui <- page_sidebar(
  title = "CSF proteome in Alzheimer's disease: interactive viewer",
  theme = bs_theme(version = 5, primary = "#2a78d6"),
  sidebar = sidebar(
    width = 300,
    selectInput("cohort", "Cohort", COHORTS),
    sliderInput("fdr", "FDR threshold", min = 0.01, max = 0.2, value = 0.05, step = 0.01),
    checkboxInput("core_only", "Highlight the 50-protein core", FALSE),
    selectizeInput("protein", "Protein", choices = choices,
                   selected = prot$protein_group[prot$label == "YWHAZ"][1]),
    helpText("Click a point in the volcano plot to select a protein."),
    hr(),
    helpText(
      "Main model: proDA, ~ group + site + age + sex (88 AD, 109 controls).",
      "Cohorts: the same model within each site, without site.",
      "Cohorts have far less power; read the direction and size of effects,",
      "not significance. Kiel has no AD cases."
    ),
    tags$a("Repository and methods",
           href = "https://github.com/sofiamucci/csf-proteomics-alzheimer-biomarkers",
           target = "_blank")
  ),
  layout_columns(
    col_widths = c(7, 5),
    card(card_header(textOutput("volcano_title")),
         plotOutput("volcano", click = "volcano_click", height = "420px")),
    card(card_header("Top protein groups"),
         tableOutput("top"))
  ),
  card(
    card_header(uiOutput("protein_header")),
    layout_columns(
      col_widths = c(4, 4, 4),
      plotOutput("forest", height = "300px"),
      plotOutput("distribution", height = "300px"),
      plotOutput("tau", height = "300px")
    )
  )
)

# ---- Server -------------------------------------------------------------------
server <- function(input, output, session) {

  res <- reactive({
    r <- d$results[d$results$cohort == input$cohort, ]
    r <- merge(r, prot[, c("protein_group", "label", "tier", "core")], by = "protein_group")
    r$significant <- !is.na(r$p_adj) & r$p_adj < input$fdr
    r
  })

  observeEvent(input$volcano_click, {
    r <- res()
    r$neglog <- -log10(pmax(r$p_adj, 1e-16))
    hit <- nearPoints(r, input$volcano_click, xvar = "log2FC", yvar = "neglog",
                      maxpoints = 1, threshold = 10)
    if (nrow(hit) == 1) updateSelectizeInput(session, "protein", selected = hit$protein_group)
  })

  output$volcano_title <- renderText({
    r <- res()
    sprintf("%s: %d of %d protein groups at FDR < %.2f (%d higher, %d lower in AD)",
            input$cohort, sum(r$significant), nrow(r), input$fdr,
            sum(r$significant & r$log2FC > 0), sum(r$significant & r$log2FC < 0))
  })

  output$volcano <- renderPlot({
    r <- res()
    r$neglog <- -log10(pmax(r$p_adj, 1e-16))
    r$status <- ifelse(!r$significant, "not significant",
                       ifelse(r$log2FC > 0, "higher in AD", "lower in AD"))
    sel <- r[r$protein_group == input$protein, ]
    p <- ggplot(r, aes(log2FC, neglog)) +
      geom_hline(yintercept = -log10(input$fdr), linetype = "dashed", colour = "grey50") +
      geom_vline(xintercept = 0, colour = "grey85") +
      geom_point(aes(colour = status), size = 1.5, alpha = 0.7) +
      scale_colour_manual(values = c("higher in AD" = "#eb6834", "lower in AD" = "#2a78d6",
                                     "not significant" = "grey75"), name = NULL) +
      labs(x = "log2 fold change, AD vs control", y = "−log10 adjusted p-value") +
      theme(legend.position = "top")
    if (input$core_only) {
      p <- p + geom_point(data = r[r$core, ], shape = 21, size = 2.6,
                          colour = "black", stroke = 0.4)
    }
    p + geom_point(data = sel, shape = 21, size = 4.5, colour = "black", stroke = 1.2) +
      geom_label(data = sel, aes(label = label), vjust = -0.8, size = 3.5)
  })

  output$top <- renderTable({
    r <- res()
    if (input$core_only) r <- r[r$core, ]
    r <- r[order(r$p_adj), ][seq_len(min(15, nrow(r))), ]
    data.frame(Gene = r$label,
               log2FC = sprintf("%.2f", r$log2FC),
               `adj. p` = formatC(r$p_adj, format = "e", digits = 1),
               Tier = r$tier, check.names = FALSE)
  }, striped = TRUE, spacing = "s")

  output$protein_header <- renderUI({
    p <- prot[prot$protein_group == input$protein, ]
    tagList(
      tags$strong(p$label), " · ", p$protein_group, " · tier: ", p$tier,
      if (p$circular) " · circular (defines the diagnosis)",
      " · ", tags$a("UniProt", href = uniprot_url(p$protein_group), target = "_blank")
    )
  })

  output$forest <- renderPlot({
    f <- d$results[d$results$protein_group == input$protein, ]
    f$cohort <- factor(f$cohort, levels = rev(COHORTS))
    f$label  <- sprintf("%s (n = %d)", f$cohort, f$n_control + f$n_AD)
    f$label  <- factor(f$label, levels = f$label[order(f$cohort)])
    ggplot(f, aes(log2FC, label, colour = cohort == MAIN)) +
      geom_vline(xintercept = 0, colour = "grey70") +
      geom_errorbar(aes(xmin = log2FC - 1.96 * se, xmax = log2FC + 1.96 * se), width = 0.2) +
      geom_point(size = 2.5) +
      scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey45"), guide = "none") +
      labs(title = "Effect by cohort (95% CI)", x = "log2 fold change, AD vs control", y = NULL)
  })

  output$distribution <- renderPlot({
    s <- d$samples
    s$value <- d$intensity[input$protein, s$run]
    ggplot(s[!is.na(s$value), ], aes(site, value, colour = group)) +
      geom_boxplot(outlier.shape = NA, position = position_dodge(0.75), width = 0.6) +
      geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
                 size = 1, alpha = 0.6) +
      scale_colour_manual(values = GRP_COLS, name = NULL) +
      labs(title = "Abundance by site", x = NULL, y = "log2 intensity") +
      theme(legend.position = "top")
  })

  output$tau <- renderPlot({
    s <- d$samples
    s$value <- d$intensity[input$protein, s$run]
    s <- s[!is.na(s$value) & !is.na(s$t_tau), ]
    ggplot(s, aes(t_tau, value, colour = group)) +
      geom_point(size = 1.2, alpha = 0.6) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      scale_x_log10() +
      scale_colour_manual(values = GRP_COLS, name = NULL) +
      labs(title = "Relation to CSF t-tau", x = "t-tau (pg/ml, log scale)", y = "log2 intensity") +
      theme(legend.position = "top")
  })
}

shinyApp(ui, server)
