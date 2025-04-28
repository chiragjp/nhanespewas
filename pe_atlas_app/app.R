# PE Browser App
# Developed using Posit Shiny
# Chirag Patel
# chirag@hms.harvard
# 12/30/24



library(tidyverse)
library(DBI)
library(ggsci)
library(DT)
library(ggrepel)
library(cowplot)
library(shiny)
library(corrr)
library(gplots)
library(pool)
library(reactable)
library(igraph)
library(ggridges)
library(ggpubr)

## get data
pool <- dbPool(drv = RSQLite::SQLite(), dbname='pe_shiny_2.sqlite')
##

p_variables <- tbl(pool, 'adjusted_meta_2') |> group_by(pvarname) |> count() |> left_join(tbl(pool, "p_variable_domain"), by=c("pvarname"="Variable.Name")) |> collect()
p_variables <- p_variables |> mutate(cat_subcat=ifelse(!is.na(psubcategory), paste(pcategory, psubcategory, sep="-"), pcategory ))
p_variables <- p_variables |> mutate(pvardesc_selector = sprintf("%s-(%s)", pvardesc, pvarname))

e_variables <- tbl(pool, 'adjusted_meta_2') |> group_by(evarname) |> count() |> left_join(tbl(pool, "e_variable_domain"), by=c("evarname"="Variable.Name")) |> collect()
e_variables <- e_variables |> mutate(cat_subcat=ifelse(!is.na(esubcategory), paste(ecategory, esubcategory, sep="-"), ecategory ))
e_variables <- e_variables |> mutate(evardesc_selector = sprintf("%s-(%s)", evardesc, evarname))
e_category_strs <- e_variables |> ungroup() |> select(ecategory, esubcategory) |> group_by(ecategory, esubcategory) |> count()
e_category_strs <- e_category_strs |> ungroup() |> mutate(cat_subcat=ifelse(!is.na(esubcategory), paste(ecategory, esubcategory, sep="-"), ecategory ))



## ui
ui <- navbarPage("Phenome-Exposome Atlas",
    tabPanel("About",
      mainPanel(h1("Welcome to the Phenome-Exposome Atlas!"),
                p("How much variation do exposures explain in phenotype?"),
                p("This is an atlas of 127K correlations between 278 phenotypes (e.g., body mass index, glucose, height, creatinine) and 651 exposures and behaviors (e.g., self-reported nutrient intake; blood lead levels; urinary phthalates; blood PFOA) across the National Health and Nutrition Examination Surveys (NHANES) in individuals from 1999-2017. For each survey, we associate the phenotype and exposure, and summarize the association size across the surveys. We assess the robustness of associations by estimating the false discovery rate, consistency with adjustments, and concordance across multiple waves of the surveys."),
                img(src="pe_fig1.png",height="650px"),
                h3("GitHub Repository:"),
                p(a("nhanes_pewas", href="https://github.com/chiragjp/nhanes_pewas"))
                )
    ),
    tabPanel("Browse By Phenotype",
        sidebarPanel(
            selectInput("pvarname", h3("By Phenotype"), choices=setNames(p_variables$pvarname,p_variables$pvardesc_selector), selected="LBXGLU")
        ),
        mainPanel(
            h1(textOutput("pvariable_description")),
            h2("Category:", textOutput("pvariable_subcategory")),
            tabsetPanel(type = "tabs",
                tabPanel("-log(P) vs. Exposure Group", plotOutput("by_phenotype_manhattan_plot")),
                tabPanel("Effect Sizes", plotOutput("by_phenotype_effect_size")),
                tabPanel("R^2", plotOutput("by_phenotype_r2_plot")),
                tabPanel("-log(P) vs. Effect Size", plotOutput("by_phenotype_volcano_plot"),
                         p("Points in the background depict summary stats for the entire category")),
                tabPanel("Globe", plotOutput("by_phenotype_globe_plot"),
                        p("Nodes are exposures associated with phenotype"),
                        p("Grey lines depict positive and significant correlations"),
                        p("Black lines depict negative and significant correlations"),
                        h3("Globe Abbreviation Legend"),
                         DTOutput("by_phenotype_globe_legend_table")),
                tabPanel("Specification Bias", plotOutput("by_phenotype_estimate_plot"))
            ),
            h3("Summary Statistics"),
            reactableOutput("by_phenotype_table"),
            ## download button
            downloadButton("by_phenotype_table_download", "Download Summary Statistics"),
            ## put in correlations here
            h3("Similar Phenotypes"),
            ## put in multivariate model
            reactableOutput("exposure_corr_table")
        )
    ),

    tabPanel("Browse By Exposure",
             sidebarPanel(
                 selectInput("evarname", h3("By Exposure"), choices=setNames(e_variables$evarname,e_variables$evardesc_selector), selected="LBXBCD")
             ),
             mainPanel(
               h1(textOutput("evariable_description")),
               h2("Category:", textOutput("evariable_subcategory")),
               tabsetPanel(type = "tabs",
                           tabPanel("-log(P) vs. Effect Size", plotOutput("by_exposure_volcano_plot"),p("Points in the background depict summary stats for the entire category")),
                           tabPanel("ES Histogram", plotOutput("by_exposure_context_plot"), p("Displaying effect sizes for FDR less than 5%")),
                           tabPanel("R^2", plotOutput("by_exposure_r2_plot")),
                           tabPanel("Specification Bias", plotOutput("by_exposure_estimate_plot"))
               ),
               h3("Summary Statistics"),
               reactableOutput("by_exposure_table"),
               ## download button
               downloadButton("by_exposure_table_download", "Download Summary Statistics"),
               ## put in correlations here
               h3("Similar Exposures"),
               ## put in multivariate model
               reactableOutput("phenotype_corr_table")
             )
    ),

    tabPanel("Browse by Exposure Group",
               sidebarPanel(
               selectInput("egroup", h3("By Group"), choices=setNames(e_category_strs$cat_subcat,e_category_strs$cat_subcat), selected="pollutant-hydrocarbon")
             ),
             mainPanel(
               h1(textOutput("egroup")),
               tabsetPanel(type="tabs",
                           tabPanel("-log(P) vs. Effect Size", plotOutput("by_exposure_group_volcano_plot")),
                           tabPanel("ES Histogram", plotOutput("by_exposure_group_context_plot")),
                           tabPanel("R^2", plotOutput("by_exposure_group_r2_plot")),
                           tabPanel("Specification Bias", plotOutput("by_exposure_group_estimate_plot"))
                ),
               h3("Summary Statistics for each E"),
               reactableOutput("by_exposure_group_table"),
               downloadButton("by_exposure_group_table_download", "Download Summary Statistics")
             )
    ),
  tabPanel("E Catalog",
    mainPanel(
      h1("Exposome Catalog"),
      reactableOutput("e_catalog_table")
    )
  ),
  tabPanel("P Catalog",
           mainPanel(
             h1("Phenome Catalog"),
             reactableOutput("p_catalog_table")
           )
  )

)

### table and figure functions

## for exposure or phenotypes
plot_adjusted_r2_cdf <- function(to_plot) {
  p <- ggplot(to_plot, aes(r2_adjusted_vs_base, color=sig_levels))
  p <- p + stat_ecdf() + scale_x_continuous(limits=c(0, .05)) + scale_color_tron()
  p <- p + theme_bw() + xlab("R-squared")
  p <- p + theme(legend.position = "bottom")
  p
}


plot_volcano <- function(to_plot, bg=NULL) {
  p <- ggplot(to_plot, aes(estimate, -log10(p.value), color=sig_levels))
  p <- p + geom_point(shape = 1, size = 1) # color = "black")
  if(!is.null(bg)) {
    p <- p + geom_point(data=bg, aes(estimate, -log10(p.value)), alpha=.01)
  }
  p <- p + scale_x_continuous(limits=c(-.2, .2)) + scale_color_tron()
  p <- p + theme_bw() + xlab("Estimate (per 1 SD of scale(E) on scale(P))") + ylab("-log10(overall p.value)")
  p
}

plot_estimates <- function(to_plot) {
  p <- ggplot(to_plot, aes(x=estimate_unadjusted, y=estimate))
  p <- p + geom_point(shape = 1, size = 1, color = "red") + scale_x_continuous(limits=c(-.5, .5)) + scale_y_continuous(limits=c(-.5, .5))
  p <- p + geom_abline() + geom_smooth(se = FALSE, method = lm)
  p <- p + xlab("Estimate (Univariate Model)") + ylab("Estimate (Demographics Adjusted)")
  p <- p + theme_bw()
  p
}


reactable_tibble <- function(meta_obj) {
  to_react <- meta_obj |>
    select(pvardesc, evardesc,estimate,expoq2, expoq3, expoq4, expoq5,rsq_adjusted_base_diff, rsq_adjusted_diff, p.value, sig_levels, total_n, i.squared.uwls, nobs) |> collect() |>
    mutate(neglogpvalue=-log10(p.value)) |> select(-p.value) |>
    rename(effect.size=estimate,r.squared=rsq_adjusted_base_diff,n_surveys=nobs, total_sample_size=total_n, p_description=pvardesc, e_description=evardesc) |>
    rename(es_25_50=expoq2, es_50_75=expoq3, es_75_90=expoq4, es_90_100=expoq5)
  to_react <- to_react |> select(p_description, e_description, effect.size, es_25_50, es_50_75, es_75_90, es_90_100, neglogpvalue, sig_levels, r.squared, i.squared.uwls, total_sample_size, n_surveys)
}

reactable_summary_stats <- function(meta_obj) {
  to_react <- reactable_tibble(meta_obj)
  reactable(to_react |> arrange(-neglogpvalue), columns = list(effect.size=colDef(format = colFormat(digits=3)),
                                                               es_25_50=colDef(format = colFormat(digits=3)),
                                                               es_50_75=colDef(format = colFormat(digits=3)),
                                                               es_75_90=colDef(format = colFormat(digits=3)),
                                                               es_90_100=colDef(format = colFormat(digits=3)),
                                                               neglogpvalue=colDef(format=colFormat(digits=1)),
                                                               i.squared.uwls=colDef(format = colFormat(digits=2)),
                                                               r.squared=colDef(format=colFormat(digits=4))

  ),searchable = TRUE, filterable=FALSE) # end reactable
}



plot_manhattan <- function(to_plot) {
  PCAP <- 20
  ## order in terms of number found; or R2.
  enewsubcategory_order <- to_plot |> filter(sig_levels == 'Bonf.<0.05') |> group_by(enewsubcategory) |> summarize(med_r2=mean(rsq_adjusted_base_diff)) |> arrange(-med_r2) |> mutate(group_number = row_number())

  to_plot <- to_plot |> left_join(enewsubcategory_order |> select(enewsubcategory, group_number))
  eplot_order <- to_plot |> select(evarname ,enewsubcategory, group_number) |> unique() |> arrange(group_number) |> mutate(order_number = row_number())

  xaxis_labels <- eplot_order |> group_by(group_number,enewsubcategory) |> summarize(med_x=round(median(order_number)), min_x=min(order_number), max_x=max(order_number))
  to_plot <- to_plot |> left_join(eplot_order |> select(evarname,order_number))
  p <- ggplot(to_plot |> arrange(order_number), aes(order_number, -log10(p.value)))
  p <- p + geom_point(aes(color=sig_levels), shape=20) + scale_color_tron()
  p <- p + scale_y_continuous(limits=c(0, PCAP+1))
  p <- p + scale_x_continuous(breaks = xaxis_labels$med_x,  labels = xaxis_labels$enewsubcategory)
  p <- p + geom_vline(data=xaxis_labels, aes(xintercept=max_x), linetype='dotted', col = 'black') + theme_bw()
  #p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom")
  p <- p + theme(axis.text.x = element_blank(), legend.position = "bottom")
  p1 <- p + xlab("") + ylab("-log10(pvalue)")

  p <- ggplot(to_plot |> arrange(order_number), aes(order_number, rsq_adjusted_base_diff))
  p <- p + geom_point(aes(color=sig_levels), shape=20) + scale_color_tron()
  p <- p + scale_x_continuous(breaks = xaxis_labels$med_x,  labels = xaxis_labels$enewsubcategory)
  p <- p + geom_vline(data=xaxis_labels, aes(xintercept=max_x), linetype='dotted', col = 'black') + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom")
  p2 <- p + xlab("") + ylab("r.sq")
  ggarrange(p1, p2,  nrow = 2, common.legend = TRUE)
}

plot_effect_size <- function(to_plot) {
  # distribution
  pvarname_lcl <- unique(to_plot$pvarname)[1]
  p <- ggplot(to_plot |> filter(!is.na(sig_levels)) |> mutate(sig_levels = factor(sig_levels, levels=c("> BY & Bonf.", "BY<0.05", "Bonf.<0.05"))),
              aes(y=exposure_level, x=I(abs(estimate)), fill=sig_levels))
  p <- p + geom_violin()
  p <- p + xlab("Absolute Value Effect Size [1 Unit of Phenotype (original units)]") + ylab("Quantile of Exposure") + ylab("Percentile of Exposure")
  p <- p + scale_fill_tron()
  p <- p + theme_bw() + theme(legend.position = "bottom")
  p
}

exposome_globe_plot <- function(exposure_tibble) {
  E_CORR_THRESHOLD_NEG <- -0.2
  E_CORR_THRESHOLD_POS <- 0.2
  MAX_NUMBER_NODES <- 100
  corr_tibble <- tbl(pool, "ee")
  ff <- corr_tibble |> right_join(exposure_tibble , by=c("xvarname"="evarname"))
  sig <- ff |> right_join(exposure_tibble, by=c("yvarname"="evarname"))

  sig_to_plot <- sig |> filter(estimate > E_CORR_THRESHOLD_POS | estimate < E_CORR_THRESHOLD_NEG) ## filter among the highest correlations -- choose some background distribution
  sig_to_plot <- sig_to_plot |> collect() |> arrange(desc(abs(estimate))) |> filter(row_number() <= MAX_NUMBER_NODES)

  sig_graph <- sig_to_plot |> select(xvarname, yvarname, everything())
  sig_vertex <- sig_graph |> select(xvarname, xcategory) |> unique() |> rename(varname = xvarname, category=xcategory)
  sig_vertex2 <- sig_graph |> select(yvarname, ycategory) |> unique() |> rename(varname = yvarname, category=ycategory)
  sig_vertex <- sig_vertex |> rbind(sig_vertex2) |> unique()
  sig_graph <- sig_graph |> mutate(estimate = ifelse(estimate > 1, 1, estimate)) |> mutate(estimate = ifelse(estimate < -1, -1, estimate))

  g <- graph_from_data_frame(sig_graph, directed = FALSE)
  E(g)$color <- ifelse(E(g)$estimate > 0, "black", "grey")
  V(g)$category <- sig_vertex$category[match(V(g)$name, sig_vertex$varname)]
  categories <- unique(V(g)$category)
  # Assign colors to each category
  color_pal <- ggsci::pal_tron("legacy", alpha = 0.7)
  category_colors <- setNames(color_pal(length(categories)), categories)
  # Map the colors to vertices based on category
  vertex_colors <- category_colors[V(g)$category]
  coords <- layout_in_circle(g)
  scaling_factor <- 5
  par(mar = c(0, 0, 0, 0))
  plot(g, vertex.color = vertex_colors, layout=coords,
       vertex.size = 10,
       vertex.color = vertex_colors,
       vertex.label = substr(V(g)$name, nchar(V(g)$name)-2, nchar(V(g)$name)),      # Label each vertex with its name
       vertex.label.cex = 0.7,         # Control the size of the labels
       vertex.label.color = "black",   # Color of the labels
       vertex.label.family = "sans",
       edge.width = abs(sqrt(E(g)$estimate)) * scaling_factor
  )
}


# server logic
server <- function(input, output, session) {
    output$by_phenotype_volcano_plot <- renderPlot({
        pvarn <- input$pvarname
        pcat <- p_variables |> filter(pvarname == pvarn) |> pull(pcategory)
        to_plot <- tbl(pool, "adjusted_meta_2") |> filter(pvarname == pvarn) |> select(estimate, p.value, sig_levels) |> collect()
        to_plot <- to_plot |> mutate(p.value = ifelse(p.value < 1e-20, 1e-20, p.value))
        bg <- tbl(pool, "adjusted_meta_2") |> filter(pcategory == pcat) |> select(estimate, p.value, sig_levels) |> collect()
        bg <- bg |> mutate(p.value = ifelse(p.value < 1e-20, 1e-20, p.value))
        plot_volcano(to_plot, bg)
    })

    output$by_phenotype_manhattan_plot <- renderPlot({
      pvarn <-input$pvarname
      pcat <- p_variables |> filter(pvarname == pvarn) |> pull(pcategory)
      to_plot <- tbl(pool, "adjusted_meta_2") |> filter(pvarname == pvarn) |> collect()
      to_plot <- to_plot |> mutate(p.value = ifelse(p.value < 1e-20, 1e-20, p.value))
      plot_manhattan(to_plot)
    })


    output$by_phenotype_globe_plot <- renderPlot({
      pvarn <-input$pvarname
      pcat <- p_variables |> filter(pvarname == pvarn) |> pull(pcategory)
      exposures_to_plot <- tbl(pool, "adjusted_meta_2") |> filter(pvarname == pvarn, sig_levels == "Bonf.<0.05")|> select(evarname)
      exposome_globe_plot(exposure_tibble = exposures_to_plot)
    })
    output$by_phenotype_globe_legend_table <- renderDT({
      pvarn <- input$pvarname
      exposures_to_plot <- tbl(pool, "adjusted_meta_2") |> filter(pvarname == pvarn, sig_levels == "Bonf.<0.05")|> select(evarname, evardesc)
      to_table <- exposures_to_plot |> collect()
      to_table <- to_table |> mutate(evarname = substr(evarname, nchar(evarname)-2, nchar(evarname))) |> rename(abbr=evarname, exposure_description=evardesc)
    })
    output$by_phenotype_effect_size <-renderPlot({
      pvarn <-input$pvarname
      pcat <- p_variables |> filter(pvarname == pvarn) |> pull(pcategory)
      es <- tbl(pool, "pe_quantile_ns") |> filter(pvarname == pvarn) |> collect()
      plot_effect_size(es)
    })

    output$by_phenotype_r2_plot <- renderPlot({
        pvarn <- input$pvarname
        to_plot <- tbl(pool, "adjusted_meta_2") |> filter(pvarname==pvarn) |> select(rsq_adjusted_base_diff,rsq_adjusted_diff, sig_levels) |> collect()
        to_plot <- to_plot |> rename(r2_adjusted_vs_base=rsq_adjusted_base_diff, r2_adjusted_vs_unadjusted=rsq_adjusted_diff)
        plot_adjusted_r2_cdf(to_plot)
    })

    output$by_phenotype_estimate_plot <- renderPlot({
        pvarn <- input$pvarname
        to_plot <- tbl(pool, "expos_wide") |> filter(pvarname == pvarn) |> select(estimate_2,estimate_1) |> rename(estimate=estimate_2, estimate_unadjusted=estimate_1) |> collect()
        plot_estimates(to_plot)
    })

    output$by_phenotype_table <- renderReactable({
        pvarn <- input$pvarname
        to_react <- tbl(pool, "adjusted_meta_2") |> filter(pvarname == pvarn)
        qns <- tbl(pool, "pe_quantile_ns") |> filter(pvarname == pvarn) |>
          select(estimate, evarname, pvarname, term) |> pivot_wider(names_from =c("term"), values_from=c("estimate"))
        to_react <- to_react |> left_join(qns, by=c("evarname", "pvarname"))
        reactable_summary_stats(to_react)
    })

    output$by_phenotype_table_download <- downloadHandler(
      ## to implement
      filename = function() {
        pvarn <- input$pvarname
        paste("by_pheno-", pvarn, '-', Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        pvarn <- input$pvarname
        to_react <- tbl(pool, "adjusted_meta_2") |> filter(pvarname == pvarn) |> reactable_tibble()
        write_csv(x=to_react, file)
      }
    )

    output$phenotype_corr_table <- renderReactable({
        pvarn <- input$pvarname
        to_react <- exposure_correlation |> filter(x == pvarn) |> select(-x) |> rename(corr=r)
        to_react <- to_react |> left_join(p_variable_domain, by=c("y"="Variable.Name")) |> select("pvardesc", "pcategory", "corr")
        reactable(to_react |> arrange(-corr), columns=list(corr=colDef(format=colFormat(digits=3))))
    })

    output$pvariable_description <- renderText({
        pvarn <- input$pvarname
        ret <- tbl(pool, "p_variable_domain") |> filter(Variable.Name == pvarn) |> pull(pvardesc)
        ret[1]
    })

    output$pvariable_category <- renderText({
      pvarn <- input$pvarname
      ret <- tbl(pool, "p_variable_domain") |> filter(Variable.Name == pvarn) |> pull(pcategory)
      ret[1]
    })

    output$evariable_description <- renderText({
        evarn <- input$evarname
        ret <- tbl(pool, "e_variable_domain") |> filter(Variable.Name == evarn) |> pull(evardesc)
        ret[1]
    })

    output$evariable_category <- renderText({
      evarn <- input$evarname
      ret <- tbl(pool, "e_variable_domain") |> filter(Variable.Name == evarn) |> pull(ecategory)
      ret[1]
    })

    output$pvariable_subcategory <- renderText({
      pvarn <- input$pvarname
      ret <- p_variables |> filter(pvarname == pvarn) |> pull(cat_subcat)
      ret[1]
    })

    output$evariable_subcategory <- renderText({
      evarn <- input$evarname
      ret <- e_variables |> filter(evarname == evarn) |> pull(cat_subcat)
      ret[1]
    })

    output$by_exposure_volcano_plot <- renderPlot({
        evarn <- input$evarname
        ecat <- e_variables |> filter(evarname == evarn) |> pull(ecategory)
        to_plot <- tbl(pool, "adjusted_meta_2") |> filter(evarname == evarn) |> select(estimate, p.value, sig_levels) |> collect()
        to_plot <- to_plot |> mutate(p.value = ifelse(p.value < 1e-20, 1e-20, p.value))
        bg <- tbl(pool, "adjusted_meta_2") |> filter(ecategory == ecat) |> select(estimate, p.value, sig_levels) |> collect()
        bg <- bg |> mutate(p.value = ifelse(p.value < 1e-20, 1e-20, p.value))
        plot_volcano(to_plot, bg)
    })

    output$by_exposure_r2_plot <- renderPlot({
      evarn <- input$evarname
      to_plot <- tbl(pool, "adjusted_meta_2") |> filter(evarname==evarn) |> select(rsq_adjusted_base_diff,rsq_adjusted_diff, sig_levels) |> collect()
      to_plot <- to_plot |> rename(r2_adjusted_vs_base=rsq_adjusted_base_diff, r2_adjusted_vs_unadjusted=rsq_adjusted_diff)
      plot_adjusted_r2_cdf(to_plot)
    })


    output$by_exposure_table <- renderReactable({
      evarn <- input$evarname
      to_react <- tbl(pool, "adjusted_meta_2") |> filter(evarname == evarn)
      qns <- tbl(pool, "pe_quantile_ns") |> filter(evarname == evarn) |>
        select(estimate, evarname, pvarname, term) |> pivot_wider(names_from =c("term"), values_from=c("estimate"))
      to_react <- to_react |> left_join(qns, by=c("evarname", "pvarname"))
      reactable_summary_stats(to_react)
    })


    output$by_exposure_table_download <- downloadHandler(
      ## to implement
      filename = function() {
        evarn <- input$evarname
        paste("by_expo-", evarn, '-', Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        evarn <- input$evarname
        to_react <- tbl(pool, "adjusted_meta_2") |> filter(evarname == evarn) |> reactable_tibble()
        write_csv(x=to_react, file)
      }
    )

    output$by_exposure_estimate_plot <- renderPlot({
      evarn <- input$evarname
      to_plot <- tbl(pool, "expos_wide") |> filter(evarname == evarn) |> select(estimate_2,estimate_1) |> rename(estimate=estimate_2, estimate_unadjusted=estimate_1) |> collect()
      plot_estimates(to_plot)
    })

    output$by_exposure_context_plot <- renderPlot({
      evarn <- input$evarname
      m <- tbl(pool, "adjusted_meta_2") |> filter(pval_BY < 0.05) |> filter(ecategory != 'allergy')
      e_estimate <- m |> filter(evarname == evarn) |> select(estimate,sig_levels) |> collect()
      vardesc <- e_variables |> filter(evarname == evarn)  |> pull(evardesc)
      vardesc <- vardesc[1]
      other_estimates <- m |> filter(evarname != evarn) |> select(ecategory, esubcategory,estimate, sig_levels) |> collect()
      other_estimates <- other_estimates |> mutate(cat_subcat=ifelse(!is.na(esubcategory), paste(ecategory, esubcategory, sep="-"), ecategory ))
      other_estimates <- other_estimates |> select(-c(ecategory, esubcategory))
      to_plot <- e_estimate |> mutate(cat_subcat=vardesc) |> rbind(other_estimates)
      to_plot <- to_plot |> mutate(col=ifelse(cat_subcat == vardesc, 'red', 'black'))
      to_plot <- to_plot |> mutate(cat_subcat=fct_relevel(cat_subcat, vardesc, after = 0 ))
      p <- ggplot(to_plot, aes(y=cat_subcat, x=estimate, fill=col))
      p <- p + geom_density_ridges(scale = 3) + scale_x_continuous(limits=c(-.5, .5))
      p <- p + theme_bw() + ylab("") + xlab("Adjusted Estimate") + scale_fill_tron()
      p <- p + theme(legend.position = "none")
      p
    })

    output$exposure_corr_table <- renderReactable({
      pvarn <- input$pvarname
      to_react <- tbl(pool, "exposure_correlation") |> filter(x == pvarn) |> select(-x) |> collect() |> rename(corr=r)
      to_react <- to_react |> left_join(p_variables, by=c("y"="pvarname")) |> select("pvardesc", "pcategory", "corr")
      reactable(to_react |> arrange(-corr), columns=list(corr=colDef(format=colFormat(digits=3))))
    })

    output$phenotype_corr_table <- renderReactable({
      evarn <- input$evarname
      to_react <- tbl(pool, "phenotype_correlation") |> filter(x == evarn) |> select(-x) |> collect() |> rename(corr=r)
      to_react <- to_react |> left_join(e_variables, by=c("y"="evarname")) |> select("evardesc", "ecategory", "corr")
      reactable(to_react |> arrange(-corr), columns=list(corr=colDef(format=colFormat(digits=3))))
    })

    output$egroup <- renderText({
      input$egroup
    })

    output$by_exposure_group_volcano_plot <- renderPlot({
      cat_sub <- input$egroup
      category_and_sub <- e_category_strs |> filter(cat_subcat == cat_sub)
      subcategory_str <- category_and_sub |> pull(esubcategory)
      ecategory_str <- category_and_sub |> pull(ecategory)
      to_plot <- tbl(pool, "adjusted_meta_2") |> filter(ecategory == ecategory_str)
      if(!is.na(subcategory_str)) {
        to_plot <- to_plot |> filter(esubcategory == subcategory_str)
      }
      to_plot <- to_plot |> select(estimate, p.value, sig_levels) |> collect()
      to_plot <- to_plot |> mutate(p.value = ifelse(p.value < 1e-20, 1e-20, p.value))
      plot_volcano(to_plot)
    })


    output$by_exposure_group_table <- renderReactable({
      cat_sub <- input$egroup
      category_and_sub <- e_category_strs |> filter(cat_subcat == cat_sub)
      subcategory_str <- category_and_sub |> pull(esubcategory)
      ecategory_str <- category_and_sub |> pull(ecategory)
      to_react <- tbl(pool, "adjusted_meta_2") |> filter(ecategory == ecategory_str)
      qns <- tbl(pool, "pe_quantile_ns") |> filter(ecategory == ecategory_str)
      if(!is.na(subcategory_str)) {
        to_react <- to_react |> filter(esubcategory == subcategory_str)
        qns <- tbl(pool, "pe_quantile_ns") |> filter(esubcategory == subcategory_str)
      }
      qns <- qns |> select(estimate, evarname, pvarname, term) |> pivot_wider(names_from =c("term"), values_from=c("estimate"))
      to_react <- to_react |> left_join(qns, by=c("evarname", "pvarname"))
      reactable_summary_stats(to_react)
    })

    output$by_exposure_group_table_download <- downloadHandler(
      filename = function() {
        cat_sub <- input$egroup
        paste("by_expo_group-", cat_sub, '-', Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        cat_sub <- input$egroup
        category_and_sub <- e_category_strs |> filter(cat_subcat == cat_sub)
        subcategory_str <- category_and_sub |> pull(esubcategory)
        ecategory_str <- category_and_sub |> pull(ecategory)
        to_react <- tbl(pool, "adjusted_meta_2") |> filter(ecategory == ecategory_str)
        if(!is.na(subcategory_str)) {
          to_react <- to_react |> filter(esubcategory == subcategory_str)
        }

        write_csv(x=to_react |> reactable_tibble(), file)
      }
    )

    output$by_exposure_group_r2_plot <- renderPlot({
      cat_sub <- input$egroup
      category_and_sub <- e_category_strs |> filter(cat_subcat == cat_sub)
      subcategory_str <- category_and_sub |> pull(esubcategory)
      ecategory_str <- category_and_sub |> pull(ecategory)
      to_plot <- tbl(pool, "adjusted_meta_2") |> filter(ecategory == ecategory_str)
      if(!is.na(subcategory_str)) {
        to_plot <- to_plot |> filter(esubcategory == subcategory_str)
      }
      to_plot <- to_plot |> select(rsq_adjusted_base_diff, rsq_adjusted_diff, sig_levels) |> collect()
      to_plot <- to_plot |> rename(r2_adjusted_vs_base=rsq_adjusted_base_diff, r2_adjusted_vs_unadjusted=rsq_adjusted_diff)
      plot_adjusted_r2_cdf(to_plot)
    })

    output$by_exposure_group_estimate_plot <- renderPlot({
      cat_sub <- input$egroup
      category_and_sub <- e_category_strs |> filter(cat_subcat == cat_sub)
      subcategory_str <- category_and_sub |> pull(esubcategory)
      ecategory_str <- category_and_sub |> pull(ecategory)
      to_plot <- tbl(pool, "expos_wide") |> filter(ecategory == ecategory_str)
      if(!is.na(subcategory_str)) {
        to_plot <- to_plot |> filter(esubcategory == subcategory_str)
      }
      to_plot <- to_plot |> select(estimate_2,estimate_1) |> rename(estimate_unadjusted=estimate_1, estimate=estimate_2) |> collect()
      plot_estimates(to_plot)
    })

    output$by_exposure_group_context_plot <- renderPlot({
      cat_sub <- input$egroup
      to_plot <- tbl(pool, "adjusted_meta_2") |> filter(pval_BY < 0.05) |> collect() |>
        mutate(cat_subcat=ifelse(!is.na(esubcategory), paste(ecategory, esubcategory, sep="-"), ecategory ))

      to_plot <- to_plot |> mutate(col=ifelse(cat_subcat == cat_sub, 'red', 'black'))
      to_plot <- to_plot |> mutate(cat_subcat=fct_relevel(cat_subcat, cat_sub, after = 0 ))
      p <- ggplot(to_plot, aes(y=cat_subcat, x=estimate, fill=col))
      p <- p + geom_density_ridges(scale = 3) + scale_x_continuous(limits=c(-.5, .5))
      p <- p + theme_bw() + ylab("") + xlab("Adjusted Estimate") + scale_fill_tron()
      p <- p + theme(legend.position = "none")
      p

    })

  output$e_catalog_table <- renderReactable({
    am2 <- tbl(pool, "adjusted_meta_2") |> select(evarname, enewsubcategory, evardesc) |> collect() |> unique()
    reactable(am2 |> rename(ID=evarname, Category=enewsubcategory, Description=evardesc), filterable = TRUE, searchable = TRUE, pageSizeOptions = c(50, 100, 200), showPageSizeOptions = TRUE)
  })

  output$p_catalog_table <- renderReactable({
    am2 <- tbl(pool, "adjusted_meta_2") |> select(pvarname, pnewsubcategory, pvardesc) |> collect() |> unique()
    reactable(am2 |> rename(ID=pvarname, Category=pnewsubcategory, Description=pvardesc), filterable = TRUE, searchable = TRUE, pageSizeOptions = c(50, 100, 200), showPageSizeOptions = TRUE)
  })

}

shinyApp(ui, server)



