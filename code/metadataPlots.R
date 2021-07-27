# Title     : SARS-COV-2 Analysis
# Objective : Plots of SARS-COV-2 samples classified to variants
# Created by: Dana Bar-Ilan
# Created on: 5/24/2021

library(ggplot2); library(dplyr); library(reshape); library(viridis); library(hrbrthemes)
library(plotly); library(tidyverse); library(htmlwidgets); library(RColorBrewer)

metadata <- read.csv("C:/Users/User/Documents/NGS/metadata.tsv", sep = '\t', quote = "")
metadata <- metadata[!duplicated(metadata$strain),]
metadata$date <- as.Date(metadata$date)
metadata <- metadata[(metadata$strain != "REF_NC_045512.2") & (metadata$knownVariant != "QC fail"),] # remove refseq row
# gather similare variants together:
metadata <- metadata[!metadata$knownVariant %in% c("suspect", "suspect + 484","suspect + 501"),]
metadata[metadata$knownVariant %in% c("?", "no variant", ""),]$knownVariant <- "no monitored variant"
metadata[metadata$knownVariant %in% c("C.37 - Chile suspect", "C.37 - suspect"),]$knownVariant <- "C.37"
metadata[metadata$knownVariant %in% c("B.1.525 - Global suspect", "B.1.525 - suspect"),]$knownVariant <- "B.1.525"
metadata[metadata$knownVariant %in% c("B.1.1.7 - UK suspect", "B.1.1.7 - UK + 484", "B.1.1.7 - UK",
                                      "B.1.1.7 - UK suspect + 484", "B.1.1.7 - UK + E484K", 
                                      "B.1.1.7 - suspect", "B.1.1.7 - suspect + 484", "B.1.1.7 + 484"),]$knownVariant <- "B.1.1.7"
metadata[metadata$knownVariant == "B.1.351 - suspect",]$knownVariant <- "B.1.351"
metadata[metadata$knownVariant == "P.1 - suspect",]$knownVariant <- "P.1"
metadata[metadata$knownVariant == "B.1.617.2 - suspect",]$knownVariant <- "B.1.617.2"
metadata[metadata$knownVariant %in% c("C.36.3 - suspect", "C.36.3 suspect"),]$knownVariant <- "C.36.3"
metadata[metadata$knownVariant == "VUI_L452R/L1063F_Israel",]$knownVariant <- "VUI_L452R_Israel"

# data frame -> select only date, known vairnat and pangolin clade to plot
df <- metadata %>% select(c("date", "knownVariant", "pangolinClade"))

knownVariant = unique(df$knownVariant) # list of unique variants
vars_dates <- data.frame(date = unique(df$date)) # list of all dates
vars_dates <- crossing(vars_dates, knownVariant) # all variants X all dates

## known variants - Area Plot ##

by_known_var <- df%>% group_by(date, knownVariant) %>% summarise(n=n()) %>% ungroup()  # counts of variants per date 
by_known_var <- vars_dates %>% left_join(by_known_var, by=c("date","knownVariant")) # make sure all dates and variants apear
by_known_var$n[is.na(by_known_var$n)] <- 0  # fill nas with 0 (variant in date count = 0)
by_known_var <- by_known_var %>% subset(date >= "2020-10-01")  # remove older dates
by_known_var <- by_known_var %>% mutate(week=as.Date(cut.Date(date, breaks = "2 weeks", start.on.monday = FALSE))) # divide to weeks
twoWeeks_count <- by_known_var %>% group_by(week, knownVariant) %>% summarise(count=sum(n)) %>% ungroup() # get sum in date for ratio
twoWeeks_count <- twoWeeks_count %>% group_by(week) %>% summarise(count=count, knownVariant=knownVariant,
                                                                  ratio = count/sum(count)) %>% ungroup()
# export plot
pdf(file = "areaPlot_ngs102.pdf", width = 10, height = 7)

area_plot <- ggplot(twoWeeks_count, aes(x=week, y=ratio, fill=knownVariant)) +
  geom_area() + ggtitle("SARS-CoV-2 Israel") + scale_x_date(date_breaks = "2 week") + theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0.5),
                                                                                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_color_brewer(palette = "Set2")
print(area_plot)
dev.off()

## known variants - Lines Plot ##
# instead of ratios, only simple count
twoWeeks_count <- by_known_var %>% group_by(week, knownVariant) %>% summarise(count=sum(n)) %>% ungroup()
twoWeeks_count <- twoWeeks_count %>% group_by(week) %>% summarise(count=count, knownVariant=knownVariant) %>% ungroup()

# export plot
pdf(file = "linesPlot_ngs102.pdf", width = 10, height = 6)

lines_plot <- ggplot(twoWeeks_count, aes(x=week, y=count)) +
  geom_line(aes(color=knownVariant)) + ggtitle("SARS-CoV-2 Israel") + scale_x_date(date_breaks = "2 week") + theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0.5),
                                                                                                                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
print(lines_plot)
dev.off()

## known variants - Stacked Bar Plot ##
# exactly like Lines Plot (count of data)

# export plot
pdf(file = "stackedbarPlot_ngs102.pdf", width = 10, height = 7)
stacekd_bar <- ggplot(twoWeeks_count, aes(fill = knownVariant,x=week, y=count)) +
  geom_bar(position = "stack", stat = "identity") + ggtitle("SARS-CoV-2 Israel") + scale_x_date(date_breaks = "2 week") + theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0.5),
                                                                                                                                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
print(stacekd_bar)
dev.off()

###### code for interactive plotly html ######
p <- plotly_build(area_plot) %>%
  layout(updatemenus = list(
    list(
      type = "buttons",
      direction = "right",
      xanchor = "center",
      yanchor = "top",
      showactive = FALSE,
      x = 1.05,
      y = -0.25,
      buttons = list(
        list(method = "restyle",
             args = list("visible", "all"),
             label = "show all"),
        list(method = "restyle",
             args = list("visible", "legendonly"),
             label = "hide all")
      )
    )
  ))
htmlwidgets::saveWidget(as_widget(p), "VariantsIsrael.html")

