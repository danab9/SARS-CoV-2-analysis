# Title     : SARS-COV-2 Analysis
# Objective : TODO
# Created by: Dana Bar-Ilan
# Created on: 5/24/2021

library(ggplot2); library(dplyr); library(reshape); library(viridis); library(hrbrthemes)
library(plotly); library(tidyverse); library(htmlwidgets); library(RColorBrewer)

nb.cols <- 17
my_colors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

# load metadata tsv
metadata <- read.csv("C:/Users/User/Documents/NGS/metadata.tsv", sep = '\t')
metadata$date <- as.Date(metadata$date)
metadata <- metadata[(metadata$strain != "REF_NC_045512.2") & (metadata$knownVariant != "QC fail"),] # remove refseq row
metadata[metadata$knownVariant %in% c("?", "no variant", ""),]$knownVariant <- "no monitored variant"
metadata[metadata$knownVariant %in% c("C.37", "C.37 - Chile suspect"),]$knownVariant <- "C.37 - Chile"
metadata[metadata$knownVariant == "B.1.429- California",]$knownVariant <- "B.1.429 - California"
metadata[metadata$knownVariant %in% c("B.1.525", "B.1.525 - Global suspect"),]$knownVariant <- "B.1.525 - Global"
metadata[metadata$knownVariant %in% c("B.1.1.7", "B.1.1.7 - UK suspect"),]$knownVariant <- "B.1.1.7 - UK"
metadata[metadata$knownVariant == "B.1.351 - SA suspect",]$knownVariant <- "B.1.351 - SA"
metadata[metadata$knownVariant == "P.1 - Manaus suspect",]$knownVariant <- "P.1 - Manaus"
metadata[metadata$knownVariant == "B.1.1.7 - UK suspect + 484",]$knownVariant <- "B.1.1.7 - UK + 484"

df <- metadata %>% select(c("date", "knownVariant", "pangolinClade"))
tmp <- data.frame(date = unique(df$date))
knownVariant = unique(df$knownVariant)
tmp <- crossing(tmp, knownVariant)


by_known_var <- df%>% group_by(date, knownVariant) %>% summarise(n=n())# only known variant column, not pangoclade
by_known_var <- tmp %>% left_join(by_known_var, by=c("date","knownVariant"))
by_known_var[is.na(by_known_var)] <- 0
by_known_var <- by_known_var %>% subset(date >= "2020-10-01")
by_known_var <- by_known_var %>% mutate(week=as.Date(cut.Date(date, breaks = "2 weeks", start.on.monday = FALSE)))
twoWeeks_count <- by_known_var %>% group_by(week, knownVariant) %>% summarise(count=sum(n)) %>% ungroup()
twoWeeks_count <- twoWeeks_count %>% group_by(week) %>% summarise(count=count, knownVariant=knownVariant,
                                                                  ratio = count/sum(count)) %>% ungroup()
area_plot <- ggplot(twoWeeks_count, aes(x=week, y=ratio, fill=knownVariant)) +
  geom_area() + ggtitle("SARS-CoV-2 Israel") + scale_x_date(date_breaks = "2 week") + theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0.5),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_color_brewer(palette = "Se")
print(area_plot)


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

# without british variant
# by_known_var_nouk <- df %>% group_by(date, knownVariant) %>% summarise(n=n())
# by_known_var_nouk <- tmp %>% left_join(by_known_var_nouk, by=c("date","knownVariant"))
# by_known_var_nouk[is.na(by_known_var_nouk)] <- 0
# by_known_var_nouk <- by_known_var_nouk[by_known_var_nouk$knownVariant != "B.1.1.7 - UK",]
# by_known_var_nouk <- by_known_var_nouk %>% group_by(date) %>% summarise(ratio = n/sum(n), knownVariant=knownVariant)
# by_known_var_nouk <- by_known_var_nouk %>% subset(date >= "2020-11-01")
#













#
# area_plot_nouk <- ggplot(by_known_var_nouk, aes(x=date, y=ratio, fill=knownVariant)) +
#   geom_area(size=0.2, colour="black") + ggtitle("SARS-CoV-2 Israel") + scale_x_date(date_breaks = "2 week")
# p_nouk <- plotly_build(area_plot_nouk) %>%
# layout(updatemenus = list(
#     list(
#         type = "buttons",
#         direction = "right",
#         xanchor = "center",
#         yanchor = "top",
#         showactive = FALSE,
#         x = 1.05,
#         y = -0.25,
#         buttons = list(
#             list(method = "restyle",
#                  args = list("visible", "all"),
#                  label = "show all"),
#             list(method = "restyle",
#                  args = list("visible", "legendonly"),
#                  label = "hide all")
#         )
#     )
# ))
# htmlwidgets::saveWidget(as_widget(p_nouk), "VariantsIsrael_noUK.html")