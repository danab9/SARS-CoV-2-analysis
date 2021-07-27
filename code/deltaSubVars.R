# Title     : Delta Sub Vairants - analysis
# Objective : Plots of SARS-COV-2 Delta sub variats
# Created by: Dana Bar-Ilan
# Created on: July 8th 2021
library(dplyr); library(ggplot2); library(tidyverse)

setwd("~/sars-cov-2-analysis")
table = read.table("data/delta_subvars.txt", sep = '\t', header = TRUE)
table$date <- as.Date(table$date, format="%d/%m/%y")
table <- na.omit(table)

delta.subtype = unique(table$delta.subtype)
vars_dates <- data.frame(date = unique(table$date)) %>% crossing(delta.subtype)


by_variant <- table%>% group_by(date, delta.subtype) %>% summarise(n=n()) %>% ungroup()  # counts of variants per date 
by_variant <- vars_dates %>% left_join(by_variant, by=c("date","delta.subtype"))
by_variant$n[is.na(by_variant$n)] <- 0
by_variant <- by_variant %>% mutate(week=as.Date(cut.Date(date, breaks = "2 weeks", start.on.monday = FALSE)))

twoWeeks_count <- by_variant %>% group_by(week, delta.subtype) %>% summarise(count=sum(n)) %>% ungroup()
twoWeeks_count <- twoWeeks_count %>% group_by(week) %>% summarise(count=count, delta.subtype=delta.subtype) %>% ungroup()

# stacked bar:
pdf(file = "stackedBar_DeltaSubtypes.pdf", width = 10, height = 7)
stacekd_bar <- ggplot(twoWeeks_count, aes(fill = delta.subtype,x=week, y=count)) +
  geom_bar(position = "stack", stat = "identity") + ggtitle("Delta sub-types") + scale_x_date(date_breaks = "2 week") + theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0.5),
                                                                                                                                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
print(stacekd_bar)
dev.off()

# area:


# known variants - Area Plot
by_variant <- table%>% group_by(date, delta.subtype) %>% summarise(n=n()) %>% ungroup()  # counts of variants per date 
by_variant <- vars_dates %>% left_join(by_variant, by=c("date","delta.subtype"))
by_variant$n[is.na(by_variant$n)] <- 0

by_variant <- by_variant %>% mutate(week=as.Date(cut.Date(date, breaks = "2 weeks", start.on.monday = FALSE)))
twoWeeks_count <- by_variant %>% group_by(week, delta.subtype) %>% summarise(count=sum(n)) %>% ungroup()
twoWeeks_count <- twoWeeks_count %>% group_by(week) %>% summarise(count=count, delta.subtype=delta.subtype,
                                                                  ratio = count/sum(count)) %>% ungroup()

pdf(file = "areaPlot_DeltaSubtype.pdf", width = 10, height = 7)

area_plot <- ggplot(twoWeeks_count, aes(x=week, y=ratio, fill=delta.subtype)) +
  geom_area() + ggtitle("Delta subtypes Israel") + scale_x_date(date_breaks = "2 week") + theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 0.5),
                                                                                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_color_brewer(palette = "Set2")
print(area_plot)
dev.off()

