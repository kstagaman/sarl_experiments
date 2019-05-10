<- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- <- # chorionate vs dechorionated example

library(data.table)
library(ggplot2)
library(cowplot)

setwd("~/Google Drive/OSU/OSU_Home_Sync/Zf_BaP/GF-vs-CVZ-wChorions/")

movement <- as.data.table(read.csv("behavior_24h_FOR_GUI_HMAT_With_Removal_2019-03-07_Single_Tab_File  03-07-2019 145240.csv"))
View(movement)
idVars <- names(movement)[!(grepl("^t", names(movement)))]
movement.melt <- melt(movement, id.vars = idVars)
movement.melt[, variable := as.numeric(gsub("t", "", variable))]
View(movement.melt)
class(movement.melt$variable)
movement.melt[, conc := factor(conc, levels = c(0, 1))]
dechor_color_scale <- scale_color_manual(
  name = "",
  values = c("#000000", "#73FAB9"),
  labels = c("dechorionated", "chorionated")
)

dechor_fill_scale <- scale_fill_manual(
  name = "",
  values = c("#000000", "#73FAB9"),
  labels = c("dechorionated", "chorionated")
)

ggplot(movement.melt, aes(x = variable, y = value, color = conc)) +
  geom_vline(xintercept = c(30, 40), color = "red", size = 0.3) +
  stat_summary(fun.y = mean, geom = "line", aes(group = conc)) +
  coord_cartesian(xlim = c(21, 48)) +
  dechor_color_scale +
  labs(x = "Time", y = "Movement")

window.data <- movement.melt[variable > 29 & variable < 41]
aucs <- window.data[
  ,
  .(
    conc = conc[1],
    AUC = DescTools::AUC(x = variable, y = value)
    ),
  by = well
  ]

wilcox.test(AUC ~ conc, data = aucs)
wilcox.test(x = aucs[conc == 0]$AUC, y = aucs[conc == 1]$AUC)

ggplot(movement.melt, aes(x = variable, y = value, color = conc)) +
  geom_vline(xintercept = c(30, 40), color = "red", size = 0.3) +
  geom_line(aes(group = well), alpha = 0.2) +
  stat_summary(fun.y = mean, geom = "line", aes(group = conc)) +
  coord_cartesian(xlim = c(21, 48)) +
  dechor_color_scale +
  labs(x = "Time", y = "Movement")

ggplot(aucs, aes(x = AUC)) +
  geom_freqpoly(aes(color = conc), binwidth = 50) +
  dechor_color_scale

ggplot(aucs, aes(x = AUC)) +
  geom_histogram(aes(fill = conc), binwidth = 50, position = "dodge") +
  dechor_fill_scale

# ggplot(aucs, aes(x = AUC)) +
#   geom_histogram(aes(fill = conc), binwidth = 50, position = "dodge", alpha = 0.4) +
#   geom_freqpoly(aes(color = conc), binwidth = 50) +
#   dechor_fill_scale +
#   dechor_color_scale

######

development <- as.data.table(
  read.csv("morphology_FOR_GUI_2019-03-07_Single_Tab_File  03-07-2019 145240.csv")
  )
development <- development[DNC_ == 0]
development <- development[, DNC_ := NULL]
View(development)

dev.melt <- melt(development, id.vars = names(development)[1:6])
dev.melt[, conc := factor(conc, levels = c(0, 1))]

ggplot(dev.melt, aes(x = variable, fill = conc)) +
  geom_bar(aes(weight = value), position = "dodge") +
  scale_fill_manual(
    name = "",
    values = c("#000000", "#73FAB9"),
    labels = c("dechorionated", "chorionated")
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Developmental marker")

cont.tbl <- xtabs(value ~ variable + conc, data = dev.melt)
cont.tbl <- cont.tbl[rowSums(cont.tbl) > 0, ]
chisq.test(cont.tbl)
BayesFactor::contingencyTableBF(cont.tbl, sampleType = "jointMulti")
names(development)
deform.cols <- c("DP24", names(development)[12:length(names(development))])
merged.dev <- development[
  ,
  .(
    conc,
    All.mort = ifelse(sum(MO24, MORT, na.rm = T) > 0, "dead", "alive"),
    All.deform = ifelse(sum(.SD, na.rm = T) > 0, "deformed", "normal")),
  .SDcols = deform.cols,
  by = well]

(cont.mort <- table(merged.dev[, .(conc, All.mort)]))
chisq.test(cont.mort)
BayesFactor::contingencyTableBF(cont.mort, sampleType = "jointMulti")

(cont.deform <- table(merged.dev[, .(conc, All.deform)]))
chisq.test(cont.deform)
BayesFactor::contingencyTableBF(cont.deform, sampleType = "jointMulti")
