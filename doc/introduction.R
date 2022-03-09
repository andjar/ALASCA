## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%", 
  fig.asp = 0.7,
  fig.width = 12,
  fig.align = "center"
)

## ----setup--------------------------------------------------------------------
library(RMASCA)

## ----dummy_data_configuration-------------------------------------------------
nPart <- 600
nGroups <- c("Controls", "Chocolate", "Salad")
nTime <- c(1, 2, 3, 4)
variables <- c("BMI", "Glucose", "VLDL", "LDL", "HDL", "ferritin", "CRP", "Happiness", "Anger", "Age")

## ----dummy_data---------------------------------------------------------------
df <- data.frame(
  partid = c(1:nPart),
  group = nGroups[sample(c(1:length(unique(nGroups))), nPart, replace = TRUE)]
)
df$time <- nTime[1]
df_temp_temp <- df
for(i in nTime[2:length(nTime)]){
  df_temp <- df_temp_temp
  df_temp$time <- i
  df <- rbind(df, df_temp)
}
for(i in 1:length(variables)){
  df[,3+i] <- rnorm(nrow(df), mean = 10, sd = 3)
}
colnames(df) <- c("partid", "group", "time", variables)

## ----dummy_data_trends--------------------------------------------------------
df$Anger[df$group == "Chocolate" & df$time == 2] <- df$Anger[df$group == "Chocolate" & df$time == 1]/2
df$Anger[df$group == "Chocolate" & df$time == 3] <- df$Anger[df$group == "Chocolate" & df$time == 1]/1.5
df$Happiness[df$group == "Chocolate" & df$time == 2] <- df$Happiness[df$group == "Chocolate" & df$time == 1]*2
df$Happiness[df$group == "Chocolate" & df$time == 3] <- df$Happiness[df$group == "Chocolate" & df$time == 1]*1.5
df$BMI[df$group == "Chocolate"] <- df$BMI[df$group == "Chocolate"]*df$time[df$group == "Chocolate"]
df$Glucose[df$group == "Chocolate"] <- df$Glucose[df$group == "Chocolate"]*df$time[df$group == "Chocolate"]

df$Anger[df$group == "Salad" & df$time == 2] <- df$Anger[df$group == "Salad" & df$time == 1]*2
df$Anger[df$group == "Salad" & df$time == 3] <- df$Anger[df$group == "Salad" & df$time == 1]*1.5
df$Happiness[df$group == "Salad" & df$time == 2] <- df$Happiness[df$group == "Salad" & df$time == 1]/2
df$Happiness[df$group == "Salad" & df$time == 3] <- df$Happiness[df$group == "Salad" & df$time == 1]/1.5
df$BMI[df$group == "Salad"] <- df$BMI[df$group == "Salad"]/df$time[df$group == "Salad"]
df$Glucose[df$group == "Salad"] <- df$Glucose[df$group == "Salad"]/df$time[df$group == "Salad"]

df$Age <- df$Age*df$time
df$HDL <- df$HDL*df$time

## ----data_dummy_melt----------------------------------------------------------
df$group <- factor(df$group)
df$group <- relevel(df$group, ref = "Controls")
df <- reshape2::melt(df, id.vars = c("partid", "group", "time"))
for(i in unique(df$variables)){
  for(j in unique(df$partid)){
    df$value[df$variables == i & df$partid == j] <- df$value[df$variables == i & df$partid == j] + sample(0:5,1)
  }
}

## ----data_dummy_plot----------------------------------------------------------
ggplot2::ggplot(data = df,
                ggplot2::aes(x = factor(time), y = value, color = group)) +
  ggplot2::geom_boxplot() +
  ggplot2::facet_wrap(~variable, scales = "free_y") +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(x = "Time", y = "Value", color = "Group")

## ----data_dummy_plot_2--------------------------------------------------------
ggplot2::ggplot(data = subset(df, partid %in% sample(unique(partid), 5)),
     ggplot2::aes(x = factor(time), y = value, group = partid, color = group)) +
     ggplot2::geom_point() + ggplot2::geom_line() +
     ggplot2::facet_wrap(~variable, scales = "free_y") +
     ggplot2::theme(legend.position = "bottom") +
     ggplot2::labs(x = "Time", y = "Value", color = "Group")

## ----RMASCA_simple_model------------------------------------------------------
model.formula <- value ~ time*group + (1|partid)
res.simple <- RMASCA(df = df, formula = model.formula)

## ----RMASCA_simple_model_plot_PC1---------------------------------------------
plot(res.simple)

## ----RMASCA_simple_model_plot_PC2---------------------------------------------
plot(res.simple, component = "PC2")

## ----RMASCA_screeplot---------------------------------------------------------
screeplot(res.simple)

## ----RMASCA_custom_plot-------------------------------------------------------
scores <- get_scores(res.simple)
loadings <- get_loadings(res.simple)

gst_1 <- ggplot2::ggplot(scores$time, ggplot2::aes(x = time, y = PC1, group = NA)) + 
  ggplot2::geom_point() + 
  ggplot2::geom_line() + 
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Time", y = paste0("Pr Comp. 1 (",round(100*scores$explained$time[1], 2),"%)"))

gst_2 <- ggplot2::ggplot(scores$time, ggplot2::aes(x = time, y = PC2, group = NA)) + 
  ggplot2::geom_point() + 
  ggplot2::geom_line() + 
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Time", y = paste0("Pr Comp. 2 (",round(100*scores$explained$time[2], 2),"%)"))

gsg_1 <- ggplot2::ggplot(scores$group, ggplot2::aes(x = time, y = PC1, group = group, color = group)) + 
  ggplot2::geom_point() + 
  ggplot2::geom_line() + 
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Time", y = paste0("Pr Comp. 1 (",round(100*scores$explained$group[1], 2),"%)"))

gsg_2 <- ggplot2::ggplot(scores$group, ggplot2::aes(x = time, y = PC2, group = group, color = group)) + 
  ggplot2::geom_point() + 
  ggplot2::geom_line() + 
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Time", y = paste0("Pr Comp. 2 (",round(100*scores$explained$group[2], 2),"%)"))

glt_1 <- ggplot2::ggplot(loadings$time, ggplot2::aes(x = covars, y = PC1)) + 
  ggplot2::geom_point() + 
  ggplot2::labs(x = "Variable", y = paste0("Pr Comp. 1 (",round(100*loadings$explained$time[1], 2),"%)")) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.4, hjust=1))

glt_2 <- ggplot2::ggplot(loadings$time, ggplot2::aes(x = covars, y = PC2)) + 
  ggplot2::geom_point() + 
  ggplot2::labs(x = "Variable", y = paste0("Pr Comp. 2 (",round(100*loadings$explained$time[2], 2),"%)")) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.4, hjust=1))

glg_1 <- ggplot2::ggplot(loadings$group, ggplot2::aes(x = covars, y = PC1)) + 
  ggplot2::geom_point() + 
  ggplot2::labs(x = "Variable", y = paste0("Pr Comp. 1 (",round(100*loadings$explained$group[1], 2),"%)")) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.4, hjust=1))

glg_2 <- ggplot2::ggplot(loadings$group, ggplot2::aes(x = covars, y = PC2)) + 
  ggplot2::geom_point() + 
  ggplot2::labs(x = "Variable", y = paste0("Pr. Comp. 2 (",round(100*loadings$explained$group[2], 2),"%)")) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.4, hjust=1))

ggpubr::ggarrange(
  ggpubr::ggarrange(gst_1, glt_1, gsg_1, glg_1, nrow = 1, widths = c(2,3,2,3), common.legend = TRUE, legend = "none"),
  ggpubr::ggarrange(gst_2, glt_2, gsg_2, glg_2, nrow = 1, widths = c(2,3,2,3), common.legend = TRUE, legend = "bottom"),
  nrow = 2
)

## ----robustness---------------------------------------------------------------
res.simple$n_validation_runs <- 50 # the more, the better

res.simple <- validate(res.simple, participant_column = "partid")

plot(res.simple)

plot(res.simple, component = "PC2")

## ----RMASCA_model_2-----------------------------------------------------------
model.formula <- value ~ time + time:group + (1|partid)
res.mod2 <- RMASCA(df = df, formula = model.formula)
plot(res.mod2)
plot(res.mod2, component = "PC2")

## ----RMASCA_model_3-----------------------------------------------------------
model.formula <- value ~ time + time:group + (1|partid)
res.mod3 <- RMASCA(df = df, formula = model.formula, separate_time_and_group = FALSE)
plot(res.mod3)
plot(res.mod3, component = "PC2")

## ----RMASCA_model_4-----------------------------------------------------------
age <- subset(df, variable == "Age")$value
df <- subset(df, variable != "Age")
df$age <- age
model.formula <- value ~ time + time:group + age + (1|partid)
res.mod4 <- RMASCA(df = df, formula = model.formula)
plot(res.mod4)
plot(res.mod4, component = "PC2")

## ----lmm_info-----------------------------------------------------------------
summary(res.simple$lmer.models[[1]])

## ----lmm_residuals------------------------------------------------------------
plot(density(residuals(res.simple$lmer.models[[1]])), main = unique(df$variable)[1])

qqnorm(residuals(res.simple$lmer.models[[1]]), main = unique(df$variable)[1]) 
qqline(residuals(res.simple$lmer.models[[1]]))

## ----lmm_pvalues--------------------------------------------------------------
head(res.simple$LMM.coefficient)

