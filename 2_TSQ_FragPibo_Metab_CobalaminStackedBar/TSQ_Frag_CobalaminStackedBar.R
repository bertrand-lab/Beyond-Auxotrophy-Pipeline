


# Setup -------------------------------------------------------------------


# Load required packages
library(tidyverse)
library(gplots)
library(tibble)
library(ggpubr)




# Frag Cobalamin Stacked Barplots ---------------------------------------



# Working directory 
setwd("~/Dropbox (Bertrand Lab)/Bertrand Lab's shared workspace/Catalina/Summer_2022/1_FracyPibo_Beyond_Auxotrophy/BA_Manuscript/BA_FragPibo_Quotas_Pipeline")

# Go to cal curve folder 
setwd("1_TSQ_Frag_Metab_CalCurveQuant")


# Frag Data Cleanup ------------------------------------------------------------


# Load and clean matab quant export dataset
frag_cobal_quant_data <-
  read.csv("metab_quant_export_df.csv") %>% # Load csv
  filter(grepl("B12-", Molecule.Name)) %>%
  dplyr::select(Replicate.Name,
                # select only needed columns
                Molecule.Name,
                B12.Treatment,
                fmol_mgC,
                Final_Peak,
                fmol_cell) %>%
  mutate(nanomol_molC = ((fmol_mgC * 1000) * 12.0107) / 1e6,
         attomol_cell = fmol_cell * 1000) %>% # create a column with converted nanomoles and zeptomoles
  # filter(!grepl("B12-CN", Molecule.Name)) %>%
  mutate(B12.Treatment = factor(B12.Treatment, levels = c("+", "-"))) %>% # convert to factors
  mutate(Replicate.Name = factor(
    Replicate.Name,
    levels = c("6", "24", "28", "15", "16", "27"),
    
    labels = c("a", "b", "c", "d", "e", "f") # give samples pleasing labels
  )) %>%
  group_by(Molecule.Name, Replicate.Name, B12.Treatment) %>%
  summarise(mean_nmol_molC = mean(nanomol_molC),
            # calculate mean nmol_mol c by replicate
            mean_Peak = mean(Final_Peak),
            mean_attomoles_cell = mean(attomol_cell)) # %>%
# filter(!grepl("-", B12.Treatment))




# Bring in reccomended norm types to calculate total cobalamin 
recc_norm <- read.csv("reccomended_norm_df.csv") %>%
  filter(grepl("B12-", Molecule.Name)) 

# Df for determining total cobalamin
total_cobalamin_df <- frag_cobal_quant_data %>%
  left_join(recc_norm, by = "Molecule.Name") %>%
  filter(reccomended_norm == "abs quant") %>%
  group_by(Replicate.Name, B12.Treatment) %>%
  summarise(total_nmol_molC = sum(mean_nmol_molC),
            total_attomoles_cell = sum(mean_attomoles_cell)) %>%
  group_by(B12.Treatment) %>%
  summarise(mean_total_nmol_molC = mean(total_nmol_molC),
            sd_total_nmol_molC = sd(total_nmol_molC),
            mean_total_attomoles_cell = mean(total_attomoles_cell),
            sd_total_attomoles_cell = sd(total_attomoles_cell))


# Frag Bring in LOD/LOQ Info ---------------------------------------------------


# Loop for checking LOD/LOQ
setwd("../0_TSQ_Frag_Metab_BMIS")

# load in df with LOD/LOQ's from BMIS script
LODQ_df <- read.csv("LODQ_export.csv") %>%
  dplyr::select(Molecule.Name, LOD, LOQ) 

# Add lod/loq values to the quant data 
frag_cobal_quant_data_LODQ <- 
  left_join(frag_cobal_quant_data, LODQ_df, by = "Molecule.Name") %>%
  mutate(Below_LOD = NA) %>%
  mutate(Below_LOQ = NA) %>%
  mutate(annotate = NA)


# Loop to notify if below LOD
for (i in 1:nrow(frag_cobal_quant_data_LODQ)){
  if (frag_cobal_quant_data_LODQ$mean_Peak[i] < frag_cobal_quant_data_LODQ$LOD[i]){
    frag_cobal_quant_data_LODQ$Below_LOD[i] <- TRUE
    frag_cobal_quant_data_LODQ$mean_nmol_molC[i] <- 0
    frag_cobal_quant_data_LODQ$mean_attomoles_cell[i] <- 0
  }
  else{
    frag_cobal_quant_data_LODQ$Below_LOD[i] <- FALSE
  }
}
  
# Loop to notify if sample is below LOQ
for (i in 1:nrow(frag_cobal_quant_data_LODQ)){
  if (frag_cobal_quant_data_LODQ$mean_Peak[i] < frag_cobal_quant_data_LODQ$LOQ[i]){
    frag_cobal_quant_data_LODQ$Below_LOQ[i] <- TRUE
  }
  else{
    frag_cobal_quant_data_LODQ$Below_LOQ[i] <- FALSE
  }
}


# If above LOD but below LOQ, annotate with a star
for (i in 1:nrow(frag_cobal_quant_data_LODQ)) {
  if (frag_cobal_quant_data_LODQ$Below_LOD[i] == "FALSE" &
      frag_cobal_quant_data_LODQ$Below_LOQ[i] == "TRUE") {
     frag_cobal_quant_data_LODQ$annotate[i] <- "*"
  }
  else{}
}


# Frag Plotting ----------------------------------------------------------------


# color palettes from metbrewer
tiep_palette <- met.brewer("Tiepolo", n=6)
b12_palette <- c(tiep_palette[2], tiep_palette[1], tiep_palette[3], tiep_palette[5])

# prep data 
frag_plot_data <- frag_cobal_quant_data_LODQ %>% 
    filter(is.na(annotate)) %>%
    mutate(Molecule.Name = factor(Molecule.Name, 
                                  levels = c("B12-Ado",
                                             "B12-CN",
                                             "B12-Me",
                                             "B12-OH")))

# Set heights for trace signifying asterisks
yloc <- c(120, 0)

# Create df for trace labels 
frag_trace_label <- frag_cobal_quant_data_LODQ %>% 
  filter(!is.na(annotate)) 

# Add trace labels to df
frag_trace_label$yloc <- yloc
  


# Plot for stacked bar with nanomoles B12/mol C in Frag 
frag_CobalStackedbar_perC <- ggplot() +
  geom_bar(
    data = frag_plot_data,
    aes(
      fill = Molecule.Name,
      y = as.numeric(mean_nmol_molC),
      x = Replicate.Name
    ),
    position = "stack",
    stat = "identity"
  ) +
    theme_classic() +
  ylab(expression(paste("Nanomoles B"["12"] ~ "mole C" ^ "-1"))) +
  xlab (expression("Biological Replicate")) +
  scale_fill_manual(
    name = expression(paste("B"["12"] * " Analog")),
    values = b12_palette,
    labels = c(expression(paste("Ado-B"["12"])), 
               expression(paste("CN-B"["12"])),
               expression(paste("Me-B"["12"])),
               expression(paste("OH-B"["12"])))) +
  theme(panel.background = element_rect(fill = "transparent"),
    # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA),
    # bg of the plot
    panel.grid.major = element_blank(),
    # get rid of major grid
    panel.grid.minor = element_blank(),
    # get rid of minor grid
    legend.background = element_rect(fill = "transparent"),
    # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    # get rid of legend panel bg
    text = element_text(size = 24, color = "black")) +
  facet_grid(. ~ B12.Treatment, 
             scales = "free_x", drop = TRUE) +
  geom_text(data = frag_trace_label, 
            aes(x = Replicate.Name,
                y = yloc,
                label = "*"
                ),
            size = 25,
            color = b12_palette[1]) +
  ylim(0,300)




# New heights for trace signifying asterisks
yloc <- c(.12, 0)

# Add trace labels to df
frag_trace_label$yloc <- yloc


# Plot for stacked bar with attomoles B12/cell in Frag 
frag_CobalStackedbar_percell <- ggplot() +
  geom_bar(
    data = frag_plot_data,
    aes(
      fill = Molecule.Name,
      y = as.numeric(mean_attomoles_cell),
      x = Replicate.Name
    ),
    position = "stack",
    stat = "identity"
  ) +
  theme_classic() +
  ylab(expression(paste("Attomoles B"["12"] ~ "Cell" ^ "-1"))) +
  xlab (expression("Biological Replicate")) +
  scale_fill_manual(
    name = expression(paste("B"["12"] * " Analog")),
    values = b12_palette,
    labels = c(expression(paste("Ado-B"["12"])), 
               expression(paste("CN-B"["12"])),
               expression(paste("Me-B"["12"])),
               expression(paste("OH-B"["12"])))) +
  theme(panel.background = element_rect(fill = "transparent"),
        # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        # bg of the plot
        panel.grid.major = element_blank(),
        # get rid of major grid
        panel.grid.minor = element_blank(),
        # get rid of minor grid
        legend.background = element_rect(fill = "transparent"),
        # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        # get rid of legend panel bg
        text = element_text(size = 24, color = "black")) +
  facet_grid(. ~ B12.Treatment,
             scales = "free_x", drop = TRUE) +
  geom_text(data = frag_trace_label,
            aes(x = Replicate.Name,
                y = yloc,
                label = "*"
            ),
            size = 25,
            color = b12_palette[1]) 
  
  
# Change wd to cobalamin folder
setwd("../2_TSQ_Frag_Metab_CobalaminStackedBar/")


# Open a pdf file
pdf("Frag_cobalstackedbar_nmol_molC.pdf") 

# print plot 
frag_CobalStackedbar_perC

# Close the pdf file
dev.off() 

# Open a pdf file
pdf("Frag_cobalstackedbar_attomol_cell.pdf") 

# print plot 
frag_CobalStackedbar_percell

# Close the pdf file
dev.off() 




# Pibo Cobalamin Stacked Bar Plots  --------------------------------------


# Go to cal curve folder 
setwd("../TSQ_Pibo_Metab_Sumdata")


# Pibo Data Cleanup ------------------------------------------------------------

pibo_cobal_quant_data <- read.csv("BA_B12Quoatas_AG.csv") %>%
  # filter(!grepl("B12-CN", Molecule.Name)) %>% #Remove cyanocobalamin
  filter(!grepl("cobalamin", Molecule.Name)) %>% # remove total cobalamin entries
  dplyr::select(
    Bioreplicate,
    Molecule.Name,
    Treatment,
    fmol.oncolumn,
    LOD.fmol.column,
    LOQ.fmol.column,
    molecules.per.cell,
    nmol.cobalamin.per.bio.rep,
    mol.carbon.per.sample
  ) %>% # keep only relevant columns
  mutate(rep.label = paste(Bioreplicate, 
                                Treatment)) %>% # put together biorep and treatment
  mutate(LOQ.fmol.column = LOD.fmol.column/3*5) %>% # Change LOQ to 5 times sd rather than 10

  mutate(rep.label = factor(
    rep.label,
    levels = c("A +", "B +", "C +", "A -", "B -", "C -"),
    labels = c("a", "b", "c", "d", "e", "f")
  )) %>%
  rownames_to_column(var = "index")




# Pibo Bring in LOD/LOQ Info ---------------------------------------------------

# load in df with LOD/LOQ's from BMIS script
LODQ_df_pibo <- read.csv("BA_B12Quoatas_AG.csv") %>%  
  filter(!grepl("cobalamin", Molecule.Name)) %>% # remove total cobalamin
  dplyr::select(fmol.oncolumn,
                LOD.fmol.column,
                LOQ.fmol.column) %>%
  mutate(LOQ.fmol.column = LOD.fmol.column/3*5) %>% # swap LOQ to 5x sd rather than 10
  rownames_to_column(var = "index")

# Loop to notify if below LOD
for (i in 1:nrow(LODQ_df_pibo)){
  if (LODQ_df_pibo$fmol.oncolumn[i] < LODQ_df_pibo$LOD.fmol.column[i]){
    LODQ_df_pibo$Below_LOD[i] <- TRUE
  }
  else{
    LODQ_df_pibo$Below_LOD[i] <- FALSE
  }
}

# Loop to notify if sample is below LOQ
for (i in 1:nrow(LODQ_df_pibo)){
  if (LODQ_df_pibo$fmol.oncolumn[i] < LODQ_df_pibo$LOQ.fmol.column[i]){
    LODQ_df_pibo$Below_LOQ[i] <- TRUE
  }
  else{
    LODQ_df_pibo$Below_LOQ[i] <- FALSE
  }
}

# If above LOD but below LOQ, annotate with a star
for (i in 1:nrow(LODQ_df_pibo)) {
  if (LODQ_df_pibo$Below_LOD[i] == "FALSE" &
      LODQ_df_pibo$Below_LOQ[i] == "TRUE") {
    LODQ_df_pibo$annotate[i] <- "*"
  }
  else{LODQ_df_pibo$annotate[i] <- "NA"}
}

# Clean up for only needed info 
LODQ_df_pibo_clean <- LODQ_df_pibo %>%
  dplyr::select(index, 
                Below_LOD, 
                Below_LOQ, 
                annotate)

# Merge back with original dataset 
pibo_cobal_quant_data_LODQ <-
  left_join(pibo_cobal_quant_data, 
            LODQ_df_pibo_clean, 
            by = "index") %>%
  group_by(Molecule.Name, 
           rep.label, 
           Treatment) %>% # group by
  summarise(
    nmolcobal = first(nmol.cobalamin.per.bio.rep),
    molC = first(mol.carbon.per.sample),
    molecules.per.cell = mean(molecules.per.cell),
    Below_LOD = first(Below_LOD), 
    Below_LOQ = first(Below_LOQ),
    annotate = first(annotate)
  ) %>% # create rows with nmolcobalamin and mol C
  mutate(mean_nmol_molC = nmolcobal / molC,
         mean_molecules.per.cell = mean(molecules.per.cell)) %>% # calculate nmolcobalamin per mole Carbon
  mutate(Treatment = factor(Treatment, levels = c("+", "-")), # Change treatment to factor
         mean_nmol_molC = ifelse(Below_LOD == TRUE, 0, mean_nmol_molC), # if below LOD, set values to 0
         mean_molecules.per.cell = ifelse(Below_LOD == TRUE, 0, mean_molecules.per.cell))
  

# FIXME plot pibo with same units as frag and in "best units" as well 
# Pibo Plotting ----------------------------------------------------------------

# Prep plot data
pibo_plot_data <- pibo_cobal_quant_data_LODQ %>%
  mutate(Molecule.Name = factor(Molecule.Name, 
                                levels = c("B12-Ado",
                                           "B12-CN",
                                           "B12-Me",
                                           "B12-OH"))) %>%
  rownames_to_column(var = "index") %>%
  filter(!(annotate == "*")) # remove trace elements 

# Create df with labels for trace entries 
trace_label <- pibo_cobal_quant_data_LODQ  %>%
  rownames_to_column(var = "index") %>%
  filter(index == 14) # Keep only entry 14 to use as trace 

# Stacked barplot for pibo (per C)
pibo_CobalStackedbar_perC <- ggplot() +
  geom_bar(
    data = pibo_plot_data,
    aes(fill = Molecule.Name, 
        y = mean_nmol_molC, 
        x = rep.label),
    position = "stack",
    stat = "identity"
  ) +
  theme_classic() +
  ylab(expression(paste("Nanomoles B"["12"] ~ "mole C" ^ "-1"))) +
  xlab (expression("Biological Replicate")) +
  scale_fill_manual(
    name = expression(paste("B"["12"] * " Analog")),
    values = b12_palette,
    labels = c(expression(paste("Ado-B"["12"])), 
               expression(paste("CN-B"["12"])), 
               expression(paste("Me-B"["12"])), 
               expression(paste("OH-B"["12"])))
  ) +
  theme(
    panel.background = element_rect(fill = "transparent"),
    # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA),
    # bg of the plot
    panel.grid.major = element_blank(),
    # get rid of major grid
    panel.grid.minor = element_blank(),
    # get rid of minor grid
    legend.background = element_rect(fill = "transparent"),
    # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    # get rid of legend panel bg
    text = element_text(size = 23, color = "black")
  ) +
  facet_grid(. ~ Treatment, scales = "free_x", drop = TRUE) +
  geom_text(data = trace_label, 
            aes(
    x = "b", 
    y = 110, 
    label = "*"),
    size = 25,
    color = b12_palette[3]) +
  ylim(0, 300)

# Stacked barplot for pibo (attomoles per cell)
pibo_CobalStackedbar_atto_percell <- ggplot() +
  geom_bar(
    data = pibo_plot_data,
    aes(fill = Molecule.Name, 
        y = (mean_molecules.per.cell / 6.02E23)* 1E18 , # convert molecules to attomoles
        x = rep.label),
    position = "stack",
    stat = "identity"
  ) +
  theme_classic() +
  ylab(expression(paste("Attomoles B"["12"] ~ "Cell" ^ "-1"))) +
  xlab (expression("Biological Replicate")) +
  scale_fill_manual(
    name = expression(paste("B"["12"] * " Analog")),
    values = b12_palette,
    labels = c(expression(paste("Ado-B"["12"])), 
               expression(paste("CN-B"["12"])), 
               expression(paste("Me-B"["12"])), 
               expression(paste("OH-B"["12"])))
  ) +
  theme(
    panel.background = element_rect(fill = "transparent"),
    # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA),
    # bg of the plot
    panel.grid.major = element_blank(),
    # get rid of major grid
    panel.grid.minor = element_blank(),
    # get rid of minor grid
    legend.background = element_rect(fill = "transparent"),
    # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    # get rid of legend panel bg
    text = element_text(size = 23, color = "black")
  ) +
  facet_grid(. ~ Treatment, scales = "free_x", drop = TRUE) +
  geom_text(data = trace_label, 
            aes(
              x = "b", 
              y = .28E-03, 
              label = "*"),
            size = 25,
            color = b12_palette[3]) +
  ylim(0, 5E-04)

# Stacked barplot for pibo (per cell)
pibo_CobalStackedbar_zepto_percell <- ggplot() +
  geom_bar(
    data = pibo_plot_data,
    aes(fill = Molecule.Name, 
        y = (mean_molecules.per.cell / 6.02E23)* 1E21 , # convert molecules to attomoles
        x = rep.label),
    position = "stack",
    stat = "identity"
  ) +
  theme_classic() +
  ylab(expression(paste("Zeptomoles B"["12"] ~ "Cell" ^ "-1"))) +
  xlab (expression("Biological Replicate")) +
  scale_fill_manual(
    name = expression(paste("B"["12"] * " Analog")),
    values = b12_palette,
    labels = c(expression(paste("Ado-B"["12"])), 
               expression(paste("CN-B"["12"])), 
               expression(paste("Me-B"["12"])), 
               expression(paste("OH-B"["12"])))
  ) +
  theme(
    panel.background = element_rect(fill = "transparent"),
    # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA),
    # bg of the plot
    panel.grid.major = element_blank(),
    # get rid of major grid
    panel.grid.minor = element_blank(),
    # get rid of minor grid
    legend.background = element_rect(fill = "transparent"),
    # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    # get rid of legend panel bg
    text = element_text(size = 23, color = "black")
  ) +
  facet_grid(. ~ Treatment, scales = "free_x", drop = TRUE) +
  geom_text(data = trace_label, 
            aes(
              x = "b", 
              y = .28, 
              label = "*"),
            size = 25,
            color = b12_palette[3]) +
  ylim(0, 5E-01)

# Change wd to cobalamin folder
setwd("../2_TSQ_Frag_Metab_CobalaminStackedBar/")


# Open a pdf file
pdf("Pibo_cobalstackedbar_nmol_molC.pdf") 

# print plot 
pibo_CobalStackedbar

# Close the pdf file
dev.off() 


# Open a pdf file
pdf("PiboFrag_faceted_cobalstackedbar_nmol_molC.pdf") 


ggarrange(
  frag_CobalStackedbar,
  pibo_CobalStackedbar,
  labels = c("a)", "b)"),
  common.legend = TRUE
)

# Close the pdf file
dev.off() 


# Boxplots ----------------------------------------------------------------


# Boxplot data
frag_boxplot_df <- frag_cobal_quant_data %>% 
  filter(!grepl("-", B12.Treatment))


# Frag boxplot
ggplot(frag_cobal_quant_data, aes(x = Molecule.Name, y = mean_nmol_molC, fill = Molecule.Name)) + 
  geom_boxplot(notch=FALSE) +
  ylab(expression(paste("Nanomoles B"["12"] ~ "mole C" ^ "-1"))) 
  


# Boxplot data
frag_boxplot_df <- frag_cobal_quant_data %>% 
  filter(!grepl("-", B12.Treatment))


# Frag boxplot
frag_boxplot <- ggplot(frag_cobal_quant_data, aes(x = Molecule.Name, y = mean_nmol_molC, fill = Molecule.Name)) + 
  geom_boxplot(notch=FALSE) +
  ylab(expression(paste("Nanomoles B"["12"] ~ "mole C" ^ "-1"))) +
  theme_classic() +
  scale_fill_manual(
    name = expression(paste("B"["12"] * " Analog")),
    values = b12_palette,
    labels = c(expression(paste("Ado-B"["12"])), expression(paste("Me-B"["12"])), expression(paste("OH-B"["12"])))
  ) +
  theme(
    panel.background = element_rect(fill = "transparent"),
    # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA),
    # bg of the plot
    panel.grid.major = element_blank(),
    # get rid of major grid
    panel.grid.minor = element_blank(),
    # get rid of minor grid
    legend.background = element_rect(fill = "transparent"),
    # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    # get rid of legend panel bg
    text = element_text(size = 23, color = "black"),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank()
  ) 
  
# Boxplot data
pibo_boxplot_df <- pibo_cobal_quant_data %>% 
  filter(!grepl("-", Treatment))


# Frag boxplot
pibo_boxplot <- ggplot(pibo_cobal_quant_data, aes(x = Molecule.Name, y = mean_nmol_molC, fill = Molecule.Name)) + 
  geom_boxplot(notch=FALSE) +
  ylab(expression(paste("Nanomoles B"["12"] ~ "mole C" ^ "-1"))) +
  theme_classic() +
  scale_fill_manual(
    name = expression(paste("B"["12"] * " Analog")),
    values = b12_palette,
    labels = c(expression(paste("Ado-B"["12"])), expression(paste("Me-B"["12"])), expression(paste("OH-B"["12"])))
  ) +
  theme(
    panel.background = element_rect(fill = "transparent"),
    # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA),
    # bg of the plot
    panel.grid.major = element_blank(),
    # get rid of major grid
    panel.grid.minor = element_blank(),
    # get rid of minor grid
    legend.background = element_rect(fill = "transparent"),
    # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"),
    # get rid of legend panel bg
    text = element_text(size = 23, color = "black"),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank()
  ) 

ggarrange(
  frag_boxplot,
  pibo_boxplot,
  labels = c("a)", "b)"),
  common.legend = TRUE
)

