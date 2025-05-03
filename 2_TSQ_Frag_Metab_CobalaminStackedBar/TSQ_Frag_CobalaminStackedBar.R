


# Setup -------------------------------------------------------------------


# Load required packages
library(tidyverse)
library(gplots)
library(tibble)
library(ggpubr)
library(MetBrewer)
library(janitor)
library(flextable)
library(here)



# Frag Data Cleanup ------------------------------------------------------------


# Load and clean matab quant export dataset
frag_cobal_quant_data <-
  read.csv(here("1_TSQ_Frag_Metab_CalCurveQuant", "metab_quant_export_df.csv")) |> # Load csv
  filter(grepl("B12-", Molecule.Name)) |>
  dplyr::select(Replicate.Name,
                # select only needed columns
                Molecule.Name,
                B12.Treatment,
                fmol_mgC,
                Final_Peak,
                fmol_cell) |>
  mutate(nanomol_molC = ((fmol_mgC * 1000) * 12.0107) / 1e6,
         attomoles_molC = nanomol_molC * 1e9,
         attomol_cell = fmol_cell * 1000) |> # create a column with converted nanomoles and zeptomoles
  # filter(!grepl("B12-CN", Molecule.Name)) |>
  mutate(B12.Treatment = factor(B12.Treatment, levels = c("+", "-"))) |> # convert to factors
  mutate(Replicate.Name = factor(
    Replicate.Name,
    levels = c("6", "24", "28", "15", "16", "27"),
    
    labels = c("a", "b", "c", "d", "e", "f") # give samples pleasing labels
  )) |>
  group_by(Molecule.Name, Replicate.Name, B12.Treatment) |>
  summarise(mean_nmol_molC = mean(nanomol_molC),
            mean_attomol_molC = mean(attomoles_molC),
            # calculate mean nmol_mol c by replicate
            mean_Peak = mean(Final_Peak),
            mean_attomoles_cell = mean(attomol_cell)) # |>
# filter(!grepl("-", B12.Treatment))




# Bring in recommended norm types to calculate total cobalamin 
recc_norm <- read.csv(here("1_TSQ_Frag_Metab_CalCurveQuant", "recommended_norm_df.csv")) |>
  filter(grepl("B12-", Molecule.Name)) 

# Df for determining total cobalamin
total_cobalamin_df <- frag_cobal_quant_data |>
  left_join(recc_norm, by = "Molecule.Name") |>
  filter(reccomended_norm == "abs quant") |>
  group_by(Replicate.Name, B12.Treatment) |>
  summarise(total_nmol_molC = sum(mean_nmol_molC),
            total_attomoles_molC = sum(mean_attomol_molC),
            total_attomoles_cell = sum(mean_attomoles_cell)) |>
  group_by(B12.Treatment) |>
  summarise(mean_total_nmol_molC = mean(total_nmol_molC),
            sd_total_nmol_molC = sd(total_nmol_molC),
            mean_total_attomol_molC = mean(total_attomoles_molC),
            sd_total_attomol_molC = sd(total_attomoles_molC),
            mean_total_attomoles_cell = mean(total_attomoles_cell),
            sd_total_attomoles_cell = sd(total_attomoles_cell))


# Frag Bring in LOD/LOQ Info ---------------------------------------------------


# load in df with LOD/LOQ's from BMIS script
LODQ_df <- read.csv(here("1_TSQ_Frag_Metab_CalCurveQuant", "LODQ_export.csv")) |>
  dplyr::select(Molecule.Name, LOD, LOQ) 

# Add lod/loq values to the quant data 
frag_cobal_quant_data_LODQ <- 
  left_join(frag_cobal_quant_data, LODQ_df, by = "Molecule.Name") |>
  mutate(Below_LOD = NA) |>
  mutate(Below_LOQ = NA) |>
  mutate(annotate = NA)


# Loop to notify if below LOD
for (i in 1:nrow(frag_cobal_quant_data_LODQ)){
  if (frag_cobal_quant_data_LODQ$mean_Peak[i] < frag_cobal_quant_data_LODQ$LOD[i]){
    frag_cobal_quant_data_LODQ$Below_LOD[i] <- TRUE
    frag_cobal_quant_data_LODQ$mean_nmol_molC[i] <- 0
    frag_cobal_quant_data_LODQ$mean_attomol_molC[i] <- 0
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
b12_palette <- c(tiep_palette[2], 
                 tiep_palette[1], 
                 tiep_palette[3], 
                 tiep_palette[5])

# prep data 
frag_plot_df <- frag_cobal_quant_data_LODQ |> 
    filter(is.na(annotate)) |>
    mutate(Molecule.Name = factor(Molecule.Name, 
                                  levels = c("B12-Ado",
                                             "B12-CN",
                                             "B12-Me",
                                             "B12-OH")))

# Set heights for trace signifying asterisks
yloc <- c(120, 0)

# Create df for trace labels 
frag_trace_label <- frag_cobal_quant_data_LODQ |> 
  filter(!is.na(annotate)) 

# Add trace labels to df
frag_trace_label$yloc <- yloc
  
# labels for facets
levels(frag_plot_df$B12.Treatment) <- c("+B[12]", "-B[12]")

# labels for facets
levels(frag_trace_label$B12.Treatment) <- c("+B[12]", "-B[12]")


# Plot for stacked bar with nanomoles B12/mol C in Frag 
frag_CobalStackedbar_perC <- ggplot() +
  geom_bar(
    data = frag_plot_df,
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
    name = expression(paste("B"["12"] * " Form")),
    values = b12_palette,
    labels = c(expression(paste("Ado-B"["12"])), 
               expression(paste("CN-B"["12"])),
               expression(paste("Me-B"["12"])),
               expression(paste("OH-B"["12"])))) +
  theme(panel.background = element_rect(fill = "transparent"
                                        ),
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
    text = element_text(size = 24, color = "black"),
    
    # Remove axis text
    axis.text.x=element_blank(), 
    
    # Remove ticks
    axis.ticks.x=element_blank()) +
  geom_text(data = frag_trace_label,
            aes(x = Replicate.Name,
                y = yloc,
                label = "*"),
            size = 25,
            color = b12_palette[1]) +
  facet_grid(. ~ B12.Treatment, 
             scales = "free_x", 
             drop = TRUE,
             labeller = label_parsed) +
  ylim(0,300)


# New heights for trace signifying asterisks
yloc <- c(1.2*1e11, 0)

# Add trace labels to df
frag_trace_label$yloc <- yloc



# Plot for stacked bar with attomoles B12/mol C in Frag 
frag_CobalStackedbar_attomoles_perC <- ggplot() +
  geom_bar(
    data = frag_plot_df,
    aes(
      fill = Molecule.Name,
      y = as.numeric(mean_attomol_molC),
      x = Replicate.Name
    ),
    position = "stack",
    stat = "identity"
  ) +
  theme_classic() +
  ylab(expression(paste("Attomoles B"["12"] ~ "mole C" ^ "-1"))) +
  xlab (expression("Biological Replicate")) +
  scale_fill_manual(
    name = expression(paste("B"["12"] * " Form")),
    values = b12_palette,
    labels = c(expression(paste("Ado-B"["12"])), 
               expression(paste("CN-B"["12"])),
               expression(paste("Me-B"["12"])),
               expression(paste("OH-B"["12"])))) +
  theme(panel.background = element_rect(fill = "transparent"
  ),
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
  text = element_text(size = 24, color = "black"),
  
  # Remove axis text
  axis.text.x=element_blank(), 
  
  # Remove ticks
  axis.ticks.x=element_blank()) +
  geom_text(data = frag_trace_label,
            aes(x = Replicate.Name,
                y = yloc,
                label = "*"),
            size = 25,
            color = b12_palette[1]) +
  facet_grid(. ~ B12.Treatment, 
             scales = "free_x", 
             drop = TRUE,
             labeller = label_parsed) +
  ylim(0,300*1e9)

# New heights for trace signifying asterisks
yloc <- c(.12, 0)

# Add trace labels to df
frag_trace_label$yloc <- yloc


# Plot for stacked bar with attomoles B12/cell in Frag 
frag_CobalStackedbar_percell <- ggplot() +
  geom_bar(
    data = frag_plot_df,
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
    name = expression(paste("B"["12"] * " Form")),
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
        text = element_text(size = 24, color = "black"),
        
        # Remove axis text
        axis.text.x=element_blank(), 
        
        # Remove ticks
        axis.ticks.x=element_blank()) +

  geom_text(data = frag_trace_label,
            aes(x = Replicate.Name,
                y = yloc,
                label = "*"
            ),
            size = 25,
            color = b12_palette[1]) +
  facet_grid(. ~ B12.Treatment,
             scales = "free_x", 
             drop = TRUE,
             labeller = label_parsed)
  


# Open a pdf file
pdf(here("2_TSQ_Frag_Metab_CobalaminStackedBar", "Frag_cobalstackedbar_nmol_molC.pdf"))

# print plot 
frag_CobalStackedbar_perC

# Close the pdf file
dev.off() 

# Open a pdf file
pdf(here("2_TSQ_Frag_Metab_CobalaminStackedBar","Frag_cobalstackedbar_attomol_cell.pdf"))

# print plot 
frag_CobalStackedbar_percell

# Close the pdf file
dev.off() 




# Faceted Stacked Bar Plot ------------------------------------------------


# Open a pdf file
pdf("Frag_faceted_cobalstackedbar_nmol_molC.pdf", width = 15, height = 13) 


ggarrange(
  frag_CobalStackedbar_perC,
  frag_CobalStackedbar_percell,
  labels = c("a)", "b)"),
  common.legend = TRUE,
  legend = "bottom",
  font.label = list(size = 25)
)

# Close the pdf file
dev.off() 

# B12 Totals Table (Table S3) -------------------------------------------------------------------

# Create a table for B12 data
frag_totals_df <- frag_plot_df |>
  
  # grab only B12 treatments 
  filter(B12.Treatment == "+B[12]") |>
  dplyr::select(Replicate.Name,
         Molecule.Name,
         mean_attomol_molC,
         mean_attomoles_cell) |>
  
  # Conver mol per cell to attomoles per cell
  mutate(mean_mol_cell = mean_attomoles_cell * 602214.15) |>
  group_by(Molecule.Name) |>
  summarise(
    mean_rep_mol_cell = round(mean(mean_mol_cell)),
    sd_rep_mol_cell = round(sd(mean_mol_cell)),
    mean_rep_nmol_molC = round(mean(mean_attomol_molC) / 1E9, digits = 2),
    sd_rep_nmol_molC = round(sd(mean_attomol_molC) / 1E9, digits =
                                  2),
    mean_rep_attomoles_cell = round(mean(mean_attomoles_cell), digits =
                                      3),
    sd_rep_attomoles_cell = round(sd(mean_attomoles_cell), digits =
                                    3)
  ) |>
  
  # Add totals
  adorn_totals() |>
  
  # Change to scientific notation 
  # Not using for now because looks gross with the sd's adorned as happens later 
  # mutate(across(!Molecule.Name, ~ formatC(.x, format = "e", digits = 3))) |> 
  
  # Add sd's to values
  unite("Nmol per Mole C",
        mean_rep_nmol_molC:sd_rep_nmol_molC,
        sep = "±") |>

  unite("Attomoles per Cell",
        mean_rep_attomoles_cell:sd_rep_attomoles_cell,
        sep = "±") |>
  unite("Molecules per Cell", mean_rep_mol_cell:sd_rep_mol_cell, sep = "±") |>

  
  # display as html element 
  flextable() |> 
  
  set_header_labels(values = 
                      list(Molecule.Name = "Cobamide")) |> 
  
  theme_zebra()
  
  
  frag_table




    




