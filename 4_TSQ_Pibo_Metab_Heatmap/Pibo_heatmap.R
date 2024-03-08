
# Load required packages
library(tidyverse)
library(gplots)
library(tibble)
library(MetBrewer)


# Heatmap for frag quotas

# Working directory 
setwd("~/Dropbox (Bertrand Lab)/Bertrand Lab's shared workspace/Catalina/Summer_2022/1_FracyPibo_Beyond_Auxotrophy/BA_Manuscript/BA_FragPibo_Quotas_Pipeline/TSQ_Pibo_Metab_Sumdata")



# Data cleanup
quant_data <-
  read.csv("Pibo_BA_heatmaps2.csv") %>% # Read in export file
  dplyr::select(Molecule.Name, Treatment, Bioreplicate, Area.per.cell, Molecules.per.cell)


heatmap_data <- quant_data %>% # Filter for only required cols
  unite("Treatment_Biorep", Treatment:Bioreplicate) %>% # make a column with both biorep numbers and B12 treatments
  mutate(
    Molecule.Name = replace(
      Molecule.Name,
      Molecule.Name == "Dimethyl-benzimidazole(DMB)",
      "DMB"
    )
  )  %>%
  mutate(
    Molecule.Name = replace(
      Molecule.Name,
      Molecule.Name == "FAMP",
      "FAMP*"
    )
  ) %>%
mutate(
  Molecule.Name = replace(
    Molecule.Name,
    Molecule.Name == "Homarine",
    "Homarine*"
  )
) %>%
  group_by(Molecule.Name, Treatment_Biorep) %>%
  dplyr::summarise(Area.per.cell = mean(Area.per.cell, na.rm = TRUE)) %>%
  dplyr::filter(!grepl('D', Treatment_Biorep)) %>%
  pivot_wider(names_from = Molecule.Name, values_from = Area.per.cell) %>%
  column_to_rownames("Treatment_Biorep") %>% # Change to rowname
  t() %>% # transpose df
  as.matrix() # change to matrix

# use metbrewer to generate color palette
palette <- rev(met.brewer("Hiroshige", n=100))

# Create label vector
labvec <- c(NA, expression("+B"[12]), NA, NA,expression("-B"[12]),  NA)

setwd("../4_TSQ_Pibo_Metab_Heatmap")

# Open a pdf file
pdf("Pibo_heatmap.pdf") 

# Create heatmap
heatmap.2(heatmap_data,
          Colv = FALSE,
          trace = "none", # remove trace
          density = "none", # remove density label
          col = palette, # colors
          labCol = labvec,
          cexCol = 2,
          srtCol = 0,
          adjCol = .5,
          scale = "row", # calculate z scores by row
          sepwidth = c(.1, .1), # cell border thickness
          colsep = 3, # space between the 3rd column
          sepcolor = "white",
          margins = c(7, 10))  # cell border color
dev.off()


# T-tests -----------------------------------------------------------------

# Data prep 

ttest_df <- quant_data %>%
  group_by(Molecule.Name, Treatment, Bioreplicate) %>%
  summarise(mean_peak_cell = mean(Area.per.cell),
            mean_Molecules.per.cell = mean(as.numeric(Molecules.per.cell))) %>%
  unite("Treatment_Biorep", Treatment:Bioreplicate, remove = FALSE) %>%
  mutate(rep.label = factor(
    Treatment_Biorep,
    levels = c("+_A", "+_B", "+_C", "-_A", "-_B", "-_C"),
    labels = c("a", "b", "c", "d", "e", "f"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("+", "-")))


# Create a function that checks metab differences with student's t-test
metab_ttests <- function(molecule_name){
  t.test(mean_peak_cell ~ Treatment, 
         data = ttest_df %>%
           filter(Molecule.Name == paste(molecule_name)),  
         var.equal = TRUE)
} 

# Create list of molecules to complete ttest on 
metabs_list <- unique(ttest_df$Molecule.Name)

# Apply function over list of metabs
lapply(metabs_list, metab_ttests)



# Individual Metabolite Plots ---------------------------------------------

FAMP_plot_data <- ttest_df %>%
  filter(Molecule.Name == "FAMP") 

# Check difference with student's t-test
t.test(mean_peak_cell ~ Treatment, data = FAMP_plot_data, var.equal = TRUE)
# P = 0.08027

FAMP_plot <- ggplot(data = FAMP_plot_data,
       aes(x = Treatment,
           y = mean_peak_cell)) +
  geom_point(stat = "identity",
             size = 10,
             color = "darkgrey") +
  ylab(expression(paste("Mean Peak Area Cell" ^ "-1" ~ "FAMP"))) +
  xlab("Treatment") + 
  theme_classic() +
  theme(
    text = element_text(size = 23, color = "black")
  ) +
  ylim(0,0.040) +
  geom_signif(comparisons = list(c("+", "-")),
              y_position = .034,
              annotation = ".",
              textsize = 10)

Homarine_plot_data <- ttest_df %>%
  filter(Molecule.Name == "Homarine") %>%
  group_by(Treatment, Bioreplicate) %>%
  mutate(zeptomoles_cell = (mean_Molecules.per.cell / 6.02E23)* 1e21)
  
  

# Check difference with student's t-test
t.test(zeptomoles_cell ~ Treatment, data = Homarine_plot_data, var.equal = TRUE)
# P = 0.03848 (relative peak)
# P = 0.03179 (zeptomoles)


Homarine_plot <- ggplot(data = Homarine_plot_data,
                    aes(x = Treatment,
                        y = mean_peak_cell)) +
  geom_point(stat = "identity",
             size = 10,
             color = "darkgrey") +
  ylab(expression(paste("Mean Peak Cell" ^ "-1" ~ "Homarine"))) +
  xlab("Treatment") + 
  theme_classic() +
  theme(
    text = element_text(size = 23, color = "black")
  ) 
   ylim(0,0.00030) +
  geom_signif(comparisons = list(c("+", "-")),
              map_signif_level = TRUE, 
               y_position = .00028,
              test = "t.test")


Homarine_plot_zeptomoles_cell <- ggplot(data = Homarine_plot_data,
                        aes(x = Treatment,
                            y = zeptomoles_cell)) +
  geom_point(stat = "identity",
             size = 10,
             color = "darkgrey") +
  ylab(expression(paste("Zeptomoles Cell" ^ "-1" ~ "Homarine"))) +
  xlab("Treatment") + 
  theme_classic() +
  theme(
    text = element_text(size = 23, color = "black")
  ) +
  ylim(0,10) +
  geom_signif(comparisons = list(c("+", "-")),
              annotations = "*",
              textsize = 10,
              y_position = 8.5)


ggarrange(Homarine_plot_zeptomoles_cell, FAMP_plot) 




