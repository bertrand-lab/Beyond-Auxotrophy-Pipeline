# This script creates a heatmap from metabolite data for F. cylindrus generated from the TSQ

# Load required packages
library(tidyverse)
library(gplots)
library(tibble)
library(MetBrewer)
library(ggsignif)
library(ggpubr)


# Go to working directory 
setwd("~/Bertrand Lab Dropbox/Bertrand Lab shared workspace/Catalina/Summer_2022/1_FracyPibo_Beyond_Auxotrophy/BA_Manuscript/BA_FragPibo_Quotas_Pipeline")

# Go to cal curve folder 
setwd("./1_TSQ_Frag_Metab_CalCurveQuant")

# Define sample names by order
samples <- c("6_+",
             "24_+",
             "28_+",
             "15_-",
             "16_-",
             "27_-")

# Data cleanup
quant_data <-
  read.csv("metabs_data_quant.csv") %>% # Read in export file
  filter(!grepl("B12", Molecule.Name)) %>% # filter out B12 analogs
  filter(reccomended_norm != "no quant") # filter out no-quant metabs
  

# Clean up data
heatmap_data <-
  quant_data %>% dplyr::select(Replicate.Name, Molecule.Name, B12.Treatment, peak_cell) %>% # Filter for only required cols
  group_by(Molecule.Name, Replicate.Name, B12.Treatment) %>%
  dplyr::summarise(mean_peak_cell = mean(peak_cell,
                                         na.rm = TRUE)) %>% # create mean by bioreplicate from techreps
  filter(Molecule.Name != "HMP") %>% #HMP peak found to be unreliable
  mutate(
    Molecule.Name = replace(
      Molecule.Name,
      Molecule.Name == "DMSP",
      "DMSP*"
    ),
    Molecule.Name = replace(
      Molecule.Name,
      Molecule.Name == "cHET",
      "cHET*")
  ) %>%
  pivot_wider(names_from = Molecule.Name,
              values_from = mean_peak_cell) %>%
  unite("Rep_B12",
        Replicate.Name:B12.Treatment) %>% # Create a column with biorep numbers and B12 treatment5
  mutate(Rep_B12 = factor(Rep_B12,
                          levels = samples)) %>%
  slice(match(samples, Rep_B12)) %>%
  column_to_rownames("Rep_B12") %>% # Change to rowname
  t() %>% # transpose df
  as.matrix() # change to matrix



# use metbrewer to generate color palette
palette <- rev(met.brewer("Hiroshige", n=100))

# Create label vector
labvec <- c(NA, expression("+B"[12]), NA, NA,expression("-B"[12]),  NA)
labvec_rows <- c("DMSP*", "cHET*", expression("B"[7]), "SAM", "SAH", "Proline", "Homarine", "GBT", expression("B"[1])) 

# Go to cal curve folder 
setwd("../3_TSQ_Frag_Metab_Heatmap")

# Open a pdf file
pdf("Frag_heatmap.pdf") 

# Create heatmap
heatmap.2(heatmap_data,
          Colv = FALSE,
          trace = "none", # remove trace
          density = "none", # remove density label
          col = palette, # colors
          scale = "row", # calculate z scores by row
          sepwidth = c(.1, .1), # cell border thickness
          colsep = 3, # space between the 3rd column
          sepcolor = "white",
          labCol = labvec,
          #labRow = labvec_rows,
          cexCol = 2,
          cexRow = 2,
          srtCol = 0,
          adjCol = .5,
          offsetCol = .75,
          margins = c(7, 10))  # cell border color


# Close the pdf file
dev.off() 



cHET_plot_data <- quant_data %>%
  filter(Molecule.Name == "cHET") %>%
  group_by(B12.Treatment, Replicate.Name) %>%
  summarise(mean_peak_cell = mean(peak_cell),
            mean_attomol_cell = mean(fmol_cell) * 1000)  %>%
  mutate(B12.Treatment = factor(B12.Treatment, levels = c("+", "-")))


# Check difference with student's t-test
t.test(mean_peak_cell ~ B12.Treatment, data = cHET_plot_data)
# P = .04923

cHET_plot_peak_cell <- 
  ggplot(data = cHET_plot_data,
       aes(x = B12.Treatment,
           y = mean_peak_cell)) +
  geom_point(stat = "identity",
             size = 10,
             color = "darkgrey") +
  ylab(expression(paste("Mean Peak Cell" ^ "-1" ~ "cHET"))) +
  xlab("Treatment") + 
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
   ylim(0,45) +
  geom_signif(comparisons = list(c("+", "-")),
              map_signif_level = FALSE, 
              y_position = 39, 
              test = "t.test",
              textsize = 10) + 
    scale_x_discrete(labels = c("+" = expression("+B"[12]),
                              "-" = expression("-B"[12])))
    

# attomol plot 
cHET_plot_peak_cell<- ggplot(data = cHET_plot_data,
       aes(x = B12.Treatment,
           y = mean_attomol_cell)) +
  geom_point(stat = "identity",
             size = 10,
             color = "darkgrey") +
  ylab(expression(paste("Attomoles cHET Cell" ^ "-1"))) +
  xlab("Treatment") + 
  theme_classic() +
  theme(
    text = element_text(size = 23, color = "black")
  ) +
   ylim(0,7) +
  geom_signif(comparisons = list(c("+", "-")),
              map_signif_level = FALSE, 
              y_position = 5, 
              test = "t.test",
              textsize = 12) + 
  scale_x_discrete(labels = c("+" = expression("+B"[12]),
                              "-" = expression("-B"[12])))


DMSP_plot_data <- quant_data %>%
  filter(Molecule.Name == "DMSP") %>%
  group_by(B12.Treatment, Replicate.Name) %>%
  summarise(mean_peak_cell = mean(peak_cell))  %>%
  mutate(B12.Treatment = factor(B12.Treatment, levels = c("+", "-")))


# Check difference with student's t-test
t.test(mean_peak_cell ~ B12.Treatment, data = DMSP_plot_data)
# P = 0.03828

DMSP_plot <- ggplot(data = DMSP_plot_data,
       aes(x = B12.Treatment,
           y = mean_peak_cell)) +
  geom_point(stat = "identity",
             size = 10,
             color = "darkgrey") +
  ylab(expression(paste("Mean DMSP Peak Area Cell" ^ "-1"))) +
  xlab("Treatment") + 
  theme_classic() +
  theme(
    text = element_text(size = 23, color = "black")
  ) +
  ylim(0,50) +
  geom_signif(comparisons = list(c("+", "-")),
              map_signif_level = FALSE, 
              y_position = 43, 
              test = "t.test",
              textsize = 12) + 
  scale_x_discrete(labels = c("+" = expression("+B"[12]),
                              "-" = expression("-B"[12])))


# Open a pdf file
pdf("Frag_faceted_metab.pdf", width = 12, height = 6.5) 

ggarrange(DMSP_plot, cHET_plot_peak_cell) 

dev.off()
