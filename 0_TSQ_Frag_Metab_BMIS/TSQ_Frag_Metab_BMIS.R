# This script is for deciding on best matched internal standards (BMIS) in quality control (QC) samples run on TSQ for Global Metabolite F. cylindrus experiment 

# Load required packages into environment
library(ggplot2)
library(dplyr)
library(raster)
library(MetBrewer)
library(gridExtra)
library(tidyr)

# set working directory 
setwd("~/Dropbox (Bertrand Lab)/Bertrand Lab's shared workspace/Catalina/Summer_2022/1_FracyPibo_Beyond_Auxotrophy/BA_Manuscript/BA_FragPibo_Quotas_Pipeline/TSQ_Frag_Metab_RAW")


# load QC data and only include QC treatments
QC.Data <-
  read.csv("ER3_149_Frag_Catalina_01_Cal03022022_output.csv") %>% filter(grepl("QC", Replicate.Name))


# Plot QC peaks per compound
QC.Data$Replicate.Name <-
  factor(
    QC.Data$Replicate.Name,
    levels = c(
      "ER3_149_2fdQC_03",
      "ER3_149_2fdQC_04" ,
      "ER3_149_2fdQC_08" ,
      "ER3_149_2fdQC_09",
      "ER3_149_2fdQC_44" ,
      "ER3_149_2fdQC_79" ,
      "ER3_149_2fdQC_114",
      "ER3_149_2fdQC_142"
    )
  )

# Fix DMB name 
QC.Data <- QC.Data %>% 
  mutate(Molecule.Name = replace(Molecule.Name, Molecule.Name == "Dimethyl-benzimidazole (DMB)", "DMB"))

QC.Peaks <- ggplot() + 
  geom_bar(data = QC.Data, aes(y = Total.Area, x = Replicate.Name), stat = "identity") + 
  ggtitle("QC Peak Areas (No Normalization)") + 
  ylab("Peak Area") +
  xlab("QC Injection") +
  facet_wrap( ~ Molecule.Name, scales = "free") +
  theme(axis.text.x = element_text(angle = 90,   hjust = 1)) + 
  scale_x_discrete(breaks = levels(QC.Data$Replicate.Name), 
                   labels=c("03", "04", "08", "09", "44", "79", "114", "142" ))
  

# Get mean, standard deviation, and cv of all QC's per compound
QC.Data.SumStats <- QC.Data %>%
  group_by(Molecule.Name) %>%
  dplyr::summarise(Mean.Peak.Area = mean(Total.Area), SD.Peak.Area = sd(Total.Area), CV.Peak.Area = raster::cv(Total.Area))



# Get Normalized Peaks (divide each light QC by each heavy QC) ------------

# Create df with only non-heavy entries
BMIS.Data <- QC.Data %>% filter(!grepl('heavy', Molecule.Name))

# Create a list of heavy compounds to normalize to
heavy_compounds <- c("B1", "B2", "B12-CN", "B7")

# Create a function that does the normalization
norm_peaks <- function(compound_name) {
  
  # Index for heavy compound from QC data
  heavy <- QC.Data %>% filter(Molecule.Name == paste(compound_name, "-heavy", sep = ""))
  
  # Amend areas of the heavy compound to the CalCurve_cals df
  BMIS.Data[paste("Heavy", compound_name, "Peaks", sep = ".")] <-
    heavy$Total.Area
  
  # Calculate normalized peaks by dividing each compound by the heavy peak
  BMIS.Data[paste("Heavy", compound_name, "Norm", sep = ".")] <-
    BMIS.Data$Total.Area / BMIS.Data[paste("Heavy", compound_name, "Peaks", sep = ".")]
 
  # save new df with normalization to global environment 
  assign("BMIS.Data", BMIS.Data, envir = globalenv()) 
  
}


# Apply with list of heavy compounds
lapply(heavy_compounds, norm_peaks)


# Compare Heavy-Normalized QC to non-normed QC ----------------------------

# Get summary stats for heavy compound normalization (sd and cv)
HeavyNormQC.Data.SumStats <- BMIS.Data %>%
  group_by(Molecule.Name) %>%
  dplyr::summarise(
    Mean.Peak.Area = mean(Total.Area),
    SD.Peak.Area = sd(Total.Area),
    CV.Peak.Area = raster::cv(Total.Area),
    Mean.Peak.Area.B1.Norm = mean(Heavy.B1.Norm),
    SD.Peak.Area.B1.Norm = sd(Heavy.B1.Norm),
    CV.Peak.Area.B1.Norm = raster::cv(Heavy.B1.Norm),
    Mean.Peak.Area.B2.Norm = mean(Heavy.B2.Norm),
    SD.Peak.Area.B2.Norm = sd(Heavy.B2.Norm),
    CV.Peak.Area.B2.Norm = raster::cv(Heavy.B2.Norm),
    `Mean.Peak.Area.B12-CN.Norm` = mean(`Heavy.B12-CN.Norm`),
    `SD.Peak.Area.B12-CN.Norm` = sd(`Heavy.B12-CN.Norm`),
    `CV.Peak.Area.B12-CN.Norm` = raster::cv(`Heavy.B12-CN.Norm`),
    Mean.Peak.Area.B7.Norm = mean(Heavy.B7.Norm),
    SD.Peak.Area.B7.Norm = sd(Heavy.B7.Norm),
    CV.Peak.Area.B7.Norm = raster::cv(Heavy.B7.Norm)
  )


# Compare CV's per compound and normalization type -------------------------------
HeavyNormQC.Data.SumStats1 <-
  HeavyNormQC.Data.SumStats %>% dplyr::select(
    Molecule.Name,
    CV.Peak.Area,
    CV.Peak.Area.B1.Norm,
    CV.Peak.Area.B2.Norm,
    `CV.Peak.Area.B12-CN.Norm`,
    CV.Peak.Area.B7.Norm
  )

colnames(HeavyNormQC.Data.SumStats1)[2:6] <-
  c("No Normalization",
    "Heavy B1",
    "Heavy B2",
    "Heavy B12-CN",
    "Heavy B7")


gather_cols <-
  c("No Normalization",
    "Heavy B1",
    "Heavy B2",
    "Heavy B12-CN",
    "Heavy B7")


# change to long format for plotting
HeavyNormQC.Data.Comp <-
  gather(HeavyNormQC.Data.SumStats1,
         Norm.Compound,
         Norm.Peak.CV,
         gather_cols,
         factor_key = TRUE)

ggplot() + 
  geom_bar(data = HeavyNormQC.Data.Comp, aes(y = Norm.Peak.CV, x = Norm.Compound, fill = Norm.Compound), stat = "identity") + 
  ylab("CV") +
  xlab("Normalization Compound") +
  facet_wrap( ~ Molecule.Name, scales = "fixed") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15)) +
  scale_fill_manual(values=met.brewer("Cross", 5), name = "Normalization Compound")



# Choose BMIS For Each Light Compound -------------------------------------

# Create a function that calculates the differences in cv before and after normalization to each heavy compound
calc_deltas <- function(compound_name) {
  
  # Subtract original compound cv from new with normalization
  HeavyNormQC.Data.SumStats[paste(compound_name)] <- HeavyNormQC.Data.SumStats$CV.Peak.Area - HeavyNormQC.Data.SumStats[paste("CV.Peak.Area", compound_name, "Norm", sep = ".")]

  assign("HeavyNormQC.Data.SumStats", HeavyNormQC.Data.SumStats, envir = globalenv())
}


# Apply with list of heavy compounds
lapply(heavy_compounds, calc_deltas)

# create a list of highest changes and which heavy standard norm led to them
HeavyNormQC.Data.SumStats.deltaCV <-
  HeavyNormQC.Data.SumStats %>%
  dplyr::select(B1, B2, `B12-CN`, B7)
  
# Pull largest cv's
BMIS_list <-
  data.frame(colnames(HeavyNormQC.Data.SumStats.deltaCV)[apply(HeavyNormQC.Data.SumStats.deltaCV, 1, which.max)]) 

# Update column name 
colnames(BMIS_list) <- "BMIS"

# Create deltaCV column with max
BMIS_list$deltaCV <-
  apply(HeavyNormQC.Data.SumStats.deltaCV, 1, max)

# Add back to dataframe
HeavyNormQC.Data.BMIS <- cbind(HeavyNormQC.Data.SumStats, BMIS_list)

# Change column namE
HeavyNormQC.Data.BMIS$BMIS_col_name <- paste("Mean.Peak.Area.",HeavyNormQC.Data.BMIS$BMIS,".Norm", sep = "")




# Make loop that sees if max change is larger than 30%
# If so, make final norm column show the normalization that caused improvement 
# Start indexing vector at 7 to ignore B1, the B12's, and B2 for BMIS selection (they have corresponding internal standards). B7 not included because B7:heavy B7 norm increases cv rather than reduces

# Make df with compounds with corresponding internal standards
HeavyNormQC.Data.BMIS.corresp <- HeavyNormQC.Data.BMIS[1:6,]
for (i in 1:nrow(HeavyNormQC.Data.BMIS.corresp)) {
  HeavyNormQC.Data.BMIS.corresp$Final.BMIS.Norm.Peak[i] <-
    HeavyNormQC.Data.BMIS.corresp[i, as.character(HeavyNormQC.Data.BMIS.corresp$BMIS_col_name[i])]
  HeavyNormQC.Data.BMIS.corresp$BMIS_used[[i]] <-
    HeavyNormQC.Data.BMIS.corresp$BMIS[i]
}


# Make df for other compounds 
HeavyNormQC.Data.BMIS.NOcorresp <- HeavyNormQC.Data.BMIS[7:26,]

# Do loop to choose BMIS
for (i in 1:nrow(HeavyNormQC.Data.BMIS.NOcorresp)){
  if (HeavyNormQC.Data.BMIS.NOcorresp$deltaCV[i] < 30){
    HeavyNormQC.Data.BMIS.NOcorresp$Final.BMIS.Norm.Peak[i] <- HeavyNormQC.Data.BMIS.NOcorresp$Mean.Peak.Area[i]
    HeavyNormQC.Data.BMIS.NOcorresp$BMIS_used[i] <- "none"
  } else {
    HeavyNormQC.Data.BMIS.NOcorresp$Final.BMIS.Norm.Peak[i] <- HeavyNormQC.Data.BMIS.NOcorresp[i,as.character(HeavyNormQC.Data.BMIS.NOcorresp$BMIS_col_name[i])]
    HeavyNormQC.Data.BMIS.NOcorresp$BMIS_used[[i]] <- HeavyNormQC.Data.BMIS.NOcorresp$BMIS[i]
  }
}

HeavyNormQC.Data.BMIS <- rbind(HeavyNormQC.Data.BMIS.corresp, HeavyNormQC.Data.BMIS.NOcorresp)

# Create dataframe for exporting df's
QC_Norm_Export <-
  HeavyNormQC.Data.BMIS %>% dplyr::select(
    Molecule.Name,
    Mean.Peak.Area,
    SD.Peak.Area,
    CV.Peak.Area,
    Mean.Peak.Area.B1.Norm,
    SD.Peak.Area.B1.Norm,
    CV.Peak.Area.B1.Norm,
    Mean.Peak.Area.B2.Norm,
    SD.Peak.Area.B2.Norm,
    CV.Peak.Area.B2.Norm,
    `Mean.Peak.Area.B12-CN.Norm`,
    `SD.Peak.Area.B12-CN.Norm`,
    `CV.Peak.Area.B12-CN.Norm`,
    Mean.Peak.Area.B7.Norm,
    SD.Peak.Area.B7.Norm,
    CV.Peak.Area.B7.Norm,
    BMIS,
    deltaCV,
    BMIS_used,
    Final.BMIS.Norm.Peak
  )

QC_Norm_Export_Sum <-
  QC_Norm_Export %>% dplyr::select(Molecule.Name, BMIS, deltaCV, BMIS_used, Final.BMIS.Norm.Peak)

# Fix weird encoding error
QC_Norm_Export <- apply(QC_Norm_Export,2,as.character)
QC_Norm_Export_Sum <- apply(QC_Norm_Export_Sum,2,as.character)


# Change wd
setwd("../0_TSQ_Frag_Metab_BMIS")

# Export results and summary from BMIS analysis
write.csv(QC_Norm_Export, file = "QC_BMIS_results.csv")
write.csv(QC_Norm_Export_Sum, file = "QC_BMIS_results_sum.csv")

# Create a table from export
grid.table(QC_Norm_Export_Sum)

# Plot QC's before and after normalization for those that get a normalization

# Get before and after CV's
QC_norm_comp_df <- HeavyNormQC.Data.BMIS %>% filter(BMIS_used != "none") %>% dplyr::select(Molecule.Name, CV.Peak.Area, deltaCV) 

# get columns for before and after normalization
QC_norm_comp_df$After_norm <- QC_norm_comp_df$CV.Peak.Area - QC_norm_comp_df$deltaCV

# change to long form data with pivot (note: use of amend gather above)
colnames(QC_norm_comp_df)[c(2,4)] <- c("Before Normalization", "After Normalization") 

QC_norm_comp_df <- QC_norm_comp_df %>% dplyr::select(Molecule.Name, "Before Normalization", "After Normalization")

HeavyNormQC.Data.Comp <- QC_norm_comp_df %>%
  pivot_longer(!Molecule.Name, names_to = "Norm", values_to = "CV")

# filter out methinonine (cv too high)
HeavyNormQC.Data.Comp<- HeavyNormQC.Data.Comp %>% dplyr::filter(Molecule.Name != "Methionine")


HeavyNormQC.Data.Comp$Norm <- factor(HeavyNormQC.Data.Comp$Norm, levels = c("Before Normalization", "After Normalization") )

# Open a pdf file
pdf("Frag_BMIS.pdf") 

# make plot comparing effects of BMIS normalization
ggplot(HeavyNormQC.Data.Comp, aes(fill=Norm, y=CV, x=Molecule.Name)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("Molecule Name") +
  ylab("Quality Control CV") +
  scale_fill_manual(values=met.brewer("Cross", 2), name = "Normalization")  + 
  theme_classic() +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90))

# Close the pdf file
dev.off() 





  
  
