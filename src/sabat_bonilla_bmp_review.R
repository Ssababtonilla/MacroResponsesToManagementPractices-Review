###Analysis script for: Meta-analysis of Stream Macroinvertebrate Responses 
#to Management Practices in Selected Regions of the Chesapeake Bay Watershed 
#Vary with Region, Land Use, Practice Type and Assessment Metrics. 

# -------------------------------------------------------------------------------------
# 0. Install & load packages
# -------------------------------------------------------------------------------------
required_pkgs <- c(
  "ggplot2","dplyr","tidyr","readr","scales","sf","maps","cowplot","knitr","kableExtra",
  "ggpubr","MuMIn","lme4","lmerTest","emmeans","ggh4x","ggspatial","stringr","forcats"
)
new_pkgs <- setdiff(required_pkgs, installed.packages()[,"Package"])

if(length(new_pkgs)) install.packages(new_pkgs)
lapply(required_pkgs, library, character.only=TRUE)

# -------------------------------------------------------------------------------------
# 1. Read & preliminarily clean raw data
# -------------------------------------------------------------------------------------
Data <- read_csv("Regional BMP Database - Data.csv", show_col_types = FALSE)

# Standardize column names
names(Data) <- make.names(names(Data), unique = TRUE)
names(Data)[names(Data) == "BMP.Group..NRCS."] <- "BMP.Group"
names(Data)[names(Data) == "Effect"]              <- "Effect.Size"

# Fix known typos / unify labels
Data <- Data %>%
  # 1a) correct a handful of BMP.Group & Measurement.Type typos
  mutate(
    BMP.Group = recode(
      BMP.Group,
      "Stream Habitat Improvement and management" = "Stream Habitat Improvement and Management",
      .default = BMP.Group
    ),
    Measurement.Type = recode(
      Measurement.Type,
      "Isotopes" = "Isotopic Measurements",
      "P/B"      = "Production/Biomass Ratio",
      .default   = Measurement.Type
    ),
    # 1b) turn blank Effect.Size into NA
  #  Effect.Size = na_if(Effect.Size, ""),
    # 1c) collapse stray blanks in Direction.of.Effect
    Direction.of.Effect = if_else(
      Direction.of.Effect %in% c("", "Positive "),
      "Positive",
      Direction.of.Effect
    ),
    # 1d) turn empty Threat.Type into NA
    Threat.Type = na_if(Threat.Type, ""),
    # 1e) unify trailing‐space variant of "Agriculture "
    Threat.Type = if_else(Threat.Type == "Agriculture ", "Agriculture", Threat.Type),
    # 1f) unify Response.Metric EPT label
    Response.Metric = recode(
      Response.Metric,
      "% EPT " = "% EPT",
      "%EPT"   = "% EPT",
      .default = Response.Metric
    ),
    # 1g) standardize ecoregion names
    Ecoregion = recode(
      Ecoregion,
      "Ridge & Valley"             = "Ridge and Valley",
      "Peidmont"                   = "Piedmont",
      "Nothern Allegheny Plateau"  = "Northern Allegheny Plateau",
      .default                      = Ecoregion
    )
  ) %>%
  # 1h) unify land‐use codes
  mutate(
    Threat.Type = case_when(
      Threat.Type == "(A) Agriculture"  ~ "Agriculture",
      Threat.Type == "(U) Urbanization" ~ "Urbanization",
      Threat.Type == "Mixed"            ~ "Agriculture/Urban",
      TRUE                               ~ Threat.Type
    ),
    # 1i) create BMP_Acronym and drop any that failed to match
    BMP_Acronym = case_when(
      BMP.Group == "Access Control"                           ~ "AC",
      BMP.Group == "Stream Habitat Improvement and Management"~ "SHIM",
      BMP.Group == "Riparian Forest Buffer"                   ~ "RFB",
      BMP.Group == "Streambank and Shoreline Protection"      ~ "SSP",
      BMP.Group == "Structures for Wildlife"                  ~ "SW",
      BMP.Group == "Riparian Herbaceous Cover"                ~ "RHC",
      BMP.Group == "Open Channel"                             ~ "OC",
      BMP.Group == "Waste Treatment"                          ~ "WT",
      BMP.Group == "Stormwater Runoff Control"                ~ "SRC",
      TRUE                                                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(BMP_Acronym))

# -------------------------------------------------------------------------------------
# 1.1 Map physiographic codes → names
# -------------------------------------------------------------------------------------
province_lookup <- tibble(
  Province_Code = c("CP","AP","RV","PE","BR","NE"),
  PROVINCE      = c(
    "Coastal Plain","Appalachian Plateaus",
    "Valley and Ridge","Piedmont","Blue Ridge","New England"
  )
)
Data <- Data %>%
  mutate(Physiographic.Province = toupper(Physiographic.Province)) %>%
  left_join(province_lookup,
            by = c("Physiographic.Province" = "Province_Code")) %>%
  filter(!is.na(PROVINCE)) %>%
  mutate(PROVINCE = toupper(PROVINCE))

# -------------------------------------------------------------------------------------
# 1.2 Convert sample sizes & SDs to numeric
# -------------------------------------------------------------------------------------
Data <- Data %>%
  mutate(across(
    c(
      Mean.Values.for.BMP, Mean.values.for.non.BMP,
      Standard.Deviation, Reference.SD,
      Total.Samples.treatment, Total.Sites
    ),
    as.numeric
  )) %>%
  mutate(
    n_trt = as.integer(Total.Samples.treatment),
    n_ctl = as.integer(Total.Sites)
  )


# -------------------------------
# 2. Compute effect sizes: log-response ratio (L) & Hedges' g
# -------------------------------
Data <- Data %>%
  mutate(
    # log-response ratio
    L = if_else(
      !is.na(Mean.Values.for.BMP) & !is.na(Mean.values.for.non.BMP) &
        Mean.values.for.non.BMP != 0,
      log(Mean.Values.for.BMP / Mean.values.for.non.BMP),
      NA_real_
    ),
    # pooled standard deviation
    Sp = if_else(
      !is.na(Standard.Deviation) & !is.na(Reference.SD) &
        !is.na(n_trt) & !is.na(n_ctl) & (n_trt + n_ctl - 2) > 0,
      sqrt(((n_trt - 1) * Standard.Deviation^2 +
              (n_ctl - 1) * Reference.SD^2) /
             (n_trt + n_ctl - 2)),
      NA_real_
    ),
    # Hedges' d
    Hedges_d = if_else(
      !is.na(Mean.Values.for.BMP) & !is.na(Mean.values.for.non.BMP) &
        !is.na(Sp) & Sp != 0,
      (Mean.Values.for.BMP - Mean.values.for.non.BMP) / Sp,
      NA_real_
    ),
    # small-sample correction J
    J = if_else(
      !is.na(n_trt) & !is.na(n_ctl) & (n_trt + n_ctl - 2) > 0,
      1 - (3 / (4 * (n_trt + n_ctl - 2) - 1)),
      NA_real_
    ),
    # Hedges' g
    Hedges_g = Hedges_d * J
  )

# -------------------------------
# 3. Define mapping: Response.Metric → Metric.Type & Metric.Group
# -------------------------------

# (a) Full structural metrics list
structural_metrics <- c(
  # Abundance Metrics (35)
  "Total Abundance", "Average Abundance Across Reps",
  "Chironomidae Abundance", "Baetis Abundance", "Stenonema Abundance",
  "Paraleptophlebia Abundance", "Ephemerella Abundance", "Ameletus Abundance",
  "Hexatoma Abundance", "Stenacron Abundance", "Caenis Abundance",
  "Neoperla Abundance", "Cheumatopsyche Abundance", "Chimarra Abundance",
  "Stenelmis Abundance", "Psephenus Abundance", "Clinger Abundance",
  "Ept Abundance", "Diptera Abundance", "Plecoptera Abundance",
  "Coleoptera Abundance", "Ephemeroptera Abundance", "Trichoptera Abundance",
  "Trichoptera", "Oligochaeta", "Decapoda", "Collector Abundance",
  "Filterer Abundance", "Scraper Abundance", "Shredder Abundance",
  "Predator Abundance", "Burrower Abundance", "Climber Abundance",
  "Swimmer Abundance", "Sprawler Abundance", "# Egg Masses In A Stream",
  
  # Additional Abundance (6)
  "Cote Abundance", "Cch Abundance", "Odonata Abundance",
  "Gastropoda Abundance", "Abundance Of Dominant Taxa",
  "Chironomidae/Oligochaeta Abundance",
  
  # Relative & Mean Relative Abundance (15)
  "Emergence Synchronization (Poor) Relative Abundance",
  "Emergence Synchronization (Well) Relative Abundance",
  "Occurrence In Drift (Rare) Relative Abundance",
  "Occurrence In Drift (Abundant) Relative Abundance",
  "Swimming Ability (Strong) Relative Abundance",
  "Swimming Ability (Weak) Relative Abundance",
  "Swimmer Relative Abundance", "Collector Relative Abundance",
  "Herbivore Relative Abundance", "Plecoptera Mean Relative Abundance",
  "Amphipoda Mean Relative Abundance", "Isopoda Mean Relative Abundance",
  "Coleoptera Mean Relative Abundance", "Diptera Mean Relative Abundance",
  "Ephemeroptera Mean Relative Abundance", "Trichoptera Mean Relative Abundance",
  "Oligochaeta Mean Relative Abundance", "Turbellaria Mean Relative Abundance",
  "Total Mean Relative Abundance",
  
  # Density Metrics (2)
  "Density", "Total Density",
  
  # Richness Metrics (20)
  "Family Richness", "Genus Richness", "Species Richness", "EPT Richness",
  "Collector Richness", "Filterer Richness", "Scraper Richness",
  "Shredder Richness", "Predator Richness", "Ephemeroptera Richness",
  "Plecoptera Richness", "Trichoptera Richness", "Total Richness",
  "# Indicator Taxa", "Average Family Richness Across Reps",
  "# Dominant Genera", "Cote Richness", "Cch Richness", "Odonata Richness",
  "Gastropoda Richness",
  
  # Diversity Metrics (6)
  "Shannon", "Simpson", "Jaccard", "Total Diversity", "Total Evenness",
  "Shannon Functional Trait-State Diversity",
  
  # Composition Metrics (remaining)
  "% EPT", "% E", "% P", "% T", "% Diptera", "% Non-Insects",
  "% Dominant Taxa", "% Dominant Families", "% Dominant Genera",
  "% Top 2 Dominant Taxa", "% 2 Dominant Taxa", "% Chironomidae",
  "% Chironomidae And Oligochaeta", "% Oligochaeta", "% Coleoptera",
  "% Elmidae", "% Crustacea", "% Isopoda", "% Insects",
  "% EPT-Hydropsychidae", "% EPT/Chironomidae", "% EPT-H",
  "% EPT-Cheumatopsyche Abundance", "% Clingers", "% Scrapers",
  "% D", "% Cheumatopsyche", "% Hydropsychidae", "% Stenelmis",
  "% Baetis", "% Calopteryx", "% Tipula", "% Chimarra",
  "% Leuctra", "% Strophopteryx", "% Prosimulium",
  "% Paraleptophlebia", "% Acroneuria", "% Epeorus", "% Dominance",
  "% Captures Chironomidae", "% Captures Baetidae",
  "% Captures Hydropsychidae", "% Captures Annelida",
  "% Captures Heptageniidae", "% Captures Simuliidae",
  "% Collectors", "% Filterers", "% Shredders", "% Predators",
  "% Gatherers", "# FFG",
  "% Captures Of FFG'S Composed Of Chironomidae (Filterers)",
  "% Captures Of FFG'S Composed Of Chironomidae (Collectors)",
  "% Captures Of FFG'S Composed Of Chironomidae (Shredders)",
  "% Captures Of FFG'S Composed Of Chironomidae (Scrapers)",
  "% Captures Of FFG'S Composed Of Chironomidae (Predators)",
  "% Captures Of FFG'S Composed Of Chironomidae (Parasites)",
  "Ratio SC/CF", "Ratio SH/Total Number", "Ratio EPT/Chironomidae",
  "Egg Mass Distribution", "Univoltine Relative Abundance",
  "Multivoltine Relative Abundance", "Emergence Synchronization",
  "Swimming Ability Relative Abundance", "Clinger Relative Abundance",
  "% 5 Dominant Taxa", "% Model Affinity", "Dominant 3 Taxa"
)

# (b) Full functional metrics list
functional_metrics <- c(
  # Biomass (14)
  "Prey Biomass", "Predator Biomass", "Shredder Biomass",
  "Scraper Biomass", "Collector Biomass", "Filterer Biomass",
  "Total Biomass", "Non-Insect Biomass", "Ephemeroptera Biomass",
  "Plecoptera Biomass", "Trichoptera Biomass", "Diptera Biomass",
  "Odonata Biomass", "Coleoptera Biomass",
  
  # Secondary Production (12)
  "Scraper Secondary Production", "Shredder Secondary Production",
  "Collector Secondary Production", "Filterer Secondary Production",
  "Predator Secondary Production", "Total Secondary Production",
  "Scraper P/B", "Shredder P/B", "Collector P/B", "Filterer P/B",
  "Predator P/B", "Total P/B",
  
  # Decomposition (1)
  "Decomposition Rate",
  
  # Isotopic (2)
  "Isotopic N", "Isotopic C",
  
  # Habitat Preferences (4)
  "% Captures Of Habitat Type (Leaf Packs)",
  "% Captures Of Habitat Type (Sediment)",
  "% Captures Of Habitat Type (Wood)",
  "Temperature Preference Metric"
)

# (c) Full tolerance metrics list
tolerance_metrics <- c(
  "Average Tolerance", "% Tolerant", "% Intolerant",
  "# Intolerant Families", "% TN Nutrient Tolerant Organisms",
  "% T-Hydropsychidae", "HBI"
)

# (d) Full biotic indices list
biotic_indices <- c(
  "HBI", "IBI", "BIBI", "EPT Index", "SCI", "VASCI", "KIBI",
  "NCBI", "FSBI", "EPTBI", "FBI", "Bioclassification", "Average IBI",
  "Fine Sediment Biotic Index"
)

# (e) Combine all response metrics
response_metrics <- c(
  structural_metrics,
  functional_metrics,
  tolerance_metrics,
  biotic_indices
)

# (f) Define Metric.Type categories (using original counts)
metric_types <- c(
  # Structural portion
  rep("Abundance", 35),    # first 35 abundance
  rep("Abundance", 6),     # next 6 abundance
  rep("Abundance", 15),    # 15 relative/mean relative
  rep("Density", 2),       # 2 density
  rep("Richness", 20),     # 20 richness
  rep("Diversity", 6),     # 6 diversity (including Simpson)
  rep("Composition Metrics", 
      length(structural_metrics) - (35 + 6 + 15 + 2 + 20 + 6)
  ),
  # Functional portion
  rep("Biomass", 14),           # 14 biomass
  rep("Secondary Production", 12), # 12 secondary
  "Decomposition",             # 1 decomposition
  rep("Isotopic Measurements", 2), # 2 isotopic
  rep("Habitat Preferences", 4),  # 4 habitat
  # Tolerance portion
  rep("Tolerance Metrics", length(tolerance_metrics)),
  # Biotic Indices portion
  rep("Biotic Indices", length(biotic_indices))
)

# (g) Define Metric.Group categories
metric_groups <- c(
  rep("Structural", length(structural_metrics)),
  rep("Functional", length(functional_metrics)),
  rep("Tolerance Metrics", length(tolerance_metrics)),
  rep("Biotic Indices", length(biotic_indices))
)

# (h) Build lookup table (uppercase) and dedupe
metric_mapping <- tibble(
  R.M          = toupper(response_metrics),
  Metric.Type  = metric_types,
  Metric.Group = metric_groups
) %>%
  distinct(R.M, .keep_all = TRUE)

# -------------------------------
# 4. Merge metric mapping into Data & fix variants
# -------------------------------
Data <- Data %>%
  mutate(
    Response.Metric = trimws(Response.Metric),
    Response.Metric = str_replace_all(Response.Metric, c(
      "Oligochoeata Mean Relative Abundance" = "Oligochaeta Mean Relative Abundance",
      "Hilsenhoff"                          = "HBI",
      "Sychronization"                     = "Synchronization",
      "Occurence"                          = "Occurrence",
      "^%Scrapers"                         = "% Scrapers",
      "Captures Heptagenaidae"             = "Captures Heptageniidae",
      "% Caleopteryx"                      = "% Calopteryx",
      "% Hydropsyche"                      = "% Hydropsychidae"
    )),
    R.M = toupper(Response.Metric)
  ) %>%
  left_join(metric_mapping, by = "R.M") %>%
  select(-R.M) %>%
  # force Simpson → Diversity
  mutate(
    Metric.Type  = if_else(Response.Metric == "Simpson", "Diversity", Metric.Type),
    Metric.Group = if_else(Response.Metric == "Simpson", "Structural", Metric.Group)
  ) %>%
  distinct()

# -------------------------------
# 4c. Sanity check – no missing Metric.Type & Simpson is Diversity
# -------------------------------
unmatched <- Data %>%
  filter(
    is.na(Metric.Type) |
      (Response.Metric == "Simpson" & Metric.Type != "Diversity")
  ) %>%
  distinct(Response.Metric, Metric.Type)

if (nrow(unmatched) > 0) {
  message("⚠️  Issue with metric mapping:\n",
          paste0(unmatched$Response.Metric, " -> ", unmatched$Metric.Type, collapse = "\n"))
} else {
  message("✅  All metrics mapped correctly, and Simpson is Diversity.")
}

# -------------------------------------------------------------------------------------
# 5. Build single cleaned df
# -------------------------------------------------------------------------------------
bmps_with_data <- Data %>%
  group_by(BMP_Acronym) %>%
  summarise(
    has_ES = any(!is.na(Effect.Size)),
    has_g  = any(!is.na(Hedges_g)),
    .groups="drop"
  ) %>%
  filter(has_ES|has_g) %>%
  pull(BMP_Acronym)

df <- Data %>%
  filter(BMP_Acronym %in% bmps_with_data) %>%
  mutate(
    Measurement.Type = case_when(
      Metric.Type=="Biomass" ~ "Biomass",
      Measurement.Type=="Production/Biomass Ratio" ~ "P/B Ratio",
      Measurement.Type=="Relative Abundance"      ~ "Abundance",
      Measurement.Type=="Decomposition"           ~ "Decomposition",
      Measurement.Type=="Isotopic Measurements"   ~ "Isotopic",
      TRUE                                         ~ Measurement.Type
    ),
    Direction = case_when(
      Effect.Size >  0.2               ~ "Positive",
      Effect.Size < -0.2               ~ "Negative",
      between(Effect.Size,-0.2,0.2)    ~ "No Response",
      TRUE                              ~ NA_character_
    )
  ) %>%
  filter(!is.na(Effect.Size), !is.nan(Effect.Size))

total_studies <- n_distinct(df$RefID)


# -------------------------------------------------------------------------------------
# 6. Q1: Most common BMPs
# -------------------------------------------------------------------------------------
most_common_bmps <- df %>%
  group_by(BMP_Acronym,BMP.Group) %>%
  summarise(Studies=n_distinct(RefID),.groups="drop") %>%
  arrange(desc(Studies))

kable(
  most_common_bmps,
  caption = "Most Common BMPs"
) %>%
  kable_styling(full_width = FALSE)
write_csv(most_common_bmps,"Most_Common_BMPs.csv")


# -------------------------------------------------------------------------------------
# 7. Q2: Macroinvertebrate responses by BMP
# -------------------------------------------------------------------------------------
resp_summary <- df %>%
  mutate(Dir.of.Effect = case_when(
    Hedges_g >= 0.8              ~ "Large Positive",
    Hedges_g >= 0.5              ~ "Medium Positive",
    Hedges_g >= 0.2              ~ "Small Positive",
    Hedges_g <= -0.8             ~ "Large Negative",
    Hedges_g <= -0.5             ~ "Medium Negative",
    Hedges_g <= -0.2             ~ "Small Negative",
    TRUE                          ~ "No Response"
  )) %>%
  group_by(BMP_Acronym,Dir.of.Effect) %>%
  summarise(
    Mean_g = mean(Hedges_g,na.rm=TRUE),
    Mean_ES= mean(Effect.Size,na.rm=TRUE),
    Studies= n_distinct(RefID),
    .groups="drop"
  )

kable(
  resp_summary,
  caption ="Responses by BMP"
) %>% 
  kable_styling(full_width=FALSE)
write_csv(resp_summary,"Responses_by_BMP.csv")

# boxplot
custom_cols <- c(
  AC="#1f78b4", SHIM="#e31a1c", RFB="#33a02c", SSP="#ff7f00",
  SW="#6a3d9a", RHC="#b2df8a", OC="#a6cee3", WT="#b15928", SRC="#fb9a99"
)

ggplot(df, aes(BMP_Acronym, Effect.Size, fill = BMP_Acronym)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.3, color = "black") +
  stat_summary(fun = mean, geom = "point", color = "blue", shape = 18, size = 3) +
  scale_fill_manual(values = custom_cols[names(custom_cols) %in% bmps_with_data]) +
  labs(
    x        = "Management Practice",
    y        = "Effect Size",
    fill     = "Practice"
  ) +
  theme_minimal() +
  theme(
    axis.title.x     = element_text(size = 22, margin = margin(t = 10)),
    axis.title.y     = element_text(size = 22, margin = margin(r = 10)),
    axis.text.x      = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 18),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin      = margin(15, 15, 15, 15)
  ) +
  geom_hline(yintercept = c(-0.2, 0.2),
             linetype = "dashed",
             color    = "red")


# -------------------------------------------------------------------------------------
# 8. Q3: Responses by Land‐Use × BMP
# -------------------------------------------------------------------------------------
df_plot <- df %>%
  mutate(
    Threat.Type = fct_relevel(
      if_else(Threat.Type=="A/U","Mixed",Threat.Type),
      "Agriculture","Urbanization","Mixed"
    ),
    combo = paste(Threat.Type,BMP_Acronym,sep="_")
  )

df_x <- df_plot %>%
  distinct(Threat.Type,BMP_Acronym,combo) %>%
  arrange(Threat.Type,BMP_Acronym) %>%
  mutate(x_pos = row_number(), x_label = BMP_Acronym)

df_plot <- df_plot %>%
  left_join(df_x,by="combo") %>%
  mutate(combo_factor = factor(combo,levels=df_x$combo))

group_labels <- df_x %>%
  group_by(Threat.Type) %>%
  summarise(min_x=min(x_pos),max_x=max(x_pos),
            mid_x=(min_x+max_x)/2,.groups="drop")

ggplot(df_plot, aes(combo_factor, Effect.Size)) +
  geom_boxplot(outlier.shape = NA, fill = "gray80", alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, color = "black") +
  stat_summary(fun = mean, geom = "point", color = "blue", shape = 18, size = 3) +
  scale_x_discrete(labels = df_x$x_label) +
  labs(
    x        = "Management Practice by Land Use",
    y        = "Effect Size"
  ) +
  theme_minimal() +
  theme(
    axis.title.x     = element_text(size = 24, margin = margin(t = 10)),
    axis.title.y     = element_text(size = 24, margin = margin(r = 10)),
    axis.text.x      = element_text(size = 20),
    axis.text.y      = element_text(size = 20),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin      = margin(15, 15, 15, 15),
    legend.position  = "none"
  ) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "red") +
  coord_cartesian(ylim = c(NA, max(df_plot$Effect.Size) + 1), clip = "off") +
  annotate(
    "text",
    x = group_labels$mid_x,
    y = max(df_plot$Effect.Size) + 0.7,
    label = group_labels$Threat.Type,
    size = 6,
    fontface = "bold"
  )

# -------------------------------------------------------------------------------------
# 9. Q4: Proportional BMP stacks by Province
# -------------------------------------------------------------------------------------
prov_chart <- df %>%
  filter(!is.na(PROVINCE)&PROVINCE!="") %>%
  group_by(PROVINCE,BMP_Acronym) %>%
  summarise(Studies=n_distinct(RefID),.groups="drop") %>%
  group_by(PROVINCE) %>%
  mutate(Proportion=Studies/sum(Studies)) %>%
  ungroup()

ggplot(prov_chart,aes(PROVINCE,Proportion,fill=BMP_Acronym))+
  geom_col(color="black")+
  scale_y_continuous(labels=percent_format(1))+
  labs(x="Province",y="Proportion of Studies",fill="BMP")+
  theme_minimal()+theme(
    axis.text.x=element_text(size=14,angle=45,hjust=1),
    axis.text.y=element_text(size=14),
    axis.title=element_text(size=16),
    panel.border=element_rect(fill=NA,colour="black")
  )+
  scale_fill_manual(values=custom_cols)

# -------------------------------------------------------------------------------------
# 10. Q5: Proportional BMP stacks by Land‐Use
# -------------------------------------------------------------------------------------
land_chart <- df %>%
  group_by(Threat.Type,BMP_Acronym) %>%
  summarise(Studies=n_distinct(RefID),.groups="drop") %>%
  group_by(Threat.Type) %>%
  mutate(Proportion=Studies/sum(Studies)) %>%
  ungroup()

ggplot(land_chart,aes(Threat.Type,Proportion,fill=BMP_Acronym))+
  geom_col(color="black")+
  scale_y_continuous(labels=percent_format(1))+
  labs(x="Land Use",y="Proportion of Studies",fill="BMP")+
  theme_minimal()+theme(
    axis.text=element_text(size=14),
    axis.title=element_text(size=16),
    panel.border=element_rect(fill=NA,colour="black")
  )+
  scale_fill_manual(values=custom_cols)


# -------------------------------------------------------------------------------------
# 11. Q6: Metric.Type summary & boxplot
# -------------------------------------------------------------------------------------
metric_summary <- df %>%
  group_by(Metric.Type) %>%
  summarise(
    Studies=n_distinct(RefID),
    Mean_ES=mean(Effect.Size,na.rm=TRUE),
    Median_ES=median(Effect.Size,na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(Percent=Studies/total_studies*100)

kable(metric_summary,
      caption = "Metric.Type Summary"
)%>%
  kable_styling(full_width=FALSE)
write_csv(metric_summary,"Metric_Type_Summary.csv")

# 1) Reverse the Metric.Type ordering
metric_levels <- rev(sort(unique(df$Metric.Type)))

# 2) Plot with the same formatting as your other figures
cr <- range(df$Effect.Size, na.rm = TRUE)

ggplot(df, aes(
  x = Effect.Size,
  y = factor(Metric.Type, levels = metric_levels)
)) +
  geom_boxplot(
    outlier.shape = NA,
    fill          = "gray80",
    alpha         = 0.7
  ) +
  geom_jitter(
    height = 0.15,    # ← jitter vertically around each category
    alpha  = 0.4,
    color  = "black"
  ) +
  stat_summary(
    fun  = mean,
    geom = "point",
    color= "blue",
    size = 3,
    shape= 18
  ) +
  labs(
    x = "Effect Size",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.title.x     = element_text(size  = 24, margin = margin(t = 10)),
    axis.text.x      = element_text(size  = 20),
    axis.text.y      = element_text(size  = 20),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin      = margin(15, 15, 15, 15),
    legend.position  = "none"
  ) +
  geom_vline(         # now a vertical reference line
    xintercept = c(-0.2, 0.2),
    linetype   = "dashed",
    color      = "red", 
    linewidth = 1.5
  ) +
  coord_cartesian(    # apply the same x-limits
    xlim = cr,
    clip = "off"
  )

# -------------------------------------------------------------------------------------
# 12. Structural vs. Functional boxplots
# -------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# 1) Identify the “top” Structural and Functional metrics from the full Data
# ----------------------------------------------------------------------------
top_str_metrics <- Data %>%
  filter(Metric.Group == "Structural", !is.na(Effect.Size)) %>%
  count(Response.Metric) %>%
  filter(n > 8) %>%
  pull(Response.Metric)

top_fun_metrics <- Data %>%
  filter(Metric.Group == "Functional", !is.na(Effect.Size)) %>%
  count(Response.Metric) %>%
  filter(n > 2) %>%
  pull(Response.Metric)

# ----------------------------------------------------------------------------
# 2) Subset your “df” (the cleaned & BMP‐filtered data) to only those top metrics
# ----------------------------------------------------------------------------
Data_structural <- df %>%
  filter(Response.Metric %in% top_str_metrics)

Data_functional <- df %>%
  filter(Response.Metric %in% top_fun_metrics)

# Optional: re‐export summaries of those top metrics
most_represented_structural <- Data_structural %>%
  group_by(Response.Metric) %>%
  summarise(
    Number_of_Studies = n_distinct(RefID),
    Mean_ES           = mean(Effect.Size, na.rm=TRUE),
    Median_ES         = median(Effect.Size, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(Percent = Number_of_Studies/total_studies*100)

most_represented_functional <- Data_functional %>%
  group_by(Response.Metric) %>%
  summarise(
    Number_of_Studies = n_distinct(RefID),
    Mean_ES           = mean(Effect.Size, na.rm=TRUE),
    Median_ES         = median(Effect.Size, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(Percent = Number_of_Studies/total_studies*100)

write_csv(most_represented_structural, "most_represented_structural.csv")
write_csv(most_represented_functional, "most_represented_functional.csv")

# ----------------------------------------------------------------------------
# 3) Compute common x‐limits (and force lower bound ≤ –1)
# ----------------------------------------------------------------------------
common_range <- range(
  c(Data_structural$Effect.Size, Data_functional$Effect.Size),
  na.rm = TRUE
)
common_range[1] <- min(common_range[1], -1)


# ----------------------------------------------------------------------------
# 4) Structural metrics boxplot, metrics in alphabetical order
# ----------------------------------------------------------------------------
metric_levels_str <- sort(unique(Data_structural$Response.Metric))

plot_structural <- ggplot(
  Data_structural,
  aes(
    x = Effect.Size,
    y = factor(Response.Metric, levels = metric_levels_str)
  )
) +
  geom_boxplot(outlier.shape = NA, fill = "gray80", alpha = 0.7) +
  geom_jitter(
    height = 0.15, alpha = 0.4, color = "black"
  ) +
  stat_summary(
    fun   = mean, geom = "point",
    color = "blue", size = 3, shape = 18
  ) +
  labs(x = "Effect Size", y = NULL) +
  theme_minimal() +
  theme(
    axis.title.x    = element_text(size = 28, margin = margin(t = 10)),
    axis.text.x     = element_text(size = 24),
    axis.text.y     = element_text(size = 24),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin     = margin(15, 15, 15, 15),
    legend.position = "none"
  ) +
  geom_vline(
    xintercept = c(-0.2, 0.2),
    linetype   = "dashed",
    color      = "red",
    linewidth  = 1.5
  ) +
  coord_cartesian(
    xlim = common_range,
    clip = "off"
  )

print(plot_structural)

# ----------------------------------------------------------------------------
# 5) Functional metrics boxplot, metrics alphabetical
# ----------------------------------------------------------------------------
metric_levels_fun <- sort(unique(Data_functional$Response.Metric))

plot_functional <- ggplot(
  Data_functional,
  aes(
    x = Effect.Size,
    y = factor(Response.Metric, levels = metric_levels_fun)
  )
) +
  geom_boxplot(outlier.shape = NA, fill = "gray80", alpha = 0.7) +
  geom_jitter(
    height = 0.15, alpha = 0.4, color = "black"
  ) +
  stat_summary(
    fun   = mean, geom = "point",
    color = "blue", size = 3, shape = 18
  ) +
  labs(x = "Effect Size", y = NULL) +
  theme_minimal() +
  theme(
    axis.title.x    = element_text(size = 28, margin = margin(t = 10)),
    axis.text.x     = element_text(size = 24),
    axis.text.y     = element_text(size = 24),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin     = margin(15, 15, 15, 15),
    legend.position = "none"
  ) +
  geom_vline(
    xintercept = c(-0.2, 0.2),
    linetype   = "dashed",
    color      = "red",
    linewidth  = 1.5
  ) +
  coord_cartesian(
    xlim = common_range,
    clip = "off"
  )

print(plot_functional)

combined <- plot_grid(
  plot_structural,
  plot_functional,
  ncol        = 1,
  align       = "v",    # align vertically
  axis        = "l",    # line up left (y) axes
  rel_heights = c(1,1)  # equal heights
)

ggsave("Fig_Structural_vs_Functional.png",
       combined,
       width  = 12,
       height = 16,  # make it tall enough for 2 panels
       dpi    = 300)
# -------------------------------------------------------------------------------------
# 13. Model fitting, diagnostics, ANOVA, AIC/AICc, and coefficient plots
# -------------------------------------------------------------------------------------

# 1) Fit your full model
full_mod <- lm(
  Effect.Size ~ PROVINCE * Threat.Type * BMP_Acronym * Metric.Type,
  data = df
)

# 1a) Diagnostic plots
tryCatch({
  par(mfrow = c(2, 2))
  plot(full_mod)
  par(mfrow = c(1, 1))
}, error = function(e) {
  message("Plotting window might be too small, please resize and try again.")
})

# 1b) ANOVA and summary of the full model
anova_full <- anova(full_mod)
summary_full <- summary(full_mod)

# 1c) Stepwise selection using AIC and then ANOVA/summary of reduced model
red_mod <- step(full_mod, direction = "both", trace = FALSE)
anova_red <- anova(red_mod)
summary_red <- summary(red_mod)


# 2) AIC/AICc table
aic_tab <- tibble(
  Model = c("Full", "Reduced"),
  AIC   = c(AIC(full_mod),  AIC(red_mod)),
  AICc  = c(AICc(full_mod), AICc(red_mod))
)
write_csv(aic_tab, "Model_AIC_Comparison.csv")
kable(aic_tab, caption = "AIC & AICc") %>% 
  kable_styling(full_width = FALSE)

# 2a)  Overall model fit (F & p) for Table 4
library(broom)

# helper: pull overall F & p from lm objects
model_overall <- function(mod){
  s <- glance(mod)
  tibble(
    Df          = s$df.residual + s$df,          # total df
    Resid_Df    = s$df.residual,
    F_statistic = s$statistic,
    p_value     = s$p.value
  )
}

overall_full    <- model_overall(full_mod)    %>% mutate(Model = "Full")
overall_reduced <- model_overall(red_mod)     %>% mutate(Model = "Reduced")

overall_tbl <- bind_rows(overall_full, overall_reduced) %>%
  select(Model, Df, Resid_Df, F_statistic, p_value)

# 2b  Append to existing AIC/AICc table
table4 <- aic_tab %>%                       # created earlier in your script
  left_join(overall_tbl, by = "Model") %>%
  relocate(AIC:AICc, .after = Model) %>%    # order columns neatly
  mutate(across(where(is.numeric), ~round(., 2)))

# 3) Helper to clean up each term into a “friendly” name
rename_term <- function(term) {
  if (term == "(Intercept)") return("Intercept")
  term <- gsub("`", "", term)
  term <- gsub("_", " ", term)
  term <- gsub("^PROVINCE", "", term, ignore.case = TRUE)
  term <- gsub("Threat\\.Type", "", term, ignore.case = TRUE)
  term <- gsub("Metric\\.Type|MetricType", "", term, ignore.case = TRUE)
  term <- trimws(term)
  
  province_map <- c(
    "COASTAL PLAIN"    = "Coastal Plain",
    "VALLEY AND RIDGE" = "Valley & Ridge",
    "PIEDMONT"         = "Piedmont"
  )
  land_use_map <- c(
    "AGRICULTURE"   = "Agriculture",
    "URBANIZATION"  = "Urbanization",
    "MIXED"         = "Mixed"
  )
  metric_map <- c(
    "RICHNESS"             = "Richness",
    "COMPOSITION METRICS"  = "Composition",
    "TOLERANCE METRICS"    = "Tolerance",
    "BIOTIC INDICES"       = "Biotic Indices",
    "HABITAT PREFERENCES"  = "Habitat Prefs",
    "BIOMASS"              = "Biomass",
    "DENSITY"              = "Density",
    "ISOTOPIC MEASUREMENTS"= "Isotopes"
  )
  all_map <- c(province_map, land_use_map, metric_map)
  
  if (grepl(":", term, fixed = TRUE)) {
    parts <- strsplit(term, ":")[[1]]
    if (length(parts) == 2) {
      p1 <- toupper(trimws(parts[1])); p2 <- toupper(trimws(parts[2]))
      p1c <- if (p1 %in% names(all_map)) all_map[p1] else tools::toTitleCase(parts[1])
      p2c <- if (p2 %in% names(all_map)) all_map[p2] else tools::toTitleCase(parts[2])
      return(paste0(p1c, " × ", p2c))
    }
  }
  
  tu <- toupper(term)
  if (tu %in% names(all_map)) return(all_map[tu])
  tools::toTitleCase(term)
}


# 4) Build a reusable plotting function
coef_plot <- function(model, ylim = c(-4, 4)) {
  # extract raw coefficients
  df_coef <- as.data.frame(summary(model)$coefficients)
  colnames(df_coef) <- c("Estimate","StdError","tValue","pValue")
  df_coef$Term      <- rownames(df_coef)
  df_coef$CleanTerm <- sapply(df_coef$Term, rename_term)
  
  # order by Estimate
  df_coef$CleanTerm <- factor(
    df_coef$CleanTerm,
    levels = df_coef$CleanTerm[order(df_coef$Estimate)]
  )
  
  ggplot(df_coef, aes(x = CleanTerm, y = Estimate)) +
    geom_point(color = "blue", size = 3) +
    geom_errorbar(aes(
      ymin = Estimate - StdError,
      ymax = Estimate + StdError
    ), width = 0.2, color = "black") +
    geom_hline(yintercept = c(-0.2, 0.2),
               linetype   = "dashed",
               color      = "red",
               linewidth  = 1) +
    coord_flip(ylim = ylim) +
    labs(x = NULL, y = "Coefficient Estimate") +
    theme_minimal() +
    theme(
      axis.title.x    = element_text(size = 20, margin = margin(t = 10)),
      axis.text.x     = element_text(size = 16),
      axis.text.y     = element_text(size = 16),
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 1),
      plot.margin     = margin(15, 15, 15, 15),
      legend.position = "none"
    )
}


# 5) Draw and save coefficient‐summary plots
cr <- c(-4, 4)

p_full <- coef_plot(full_mod, ylim = cr)
p_red  <- coef_plot(red_mod,  ylim = cr)

print(p_full)
print(p_red)

ggsave("Full_Model_Coefficients.png",    p_full, width = 10, height = 6, dpi = 300)
ggsave("Reduced_Model_Coefficients.png", p_red,  width = 10, height = 6, dpi = 300)

# 6) Regions–only main effects
coef_data <- broom::tidy(full_mod) %>%
  rename(Estimate = estimate, StdError = std.error, Term = term)
prov_main <- coef_data %>%
  filter(grepl("^PROVINCE", Term) & !grepl(":", Term)) %>%
  mutate(CleanTerm = vapply(Term, rename_term, ""))

# ensure reference level appears
if (!any(prov_main$CleanTerm == "Appalachian Plateaus")) {
  ref <- coef_data %>% filter(Term == "(Intercept)")
  prov_main <- bind_rows(
    tibble(
      Estimate  = ref$Estimate,
      StdError  = ref$StdError,
      Term      = "PROVINCEAPPALACHIAN PLATEAUS",
      CleanTerm = "Appalachian Plateaus"
    ),
    prov_main
  )
}

p_region <- ggplot(prov_main, aes(
  x = reorder(CleanTerm, Estimate),
  y = Estimate
)) +
  geom_point(color = "blue", size = 3) +
  geom_errorbar(aes(
    ymin = Estimate - StdError,
    ymax = Estimate + StdError
  ), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1.5) +
  coord_flip() +
  labs(x = "Physiographic Province", y = "Coefficient Estimate") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 24, margin = margin(t = 10)),
    axis.title.y = element_text(size = 24, margin = margin(r = 10)),
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin  = margin(15, 15, 15, 15),
    legend.position = "none"
  )

print(p_region)
ggsave("Regions_Only_Coefficients.png", p_region, width = 10, height = 6, dpi = 300)

# 7) Province × Metric.Type interactions
prov_met <- coef_data %>%
  filter(grepl("^PROVINCE.*Metric\\.Type", Term)) %>%
  mutate(CleanTerm = vapply(Term, rename_term, ""))

p_region_metric <- ggplot(prov_met, aes(
  x = reorder(CleanTerm, Estimate),
  y = Estimate
)) +
  geom_point(color = "blue", size = 3) +
  geom_errorbar(aes(
    ymin = Estimate - StdError,
    ymax = Estimate + StdError
  ), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1.5) +
  coord_flip() +
  labs(x = "Province × Metric Type", y = "Coefficient Estimate") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 24, margin = margin(t = 10)),
    axis.title.y = element_text(size = 24, margin = margin(r = 10)),
    axis.text.x  = element_text(size = 20),
    axis.text.y  = element_text(size = 20),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin  = margin(15, 15, 15, 15),
    legend.position = "none"
  )

print(p_region_metric)
ggsave("Province_Metric_Interactions.png", p_region_metric, width = 10, height = 6, dpi = 300)


# -------------------------------------------------------------------------------------
# 14. Generic summary‐table exports
# -------------------------------------------------------------------------------------
summarize_table <- function(df, group_vars, file_name){
  out <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      Studies   = n_distinct(RefID),
      Mean_ES   = mean(Effect.Size, na.rm = TRUE),
      Median_ES = median(Effect.Size, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    mutate(
      Percent = Studies / total_studies * 100
    ) %>%
    arrange(desc(Percent))
  
  write_csv(out, paste0(file_name, ".csv"))
  return(out)
}

tbl1 <- summarize_table(df, c("Year", "BMP_Acronym"),          "Summary_By_Year_BMP")
tbl2 <- summarize_table(df, c("Threat.Type", "PROVINCE"),     "Summary_By_LandUse_Province")
tbl3 <- summarize_table(df, c("Metric.Group", "Metric.Type"), "Summary_By_MetricGroup_Type")
tbl4 <- summarize_table(df, c("PROVINCE", "Threat.Type", "BMP_Acronym",
                              "Metric.Group", "Metric.Type", "Response.Metric"),
                        "Summary_By_MetricGroup_Type")
tbl5 <- summarize_table(df, c("PROVINCE", "Threat.Type", "BMP_Acronym", "Metric.Type"),
                        "Summary_By_MetricType")
tbl6 <- summarize_table(df, c("Threat.Type", "BMP_Acronym"),
                        "Summary_By_MetricType")
# 15) Print with kable, naming the caption explicitly:
kable(tbl4, caption = "By Metric Group / Response Metric") %>%
  kable_styling(full_width = FALSE)

kable(tbl5, caption = "By Province, Threat, BMP & Metric Type") %>%
  kable_styling(full_width = FALSE)

write_csv(tbl5, "Province_Threat_BMP_Metric Type.csv")

# -------------------------------------------------------------------------------------
# 15. Create Literature Review Table
# -------------------------------------------------------------------------------------
lit_review <- df %>%
  select(Paper.Title, Physiographic.Province, Effect.Size) %>%
  group_by(Paper.Title, Physiographic.Province) %>%
  summarise(
    Mean_Effect_Size   = mean(Effect.Size, na.rm = TRUE),
    Median_Effect_Size = median(Effect.Size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  distinct()

kable(lit_review, caption = "List of Papers Used in the Literature Review") %>%
  kable_styling(full_width = FALSE)

write_csv(lit_review, "Literature_Review_Papers.csv")

lit_review_metrics <- df %>%
  select(Paper.Title, Metric.Type, Effect.Size) %>%
  group_by(Paper.Title, Metric.Type) %>%
  summarise(
    Mean_Effect_Size   = mean(Effect.Size, na.rm = TRUE),
    Median_Effect_Size = median(Effect.Size, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  distinct()
# -------------------------------------------------------------------------------------
# End of script — all steps in one file, fully reproducible.
# -------------------------------------------------------------------------------------
