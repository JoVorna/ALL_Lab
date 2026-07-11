# Survival analysis

# 2026-06-11 Joonas Vornanen


# Start with clin file
library(data.table)
library(survival)
library(survminer)
setwd("/Users/joonasvornanen/kesatyo R/TCGA")

clin <- fread("luad_lusc_hnsc_clin.tsv")


colnames(clin)

cols_to_keep <- c("Patient ID", "Sample ID", "Overall Survival (Months)", "Overall Survival Status",
                  "TCGA PanCanAtlas Cancer Type Acronym")

clin <- clin[, ..cols_to_keep]

# Remove NA's from clin file 

clin <- na.omit(clin) #Approx 70 samples gone

# Prepare LUAD

luad <- fread('/Users/joonasvornanen/kesatyo R/TCGA/luad_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt')
luad <- as.data.frame(luad)
# Check gene duplicates
luad_duplicated <- luad[duplicated(luad$Hugo_Symbol), ]
luad_duplicated$Hugo_Symbol
#FGF13,ELMOD1,NKAIN3,PALM2AKAP2,QSOX1,SNAP47,TMEM8B only duplicated genes
rm(luad_duplicated)

luad <- luad[!duplicated(luad$Hugo_Symbol), ]


# Add genes as rownames

rownames(luad) <- luad$Hugo_Symbol

luad_cols_to_keep <- intersect(clin$`Sample ID`, colnames(luad))

luad <- luad[, luad_cols_to_keep]

# LUSC

lusc <- fread('/Users/joonasvornanen/kesatyo R/TCGA/lusc_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt')
lusc <- as.data.frame(lusc)
# Check gene duplicates
lusc_duplicated <- lusc[duplicated(lusc$Hugo_Symbol), ]
lusc_duplicated$Hugo_Symbol
#Same genes as before duplicated also here
rm(lusc_duplicated)

lusc <- lusc[!duplicated(lusc$Hugo_Symbol), ]

# Add genes as rownames

rownames(lusc) <- lusc$Hugo_Symbol

lusc_cols_to_keep <- intersect(clin$`Sample ID`, colnames(lusc))

lusc <- lusc[, lusc_cols_to_keep]

# HNSC

hnsc <- fread('/Users/joonasvornanen/kesatyo R/TCGA/hnsc_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt')
hnsc <- as.data.frame(hnsc)
# Check gene duplicates
hnsc_duplicated <- hnsc[duplicated(hnsc$Hugo_Symbol), ]
hnsc_duplicated$Hugo_Symbol
#Same as before
rm(hnsc_duplicated)

hnsc <- hnsc[!duplicated(hnsc$Hugo_Symbol), ]

# Add genes as rownames

rownames(hnsc) <- hnsc$Hugo_Symbol

hnsc_cols_to_keep <- intersect(clin$`Sample ID`, colnames(hnsc))

hnsc <- hnsc[, hnsc_cols_to_keep]



# Pre processing done

# Get common genes

common_genes <- Reduce(intersect, list(
  rownames(luad),
  rownames(lusc),
  rownames(hnsc)
))

luad <- luad[common_genes, ]
lusc <- lusc[common_genes, ]
hnsc <- hnsc[common_genes, ]

# Combine into 1 var

expr <- cbind(luad, lusc, hnsc)
expr <- na.omit(expr) # Remove NA's
# Sanity check -- Match expr file colnames with clin file Sample ID's

head(colnames(expr))
head(clin$`Sample ID`) # OK!

# Find overlap
samples <- intersect(colnames(expr), clin$`Sample ID`)

length(samples)


clin <- clin[match(colnames(expr), clin$`Sample ID`), ]

# Final sanity check

all(clin$`Sample ID` == colnames(expr)) # TRUE, OK

# Now process survival col

clin$status <- ifelse(
  grepl("DECEASED", clin$`Overall Survival Status`),
  1,
  0
)

clin$time <- as.numeric(clin$`Overall Survival (Months)`)
##############################################################
# Survival analysis###########################################
##############################################################

# Because expr and clin file sample indexes match, we can get gene expr like this:

gene <- "TNFRSF1A"

expr_gene <- as.numeric(expr[gene, ])

clin$expression <- expr_gene

#Optional way to make this data frame work
##############################################################
# 
# expr_df <- as.data.frame(t(expr))
# 
# 
# 
# 
# 
# 
# # Add sample IDs as a column
# expr_df$`Sample ID` <- rownames(expr_df)
# 
# 
# clin_merged <- merge(
#   clin,
#   expr_df,
#   by = "Sample ID"
# )
# clin_merged[[gene]]
# 
# 
# clin <- clin_merged
# ##############################################################


# Get the diseases you want:

#clin <- clin[grepl("LUSC", clin$`TCGA PanCanAtlas Cancer Type Acronym`), ]

# Get cutoff, quantile 75% for clin expression


# cutoff <- quantile(clin$expression, probs = 0.75)
# cutoff2 <- quantile(clin$expression, probs = 0.25)
# clin$group <- ifelse(
#   clin$expression >= cutoff,
#   "High",ifelse(
#     clin$expression <= cutoff2,
#     "Low",
#     NA
#   )
# )

# Cutoff at mean.

clin$group <- ifelse(
  clin$expression >= mean(clin$expression),
  "High",
  "Low"
)

# Fit survival model

fit <- survfit(
  Surv(time, status) ~ group,
  data = clin
)

# Plot

ggsurvplot(
  fit,
  data = clin,
  pval = TRUE,
  risk.table = TRUE
)

# Cox reg so we get hazard ratios

cox <- coxph(
  Surv(time, status) ~ expression,
  data = clin
)

summary(cox)





# TODO

# 2 wasy to make expr df, either like now, expr is a col
# or transpose whole df and merge wit clin and you get genes as cols also


