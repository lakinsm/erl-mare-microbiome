## Code to prepare publication figures for the Lakin et al. JEVS manuscript

library(vegan)
library(data.table)
library(ggplot2)
library(metagenomeSeq)

setwd('/home/lakinsm/Documents/ferrisBioinformatics/JEVS_manuscript/')

source('publication_functions.R')

set.seed(154)  # Seed the RNG, necessary for reproducibility


metadata_filepath = 'ERL_microbiome_metadata.csv'


# In which column of the metadata file are the sample IDs stored?
sample_column_id = 'ID'

##########################
## Import & Format Data ##
##########################
## These files should be standard for all analyses, as they are
## the output matrices from AMR++ nextflow.  Additionally,
## you will need to obtain the most recent megares annotations file
## from megares.meglab.org

# Load the data, MEGARes annotations, and metadata
temp_kraken <- read.table('analytic_data/ERL_microbiome_analytic_data.csv', header=T, row.names=1, sep=',')
num_features <- apply(temp_kraken, 2, function(x) sum(x > 0))
temp_kraken <- temp_kraken[, num_features > 1]
temp_kraken <- temp_kraken[!grepl('Unassigned', rownames(temp_kraken)), ]
kraken <- newMRexperiment(temp_kraken[rowSums(temp_kraken) > 0, ])

metadata <- read.csv(metadata_filepath, header=T)
metadata[, sample_column_id] <- make.names(metadata[, sample_column_id])


# Calculate normalization factors on the analytic data.
# We use Cumulative Sum Scaling as implemented in metagenomeSeq.
# You will likely get a warning about this step, but it's safe to ignore
cumNorm(kraken)


# Extract the normalized counts into data tables for aggregation
kraken_norm <- data.table(MRcounts(kraken, norm=T))
kraken_raw <- data.table(MRcounts(kraken, norm=F))


# Aggregate the kraken data using the rownames:
# this set of commands splits the rownames into their taxonomic levels and
# fills empty values with NA.  We then join that taxonomy data table with
# the actual data and aggregate using lapply as before.
kraken_taxonomy <- data.table(id=rownames(kraken))
setDT(kraken_taxonomy)[, c('Domain',
                           'Phylum',
                           'Class',
                           'Order',
                           'Family',
                           'Genus',
                           'Species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]
setkey(kraken_taxonomy, id)
kraken_norm[, id :=(rownames(kraken)), ]
setkey(kraken_norm, id)
kraken_norm <- kraken_taxonomy[kraken_norm]  # left outer join

kraken_raw[, id :=(rownames(kraken)), ]
setkey(kraken_raw, id)
kraken_raw <- kraken_taxonomy[kraken_raw]  # left outer join


# Group the kraken data by level for analysis, removing NA entries
kraken_domain <- kraken_norm[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:8]
kraken_domain_analytic <- newMRexperiment(counts=kraken_domain[, .SD, .SDcols=!'Domain'])
rownames(kraken_domain_analytic) <- kraken_domain$Domain

kraken_domain_raw <- kraken_raw[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:8]
kraken_domain_raw_analytic <- newMRexperiment(counts=kraken_domain_raw[, .SD, .SDcols=!'Domain'])
rownames(kraken_domain_raw_analytic) <- kraken_domain_raw$Domain

kraken_phylum <- kraken_norm[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:8]
kraken_phylum_analytic <- newMRexperiment(counts=kraken_phylum[, .SD, .SDcols=!'Phylum'])
rownames(kraken_phylum_analytic) <- kraken_phylum$Phylum

kraken_phylum_raw <- kraken_raw[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:8]
kraken_phylum_raw_analytic <- newMRexperiment(counts=kraken_phylum_raw[, .SD, .SDcols=!'Phylum'])
rownames(kraken_phylum_raw_analytic) <- kraken_phylum_raw$Phylum

kraken_class <- kraken_norm[!is.na(Class) & Class != 'NA', lapply(.SD, sum), by='Class', .SDcols=!1:8]
kraken_class_analytic <- newMRexperiment(counts=kraken_class[, .SD, .SDcols=!'Class'])
rownames(kraken_class_analytic) <- kraken_class$Class

kraken_class_raw <- kraken_raw[!is.na(Class) & Class != 'NA', lapply(.SD, sum), by='Class', .SDcols=!1:8]
kraken_class_raw_analytic <- newMRexperiment(counts=kraken_class_raw[, .SD, .SDcols=!'Class'])
rownames(kraken_class_raw_analytic) <- kraken_class_raw$Class

kraken_order <- kraken_norm[!is.na(Order) & Order != 'NA', lapply(.SD, sum), by='Order', .SDcols=!1:8]
kraken_order_analytic <- newMRexperiment(counts=kraken_order[, .SD, .SDcols=!'Order'])
rownames(kraken_order_analytic) <- kraken_order$Order

kraken_order_raw <- kraken_raw[!is.na(Order) & Order != 'NA', lapply(.SD, sum), by='Order', .SDcols=!1:8]
kraken_order_raw_analytic <- newMRexperiment(counts=kraken_order_raw[, .SD, .SDcols=!'Order'])
rownames(kraken_order_raw_analytic) <- kraken_order_raw$Order

kraken_family <- kraken_norm[!is.na(Family) & Family != 'NA', lapply(.SD, sum), by='Family', .SDcols=!1:8]
kraken_family_analytic <- newMRexperiment(counts=kraken_family[, .SD, .SDcols=!'Family'])
rownames(kraken_family_analytic) <- kraken_family$Family

kraken_family_raw <- kraken_raw[!is.na(Family) & Family != 'NA', lapply(.SD, sum), by='Family', .SDcols=!1:8]
kraken_family_raw_analytic <- newMRexperiment(counts=kraken_family_raw[, .SD, .SDcols=!'Family'])
rownames(kraken_family_raw_analytic) <- kraken_family_raw$Family

kraken_genus <- kraken_norm[!is.na(Genus) & Genus != 'NA', lapply(.SD, sum), by='Genus', .SDcols=!1:8]
kraken_genus_analytic <- newMRexperiment(counts=kraken_genus[, .SD, .SDcols=!'Genus'])
rownames(kraken_genus_analytic) <- kraken_genus$Genus

kraken_genus_raw <- kraken_raw[!is.na(Genus) & Genus != 'NA', lapply(.SD, sum), by='Genus', .SDcols=!1:8]
kraken_genus_raw_analytic <- newMRexperiment(counts=kraken_genus_raw[, .SD, .SDcols=!'Genus'])
rownames(kraken_genus_raw_analytic) <- kraken_genus_raw$Genus

kraken_species <- kraken_norm[!is.na(Species) & Species != 'NA', lapply(.SD, sum), by='Species', .SDcols=!1:8]
kraken_species_analytic <- newMRexperiment(counts=kraken_species[, .SD, .SDcols=!'Species'])
rownames(kraken_species_analytic) <- kraken_species$Species

kraken_species_raw <- kraken_raw[!is.na(Species) & Species != 'NA', lapply(.SD, sum), by='Species', .SDcols=!1:8]
kraken_species_raw_analytic <- newMRexperiment(counts=kraken_species_raw[, .SD, .SDcols=!'Species'])
rownames(kraken_species_raw_analytic) <- kraken_species_raw$Species


# Make long data frame for plotting with ggplot2
kraken_melted_analytic <- rbind(melt_dt(MRcounts(kraken_domain_analytic), 'Domain'),
                                melt_dt(MRcounts(kraken_phylum_analytic), 'Phylum'),
                                melt_dt(MRcounts(kraken_class_analytic), 'Class'),
                                melt_dt(MRcounts(kraken_order_analytic), 'Order'),
                                melt_dt(MRcounts(kraken_family_analytic), 'Family'),
                                melt_dt(MRcounts(kraken_genus_analytic), 'Genus'),
                                melt_dt(MRcounts(kraken_species_analytic), 'Species'))
kraken_melted_raw_analytic <- rbind(melt_dt(MRcounts(kraken_domain_raw_analytic), 'Domain'),
                                    melt_dt(MRcounts(kraken_phylum_raw_analytic), 'Phylum'),
                                    melt_dt(MRcounts(kraken_class_raw_analytic), 'Class'),
                                    melt_dt(MRcounts(kraken_order_raw_analytic), 'Order'),
                                    melt_dt(MRcounts(kraken_family_raw_analytic), 'Family'),
                                    melt_dt(MRcounts(kraken_genus_raw_analytic), 'Genus'),
                                    melt_dt(MRcounts(kraken_species_raw_analytic), 'Species'))

# Ensure that the metadata entries match the factor order of the MRexperiments
metadata <- data.table(metadata[match(colnames(MRcounts(kraken_domain_analytic)), metadata[, sample_column_id]), ])
setkeyv(metadata, sample_column_id)


# Vector of objects for iteration and their names
kraken_analytic_data <- c(kraken_domain_analytic,
                          kraken_phylum_analytic,
                          kraken_class_analytic,
                          kraken_order_analytic,
                          kraken_family_analytic,
                          kraken_genus_analytic,
                          kraken_species_analytic)
kraken_analytic_names <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
kraken_raw_analytic_data <- c(kraken_domain_raw_analytic,
                              kraken_phylum_raw_analytic,
                              kraken_class_raw_analytic,
                              kraken_order_raw_analytic,
                              kraken_family_raw_analytic,
                              kraken_genus_raw_analytic,
                              kraken_species_raw_analytic)
kraken_raw_analytic_names <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

for( l in 1:length(kraken_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(kraken_analytic_data[[l]])), metadata[[sample_column_id]])
    pData(kraken_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(kraken_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(kraken_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_analytic_data[[l]])))
    rownames(fData(kraken_analytic_data[[l]])) <- rownames(MRcounts(kraken_analytic_data[[l]]))
}

for( l in 1:length(kraken_raw_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(kraken_raw_analytic_data[[l]])), metadata[[sample_column_id]])
    pData(kraken_raw_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(kraken_raw_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(kraken_raw_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_raw_analytic_data[[l]])))
    rownames(fData(kraken_raw_analytic_data[[l]])) <- rownames(MRcounts(kraken_raw_analytic_data[[l]]))
}


###############################
## Alpha rarefaction phylum ###
###############################
publication_alpha_rarefaction(data_list=list(kraken_raw_analytic_data[[2]]),
                      data_names=c('Phylum'),
                      metadata=metadata,
                      sample_var=sample_column_id,
                      group_var='SamplingLocation',
                      analysis_subset=list('Horse != HBlank'),
                      outdir='figures/',
                      data_type='Microbiome')

##################
## NMDS phylum ###
##################
publication_ordination(data_list=list(kraken_raw_analytic_data[[2]]),
                              data_names=c('Phylum'),
                              metadata=metadata,
                              sample_var=sample_column_id,
                              hull_var='SamplingLocation',
                              analysis_subset=list('Horse != HBlank'),
                              outdir='figures/',
                              data_type='Microbiome',
                              method='NMDS')


###################################
## Pathogen Bar Graph and Tables ##
###################################
dat <- data.table(
    read.csv('/home/lakinsm/Documents/ferrisBioinformatics/2017Feb14_microbiome/stats/Microbiome/LocationDzFixedHorseRandom/LocationDzFixedHorseRandom_Microbiome_Genus_SamplingLocationClitoris - SamplingLocationVestibule_Model_Contrasts.csv')
)

dat[dat$Node.Name %in% c('Klebsiella', 'Escherichia', 'Pseudomonas', 'Streptococcus')] -> strep
dat[dat$Node.Name %in% c('Klebsiella', 'Escherichia', 'Pseudomonas', 'Streptococcus') & grepl('Uterus', dat$Contrast) & dat$logFC < 0] -> uterine

setkeyv(strep, c('Node.Name', 'Contrast'))
strep <- strep[, .SD, .SDcols=c('Node.Name', 'Contrast', 'logFC', 'CI.L', 'CI.R', 'adj.P.Val')]
names(strep) <- c('Genus', 'Comparison', 'Log2 Fold Change', 'Confidence Interval Left', 'Confidence Interval Right', 'Adjusted P Value')
strep$Comparison <- gsub('SamplingLocation', '', strep$Comparison)
strep$Comparison <- gsub('Clitoris', 'Clitoral Fossa', strep$Comparison)
strep$Comparison <- gsub('Vagina', 'Vaginal Vault', strep$Comparison)
strep$Comparison <- gsub('Vestibule', 'Vaginal Vestibule', strep$Comparison)


setkeyv(uterine, c('Node.Name', 'Contrast'))
uterine <- uterine[, .SD, .SDcols=c('Node.Name', 'Contrast', 'logFC', 'CI.L', 'CI.R', 'adj.P.Val')]
names(uterine) <- c('Genus', 'Comparison', 'Log2 Fold Change', 'Confidence Interval Left', 'Confidence Interval Right', 'Adjusted P Value')
uterine$Comparison <- gsub('SamplingLocation', '', uterine$Comparison)
uterine$Comparison <- gsub('Clitoris', 'Clitoral Fossa', uterine$Comparison)
uterine$Comparison <- gsub('Vagina', 'Vaginal Vault', uterine$Comparison)
uterine$Comparison <- gsub('Vestibule', 'Vaginal Vestibule', uterine$Comparison)

write.csv(uterine, 'figures/Table2.csv', row.names=F)

write.csv(strep, 'figures/Table3.csv', row.names=F)


pathogens <- data.table(cbind(rownames(kraken_genus_analytic), MRcounts(kraken_genus_analytic)))
names(pathogens)[1] <- 'Name'

mdt <- pathogens[Name %in% c('Klebsiella', 'Escherichia', 'Pseudomonas', 'Streptococcus')]
n <- mdt$Name
mdt <- mdt[, lapply(.SD, as.numeric), .SDcols=!'Name']
subset <- colnames(mdt)[colSums(mdt) > 0]
mdt <- mdt[, .SD, .SDcols=subset]
mdt <- cbind(n, mdt)
names(mdt)[1] <- 'Name'

mdt <- melt(mdt, id.vars='Name', variable.name='Sample', value.name='Abundance')
mdt <- mdt[Sample != 'HBlank']
setkey(mdt, Sample)

mdt <- metadata[mdt]

temp_mdt_labels <- as.character(mdt$SamplingLocation)
temp_mdt_labels <- gsub('Clitoris', 'Clitoral Fossa', temp_mdt_labels)
temp_mdt_labels <- gsub('Vagina', 'Vaginal Vault', temp_mdt_labels)
temp_mdt_labels <- gsub('Vestibule', 'Vaginal Vestibule', temp_mdt_labels)

mdt$SamplingLocation <- factor(temp_mdt_labels, levels=c('Clitoral Fossa', 'Vaginal Vestibule', 'Vaginal Vault', 'Uterus'),
                               ordered=T)

temp_meta <- data.table(ID=colnames(temp_kraken))
setkey(temp_meta, ID)
temp_meta <- metadata[temp_meta]
temp_meta_labels <- temp_meta$SamplingLocation
temp_meta_labels <- gsub('Clitoris', 'Clitoral Fossa', temp_meta_labels)
temp_meta_labels <- gsub('Vagina', 'Vaginal Vault', temp_meta_labels)
temp_meta_labels <- gsub('Vestibule', 'Vaginal Vestibule', temp_meta_labels)
temp_meta$SamplingLocation <- factor(temp_meta_labels, levels=c('Clitoral Fossa', 'Vaginal Vestibule', 'Vaginal Vault', 'Uterus'),
                                     ordered=T)
sample_nums <- table(temp_meta$SamplingLocation)

mdt[, MeanAbundance := ( as.numeric(Abundance / sample_nums[SamplingLocation]) )]

png('figures/Figure3.png', width=1200, height=900)
g <- ggplot(mdt, aes(x=SamplingLocation, y=MeanAbundance)) + geom_bar(stat='identity') +
    facet_wrap(~Name, ncol=1, scales='free_y') +
    xlab('\nSampling Location') + ylab('Mean Normalized Abundance\n') +
    theme(panel.grid.minor = element_blank(),
          strip.text.x=element_text(size=20),
          strip.text.y=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=20, angle=45, hjust=1),
          axis.title.x=element_text(size=30),
          axis.title.y=element_text(size=30),
          legend.position="bottom",
          panel.margin=unit(0.1, "lines"),
          plot.title=element_text(size=30, hjust=0.5),
          legend.text=element_text(size=18),
          legend.title=element_blank())
print(g)
dev.off()


############################
## Lactobacillaceae ########
############################
lact <- data.table(cbind(rownames(kraken_order_analytic), MRcounts(kraken_order_analytic)))
names(lact)[1] <- 'Name'

mdt <- lact[Name == 'Lactobacillales']
n <- mdt$Name
mdt <- mdt[, lapply(.SD, as.numeric), .SDcols=!'Name']
subset <- colnames(mdt)[colSums(mdt) > 0]
mdt <- mdt[, .SD, .SDcols=subset]
mdt <- cbind(n, mdt)
names(mdt)[1] <- 'Name'

mdt <- melt(mdt, id.vars='Name', variable.name='Sample', value.name='Abundance')
mdt <- mdt[Sample != 'HBlank']
setkey(mdt, Sample)

mdt <- metadata[mdt]

temp_mdt_labels <- as.character(mdt$SamplingLocation)
temp_mdt_labels <- gsub('Clitoris', 'Clitoral Fossa', temp_mdt_labels)
temp_mdt_labels <- gsub('Vagina', 'Vaginal Vault', temp_mdt_labels)
temp_mdt_labels <- gsub('Vestibule', 'Vaginal Vestibule', temp_mdt_labels)

mdt$SamplingLocation <- factor(temp_mdt_labels, levels=c('Clitoral Fossa', 'Vaginal Vestibule', 'Vaginal Vault', 'Uterus'),
                               ordered=T)

temp_meta <- data.table(ID=colnames(temp_kraken))
setkey(temp_meta, ID)
temp_meta <- metadata[temp_meta]
temp_meta_labels <- temp_meta$SamplingLocation
temp_meta_labels <- gsub('Clitoris', 'Clitoral Fossa', temp_meta_labels)
temp_meta_labels <- gsub('Vagina', 'Vaginal Vault', temp_meta_labels)
temp_meta_labels <- gsub('Vestibule', 'Vaginal Vestibule', temp_meta_labels)
temp_meta$SamplingLocation <- factor(temp_meta_labels, levels=c('Clitoral Fossa', 'Vaginal Vestibule', 'Vaginal Vault', 'Uterus'),
                                     ordered=T)
sample_nums <- table(temp_meta$SamplingLocation)

mdt[, MeanAbundance := ( as.numeric(Abundance / sample_nums[SamplingLocation]) )]

png('figures/Lactobacillales.png', width=1200, height=900)
g <- ggplot(mdt, aes(x=SamplingLocation, y=MeanAbundance)) + geom_bar(stat='identity') +
    xlab('\nSampling Location') + ylab('Mean Normalized Abundance\n') +
    theme(panel.grid.minor = element_blank(),
          strip.text.x=element_text(size=20),
          strip.text.y=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=20, angle=45, hjust=1),
          axis.title.x=element_text(size=30),
          axis.title.y=element_text(size=30),
          legend.position="bottom",
          panel.margin=unit(0.1, "lines"),
          plot.title=element_text(size=30, hjust=0.5),
          legend.text=element_text(size=18),
          legend.title=element_blank())
print(g)
dev.off()
