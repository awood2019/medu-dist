#analyzing medusozoan eDNA sequences from GoM

#load packages
library(dada2)

#--------code from GoMx phyla distribution script to make phyloseq object-------

#load packages
library(tidyverse) ; packageVersion("tidyverse") 
library(phyloseq) ; packageVersion("phyloseq") 
library(vegan) ; packageVersion("vegan") 
library(DESeq2) ; packageVersion("DESeq2") 
library(dendextend) ; packageVersion("dendextend") 
library(viridis) ; packageVersion("viridis") 
library("ggplot2")

#set working directory 
setwd("~/Desktop/AW RStudio/data/gomx-phy-dist")

#load components of phyloseq object: taxonomy table, count table, and sample data table
tax_tab <- read_csv("rep-seqs-phylum.csv", show_col_types = FALSE) #loading taxonomy table w/ ASVs, sequence, & phyla
count_tab <- read_delim("table.tsv") #loading count table w/ ASV counts for each sample
sample_info_tab <- read_csv("anth-28S-sampledata_20231016.csv") #loading sample data table w/ sample metadata

#coerce tables into proper format to make phyloseq object
#tax_tab_phy: includes taxonomic information for each representative (ASV) sequence
phylum <- tax_tab$Phylum #pulling out phylum column from taxonomy table
tax_tab_phy <- tibble(phylum) #making phyla into a tibble containing phylum for each sequence
tax_tab_phy <- as.matrix(tax_tab_phy) #make tibble into matrix
row.names(tax_tab_phy) <- tax_tab$Sequence #make sequence column the row names
tax_tab_phy[is.na(tax_tab_phy)] <- "< 85% similarity to top BLAST hit" #change NA values to more accurate description

#count_tab_phy: includes all ASVs and their abundances in each sample (row.names must match row.names of tax_tab_phy)
count_tab_phy <- select(count_tab, -"...1") #delete this weird column
row.names(count_tab_phy) <- count_tab$...1 #make sequences the row names (ignore warning message)

#sample_info_tab_phy: table that includes sample information for all samples (row.names must equal col.names in count table)
sample_info_tab <- sample_info_tab %>% mutate(depth_bin = cut_width(sample_info_tab$Depth, width = 10, boundary = 0)) #create column for depth range as a factor
sample_info_tab_phy <- sample_info_tab
sample_info_tab_phy <- sample_info_tab_phy[-c(55,56),] #delete the last 2 rows because they have NAs across the board
sample_data <- sample_data(sample_info_tab_phy) #convert to phyloseq component now because row names get changed by sample_data command
row.names(sample_data) <- sample_data$File.name #change row names to match file name

#make phyloseq object 
ASV_physeq <- phyloseq(otu_table(count_tab_phy, taxa_are_rows = TRUE), tax_table(tax_tab_phy), sample_data)
ASV_physeq <- prune_taxa(taxa_sums(ASV_physeq) > 0, ASV_physeq) #pruning out ASVs with zero counts
saveRDS(ASV_physeq, 'allphy_physeq.rds') #save phyloseq object

#transform phyloseq object to dataframe 
df_ASV_physeq <- ASV_physeq %>% psmelt() #melt phyloseq object to long dataframe
head(df_ASV_physeq)

#----------end of the code from GoM phyla distribution script--------

#set working directory 
setwd("~/Desktop/AW RStudio/data/gomx-medu-dist")

#import files
genbank_seqs <- read_csv("28S-Hyd-Scy-Cub-Staur-Genbank-barcode-nodecimal.csv") #file of 28S barcodes from sequences identified as Hydrozoa, Scyphozoa, Cubozoa, or Staurozoa in Genbank & their accession ID without a decimal
accid_taxid <- read.delim("accid_taxid.txt") #file of accession IDs and tax IDs
unique_taxid <- read.delim("unique_taxid.txt") #file of unique tax IDs
blast_tax <- read.delim("blastn_hits_taxonomy.txt") #file of tax ID and BLAST taxonomy

#create table of medusozoan sequences and their taxonomy by combining accesssion ID, sequence, and taxonomy
accid_taxid_unique <- accid_taxid %>% distinct(TaxID, .keep_all = TRUE) #some accession IDs have the same TaxID, so select unique TaxIDs that are the first occurrence in df & corresponding accession ID
medu_classifier <- left_join(genbank_seqs, accid_taxid_unique, by = "Name") #join unique TaxIDs with Genbank sequences by their accession ID
blast_tax_unique <- blast_tax %>% distinct(TaxID, .keep_all=TRUE) #multiple occurrences of the same TaxID, so select unique TaxIDs that are the first occurrence in df & corresponding tax classification
medu_classifier <- right_join(medu_classifier, blast_tax, by = "TaxID") #join df of sequences w/ taxonomic classification
length(unique(medu_classifier$Species)) #how many species-level sequences are in the classifier?

#save table to make assignTaxa and assignSpecies classifiers in CLI
write.csv(medu_classifier, "~/Desktop/AW RStudio/results/gomx-medu-dist/medu_preclassifier.csv")

#select only the cnidarians from phyloseq object
cnid.table <- df_ASV_physeq %>% 
  subset(phylum == "Cnidaria") %>%
  filter(Abundance >0) #only the samples that contained cnidarians

#create vector of only unique sequences from cnid.table so we can compare them to FASTA files
seqs <- unique(cnid.table$OTU)

#assignTaxonomy using FASTA file as reference
taxa <- assignTaxonomy(seqs, "28S-Medu_assignTaxonomy.fasta", multi=TRUE, minBoot = 80,
                       taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) #classifying to genus level with assignTaxonomy
unname(taxa)
unique(taxa [, 3]) #which classes were identified?

#addSpecies using FASTA file as reference
taxa.species <- addSpecies(taxa, "28S-Medu_assignSpecies.fasta", allowMultiple = TRUE)  #finding 100% matches to our reference database of Gulf of Mexico ctenos with assignSpecies
taxa.species <- subset(taxa.species, select =c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) #get rid of weird NA column
unique(taxa.species [, 6]) #which genera were identified?
unique(taxa.species [, 7]) #which species were identified?

#convert from vector of ASV taxonomy to dataframe
taxa.species.df <- as.data.frame(taxa.species)
taxa.species.df$seq <- row.names(taxa.species.df) #make sequence column from the row names
length(unique(taxa.species.df[["seq"]])) #how many unique ASVs are in table?

#save tables of ASVs identified to Genus level
setwd("~/Desktop/AW RStudio/results/gomx-medu-dist")
classified.genera <- taxa.species.df[!is.na(taxa.species.df$Genus), ] #select all ASVs identified to the Genus level
classified.genera.summary <- as.data.frame(table(taxa.species.df$Genus)) #create table of summed occurrences of identified genera
write.table(classified.genera, file = 'medu_classified_genera.tsv', sep = "\t", row.names = FALSE, quote=FALSE) #save table of ASVs identified to Genus level

#make table of abundance counts for each unique ctenophore ASV
counts_per_ASV <- cnid.table %>% 
  group_by(OTU) %>%
  summarize(totalcount = sum(Abundance)) 
names(counts_per_ASV)[names(counts_per_ASV) == "OTU"] <- "seq" #rename this column seq so we can combine it with taxa.species.df
names(counts_per_ASV)[names(counts_per_ASV) == "totalcount"] <- "unique abundance" #name this column so we know it represents the counts of unique ASVs present in dataset

#combine abundance counts table of cnidarian ASVs w/ classifier table
taxa.species.df <- left_join(taxa.species.df, counts_per_ASV, by = "seq")

#save as a table
write.table(taxa.species.df, file = 'Medu_classifier_results.tsv', sep = "\t", row.names = FALSE, quote=FALSE) #writing the ASV counts table with the taxonomic classifications of each cteno ASV

#--------------plotting bar plots of abundance with depth------------------

#create table from cnid.table with only depth, method, & abundance 
cnid.depth <- cnid.table %>%
  select(OTU, Abundance, depth_bin, CTD.ROV)

#join with table of medu species taxonomy
species.depth <- right_join(taxa.species.df, cnid.depth, by = c("seq" = "OTU")) 
species.depth <- species.depth %>% #take out all samples that were blanks or controls
  filter(CTD.ROV != "Negative" & CTD.ROV != "NTC")

#figure out which depth ranges were sampled at so we can make depth bins for plots
unique(species.depth$depth_bin)

#plot by species
ggplot(species.depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Species)) + #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  ggtitle("Depth Distribution of Hydrozoa, Scyphozoa, Cubozoa, and Staurozoa Genera in Gulf of Mexico")+
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

#plot by genus
ggplot(species.depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Genus)) + #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  ggtitle("Depth Distribution of Hydrozoa, Scyphozoa, Cubozoa, and Staurozoa Genera in Gulf of Mexico")+
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

#plot by family
ggplot(species.depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Family)) + #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") + #facet by sampling method
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  ggtitle("Depth Distribution of Hydrozoa, Scyphozoa, Cubozoa, and Staurozoa Families in Gulf of Mexico")+
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

#plot by order
ggplot(species.depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Order)) + #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  ggtitle("Depth Distribution of Hydrozoa, Scyphozoa, Cubozoa, and Staurozoa Orders in Gulf of Mexico")+
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

#---------------plotting cnidarian ASVs in ordination plot with UniFrac distance matrix---------------

#load packages needed for PCoA ordination plot
library(vegan)
library(dplyr)
library(ape)
library(phyloseq)
library(Biostrings)
library(ShortRead)
library(tidysq)
library(ggalt)
library(phytools)

#set working directory
setwd("~/Desktop/AW RStudio/data/gomx-medu-dist")

#make phyloseq object with only cnidarians
cnid_physeq = subset_taxa(ASV_physeq, phylum == "Cnidaria") 

#remove rows with 0 cnidarians so the distance matrix can be calculated
cnid_physeq_0 <- prune_samples(sample_sums(cnid_physeq) >0, cnid_physeq) 

#also remove rows with less than 10 and 100 ctenophores to see how ordination changes 
cnid_physeq_10 <- prune_samples(sample_sums(cnid_physeq) >10, cnid_physeq) 
cnid_physeq_100 <- prune_samples(sample_sums(cnid_physeq) >100, cnid_physeq) 

#export table of ASVs so we can put it into Geneious to make a phylogenetic tree
tax_cnid <- tax_table(cnid_physeq_0) #take component of ASVs out of phyloseq object
write.csv(tax_cnid, "~/Desktop/AW RStudio/data/gomx-medu-dist/tax_cnid.csv") #save as csv file

#add tree of ctenophores to each phyloseq object
cnid_tree <- read.nexus("~/Desktop/AW RStudio/data/gomx-medu-dist/tax_cnid alignment FastTree Tree.nex") #import Geneious tree of aligned ctenophores
rooted_tree <- midpoint.root(cnid_tree)
physeq_tree_abdndce0 <- merge_phyloseq(cnid_physeq_0, rooted_tree) #merge tree with existing phyloseq object for abundance >0
physeq_tree_abdndce10 <- merge_phyloseq(cnid_physeq_10, rooted_tree) #merge tree with existing phyloseq object for abundance >10
physeq_tree_abdndce100 <- merge_phyloseq(cnid_physeq_100, rooted_tree) #merge tree with existing phyloseq object for abundance >100

#calculate Unifrac distance matrix to each phyloseq object
uni_matrix_0 <- UniFrac(physeq_tree_abdndce0, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE) #distance matrix for phyloseq object with abundance >0
uni_matrix_10 <- UniFrac(physeq_tree_abdndce10, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE) #distance matrix for phyloseq object with abundance >10
uni_matrix_100 <- UniFrac(physeq_tree_abdndce100, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE) #distance matrix for phyloseq object with abundance >100

#perform ordination for each phyloseq object
cnid_ordination_0 <- ordinate(physeq_tree_abdndce0, method = "PCoA", uni_matrix_0, weighted=TRUE) #ordination for phyloseq object with abundance >0
cnid_ordination_10 <- ordinate(physeq_tree_abdndce10, method = "PCoA", uni_matrix_10, weighted=TRUE) #ordination for phyloseq object with abundance >10
cnid_ordination_100 <- ordinate(physeq_tree_abdndce100, method = "PCoA", uni_matrix_100, weighted=TRUE) #ordination for phyloseq object with abundance >100

#make PCOA plot for each phyloseq object
#plot for phyloseq object with abundance >0
plot_ordination(physeq_tree_abdndce0, cnid_ordination_0, color = "Depth", shape="CTD.ROV") + #define the point color by depth and point shape by sampling method 
  geom_encircle(aes(group=Site, fill=Site), alpha=0.2, s_shape=1, expand=0) + #add polygon shapes to enclose points by site
  theme_classic()+
  scale_color_continuous(type="viridis", option="D", direction=-1) +
  ggtitle("PCoA of Cnidarian ASVs with Abundance >0")

#plot for phyloseq object with abundance >10
plot_ordination(physeq_tree_abdndce10, cnid_ordination_10, color = "Depth", shape="CTD.ROV") + #define the point color by depth and point shape by sampling method 
  geom_encircle(aes(group=Site, fill=Site), alpha=0.2, s_shape=1, expand=0) + #add polygon shapes to enclose points by site
  theme_classic()+
  scale_color_continuous(type="viridis", option="D", direction=-1) +
  ggtitle("PCoA of Cnidarian ASVs with Abundance >10")

#plot for phyloseq object with abundance >100
plot_ordination(physeq_tree_abdndce100, cnid_ordination_100, color = "Depth", shape="CTD.ROV") + #define the point color by depth and point shape by sampling method 
  geom_encircle(aes(group=Site, fill=Site), alpha=0.2, s_shape=1, expand=0) + #add polygon shapes to enclose points by site
  theme_classic()+
  scale_color_continuous(type="viridis", option="D", direction=-1) +
  ggtitle("PCoA of Cnidarian ASVs with Abundance >100")

