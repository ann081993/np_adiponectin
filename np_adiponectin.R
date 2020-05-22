#### 2020-05-21 np_adiponectin
# load organism list
setwd("C:/Users/ASC/Desktop/projects/NP_adiponectin")
np_organism_list <- read.csv("np_organism_list.csv", stringsAsFactors = FALSE)
{
np_organism_list$organism <- gsub("\\(", " \\(", trimws(np_organism_list$organism))

genus <- sapply(strsplit(np_organism_list$organism, split = "\\s+"), "[[", 1)
species <- sapply(strsplit(np_organism_list$organism, split = "\\s+"), "[[", 2)

head(genus)
head(species)
organism <- unique(paste(genus, species))
}
organism <- np_organism_list$corrected_name

# organism to compounds (PubChem)
organism[1]
organism[2]
"https://pubchem.ncbi.nlm.nih.gov/#query=[{%22query%22:%22Ribes%20fasciculatum%22,%22collection%22:%22pubmed%22,%22id_type%22:%22cid%22}]&collection=compound"
"https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22compound%22,%22where%22:{%22ands%22:[{%22input%22:{%22type%22:%22netcachekey%22,%22key%22:%22g8Qhv7Z008jk5lv_2YcS7-aTgPPOIIA4-h2bdOEMiXXhFbU%22}}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_compound_composite_query_Ribes%20fasciculatum%22}"
"https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22compound%22,%22where%22:{%22ands%22:[{%22input%22:{%22type%22:%22netcachekey%22,%22key%22:%22g8Qhv7Z008jk5lv_2YcS7-aTgPPOIIA4-h2bdOEMiXXhFbU%22}}]},%22order%22:[%22relevancescore,desc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PubChem_compound_composite_query_Lespedeza%20bicolor%22}"

# organism to compounds (NPASS)
# load NPASS
readLines("NPASSv1.0_download_naturalProducts_species_pair.txt", 10)
sp_info <- read.csv("NPASSv1.0_download_naturalProducts_speciesInfo.txt", stringsAsFactors = FALSE, sep = "\t")
sp_np_pair <- read.csv("NPASSv1.0_download_naturalProducts_species_pair.txt", stringsAsFactors = FALSE, sep = "\t")
np_info <- read.csv("NPASSv1.0_download_naturalProducts_generalInfo.txt", stringsAsFactors = FALSE, sep = "\t")
np_prop <- read.csv("NPASSv1.0_download_naturalProducts_properties.txt", stringsAsFactors = FALSE, sep = "\t")

sp_info[grep("Chrysanthemum indicum", sp_info$org_name), ]

conversion <- function(df, from, from_list, to) {
        df <- df[, c(from, to)]
        return(df[df[, 1] %in% from_list, 2])
}

organism[!(organism %in% sp_info$org_name)]
table(organism %in% sp_info$org_name)
table(organism %in% sp_info$species_name)

org_id <- conversion(sp_info, "org_name", organism, "org_id")
conversion(sp_info, "org_id", org_id, "org_name")

np_id <- conversion(sp_np_pair, "org_id", org_id, "np_id")
np_cid <- conversion(np_info, "np_id", np_id, "pubchem_cid")
np_name <- conversion(np_info, "np_id", np_id, "pref_name")

cids <- NULL
smiles <- NULL
for(o in org_id) {
        np_id <- conversion(sp_np_pair, "org_id", o, "np_id")
        #smi <- conversion(np_prop, "np_id", np_id, "canonical_smiles")
        #smiles <- c(smiles, smi)
        cids <- c(cids, length(np_id))
}

source("https://raw.githubusercontent.com/ann081993/fp_util/master/fp_util_src.R") 
smi2fp <- function(smi, method = "pubchem") {
        result <- list()
        for (i in smi) {
                cat(i)
                parsed_smi <- suppressWarnings(parse.smiles(i)[[1]])
                fp <- try(get.fingerprint(parsed_smi, method), silent = TRUE)
                if(!is(fp, 'try-error') | !is.null(parsed_smi)) {
                        result <- append(result, fp)
                        cat("\tO", "\n")
                } else cat("\n")
        }
        if(length(result) > 0) result <- fp.to.matrix(result)
        return(result)
}

parse.smiles("O=c1cc2CCn3c2c(c1)c1cc2OCOc2cc1c3") # error example
parse.smiles("[O+]=c1cc2CCn3c2c(c1)c1cc2OCOc2cc1c3")
isvalid.formula("O=c1cc2CCn3c2c(c1)c1cc2OCOc2cc1c3")

smiles <- smiles[!grepl("[.]", smiles)]
smiles <- smiles[!grepl("[n+]", smiles)]
np_fp <- smi2fp(smiles)
fppca <- prcomp(np_fp)

# visualize
library(ggplot2)
library(ggrepel)
library(gplots)
library(viridis)
library(cluster)
library(ggpubr)

df_plot <- as.data.frame(cbind(fppca$x[, 1:2]))
ggscatter(df_plot, x = "PC1", y = "PC2")

# phylogenic analysis
# https://www.molecularecologist.com/2017/02/phylogenetic-trees-in-r-using-ggtree/
library(easyPubMed)
library(ape) ###########
library(msa) # including Biostrings
library(genbankr) # trash
library(XML)
# GenBank 
# https://peerj.com/articles/2546/ rbcL gene, GeneBank, BLASTN
# https://www.pnas.org/content/pnas/91/12/5730.full.pdf
# https://en.wikipedia.org/wiki/RuBisCO
# https://spectrum.library.concordia.ca/6741/1/Dayanandan_AnnalsMissouriBotanicalGardens_1993.pdf
# Biostrings: https://bioconductor.org/packages/release/bioc/vignettes/Biostrings/inst/doc/MultipleAlignments.pdf
get_genebank <- function (pubmed_query_string, api_key = NULL) {
        old_warn <- options()$warn
        options(warn = -1)
        t_0 <- Sys.time()
        myQuery <- as.character(pubmed_query_string)
        myQuery <- gsub(" ", "+", myQuery, fixed = TRUE)
        myPubmedURL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?", 
                             "db=nuccore&term=", myQuery, "&usehistory=y", 
                             sep = "")
        if (!is.null(api_key)) {
                myPubmedURL <- paste(myPubmedURL, "&api_key=", 
                                     api_key, sep = "")
        }
        idXML <- NULL
        try_num <- 1
        while (is.null(idXML)) {
                if (try_num > 1) 
                        Sys.sleep(time = 2)
                t_1 <- Sys.time()
                if (as.numeric(difftime(t_1, t_0, units = "mins")) > 
                    2) {
                        message("Killing the request! Something is not working. Please, try again later")
                        return()
                }
                idXML <- tryCatch({
                        IDconnect <- suppressWarnings(url(myPubmedURL, open = "rb", 
                                                          encoding = "UTF8"))
                        idXML <- suppressWarnings(readLines(IDconnect, warn = FALSE, 
                                                            encoding = "UTF8"))
                        idXML <- paste(idXML, collapse = "")
                        if (grepl("<ERROR>", substr(idXML, 1, 250))) {
                                NULL
                        }
                        else {
                                idXML
                        }
                }, error = function(e) {
                        NULL
                }, finally = {
                        try(suppressWarnings(close(IDconnect)), silent = TRUE)
                })
                myIDlist <- NULL
                if (!is.null(idXML)) {
                        tryCatch({
                                myIDlist <- list()
                                my_tags <- c("Count", "RetMax", "RetStart", 
                                             "QueryKey", "WebEnv", "IdList", 
                                             "TranslationSet", "QueryTranslation")
                                for (j in 1:length(my_tags)) {
                                        ttag <- my_tags[j]
                                        xx <- custom_grep(idXML, tag = ttag, "char")
                                        myIDlist[[ttag]] <- xx[1]
                                }
                                nutag <- "Id"
                                xx <- myIDlist[["IdList"]]
                                xx <- custom_grep(xx, "Id", format = "list")
                                names(xx) <- rep("Id", length(xx))
                                myIDlist[["IdList"]] <- xx
                                xx <- myIDlist[["TranslationSet"]]
                                myIDlist[["TranslationSet"]] <- list()
                                nutag <- c("From", "To")
                                for (z in nutag) {
                                        yy <- custom_grep(xx, z, format = "char")
                                        myIDlist[["TranslationSet"]][[z]] <- yy[1]
                                }
                        }, error = function(e) {
                                idXML <- NULL
                        })
                }
                if (!is.list(myIDlist)) {
                        idXML <- NULL
                }
                try_num <- try_num + 1
        }
        myIDlist[["OriginalQuery"]] <- myQuery
        myIDlist[["APIkey"]] <- api_key
        options(warn = old_warn)
        return(myIDlist)
}

gb_total <- NULL
labels <- NULL
n = 1

gb_total <- gb_total[1]
o <- organism[n]; n <- n + 1
gb_list <- get_genebank(paste0(o, " rbcL"))
gb_gb <- read.GenBank(unlist(gb_list$IdList), species.names = T)
gb_gb <- gb_gb[nchar(gb_gb) / 6 < 2000 & nchar(gb_gb) / 6 > 800][1]

if(length(gb_gb) == 1) {
        gb_total <- c(gb_total, gb_gb)
        labels <- c(labels, o)
}

######### filter sequence with NA or sequence with unknown nucleotide

write.dna(gb_total, file ="gb.fasta", format = "fasta")
seqs <- readDNAStringSet("gb.fasta") # names(seqs) <- labels
seqs_msa <- msaMuscle(seqs) # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0153008
seqs_msa <- msaClustalOmega(seqs)
seqs_msa <- msaConvert(seqs_msa, type = "ape::DNAbin")

dist.aa(seqs_msa)
labels(seqs_msa)
image(seqs_msa)
as.matrix(seqs_msa)

plot(hclust(dist.aa(seqs_msa)))
plot(hclust(dist.aa(seqs_msa)),
     xlab = "", ylab = "Number of differet nt in rbcL gene", cex = 1)

plot(as.dendrogram(hclust(dist.aa(seqs_msa)), hang=0.02),
     xlab = "Number of differet nt in rbcL gene", horiz = TRUE)

image(woodmouse)
dist.aa(woodmouse)

#### 2020-05-22 total analysis of NPASS
sp_info <- read.csv("NPASSv1.0_download_naturalProducts_speciesInfo.txt", stringsAsFactors = FALSE, sep = "\t")
sp_np_pair <- read.csv("NPASSv1.0_download_naturalProducts_species_pair.txt", stringsAsFactors = FALSE, sep = "\t")
np_info <- read.csv("NPASSv1.0_download_naturalProducts_generalInfo.txt", stringsAsFactors = FALSE, sep = "\t")
np_prop <- read.csv("NPASSv1.0_download_naturalProducts_properties.txt", stringsAsFactors = FALSE, sep = "\t")

t <- table(sp_info$genus_name)
hist(ifelse(t > 30, 30, t), breaks = 20)

max_sim <- function(dist_mat) { # max similarity for each compounds from fpjacdis
        c <- ncol(dist_mat)
        r <- nrow(dist_mat)
        if (c == r) dist_mat[!upper.tri(dist_mat) & !lower.tri(dist_mat)] <- 0
        return(do.call(pmax, data.frame(dist_mat)))
}

sp_info <- sp_info[sp_info$kingdom_name == "Viridiplantae", ]
# filter out species producing < 2 compounds or > 91 compounds (99%)
num_cpds <- sapply(sp_info$org_id, function(x) length(conversion(sp_np_pair, "org_id", x, "np_id")))
summary(num_cpds); quantile(num_cpds, probs = seq(0, 1, 0.01))
head(num_cpds[order(num_cpds, decreasing = T)])
conversion(sp_info, "org_id", "NPO1797", "org_name")
conversion(sp_np_pair, "org_id", "NPO1797", "np_id")

conversion(np_info, "np_id", conversion(sp_np_pair, "org_id", "NPO1797", "np_id"), "pref_name")
ggboxplot(num_cpds, add = "jitter") + scale_y_log10()

sp_info <- sp_info[num_cpds >= 2 & num_cpds <= 91, ]

# filter out species alone in its genus
num_fam <- table(sp_info$family_name)
summary(as.numeric(num_fam)); quantile(num_fam, probs = seq(0, 1, 0.05))
ggboxplot(as.numeric(num_fam), add = "jitter") + scale_y_log10()

num_gen <- table(sp_info$genus_name)
summary(as.numeric(num_gen)); quantile(num_gen, probs = seq(0, 1, 0.05))
ggboxplot(as.numeric(num_gen), add = "jitter") + scale_y_log10()

sp_info <- sp_info[duplicated(sp_info$genus_name), ]
dim(sp_info)

# chemical diversity comparison between to species
chem_diver <- function(sp1, sp2, visualize = FALSE) {
        nps1 <- conversion(sp_np_pair, "org_id", sp1, "np_id")
        nps2 <- conversion(sp_np_pair, "org_id", sp2, "np_id")
        
        smi1 <- conversion(np_prop, "np_id", nps1, "canonical_smiles"); cat(length(smi1), "\n")
        smi2 <- conversion(np_prop, "np_id", nps2, "canonical_smiles"); cat(length(smi2), "\n")
        
        fpmatrix1 <- smi2fp(smi1); l1 <- nrow(fpmatrix1)
        fpmatrix2 <- smi2fp(smi2); l2 <- nrow(fpmatrix2)
        if(l1 <= 1 | l2 <= 1 | is.null(l1) | is.null(l2)) return()
        
        fpjacdis <- fp2jacdis(rbind(fpmatrix1, fpmatrix2))
        in1 <- fpjacdis[1:l1, 1:l1]
        in1 <- max_sim(in1)
        in2 <- fpjacdis[l1 + 1:l2, l1 + 1:l2]
        in2 <- max_sim(in2)
        between <- fpjacdis[1:l1, l1 + 1:l2]
        between <- c(max_sim(between), max_sim(t(between)))

        if (visualize) {
                df_plot <- data.frame(Similarity = as.numeric(c(in1, in2, between)),
                                      Pair = c(rep("sp1-sp1", length(in1)),
                                               rep("sp2-sp2", length(in2)),
                                               rep("sp1-sp2", length(between))))
                ggboxplot(df_plot, x = "Pair", y = "Similarity", add = "jitter")
        } else {
                return(c(median(in1), median(in2), median(between)))
        }
}
        
n = 1
n = n + 1; chem_diver(sp_info$org_id[n], sp_info$org_id[nrow(sp_info) - n])
        
genus <- unique(sp_info$genus_name)
intra <- NULL
inter <- NULL
for(g in genus[1:10]) {
        ind <- (sp_info$genus_name == g)
        sp_count <- table(ind)[2]
        sp_in <- names(table(sp_info$org_id[ind]))
        sp_out <- names(table(sp_info$org_id[!ind]))
        for(n in 1:sp_count) {
             candidate <- sample(sp_in, size = 2)
             intra <- rbind(intra, c(chem_diver(candidate[1], candidate[2]), g))
             inter <- rbind(inter, c(chem_diver(candidate[1], sample(sp_out, size = 1)), g))
             inter <- rbind(inter, c(chem_diver(candidate[2], sample(sp_out, size = 1)), g))
        }
}

df_plot <- data.frame(Similarity = as.numeric(c(inter[, 3], intra[, 3])),
                      Sp = c(inter[, 4], intra[, 4]),
                      Type = c(rep("inter", nrow(inter)),
                               rep("intra", nrow(intra))))
ggboxplot(df_plot, x = "Sp", y = "Similarity", add = "jitter", color = "Type")


