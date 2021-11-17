#### BINF 6210 - Assignment 2 - Due Friday October 29, 2021 by 5 pm ----
# By Shalvi Chirmade
# Secondary author: David Jamieson
# Tertiary authors: Emily Maier, Eden Blanchard

### HOW DOES GEOGRAPHY RELATE TO DIVERSIFICATION OF THE GENUS CANIS?



#### 1- INTRODUCTION ----

#I have always been fascinated with dogs and wolves for as long as I can remember. As a child, our family never had pets and the thought of a dog was exhilarating. The variety of breeds was beyond my comprehension but all I wanted was to learn more. Little did my young mind know, all breeds of household dogs encompass not only one species but a single subspecies of Canis lupus, Canis lupus familiaris. There are 32 subspecies of Canis lupus found worldwide (Snyder, 1991); they differ based on genetic and phenotypic attributes such as size, color and geographical distribution (Snyder, 1991). This genus is one of the most common carnivore taxa found globally with at least one species found in every continent, with an exception of Antarctica (Padilla & Hilton, 2015). Due to their extensive geographical habitat, researchers have observed a wide diversity throughout the taxonomic group. Geographical locations include deserts, tropical forests, coastal areas and even snow-topped mountain peaks (Padilla & Hilton, 2015)! However, over the past few years, experts have noticed a change in habitat for a few species; for example, Canis latrans were uniquely found in western North America and are now seen all over the continent (Sillero-Zubiri et al., 2004). The study of this family is of great importance as several taxa within this genus are considered critically endangered and reaching extinction (Padilla & Hilton, 2015). Unfortunately, as we will see in my data analysis, there are only a handful of sequences found in NCBI for these particular threatened species. A higher amount of studies are being conducted on commonly-found species, Canis latrans, and very few on these endangered species. As seen in the literature, these species are difficult to locate and are found in isolated areas (Gotelli et al., 1994). There are currently various conservation plans carried out by IUCN which includes The Canid Action Plan encompassing the whole family of canids (Sillero-Zubiri et al., 2004). This plan was authored by 88 experts in 2004 and has seen a great improvement over the past decade.

#The Canis genus has an expansive geographic distribution, but relatively low speciation. Will clustering by genetic sequences yield trends of varied allopatric evolution based on geographic isolation, or will an overall sympatric trend be present among the global population of this genus?




#### 2- DATA ACQUISITION ----


##Install and/or load packages required for this script.

#Packages from CRAN

#install.packages("ape")
library(ape)
#install.packages("maps")
library(maps) #Must be loaded before mapdata.
#install.packages("mapdata")
library(mapdata)
#install.packages("phangorn")
library(phangorn)
#install.packages("phytools")
library(phytools)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("rentrez")
library(rentrez)
#install.packages("rgbif")
library(rgbif)
#install.packages("RSQLite")
library(RSQLite)
#install.packages("stringi")
library(stringi)
#install.packages("seqinr")
library(seqinr)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("viridis")
library(viridis)

#Packages from Bioconductor

#Install Bioconductor installation system, if needed.
#install.packages("BiocManager")

#Can install in one command, using the next line or install individually if needed.
#BiocManager::install(c("Biostrings", "muscle", "DECIPHER"))

#library(BiocManager)
library(Biostrings)
#BiocManager::install("muscle")
library(muscle)
#BiocManager::install("DECIPHER")
library(DECIPHER)



#Set my working directory
getwd()
#setwd("/Users/shalvichirmade/Documents/MBinf/BINF 6210/Assignments/2")



##Extract data needed from NCBI using rentrez package


#The genus I have chosen is, Canis. The gene I have chosen is cytochrome B; the literature shows this particular gene from mitochondrial DNA to be of importance in the study of canid phylogeny (Gaubert et al., 2012; Okumura et al., 1998 ;Rumman et al, 2019).


#I am going to use the nucleotide database, nuccore, to obtain DNA sequence data for this study. I need to look at the fields I can search for under this database.
Fields_Nuccore <- entrez_db_searchable("nuccore")
Fields_Nuccore

#Need to use ORGN for the taxonomic name, GENE for the gene I am searching for and SLEN to filter for the proper sequence length.
#Searching the database to see how many hits I get 
cytb_search <- entrez_search(db = "nuccore", term = "Canis[ORGN] AND (Cyt B[Gene] OR CytB[Gene])")
cytb_search

#We have 2876 hits and only 20 IDs because I did not set the retmax. The number of hits is based on my search on October 17, 2021; this number has increased now.

#According to the UniProt database, the average sequence length of cytB is 379 bp; here is the link for your reference. https://www.uniprot.org/uniprot/?query=cyb+canis&sort=score
#I will set the SLEN length accordingly to accurately capture the correct sequences.
#In my previous searches, I received a hit count of 567 which was too large for entrez_fetch. I will create a smaller dataset by further restricting sequence length.

cytb_length_search <- entrez_search(db = "nuccore", term = "Canis[ORGN] AND (Cyt B[Gene] OR CytB[Gene]) AND 350:400[SLEN]")
cytb_length_search

#Now I have 449 hits. How many IDs have been found?
summary(cytb_length_search) #20


#I had previously limited my retmax search criteria to prevent receiving a fail error from extracting too many records. After speaking to Jacqueline, she helped me extract my data using the web history tool which allowed for a larger sequence accumulation. Here are the lines I used for entrez_fetch; most lines have been deleted but I left a few to visualize the difference.


#Rerun the database search using the retmax criteria to receive all sequence data possible without receiving a fail error.
#cytb_length_search <- entrez_search(db = "nuccore", term = "Canis[ORGN] AND (Cyt B[Gene] OR CytB[Gene]) AND 350:400[SLEN]", retmax = 300)
#cytb_length_search
#Now the number of IDs should be 300.
#summary(cytb_length_search)
#Can see the change in the number of IDs.

##Next, I need to acquire the nucleotide sequences of my search criteria. This is done by using the IDs and creating a fasta file.
#cytb_fetch <- entrez_fetch(db = "nuccore", id = cytb_length_search$ids, rettype = "fasta")

#Inspect the data created.
#summary(cytb_fetch)
#head(cytb_fetch)
#Can see that all the sequences are in one long character form. We can separate out each sequence by using the "/n" deliminator.

#Write the file to my current working directory.
#write(cytb_fetch, "cytb_fetch.fasta", sep = "\n")
#Checked file in TextEdit to make sure all sequences fit the requirements I specified. All sequences have been separated by the "new line" deliminator.

#Read the fasta file into DNAStringSet to be able to further analyze the nucleotide sequences.
#cytb_stringSet <- readDNAStringSet("cytb_fetch.fasta")

#view(cytb_stringSet)
#We can see that the row names are the title of each nucleotide sequence and the sequence itself is in a column called "x".

#Interpolate the data into a data frame (cannot use tibble to later convert to DNAStringSet).
#dfCytB <- data.frame(cytB_Title = names(cytb_stringSet), cytB_Sequence = paste(cytb_stringSet))
#View(dfCytB)



##Here is the code I used to obtain the data I currently have. I used Jacqueline's functions to extract a larger data set without having to restrict retmax.

#During the TA session with Jacqueline, she requested for us to send our data files along with our R script. I have commented the lines of code used to obtain this data and left the code used to extract the data file for your convenience.

#Source Jacqueline's functions file.
#source("Entrez_Functions.R")

#As the function states, this line of code allows the user to copy their search output into fasta files with 100 sequences in each.
#FetchFastaFiles(searchTerm = "Canis[ORGN] AND (Cyt B[Gene] OR CytB[Gene]) AND 300:400[SLEN]", seqsPerFile = 100, fastaFileName = "Canis_cytB")

#Merging the sequences into a data frame. I am using generic names while extracting data to allow for reproducible analysis.
#dfNCBI <- MergeFastaFiles(filePattern = "Canis_cytB*")


#Write the file into my working directory so I don't have to download the data every time I run this script.
#write_tsv(dfNCBI, "CanisCytBDataNCBI.tsv")

#Download the data from my working directory into the working environment.
dfNCBI <- read_tsv(file = "CanisCytBDataNCBI.tsv")


#Remove the variables we do not require.
rm(Fields_Nuccore,cytb_search, cytb_length_search)



##Extract data needed from GBIF using the rgbif package.


#As before, I have left the lines of code used to extract the data in comment form. 

#Search for occurrences of the genus Canis in the GBIF database. Using occ_data() instead of occ_search() as the extra columns are not needed for this analysis. I am using generic names while extracting data to allow for reproducible analysis.
#Taxon_GBIF <- occ_data(scientificName = "Canis")

#There's a list list of metadata and data. We need the data tibble.
#dfGBIF <- Taxon_GBIF$data
#View(dfGBIF)


#Write the file into my working directory so I don't have to download the data every time I run this script.
#write_tsv(dfGBIF, "CanisDataGBIF.tsv")

#Download the data from my working directory into the working environment.
dfGBIF <- read_tsv(file = "CanisDataGBIF.tsv")

#Analysis and filtering will be done in the next section.




#### 3- NCBI DATA FILTERATION AND QUALITY CONTROL ----


##Analyzing and filtering the data obtained from NCBI.

##This data is from the web history objects.

View(dfNCBI)
#I have 469 observations

#Evaluate the data merged. What class is it?
class(dfNCBI)

#Create a new column to extract the scientific name of the organism for each sequence. Some rows have a subspecies name in the fourth word.
dfNCBI$Species_Name <- word(dfNCBI$Title, 2L, 4L)
View(dfNCBI)

#Rearrange the columns to allow for a neater view.
dfNCBI <- dfNCBI[, c("Title", "Species_Name", "Sequence")]
View(dfNCBI)


#How many species are in my data frame?
length(unique(dfNCBI$Species_Name)) #19

#What are they?
unique(dfNCBI$Species_Name)


#Create a data frame with the unfamiliar Species_Name, to investigate if they should be removed.
dfUnfamiliar <- dfNCBI %>%
  filter(Species_Name == "Canis environmental sample")
View(dfUnfamiliar)

#Decide to removed the environmental sample as this will not be useful for our further phylogeny and geographical analysis.
dfNCBI_CanisCytB <- dfNCBI %>%
  filter(!Species_Name == "Canis environmental sample")

#See if this worked.
length(unique(dfNCBI_CanisCytB$Species_Name)) #18
unique(dfNCBI_CanisCytB$Species_Name)
#It has been removed.

#Remove the third word in Species_Name if it is not part of the taxonomic identity.
dfNCBI_CanisCytB$Species_Name <- str_remove(dfNCBI_CanisCytB$Species_Name, " isolate| haplotype| cytochrome| mitochondrion| voucher")

#dfNCBI_CanisCytB <- dfNCBI_CanisCytB %>%
#  mutate(Species_Name = str_remove(Species_Name, " isolate")) %>%
#  mutate(Species_Name = str_remove(Species_Name, " haplotype")) %>%
#  mutate(Species_Name = str_remove(Species_Name, " cytochrome")) %>%
#  mutate(Species_Name = str_remove(Species_Name, " mitochondrion")) %>%
#  mutate(Species_Name = str_remove(Species_Name, " voucher"))
#Added a space before the word so it takes away the end space in the column.

#Find a way to do this all in one line. The next few lines show what doesn't work.

#dfTest <- dfNCBI_CanisCytB %>%
#  mutate(Species_Name = str_remove(Species_Name, " isolate" | " haplotype" | " cytochrome" | " mitochondrion" | " voucher"))

#dfTest <- dfNCBI_CanisCytB %>%
# mutate(Species_Name = str_remove_all(Species_Name, c("isolate", "haplotype", "cytochrome", "mitochondrion", "partial", "voucher")))

#Create a character vector with the words I need to remove.
#WordsToRemove <- c("isolate", "haplotype", "cytochrome", "mitochondrion", "partial", "voucher")

#dfTest <- dfNCBI_CanisCytB %>%
#mutate(Species_Name = str_remove(Species_Name, WordsToRemove))
#Removes some but not all..



#Check to see if it has worked.
length(unique(dfNCBI_CanisCytB$Species_Name)) #12
sort(unique(dfNCBI_CanisCytB$Species_Name))

#Take a closer look at the data frame.
dim(dfNCBI_CanisCytB)
class(dfNCBI_CanisCytB)

#Some entries use the subspecies only and do not specify species. I will leave it in for now and decide to add the species name later if required.


#Count the number of entries per species. 
dfSpecies <- (dfNCBI_CanisCytB) %>%
  count(Species_Name)
View(dfSpecies)

#There are 411 sequences of Canis latrans. After aligning all these sequences, I noticed that they were very similar. I checked a few of the accession numbers on the web browser to make sure they match to Canis latrans. To have a smaller data set and to allow for better alignment, I will be randomly selecting 20 records from this species to better align my sequences for the final analysis. 

#Decided to write this as a function. This adds predictability to the code allowing the user to easily change the species or sample size number. Could also use this function for any data frames with a column named "Species_Name".
sample.species <- function(df, species, n){
  #Filter only the species we're interested in.
  dfsample <- df %>%
    filter(Species_Name == species)
  #Sample for the number of records we want.
  dfspec <- sample_n(tbl = dfsample, size = n)
  #Create a data frame with all other records
  dfrest <- df %>%
    filter(!Species_Name == species)
  #Rbind the two data frames together
  df <- rbind(dfspec, dfrest)
  #Return the new data frame.
  return(df)
  
}
set.seed(2021)
#Limit the number of records to 20 for the species Canis latrans
dfNCBI_CanisCytB <- sample.species(df = dfNCBI_CanisCytB, species = "Canis latrans", n = 20)


#First create a data frame with just Canis latrans and then randomly select 20 records. Setting seed so the records will be reproducible.
#dfCanislatrans <- dfNCBI_CanisCytB %>%
#  filter(Species_Name == "Canis latrans")
#
#set.seed(2021)
#dfCanislatrans <- dfCanislatrans %>%
#  sample_n(20)
#
#dim(dfCanislatrans)
#unique(dfCanislatrans$Species_Name)
#
##Next, create a second data frame with all the other species records.
#dfOthers <- dfNCBI_CanisCytB %>%
#  filter(!Species_Name == "Canis latrans")
#View(dfOthers)

#Now merge these two data frames and overwrite previously created dfNCBI_CanisCytB.
dfNCBI_CanisCytB <- base::rbind(dfCanislatrans, dfOthers)

#Check the data frame created.
dim(dfNCBI_CanisCytB)
class(dfNCBI_CanisCytB)

#Checking the number of records per species again.
dfSpecies <- (dfNCBI_CanisCytB) %>%
  count(Species_Name)
View(dfSpecies)
#Can now see that the number of records for Canis larans is 20.
#I re-ran this set of code to make sure the set seed was working and was choosing the same records.


#Check sequence lengths of data downloaded from NCBI
summary(str_count(dfNCBI_CanisCytB$Sequence))

#Create a histogram to showcase the sequence length seen in the data. Credit for italics in the title goes to Nykole Crevits and rischan from stackoverflow https://stackoverflow.com/questions/23969726/italics-and-normal-text-in-a-main-plot-title .

#Build a color palette for the histogram; just for fun. I'm using the same color palette as my last assignment, it is colorblind-friendly. The colors have no meaning in this graph, it just looked nicer than plain grey.
colors <- brewer.pal(8, "BrBG")

hist(str_count(dfNCBI_CanisCytB$Sequence), 
     main = substitute(paste(italic('Canis'), " Sequence Distribution")), 
     xlab = "Sequence Length (bp)",
     col = colors)
#As previously researched, most of the sequences are in the range found on UniProt.


#Check to see if there are any N's or gaps in the sequences.
mean(str_count(string = dfNCBI_CanisCytB$Sequence, pattern = "N")) #0.05194805
mean(str_count(string = dfNCBI_CanisCytB$Sequence, pattern = "-")) #0
#There are no gaps and barely any N's.

#As my sequence data has been downloaded from NCBI, there is an insignificant amount of N's in the entirety of this data set. Due to this, sequence manipulation is not needed.

#Remove items from the workspace that are no longer needed.
rm(dfCanislatrans, dfOthers, dfSpecies, dfUnfamiliar, colors, dfNCBI)




#### 4- GBIF DATA FILTERATION AND QUALITY CONTROL ----


#Interpret the available columns and data provided.
names(dfGBIF)
length(dfGBIF)

#Subset the data frame created for columns needed for analyzing the data.
dfGBIF_Subset <- dfGBIF [ , c("gbifID", "kingdom", "phylum", "class", "order", "family", "genus", "species", "country", "countryCode", "decimalLatitude", "decimalLongitude", "taxonomicStatus", "institutionCode")]
View(dfGBIF_Subset)


#Check the data frame.
dim(dfGBIF_Subset)
class(dfGBIF_Subset)

#Check to see if all records have the same kingdom, phylum, class, order, family and genus. Comment on anomalies.
apply(X = dfGBIF_Subset[, c("kingdom", "phylum", "class", "order", "family", "genus")], MARGIN = 2, FUN = unique)

#unique(dfGBIF_Subset$kingdom)
#unique(dfGBIF_Subset$phylum)
#unique(dfGBIF_Subset$class)
#unique(dfGBIF_Subset$order)
#unique(dfGBIF_Subset$family)
#unique(dfGBIF_Subset$genus)
#There are no anomalies or even NAs!


#What are the unique species in this database?
length(unique(dfGBIF_Subset$species)) #5
unique(dfGBIF_Subset$species)
#There is a NA value. I want to view the record for this one and analyse if it can be discarded.

dfSpeciesNA <- dfGBIF_Subset %>%
  filter(is.na(species))
View(dfSpeciesNA)
#I will delete these rows as they are not needed for future analysis.
rm(dfSpeciesNA)


#The taxonomic status column tells us the validity of the taxonomic classification https://rs.gbif.org/vocabulary/gbif/taxonomic_status.xml#.
unique(dfGBIF_Subset$taxonomicStatus)
#There is ACCEPTED and SYNONYM.

#How many records are there of each?
dfGBIF_Subset %>%
  group_by(taxonomicStatus) %>%
  count()
#There are 460 ACCEPTED and 40 SYNONYM.
#SYNONYM records have an unclear nomenclature status; I will remove these records are they are not needed for our future analysis.


#What countries are these records from? 
length(unique(dfGBIF_Subset$country)) #24
unique(dfGBIF_Subset$country)

#As we are doing a geographical analysis, unspecified countries are not useful to us. Check to see if they have other geographical data, like latitude and longitude. If so, I will need to specify the country if possible.
dfNACountry <- dfGBIF_Subset %>%
  filter(is.na(country))
View(dfNACountry)
#There are two records with NA in country and both records do not have other geographical data in either latitude, longitude or country code. Will be deleted.
rm(dfNACountry)


##Filter the data frame with only the records needed for further analysis.
dfGBIF_Subset <- dfGBIF_Subset %>%
  filter(!is.na(species)) %>%
  filter(!is.na(country)) %>%
  filter(taxonomicStatus == "ACCEPTED")

#Check newly created data frame.
dim(dfGBIF_Subset)


#How many records are they for each species?

dfGBIF_Subset %>%
  group_by(species) %>%
  count()

#How many records are they for each country?

dfGBIF_Subset %>%
  group_by(country) %>%
  count()


#How many records of unique species in each country?

dfSpeciesbyCountry <- dfGBIF_Subset %>%
  group_by(country, species) %>%
  count(species)
View(dfSpeciesbyCountry)
#Can see that no country has records of all four species. USA seems to have done quite a lot of research on canids in their geographical regions. The geophylogeny image at the end of the script will allow us to see if they are spread around the country or only from an isolated region.

#Remove data frames that are not necessary for further analysis.
rm(dfGBIF, dfSpeciesbyCountry)



#### 5- SEQUENCE ALIGNMNET ----


#Convert tibble to a data frame for further analysis.
dfNCBI_CanisCytB <- as.data.frame(dfNCBI_CanisCytB)
class(dfNCBI_CanisCytB)


#For the remainder of our analysis, we want to name each of these sequences to be able to re-identify the record they belong to. I will be creating a new column with just the Accession numbers for each record and then assigning this to the sequence.
dfNCBI_CanisCytB$Accession_Number <- word(dfNCBI_CanisCytB$Title, 1L)

#Check to make sure they are all unique
length(unique(dfNCBI_CanisCytB$Accession_Number)) == length(dfNCBI_CanisCytB$Title)


#Here I am making a new subset data frame for my alternative alignment strategy. I have to do this step before converting the sequence column to a DNAStringSet. I will give a more detailed explanation after you view the first alignment. This second alignment will come afterwards.

#Create a vector containing the records I want to remove for the alternative alignment.
accessionremove <- c("AF028165.1", "AF028168.1", "AF028167.1", "AF028166.1", "AF028162.1", "AF028160.1", "JX849653.1", "JX849648.1", "L29416.1", "L29415.1")

#Create the subset data frame.
dfAlignment2 <- dfNCBI_CanisCytB %>%
  filter(!Accession_Number %in% accessionremove)
View(dfAlignment2)


##For the sequences to be aligned using the Biostrings packages, the column must be converted to a DNAStringSet class. I am going to create a new column to carry this out. 
dfNCBI_CanisCytB$SequenceDSS <- DNAStringSet(dfNCBI_CanisCytB$Sequence)

#Check to see if it has worked.
class(dfNCBI_CanisCytB$SequenceDSS)

#Reorder the data frame for easier viewing.
dfNCBI_CanisCytB <- dfNCBI_CanisCytB[, c("Title", "Accession_Number" ,"Species_Name", "Sequence", "SequenceDSS")]
View(dfNCBI_CanisCytB)


#Assign accession number to each sequence so the sequence will have a unique name associated to it. I used this when I was analyzing my sequence alignment and trying to find the incongruous records. It is commented out now as in my later analysis, the name of my species is more important for answering my question than the accession number. I had a talk with Sally about this and she agreed that using species names for my final analysis was appropriate for my question.
#names(dfNCBI_CanisCytB$SequenceDSS) <- dfNCBI_CanisCytB$Accession_Number


#I'm going to assign the corresponding species name to each sequence; due to the nature of my question, the species name is most important. As the NCBI data does not have geographical data associated with each sequence, there is no reason, at this moment, for me to categorize each unique sequence to its accession number.
names(dfNCBI_CanisCytB$SequenceDSS) <- dfNCBI_CanisCytB$Species_Name

#Check to see if it worked.
names(dfNCBI_CanisCytB$SequenceDSS)

#Now that we have our sequences ready for alignment, lets now view them neatly in our browser.
BrowseSeqs(dfNCBI_CanisCytB$SequenceDSS)
#We can see in this visualization, that some sequences are quite dissimilar.


#I am going to do my alignment based on default methods as this gene is highly conserved in my taxa. According to the results, I will edit my methods accordingly to create a better alignment. After multiple tries of the alignment, I found this combination of arguments to produce the "best" alignment. I am using the default setting for maxiters as my number of sequences is low and have set gap penalty high due to some species having variable sequences in comparison to the others, which causes a lot of gaps being formed without a high penalty.
CytB_alignment <- DNAStringSet(muscle::muscle(dfNCBI_CanisCytB$SequenceDSS, gapopen = -1000), use.names = T)

#Check class of alignment made.
class(CytB_alignment)

#View alignment in browser.
BrowseSeqs(CytB_alignment)
#Can see that sequences from different species align noticeably..

#Write the alignment file to my working directory to analyze further in MEGA. Commented so it will not re-write the file.
#writeXStringSet(CytB_alignment, file = "Canis_Cytb_Alignment.fas", format = "fasta")

#Time to evaluate our alignment.
CytB_alignment[1] #shows the first sequence
length(CytB_alignment[[1]])
#After the alignment, this particular sequence went from 369 sequences to 566. As seen in the visual alignment, there were multiple gaps added at the beginning, middle and end of the sequence to allow for species gene variability.

#Let's look at the mean number of gaps in the alignment.
mean(unlist(lapply(X = CytB_alignment, FUN = str_count, pattern = "-"))) #203.026
#We can see that the alignment has caused an average of 203 gaps being added to each sequence. 

#I also viewed this sequence alignment in MEGA; I can see a lot of variability in about ten sequences out of the total. Even after translation, a huge variety of amino acids are seen for those particular sequences.



#After analyzing this alignment result with Sally, she helped me realize that this particular alignment cannot be used for my further analysis. There are a few sequences in this data set that are misaligning. Sally's recommendations were to look at these specific accession numbers and see either if they are from a different section from the gene, the reverse compliment, or just poor sequences. Here is a list of the accession numbers I analysed. 

#AF028165.1 - Canis lupus - Study A
#AF028168.1 - Canis simensis - Study A
#AF028167.1 - Canis mesomelas elongae - Study A
#AF028166.1 - Canis mesomelas elongae - Study A
#AF028162.1 - Canis aureus - Study A
#AF028160.1 - Canis adustus - Study A
#JX849653.1 - Canis lupus familiaris - isolate Taiwan12 - Study B
#JX849648.1 - Canis lupus familiaris - isolate Taiwan23 - Study B
#L29416.1 - Canis simensis - Study C
#L29415.1 - Canis simensis - Study C

#Study A - Wayne et al., 1997 - As seen in the paper and the GenBank information, most of these sequences are from separate regions of the cytochrome B gene. The primers used for PCR analysis encompassed 2001 bp region of mitochondrial DNA. This region includes three protein coding regions: cytochrome B, COI and COII.
#Study B - This study says "unpublished" so I was unable to find more information about these sequences.
#Study C - Gottelli et al., 1994 - This study was similar in their sequence extraction as Study A; they used primers to encompass the same 2001 bp region of mtDNA. Both papers actually used the same primer set for cytochrome B! This could be the reason why these sequences are similar to each other but dissimilar to the rest!

#After analyzing each paper and researching deeper about the sequences that do align well, I have decided to remove these ten sequences and perform another sequence alignment with the remainder. The dendrogram that was created using this first alignment had a distance scale of 1.2! I hope to create a better alignment; let's see!



##Now, let's use the updated data frame (line 453) which excludes these misaligned sequences. I'm nervous and excited to see what becomes of it and how it changes the dendrogram! I'm using the same order of steps as before, so I will only provide brief comments for every line of code.

#Make new column containing sequence data in DNASTringSet class.
dfAlignment2$SequenceDSS <- DNAStringSet(dfAlignment2$Sequence)

#Check to see if it has worked.
class(dfAlignment2$SequenceDSS)

#Reorder the data frame for easier viewing.
dfAlignment2 <- dfAlignment2[, c("Title", "Accession_Number" ,"Species_Name", "Sequence", "SequenceDSS")]
View(dfAlignment2)

#Use species names to name each sequence.
names(dfAlignment2$SequenceDSS) <- dfAlignment2$Species_Name
names(dfAlignment2$SequenceDSS)

#Carry out alignment.
CytB_alignment2 <- DNAStringSet(muscle::muscle(dfAlignment2$SequenceDSS, gapopen = -100), use.names = T)

#View alignment in browser.
BrowseSeqs(CytB_alignment2)

#I DID IT! I FINALLY DID IT! It has taken me so long to figure out this alignment issue and I'm so happy it worked!! Before committing to this alignment, I carried out clustering to make sure it was creating a correct dendrogram as well! The first alignment caused my dendrogram scale to go to 1.2, so an evolutionary distance of 120%, and this alternate alignment has brought my scale bar down to 0.15! I will explain further in the next section.

#We can still see a few anomalies in this alignment but as we will notice later on, they belong to species that form distinct clusters.

#Remove first alignment data.
rm(CytB_alignment, accessionremove)



#### 6- CLUSTERING AND PHYLOGENY ----


#The model I have chosen is based on my research from the literature. I found a few papers that use K80 as their model for the Canis and cytochrome B (Gaubert et al., 2012; Okumura et al., 1998). K80, also known as K2P, is the Kimura 2-parameter model (Kimura, 1980) which distinguishes the two differences in nucleotide change, transitions or transversions.

#I first used modelTest() to determine the best fit model for my sequences; the lowest AIC score was allotted to GTR+G+I. I ran this function when I was initially using my whole data set, about 900 sequences. After I made the decision to cut down my number of sequences (first alignment), modelTest() then returned GTR+G as the best suitable model. I found one paper that used GTR+G+I for Canis lupus familiaris (Guo et al., 2015) only because "Model Selection" on MEGA suggested it to the authors. The authors were also conducting the analysis using other genera along with Canis. As the rest of the papers I examined used K80 for this genus, I continued to use the same model as it is available in dist.dna(). I hope this is something I can address with my grouep members in Assignment 3 and compare the difference in alignment using a more complicated model.
#CytB_alignment_phyDat <- as.phyDat(CytB_alignment_DNAbin)
#model <- modelTest(CytB_alignment_phyDat)


#First convert the alignment into DNAbin format. 
CytB_alignment2_DNAbin <- as.DNAbin(CytB_alignment2)

#Create a distance matrix using the chosen model. 
CytB_distMatrix <- dist.dna(CytB_alignment2_DNAbin, model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE) #pairwise deletion - missing data will be ignored in a pairwise fashion, can lose a lot of information if you use complete deletions

class(CytB_distMatrix) #matrix array
summary(as.vector(CytB_distMatrix))

#I had to thoroughly analyze my distance matrix when using my first alignment as my highest number was 2.66! I used the "raw" model in dist.dna() to help understand further but I didn't get any new information. Now that I'm using my second alignment, the highest value is 0.36 which is an expected number when comparing evolutionary similar sequences.


#Finally create the dendrogram based on a selected clustering method and threshold. I have chosen this method as Neighbor-Joining trees as it considers the variation of evolutionary rates while creating the phylogenetic tree. As this gene is supposed to be well conserved and showcase accurate evolutionary change in the genus Canis, I believe a Neighbor-Joining tree will create the most accurate phylogenetic tree for this data.
par(mar = c(8,5,2,2)) #The bottom part of the dendrogram was cut off, had to change margins.
Cytb_clusters <- IdClusters(CytB_distMatrix,
                           method = "NJ",
                           cutoff = 0.1,
                           showPlot = TRUE,
                           type = "both",
                           verbose = TRUE)
par(mar = c(2,2,2,2)) #Setting back to default.


#Comment on error
#Duplicated labels in myDistMatrix appended with index.
#This error is created because I have multiple records with the same name. The dendrogram adds numbers to each record in the image but does not add the numbers to the actual record name.

#Here we can see that the dendrogram has definitely proven our assumption right during alignment! Canis lupus and its subspecies form distinct clusters. So, hypothesizing that these sequences are evolutionarily diverse based on the other species, such as Canis latrans, was correct. The result of this phylogenetic tree is also seen in the literature. Figure 5b in Rumman et al. shows similar phylogenetic evolution for this genus. I will talk more about this in the discussion. 



#How many number of clusters do I have?
length(unique(unlist(Cytb_clusters[[1]][1]))) #7

class(Cytb_clusters)
Cytb_clusters[1] #Tells you which cluster each sequence belongs to.
class(Cytb_clusters[[1]])
View(Cytb_clusters[[1]])
#We can see that some clusters only carry a single sequence. This can mean either that these sequences are associated with species that are very different than the others or it can mean that the sequence quality was poor and hence did not fall under another related species cluster. However, as seen in the data frame, these sequences belong to the evolutionarily diverse species; this is accurate based on what was seen earlier. Will comment on this during the discussion section.




#### 7- GEOPHYLOGENY ----

#Earlier in this script, during data wrangling, we realized that NCBI provided us with 12 different species names while GBIF only provided four. After further investigation, I found that some records from NCBI used the subspecies name as the species, which creates some confusion during analysis. As the data from GBIF does not provide any subspecies names (other than Canis lycaon), I will have to manipulate the NCBI species name to account for the four species found on GBIF.

#Canis aureus, Canis latrans and Canis lycaon have matching records in NCBI; they will not be changed.

#Canis lupus encompasses six species names from NCBI: Canis familiaris, Canis lupus, Canis lupus chanco, Canis lupus familiaris, Canis lupus pallipes, and Canis rufus.

#Three species from NCBI are not accounted for in GBIF and cannot be used for this analysis: Canis adustus, Canis mesmolas elongae and Canis simensis.

#For assignment 3, I would like to pull records from BOLD as well and incorporate new geographical data into this matrix. I hope these species that cannot be used today will be able to be incorporated later on.


#I will start by making a crude phylogeny tree based on the dendrogram created above. This has to be done due to the varied species names found from both databases (I checked my reasoning with Sally and she agreed that this is the best path to take right now). I am just using the species name and leaving out the genus name as read.tree conjuncted the names together; I find it easier to read the tree with just the species.
CanisNames <- "(lupus, (aureus, (lycaon, latrans)));"

class(CanisNames) #character


#Create the crude tree using these names found in GBIF.
CanisTree <- read.tree(text = CanisNames)

#Check the class of this newly formed tree.
class(CanisTree) #phylo

#Plot the tree and see if CanisNames needs to be changed to form the correct tree. Commented out so it does not create an additional image for my assignment.
#plot(CanisTree)


#View the phylo object.
CanisTree
#Tree has 4 tips (one for every species in GBIF), and 3 internal nodes. The tip labels correspond to the species names and the tree is rooted with no branch lengths.

#I will create branch lengths to approximately match the dendrogram created earlier.
CanisTree$edge.lengths <- c(4, 1, 1, 1)

#Make sure the edge lengths have been incorporated into the phylo object.
CanisTree[[4]]


#Now I need to create the matrix of geographical data. I have 458 records in my GBIF data frame. I will see how crowded the geophylogeny map becomes; if it is hard to read, I will cut down the number of records by limiting a few per country.

#Create a small data frame with just the columns required for this matrix.
dfGBIF <- dfGBIF_Subset %>%
  select(species, decimalLatitude, decimalLongitude)

#Analyze and make sure there are no NAs.
unique(dfGBIF$species) #"Canis lupus"   "Canis latrans" "Canis aureus"  "Canis lycaon" 
sum(is.na(dfGBIF$decimalLatitude)) #0
sum(is.na(dfGBIF$decimalLongitude)) #0

#When creating matrix for GBIF, remove genus name. This way, the species name matches the tip names from phylo object, CanisTree.
dfGBIF$species <- word(dfGBIF$species, 2L)

#Make sure it worked.
unique(dfGBIF$species)

#Make the species column into the rownames.
#rownames(dfGBIF)
#rownames(dfGBIF) <- dfGBIF$species 
#The error message states, "duplicate 'row.names' are not allowed". After researching this error, it is found that unfortunately, you cannot have duplicate rownames in a data frame, however, you can in a matrix. Let me try doing this after converting the data frame to a matrix.


#Convert the data frame into a matrix
matrixGBIF <- as.matrix(dfGBIF)

#Check to see if it worked.
class(matrixGBIF)
View(matrixGBIF)


#Remove the species column as it is no longer needed.
matrixGBIF <- matrixGBIF[ , -1]

#Check if it worked.
View(matrixGBIF)


#Convert the character matrix into a numeric matrix; a numerical matrix is needed for the use of the function phylo.to.map(). Credit for using the apply() function goes to Jacqueline May. I had forgotten to realize that a matrix can only hold one type of vector. When I was just using as.numeric(), the matrix went from having a dimension of 1:458, 1:2 to 1:966. After talking to Jacqueline and reading the documentation of the function, I now know that apply() can be used for matrix just like laaply() can be used for a data frame.
matrixGBIF <- apply(X = matrixGBIF, MARGIN = 2, FUN = as.numeric) 

#Convert rownames into the species name for each record.
rownames(matrixGBIF) <- dfGBIF$species
#https://stackoverflow.com/questions/31707654/r-set-duplicate-row-names-to-a-numeric-data-frame

#Check if it worked.
rownames(matrixGBIF)


#As mentioned in the sample script, I have to make my CanisTree ultrametric. This means that the tree to be used has to be rooted and have equal distance from all terminal nodes to the root. Setting the edge lengths earlier was redundant.
CanisUltra <- phytools::force.ultrametric(CanisTree)

#Make sure the tree is rooted and bifurcating.
ape::is.rooted(CanisUltra) #TRUE
ape::is.binary(CanisUltra) #TRUE

#Set the color scale to use for our geophylogeny tree. Using the viridis package for a nice representation of the different species.
CanisColors <- setNames(sample(plasma(n = Ntip(CanisUltra))),
                 CanisTree$tip.label)

#Check object made.
class(CanisColors)
CanisColors
#Associated a color to each species. Can change this later on if the colors are not well distinguished.

#Create a geophylogeny tree.
CanisPhyloMap <- phytools::phylo.to.map(tree = CanisUltra, coords = matrixGBIF, rotate = TRUE, type = "phylogram", fsize = 0.7, plot = F)


#Plotting the map using the associated colors.

phytools::plot.phylo.to.map(x = CanisPhyloMap, rotate = TRUE, type = "phylogram", colors = CanisColors, from.tip = T, lwd = 2, psize = 2, ftype = "b")
#For some reason, the full phylogram does not show the first time you run this line; have to run it again for full coverage. I can't find an explanation as to why this happens.


#I tried to remove the labels and create a legend instead but the legend is cut off when the image is zoomed. Also had to create a new character vector for the species names because in phylo.to.map, the species names on the lowest branches reverse in order so the legend associates the wrong color with the species. I hope this is also something I can work on for Assignment 3.
#labels <- c("Canis lupus", "Canis aureus", "Canis Lycaon", "Canis latrans")
#with(phytools::plot.phylo.to.map(x = CanisPhyloMap, rotate = TRUE, type = "phylogram", colors = CanisColors, from.tip = F, lwd = 2, psize = 2, ftype = "off"), 
#     legend(x = -150, y = 200, legend = labels, fill = CanisColors))





#### 8- RESULTS AND DISCUSSION


#As we have seen during this assignment, there have been numerous studies conducted on commonly-found Canis species and a low amount carried out on remote and endangered species. This was seen in both NCBI and GBIF data; NCBI had 411 records on Canis latrans, a common species found in North America, and GBIF had 414 records for the same species. The remaining records in both databases were minimal per species in comparison. I hope with the advancement of bioinformatics and more open source data, this genus will be studied at a higher level, especially in respect to IUCN's conservation plans (Sillero-Zubiri et al., 2004). With the records currently available, let me continue on with my analysis. As seen in the literature, cytochrome b is a common gene studied in respect to Canis evolution (Gaubert et al., 2012; Okumura et al., 1998 ;Rumman et al, 2019). According to these studies and the UniProt database, I chose a sequence length between 300 and 400bp. Figure 1 shows the average base pair length to be around 369 which closely relates to the data found in UniProt. During the process of narrowing down the optimum alignment, I noticed a clear distinction in a few species sequences in comparison to the others. Canis mesmomelas and Canis simensis, for example, were noticeably different in comparison to Canis latrans or Canis rufus. As mentioned during this section of the analysis, a few records were removed from the final alignment due to their sequences belonging to varied regions of the cytochrome b gene. This was confirmed by reading their respective papers and noticing the primers used for sequence amplification (Gottelli et al., 1994; Wayne et al., 1997). 

#Now, let's get to analyzing the results and answering my initial question! Figure 2 shows the dendrogram created using the aligned sequences. Here we can see that the clustering result helps answer my speculation during alignment; Canis latrans and Canis rufus were analogous in sequence similarity and show close clustering on the tree. Canis mesmomelas and Canis simensis, the diverse sequences, show separate clustering and higher evolutionary distance from the other species. Canis lupus, on the other hand, shows a strong evolutionary distance between other species as well as from itself! As I mentioned in the introduction, Canis lupus encompasses a wide range of genetically diverse canids (Snyder, 1991) which is shown in the cluster distribution in Figure 2 as well as their wide geographical distribution in Figure 3. Figure 5b from Rumman et al., also shows evolutionary distinction of Canis lupus in regards to the other species and shows a close relation between Canis latrans and Canis rufus. The geophylogeny image, Figure 3, shows the geographical distribution of four species in respect to the data available on GBIF. In comparison to Figure 2, Canis latrans and Canis lycaon are part of the same clade and are both found in North America; this shows sympatric speciation. Allpatric speciation is seen in respect to Canis lupus from their immensely diverse phylogeny and their distribution across the world (Sillero-Zubiri et al., 2004). Canis aureus, in Figure 3, shows geographical distribution in the northern parts of Africa and south Asia; this is confirmed by the research conducted by Siller-Zubiri et al. and described in Figure 2 by forming its own cluster. Canis simensis geographical data is not found in GBIF however in the literature, it shows an isolated habitat of only Ethiopia (Gottelli et al., 1994). As it forms a separate branch in the phylogeny, can we say it shows allopatric speciation? Due to to wide expanse of canid distribution, it has been difficult to exactly state whether evolution of this genus is based on allopatric or sympatric speciation. The literature, and my own analysis, show a varied evolutionary trend for this genus. However, overall, I believe this genus embodies a sympatric trend of evolution.






#### 9- ACKNOWLEDGMENTS -----

#Nykole Crevits
#Idea for using italics in my histogram plot title.
#Helped find function to change my margins on the dendrogram image.
#Talk through my sequence alignment to find the reason behind misalignment.

#Jacqueline May
#Explained the use of the function apply() to change a character vector into a numeric vector in a matrix.

#Dr. Sarah Adamowicz
#Discussed the reasons behind my misalignment and the different ways I could work through the issue.
#Descriptive lectures and detailed scripts which helped guide my understanding in new R packages and find a pathway through this assignment.

#My study group --> Nykole Crevits, Patricia Balbon and Emily Maier
#Continuous guidance and support while working through this assignment.

#Omar Gonzalez, 2019. R:Set duplicate 'row.names' to a numeric data frame.https://stackoverflow.com/questions/31707654/r-set-duplicate-row-names-to-a-numeric-data-frame

#Rischan, 2014. Italics and normal text in a main plot title. https://stackoverflow.com/questions/23969726/italics-and-normal-text-in-a-main-plot-title

#Tai Galili, 2018. CLuster labels are cut off on horizontal hclust dendrogram. https://stackoverflow.com/questions/51429782/cluster-labels-are-cut-off-on-horizontal-hclust-dendrogram


#### 10- REFERENCES ----

 #1.	Gaubert P, Bloch C, Benyacoub S, Abdelhamid A, Pagani P, Djagoun CAMS, et al. (2012) Reviving the African Wolf Canis lupus lupaster in North and West Africa: A Mitochondrial Lineage Ranging More than 6,000 km Wide. PLoS ONE 7(8): e42740.
# 2.	Gottelli, D., Sillero-Zubiri, C., Applebaum, G. D., Roy, M. S., Girman, D. J., Garcia-Moreno, J., Ostrander, E. A., & Wayne, R. K. (1994). Molecular genetics of the most endangered canid: the Ethiopian wolf Canis simensis. Molecular ecology, 3(4), 301–312.
# 3.	Guo, X., Pei, J., Bao, P., Yan, P., & Lu, D. (2016). The complete mitochondrial genome of Hequ Tibetan Mastiff Canis lupus familiaris (Carnivora: Canidae). Mitochondrial DNA. Part A, DNA mapping, sequencing, and analysis, 27(6), 4659–4660.
# 4.	Kimura M. (1980). A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. Journal of molecular evolution, 16(2), 111–120.
# 5.	 Naohiko Okumura, Naotaka Ishiguro, Masuo Nakano, Akira Matsui, Makoto Sahara (1998). Genetic Variation of the Mitochondrial DNA Cytochrome b Region in Japanese Native Dog Breeds (Canis familiaris). Zoological Science, 15(5), 699-701
# 6.	Padilla, L. R., & Hilton, C. D. (2015). Canidae. Fowler's Zoo and Wild Animal Medicine, Volume 8, 457–467. 
# 7.	Robert K. Wayne, Eli Geffen, Derek J. Girman, Klaus P. Koepfli, Lisa M. Lau, Charles R. Marshall, Molecular Systematics of the Canidae, Systematic Biology, Volume 46, Issue 4, December 1997, Pages 622–653
# 8.	Rumman, U., Sarwar, G., Janjua, S., Khan, F. M., & Nazir, F. (2019). Comparative Mito-Genomic Analysis of Different Species of Genus Canis by Using Different Bioinformatics Tools, Journal of Bioresource Management, 6 (1).
# 9.	Sillero-Zubiri, C., Hoffmann, M. and Macdonald, D.W. (eds). 2004. Canids: Foxes, Wolves, Jackals and Dogs. Status Surveyand Conservation Action Plan. IUCN/SSC Canid Specialist Group. Gland, Switzerland and Cambridge, UK
# 10.	Snyder, S. (n.d.). Wildlife Species: Canis lupus. Canis lupus. Retrieved October 25, 2021, from https://www.fs.fed.us/database/feis/animals/mammal/calu/all.html. 


##R Package Citations

# 1.	Chamberlain S, Barve V, Mcglinn D, Oldoni D, Desmet P, Geffert L, Ram K (2021). _rgbif: Interface to the Global Biodiversity Information Facility API_. R package version 3.6.0
# 2.	Charif D, Lobry J (2007). “SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.” In Bastolla U, Porto M, Roman H, Vendruscolo M (eds.), Structural approaches to sequence evolution: Molecules, networks, populations, series Biological and Medical Physics, Biomedical Engineering, 207-232. Springer Verlag, New York. ISBN : 978-3-540-35305-8.
# 3.	Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res 32:1792-1797.
# 4.	  Enhancements by Thomas P Minka and Alex Deckmyn. (2021). maps: Draw Geographical Maps. R package version 3.4.0
# 5.	Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2.
# 6.	Gagolewski M (2021). _stringi: Fast and portable character string processing in R_. R package version 1.7.5
# 7.	H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2021). Biostrings: Efficient manipulation of biological strings. R package version 2.60.2.
# 8.	Kirill Müller, Hadley Wickham, David A. James and Seth Falcon (2021). RSQLite: 'SQLite' Interface for R. R package version 2.2.8.
# 9.	Original S code by Richard A. Becker and Allan R. Wilks. R version by Ray Brownrigg. (2018). mapdata: Extra Map Databases. R package version 2.3.0.
# 10.	Original S code by Richard A. Becker, Allan R. Wilks. R version by Ray Brownrigg.
# 11.	Paradis E. & Schliep K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35: 526-528.
# 12.	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3 217-223.
# 13.	Schliep, K., Potts, A. J., Morrison, D. A., Grimm, G. W. (2017), Intertwining phylogenetic trees and networks. Methods in Ecology and Evolution, 8: 1212--1220.
# 14.	Simon Garnier, Noam Ross, Robert Rudis, Antônio P. Camargo, Marco Sciaini, and Cédric Scherer (2021). Rvision - Colorblind-Friendly Color Maps for R. R package version 0.6.2.
# 15.	Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686
# 16.	Winter, D. J. (2017) rentrez: an R package for the NCBI eUtils API The R Journal 9(2):520-526
# 17.	Wright ES (2016). “Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R.” _The R Journal_, *8*(1), 352-359.
