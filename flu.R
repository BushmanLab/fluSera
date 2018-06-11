library(dnar)
library(Biostrings)
library(dplyr)
library(tidyr)
devtools::install_github("kevinmcc21/merf")

# Update these variables as appropriate
projectDir <- "/home/kevin/projects/flu/"
R_image <- paste0(projectDir,"Rimages/mutations.RData")


# Functions ----

processVCF <- function(vcfLine, fastaSeq) {
  pos <- as.numeric(vcfLine[3])
  var <- as.character(vcfLine[6])
  ref <- as.character(vcfLine[5])
  fields <- strsplit(vcfLine[11],":")[[1]]
  if(nchar(ref)>1 | nchar(var)>1) {
    return(unlist(unname(c(vcfLine[1:3],vcfLine[5],var,"n/a","n/a","n/a","n/a",
                           fields[2],fields[5],fields[6],fields[7],fields[8]))))
  }
  codonStart <- 3*floor((pos - 1) / 3) + 1
  codonStop <- codonStart + 2
  oldCodon <- substr(fastaSeq,codonStart,codonStop)
  newCodon <- unlist(strsplit(oldCodon,""))
  newCodon[((pos - 1) %% 3) + 1] <- var
  newCodon <- paste0(newCodon, collapse="")
  return(unlist(unname(c(vcfLine[1:3],ref,var,oldCodon,newCodon,
                         as.character(Biostrings::translate(Biostrings::DNAString(oldCodon))),
                         as.character(Biostrings::translate(Biostrings::DNAString(newCodon))),
                         fields[2],fields[5],fields[6],fields[7],fields[8]))))
}

getSampleType <- function(X) { strsplit(X,"[0-9]")[[1]][1] }
getSampleNumber <- function(X) {  tail(strsplit(X,"[A-Z]")[[1]],1) }

# Sample pre processing ----

system(paste0("wget -O ",
              paste0(projectDir, "sampleGroupMetadata.tsv"),
              " \"https://docs.google.com/spreadsheets/d/1tqRTMY2kXH5DZ5GJtq4kdDaAP5sTvmeGITiWd6JZXqM/pub?gid=0&single=true&output=tsv\""))

sampleMetadata <- read.table(paste0(projectDir, "sampleGroupMetadata.tsv",
                             sep="\t",header = TRUE)
sampleMetadata <- sampleMetadata[with(sampleMetadata, order(SampleType, SampleGroup)),]
save.image(file=R_image)

fa <- paste0(projectDir, "refFinal/allSegs.fa")
outDir <- paste0(projectDir, "allSegs/")

sysPreprocess <- paste0("for R1 in /home/kevin/projects/flu/fastq/*/*R1*; do
  R2=${R1/R1/R2};
  /home/kevin/dev/general/varScanPreprocess.sh ", fa, " ${R1} ${R2} ", outDir,";
  done")
system(paste0("bash -c ", shQuote(sysPreprocess)))

fastq.lines <- system("for fastq in /home/kevin/projects/flu/fastq/*/*R1*; do
                      sample=${fastq##*/}; sample=${sample%_*}; 
                      fl=$(zcat $fastq | wc -l); echo $sample $fl; done",
                      intern=TRUE)
fastq.lines <- t(sapply(strsplit(fastq.lines, " "), function(x) c(x[[1]], x[[2]])))

bam.lines <- system("for bam in /home/kevin/projects/flu/allSegs/bam/*.bam; do
                      sample=${bam##*/}; sample=${sample%%.*};
                      bl=$(samtools view $bam | wc -l); echo $sample $bl; done",
                    intern=TRUE)
bam.lines <- t(sapply(strsplit(bam.lines, " "), function(x) c(x[[1]], x[[2]])))
bam.lines.df <- data.frame(library=bam.lines[,1],
                           numMappedReads=as.numeric(bam.lines[,2]))


pp.summary <- data.frame(library=fastq.lines[,1],
                         numReads=as.numeric(fastq.lines[,2])/2,
                         check.names = FALSE) %>%
  separate(library,c("sampleID","replicateNumber"),remove=FALSE,sep="-",fill="right") %>%
  left_join(bam.lines.df) %>%
  rowwise() %>%
  mutate(sampleType=getSampleType(sampleID),
         sampleNumber=getSampleNumber(sampleID)) %>%
  left_join(sampleMetadata %>% select(SampleID,sampleGroup=SampleGroup),by=c("sampleID"="SampleID")) %>%
  arrange(sampleType,sampleGroup,sampleNumber,replicateNumber)

save.image(file=R_image)

# Call mutations ----

ref <- paste0(projectDir, "refFinal/allSegs.fa")
dir <- paste0(projectDir, "allSegs/")

cmd <- paste0("for PILEUP in ",dir,"pileup/*.pileup;",
              "do SAMPLE=${PILEUP##*/}; SAMPLE=${SAMPLE%%_*};",
              "varscan mpileup2cns ${PILEUP} --variants --output-vcf --min-coverage 30 ",
              "--min-var-freq 0.01 > ",dir,"vcf/${SAMPLE}.vcf; done");
system(cmd)

sysCombineVariants <- paste0("for v in ", dir, "/vcf/*.vcf;
                          do sample=${v##*/}; sample=${sample%%.*}; 
                          grep -v ",'"#"'," $v |
                          awk -v sample=",'"',"$sample",'"'," '{print sample ",'"',"\t",'"'," $0;}';
                          done > ", dir, "/variants.vcf")
system(paste0("bash -c ", shQuote(sysCombineVariants)))

save.image(file=R_image)

# Read in variants ----

fastaOther <- read.fa(paste0(projectDir, "refFinal/allSegs.fa"))
vcfOther <- paste0(projectDir, "allSegs/variants.vcf")
allVariants <-
  as.data.frame(t(apply(read.table(vcfOther, stringsAsFactors = FALSE) %>%
                          mutate(V2=ifelse(is.na(V2),"NA",V2)), 
                        1, function(x) {
                          processVCF(x, fastaOther$seq[fastaOther$name==x["V2"]])
                          })),
                stringsAsFactors=FALSE)

# Filter variants based on segment + position----

allVariants <- allVariants %>%
  filter((V2!="NP" | between(as.numeric(V3),5,1548))) %>%
  filter((V2!="NS1" | between(as.numeric(V3),5,874))) %>%
  filter((V2!="M" | between(as.numeric(V3),5,1011))) %>%
  filter((V2!="PB2" | between(as.numeric(V3),5,2325))) %>%
  filter((V2!="PB1" | between(as.numeric(V3),14,2326))) %>%
  filter((V2!="PA" | between(as.numeric(V3),5,2217))) %>%
  filter((V2!="HA" | as.numeric(V3) < 1667))

allVariants$V3 <- as.numeric(as.character(allVariants$V3))
allVariants$V10 <- as.numeric(as.character(allVariants$V10))
allVariants$V11 <- as.numeric(as.character(allVariants$V11))
allVariants$V12 <- as.numeric(as.character(allVariants$V12))
allVariants$V13 <- as.numeric(sub("%","",as.character(allVariants$V13)))
allVariants$V14 <- as.numeric(as.character(allVariants$V14))

colnames(allVariants) <- c("library", "segment", "position", "reference base", "variant base", 
                           "reference codon", "variant codon", "reference AA", "variant AA",
                           "genotyping quality","reference supporting reads","variant supporting reads",
                           "VAF","p-value")

# Apply other filters & reshape data frame

sampleVariants <- allVariants %>%
  mutate(synonymous = as.character(`reference AA`)==as.character(`variant AA`)) %>%
  separate(library,c("sampleID","replicateNumber"),remove=FALSE,sep="-",fill="right") %>%
  left_join(pp.summary, by=c("library","sampleID")) %>%
  # remove variants from samples with low read count
  # also remove sample E14 (sample preparation errors)
  filter(sampleID!="E14" & numMappedReads > 5000) %>%
  mutate(sampleID=factor(sampleID,levels=unique(sampleID))) %>%
  mutate(AAposition=ceiling(position/3)) %>%
  # adjust amino acid position for HA variants beyond 133
  mutate(AAposition=ifelse(segment=="HA" & AAposition>133,AAposition-1,AAposition)) %>%
  group_by(sampleID, segment, position, AAposition, `reference base`, `variant base`,
           `reference codon`, `variant codon`, `reference AA`, `variant AA`,
           synonymous, sampleGroup, sampleType, sampleNumber) %>%
  # remove variants appearing in only one replicate
  filter(n()>1) %>%
  summarize(VAFmean=round(mean(VAF),1),VAFsd=sd(VAF)) %>% ungroup() %>%
  arrange(factor(segment,levels=c("HA","NA","M","NP","NS1","PA","PB1","PB2"))) %>%
  # remove variants with VAF < 2.5%
  filter(VAFmean>2.5)

variants <- sampleVariants[sampleVariants$synonymous==FALSE,
                           c("segment", "sampleID","sampleGroup",
                             "sampleType", "position", "AAposition",
                             "reference base", "variant base", "reference AA",
                             "variant AA", "VAFmean")]
variants$sampleGroup <- 
  factor(variants$sampleGroup,levels=unique(variants$sampleGroup))
rownames(variants) <- NULL

save.image(file=R_image)

# Variant plot data ----

numSamples <- length(unique(sampleVariants$sampleID))

HA.xMin <- 100
HA.xMax <- 250
HA.xRange <- HA.xMax - HA.xMin + 1


plotPoints <- data.frame(AAposition=rep(seq(HA.xMin,HA.xMax),numSamples))

library(dplyr)

HA.plotData <- sampleVariants %>%
  filter(segment=="HA" & between(AAposition,HA.xMin,HA.xMax)) %>%
  left_join(sampleMetadata %>% select(SampleID,SampleName), by=c("sampleID"="SampleID")) %>%
  select(AAposition, VAFmean, VAFsd, sampleID, sampleGroup, sampleType, SampleName) %>%
  group_by(AAposition,sampleID) %>%
  mutate(varNo=as.character(seq(1,n()))) %>% ungroup()

HA.plotData1 <- HA.plotData %>% filter(sampleID %in% c("E2","E4","E6")) %>% droplevels()
HA.plotData2 <- HA.plotData %>% filter(sampleID %in% c("E15","E16","E18","E19","E20")) %>% droplevels()

save.image(file=R_image)




