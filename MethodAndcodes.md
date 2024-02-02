# Method-zhangkuixing  *2024/02/01*

## Software and Algorithms

| fastp     | Chen et al., 2018                                            | https://github.com/OpenGene/fastp                     |
| --------- | ------------------------------------------------------------ | ----------------------------------------------------- |
| Bowtie2   | [Langmead and Salzberg, 2012](https://www.sciencedirect.com/science/article/pii/S1097276517307578?via%3Dihub#bib23) | http://bowtie-bio.sourceforge.net/bowtie2/index.shtml |
| Samtools  | [Li et al., 2009](https://www.sciencedirect.com/science/article/pii/S1097276517307578?via%3Dihub#bib25) | http://samtools.sourceforge.net/                      |
| Bedtools  | [Quinlan and Hall, 2010](https://www.sciencedirect.com/science/article/pii/S1097276517307578?via%3Dihub#bib35) | http://bedtools.readthedocs.io/en/latest/             |
| sambamba  | Tarasov et al., 2015                                         | https://www.open-bio.org/wiki/Sambamba                |
| picard    | Broad Institute, 2019                                        | https://broadinstitute.github.io/picard/              |
| MACS2     | Zhang et al., 2015                                           | https://github.com/macs3-project/MACS                 |
| MEME      | Timothy L et al., 2015                                       | https://meme-suite.org/meme/tools/meme                |
| R         | N/A                                                          | https://www.r-project.org/                            |
| Deeptools | [Ramírez et al., 2014](https://www.sciencedirect.com/science/article/pii/S1097276517307578?via%3Dihub#bib36) | http://deeptools.readthedocs.io/en/latest/            |

## ChIP-seq and data analysis

ChIP-seq raw FASTQ reads were trimmed by fastp (version 0.12.4) and aligned to hg38 reference genome using bowtie2 aligner (version 2.4.4) with the default setting. BAM files were sorted and filterd by samtools (version 1.9). Enriched peaks were determined using MACS2 (version 2.2.7.1) with -*P-*value < 0.01 cutoff.

## CUT&Tag and data analysis

CUT&Tag raw data was processed consistent with ChIP-seq with modifications: For spike-in normalization, the reads were aligned to reference of hg38 along with the *E. coli* genome by Bowtie2 with the options *--end-to-end --very-sensitive --no-mixed --no-discordant*. And reads of spike-in or not were separated by samtools (version 1.9) option *-b -L*. Results BAM files were filtered using samtools (version 1.9) for a MAPQ score of 20 and deduplicated with picard (version 2.27.1) MarkDuplicates. Blacklisted regions were removed from the BAM file with bedtools (v2.30.0) intersect using ENCODE blacklist bed file (2019, Amemiya, Sci Rep).  Data were normalized for spike-in and library-size with bedtools genomecov (version v2.30.0). Resulted bedgraph file was adapted for peak-calling and converted to a bigwig using ucsc bedGraphToBigWig (version 4) for data visualization and metaplot.

Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z

## Motif analysis for public ChIP-seq peaks

ChIP-seq peaks for 955 proteins were directly downloaded from ENCODE projects. Additional ChIP-seq datasets for 13 proteins were downloaded from GEO, and mapped via bowtie2 (version 2.4.4). PCR duplicates were removed via sambamba (version 0.7.1). Resulted bam files were converted to bigWig files using deeptools (version 3.5.1) with the options *-e 200 --binSize 10 --normalizeUsing CPM*. Reads were extended to the average predicted fragment length, 200bp, for peak calling using MACS2 (version 2.2.7.1) with options *--shift 0 --nolambda --nomodel --keep-dup = all*. Sequence motifs of different lengths (6-8 bp, 6-10 bp and 6-12 bp) were analyzed with MEME (version 5.5.3) for ChIP-seq peaks. The top two motifs for each ChIP-seq dataset were visualized using the motifStack package (version 1.42.0).

## S1-END-seq data analysis

S1-END-seq data in five cell lines (RKO, KM12, SW48, SW620 and SW837) were aligned to reference genome (hg38) using bowtie2 (version 2.4.4) with default parameters. Non-redundant uniquely-mapped reads with kept with samtools (version 1.9). MACS2 (version 2.2.7.1) was used with options *-p 1e-5 --nolambda --nomodel --keep-dup = all* to call peaks. Peaks overlapped with homopyrimidine:homopurine repeats with mirror symmetry were considered as H-DNA peaks. H-DNA peaks in five cell lines were combined as a reference set of H-DNA-forming regions. ChIP-seq peaks within ±500 bp of a H-DNA peak were considered as interacting peaks. If multiple ChIP-seq datasets were available for one protein, the one with highest number of interacting peaks was chosen. Bigwig files were generated in terms of CPM with the bin size = 10 for heatmap and aggregation visualizations with deeptools (version 3.5.1).

## END-seq data analysis

END-seq dataset of doxycycline-induced shWRN and its control (2020, André Nussenzweig, Nature) were downloaded and produced consistently with S1-END-seq analysis we described. Shared peaks of replicates were taken and filtered by a 20-fold enrichment over control. Peaks were merged using bedtools (version 2.30.0) with optin *-d 1000* and set a minimum size of 1 kb. Peaks overlapped with homopyrimidine:homopurine repeats with mirror symmetry were considered as H-DNA peak. 

van Wietmarschen, N., Sridharan, S., Nathan, W.J. *et al.* Repeat expansions confer WRN dependence in microsatellite-unstable cancers. *Nature* **586**, 292–298 (2020). https://doi.org/10.1038/s41586-020-2769-8

## Jaccard similarity coefficients

Homopyrimidine:homopurine regions of any two ChIP-seq peaks were intersected to calculate jaccard similarity coefficient by bedtools *jaccard* (version 2.30.0) and visualized by package pheatmap (version 4.3.2).



# code

## Figure 3

Analysis procedure

We downloaded S1-END-Seq raw data with two replicates of five cell lines (SRR19364482 to SRR19364491), filtered with low quality bases and adapters, and aligned to hg38 reference genome. We keep uniquely mapped reads with high mapping quality and remove PCR duplicates. Bam files were then converted into bigwig files with normalization of CPM, binszie of 10. And MACS2 was used with -nomodel parameter to call peaks, followed with Gabriel et al. For narrow peaks, replicates of samples were intersected and merged with other samples. For merged five cell line peaks (14674), we predicted HDNA mirror repeats loci on these peaks, with total number of 10696. And we seperated the mirror repeats into homopyrimidine and homopurine strand (5461 and 5237, respectively). We then plot signal heatmap and profile of selected mass spectrometry (MS) proteins with their ChIP-seq data, based on the homopyrimidine strand direction, the representive cell line KM12 signals showed the pattern of homopyrimidine enrichment flanking the repeat downstream on the plus strand. And found three patterns that HNRNPK, PCBP1 and PCBP2 showed enriched signals among the upstream of repeat center, DDX5, SMARCA5 and PSIP1 showed features that are consistent with the KM12 HDNA signal, DNMT1, TOP2A and MATR3 showed enriched signal among the downstream of KM12 HDNA signal, with a break point region.
We also compared MS proteins, known triplex proteins and all transcription factor (TF) ChIP-Seq peaks overlapped with S1-END-Seq HDNA mirror peaks, summarized and sorted the number of overlapping. Briefly, we downloaded public TF ChIP-seq data and analyzed MS and known triplex proteins raw data as our method described, and promoter region was removed to avoid transcription hot spot infulence. Then  intersected  with S1-END-seq HDNA peaks. We extended the S1-END-seq HDNA peaks to 500 bp from 5ʹ to 3ʹ, as consideration of some patterns that bind far from mirror center. The results were shown as bar and CDF plot, and heatmap signal profiles were also presented in supplementary materials. What's more , 20 candidates protein peaks were loaded into MEME to search for the motif feature.

```shell
--------------------  This is a shell script  -------------------- 
## Fastqc and align
#Download the raw data of S1-END-Seq for five cell lines
#Download sra data (SRR_Acc_List_S1_END-Seq_FiveCellines.txt)
mkdir SRA_FiveCelline && cd SRA_FiveCelline
mkdir {raw, clean, align, rmdup, macs2}
cat SRR_Acc_List_S1_END-Seq_FiveCellines.txt | while read id; do prefetch ${id}; done
#Convert sra to fastq
cat SRR_Acc_List_S1_END-Seq_FiveCellines.txt | while read id; \
do \
  fastq-dump ./${id}/${id}.sra --split-files --gzip --defline-qual '+' -A ${id} \
  -O ./raw; \
done
#Use fastp to perform quality control and remove low quality bases and adapters
cat SRR_Acc_List_S1_END-Seq_FiveCellines.txt| while read id; \
do \
  fastp -i ./raw/${id}/${id}.fastq.gz \
  -o ./clean/${id}.fq.gz \
  -h ./clean/${id}.fastp.html; \
done

## Align to hg38
#Index
mkdir /bowtie2_index/hg38 && cd /bowtie2_index/hg38
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bowtie2-build hg38.fa hg38 
cd ../
#Align
bowtie2_index = "/bowtie2_index/hg38/hg38"
cat SRR_Acc_List_S1_END-seq_FiveCellines.txt | while read id; \ 
do \
  bowtie2 -p 8 -x ${bowtie2_index} -U ./clean/${id}.fq.gz | \
  samtools view -bS - | \
  samtools view -bh -q 3 - > ./align/${id}.flt.bam \
done
#output minus strand alignments
cd ./align
ls *.flt.bam | while read id; \
do \
  samtools view -f 16 -h ${id} >${id%.*}.reverse.bam;
done
#output plus strand alignments
ls *.flt.bam | while read id; \
do \
  samtools view -F 16 -h ${id} >${id%.*}.forward.bam;
done
#sort
ls *.bam | while read id; \
do \
  samtools sort -O bam -@ 8 ${id} > ./align/${id%.*}.sorted.bam; \
done
#index
ls *.sorted.bam | while read id; do samtools index ${id}; done
#convet to bigwig
ls *.sorted.bam | while read id; \
do \
  bamCoverage --normalizeUsing CPM \
  -b ${id} -o ${id%.*}.bw --binSize 10 -p 6; \
done

## Rmdup
#samtools markdup
ls *.sorted.bam | while read id; \
do \
  samtools markdup -r ${id} ../rmdup/${id%.*}.rmdup.bam; \
done

## peak calling
#sample rename
cd ../
mkdir rename && cd rename
cat rename.txt 
#
SRR19364482	KM12_replicate2
SRR19364483	SW837_replicate2
SRR19364484	RKO_replicate2
SRR19364485	SW48_replicate2
SRR19364486	SW837_replicate1
SRR19364487	SW620_replicate2
SRR19364488	SW620_replicate1
SRR19364489	RKO_replicate1
SRR19364490	SW48_replicate1
SRR19364491	KM12_replicate1
#
for i in {1..11} \
do \
  id=$(cat rename.txt | awk '{print $1}'| head -n $i | tail -1) \
  filename=$(cat rename.txt | awk '{print $2}'| head -n $i | tail -1) \
  mv ../rmdup/${id}.sorted.rmdup.bam ./${filename}.sorted.rmdup.bam | \
done

ls *.sorted.rmdup.bam | while read id; do samtools index ${id}; done

#call peak
ls *.sorted.rmdup.bam | while read id; \
do \
  macs2 callpeak -t ${id} -f BAM -g hs -n ${id%%.*} \
  -p 1e-5 --nolambda --nomodel -B \
  --keep-dup=all \
  --outdir ../macs2/; \
done

cd ../macs2
## Peak merge
#intersected between replicates and merged between samples
mkdir Intersect_replicates && cp *.narrowPeak Intersect_replicates
cd Intersect_replicates
#sort
ls *.narrowPeak | while read id; do cut -f 1-3 ${id} |\
sort -k1,1 -k2,2n -k3,3n > ${id%.*}.sorted.bed; done
#intersect
ls *.sorted.bed | while read id; \
do \
  bedtools intersect -a ${id%%_*}_replicate1_peaks.sorted.bed \
  -b ${id%%_*}_replicate2_peaks.sorted.bed \
  > ${id%%_*}_intersect_peaks.sorted.bed; \
done
#merge
cat *intersect_peaks.sorted.bed > FiveCelline_merge.bed
cat FiveCelline_merge.bed | sort -k1,1 -k2,2n -k 3,3n | sort -u \
bedtools merge -i - > FiveCelline_merge.sorted.merge.bed

## Define H-DNA region
#predict HNDA mirror repeats
cat FiveCelline_merge.sorted.merge.bed | \
bedtools getfasta -fi GRCh38.primary_assembly.genome.fa \
-bed - -fo FiveCelline_merge.sorted.merge.fa
python3 HDNA_Finder.py -d FiveCelline_merge.sorted.merge.fa FiveCelline_merge.sorted.merge.txt (with error rate of 0)
#extract mirror repeat
cat FiveCelline_merge.sorted.merge.txt | grep "mirror"  | cut -f 1,4,5 | \
awk -F"[:-]" 'BEGIN{ OFS="\t"; }{print $1,$2,$3,$4,$5}' | \
awk 'BEGIN{ OFS="\t"; }{print $1,$2+$4-1,$2+$5}' | sort -k1,1 -k2,2n -k3,3n | \
bedtools merge -i - | sort -u > FiveCelline_merge.sorted.merge.mirror.bed

#seperate homopyrimidine and homopurine strand
cat FiveCelline_merge.sorted.merge.bed | \
bedtools getfasta -fi GRCh38.primary_assembly.genome.fa \
-bed - -fo FiveCelline_merge.sorted.merge.mirror.fa

seqkit fx2tab -H -n -i -B a -B t -B c -B g -B ct -B ag \
FiveCelline_merge.sorted.merge.mirror.fa | awk '$6>50' | cut -f 1 | sed '1d' | \
awk -F"[:-]" 'BEGIN{ OFS="\t"; }{print $1,$2,$3}'  > FiveCelline_merge.sorted.HDNA_mirror_CT_rich.bed  #AG: '$7>50'
#turn to bed6 file to prepare for plot by deeptools
cat FiveCelline_merge.sorted.HDNA_mirror_CT_rich.bed | awk 'BEGIN {OFS="\t";} {print $0,"1","2","+"}' > FiveCelline_merge.sorted.HDNA_mirror_CT_rich.strand.bed
cat FiveCelline_merge.sorted.HDNA_mirror_AG_rich.bed | awk 'BEGIN {OFS="\t";} {print $0,"1","2","-"}' > FiveCelline_merge.sorted.HDNA_mirror_AG_rich.strand.bed
cat FiveCelline_merge.sorted.HDNA_mirror_CT_rich.strand.bed FiveCelline_merge.sorted.HDNA_mirror_AG_rich.strand.bed | sort -k1,1 -k2,2n -k3,3n \
> FiveCelline_merge.sorted.HDNA_mirror_CTandAG_rich.strand.sorted.bed


# Signal profile based on homopyrimidine direction repeat
#KM12 plus and minus strand on S1-END-Seq HDNA
computeMatrix reference-point --referencePoint center -p 5 --missingDataAsZero \
-S KM12_replicate1.forward.sorted.bw KM12_replicate1.reverse.sorted.bw \
-R FiveCelline_merge.sorted.HDNA_mirror_CT_rich.bed \ 
-a 1000 -b 1000 -o ./KM12_FiveCellineHDNAmirror_CT_center.gz
#AG strand is the same method

#MS proteins on S1-END-Seq HDNA
#input bed was bed6 with strand information in 6th coloumn, therefore all signals could be plot with one direction by deeptools.
#For KM12 stranded signal, the R script helped to plot with one direction.
computeMatrix reference-point --referencePoint center -p 5 --missingDataAsZero \
-S MS.bw -R FiveCelline_merge.sorted.HDNA_mirror_CTandAG_rich.strand.sorted.bed \
-a 2000 -b 2000 -o MS_S1ENDSeq_CTAG_center_2K.gz

plotHeatmap -m MS_S1ENDSeq_CTAG_center_2K.gz  \
-out MS_S1ENDSeq_CTAG_center_2K_Heatmap.pdf --plotFileFormat pdf \
--dpi 720 --heatmapHeight 8 --heatmapWidth 2 \
--zMin 0 --zMax 0.2 #the color range varies from different MS proteins

plotProfile -m MS_S1ENDSeq_CTAG_center_2K.gz \
-out MS_S1ENDSeq_CTAG_center_2K_Profile.pdf \
--plotFileFormat pdf --perGroup --dpi 720

# Detailed codes of profile and heatmap part in Figure 3
###Heatmap of selected MS proteins###
ls *.bw | while read id; \
do \
  computeMatrix reference-point --referencePoint center \
  -p 8 --missingDataAsZero -S ${id} \
  -R FiveCelline_merge.sorted.HDNA_mirror_CTandAG_rich.strand.sorted.bed \
  -a 2000 -b 2000 -o ${id%.*}_S1ENDSeq_CTAG_center_2K.gz; \
done
plotHeatmap -m DDX5_S1ENDSeq_CTAG_center_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o DDX5_S1ENDSeq_CTAG_center_2K_04color.pdf --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 0.4 --heatmapHeight 8 --heatmapWidth 2
plotHeatmap -m TOP2A_S1ENDSeq_CTAG_center_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o TOP2A_S1ENDSeq_CTAG_center_2K_02color.pdf --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 0.2 --heatmapHeight 8 --heatmapWidth 2
plotHeatmap -m TIF1B_S1ENDSeq_CTAG_center_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o TIF1B_S1ENDSeq_CTAG_center_2K_01color.pdf --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 0.1 --heatmapHeight 8 --heatmapWidth 2
plotHeatmap -m SMARCA5_S1ENDSeq_CTAG_center_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o SMARCA5_S1ENDSeq_CTAG_center_2K_02color.pdf --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 0.2 --heatmapHeight 8 --heatmapWidth 2
plotHeatmap -m PSIP1_S1ENDSeq_CTAG_center_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o PSIP1_S1ENDSeq_CTAG_center_2K_02color.pdf --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 0.2 --heatmapHeight 8 --heatmapWidth 2
plotHeatmap -m PHB2_S1ENDSeq_CTAG_center_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o PHB2_S1ENDSeq_CTAG_center_2K_02color.pdf --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 0.2 --heatmapHeight 8 --heatmapWidth 2
plotHeatmap -m PCBP1_S1ENDSeq_CTAG_center_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o PCBP1_S1ENDSeq_CTAG_center_2K_03color.pdf --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 0.3 --heatmapHeight 8 --heatmapWidth 2
plotHeatmap -m PCBP2_S1ENDSeq_CTAG_center_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o PCBP2_S1ENDSeq_CTAG_center_2K_03color.pdf --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 0.3 --heatmapHeight 8 --heatmapWidth 2
plotHeatmap -m HNRNPK_S1ENDSeq_CTAG_center_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o HNRNPK_S1ENDSeq_CTAG_center_2K_03color.pdf --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 0.3 --heatmapHeight 8 --heatmapWidth 2
plotHeatmap -m MATR3_S1ENDSeq_CTAG_center_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o MATR3_S1ENDSeq_CTAG_center_2K_015color.pdf --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 0.15 --heatmapHeight 8 --heatmapWidth 2
plotHeatmap -m DNMT1_S1ENDSeq_CTAG_center_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o DNMT1_S1ENDSeq_CTAG_center_2K_02color.pdf --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 0.2 --heatmapHeight 8 --heatmapWidth 2

### Profile of three patterns ###
##plot MS signal and KM12 HDNA signal seperately
computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S HNRNPK.bw PCBP1.bw PCBP2.bw -R FiveCelline_merge.sorted.HDNA_mirror_CTandAG_rich.strand.sorted.bed -a 2000 -b 2000 -o HNRNPK_PCBP1_PCBP2_S1ENDSeq_CTAG_center_2K.gz
computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S DDX5.bw SMARCA5.bw PSIP1.bw -R FiveCelline_merge.sorted.HDNA_mirror_CTandAG_rich.strand.sorted.bed -a 2000 -b 2000 -o DDX5_SMARCA5_PSIP1_S1ENDSeq_CTAG_center_2K.gz
computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S DNMT1.bw TOP2A.bw MATR3.bw -R FiveCelline_merge.sorted.HDNA_mirror_CTandAG_rich.strand.sorted.bed -a 2000 -b 2000 -o DNMT1_TOP2A_MATR3_S1ENDSeq_CTAG_center_2K.gz
computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S KM12_forward.bw KM12_reverse.bw -R FiveCelline_merge.sorted.HDNA_mirror_CT_rich.strand.bed -a 2000 -b 2000 -o KM12_S1ENDSeq_CT_center_2K.gz
computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S KM12_forward.bw KM12_reverse.bw -R FiveCelline_merge.sorted.HDNA_mirror_AG_rich.strand.bed -a 2000 -b 2000 -o KM12_S1ENDSeq_AG_center_2K.gz
gunzip *.gz  # Then turn to R script

### TF ChIP-seq overlapped with S1-END-seq HDNA ###
# Download public ChIP-seq data
mkdir TF_ChIP_K562
mkdir TF_ChIP_HepG2
#Download all TF ChIP-seq data for K562 and HepG cell lines from ENCODE, seperately in directory
cat DownloadFile.txt | while read id; do wget ${id}; done
ls *.bed.gz | while read id; do gunzip ${id}; done
cat metadata.tsv | sed '1d' | cut -f 1,23 > ChIP-seq_ID_name.txt
mkdir rename
#sample rename
#!bin/bash
i=0
declare -A dict

cat ChIP-seq_ID_name.txt | while read f
do
  if [ $i -ne 0 ]; then
        full=`echo $f | cut -d " " -f 2`
        ids=`echo $f | cut -d " " -f 1`;
        if [ -e $full ]; then
             mv ../$ids.bed $full.bed
             dict[$full]=$[${dict[$full]}+1]
        else
             dict[$full]=1
             mv ../$ids.bed $full.bed
        fi
        let i++
   else
        let i++
   fi
done

mkdir WithoutPromoter
#remove promoter region
#obtain gencode.v42.primary_assembly.annotation.gtf promoter region
sed 's/"/\t/g' gencode.v42.primary_assembly.annotation.gtf | \
awk 'BEGIN{OFS=FS="\t"}{if($3=="gene") \
  {ensn=$10; symbol=$16; if($7=="+") {start=$4-1; up=start-2000; \
  if(up<0) up=0; dw=start+2000; print $1,up, dw, ensn, symbol, $7;} \
  else if($7=="-") {start=$5-1; up=start+2000; dw=start-2000; if(dw<0) \
  dw=0; print $1,dw,up,ensn,symbol,$7}}} ' | \
  sort -k1,1 -k2,2n > gencode.v42.primary_assembly.promoter.bed
#Take the complement of TF to remove promoter regions
ls *.bed | while read id; \
do \
  bedtools subtract -a ${id} \
  -b gencode.v42.primary_assembly.promoter.bed \
  > ./WithoutPromoter/${id%.*}.withoutPromoter.bed; \
done
cd ./WithoutPromoter

# Extend HDNA by 500bp left and right
bedtools slop -i S1ENDSeq_FiveCelline_merge.sorted.merge.mirror.bed \
-g hg38.chrom.sizes -b 500 \
> FiveCelline_merge.sorted.merge.mirror_500bp.bed

#intesect with S1-END-seq HDNA
ls *.withoutPromoter.bed | while read id; \
do \
  bedtools intersect -a ${id} \
  -b FiveCelline_merge.sorted.merge.mirror_500bp.bed -wa | \
  sort -u | wc -l > ./${id%.*}_Intersect_S1ENDSeqHDNA.txt; \
done
cat *_Intersect_S1ENDSeqHDNA.txt > ChIP-seq_Intersect_S1ENDSeqHDNA.txt  #write intersected number
ls | grep ".withoutPromoter.bed$" > ID.txt #write ChIP-seq TF name
cat ID.txt | cut -d "." -f 1 > ID_cut.txt
paste ChIP-seq_Intersect_S1ENDSeqHDNA.txt ID_cut.txt > ID_peaksIntersectNumber.txt  #merge ID and number
cat ID_peaksIntersectNumber.txt | sort -rn -k1 | less -S #check rank
# Then turn to R script

#For ChIP-seq total peaks intersected with HDNA peak or S1-END-seq shared peaks, the method is similar

--------------------  This is a R script  -------------------- 

# CDF plot in Fig.3b

# profile in Fig.3c-e 

library(tidyverse)
library(stringr)
library(reshape2)
setwd("~/Desktop/TriplexProteome/Bin/S1_END_seq")

##读取数据
data_path_pattern2 <- "DDX5_SMARCA5_PSIP1_S1ENDSeq_shared_CTAG_center_2K"
data_path_pattern3 <- "DNMT1_TOP2A_MATR3_S1ENDSeq_shared_CTAG_center_2K"
data_path_pattern1 <- "HNRNPK_PCBP1_PCBP2_S1ENDSeq_shared_CTAG_center_2k"
data_path_KM12 <- "KM12_S1ENDSeq_shared_CT_center_2K"

CT_pattern1 <- read_delim(data_path_pattern1,
                        delim = "\t",
                        skip = 1,
                        col_names = F) %>%
  na.omit()
CT_pattern2 <- read_delim(data_path_pattern2,
                        delim = "\t",
                        skip = 1,
                        col_names = F) %>%
  na.omit()
CT_pattern3 <- read_delim(data_path_pattern3,
                          delim = "\t",
                          skip = 1,
                          col_names = F) %>%
  na.omit()

CT_KM12 <- read_delim(data_path_KM12,
                      delim = "\t",
                      skip = 1,
                      col_names = F) %>%
  na.omit()

data_agg_merge_CT_pattern1 <- CT_pattern1[,-c(1:3,5,6)] 
data_agg_merge_CT_pattern1_mean <- colMeans(data_agg_merge_CT_pattern1[,-1]) %>%
  as.data.frame()
data_agg_merge_CT_pattern1_mean$bin_number <- c(1:400)
data_agg_merge_CT_pattern1_mean$type[1:400] <- "HNRNPK"
data_agg_merge_CT_pattern1_mean$type[401:800] <- "PCBP1"
data_agg_merge_CT_pattern1_mean$type[801:1200] <- "PCBP2"
#data_agg_merge_CT_strand_mean$type[1201:1600] <- "KM12_fwd"
#data_agg_merge_CT_strand_mean$type[1601:2000] <- "KM12_rev"

## KM12 bw
data_agg_merge_CT_KM12 <- CT_KM12[,-c(1:3,5,6)] 
data_agg_merge_CT_KM12_mean <- colMeans(data_agg_merge_CT_KM12[,-1]) %>%
  as.data.frame()
data_agg_merge_CT_KM12_mean$bin_number <- c(1:400)
data_agg_merge_CT_KM12_mean$type[1:400] <- "Forward"
data_agg_merge_CT_KM12_mean$type[401:800] <- "Reverse"


colnames(data_agg_merge_CT_pattern1_mean) <- c("signal", "bin_number","type")
colnames(data_agg_merge_CT_KM12_mean) <- c("signal", "bin_number","type")
# plot
data_plot_pattern1 <- data_agg_merge_CT_pattern1_mean
data_plot_pattern1$type <- factor(data_plot_pattern1$type, levels  = unique(data_plot_pattern1$type))

data_plot_KM12 <- data_agg_merge_CT_KM12_mean
data_plot_KM12$type <- factor(data_plot_KM12$type, levels  = unique(data_plot_KM12$type))

y_label <- "signal (CPM)"
title_label <- ""
y_limit <- c(round(min(data_plot_KM12$signal)),
             round(max(data_plot_KM12$signal)) + 1)
y_axis <- seq(round(min(data_plot_KM12$signal)),
              round(max(data_plot_KM12$signal))+1,
              2)
x_label <- "Repeat Center (kb)"
legend_label <- ""
fill_value <- c("#AAB1FF","#004CFF","#29FFCE","#000000","#808080")
size <- 4
y <- data_plot_KM12$signal
x <- data_plot_KM12$bin_number

x_limit <- c(0,max(data_plot_KM12$bin_number))
x_axis <- seq(0,max(data_plot_KM12$bin_number),
              max(data_plot_KM12$bin_number)/4)
x_axis_pattern1 <- seq(0,max(data_plot_pattern1$bin_number),
                       max(data_plot_pattern1$bin_number)/4)

legend_position <- c(0.8,0.7)
legend_direction <- "vertical"

top.mar=0.6
right.mar=0.6
bottom.mar=0.6
left.mar=0.6

mytheme<-theme_classic()+
  theme(text=element_text(family = "sans",colour ="black",size = 20),
        legend.text=element_text(colour ="black",size = 15),
        legend.title=element_text(colour ="black",size = 15),
        legend.key.size=unit(6,units = "mm"),
        legend.position="top",
        axis.line = element_line(size = 0.8,colour = "gray30"),
        axis.ticks = element_line(size = 0.8,colour = "gray30"),
        axis.ticks.length = unit(1.8,units = "mm"),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 


library(ggplot2)
library(lubridate)
p1 <- ggplot(data = data_plot_pattern1, aes(x = data_plot_pattern1$bin_number, colors=type)) +
  geom_line(aes(y = ifelse(type %in% c("HNRNPK", "PCBP1", "PCBP2"), signal, NA))) +
  scale_x_continuous(breaks = x_axis_pattern1,
                     labels = c("-2","-1","Center","1","2")) +
  #scale_color_manual(values= c("#67c2a3" ,"#29abe2", "#e889bd")) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, color = 'black'), 
        axis.text.y = element_text(color = 'blue'), axis.ticks.y = element_line(color = 'blue'), 
        axis.title.y = element_text(color = 'blue')) +
  labs(y = 'MS protein signal (CPM)')
p1
p2 <- ggplot(data = data_plot_KM12, aes(x = x,colors=type)) +
  geom_line(aes(y = ifelse(type %in% c("Forward", "Reverse"), signal, NA))) +
  scale_x_continuous(breaks = x_axis,
                     labels = c("-2","-1","Center","1","2")) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, color = 'black'), 
        axis.text.y = element_text(color = 'blue'), axis.ticks.y = element_line(color = 'blue'), 
        axis.title.y = element_text(color = 'blue')) +
  labs(y = 'HDNA signal (CPM)')
p2
y2_plot(p1, p2)

ggsave("Pattern1_shared_S1ENDSeq.pdf")
ggsave("Pattern2_shared_S1ENDSeq.pdf")
ggsave("Pattern3_shared_S1ENDSeq.pdf")

library(gtable)
library(grid)

# run before plot
# Define a function to combine ggplot2 plot results for constructing a dual-axis plot：
# Ref: https://stackoverflow.com/questions/36754891/ggplot2-adding-secondary-y-axis-on-top-of-a-plot

y2_plot <- function(p1, p2) {
  p1 <- ggplotGrob(p1)
  p2 <- ggplotGrob(p2)
  
  # Get the location of the plot panel in p1.
  # These are used later when transformed elements of p2 are put back into p1
  pp <- c(subset(p1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  p1 <- gtable_add_grob(p1, p2$grobs[[which(p2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from p2
  index <- which(p2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- p2$grobs[[index]]                # Extract that grob
  ylab <- hinvert_title_grob(ylab)         # Swap margins and fix justifications
  
  # Put the transformed label on the right side of p1
  p1 <- gtable_add_cols(p1, p2$widths[p2$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from p2 (axis line, tick marks, and tick mark labels)
  index <- which(p2$layout$name == 'axis-l')  # Which grob
  yaxis <- p2$grobs[[index]]                  # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of p1
  p1 <- gtable_add_cols(p1, p2$widths[p2$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(p1)
}


# Bar plot in Fig.S3d 
library(ggplot2)

# read data
known_withPromoter <- read.table("./known_withPromoter_Intersect_S1ENDSeqHDNA.txt",header = F,sep = "\t")
known_withPromoter_shared <- read.table("./known_withPromoter_Intersect_S1ENDSeqHDNA_shared.txt",header = F,sep = "\t")
known_withoutPromoter <- read.table("./known_withoutPromoter_Intersect_S1ENDSeqHDNA_shared.txt",header = F,sep = "\t")
ms_withPromoter <- read.table("./MS_withPromoter_Intersect_S1ENDSeqHDNA.txt",header = F,sep = "\t")
ms_withPromoter_shared<-read.table("./ms_withPromoter_Intersect_S1ENDSeqHDNA_shared.txt",header = F,sep = "\t")
ms_withoutPromoter<-read.table("./MS_withoutPromoter_Intersect_S1ENDSeqHDNA_shared.txt",header = F,sep = "\t")
k562_withPromoter<-read.table("./TF_K562_withPromoter_Intersect_S1ENDSeq.txt",header = F,sep = "\t")
k562_withPromoter_shared<-read.table("./TF_K562_withPromoter_Intersect_S1ENDSeq_shared.txt",header = F,sep = "\t")
k562_withoutPromoter<-read.table("./TF_K562_withoutPromoter_Intersect_S1ENDSeq_shared.txt",header = F,sep = "\t")
hepg2_withPromoter<-read.table("./TF_HepG2_withPromoter_Intersect_S1ENDSeq.txt",header = F,sep = "\t")
hepg2_withPromoter_shared<-read.table("./TF_HepG2_withPromoter_Intersect_S1ENDSeq_shared.txt",header = F,sep = "\t")
hepg2_withoutPromoter<-read.table("TF_HepG2_withoutPromoter_Intersect_S1ENDSeq_shared.txt",header = F,sep = "\t")


k562_withPromoter$V1 <- gsub("-human", "", k562_withPromoter$V1)
k562_withPromoter_shared$V1 <- gsub("-human", "", k562_withPromoter_shared$V1)
k562_withoutPromoter$V1 <- gsub("-human", "", k562_withoutPromoter$V1)
hepg2_withPromoter$V1 <- gsub("-human", "", hepg2_withPromoter$V1)
hepg2_withPromoter_shared$V1 <- gsub("-human", "", hepg2_withPromoter_shared$V1)
hepg2_withoutPromoter$V1 <- gsub("-human", "", hepg2_withoutPromoter$V1)


tf_withPromoter<-merge(k562_withPromoter, hepg2_withPromoter, by = "V1", all = TRUE)
tf_withPromoter[is.na(tf_withPromoter)] <- 0
tf_withPromoter_shared<-merge(k562_withPromoter_shared, hepg2_withPromoter_shared, by = "V1", all = TRUE)
tf_withPromoter_shared[is.na(tf_withPromoter_shared)] <- 0
tf_withoutPromoter<-merge(k562_withoutPromoter, hepg2_withoutPromoter, by = "V1", all = TRUE)
tf_withoutPromoter[is.na(tf_withoutPromoter)] <- 0

# calculate count and total (need total number)
tf_withPromoter$total <- ifelse(tf_withPromoter$V3.x > tf_withPromoter$V3.y, 
                                tf_withPromoter$V2.x, tf_withPromoter$V2.y) # total
tf_withPromoter_shared$total <- ifelse(tf_withPromoter_shared$V3.x > tf_withPromoter_shared$V3.y, 
                                       tf_withPromoter_shared$V2.x, tf_withPromoter_shared$V2.y) # total
tf_withoutPromoter$total <- ifelse(tf_withoutPromoter$V3.x > tf_withoutPromoter$V3.y, 
                                   tf_withoutPromoter$V2.x, tf_withoutPromoter$V2.y) # total
tf_withPromoter$count <- pmax(tf_withPromoter$V3.x, tf_withPromoter$V3.y, na.rm = TRUE) # max
tf_withPromoter_shared$count <- pmax(tf_withPromoter_shared$V3.x, tf_withPromoter_shared$V3.y, na.rm = TRUE) # max
tf_withoutPromoter$count <- pmax(tf_withoutPromoter$V3.x, tf_withoutPromoter$V3.y, na.rm = TRUE) # max


colnames(tf_withPromoter)[1]<-"protein"
colnames(tf_withPromoter_shared)[1]<-"protein"
colnames(tf_withoutPromoter)[1]<-"protein"

TF_withPromoter<-tf_withPromoter[,c(1,6,7)]
TF_withPromoter_shared<-tf_withPromoter_shared[,c(1,6,7)]
TF_withoutPromoter<-tf_withoutPromoter[,c(1,6,7)]

protein_filter<-c("CHD4","DDX17","DDX21","DDX5","DNMT1","HNRNPC","HNRNPK","MATR3","NUP93"
                  ,"PCBP1","PCBP2","PDS5A","PHB2","PSIP1", "SFPQ","SMARCA5","TIF1B","TOP2A","TOP2B","XRCC5",
                  "DDX11","BRIP1","WRN","BLM","DHX9","p53","RPA1","HNRNPL","HNRNPA2B1","PTBP1","ORC4",
                  "NONO","U2AF2","Vim", "GFAP","DES")

TF_withPromoter<-TF_withPromoter[!(TF_withPromoter$protein %in% protein_filter),]
TF_withPromoter_shared<-TF_withPromoter_shared[!(TF_withPromoter_shared$protein %in% protein_filter),]
TF_withoutPromoter<-TF_withoutPromoter[!(TF_withoutPromoter$protein %in% protein_filter),]

TF_withPromoter$type<-"TF_withPromoter"
TF_withPromoter_shared$type<-"TF_withPromoter_shared"
TF_withoutPromoter$type<-"TF_withoutPromoter"

TF_withPromoter$ratio<-TF_withPromoter$count/TF_withPromoter$total
TF_withPromoter_shared$ratio<-TF_withPromoter_shared$count/TF_withPromoter_shared$total
TF_withoutPromoter$ratio<-TF_withoutPromoter$count/TF_withoutPromoter$total

ms_withPromoter$type<-"ms_withPromoter"
ms_withPromoter_shared$type<-"ms_withPromoter_shared"
ms_withoutPromoter$type<-"ms_withoutPromoter"

ms_withPromoter <- ms_withPromoter[order(ms_withPromoter[, 3], decreasing = TRUE), ]
ms_withPromoter_shared <- ms_withPromoter_shared[order(ms_withPromoter_shared[, 3], decreasing = TRUE), ]
ms_withoutPromoter <- ms_withoutPromoter[order(ms_withoutPromoter[, 3], decreasing = TRUE), ]

ms_withPromoter$ratio<-ms_withPromoter$V3/ms_withPromoter$V2
ms_withPromoter_shared$ratio<-ms_withPromoter_shared$V3/ms_withPromoter_shared$V2
ms_withoutPromoter$ratio<-ms_withoutPromoter$V3/ms_withoutPromoter$V2

known_withPromoter$type<-"known_withPromoter"
known_withPromoter_shared$type<-"known_withPromoter_shared"
known_withoutPromoter$type<-"known_withoutPromoter"

known_withPromoter$ratio<-known_withPromoter$V3/known_withPromoter$V2
known_withPromoter_shared$ratio<-known_withPromoter_shared$V3/known_withPromoter_shared$V2
known_withoutPromoter$ratio<-known_withoutPromoter$V3/known_withoutPromoter$V2

colnames(ms_withPromoter)<-c("protein","total","count","type")
colnames(ms_withPromoter_shared)<-c("protein","total","count","type")
colnames(ms_withoutPromoter)<-c("protein","total","count","type")
colnames(known_withPromoter)<-c("protein","total","count","type")
colnames(known_withPromoter_shared)<-c("protein","total","count","type")
colnames(known_withoutPromoter)<-c("protein","total","count","type")
colnames(TF_withPromoter)<-c("protein","total","count","type")
colnames(TF_withPromoter_shared)<-c("protein","total","count","type")
colnames(TF_withoutPromoter)<-c("protein","total","count","type")


MS_withPromoter<-ms_withPromoter[order(-ms_withPromoter$count),]
MS_withPromoter_shared<-ms_withPromoter_shared[order(-ms_withPromoter_shared$count),]
MS_withoutPromoter<-ms_withoutPromoter[order(-ms_withoutPromoter$count),]

df1_withPromoter<-rbind(ms_withPromoter,known_withPromoter)
df1_withPromoter_shared<-rbind(ms_withPromoter_shared,known_withPromoter_shared)
df1_withoutPromoter<-rbind(ms_withoutPromoter,known_withoutPromoter)

df2_withPromoter<-rbind(df1_withPromoter,TF_withPromoter)
df2_withPromoter_shared<-rbind(df1_withPromoter_shared,TF_withPromoter_shared)
df2_withoutPromoter<-rbind(df1_withoutPromoter,TF_withoutPromoter)

#CDF
p1_withPromoter_shared<-ggplot(df2_withPromoter_shared,aes(x=count,color=type))+
  stat_ecdf()+
  labs(x="Number of ChIP-seq peaks \n intersected with shared HDNA",
       y="cumulative distribution funcition (CDF)")+
  scale_color_manual(values =c("#6A3D9A","#1F78B4" ,"#33A02C"))+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.title.y = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.text = element_text(color="black",size=11),
        panel.border = element_rect(linetype = "solid",colour = "black", fill = "NA",size=1),
        panel.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.spacing.x=unit(0.1,'cm'),
        legend.background=element_blank(),
        legend.position = c(1,0), legend.justification = c(1,0))

p1_withPromoter_shared + coord_cartesian(xlim = c(0, 200))

p1_withoutPromoter<-ggplot(df2_withoutPromoter,aes(x=count,color=type))+
  stat_ecdf()+
  labs(x="Number of ChIP-seq peaks \n intersected with shared HDNA",
       y="cumulative distribution funcition (CDF)")+
  scale_color_manual(values =c("#6A3D9A","#1F78B4" ,"#33A02C"))+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.title.y = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.text = element_text(color="black",size=11),
        panel.border = element_rect(linetype = "solid",colour = "black", fill = "NA",size=1),
        panel.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.spacing.x=unit(0.1,'cm'),
        legend.background=element_blank(),
        legend.position = c(1,0), legend.justification = c(1,0))

p1_withoutPromoter + coord_cartesian(xlim = c(0, 200))

# ks.test

ks.test(ms_withoutPromoter$count,tf_withPromoter$count) #p-value = 3e-6
ks.test(known_withPromoter$count,tf_withPromoter$count) #p-value = 4e-3
ks.test(ms_withoutPromoter$count,known_withPromoter$count) #p-value = p = 0.7

ks.test(ms_withoutPromote_shared$count,tf_withoutPromoter_shared$count) #p-value = 5.685e-06
ks.test(known_withoutPromoter_shared$count,tf_withoutPromoter_shared$count) #p-value = 0.007002
ks.test(ms_withoutPromoter_shared$count,known_withoutPromoter_shared$count) #p-value = 0.5259

#bar
colnames(ms_withoutPromoter)<-c("protein","total","count","type","ratio")
colnames(ms_withoutPromoter_shared)<-c("protein","total","count","type","ratio")


p2_withPromoter_shared<-ggplot(ms_withPromoter_shared,aes(x=protein,y=count,fill=type))+
  geom_bar(stat = 'identity',position = "stack")+scale_x_discrete(limits = ms_withPromoter_shared$protein)+
  labs(x=NULL)+
  scale_fill_manual(values = c("#6A3D9A"))+
  theme(panel.background = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.title.y = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.text = element_text(color="black",size=11),axis.text.x = element_text(angle = 90,
                                                                                   hjust = 1),
        panel.border = element_rect(linetype = "solid",colour = "black", fill = "NA",size=1),
        panel.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.spacing.x=unit(0.1,'cm'),
        legend.background=element_blank(),
        legend.position = c(1,1), legend.justification = c(1,1))+
  geom_text(aes(x= protein,y=count,label=count),vjust=-0.5,size=3.5,fontface='bold')+
  #geom_line(aes(x= protein,y=ratio*100*60,group=1),linetype=3,cex=1)+
  #geom_point(aes(x= protein,y=ratio*100*60),color='#33A02C',size=3.5)+
  scale_y_continuous(expand = c(0,0),limits = c(0,65),sec.axis = sec_axis(~ . / 60,
                                                                            name = "Ratio"),name = "Count")
p2_withPromoter_shared

p2_withoutPromoter<-ggplot(ms_withoutPromoter,aes(x=protein,y=count,fill=type))+
  geom_bar(stat = 'identity',position = "stack")+scale_x_discrete(limits = ms_withoutPromoter$protein)+
  labs(x=NULL)+
  scale_fill_manual(values = c("#6A3D9A"))+
  theme(panel.background = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.title.y = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.text = element_text(color="black",size=11),axis.text.x = element_text(angle = 90,
                                                                                   hjust = 1),
        panel.border = element_rect(linetype = "solid",colour = "black", fill = "NA",size=1),
        panel.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.spacing.x=unit(0.1,'cm'),
        legend.background=element_blank(),
        legend.position = c(1,1), legend.justification = c(1,1))+
  geom_text(aes(x= protein,y=count,label=count),vjust=-0.5,size=3.5,fontface='bold')+
  #geom_line(aes(x= protein,y=ratio*100*60,group=1),linetype=3,cex=1)+
  #geom_point(aes(x= protein,y=ratio*100*60),color='#33A02C',size=3.5)+
  scale_y_continuous(expand = c(0,0),limits = c(0,65),sec.axis = sec_axis(~ . / 60,
                                                                            name = "Ratio"),name = "Count")
p2_withoutPromoter
```





## Fig. Rx.

### END-Seq analysis

```shell
--------------------  This is a shell script  -------------------- 

# The mapping method ref to S1-END-seq
cat ENDseq_KM12_DOX_rep1_peaks.narrowPeak | awk '$7>20' > ENDseq_KM12_DOX_rep1_peaks_FC20.bed
cat ENDseq_KM12_DOX_rep2_peaks.narrowPeak | awk '$7>20' > ENDseq_KM12_DOX_rep2_peaks_FC20.bed

bedtools intersect -a ENDseq_KM12_DOX_rep1_peaks_FC20.bed -b ENDseq_KM12_DOX_rep2_peaks_FC20.bed | sort -k1,1 -k2,2n -k3,3n | sort -u | bedtools merge -i - -d 1000 > ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K.srt.bed
cat ENDseq_KM12_DOX_rep12_peaks_FC20.srt.bed | wc -l
6907
cat ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K.srt.bed | wc -l
4323

bedtools getfasta -fi ~/Triplex/Ref_genome/hg38/GRCh38.primary_assembly.genome.fa -bed ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K.srt.bed -fo ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K.srt.fa

# HDNA on total peak and intersect with S1-END-seq
python3 HDNA_Finder.py -d ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K.srt.fa ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA.txt (with error rate of 0)

cat ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA.txt | cut -f 1,4,5 | awk -F"[:-]" 'BEGIN{ OFS="\t"; }{print $1,$2,$3,$4,$5}' | awk 'BEGIN{ OFS="\t"; }{print $1,$2+$4-1,$2+$5}' | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - | wc -l
2096 > ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA.merge.srt.bed

bedtools subtract -a ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K.srt.bed -b ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA.merge.srt.bed -A > ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_subtractHDNA.srt.bed

computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S ~/Triplex/DNMT1.bw ~/Triplex/MATR3.bw ~/Triplex/TOP2A.bw -R ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_subtractHDNA.srt.bed -a 4000 -b 4000 -o ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_subtractHDNA_4K.gz && plotHeatmap -m ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_subtractHDNA_4K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_subtractHDNA_4K_Heatmap.pdf --plotFileFormat pdf --dpi 720 --heatmapHeight 8 --heatmapWidth 2 --zMin 0 --zMax 0.2 --yMin 0.03 --yMax 0.075

computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S ~/Triplex/DNMT1.bw ~/Triplex/MATR3.bw ~/Triplex/TOP2A.bw -R ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA_CTAG_rich_strand.bed -a 4000 -b 4000 -o ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA_CTAG_rich_strand_4K.gz && plotHeatmap -m ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA_CTAG_rich_strand_4K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o ./ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA_CTAG_rich_strand_4K_Heatmap.pdf --plotFileFormat pdf --dpi 720 --heatmapHeight 8 --heatmapWidth 2 --zMin 0 --zMax 0.2 --yMin 0.03 --yMax 0.075

computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S ~/Triplex/DNMT1.bw ~/Triplex/MATR3.bw ~/Triplex/TOP2A.bw -R ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_subtractHDNA.srt.bed -a 2000 -b 2000 -o ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_subtractHDNA_2K.gz && plotHeatmap -m ./ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_subtractHDNA_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o ./ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_subtractHDNA_2K_Heatmap.pdf --plotFileFormat pdf --dpi 720 --heatmapHeight 8 --heatmapWidth 2 --zMin 0 --zMax 0.2 --yMin 0.03 --yMax 0.075
computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S ~/Triplex/DNMT1.bw ~/Triplex/MATR3.bw ~/Triplex/TOP2A.bw -R ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA_CTAG_rich_strand.bed -a 2000 -b 2000 -o ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA_CTAG_rich_strand_2K.gz && plotHeatmap -m ./ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA_CTAG_rich_strand_2K.gz --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd -o ./ENDseq_KM12_DOX_rep12_peaks_FC20_merge1K_HDNA_CTAG_rich_strand_2K_Heatmap.pdf --plotFileFormat pdf --dpi 720 --heatmapHeight 8 --heatmapWidth 2 --zMin 0 --zMax 0.2 --yMin 0.03 --yMax 0.075

```





### S1-END-seq shared peak

```shell
# This is a shell script
cd /share/home/zhangkx/Triplex/DNA_IP_chip-seq/Afteranalysis/HNRNPK/S1_END_seq/SRA_FiveCelline/align/total/macs2/Intersect_replicates/

intervene venn -i KM12_intersect_peaks.sorted.bed RKO_intersect_peaks.sorted.bed SW48_intersect_peaks.sorted.bed SW620_intersect_peaks.sorted.bed SW837_intersect_peaks.sorted.bed --output ./ --save-overlaps




# This is a R script


```





### ChIP-seq peaks HDNA region define and ratio plot

```shell
--------------------  This is a shell script  -------------------- 
cd ~/Triplex/20peaks_candidates
ls *.bed | while read id; do bedtools shuffle -i ${id} -excl ${id} -noOverlapping -seed 123 -g ~/Triplex/Ref_genome/hg38/hg38.chrom.sizes > ./${id%.*}_shuffle.bed; done
ls *.bed | while read id; do bedtools getfasta -fi ~/Triplex/Ref_genome/hg38/GRCh38.primary_assembly.genome.fa -bed ${id} -fo ${id%.*}.fa; done
ls *.fa | while read id; do seqkit rmdup -n ${id} -o ${id%.*}.uniq.fa; done
ls *.uniq.fa | while read id; do python3 ~/Triplex/HDNA_Finder.py -d ${id} ./${id%.*}_HDNA.txt; done # with error rate of 20%
ls *_HDNA.txt | while read id; do cat ${id} | cut -f 1 | sed '1d' | awk -F"[:-]" 'BEGIN{ OFS="\t"; }{print $1,$2,$3,$4,$5}' | sort -k1,1 -k2,2n -k3,3n | sort -u > ./${id%.*}.merge.srt.bed; done

--------------------  This is a R script  -------------------- 
setwd("~/Desktop/TriplexProteome/Bin")

library(readxl)
library(tidyverse)
library(ggplot2)
library(paletteer) 
## plot 20 candidates and their negative control HDNA ratio
HDNA_data <- read_excel("./HDNA_differentError_load.xlsx",
                       sheet = 1)
HDNA_data$E20_per <- HDNA_data$E20 / HDNA_data$peak
HDNA_data$E20_random_per <- HDNA_data$E20_random / HDNA_data$peak
HDNA_data %>%
  select(1,9,10) -> HDNA_data_E20_random
E20_random_point <- 
  ggplot(data = HDNA_data_E20_random, 
         mapping = aes(x = E20_per*100, y = E20_random_per*100)) +
  geom_point(aes(colour=E20_per))+
  scale_color_paletteer_c(palette="ggthemes::Red-Gold", direction = 1, name="Error Rate")+
  # scale_colour_distiller(palette = "PuOr")+
  scale_y_continuous(limits = c(20, 80), expand=expansion(add = c(0, 10)))+ 
  scale_x_continuous(limits = c(20, 100), expand=expansion(add = c(0, 10)))+
  geom_abline(slope = 1, intercept = 0, lty="dashed") + 
  labs(x="ChIP-seq peaks (%)",
       y="Control regions (%)")+
  # theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1))+
  theme(panel.grid = element_blank(),
        legend.position = c(0.3,0.7))
E20_random_point
ggsave("E20_random_point_doubleAxis.pdf")
```



### Jaccard index

```R
--------------------  This is a shell script  -------------------- 
#!/bin/bash

files=("CHD4.withoutPromoter.mergedHDNA.srt.bed" "DDX17.withoutPromoter.mergedHDNA.srt.bed" "DDX21.withoutPromoter.mergedHDNA.srt.bed" "DDX5.withoutPromoter.mergedHDNA.srt.bed" "DNMT1.withoutPromoter.mergedHDNA.srt.bed" "HNRNPC.withoutPromoter.mergedHDNA.srt.bed" "HNRNPK.withoutPromoter.mergedHDNA.srt.bed" "MATR3.withoutPromoter.mergedHDNA.srt.bed" "NUP93.withoutPromoter.mergedHDNA.srt.bed" "PCBP1.withoutPromoter.mergedHDNA.srt.bed" "PCBP2.withoutPromoter.mergedHDNA.srt.bed" "PDS5A.withoutPromoter.mergedHDNA.srt.bed" "PHB2.withoutPromoter.mergedHDNA.srt.bed" "PHB2.withoutPromoter.mergedHDNA.srt.bed" "PSIP1.withoutPromoter.mergedHDNA.srt.bed" "SFPQ.withoutPromoter.mergedHDNA.srt.bed" "SMARCA5.withoutPromoter.mergedHDNA.srt.bed" "TIF1B.withoutPromoter.mergedHDNA.srt.bed" "TOP2A.withoutPromoter.mergedHDNA.srt.bed" "TOP2B.withoutPromoter.mergedHDNA.srt.bed" "XRCC5.withoutPromoter.mergedHDNA.srt.bed")  

matrix_file="jaccard_matrix.txt"
truncate -s 0 "$matrix_file"

for ((i = 0; i < ${#files[@]}; i++)); do
    for ((j = i + 1; j < ${#files[@]}; j++)); do
        file1="${files[$i]}"
        file2="${files[$j]}"
        prefix=$(echo "$file1" | cut -d'.' -f1)_$(echo "$file2" | cut -d'.' -f1)
        bedtools jaccard -a "$file1" -b "$file2" > "${prefix}_jaccard.txt"
        jaccard=$(awk 'NR==2 {print $3}' "${file1%%.*}_${file2%%.*}_jaccard.txt")
        echo -n "$jaccard " >> "$matrix_file"
    done 
    echo >> "$matrix_file"
done

--------------------  This is a R script  -------------------- 
setwd("~/Desktop/TriplexProteome/Bin")
library(readxl)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

jaccard_file <- "jaccard_candicates_csv.csv"

mat <- as.matrix(read.csv(jaccard_file, header = T, row.names = 1))

p1 <- pheatmap(mat, cluster_cols=F, cluster_rows=F,
         color = colorRampPalette(brewer.pal(9, "OrRd"))(100))

p2 <- pheatmap(mat[21-p2$tree_col$order,],
               cluster_cols=T, cluster_rows = F,
               treeheight_col = 0, 
               border_color = NA,
               legend_breaks = seq(0,1,0.2), 
               cellwidth = 7.5, cellheight = 7.5,
               color = colorRampPalette(brewer.pal(9, "OrRd"))(100))
```



### ChIP_seq_DDX3X_yH2AX







### BQQ_Jel466_CUTTAG



```bash
https://github.com/Boyle-Lab/Blacklist/tree/master/lists

cat "/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/4.peakCalling/SEACR/DMSO-1_seacr_top0.01.peaks.stringent.bed" | grep -v ^"chrM" | 
  bedtools intersect -v -a - -b "/share/home/zhangkx/Triplex/Ref_genome/hg38-blacklist.v2.bed" > 
```

call peak and sum coverage

```bash
SEACR_1.3.sh $rootpath/3.filter/bedgraph/DMSO_BQQ_bowtie2.merge.fragments.normalized.srt.bedgraph 0.01 non stringent /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/4.peakCalling/SEACR/DMSO_BQQ_merge_normalized_seacr_top0.01

cat /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/4.peakCalling/SEACR/DMSO_BQQ_merge_normalized_seacr_top0.01.stringent.bed | grep -v ^"chrM" | bedtools intersect -v -a - -b "/share/home/zhangkx/Triplex/Ref_genome/hg38-blacklist.v2.bed" > /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/4.peakCalling/SEACR/DMSO_BQQ_merge_normalized_seacr_top0.01_noBlacklist.bed

multiBamSummary BED-file --BED /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/4.peakCalling/SEACR/DMSO_BQQ_merge_normalized_seacr_top0.01_noBlacklist.bed --bamfiles {rootpath}/3.filter/DMSO_BQQ_bowtie2.q10.merge.srt.bam {rootpath}/3.filter/DMSO_bowtie2.q10.merge.srt.bam {rootpath}/3.filter/BQQ4_bowtie2.q10.merge.srt.bam {rootpath}/3.filter/BQQ24_bowtie2.q10.merge.srt.bam --labels DMSO_BQQ DMSO BQQ4h BQQ24h --outRawCounts /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/4.peakCalling/CDF/coverage/NotSample_multiBamSummary_all_seacr001_noBlacklist_readCounts.tab
```



FC on unified peak

```R
setwd("/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/4.peakCalling/CDF/coverage/")

library(ggplot2)
library(reshape2)

γH2AX <- read.table(file = "./NotSample_multiBamSummary_all_seacr001_noBlacklist_readCounts.tab",col.names = c("chr","start","end","yH2AX","NC","si2","si4"))
Jel466 <- read.table(file = "./SampleMergeByBam2_multiBamSummary_Jel466_SEACR001_Counts.tab",col.names = c("chr","start","end","Jel466","NC","si2","si4"))

γH2AX$NC_spike <- γH2AX$NC*.00937943764643647025
γH2AX$si2_spike <- γH2AX$si2*.01157597660726647203
γH2AX$si4_spike <- γH2AX$si4*.01566244302786348614
γH2AX$si2_ratio <- γH2AX$si2_spike/γH2AX$NC_spike
γH2AX$si4_ratio <- γH2AX$si4_spike/γH2AX$NC_spike


γH2AX$si2_ratio <- γH2AX$si2/γH2AX$NC
γH2AX$si4_ratio <- γH2AX$si4/γH2AX$NC
Jel466$si2_ratio <- Jel466$si2/Jel466$NC
Jel466$si4_ratio <- Jel466$si4/Jel466$NC


γH2AX_ratio <- cbind(γH2AX$si2_ratio, γH2AX$si4_ratio)
colnames(γH2AX_ratio) <- c('siRNA2','siRNA4')
γH2AX_ratio_melt <-
  melt(γH2AX_ratio, variable.name = "siRNA",value.name = "ratio")
γH2AX_ratio_melt_data <- as.data.frame(γH2AX_ratio_melt[,2:3])
colnames(γH2AX_ratio_melt_data) <- c('siRNA','ratio')

Jel466_ratio <- cbind(Jel466$si2_ratio, Jel466$si4_ratio)
colnames(Jel466_ratio) <- c('si2_ratio','si4_ratio')
Jel466_ratio_melt <-
  melt(Jel466_ratio, variable.name = "siRNA",value.name = "ratio")
Jel466_ratio_melt_data <- as.data.frame(Jel466_ratio_melt[,2:3])
colnames(Jel466_ratio_melt_data) <- c('siRNA','ratio')


p_γH2AX <- ggplot(γH2AX_ratio_melt_data,aes(x=log2(ratio),color=siRNA,fill=siRNA)) + 
  geom_histogram(position = "identity",bins = 40,aes(y = ..density..),alpha=0.5,linewidth=0.2) +
  geom_line(stat="density",aes(col=siRNA)) +
  #coord_cartesian(xlim=c(6.5,12)) +
  ylab("Density") +
  xlab("log2(Fold Change of Signal)") +
  geom_vline(xintercept = 0) +
  theme(panel.grid.major = element_line(colour="NA"),
        panel.grid.minor = element_line(colour="NA"),
        panel.background = element_rect(fill="NA"),
        panel.border = element_rect(colour="black", fill=NA)) +
  theme(legend.position=c(0,1), legend.justification=c(0,1),
        legend.background = element_rect(fill="transparent")) +
  scale_fill_manual(values = c("#ef8a62","#67a9cf", "#9c7cbc")) +
  scale_color_manual(values = c("#ef8a62","#67a9cf", "#9c7cbc"))
p_γH2AX
ggsave("Density_histogram_BQQ_notSample.pdf",width =3.44, height =3.25)
ggsave("CDF_line_BQQ_notSample.pdf",width =3.44, height =3.25)


p1 <-ggplot(γH2AX_ratio_melt_data,aes(x=log2(ratio),color=siRNA))+
  stat_ecdf()+
  labs(x="log2(Fold Change of Signal)",
       y="cumulative distribution funcition (CDF)")+
  scale_color_manual(values =c("#ef8a62","#67a9cf", "#9c7cbc"))+
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.title.y = element_text(margin = margin(r =5),size = 11,color="black",face="bold"),
        axis.text = element_text(color="black",size=11),
        panel.border = element_rect(linetype = "solid",colour = "black", fill = "NA",linewidth=1),
        panel.spacing.x = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.text = element_text(color="black",size=10),
        legend.spacing.x=unit(0.1,'cm'),
        legend.background=element_blank(),
        legend.position = c(1,0), legend.justification = c(1,0))
p1+ coord_cartesian(xlim = c(5, 15))
p1

```







upset

```shell
intervene upset -i DMSO_merge_normalized_seacr_top0.01_noBlacklist.bed BQQ4_merge_normalized_seacr_top0.01_noBlacklist.bed BQQ24_merge_normalized_seacr_top0.01_noBlacklist.bed DMSO_BQQ4_merge_normalized_seacr_top0.01_noBlacklist.bed DMSO_BQQ4_BQQ24_merge_normalized_seacr_top0.01_noBlacklist.bed --output ./upset/all --save-overlaps

intervene upset -i DMSO_macs2.bed BQQ4_macs2.bed BQQ24_macs2.bed DMSO_BQQ_macs2.bed --output ./all_macs2 --save-overlaps

cat 1001_DMSO_DMSO_BQQ.bed | awk '{OFS=FS="\t"}{print $0,"D"}' | less -S
cat 1001_DMSO_DMSO_BQQ.bed | awk '{OFS=FS="\t"}{print $0,"DMSO"}' > DMSO_BQQ_Group.DMSO.bed
cat 1101_DMSO_BQQ4_DMSO_BQQ.bed | awk '{OFS=FS="\t"}{print $0,"DMSO-BQQ4"}' > DMSO_BQQ_Group.DB4.bed
cat 1011_DMSO_BQQ24_DMSO_BQQ.bed | awk '{OFS=FS="\t"}{print $0,"DMSO-BQQ24"}' > DMSO_BQQ_Group.DB24.bed
cat 1111_DMSO_BQQ4_BQQ24_DMSO_BQQ.bed| awk '{OFS=FS="\t"}{print $0,"DMSO-BQQ4-BQQ24"}' > DMSO_BQQ_Group.DB4B24.bed
cat 0101_BQQ4_DMSO_BQQ.bed | awk '{OFS=FS="\t"}{print $0,"BQQ4"}' > DMSO_BQQ_Group.B4.bed
cat 0011_BQQ24_DMSO_BQQ.bed | awk '{OFS=FS="\t"}{print $0,"BQQ24"}' > DMSO_BQQ_Group.B24.bed
cat 0111_BQQ4_BQQ24_DMSO_BQQ.bed| awk '{OFS=FS="\t"}{print $0,"BQQ4-BQQ24"}' > DMSO_BQQ_Group.B4B24.bed

wc -l DMSO_BQQ_Group*
 10035 DMSO_BQQ_Group.all.srt.bed
  1046 DMSO_BQQ_Group.B24.bed
   182 DMSO_BQQ_Group.B4B24.bed
   367 DMSO_BQQ_Group.B4.bed
   362 DMSO_BQQ_Group.DB24.bed
   330 DMSO_BQQ_Group.DB4B24.bed
   153 DMSO_BQQ_Group.DB4.bed
   731 DMSO_BQQ_Group.DMSO.bed
  6864 DMSO_BQQ_Group.specific.bed


wc -l DMSO_BQQ_Group*
  1046 DMSO_BQQ_Group.B24.bed
   182 DMSO_BQQ_Group.B4B24.bed
   367 DMSO_BQQ_Group.B4.bed
   362 DMSO_BQQ_Group.DB24.bed
   330 DMSO_BQQ_Group.DB4B24.bed
   153 DMSO_BQQ_Group.DB4.bed
   731 DMSO_BQQ_Group.DMSO.bed
  6864 DMSO_BQQ_Group.specific.bed
 10035 total



# define regions with sort
computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S DMSO_bowtie2.hg38.bw BQQ4_bowtie2.hg38.bw BQQ24_bowtie2.hg38.bw -R "/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/MACS2/upset/all/sets/DMSO_BQQ.group1.bed" -a 2000 -b 2000 -o ./profile/group1_common_macs2_2K.gz --sortRegions no && plotHeatmap -m ./profile/group1_common_macs2_2K.gz -out ./profile/group1_common_macs2_2K.heatmap.pdf --whatToShow 'plot, heatmap and colorbar' --colorList '#2b348e,#60beea,#9bc97c,#f0eb56,#e06b37,#8a2724' --plotFileFormat pdf --dpi 720 --zMin 0 --zMax 1 --heatmapHeight 8 --heatmapWidth 2
# header chane row number
cat group1_common_macs2_2K.header group1_common_macs2_2K.matrix group2_down_macs2_2K.matrix group3_up_macs2_2K.matrix | gzip -c > Group_macs2_2k.gz


computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S DMSO_bowtie2.hg38.bw BQQ4_bowtie2.hg38.bw BQQ24_bowtie2.hg38.bw -R /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/MACS2/upset/all/sets/DMSO_BQQ_Group.B4.bed -a 2000 -b 2000 --sortUsingSamples 1 -o ./profile/group2_B4_sortDMSO_macs2_2k.gz

computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S DMSO_bowtie2.hg38.bw BQQ4_bowtie2.hg38.bw BQQ24_bowtie2.hg38.bw -R /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/MACS2/upset/all/sets/DMSO_BQQ_Group.B24.bed -a 2000 -b 2000 --sortUsingSamples 1 -o ./profile/group2_B24_sortDMSO_macs2_2k.gz

computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S DMSO_bowtie2.hg38.bw BQQ4_bowtie2.hg38.bw BQQ24_bowtie2.hg38.bw -R /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/MACS2/upset/all/sets/DMSO_BQQ_Group.B4B24.bed -a 2000 -b 2000 --sortUsingSamples 1 -o ./profile/group2_B4B24_sortDMSO_macs2_2k.gz

computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S DMSO_bowtie2.hg38.bw BQQ4_bowtie2.hg38.bw BQQ24_bowtie2.hg38.bw -R /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/MACS2/upset/all/sets/DMSO_BQQ_Group.DMSO.bed -a 2000 -b 2000 --sortUsingSamples 1 -o ./profile/group1_DMSO_sortDMSO_macs2_2k.gz

plotHeatmap -m "/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group1_common_macs2_2K.gz" -out /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group1_common_macs2_2K_sortDMSO.heatmp.pdf --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd --plotFileFormat pdf --dpi 720 --heatmapHeight 8 --heatmapWidth 2 --sortUsingSamples 1 --outFileNameMatrix ./profile/group1_common_sortDMSO_macs2_2k.gz


computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S DMSO_bowtie2.hg38.bw BQQ4_bowtie2.hg38.bw BQQ24_bowtie2.hg38.bw -R /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/MACS2/upset/all/sets/DMSO_BQQ_Group.B4.bed -a 2000 -b 2000 -o ./profile/group2_B4_macs2_2k.gz && plotHeatmap -m ./profile/group2_B4_macs2_2k.gz -out ./profile/group2_B4_macs2_2k.heatmap.pdf --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd --plotFileFormat pdf --dpi 720 --heatmapHeight 8 --heatmapWidth 2 --sortUsingSamples 1 --outFileNameMatrix ./profile/group2_B4_sortDMSO_macs2_2k.gz

computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S DMSO_bowtie2.hg38.bw BQQ4_bowtie2.hg38.bw BQQ24_bowtie2.hg38.bw -R /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/MACS2/upset/all/sets/DMSO_BQQ_Group.B24.bed -a 2000 -b 2000 -o ./profile/group2_B24_macs2_2k.gz && plotHeatmap -m ./profile/group2_B24_macs2_2k.gz -out ./profile/group2_B24_macs2_2k.heatmap.pdf --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd --plotFileFormat pdf --dpi 720 --heatmapHeight 8 --heatmapWidth 2 --sortUsingSamples 1 --outFileNameMatrix ./profile/group2_B24_sortDMSO_macs2_2k.gz

computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S DMSO_bowtie2.hg38.bw BQQ4_bowtie2.hg38.bw BQQ24_bowtie2.hg38.bw -R /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/MACS2/upset/all/sets/DMSO_BQQ_Group.B4B24.bed -a 2000 -b 2000 -o ./profile/group2_B4B24_macs2_2k.gz && plotHeatmap -m ./profile/group2_B4B24_macs2_2k.gz -out ./profile/group2_B4B24_macs2_2k.heatmap.pdf --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd --plotFileFormat pdf --dpi 720 --heatmapHeight 8 --heatmapWidth 2 --sortUsingSamples 1 --outFileNameMatrix ./profile/group2_B4B24_sortDMSO_macs2_2k.gz

computeMatrix reference-point --referencePoint center -p 8 --missingDataAsZero -S DMSO_bowtie2.hg38.bw BQQ4_bowtie2.hg38.bw BQQ24_bowtie2.hg38.bw -R /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/MACS2/upset/all/sets/DMSO_BQQ_Group.DMSO.bed -a 2000 -b 2000 -o ./profile/group3_DMSO_macs2_2k.gz && plotHeatmap -m ./profile/group3_DMSO_macs2_2k.gz -out ./profile/group3_DMSO_macs2_2k.heatmap.pdf --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd --plotFileFormat pdf --dpi 720 --heatmapHeight 8 --heatmapWidth 2 --sortUsingSamples 1 --outFileNameMatrix ./profile/group3_DMSO_sortDMSO_macs2_2k.gz

zcat "/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group1_DMSO_sortDMSO_macs2_2k.gz" | tail -n 731 > /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group1_DMSO_sortDMSO_macs2_2k.matrix
zcat "/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B4B24_sortDMSO_macs2_2k.gz" | wc -l 183
zcat "/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B4B24_sortDMSO_macs2_2k.gz" | tail -n 182 > /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B4B24_sortDMSO_macs2_2k.matrix
zcat "/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B24_sortDMSO_macs2_2k.gz" | wc -l
1047
zcat "/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B24_sortDMSO_macs2_2k.gz" | tail -n 1046 > /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B24_sortDMSO_macs2_2k.matrix
zcat "/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B4_sortDMSO_macs2_2k.gz" | wc -l
368
zcat "/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B4_sortDMSO_macs2_2k.gz" | tail -n 367 > /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B4_sortDMSO_macs2_2k.matrix


cat /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/Group_sort_macs2_2k.header /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group1_common_sortDMSO_macs2_2k.matrix /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group3_DMSO_sortDMSO_macs2_2k.matrix /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B4_sortDMSO_macs2_2k.matrix /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B4B24_sortDMSO_macs2_2k.matrix /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/group2_B24_sortDMSO_macs2_2k.matrix | gzip -c > "/share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/Group_sortDMSO_3_macs2_2k.gz"

plotHeatmap -m /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/Group_sortDMSO_3_macs2_2k.gz -out /share/home/zhangkx/Triplex/CUTTAG_BQQ_Jel466_Hela_20231229/New_0108/2.map/bowtie2_result/linryMethod/merge/profile/Group_sortDMSO_3_macs2_2k.heatmap.pdf --whatToShow 'plot, heatmap and colorbar' --colorMap OrRd --plotFileFormat pdf --dpi 720 --heatmapHeight 8 --heatmapWidth 2 --sortRegions no

```





### DDX3XSV40_BQQ_CUTTAG



#CUTTAG process

```shell

#BSUB -J BQQ_triplex
#BSUB -n 25
#BSUB -o BQQ_triplex_out.%J.txt
#BSUB -e BQQ_triplex_error.%J.txt
#BSUB -m node01
rootpath=/share/home/zhangkx/Triplex/CUTTAG_BQQ_DDX3XSV40_hela/
cd ${rootpath}
#mkdir -p {qc,clean,qc_clean,align,rmdup,flt,bw,peaks,stat,part}
Datapath=/share/home/zhangkx/Triplex/CUTTAG_BQQ_DDX3XSV40_hela/Raw/CP2023112100040/H101SC23121304/RSSQ00504/X101SC23121304-Z01/X101SC23121304-Z01-J006/01.RawData/

for name in mock DMSO-1 DMSO-2 DMSO-3 BQQ4-4 BQQ4-5 BQQ4-6 BQQ24-7 BQQ24-8 BQQ24-9
do
        fastqc -t 8 -o ${rootpath}/qc/ ${Datapath}/${name}/${name}_1.fq.gz
       fastqc -t 8 -o ${rootpath}/qc/ ${Datapath}/${name}/${name}_2.fq.gz &
done
wait

for name in mock DMSO-1 DMSO-2 DMSO-3 BQQ4-4 BQQ4-5 BQQ4-6 BQQ24-7 BQQ24-8 BQQ24-9
do
       fastp -i ${Datapath}/${name}/${name}_1.fq.gz -o ${rootpath}/qc_clean/${name}_R1.clean.fq.gz \
       -I ${Datapath}/${name}/${name}_2.fq.gz -O ${rootpath}/qc_clean/${name}_R2.clean.fq.gz &
done
wait

cd ${rootpath}/qc_clean
ls *R1.clean.fq.gz | while read id; do bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant -p 10 -x /share/home/zhangkx/Triplex/bowtie2_index/hg38_spike/hg38_spike -1 $id -2 ${id%_R1*}_R2.clean.fq.gz | samtools sort -@ 10 -o - > ${id%_R1*}.sorted.bam; done

mv *.bam ../align

cd ../align

ls *.bam | xargs -i samtools index {}

ls *.bam | while read id; do samtools flagstat $id > $(basename $id ".bam").stat ; done

mv *.stat ../stat

ls *sorted.bam|while read id ; do samtools view -b -L "/share/home/zhangkx/Triplex/Ref_genome/hg38/hg38.bed" $id > ${id%%.*}.hg38.bam; done

mv *.hg38.bam ../part/

ls *sorted.bam|while read id ; do samtools view -b -L "/share/home/zhangkx/Triplex/Ref_genome/hg38/CnT_SpikeIn.bed" $id > ${id%%.*}.spike.bam; done

mv *.spike.bam ../part/

cd ../part

ls *.bam|while read id; do echo $id; /share/apps/anaconda3/envs/python3/bin/picard -Xmx4g  MarkDuplicates REMOVE_DUPLICATES=true Input=$id OUTPUT=${id%.*}.rmdup.bam METRICS_FILE=${id%.*}.metrics; done

mv *.rmdup.bam ../rmdup

cd ../rmdup

ls *.rmdup.bam|while read id ; do samtools view -F 4 -f 3 -bh -q 20 $id > ${id%.rmdup*}.flt.bam ; done


mv *.flt* ../flt

cd ${rootpath}/flt

ls *.bam|while read id
do
       samtools sort -n $id > ${id%.*}.sortn.bam
done

ls *spike.flt.sortn.bam|while read id
do

        num=`samtools view $id | wc -l`
        num=`expr $num / 2`
        a=1000
        sc=`perl -e "print sprintf('%.4f',$a/$num)"`
        echo "${id}-${sc}" >> scalefactor.txt

        samtools view -H ${id%spike.flt.sortn*}hg38.flt.sortn.bam | grep ^@SQ | cut -d : -f 2- | sed 's/LN://g' | sort -k 1,1 -k 2,2n > ref
        bedtools bamtobed -i ${id%spike.flt.sortn*}hg38.flt.sortn.bam -bedpe | \
        perl -ne 'use List::Util qw(min max); @t=split; print join("\t",$t[0],min($t[1],$t[4]),max($t[2],$t[5]));print "\n"' | sort -k 1,1 -k 2,2n -k 3,3n | \
        bedtools genomecov -i - -g ref -bga -split -scale $sc | \
        sort -k 1,1 -k 2,2n > ${id%spike.flt.sortn*}hg38.bg
        ls *.bg | while read x; do LC_COLLATE=C sort -k1,1 -k2,2n $x > ${x}.sorted ; done
        bedGraphToBigWig ${id%spike.flt.sortn*}hg38.bg.sorted ref ${id%spike.flt.sortn*}hg38.bw && rm ref
done

mv *.bw ../bw

# merge replicates and treatment
samtools merge DMSO.hg38.flt.merge.bam DMSO-1.hg38.flt.bam DMSO-2.hg38.flt.bam DMSO-3.hg38.flt.bam
samtools merge DMSO.spike.flt.merge.bam DMSO-1.spike.flt.bam DMSO-2.spike.flt.bam DMSO-3.spike.flt.bam
samtools merge BQQ4.hg38.flt.merge.bam BQQ4-4.hg38.flt.bam BQQ4-5.hg38.flt.bam BQQ4-6.hg38.flt.bam
samtools merge BQQ4.spike.flt.merge.bam BQQ4-4.spike.flt.bam BQQ4-5.spike.flt.bam BQQ4-6.spike.flt.bam
samtools merge BQQ24.hg38.flt.merge.bam BQQ24-7.hg38.flt.bam BQQ24-8.hg38.flt.bam BQQ24-9.hg38.flt.bam
samtools merge BQQ24.spike.flt.merge.bam BQQ24-7.spike.flt.bam BQQ24-8.spike.flt.bam BQQ24-9.spike.flt.bam
samtools merge DMSO_BQQ4_BQQ24.hg38.flt.merge.bam DMSO.hg38.flt.merge.bam BQQ4.hg38.flt.merge.bam BQQ24.hg38.flt.merge.bam
samtools merge DMSO_BQQ4_BQQ24.spike.flt.merge.bam DMSO.spike.flt.merge.bam BQQ4.spike.flt.merge.bam BQQ24.spike.flt.merge.bam

mkdir -p merge
mv *.merge.bam ./merge
cd ./merge
ls *.merge.bam | while read id; do samtools sort -n ${id} > ${id%.*}.sortn.bam; done
ls *spike.flt.merge.sortn.bam|while read id
do

        num=`samtools view $id | wc -l`
        num=`expr $num / 2`
        a=1000
        sc=`perl -e "print sprintf('%.4f',$a/$num)"`

        samtools view -H ${id%spike.flt.merge.sortn*}hg38.flt.merge.sortn.bam | grep ^@SQ | cut -d : -f 2- | sed 's/LN://g' | sort -k 1,1 -k 2,2n > ref
        bedtools bamtobed -i ${id%spike.flt.merge.sortn*}hg38.flt.merge.sortn.bam -bedpe | \
        perl -ne 'use List::Util qw(min max); @t=split; print join("\t",$t[0],min($t[1],$t[4]),max($t[2],$t[5]));print "\n"' | sort -k 1,1 -k 2,2n -k 3,3n | \
        bedtools genomecov -i - -g ref -bga -split -scale $sc | \
        sort -k 1,1 -k 2,2n > ${id%spike.flt.merge.sortn*}hg38.bg
        ls *.bg | while read x; do LC_COLLATE=C sort -k1,1 -k2,2n $x > ${x}.sorted ; done
        bedGraphToBigWig ${id%spike.flt.merge.sortn*}hg38.bg.sorted ref ${id%spike.flt.merge.sortn*}hg38.bw && rm ref
done
```
