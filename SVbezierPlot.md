---
title: "Bezier curves-like plot for structural variants visualisation"
author: "Jacek Marzec"
date: "25/05/2018"
output: 
  html_document:
    keep_md: yes
---

Script generating a bezier curves-like plot representing structural variants (SVs) listed in specified BEDPE file. Currently the script plots take into account only the start of each segement ("START_A" and "START_B" columns in the BEDPE file) to plot the bezier curves between the different genomic locations of the event.
There is also a command line version of this script where one can specify the input file (\--bedpe), output plot name (\--output), as well as the genome build version (\--g_build) using the command line arguments.


#### Load libraries


```r
suppressMessages(library(ggplot2))
suppressMessages(library(ggforce))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
```

#### Read BEDPE file


Read the BEDPE file and extract the SVs info including start and end chromosomes, positions, as well as the variant type


```r
##### Read the BEDPE file
sv.data <- read.table(file="example_data/structural/example-manta-pass.bedpe", sep="\t", as.is=TRUE, header=TRUE, comment.char="")

##### Get the SV chromosomes...
chr1 <- paste0("chr", sv.data$X.CHROM_A)
chr2 <- paste0("chr", sv.data$CHROM_B)

##### ...positions
pos1 <- sv.data$START_A
pos2 <- sv.data$START_B

##### ... and the events type
type <- sv.data$ID
type = gsub("Manta","", type)
type = gsub(":\\d+","", type)
```

#### Prepare x-axis coordinates info for ggplot

This part of the script converts the genomic positions from hg38 (other genome build can be selected using the command line script) to coordinates that can be plotted on the ggplot x-axis.

Start with calculating the whole genome length. Here we consider chromosomes 1-22, X and Y


```r
genome.length = sum(seqlengths(Hsapiens)[1:24])
```


Now calculate fake chromosomes' start positions so that they match with the x-axis coordinates in the ggplot


```r
chrs_fake_starts <- vector("list", 24)
chrs_fake_starts  <- setNames(chrs_fake_starts,  names(Hsapiens)[1:24] )
```

Chromosome 1 has coordingate 0


```r
chrs_fake_starts[["chr1"]] <- 0
```

The coordinates for the remaining chromosomes will be calculated by adding the lengths of individual preceding chromosomes


```r
length_sum <- 0

for ( i in 2:length(chrs_fake_starts) ) {

#	cat(paste0("\nCalculations for " , names(chrs_fake_starts)[i], "...", sep=""))
	cat(paste("\nThe fake start position for " , names(chrs_fake_starts)[i], " is ", length_sum + as.numeric(seqlengths(Hsapiens)[[i-1]]), sep=""))
#	cat(paste("\nLength of " , names(chrs_fake_starts)[i-1], " = ", as.numeric(seqlengths(Hsapiens)[[i-1]]), " and the sum of the preceding chromosomes = ", length_sum, ".\n\n", sep=""))

	length_sum <- length_sum + as.numeric(seqlengths(Hsapiens)[[i-1]])
	chrs_fake_starts[[names(Hsapiens)[i]]] <- length_sum
}
```

```
## 
## The fake start position for chr2 is 248956422
## The fake start position for chr3 is 491149951
## The fake start position for chr4 is 689445510
## The fake start position for chr5 is 879660065
## The fake start position for chr6 is 1061198324
## The fake start position for chr7 is 1232004303
## The fake start position for chr8 is 1391350276
## The fake start position for chr9 is 1536488912
## The fake start position for chr10 is 1674883629
## The fake start position for chr11 is 1808681051
## The fake start position for chr12 is 1943767673
## The fake start position for chr13 is 2077042982
## The fake start position for chr14 is 2191407310
## The fake start position for chr15 is 2298451028
## The fake start position for chr16 is 2400442217
## The fake start position for chr17 is 2490780562
## The fake start position for chr18 is 2574038003
## The fake start position for chr19 is 2654411288
## The fake start position for chr20 is 2713028904
## The fake start position for chr21 is 2777473071
## The fake start position for chr22 is 2824183054
## The fake start position for chrX is 2875001522
## The fake start position for chrY is 3031042417
```

Calculate the coordinates for x-axis labels (chr1, chr2...) for ggplot by adding the half-lenght of each chrosomome to its fake start


```r
chrs_fake_label.pos <- vector("list", 24)
chrs_fake_label.pos  <- setNames(chrs_fake_label.pos,  names(Hsapiens)[1:24] )

for ( i in 1:length(chrs_fake_starts) ) {

	chrs_fake_label.pos[[names(Hsapiens)[i]]] <- seqlengths(Hsapiens)[[i]]/2 + chrs_fake_starts[[names(Hsapiens)[i]]]

	cat(paste("\nThe x-axis coordinate for " , names(chrs_fake_starts)[i], " label is ", chrs_fake_label.pos[[names(Hsapiens)[i]]], " = ",  seqlengths(Hsapiens)[[i]]/2, " (half-length) + ", chrs_fake_starts[[names(Hsapiens)[i]]]," (fake start)", sep=""))
#	cat(paste("\nLength of " , names(chrs_fake_starts)[i], " = ", seqlengths(Hsapiens)[[i]], " and the fake start of chromosome ", names(Hsapiens)[i], " = ", chrs_fake_starts[[names(Hsapiens)[i]]], ".\n\n", sep=""))
}
```

```
## 
## The x-axis coordinate for chr1 label is 124478211 = 124478211 (half-length) + 0 (fake start)
## The x-axis coordinate for chr2 label is 370053186.5 = 121096764.5 (half-length) + 248956422 (fake start)
## The x-axis coordinate for chr3 label is 590297730.5 = 99147779.5 (half-length) + 491149951 (fake start)
## The x-axis coordinate for chr4 label is 784552787.5 = 95107277.5 (half-length) + 689445510 (fake start)
## The x-axis coordinate for chr5 label is 970429194.5 = 90769129.5 (half-length) + 879660065 (fake start)
## The x-axis coordinate for chr6 label is 1146601313.5 = 85402989.5 (half-length) + 1061198324 (fake start)
## The x-axis coordinate for chr7 label is 1311677289.5 = 79672986.5 (half-length) + 1232004303 (fake start)
## The x-axis coordinate for chr8 label is 1463919594 = 72569318 (half-length) + 1391350276 (fake start)
## The x-axis coordinate for chr9 label is 1605686270.5 = 69197358.5 (half-length) + 1536488912 (fake start)
## The x-axis coordinate for chr10 label is 1741782340 = 66898711 (half-length) + 1674883629 (fake start)
## The x-axis coordinate for chr11 label is 1876224362 = 67543311 (half-length) + 1808681051 (fake start)
## The x-axis coordinate for chr12 label is 2010405327.5 = 66637654.5 (half-length) + 1943767673 (fake start)
## The x-axis coordinate for chr13 label is 2134225146 = 57182164 (half-length) + 2077042982 (fake start)
## The x-axis coordinate for chr14 label is 2244929169 = 53521859 (half-length) + 2191407310 (fake start)
## The x-axis coordinate for chr15 label is 2349446622.5 = 50995594.5 (half-length) + 2298451028 (fake start)
## The x-axis coordinate for chr16 label is 2445611389.5 = 45169172.5 (half-length) + 2400442217 (fake start)
## The x-axis coordinate for chr17 label is 2532409282.5 = 41628720.5 (half-length) + 2490780562 (fake start)
## The x-axis coordinate for chr18 label is 2614224645.5 = 40186642.5 (half-length) + 2574038003 (fake start)
## The x-axis coordinate for chr19 label is 2683720096 = 29308808 (half-length) + 2654411288 (fake start)
## The x-axis coordinate for chr20 label is 2745250987.5 = 32222083.5 (half-length) + 2713028904 (fake start)
## The x-axis coordinate for chr21 label is 2800828062.5 = 23354991.5 (half-length) + 2777473071 (fake start)
## The x-axis coordinate for chr22 label is 2849592288 = 25409234 (half-length) + 2824183054 (fake start)
## The x-axis coordinate for chrX label is 2953021969.5 = 78020447.5 (half-length) + 2875001522 (fake start)
## The x-axis coordinate for chrY label is 3059656124.5 = 28613707.5 (half-length) + 3031042417 (fake start)
```

#### Calculate ggplot x-axis coordinates for SV

Calculate the coordinates to draw bezier curves by adding the SV position info to the fake start coordinates of corresponding chromosomes


```r
pos1_fake <- vector("list", nrow(sv.data))
pos2_fake <- vector("list", nrow(sv.data))

for ( i in 1:nrow(sv.data) ) {

	cat(paste("\nCalculations for SV: " , paste( chr1[i], pos1[i], sep=" " ), "-",  paste( chr2[i], pos2[i], sep=" " ), sep=""))
	cat(paste("\nThe x-axis coordinate for position 1 is ", chrs_fake_starts[[chr1[i]]] + pos1[i], " = ",  chrs_fake_starts[[chr1[i]]], " (the fake start of ", chr1[i],") + ", pos1[i], " (the real position 1)", sep=""))
	cat(paste("\nThe x-axis coordinate for position 2 is ", chrs_fake_starts[[chr2[i]]] + pos2[i], " = ",  chrs_fake_starts[[chr2[i]]], " (the fake start of ", chr2[i],") + ", pos2[i], " (the real position 2).\n", sep=""))

	pos1_fake[[i]] <- chrs_fake_starts[[chr1[i]]] + pos1[i]
	pos2_fake[[i]] <- chrs_fake_starts[[chr2[i]]] + pos2[i]
}
```

```
## 
## Calculations for SV: chr1 29416966-chr5 157376567
## The x-axis coordinate for position 1 is 29416966 = 0 (the fake start of chr1) + 29416966 (the real position 1)
## The x-axis coordinate for position 2 is 1037036632 = 879660065 (the fake start of chr5) + 157376567 (the real position 2).
## 
## Calculations for SV: chr1 39559701-chr9 13237205
## The x-axis coordinate for position 1 is 39559701 = 0 (the fake start of chr1) + 39559701 (the real position 1)
## The x-axis coordinate for position 2 is 1549726117 = 1536488912 (the fake start of chr9) + 13237205 (the real position 2).
## 
## Calculations for SV: chr1 211768443-chr7 121001796
## The x-axis coordinate for position 1 is 211768443 = 0 (the fake start of chr1) + 211768443 (the real position 1)
## The x-axis coordinate for position 2 is 1353006099 = 1232004303 (the fake start of chr7) + 121001796 (the real position 2).
## 
## Calculations for SV: chr1 226466151-chr5 139934838
## The x-axis coordinate for position 1 is 226466151 = 0 (the fake start of chr1) + 226466151 (the real position 1)
## The x-axis coordinate for position 2 is 1019594903 = 879660065 (the fake start of chr5) + 139934838 (the real position 2).
## 
## Calculations for SV: chr1 233997283-chr3 71583534
## The x-axis coordinate for position 1 is 233997283 = 0 (the fake start of chr1) + 233997283 (the real position 1)
## The x-axis coordinate for position 2 is 562733485 = 491149951 (the fake start of chr3) + 71583534 (the real position 2).
## 
## Calculations for SV: chr2 27815956-chrX 134327508
## The x-axis coordinate for position 1 is 276772378 = 248956422 (the fake start of chr2) + 27815956 (the real position 1)
## The x-axis coordinate for position 2 is 3009329030 = 2875001522 (the fake start of chrX) + 134327508 (the real position 2).
## 
## Calculations for SV: chr2 131910338-chr3 13604077
## The x-axis coordinate for position 1 is 380866760 = 248956422 (the fake start of chr2) + 131910338 (the real position 1)
## The x-axis coordinate for position 2 is 504754028 = 491149951 (the fake start of chr3) + 13604077 (the real position 2).
## 
## Calculations for SV: chr2 154025057-chr4 66465079
## The x-axis coordinate for position 1 is 402981479 = 248956422 (the fake start of chr2) + 154025057 (the real position 1)
## The x-axis coordinate for position 2 is 755910589 = 689445510 (the fake start of chr4) + 66465079 (the real position 2).
## 
## Calculations for SV: chr3 13604073-chr2 131910342
## The x-axis coordinate for position 1 is 504754024 = 491149951 (the fake start of chr3) + 13604073 (the real position 1)
## The x-axis coordinate for position 2 is 380866764 = 248956422 (the fake start of chr2) + 131910342 (the real position 2).
## 
## Calculations for SV: chr3 58018758-chr14 73879849
## The x-axis coordinate for position 1 is 549168709 = 491149951 (the fake start of chr3) + 58018758 (the real position 1)
## The x-axis coordinate for position 2 is 2265287159 = 2191407310 (the fake start of chr14) + 73879849 (the real position 2).
## 
## Calculations for SV: chr3 71583528-chr1 233997289
## The x-axis coordinate for position 1 is 562733479 = 491149951 (the fake start of chr3) + 71583528 (the real position 1)
## The x-axis coordinate for position 2 is 233997289 = 0 (the fake start of chr1) + 233997289 (the real position 2).
## 
## Calculations for SV: chr3 89421180-chr16 27166429
## The x-axis coordinate for position 1 is 580571131 = 491149951 (the fake start of chr3) + 89421180 (the real position 1)
## The x-axis coordinate for position 2 is 2427608646 = 2400442217 (the fake start of chr16) + 27166429 (the real position 2).
## 
## Calculations for SV: chr3 174703742-chr6 3784082
## The x-axis coordinate for position 1 is 665853693 = 491149951 (the fake start of chr3) + 174703742 (the real position 1)
## The x-axis coordinate for position 2 is 1064982406 = 1061198324 (the fake start of chr6) + 3784082 (the real position 2).
## 
## Calculations for SV: chr4 21404394-chr6 29520405
## The x-axis coordinate for position 1 is 710849904 = 689445510 (the fake start of chr4) + 21404394 (the real position 1)
## The x-axis coordinate for position 2 is 1090718729 = 1061198324 (the fake start of chr6) + 29520405 (the real position 2).
## 
## Calculations for SV: chr4 66465068-chr2 154025068
## The x-axis coordinate for position 1 is 755910578 = 689445510 (the fake start of chr4) + 66465068 (the real position 1)
## The x-axis coordinate for position 2 is 402981490 = 248956422 (the fake start of chr2) + 154025068 (the real position 2).
## 
## Calculations for SV: chr4 85444150-chr10 105326115
## The x-axis coordinate for position 1 is 774889660 = 689445510 (the fake start of chr4) + 85444150 (the real position 1)
## The x-axis coordinate for position 2 is 1780209744 = 1674883629 (the fake start of chr10) + 105326115 (the real position 2).
## 
## Calculations for SV: chr5 139934831-chr1 226466158
## The x-axis coordinate for position 1 is 1019594896 = 879660065 (the fake start of chr5) + 139934831 (the real position 1)
## The x-axis coordinate for position 2 is 226466158 = 0 (the fake start of chr1) + 226466158 (the real position 2).
## 
## Calculations for SV: chr5 149927236-chr7 31619905
## The x-axis coordinate for position 1 is 1029587301 = 879660065 (the fake start of chr5) + 149927236 (the real position 1)
## The x-axis coordinate for position 2 is 1263624208 = 1232004303 (the fake start of chr7) + 31619905 (the real position 2).
## 
## Calculations for SV: chr5 157376560-chr1 29416973
## The x-axis coordinate for position 1 is 1037036625 = 879660065 (the fake start of chr5) + 157376560 (the real position 1)
## The x-axis coordinate for position 2 is 29416973 = 0 (the fake start of chr1) + 29416973 (the real position 2).
## 
## Calculations for SV: chr6 3784070-chr3 174703754
## The x-axis coordinate for position 1 is 1064982394 = 1061198324 (the fake start of chr6) + 3784070 (the real position 1)
## The x-axis coordinate for position 2 is 665853705 = 491149951 (the fake start of chr3) + 174703754 (the real position 2).
## 
## Calculations for SV: chr6 29520398-chr4 21404401
## The x-axis coordinate for position 1 is 1090718722 = 1061198324 (the fake start of chr6) + 29520398 (the real position 1)
## The x-axis coordinate for position 2 is 710849911 = 689445510 (the fake start of chr4) + 21404401 (the real position 2).
## 
## Calculations for SV: chr6 107053182-chr14 92995870
## The x-axis coordinate for position 1 is 1168251506 = 1061198324 (the fake start of chr6) + 107053182 (the real position 1)
## The x-axis coordinate for position 2 is 2284403180 = 2191407310 (the fake start of chr14) + 92995870 (the real position 2).
## 
## Calculations for SV: chr7 26951381-chr19 19107309
## The x-axis coordinate for position 1 is 1258955684 = 1232004303 (the fake start of chr7) + 26951381 (the real position 1)
## The x-axis coordinate for position 2 is 2673518597 = 2654411288 (the fake start of chr19) + 19107309 (the real position 2).
## 
## Calculations for SV: chr7 31619899-chr5 149927242
## The x-axis coordinate for position 1 is 1263624202 = 1232004303 (the fake start of chr7) + 31619899 (the real position 1)
## The x-axis coordinate for position 2 is 1029587307 = 879660065 (the fake start of chr5) + 149927242 (the real position 2).
## 
## Calculations for SV: chr7 47593877-chr22 32216386
## The x-axis coordinate for position 1 is 1279598180 = 1232004303 (the fake start of chr7) + 47593877 (the real position 1)
## The x-axis coordinate for position 2 is 2856399440 = 2824183054 (the fake start of chr22) + 32216386 (the real position 2).
## 
## Calculations for SV: chr7 121001792-chr1 211768447
## The x-axis coordinate for position 1 is 1353006095 = 1232004303 (the fake start of chr7) + 121001792 (the real position 1)
## The x-axis coordinate for position 2 is 211768447 = 0 (the fake start of chr1) + 211768447 (the real position 2).
## 
## Calculations for SV: chr8 57311101-chr8 57311403
## The x-axis coordinate for position 1 is 1448661377 = 1391350276 (the fake start of chr8) + 57311101 (the real position 1)
## The x-axis coordinate for position 2 is 1448661679 = 1391350276 (the fake start of chr8) + 57311403 (the real position 2).
## 
## Calculations for SV: chr9 13237174-chr1 39559732
## The x-axis coordinate for position 1 is 1549726086 = 1536488912 (the fake start of chr9) + 13237174 (the real position 1)
## The x-axis coordinate for position 2 is 39559732 = 0 (the fake start of chr1) + 39559732 (the real position 2).
## 
## Calculations for SV: chr9 103093603-chr16 25645303
## The x-axis coordinate for position 1 is 1639582515 = 1536488912 (the fake start of chr9) + 103093603 (the real position 1)
## The x-axis coordinate for position 2 is 2426087520 = 2400442217 (the fake start of chr16) + 25645303 (the real position 2).
## 
## Calculations for SV: chr10 52166755-chr12 57788228
## The x-axis coordinate for position 1 is 1727050384 = 1674883629 (the fake start of chr10) + 52166755 (the real position 1)
## The x-axis coordinate for position 2 is 2001555901 = 1943767673 (the fake start of chr12) + 57788228 (the real position 2).
## 
## Calculations for SV: chr10 105326110-chr4 85444155
## The x-axis coordinate for position 1 is 1780209739 = 1674883629 (the fake start of chr10) + 105326110 (the real position 1)
## The x-axis coordinate for position 2 is 774889665 = 689445510 (the fake start of chr4) + 85444155 (the real position 2).
## 
## Calculations for SV: chr12 11271298-chr19 15581760
## The x-axis coordinate for position 1 is 1955038971 = 1943767673 (the fake start of chr12) + 11271298 (the real position 1)
## The x-axis coordinate for position 2 is 2669993048 = 2654411288 (the fake start of chr19) + 15581760 (the real position 2).
## 
## Calculations for SV: chr12 49559915-chr13 48666269
## The x-axis coordinate for position 1 is 1993327588 = 1943767673 (the fake start of chr12) + 49559915 (the real position 1)
## The x-axis coordinate for position 2 is 2125709251 = 2077042982 (the fake start of chr13) + 48666269 (the real position 2).
## 
## Calculations for SV: chr12 57788225-chr10 52166758
## The x-axis coordinate for position 1 is 2001555898 = 1943767673 (the fake start of chr12) + 57788225 (the real position 1)
## The x-axis coordinate for position 2 is 1727050387 = 1674883629 (the fake start of chr10) + 52166758 (the real position 2).
## 
## Calculations for SV: chr12 120146501-chr16 14583834
## The x-axis coordinate for position 1 is 2063914174 = 1943767673 (the fake start of chr12) + 120146501 (the real position 1)
## The x-axis coordinate for position 2 is 2415026051 = 2400442217 (the fake start of chr16) + 14583834 (the real position 2).
## 
## Calculations for SV: chr13 48666260-chr12 49559924
## The x-axis coordinate for position 1 is 2125709242 = 2077042982 (the fake start of chr13) + 48666260 (the real position 1)
## The x-axis coordinate for position 2 is 1993327597 = 1943767673 (the fake start of chr12) + 49559924 (the real position 2).
## 
## Calculations for SV: chr14 73879831-chr3 58018776
## The x-axis coordinate for position 1 is 2265287141 = 2191407310 (the fake start of chr14) + 73879831 (the real position 1)
## The x-axis coordinate for position 2 is 549168727 = 491149951 (the fake start of chr3) + 58018776 (the real position 2).
## 
## Calculations for SV: chr14 92995864-chr6 107053188
## The x-axis coordinate for position 1 is 2284403174 = 2191407310 (the fake start of chr14) + 92995864 (the real position 1)
## The x-axis coordinate for position 2 is 1168251512 = 1061198324 (the fake start of chr6) + 107053188 (the real position 2).
## 
## Calculations for SV: chr16 14583823-chr12 120146512
## The x-axis coordinate for position 1 is 2415026040 = 2400442217 (the fake start of chr16) + 14583823 (the real position 1)
## The x-axis coordinate for position 2 is 2063914185 = 1943767673 (the fake start of chr12) + 120146512 (the real position 2).
## 
## Calculations for SV: chr16 25645288-chr9 103093618
## The x-axis coordinate for position 1 is 2426087505 = 2400442217 (the fake start of chr16) + 25645288 (the real position 1)
## The x-axis coordinate for position 2 is 1639582530 = 1536488912 (the fake start of chr9) + 103093618 (the real position 2).
## 
## Calculations for SV: chr16 27166424-chr3 89421185
## The x-axis coordinate for position 1 is 2427608641 = 2400442217 (the fake start of chr16) + 27166424 (the real position 1)
## The x-axis coordinate for position 2 is 580571136 = 491149951 (the fake start of chr3) + 89421185 (the real position 2).
## 
## Calculations for SV: chr19 15581754-chr12 11271304
## The x-axis coordinate for position 1 is 2669993042 = 2654411288 (the fake start of chr19) + 15581754 (the real position 1)
## The x-axis coordinate for position 2 is 1955038977 = 1943767673 (the fake start of chr12) + 11271304 (the real position 2).
## 
## Calculations for SV: chr19 19107305-chr7 26951385
## The x-axis coordinate for position 1 is 2673518593 = 2654411288 (the fake start of chr19) + 19107305 (the real position 1)
## The x-axis coordinate for position 2 is 1258955688 = 1232004303 (the fake start of chr7) + 26951385 (the real position 2).
## 
## Calculations for SV: chr22 17701854-chrX 28651228
## The x-axis coordinate for position 1 is 2841884908 = 2824183054 (the fake start of chr22) + 17701854 (the real position 1)
## The x-axis coordinate for position 2 is 2903652750 = 2875001522 (the fake start of chrX) + 28651228 (the real position 2).
## 
## Calculations for SV: chr22 32216382-chr7 47593881
## The x-axis coordinate for position 1 is 2856399436 = 2824183054 (the fake start of chr22) + 32216382 (the real position 1)
## The x-axis coordinate for position 2 is 1279598184 = 1232004303 (the fake start of chr7) + 47593881 (the real position 2).
## 
## Calculations for SV: chrX 28651217-chr22 17701865
## The x-axis coordinate for position 1 is 2903652739 = 2875001522 (the fake start of chrX) + 28651217 (the real position 1)
## The x-axis coordinate for position 2 is 2841884919 = 2824183054 (the fake start of chr22) + 17701865 (the real position 2).
## 
## Calculations for SV: chrX 134327496-chr2 27815968
## The x-axis coordinate for position 1 is 3009329018 = 2875001522 (the fake start of chrX) + 134327496 (the real position 1)
## The x-axis coordinate for position 2 is 276772390 = 248956422 (the fake start of chr2) + 27815968 (the real position 2).
```

Get random number for the bezier curves' heigths and caluclate the middle point for each bezier curve


```r
beziers.height <- runif(nrow(sv.data), 1, 2)
beziers.mid <- unlist(pos1_fake)+(unlist(pos2_fake)-unlist(pos1_fake))/2
```

Create dataframe with beziers curves info

```r
beziers <- data.frame(
    x = c(rbind( unlist(pos1_fake), beziers.mid, unlist(pos2_fake) )),
    y = c(rbind( 0.2, beziers.height, 0.2 ) ),
    group = rep( paste( chr1, pos1, chr2, pos2, sep="_" ), each=3),
		svtype = rep( type, each=3)
)
```


Generate a bezier curves-like plot representing SVs


```r
ggplot() + geom_bezier(aes(x= x, y = y, group = group, color = svtype ), data = beziers, show.legend = TRUE, size = 0.4) +

		##### Remove default axes labels and grey backgroud
		theme(axis.title.x=element_blank(), axis.text.x= element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y= element_blank(), axis.ticks.y=element_blank(),
		##### ...and the grey backgroud
					panel.background = element_rect(fill = NA),
		##### ...change the legend parameters
					legend.title=element_text(size=8), legend.text=element_text(size=7), legend.key.size = unit(0.7,"line"), legend.key= element_blank(), legend.position = c(0.97,0.7) ) +

		##### Set the axes limits
		scale_x_continuous(limits = c(1, genome.length)) +
		scale_y_continuous(limits = c(0, 2)) +

		##### Add chromosomes boundaries
		geom_segment(aes(x = c(1,unlist(chrs_fake_starts)[2:24],genome.length) , xend = c(1,unlist(chrs_fake_starts)[2:24],genome.length), y = 0, yend = 0.2), colour = 'grey', size = 0.2) +

		labs( color = "SV type") +

		##### Add chromosomes labels
		annotate(geom = 'text', label = names(chrs_fake_label.pos), x = unlist(chrs_fake_label.pos), y = 0.1, size = 2, angle = 45)
```

![](SVbezierPlot_files/figure-html/beziers_plot-1.png)<!-- -->

Print session info


```r
sessionInfo
```

```
## function (package = NULL) 
## {
##     z <- list()
##     z$R.version <- R.Version()
##     z$platform <- z$R.version$platform
##     if (nzchar(.Platform$r_arch)) 
##         z$platform <- paste(z$platform, .Platform$r_arch, sep = "/")
##     z$platform <- paste0(z$platform, " (", 8 * .Machine$sizeof.pointer, 
##         "-bit)")
##     z$locale <- Sys.getlocale()
##     if (.Platform$OS.type == "windows") {
##         z$running <- win.version()
##     }
##     else if (nzchar(Sys.which("uname"))) {
##         uname <- system("uname -a", intern = TRUE)
##         os <- sub(" .*", "", uname)
##         z$running <- switch(os, Linux = if (file.exists("/etc/os-release")) {
##             tmp <- readLines("/etc/os-release")
##             t2 <- if (any(startsWith(tmp, "PRETTY_NAME="))) sub("^PRETTY_NAME=", 
##                 "", grep("^PRETTY_NAME=", tmp, value = TRUE)[1L]) else if (any(startsWith(tmp, 
##                 "NAME"))) sub("^NAME=", "", grep("^NAME=", tmp, 
##                 value = TRUE)[1L]) else "Linux (unknown distro)"
##             sub("\"(.*)\"", "\\1", t2)
##         } else if (file.exists("/etc/system-release")) {
##             readLines("/etc/system-release")
##         }, Darwin = {
##             ver <- readLines("/System/Library/CoreServices/SystemVersion.plist")
##             ind <- grep("ProductUserVisibleVersion", ver)
##             ver <- ver[ind + 1L]
##             ver <- sub(".*<string>", "", ver)
##             ver <- sub("</string>$", "", ver)
##             ver1 <- strsplit(ver, ".", fixed = TRUE)[[1L]][2L]
##             sprintf("%s %s %s", ifelse(as.numeric(ver1) < 12, 
##                 "OS X", "macOS"), switch(ver1, `6` = "Snow Leopard", 
##                 `7` = "Lion", `8` = "Mountain Lion", `9` = "Mavericks", 
##                 `10` = "Yosemite", `11` = "El Capitan", `12` = "Sierra", 
##                 `13` = "High Sierra", ""), ver)
##         }, SunOS = {
##             ver <- system("uname -r", intern = TRUE)
##             paste("Solaris", strsplit(ver, ".", fixed = TRUE)[[1L]][2L])
##         }, uname)
##     }
##     if (is.null(package)) {
##         package <- grep("^package:", search(), value = TRUE)
##         keep <- vapply(package, function(x) x == "package:base" || 
##             !is.null(attr(as.environment(x), "path")), NA)
##         package <- .rmpkg(package[keep])
##     }
##     pkgDesc <- lapply(package, packageDescription, encoding = NA)
##     if (length(package) == 0) 
##         stop("no valid packages were specified")
##     basePkgs <- sapply(pkgDesc, function(x) !is.null(x$Priority) && 
##         x$Priority == "base")
##     z$basePkgs <- package[basePkgs]
##     if (any(!basePkgs)) {
##         z$otherPkgs <- pkgDesc[!basePkgs]
##         names(z$otherPkgs) <- package[!basePkgs]
##     }
##     loadedOnly <- loadedNamespaces()
##     loadedOnly <- loadedOnly[!(loadedOnly %in% package)]
##     if (length(loadedOnly)) {
##         names(loadedOnly) <- loadedOnly
##         pkgDesc <- c(pkgDesc, lapply(loadedOnly, packageDescription))
##         z$loadedOnly <- pkgDesc[loadedOnly]
##     }
##     z$matprod <- as.character(options("matprod"))
##     es <- extSoftVersion()
##     z$BLAS <- as.character(es["BLAS"])
##     z$LAPACK <- La_library()
##     class(z) <- "sessionInfo"
##     z
## }
## <bytecode: 0x7f851b1ff490>
## <environment: namespace:utils>
```
