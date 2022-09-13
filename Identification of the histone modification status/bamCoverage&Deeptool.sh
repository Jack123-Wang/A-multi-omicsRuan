#step1:
#The RPGC (reads per genome coverage) of ten histone markers were normalized by using the bamCoverage tool in the deepTools software (3.1.3).

ls *.bam|while read id;do bamCoverage -b $id -o $(basename $id .merged.positionsort.bam)_RPGC.bw --effectiveGenomeSize 2652783500 --normalizeUsing  RPGC ;done

#step2:
#The degree of chromatin modification for each co-occupying binding site was calculated using the computeMatrix tool in the deepTools software. Visualization was achieved by using the plotHeatmap tool in the deepTools software.
ls *.bed|while read id;do computeMatrix reference-point -p 15 --referencePoint center -b 3000 -a 3000 -R $id -S GSE25532_H3K27me3_RPGC.bw GSE29413_H3K9me3_RPGC.bw  GSE12241_H4K20me3_RPGC.bw GSE29218_H3K4me3_RPGC.bw  GSE27827_H3K4me2_RPGC.bw GSE31284_H3K9ac_RPGC.bw GSE41589_H3K36me3_RPGC.bw GSE11724_H3K79me2_RPGC.bw  GSE30203_H3K4me1_RPGC.bw  GSE24164_H3K27ac_RPGC.bw --samplesLabel "H3K27me3" "H3K9me3" "H4K20me3" "H3K4me3" "H3K4me2" "H3K9ac" "H3K36me3" "H3K79me2" "H3K4me1" "H3K27ac" -o $(basename $id .bed).gz ; done
