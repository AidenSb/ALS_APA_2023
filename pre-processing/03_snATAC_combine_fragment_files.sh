#create single fragment file for ease of barcode/celltype/disease type selection in future analyses/visualizations

mkdir /ALS_snATAC/fragments

cp /ALS_snATAC/CTRL1_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/CTRL1_fragments.tsv.gz
cp /ALS_snATAC/CTRL2_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/CTRL2_fragments.tsv.gz
cp /ALS_snATAC/CTRL3_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/CTRL3_fragments.tsv.gz
cp /ALS_snATAC/CTRL4_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/CTRL4_fragments.tsv.gz
cp /ALS_snATAC/C9ALSFTLD1_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/C9ALSFTLD1_fragments.tsv.gz
cp /ALS_snATAC/C9ALSFTLD2_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/C9ALSFTLD2_fragments.tsv.gz
cp /ALS_snATAC/C9ALSFTLD3_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/C9ALSFTLD3_fragments.tsv.gz
cp /ALS_snATAC/C9ALSFTLD5_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/C9ALSFTLD5_fragments.tsv.gz
cp /ALS_snATAC/C9ALSFTLD6_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/C9ALSFTLD6_fragments.tsv.gz
cp /ALS_snATAC/C9ALSFTLD7_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/C9ALSFTLD7_fragments.tsv.gz
cp /ALS_snATAC/C9ALSnoFTLD1_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/C9ALSnoFTLD1_fragments.tsv.gz
cp /ALS_snATAC/C9ALSnoFTLD2_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/C9ALSnoFTLD2_fragments.tsv.gz
cp /ALS_snATAC/C9ALSnoFTLD3_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/C9ALSnoFTLD3_fragments.tsv.gz
cp /ALS_snATAC/sALSnoFTLD1_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/sALSnoFTLD1_fragments.tsv.gz
cp /ALS_snATAC/sALSnoFTLD2_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/sALSnoFTLD2_fragments.tsv.gz
cp /ALS_snATAC/sALSnoFTLD3_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/sALSnoFTLD3_fragments.tsv.gz
cp /ALS_snATAC/sALSnoFTLD4_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/sALSnoFTLD4_fragments.tsv.gz
cp /ALS_snATAC/sALSnoFTLD5_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/sALSnoFTLD5_fragments.tsv.gz
cp /ALS_snATAC/sALSnoFTLD6_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/sALSnoFTLD6_fragments.tsv.gz
cp /ALS_snATAC/sALSnoFTLD7_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/sALSnoFTLD7_fragments.tsv.gz
cp /ALS_snATAC/sALSnoFTLD8_snATAC/fragments.tsv.gz /ALS_snATAC/fragments/sALSnoFTLD8_fragments.tsv.gz

# decompress files and add the same cell prefix as was added to the seurat object

cd /ALS_snATAC/fragments

gzip -dc "CTRL1_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"CTRL1_"$4,$5}' - > CTRL1_fragments.tsv
gzip -dc "CTRL2_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"CTRL2_"$4,$5}' - > CTRL2_fragments.tsv
gzip -dc "CTRL3_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"CTRL3_"$4,$5}' - > CTRL3_fragments.tsv
gzip -dc "CTRL4_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"CTRL4_"$4,$5}' - > CTRL4_fragments.tsv

gzip -dc "C9ALSFTLD1_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"C9ALSFTLD1_"$4,$5}' - > C9ALSFTLD1_fragments.tsv
gzip -dc "C9ALSFTLD2_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"C9ALSFTLD2_"$4,$5}' - > C9ALSFTLD2_fragments.tsv
gzip -dc "C9ALSFTLD3_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"C9ALSFTLD3_"$4,$5}' - > C9ALSFTLD3_fragments.tsv
gzip -dc "C9ALSFTLD5_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"C9ALSFTLD5_"$4,$5}' - > C9ALSFTLD5_fragments.tsv
gzip -dc "C9ALSFTLD6_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"C9ALSFTLD6_"$4,$5}' - > C9ALSFTLD6_fragments.tsv
gzip -dc "C9ALSFTLD7_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"C9ALSFTLD7_"$4,$5}' - > C9ALSFTLD7_fragments.tsv

gzip -dc "C9ALSnoFTLD1_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"C9ALSnoFTLD1_"$4,$5}' - > C9ALSnoFTLD1_fragments.tsv
gzip -dc "C9ALSnoFTLD2_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"C9ALSnoFTLD2_"$4,$5}' - > C9ALSnoFTLD2_fragments.tsv
gzip -dc "C9ALSnoFTLD3_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"C9ALSnoFTLD3_"$4,$5}' - > C9ALSnoFTLD3_fragments.tsv

gzip -dc "sALSnoFTLD1_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"sALSnoFTLD1_"$4,$5}' - > sALSnoFTLD1_fragments.tsv
gzip -dc "sALSnoFTLD2_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"sALSnoFTLD2_"$4,$5}' - > sALSnoFTLD2_fragments.tsv
gzip -dc "sALSnoFTLD3_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"sALSnoFTLD3_"$4,$5}' - > sALSnoFTLD3_fragments.tsv
gzip -dc "sALSnoFTLD4_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"sALSnoFTLD4_"$4,$5}' - > sALSnoFTLD4_fragments.tsv
gzip -dc "sALSnoFTLD5_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"sALSnoFTLD5_"$4,$5}' - > sALSnoFTLD5_fragments.tsv
gzip -dc "sALSnoFTLD6_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"sALSnoFTLD6_"$4,$5}' - > sALSnoFTLD6_fragments.tsv
gzip -dc "sALSnoFTLD7_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"sALSnoFTLD7_"$4,$5}' - > sALSnoFTLD7_fragments.tsv
gzip -dc "sALSnoFTLD8_fragments.tsv.gz" | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"sALSnoFTLD8_"$4,$5}' - > sALSnoFTLD8_fragments.tsv

# merge files (avoids having to re-sort)
sort -k1,1V -k2,2n -k3,3n CTRL1_fragments.tsv CTRL2_fragments.tsv CTRL3_fragments.tsv CTRL4_fragments.tsv C9ALSFTLD1_fragments.tsv C9ALSFTLD2_fragments.tsv C9ALSFTLD3_fragments.tsv C9ALSFTLD5_fragments.tsv C9ALSFTLD6_fragments.tsv C9ALSFTLD7_fragments.tsv C9ALSnoFTLD1_fragments.tsv C9ALSnoFTLD2_fragments.tsv C9ALSnoFTLD3_fragments.tsv sALSnoFTLD1_fragments.tsv sALSnoFTLD2_fragments.tsv sALSnoFTLD3_fragments.tsv sALSnoFTLD4_fragments.tsv sALSnoFTLD5_fragments.tsv sALSnoFTLD6_fragments.tsv sALSnoFTLD7_fragments.tsv sALSnoFTLD8_fragments.tsv > fragments.tsv

# block gzip compress the merged file (use -@ # for # threads)
bgzip -@ 16 fragments.tsv 

# index the bgzipped file
tabix -p bed fragments.tsv.gz


