#process snATAC-Seq fastq with Cell Ranger ATAC (version 2.0)

cellranger-atac count --id=CTRL1_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=CTRL1_snATAC &
wait;

cellranger-atac count --id=CTRL2_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=CTRL2_snATAC &
wait;

cellranger-atac count --id=CTRL3_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=CTRL3_snATAC &
wait;

cellranger-atac count --id=CTRL4_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=CTRL4_snATAC &
wait;

cellranger-atac count --id=C9ALSFTLD1_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=C9ALSFTLD1_snATAC &
wait;

cellranger-atac count --id=C9ALSFTLD2_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=C9ALSFTLD2_snATAC &
wait;

cellranger-atac count --id=C9ALSFTLD3_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=C9ALSFTLD3_snATAC &
wait;

cellranger-atac count --id=C9ALSFTLD5_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=C9ALSFTLD4_snATAC &
wait;

cellranger-atac count --id=C9ALSFTLD6_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=C9ALSFTLD1_snATAC &
wait;

cellranger-atac count --id=C9ALSFTLD7_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=C9ALSFTLD2_snATAC &
wait;

cellranger-atac count --id=C9ALSnoFTLD1_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=C9ALSnoFTLD1_snATAC &
wait;

cellranger-atac count --id=C9ALSnoFTLD2_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=C9ALSnoFTLD2_snATAC &
wait;

cellranger-atac count --id=C9ALSnoFTLD3_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=C9ALSnoFTLD3_snATAC &
wait;

cellranger-atac count --id=sALSnoFTLD1_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=sALSnoFTLD1_snATAC &
wait;

cellranger-atac count --id=sALSnoFTLD2_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=sALSnoFTLD2_snATAC &
wait;

cellranger-atac count --id=sALSnoFTLD3_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=sALSnoFTLD3_snATAC &
wait;

cellranger-atac count --id=sALSnoFTLD4_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=sALSnoFTLD4_snATAC &
wait;

cellranger-atac count --id=sALSnoFTLD5_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=sALSnoFTLD5_snATAC &
wait;

cellranger-atac count --id=sALSnoFTLD6_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=sALSnoFTLD6_snATAC &
wait;

cellranger-atac count --id=sALSnoFTLD7_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=sALSnoFTLD1_snATAC &
wait;

cellranger-atac count --id=sALSnoFTLD8_snATAC \
                      --reference=/ALS_snATAC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=/ALS_snATAC/fastq/ \
                      --sample=sALSnoFTLD2_snATAC &
wait;

