# AMULET: ATAC-seq MULtiplet Estimation Tool (https://ucarlab.com/AMULET/)
#format for: AMULET.sh /path/to/fragments.tsv.gz /path/to/singlecell.csv /path/to/human_autosomes.txt /path/to/repeatfilter.bed /path/to/output/ /path/to/shellscript/
#from /path/to/ALS_snATAC

mkdir CTRL1_snATAC_AMULET
mkdir CTRL2_snATAC_AMULET
mkdir CTRL3_snATAC_AMULET
mkdir CTRL4_snATAC_AMULET
mkdir C9ALSFTLD1_snATAC_AMULET
mkdir C9ALSFTLD2_snATAC_AMULET
mkdir C9ALSFTLD3_snATAC_AMULET
mkdir C9ALSFTLD5_snATAC_AMULET
mkdir C9ALSFTLD6_snATAC_AMULET
mkdir C9ALSFTLD7_snATAC_AMULET
mkdir C9ALSnoFTLD1_snATAC_AMULET
mkdir C9ALSnoFTLD2_snATAC_AMULET
mkdir C9ALSnoFTLD3_snATAC_AMULET
mkdir sALSnoFTLD1_snATAC_AMULET
mkdir sALSnoFTLD2_snATAC_AMULET
mkdir sALSnoFTLD3_snATAC_AMULET
mkdir sALSnoFTLD4_snATAC_AMULET
mkdir sALSnoFTLD5_snATAC_AMULET
mkdir sALSnoFTLD6_snATAC_AMULET
mkdir sALSnoFTLD7_snATAC_AMULET
mkdir sALSnoFTLD8_snATAC_AMULET

./AMULET.sh /ALS_snATAC/CTRL1_snATAC/fragments.tsv.gz /ALS_snATAC/CTRL1_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/CTRL1_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/CTRL2_snATAC/fragments.tsv.gz /ALS_snATAC/CTRL2_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/CTRL2_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/CTRL3_snATAC/fragments.tsv.gz /ALS_snATAC/CTRL3_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/CTRL3_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/CTRL4_snATAC/fragments.tsv.gz /ALS_snATAC/CTRL4_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/CTRL4_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/C9ALSFTLD1_snATAC/fragments.tsv.gz /ALS_snATAC/C9ALSFTLD1_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/C9ALSFTLD1_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/C9ALSFTLD2_snATAC/fragments.tsv.gz /ALS_snATAC/C9ALSFTLD2_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/C9ALSFTLD2_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/C9ALSFTLD3_snATAC/fragments.tsv.gz /ALS_snATAC/C9ALSFTLD3_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/C9ALSFTLD3_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/C9ALSFTLD5_snATAC/fragments.tsv.gz /ALS_snATAC/C9ALSFTLD5_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/C9ALSFTLD5_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/C9ALSFTLD6_snATAC/fragments.tsv.gz /ALS_snATAC/C9ALSFTLD6_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/C9ALSFTLD6_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/C9ALSFTLD7_snATAC/fragments.tsv.gz /ALS_snATAC/C9ALSFTLD7_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/C9ALSFTLD7_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/C9ALSnoFTLD1_snATAC/fragments.tsv.gz /ALS_snATAC/C9ALSnoFTLD1_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/C9ALSnoFTLD1_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/C9ALSnoFTLD2_snATAC/fragments.tsv.gz /ALS_snATAC/C9ALSnoFTLD2_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/C9ALSnoFTLD2_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/C9ALSnoFTLD3_snATAC/fragments.tsv.gz /ALS_snATAC/C9ALSnoFTLD3_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/C9ALSnoFTLD3_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/sALSnoFTLD1_snATAC/fragments.tsv.gz /ALS_snATAC/sALSnoFTLD1_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/sALSnoFTLD1_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/sALSnoFTLD2_snATAC/fragments.tsv.gz /ALS_snATAC/sALSnoFTLD2_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/sALSnoFTLD2_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/sALSnoFTLD3_snATAC/fragments.tsv.gz /ALS_snATAC/sALSnoFTLD3_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/sALSnoFTLD3_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/sALSnoFTLD4_snATAC/fragments.tsv.gz /ALS_snATAC/sALSnoFTLD4_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/sALSnoFTLD4_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/sALSnoFTLD5_snATAC/fragments.tsv.gz /ALS_snATAC/sALSnoFTLD5_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/sALSnoFTLD5_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/sALSnoFTLD6_snATAC/fragments.tsv.gz /ALS_snATAC/sALSnoFTLD6_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/sALSnoFTLD6_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/sALSnoFTLD7_snATAC/fragments.tsv.gz /ALS_snATAC/sALSnoFTLD7_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/sALSnoFTLD7_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

./AMULET.sh /ALS_snATAC/sALSnoFTLD8_snATAC/fragments.tsv.gz /ALS_snATAC/sALSnoFTLD8_snATAC/singlecell.csv /ALS_snATAC/AMULET-v1.1/human_autosomes.txt /ALS_snATAC/AMULET-v1.1/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed /ALS_snATAC/sALSnoFTLD8_snATAC_AMULET/ /ALS_snATAC/AMULET-v1.1/

sort "/ALS_snATAC/CTRL1_snATAC_AMULET/CTRL1_MultipletCellIds.txt" "/ALS_snATAC/CTRL1_snATAC_AMULET/CTRL2_MultipletCellIds.txt" "/ALS_snATAC/CTRL3_snATAC_AMULET/CTRL3_MultipletCellIds.txt" "/ALS_snATAC/CTRL4_snATAC_AMULET/CTRL4_MultipletCellIds.txt" "/ALS_snATAC/C9ALSFTLD1_snATAC_AMULET/C9ALSFTLD1_MultipletCellIds.txt" "/ALS_snATAC/C9ALSFTLD2_snATAC_AMULET/C9ALSFTLD2_MultipletCellIds.txt" "/ALS_snATAC/C9ALSFTLD3_snATAC_AMULET/C9ALSFTLD3_MultipletCellIds.txt" "/ALS_snATAC/C9ALSFTLD5_snATAC_AMULET/C9ALSFTLD5_MultipletCellIds.txt" "/ALS_snATAC/C9ALSFTLD6_snATAC_AMULET/C9ALSFTLD6_MultipletCellIds.txt" "/ALS_snATAC/C9ALSFTLD7_snATAC_AMULET/C9ALSFTLD7_MultipletCellIds.txt" "/ALS_snATAC/C9ALSnoFTLD1_snATAC_AMULET/C9ALSnoFTLD1_MultipletCellIds.txt" "/ALS_snATAC/C9ALSnoFTLD2_snATAC_AMULET/C9ALSnoFTLD2_MultipletCellIds.txt" "/ALS_snATAC/C9ALSnoFTLD3_snATAC_AMULET/C9ALSnoFTLD3_MultipletCellIds.txt" "/ALS_snATAC/sALSnoFTLD1_snATAC_AMULET/sALSnoFTLD1_MultipletCellIds.txt" "/ALS_snATAC/sALSnoFTLD2_snATAC_AMULET/sALSnoFTLD2_MultipletCellIds.txt" "/ALS_snATAC/sALSnoFTLD3_snATAC_AMULET/sALSnoFTLD3_MultipletCellIds.txt" "/ALS_snATAC/sALSnoFTLD4_snATAC_AMULET/sALSnoFTLD4_MultipletCellIds.txt" "/ALS_snATAC/sALSnoFTLD4_snATAC_AMULET/sALSnoFTLD4_MultipletCellIds.txt" "/ALS_snATAC/sALSnoFTLD5_snATAC_AMULET/sALSnoFTLD5_MultipletCellIds.txt" "/ALS_snATAC/sALSnoFTLD6_snATAC_AMULET/sALSnoFTLD6_MultipletCellIds.txt" "/ALS_snATAC/sALSnoFTLD7_snATAC_AMULET/sALSnoFTLD7_MultipletCellIds.txt" "/ALS_snATAC/sALSnoFTLD8_snATAC_AMULET/sALSnoFTLD8_MultipletCellIds.txt" > All_MutipletCellIds.txt
