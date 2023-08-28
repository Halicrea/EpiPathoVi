#!/usr/bin/env python3
#########################################
## 26/04/2023
## Par Elyna Bouchereau
## Fichier: Filter_PacBio_Data.py
###########################################
import pandas as ours

df_PACBIO = ours.read_csv(r'../Differences_epi_12-016_02-041_07-115/Methylation_Vaestu_alt.csv')
print(df_PACBIO)
df_12_016 = ours.DataFrame(df_PACBIO, columns=['seqid','position','context','type','strand','motif','statut','score_016','coverage_016','IPDRatio_016','fraction_016'])
print(df_12_016)

methyl_GATC_016 = ours.DataFrame()