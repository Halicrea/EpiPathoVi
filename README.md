# EpiVib

## Description
Au cours des dernières décennies, la fréquence d’apparition de maladies émergentes (MIE) a augmenté. Il a été montré que l’activité humaine, la destruction des habitats, les changements climatiques et le transport d’animaux sont d’importants facteurs favorisant leurs apparitions. Une de ces maladies émergentes affecte drastiquement les huîtres adultes et est causée par une bactérie qui a émergé récemment comme pathogène : V. aestuarianus. Les mécanismes sous-jacents sa virulence restent méconnus et nécessitent d’être étudiés. Au cours de ce stage, nous avons étudié les régulations transcriptomiques et épigénétiques de Va 12/016 in vitro dans différentes conditions.   

Ce projet à pour but de permettre de combler le manque d'outils sur l'analyse de l'épigénome bactérien en combinant des outils existant et développant des outils et des interfaces afin de faire le lien entre chaque parties. L'objectif est de pouvoir rendre l'analyse de l'épigénome accessible et intégrative afin de permettre à un utilisateur quelquonque d'identifier des sites de régulations épigénétique pouvant impacter des clusters géniques.

## Badges
Multilingual scripts

## Visuals


## Installation
Yet to be implemented.

## Usage
In construction

## Support
elyna.bouchereau@ifremer.fr

## Roadmap
- [x] Traitement de données de séquençage PacBio (SMRT-seq) et obtenir les coordonnées génomiques des marques de méthylations.
- [x] Réaliser une visualisation du méthylome avec circos.
- [x] Listé et représenté les régions hypométhylées par rapport à une condition étudiée.
- [x] Analyse de régulations transcriptomique à partir de séquençage illumina.
- [x] Enrichir les résultats épigénétiques avec les résultats transcriptomique. (Reste encore manuelle)
- [ ] Créer une interface facilitant les différentes étapes.
- [ ] Installation automatique de l'application et des dépendances.

## Ramification
- 0_bin : exectubles binaires (c++)
- 1_obj : fichiers .o pour les dépendances c++ (vide pour l'instant) 
- 2_src : *scripts* shell, R, python et c++
- - Analyse_RNAseq.Rmd : Analyse RNA-seq
- - clustering_basemode.R : Visualisation des données brutes en sortie du pipeline de traitement SMRT-seq
- - preprocess_SMRT_OUT.R : Combiner les tableaux pour identifier les différentes de méthylation entre les conditions expérimentales
- - rna-seq-test.pbs : Pour visualiser la qualité des *reads* illumina
- - run-sortmerna.pbs : Identification et nettoyage des ARN ribosomaux
- - run-trimming.pbs : Nettoyage des *reads* illumina
- - run-Bowtie2-DB-Idx.pbs : Indexage pour le génome de référence utilisé
- - run-mapping.pbs : Alignement des *reads* sur le génome de référence
- - run-htseqcount.pbs : Quantification des transcripts avec Salmon
- - pipeline_smrt_seq.pbs : Traitement des données de SMRT-seq
- - search_motif.py : Identification des motifs assoicé à des marques de méthylation dans un génome donné
- - Filter_PacBio_Data.R : Adaptation du script précédent et *wrapper* autour de l'algorithme en c++ prenant en compte le contexte des motifs
- - Finding_context_of_motif.cpp : Identificationd des coordonnées génomique de motifs dans un génome donné.
- R_results/ : Résultats des scripts R
- Transcriptomic/ : statistiques du pipeline d'analyse RNA-seq
- Epigenetic/ : Visualisation circos (s'utilise en ligne de commande: circos --conf bacteria.circos.conf)
- Differences_epi_12-016_02-041_07-115 : Résultats des données traitées par Ludovic LEGRAND

## Authors and acknowledgment
Elyna BOUCHEREAU

## License

## Project status
En développement.
