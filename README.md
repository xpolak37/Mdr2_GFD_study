# Mdr2 mice GFD study

<div style="font-size: larger;">
Monika Cahová <sup>1</sup>, Ondřej Fabián <sup>2,3,4</sup>, Michal Kraus <sup>5</sup>, Martin Leníček <sup>6</sup>, Marie Heczková <sup>1</sup>, Veronika Motúzová <sup>5</sup>, Petra Polakovičová <sup>1,7</sup>, Miriam Brátová <sup>1</sup>, Alena Bohdanecká <sup>1,8</sup>, Aleksandra Bartelik <sup>5,6</sup>, Vanda Filová <sup>9</sup>, Mojmír Hlavatý <sup>8,10</sup>, Pavel Drastich <sup>10</sup>, Miloslav Kverka <sup>5</sup>
</div>

<sup>1</sup> Experimental Medicine Centre, Institute for Clinical and Experimental Medicine (IKEM), Prague, Czech Republic<br>
<sup>2</sup> Clinical and Transplant Pathology Centre, Institute for Clinical and Experimental Medicine (IKEM), Prague, Czech Republic<br>
<sup>3</sup> Department of Pathology and Molecular Medicine, Third Faculty of Medicine, Charles University and Thomayer Hospital, Prague, Czech Republic<br>
<sup>4</sup> Department of Technical Disciplines in Health Care, University of Presov, Presov, Slovakia<br>
<sup>5</sup> Laboratory of Cellular and Molecular Immunology, Institute of Microbiology AS CR, Prague, Czech Republic<br>
<sup>6</sup> Institute of Medical Biochemistry and Laboratory Diagnostics, First Faculty of Medicine and General University Hospital in Prague, Charles University, Prague, Czech Republic<br>
<sup>7</sup> Faculty of Science, Charles University, Prague, Czech Republic<br>
<sup>8</sup> First Faculty of Medicine, Charles University, Prague, Czech Republic<br>
<sup>9</sup> Laboratory Methods Division, Institute for Clinical and Experimental Medicine (IKEM), Prague, Czech Republic<br>
<sup>10</sup> Department of Hepatogastroenterology, Institute for Clinical and Experimental Medicine (IKEM), Prague, Czech Republic

<br>

---------------------------------------------------------------------------------------------------

**General information**

This repository provides a report of the study **Sex-dependent attenuation of sclerosing cholangitis by a gluten-free diet in Mdr2 knockout mice**, asking how **diet** (gluten-free diet, GFD, vs. standard diet, SD) and **sex** jointly modifies the course of cholestatic liver disease. 

All reported results can be reproduced using the code in this repository. Feel free to contact Petra Polakovicova by [petra.polakovicova@ikem.cz](mailto:petra.polakovicova@ikem.cz) or Monika Cahova [monika.cahova@ikem.cz](mailto:monika.cahova@ikem.cz)if you have any questions about the computational part of the study.

📚 **Citation**

The study has not been published yet. If you find this code and report helpful, please cite this repository; the reference will be updated after publication.

💾 **Data Availability**
The raw liver RNA sequencing and 16S rRNA amplicon sequencing reads have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number PRJEB129687.

----------------------------------------------------------------------------------------------------

## Report info

This project analyses liver RNA-seq data from two independent mouse experiments with a Diet × Sex design:

- **P76_3** (discovery cohort): 44 animals (SD: 10F/11M; GFD: 11F/12M), used to generate the hypotheses
- **P76_4** (validation cohort): 38 animals after QC (SD: 10F/9M; GFD: 10F/9M), used only to confirm or reject the P76_3 hypotheses

Liver expression is further linked to serum bile acids and clinical liver markers (AST, ALT, AlkP, total bilirubin).

📁 **Project Structure**

Below is an overview of the folder structure:

PSC_mice/   
├── rnaseq/                          
│   │   │   ├── custom_functions_rnaseq.R   
│   │   │   ├── rnaseq_P76_3_part1.qmd   (+ .html)  
│   │   │   ├── rnaseq_P76_3_part2.qmd   (+ .html)  
│   │   │   ├── rnaseq_P76_3_part3.qmd   (+ .html)  
│   │   │   └── validation_P76_4.qmd     (+ .html)  
│   
└── phenotype_microbiome/              
    ├── README.md                       
    ├── R/  
    │   ├── functions.R     
    │   ├── 01_liver_tests.R          
    │   ├── ...     
    │   └── 07_microbiome.R 

## Methodology

🖥️ **Bioinformatics processing**

Sequencing reads were processed with a custom pipeline available at [xpolak37/RNAseq_workflow](https://github.com/xpolak37/RNAseq_workflow) and  [xpolak37/16S-rRNA-amplicon-workflow](https://github.com/xpolak37/16S-rRNA-amplicon-workflow).

📊 **Statistical analysis**

All endpoints were analysed with a single model on a four-group factor (F_SD, F_GFD, M_SD, M_GFD) with common residual variance. Four contrasts were pre-specified as primary: the diet effect (GFD vs SD) in females and in males, and the sex effect on SD and on GFD. Multiplicity was controlled by the Holm procedure, applied separately to the diet and the sex contrast families. The diet-by-sex interaction was reported per endpoint as a secondary, uncorrected test; a pooled main-effects model was not used because the sex effect depended on diet. The exploratory cohort generated hypotheses and the validation cohort tested them: confirmatory hypotheses and cut points were pre-specified at the contrast level and frozen before the validation data were accessed, together with non-evaluability rules (for example, an endpoint was not evaluated when more than 50% of a group fell below the limit of quantification). For high-dimensional datasets, multiplicity was controlled by Benjamini-Hochberg instead of Holm. 

## Results

The code for transcriptomics with reported results can be found:
- [P76_3_part1](https://xpolak37.github.io/Mdr2_GFD_study/rnaseq/rnaseq_P76_3_part1.html): data import, quality control (detected genes, sequencing depth) and exploratory analysis
- [P76_3_part2](https://xpolak37.github.io/Mdr2_GFD_study/rnaseq/rnaseq_P76_3_part2.html): global structure (PCA, PERMANOVA, clustering), differential expression, collection-order checks and curated gene panels
- [P76_3_part3](https://xpolak37.github.io/Mdr2_GFD_study/rnaseq/rnaseq_P76_3_part3.html): integration with serum bile acids and clinical liver markers, and pathway enrichment (Hallmark, KEGG, Reactome)
- [P76_4_validation](https://xpolak37.github.io/Mdr2_GFD_study/rnaseq/rnaseq_validation_P76_4.html): tests whether the main P76_3 findings replicate in the independent P76_4 cohort

**Liver phenotype, immunology and microbiota**: the results of these analyses are reported in the manuscript. The R scripts that produce them are in [`phenotype_microbiome/R/`](phenotype_microbiome/R/), and the [`phenotype_microbiome/README.md`](phenotype_microbiome/README.md) explains how to run them.

---------------------------------------------------------------------------------------------------

## Acknowledgment

This study was supported by the Ministry of Health of the Czech Republic, grant no. NU22-06-00269, DRO – IKEM, IN 00023001, and by the Ministry of Education, Youth and Sports of the Czech Republic under the Programme Johannes Amos Comenius, project number CCZ.02.01.01/00/22_008/0004597 (INTER-MICRO).
