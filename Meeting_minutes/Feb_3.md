# Agenda
- Think of research questions (likely to do with either the MICB 475 parkisons datasets or maxines newly approved dataset)
    - Gastric Cancer and Parkinsons connection (https://pmc.ncbi.nlm.nih.gov/articles/PMC7950232/#r6)
    - Found another dataset : Gut microbiome dysbiosis across early Parkinson’s disease, REM sleep behavior disorder and their first-degree relatives
(https://pmc.ncbi.nlm.nih.gov/articles/PMC10154387/#notes3)
- Discuss the best research question for our project
- Outline next steps for individuals once our research question has been confirmed

#### Barcode paper ####
- determined 15 DEGs between gastric cancer (GC) and Parkinson's disease (PD)
- 3 gene ontology (GO) biological process terms common
  - Positive regulation of proteasomal ubiquitin-dependent protein catabolic process
  - Myeloid cell differentiation
  - Leukocyte differentiation

#### Potential Research Question ####
How does gut microbiome composition and diet‑related microbial dysbiosis modulate the expression of differentially expressed genes and the shared gene‑ontology biological processes that are common between Parkinson’s disease and gastric cancer?


### Meeting Notes ###
- we found a paper that showed a potential link between gene expressions to analyze the comparison between the gut microbiome of in Parkinson's and Gastric Conditions
- we want to profile parkinsons, gastric cancer, and healthy patients
- Another Idea: could potentially also look at patients with higher BMI and REM sleep in Parkinson's patients

# Agenda
- Think of research questions (likely to do with either the MICB 475 parkisons datasets or maxines newly approved dataset)
    - Gastric Cancer and Parkinsons connection (https://pmc.ncbi.nlm.nih.gov/articles/PMC7950232/#r6)
    - Found another dataset : Gut microbiome dysbiosis across early Parkinson’s disease, REM sleep behavior disorder and their first-degree relatives
(https://pmc.ncbi.nlm.nih.gov/articles/PMC10154387/#notes3)
- Discuss the best research question for our project
- Outline next steps for individuals once our research question has been confirmed

#### Barcode paper ####
- determined 15 DEGs between gastric cancer (GC) and Parkinson's disease (PD)
- 3 gene ontology (GO) biological process terms common
  - Positive regulation of proteasomal ubiquitin-dependent protein catabolic process
  - Myeloid cell differentiation
  - Leukocyte differentiation

#### Potential Research Question ####
How does gut microbiome composition and diet‑related microbial dysbiosis modulate the expression of differentially expressed genes and the shared gene‑ontology biological processes that are common between Parkinson’s disease and gastric cancer?


### Meeting Notes ###
-Combining parkinsons and gastric cancer datasets and looking for connections between those and the microbiomes of healthy and unhelthy patients could be a strong potential direction.

-Parkinsons and REM sleep disorder correlations from maxines paper could be a potential topic (good metadata but the groupings of patients with or without disease is complex)
    -4 groups for this data set: 
    control (no diesease, n=130), 
    RBD-relative (control for people with live with people who have REM diseasen=132), 
    RBD (n=175), 
    early-PD (n=36)
-beta diversity control vs PD was significant
-correlation analysis on abundance: slight effects
-seems as if they were very comprehensive in the paper already

-Determined we will stick with the PD and GC datasets

-merging data: -what variable region?
               -single end sequencing or paried?
               -if we do them in parrallel; merge the table and rep-seq following individual import and denoising

-compare gc patients with gc controls and parkinsons patients with parkinons controls, do not compare controls

analysis:
1.Diversity metrics
2.Core microbiome (4-way ven diagram)
3.Indicator taxa
4.Deseq
5.Functional analysis (machone learing? depends on results)

Research Question: Are there any similarities of the microbiome of GC and PD patients? Which pathways are prevalent here?

-Send Ritu an email before trimming
-Save files to github along the way
-We can run both sets in parallel 

Johny: Import Parkinsons
Mixine: Find papers
Cadan: Start proposal introduction
Sina: Import GC
Aryan: Finding papers
