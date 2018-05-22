# Facial_masculinity_MHC
Data and scripts for "Facial masculinity does not appear to be a condition-dependent male ornament and does not reflect MHC heterozygosity in humans"

If you use any of the data or scripts for your own research, please cite the biorxiv paper:\
Zaidi AA, White JD, Mattern BC, Liebowitz CR, Puts DA, Claes P, Shriver MD. 2018. Facial masculinity does not appear to be a condition-dependent male ornament in humans and does not reflect MHC heterozygosity. bioRxiv [Internet]:322255. Available from: https://www.biorxiv.org/content/early/2018/05/17/322255


Details of folders:

## Scripts
  1. calc_masc_03262018.R : Script used to generate facial masculinity measures from 3D coorsinates data
  2. summary_stats_03262018.R : Script used to calculate summary statistics (e.g. degree of sexual dimorphism, cohen's D, and      Levene's test of equal variances) on facial masculinity
  3. masc_height_association_03262018.R : Script used to test for associations between facial masculinity and height                (Hypothesis 1). 
  4. height_hla_association_02272018.R : Script used to test for association between height and MHC heterozygosity 
     (Hypothesis 2).
  5. masc_hla_association_03262018.R : Script used to test for association between facial masculinity and MHC heterozygosity
     (Hypothesis 3).
  6. Plot1Face.R: Script used to visualize high-dimensional facial masculinity scores, as well as quasi-landmark based test        statistics (cohen's D, beta coefficients, p-values etc.) as facial heatmaps.
     
## Dataset
  1. _Euro_demographic_03292018.dat_\
  Dataframe containing demographic data for all 1,233 individuals used in the study. Columns are in the following order:\
  ... IID: Unique identifier for individual\
  ... phet_genome: Heterozygosity across genome-wide (LD-pruned) SNPs\
  ... phet_hla: Heterozygosity across MHC locus\
  ... Sex\
  ... Age\
  ... Height: In centimeters\
  ... Weight: In Kg\
  ... gPC1-4: genetic PC scores\
  
  2. _euro_1233_masc_het_03292018.dat_\
  Dataframe containing all the columns from Euro_demographic_03292018.dat and two additional columns:\
  a) avg.masc: Average facial masculinity calculated across all 7,150 QLs for each person\
  b) avg.masc.unit: Average facial masculinity scaled by the Euclidean distance between female and male consensus faces. 
     In other words, this is average facial masculinity divided by the average facial masculinity of the male consensus face.
     
  3. _qlmasc_1233_noprop_03292018.txt.zip_\
  Zipped text file containing high-dimensional facial masculinity scores (dimensions: 7,150 x 1,233) i.e. 3D facial masculinity calculated for 7,150 quasi-landmarks (rows) for each of 1,233 individuals (columns). Referred to in the paper as FM<sub>QL</sub>. No headers or indices included. The order of the individuals is the same as the order of individuals in Euro_demographic_03292018.dat file.
  
  4. _Refscan.obj_\
  .OBJ file containing the 3D facial template used to visualize facial heatmaps etc.
  
## Results

  #### Summary_dat
  
  1. _ql_gsd_03292018.txt_\
  ... dataframe of length 7,150 x 1. Each row is the absolute degree of sexual dimorphism for each quasi-landmark
  
  2. _ql_sex_cohenD_03292018.txt_\
  ... dataframe of length 7,150 x 1. Each row is the Cohen's D of sexual dimorphism for each quasi-landmark
  
  ### Masc_v_height
  Linear model results: facial masculinity ~ height + covariates
  
  1. _lm_overallmasc_height_results_03292018_\
  ... Results of linear model where overall facial masculinity was used as a response
  
  2. _qlmasc_height_b_noprop_sd_04142018.txt_\
  ... Beta coefficients of linear model where FM<sub>QL</sub> was used individually as a response
  
  3. _qlmasc_height_p_noprop_sd_04142018.txt_\
  ... p-value for beta coefficients
  
  4. _qlmasc_height_q_noprop_sd_04142018.txt_\
  ... indicator variable to indicate if p-value passes Bonferroni correction
  ...... 0 if p-value <0.05/7,150
  ...... 1 if p-value >= 0.05/7,150
  
  5. _qlmasc_h_bysex_noprop_sd_04142018.txt_\
  ... slope_female: beta coefficients for FM<sub>QL</sub> - females only
  ... slope_male: beta coefficients for FM<sub>QL</sub> - males only
  ... diff: slope_male - slope_female
  ... z: Z-score for difference in slope
  ... p.value: Two-sided p-value for Z-score
  ... q.value: indicator variable indicating if p-value passes Bonferroni correction
  


  
  

