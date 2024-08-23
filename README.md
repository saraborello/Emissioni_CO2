# Predict CO2 Emission
The construction industry has a significant impact on the environment, contributing substantially to global CO2 emissions. Understanding the factors behind these emissions is crucial for formulating effective reduction strategies. Besides daily energy consumption, the design, construction, and demolition of buildings also play a role in CO2 emissions.

**OBJECTIVE**
The objective of this study is to develop a robust linear model, essential for its prospective use as a predictive tool on datasets not yet examined. The model is specifically aimed at discerning and quantifying the effect of multiple variables on the quantities of CO2 emitted by a residential building.

## Table of Contents
- [Exploratory Analysis](#Exploratory-Analysis)
   - [Dataset Description ](#Dataset-Description)
   - [Target Variable ](#Target-Variable)
   - [Selection of Variables](#Selection-of-Variables)
   - [Missing Data](#Missing-Data)
   - [Elimination of Problematic Observations](#Elimination-of-Problematic-Observations)
   - [Optimal Grouping and Transformations of X](#Optimal-Grouping-and-Transformations-of-X)


## Exploratory Analysis
### Dataset Description 
The dataset used in this study was obtained from the OpenData Portal of the Lombardy Region, specifically from the CENED Database, which focuses on the Energy Certification of Buildings. This data archive, comprising Energy Performance Certificates, consists of 1.52 million observations organized into 40 variables. Due to the large volume of the dataset, systematic sampling was performed, reducing the number of observations to 202,019 units. The investigation was subsequently focused exclusively on buildings with continuous residential use (variable ` DESTINAZIONE_DI_USO = E.1(1)`  as per DPR 412/1993), thus excluding those for public use, which led to a further reduction in observations, stabilizing at a total of 16,858 units.

### Target Variable 
The selected target variable represents CO2 emissions, quantified annually. Specifically, this variable is measured in KgCO2eq/m2year, providing an indicator of greenhouse gas emissions per unit of surface area. This metric facilitates the comparison between buildings of different sizes. From Figure 1, it is observed that emissions are concentrated between 0 and 100 KgCO2eq/m2year.

<div align="center">
  
  <img width="628" alt="Target variable" src="https://github.com/user-attachments/assets/9ecc60ac-9d27-4573-8190-6d8bef98f57d">
  
  Figure 1: Target Variable Distribution

</div>



### Selection of Variables 
During the model development phase, the entire range of 40 available variables was not utilized; instead, the selection was limited to 14 variables deemed most relevant, with their descriptions detailed in the Annex section. A graphical analysis of these variables in relation to the dependent variable (Figures 2 and 3) revealed anomalies suggesting the presence of input errors, likely due to the excessively broad scale. These observations will be eliminated in the next steps.

<div align="center">
  
  <img width="695" alt="Distribution of Cov vs CO2 1 " src="https://github.com/user-attachments/assets/df201508-4c85-4d49-8d41-05165231b955">

  
  Figure 2: Distribution of Covariates vs. CO2 Emission Variable

  <img width="723" alt="Distribution of Cov vs CO2 " src="https://github.com/user-attachments/assets/c49873d6-8198-445b-a595-71a1cd6cbdee">

  Figure 3: Distribution of Covariates vs. CO2 Emission Variable
</div>

### Missing Data
Before proceeding with the actual analysis, a phase is required where the dataset needs to be cleaned by selecting only the observations that can effectively be used in the analysis.
From the examination of Figure 4, it is observed that the variables `MOTIVAZIONE_APE` and `SUPERFICIE_VETRATA_OPACA` have few missing values, 2 and 4 respectively, suggesting a possible random absence of data; therefore, these observations will be excluded.
<div align="center">
  <img width="596" alt="Missing Data" src="https://github.com/user-attachments/assets/e5f8bd75-e2d9-4b22-a607-62cc440da593">
  
  Figure 4: Missing Data


</div>

### Elimination of Problematic Observations
Initial graphical analyses revealed that the dataset contains some anomalous values. These are presumably due to errors during the data entry process by users. Tables 1 and 2 list the instances where such irregularities occur; specifically, for the variable EFER, which measures the energy contribution of renewable energy systems, negative values were detected, which are not plausible. Similarly, the CO2 variable shows incorrect orders of magnitude.

<div align="center">
  
  <img width="709" alt="Problematic Observation" src="https://github.com/user-attachments/assets/24a3d4e1-fcd3-4e77-b86b-5a4a8f63edbd">


</div>

### Optimal Grouping and Transformations of X
After addressing the missing values, before proceeding with the analysis, it is crucial to make adjustments to certain categorical variables, as they have an excessive number of levels. The variables in question are `CLASSE_ENERGETICA`, `PROVINCIA`, and `MOTIVAZIONE_APE`, as shown in Figure 5.

<div align="center">
  
  <img width="627" alt="Covariates before grouping" src="https://github.com/user-attachments/assets/129585cd-f983-4ecf-b5aa-71e274257476">
  
  Figure 5: Covariates before grouping
  
</div>


The variable `PROVINCIA`, shown in the central graph of Figure 5, consists of 12 levels. Therefore, it was decided to categorize it into 3 distinct groups based on population density, as this can better reflect the intensity of energy use and the related infrastructure needs. The new partition is shown in Figure 6.

<div align="center">
  
  <img width="673" alt="Provinces by population density" src="https://github.com/user-attachments/assets/94a20fc2-186f-4135-a5f4-9af51fc1f2fe">
  
  Figure 6: Provinces by population density
  
</div>

The variable `CLASSE_ENERGETICA`, initially divided into 8 levels as shown in the right graph of Figure 5, was reorganized into 3 distinct groups as shown in Table 3.

<div align="center">

| Group | Description   | Count |
|-------|---------------|-------|
| 1     | High Class     | 332   |
| 2     | Medium Class   | 3,678 |
| 3     | Low Class      | 12,834 |


Table 3: Energy Class

</div>

The variable `MOTIVAZIONE_APE` had 14 levels, so the optimal grouping methodology was adopted to aggregate similar categories based on their similarity in the mean with respect to the target variable, reducing unnecessary levels. Following this process, 6 clusters were identified, visible in the dendrogram shown in Figure 7 and detailed in the Table.

<div align="center">
  
  <img width="666" alt="Optimal Grouping" src="https://github.com/user-attachments/assets/ba688f8f-1be0-4534-a6df-cfc38dd7e4fa">

  Figure 7: Optimal Grouping
  
</div>


<div align="center">
   
| Group     | Description                             | Count |
|-----------|-----------------------------------------|-----------|
| 1         | New constructions                       | 1,017     |
| 2         | Major structural interventions          | 487       |
| 3         | Energy renovation                       | 950       |
| 4         | Other                                   | 2,219     |
| 5         | Transfer for consideration              | 6,976     |
| 6         | Lease agreement                         | 5,195     |

Table 4: Reasons for Opening the APE Procedure

</div>


Figure 8 shows that the optimal grouping was effectively carried out and that there is a noticeable difference between the groups with respect to the dependent variable.

<div align="center">
  
  <img width="635" alt="Covariates after grouping" src="https://github.com/user-attachments/assets/ac7928d7-5721-40f1-9e02-42a0138705d9">

  Figure 8: Covariates after grouping
  
</div>

After addressing the missing values, correcting typographical errors, and recoding the categorical variables to have a manageable number of levels, the pre-processing phase of the dataset is concluded. This conclusion allows for the first estimation of the model to proceed. Specifically, Figure 9 clearly illustrates the distribution of continuous covariates in relation to the target, offering a more intelligible representation compared to the one presented in the previous section. It is observed that some non-linear relationships exist between the covariates and the dependent variable, which could present challenges related to the assumption of linearity within the model.

<div align="center">
  
<img width="654" alt="Covariates after pre-processing" src="https://github.com/user-attachments/assets/80524499-c8be-42b7-ba6e-6c38373d6ed6">

  Figure 9: Covariates after pre-processing
  
</div>
