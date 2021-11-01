# PICCAgeing
Devine, Neumann, Levari, Wilson, &amp; Eppinger (2021), Human Ageing is Associated with More Rigid Concept Spaces. Data, task code, analysis scripts, and model code. 

## Directories
``` 
.
├── Analysis.R
├── data
│   ├── Dots
│   └── Ethics
├── Models
│   ├── HDDM
│   ├── DDM
│   └── seq
├── PICCOA_2021.RData
├── README.md
└── TaskCode
    ├── EthicPICC.py
    ├── PICCcolour.py
    └── stimuli


```

### **data**
Contains subject data. Id # < 100 is OA, > 100 is YA. Dots and Ethics task data are saved in their own directories (as .csv)

### **Models**
Contains model code and individual model fits for both the HDDM (`HDDM`, Ratcliff & McKoon, 2008; Wiecki et al., 2013) and sequential decision-making model (`seq`, Wilson, 2018). `DDM` contains code and fits for a DDM model estimated with MLE used in a previous draft of the manuscript

### **TaskCode**
PsychoPy scrips for the the experiment. PICCcolour.py is the Dots Task. EthicPICC.py is the Ethics task. /stimuli/ contains all the stimuli used in the Ethics task--i.e., the research scenarios. 

## Standalone Files
* Analysis.R is the main analysis script
* PICCOA_2021.RData is the R environment for an executed Analysis.R
* README.md is this file


