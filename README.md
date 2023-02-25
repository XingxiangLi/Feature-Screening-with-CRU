# Feature Screening with Conditional Rank Utility for Big-data Classification

## What is included?

This folder contains R code and data to reproduce the numerical results presented in the 
paper (i.e. Figures 1-2 and Tables 2-7). Note that Table 1 in the paper is purely descriptive 
and is not a result of any numerical experiment.

There are six items in the main folder, including two files and four subfolders. The Two files 
are: "proj_file.Rproj" and "readme.txt". The four subfolders are: "scripts", "data", "figures" and "tables".
  
Details of each item are described below. 
  
  - "proj_file.Rproj": a project file for opening the project with Rstudio. 

  - "readme.txt": documentation on workflow and related information about the replication code.
  
  - subfolder "scripts" contains the following files:
      - "CRU_function.R": a helper R file containing functions needed to run the replication code; 
      - "CRU_ReplicationCode_Full.R": full replication code to reproduce the numerical results;
      - "CRU_ReplicationCode_ShortRunTime.R": replication code with short runtime for testing purposes.

  - subfolder "data" contains the following 6 data files (feature sets) that are analyzed in Section 4.2; 
    see Data section of this readme file for more information about these data files.
     - fac.txt  
     - fou.txt  
     - kar.txt  
     - mor.txt  
     - pix.txt  
     - zer.txt  
  
  - subfolder "figures" is a placeholder to store the two figure outputs after running the replication code.
    The two figure outputs (in .pdf format) correspond to Figures 1-2 in the paper.
 
  - subfolder "tables" is a placeholder to store the six table outputs after running the replication code.
    The six table outputs (in .xlsx format) correspond to Tables 2-7 in the paper.
