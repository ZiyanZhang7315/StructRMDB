# StructRMDB

These files contain the simulation code for the **StructRMDB** database.

## Pre-RNA Pipeline

### 1. `preRNA_sequence.R`
**Input:**  
- Genome file  
- GTF file  
- TxDb file  

**Output:**  
- Modified and unmodified Pre RNA sequences for each site  
- Transcript name and genomic location information  

---

### 2. `Secondary_structure_prediction_rnastructure.sh` and `Secondary_structure_prediction_viennaRNA.sh`
**Input:**  
- Modified and unmodified sequences for each site  

**Output:**  
- Predicted secondary structures  
- Minimum Free Energy (MFE) values  
- RNAforester scores  
- Secondary structure CT files  

---

### 3. `RNAsmc_score_and_combine.R`
**Input:**  
- Predicted secondary structures with MFE and RNAforester scores  

**Output:**  
- RNAsmc score  
- Combined dataframe including all previous results  

---

### 4. `Secondary_structure_visualization.R`
**Input:**  
- Secondary structure CT file  
- Site location file  

**Output:**  
- Visualized image for each predicted sequence

## Mature-RNA Pipeline

### 1. `Exon_sequence.R`
**Input:**  
- Genome file  
- GTF file  
- TxDb file  

**Output:**  
- exon sequence
- exon information 

---

### 2. `matureRNA_sequence.R`
**Input:**  
- exon sequence
- exon information  

**Output:**  
- Modified and unmodified mature RNA sequences for each site  
- Transcript name and genomic location information  

---

### 3. `Secondary_structure_prediction_rnastructure.sh` and `Secondary_structure_prediction_viennaRNA.sh`
**Input:**  
- Modified and unmodified sequences for each site  

**Output:**  
- Predicted secondary structures  
- Minimum Free Energy (MFE) values  
- RNAforester scores  
- Secondary structure CT files  

---

### 4. `RNAsmc_score_and_combine.R`
**Input:**  
- Predicted secondary structures with MFE and RNAforester scores  

**Output:**  
- RNAsmc score  
- Combined dataframe including all previous results  

---

### 5. `Secondary_structure_visualization.R`
**Input:**  
- Secondary structure CT file  
- Site location file  

**Output:**  
- Visualized image for each predicted sequence 
