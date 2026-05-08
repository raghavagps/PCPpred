# PCPpred

## About PCPpred 
PCPpred is a webserver (https://webs.iiitd.edu.in/raghava/pcppred/) and standalone package (https://drive.usercontent.google.com/download?id=1kfFQRpqqluPE0jIWn2zdJg5rbpwvPVOu&export=download&authuser=0) designed for cyclic peptide researchers. It enables:
- Conversion of peptide sequences from MAP (Modifications and Annotations in Protein) format to SMILES or HELM notation.
- Permeability prediction for cyclic peptides using SMILES as input across multiple assays: PAMPA, Caco-2, RRCK, and MDCK.
---
# Zenodo
https://doi.org/10.5281/zenodo.19910197
# dataset.zip
Contain all the dataset

# Prerequisites

- Python: 3.12.6  
- Java: JRE 6+ (Required for PaDEL-Descriptor)

# Install dependencies using:
```bash
pip install -r requirements.txt
```

`requirements.txt`:
```
joblib==1.4.2
mordred==1.2.0
numpy==1.26.4
padelpy==0.1.16
pandas==2.2.3
rdkit==2025.3.1
scikit-learn==1.6.0
scipy==1.13.1
seaborn==0.13.2
tqdm==4.66.5
transformers==4.44.2
xgboost==2.1.3
torch
lightgbm==4.5.0
```
## All the data used and processed for the predicton of permeability for different assay type can be downloaded from here: Data:  https://drive.usercontent.google.com/download?id=1nhW0Qc9IGv_hFwz7OaCRjVej-JdhqeWT&export=download&authuser=0
## Chemical language models fine-tuned to predict permeability of cyclic peptides for all assay types can be downloaded from here: CLM models: https://drive.usercontent.google.com/download?id=1NR9z0x9gFNOxnUBkG8ANN62NqZqkEUpY&export=download&authuser=0
## Stacked ensemble models can be downloaded from here: Models: https://drive.usercontent.google.com/download?id=1C2GXA4nitT8wqt_Rl7WA8AJvK1HuBWVK&export=download&authuser=0

## User can directly use SMILES representation (one entry in each line) as an input in the text box or in an input file or they can choose to design their custom cyclic peptide in MAP (modification and annotation of proteins) format which can be converted into SMILES representation or can directly be used as input to predict Permeability of peptides.

# Functionalities iniclude: 

# MAP to SMILES Converter
## Overview
This Python script converts peptide sequences in MAP (Modifications and annotations in protein) format to SMILES (Simplified Molecular Input Line Entry System) notation. It supports both single MAP sequence inputs and batch processing from a file, utilizing the RDKit library for chemical structure manipulation.
Prerequisites

## Usage
The script can be run from the command line with two modes: single sequence or file input.
Single Sequence Mode
Convert a single MAP sequence to SMILES and print the result to the console.

## CLI:
```bash
python map_to_smiles.py -s "Peptide_001\n{nnr:ABU}{nnr:0OZ}{nnr:9XD}V{nnr:9XD}AA{d}{nnr:9XD}{nnr:9XD}{nnr:0Q3}{nnr:MBM}{cyc:N-C}"
```

## File Input Mode
Convert multiple MAP sequences from an input file and write the corresponding SMILES to an output file. You can specify a custom output file path (including directory) using the -o option. If not specified, the output will be written to results/smiles_output_<input_filename>.

## CLI:
```bash
python map_to_smiles.py -f input_map_sequences.txt -o /path/to/output/smiles.txt
```
or, with default output path:
```bash
python map_to_smiles.py -f input_map_sequences.txt
```
Input File Format: For each peptide there will be two line, first header line, eg: ">peptide001" and second line should contain one MAP sequence.
Output: SMILES strings are written to the specified output file (or default results/smiles_output_<input_filename>) in the same order as the input sequences. The script will create the output directory if it does not exist.

# Notes
If cyclization or SMILES generation is not possible, an error message is returned for the specific sequence.

## Acknowledgments
Original code by Charles Xu and others (2021-2024)

# MAP to HELM Converter
## Overview
This Python script converts peptide sequences in MAP (Molecular Assembly Pattern) format to HELM (Hierarchical Editing Language for Macromolecules) notation. It supports both single MAP sequence inputs with a peptide ID and batch processing from a file, utilizing the RDKit library for chemical structure manipulation and pandas for data handling.

## Usage
The script can be run from the command line with two modes: single sequence or file input.
Single Sequence Mode
Convert a single MAP sequence to HELM notation and print the result to the console. A peptide ID must be provided.

## CLI:
```bash
python map_to_helm.py -s "Peptide_001\n{nnr:ABU}{nnr:0OZ}{nnr:9XD}V{nnr:9XD}AA{d}{nnr:9XD}{nnr:9XD}{nnr:0Q3}{nnr:MBM}{cyc:N-C}" -i "001"
```
## File Input Mode
Convert multiple MAP sequences from an input file and write the corresponding HELM sequences to an output file. You can specify a custom output file path (including directory) using the -o option. If not specified, the output will be written to results/helm*output*<input_filename>.

## CLI:
```bash
python map_to_helm.py -f input_map_sequences.txt -o /path/to/output/helm_sequences.txt
```
or, with default output path:
```bash
python map_to_helm.py -f input_map_sequences.txt
```
Input File Format: For each peptide there will be two line, first header line, eg: ">peptide001" and second line in the input file should contain a MAP sequence and a peptide ID, separated by a comma (e.g., MAP_sequence,peptide_id).
Output: HELM sequences are written to the specified output file (or default results/helm_output_<input_filename>) in the same order as the input sequences. The script will create the output directory if it does not exist.

The input file must have each line formatted as MAP_sequence,peptide_id. Lines with incorrect formatting will result in an error message in the output file.

# Note: To get the SMILES for custom peptides use MAP to SMILES converter to create SMILES and then use them to predict.

```****************************************************************************************************************************```

# SMILES PAMPA Permeability Prediction
## Overview
This Python script predicts PAMPA-based permeability for cyclic peptides  based on SMILES strings. It uses a stacked ensemble machine learning architecture combining molecular descriptors, fingerprints, embeddings, and atomic features. 


## Usage

This script is command-line driven and operates in batch mode only — predicting permeability for all SMILES in the given input file.

### Command Line Format
```bash
python predict_permeability_pampa.py --input <path_to_input_smiles_file> --model <model_name> [--output <path_to_output_file>]
```

### Arguments

- `--input` (str, required):  
  Path to the input `.txt` file containing SMILES strings (one per line).

- `--model` (str, required):  
  Meta-model to use for prediction. Choose from:  
  `lgb`, `decision_tree`, `random_forest`, `gradient_boosting`, `adaboost`, `xgb`, `extra_trees`, `linear`, `knn`, `svr`, `mlp`.

- `--output` (str, optional):  
  Path to the output `.csv` file. If not provided, the output will be saved to `results/output_pampa.csv`.

### Example Usages

**Basic Prediction:**
```bash
python predict_permeability_pampa.py --input data/smiles_input.txt --model random_forest
```

**Prediction with Custom Output File:**
```bash
python predict_permeability_pampa.py --input data/smiles_input.txt --model mlp --output results/mlp_predictions.csv
```

## Input File Format

- Plain text file (`.txt`)
- Each line should contain one valid SMILES string.

**Example (smiles_input.txt):**
```
CC(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)O
CCN(CC)CCCC(C)NC(=O)c1ccc(Cl)cc1
CC(C)C[C@@H](NC(=O)[C@H](Cc1ccccc1)N)C(=O)O
```

## Output Format

- Output is a `.csv` file with two columns:
  - `SMILES`: Original input SMILES
  - `Permeability`: Predicted permeability 



# SMILES Caco-2 Permeability Prediction

## Overview
This Python script predicts **Caco-2 cell permeability** for cyclic peptides based on SMILES strings. It uses a stacked ensemble architecture combining descriptors, fingerprints, embeddings, and atomic-level features. 

## Usage

The script is command-line driven and performs **batch prediction** from a `.txt` file containing SMILES.

### Command Line Format
```bash
python predict_permeability_caco2.py --input <path_to_input_smiles_file> [--output <path_to_output_file>]
```

### Arguments
- `--input` (str, required):  
  Path to the input `.txt` file containing SMILES strings (one per line).

- `--output` (str, optional):  
  Path to the output `.csv` file. If not provided, the output will be saved to `results/output_caco2.csv`.

### Example Usages

**Basic Prediction:**
```bash
python predict_permeability_caco2.py --input data/smiles_input.txt
```

**Prediction with Custom Output File:**
```bash
python predict_permeability_caco2.py --input data/smiles_input.txt --output results/caco2_predictions.csv
```

## Input File Format

- Plain text file (`.txt`)
- Each line should contain one valid SMILES string.

**Example (smiles_input.txt):**
```
CC(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)O
CCN(CC)CCCC(C)NC(=O)c1ccc(Cl)cc1
CC(C)C[C@@H](NC(=O)[C@H](Cc1ccccc1)N)C(=O)O
```

## Output Format

- Output is a `.csv` file with two columns:
  - `SMILES`: Original input SMILES
  - `Permeability`: Predicted permeability 


# SMILES RRCK Permeability Prediction

## Overview

This Python script predicts **RRCK cell permeability** for cyclic peptides represented as SMILES strings. It uses a stacked ensemble of machine learning models based on four molecular representation types: descriptors, fingerprints, embeddings, and atomic-level features.

## Usage

The script is executed via the command line and accepts SMILES input from a `.txt` file. The results are saved to a `.csv` file.

### Command Line Format

```bash
python predict_permeability_rrck.py --input <path_to_input_smiles_file> [--output <path_to_output_file>]
```

### Arguments

- `--input` (str, required):  
  Path to the input `.txt` file containing SMILES strings, one per line.

- `--output` (str, optional):  
  Path to the output `.csv` file for saving predictions. If not provided, output is saved to `results/output_rrck.csv`.

### Example Usages

**Basic usage:**
```bash
python predict_permeability_rrck.py --input data/smiles_input.txt
```

**With custom output path:**
```bash
python predict_permeability_rrck.py --input data/smiles_input.txt --output results/rrck_predictions.csv
```

## Input Format

- Plain text file (`.txt`)
- One SMILES string per line

**Example:**
```
CC(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)O
CCN(CC)CCCC(C)NC(=O)c1ccc(Cl)cc1
CC(C)C[C@@H](NC(=O)[C@H](Cc1ccccc1)N)C(=O)O
```

## Output Format

The output is a `.csv` file with the following columns:

- `SMILES`: Input SMILES string
- `Permeability`: Predicted RRCK permeability 


# SMILES MDCK Permeability Prediction
## Overview

This Python script predicts **MDCK cell permeability** of cyclic peptides from their SMILES strings. It uses the Klekota-Roth fingerprint representation and an AdaBoost ensemble regressor model (5-fold averaging) to generate the permeability prediction.

## Usage

The script takes a text file of SMILES strings and outputs a `.csv` file with predicted permeability values. Fingerprints are calculated using PaDEL-Descriptor.

### Command Line Format

```bash
python predict_permeability_mdck.py --input <path_to_input_smiles_file> [--output <path_to_output_file>]
```

### Arguments

- `--input` (str, required):  
  Path to the input `.txt` file containing SMILES strings, one per line.

- `--output` (str, optional):  
  Path to the output `.csv` file for saving predictions. If not provided, output is saved to `results/output_mdck.csv`.

### Example Usages

**Basic usage:**
```bash
python predict_permeability_mdck.py --input data/smiles_input.txt
```

**With custom output path:**
```bash
python predict_permeability_mdck.py --input data/smiles_input.txt --output results/mdck_predictions.csv
```

## Input Format

- A plain text (`.txt`) file
- One valid SMILES string per line

**Example:**
```
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
CN1CCCC1C2=CC=CC=C2
CC(C)C1=CC=C(C=C1)O
```

## Output Format

The output `.csv` file will contain two columns:

- `SMILES`: The original SMILES input
- `Permeability`: The predicted MDCK permeability value

---


