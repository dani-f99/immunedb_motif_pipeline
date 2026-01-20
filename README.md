--------------------------------------------------------------------------------
                              SYSTEM IMMUNOLOGY LAB
                                Haifa University
                                 Daniel Fridman 
                                      2026
--------------------------------------------------------------------------------


# Identification of Genetic Motifs in SARS-CoV-2 Specific B-Cell Receptor Repertoires Repertoiry

## 1. OVERVIEW
is repository contains the complete bioinformatics pipeline used to process and analyze 
BCR repertoires from ImmuneDB. It includes the scripts for motif identification and selection 
pressure analysis used in the study: "Identification of Genetic Motifs in SARS-CoV-2 Specific 
B-Cell Receptor Repertoires."


## 2. PREREQUISITES
Please ensure the following python modules are installed:
- `Pandas`
- `NumPy` 
- `SciPy`
- `tqdm`
- `pandarallel`
- `sqlalchemy`
- `pymysql`


## 3. USAGE GUIDE
1. Congifgure the `sql_congif.json` file (see section 4 - config) by:
   - rename `sql_congif.json.example` to `sql_congif.json`.
   - Fill the sql connection and database information as illusrated.
   - Once the custom python modules will be loaded the script will initiate
   the required folders and import the `sql_congif.json` information into `config`
   variable.
2. Run the `main.py` file via cmd

`Note: Run report will be saved in the reports folder`


## 4. CONFIG.JSON CONFIGURATION
The `config.json` is in a json format and it's purpose is to configure the ImmuneDB MySQL connection
and database on which the process will be performed.
before program usage:
- `sql`: Configure the sql connection information.
- `database`: Configure the database name and subject id. 
    - subject_id in the format of "1,2,3,...,n"


## 5. Pipeline Steps & Components
1. MySQL tables import.
2. Creation of 


## 5. DIRECTORY STRUCTURE
# Project tree
 * [tree-md](./tree-md)
 * [dir2](./dir2)
   * [file21.ext](./dir2/file21.ext)
   * [file22.ext](./dir2/file22.ext)
   * [file23.ext](./dir2/file23.ext)
 * [dir1](./dir1)
   * [file11.ext](./dir1/file11.ext)
   * [file12.ext](./dir1/file12.ext)
 * [file_in_root.ext](./file_in_root.ext)
 * [README.md](./README.md)
 * [dir3](./dir3)


## 6. RESOUCES
- See origin authurs pdf tuorial at the `source\bshara` folder (readme.pdf).

