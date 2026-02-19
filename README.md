--------------------------------------------------------------------------------
                              SYSTEM IMMUNOLOGY LAB
                                Haifa University
                                 Daniel Fridman 
                                      2026
--------------------------------------------------------------------------------


# Identification of Genetic Motifs in SARS-CoV-2 Specific B-Cell Receptor Repertoires Repertoiry

## 1. OVERVIEW
This repository contains the bioinformatics program used to perform subtition mutation survival analysis
om BCR repertoires from ImmuneDB. 


## 2. PREREQUISITES
Please ensure the following python modules are installed:
- `Pandas`
- `NumPy` 
- `SciPy`
- `tqdm`
- `sqlalchemy`


## 3. USAGE GUIDE
1. Congifgure the `congif.json` file (see section 4 - config) by:
   - rename `congif.json.example` to `congif.json`.
   - Fill the sql connection and database information as illusrated.
   - Once the custom python modules will be loaded the script will initiate
   the required folders and import the `congif.json` information into `config`
   variable.
2. Follow the steps in `notebook_tutorial.ipynb` file.

`Note: Run report will be saved in the reports folder`


## 4. CONFIG.JSON CONFIGURATION
The `config.json` is in a json format and it's purpose is to configure the ImmuneDB MySQL connection
and database on which the process will be performed.
before program usage:
- `sql`: Configure the sql connection information.
  - `adress` - ip adress of the MySQL server.
  - `port` - port of the MySQL server.
  - `username` - username credentials used to connect to the MySQL server.
  - `password` - password credentials used to connect to the MySQL server.
- `database`: Configure the database name and subject id. 
    - `db_name` - Database name.
    - `subject_id` - On which subjects the analysis will be perforemed.
    - `metadata_label` - additional metdata column that will be added to the table.
- `substitution_sruvival`: Configuring the subtition survival analysis steps.
    - `req_tables` - Required tabels for this analysis. 
    - `req_metadata` - Required metadata from the metadata table (column).
    - `new_metadata_labels` - option to rename the metadata column names.
    - `clone_filter` - Filter out clones with less than 10 unique mutations.
    - `fig_remap` - changing the order of the columns.
    - `fig_order` - order of the metadata labels in the correlation figure.
    - `fig_order_change` - boolean value, used to change the labels order in the output figure.


## 5. Pipeline Steps & Components
1. MySQL tables import.
2. Creation of mutations dataframe.
3. filtring of mutations dataframe.
4. creation of correlation talbes and visualization.


## 5. DIRECTORY STRUCTURE
 * data_raw -> Folder which contains the raw tables downloaded from the sql server.
 * data_temp -> Processed tables used for the analysis.
 * reports -> Run reports, consult incase of error.
 * results_figures -> Figure output folder.
 * source -> Source code.


## 6. RESOUCES
- ImmuneDB GitHub - https://github.com/arosenfeld/immunedb

