# Import required packages
from sqlalchemy import create_engine
from datetime import datetime
import importlib.metadata
import pandas as pd
import unittest
import json
import sys
import os


##############################################################
# Custom function that cheeck if python packages are installed
def check_packages(package_names):
    status = {}

    n_total = len(package_names)
    n_installed = 0
    package_df = pd.DataFrame(index=package_names, columns=["installed","not_installed"], data=0)

    for pkg in package_names:
        try:
            # metadata.version returns the version string if installed
            dist_version = importlib.metadata.version(pkg)
            status[pkg] = f"Installed (v{dist_version})"
            package_df.loc[pkg, "installed"] = 1
            n_installed += 1

        except importlib.metadata.PackageNotFoundError:
            status[pkg] = "Not Installed !!!"
            package_df.loc[pkg, "not_installed"] = 1
    
    print(f"> {n_installed}/{n_total} packages are installed:")

    for i in status:
        print(f"{i} package is {status[i]}")

    return package_df

####################################################################
# DNA codon table - used for the translation from nt to aa sequences
protein = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
            "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
            "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
            "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
            "TCT" : "s", "CCT" : "P", "ACT" : "T", "GCT" : "A",
            "TCC" : "s", "CCC" : "P", "ACC" : "T", "GCC" : "A",
            "TCA" : "s", "CCA" : "P", "ACA" : "T", "GCA" : "A",
            "TCG" : "s", "CCG" : "P", "ACG" : "T", "GCG" : "A",
            "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
            "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
            "TAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
            "TAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
            "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
            "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
            "TGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
            "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G", 
            "---" : "-"
            }


################################################################################$$##########
# Reading information from json file. Used to extract the parameters from the `config.json`.
def read_json(path:str = "config.json") -> dict:
    """
    path : str -> path of the json file
    """

    with open(path) as config:
        config_f = json.load(config)

    return config_f


#####################################################
# Creating folder according to the and program scheme
def create_folders(req_folders : list = ["temp_data", "tms_input", "reports"]):
    """
    req_folders : str -> required folders path, if subfolder exsits input '\\' between folders.
    """

    for folder in req_folders:
        if os.path.exists(folder) is False:
            os.mkdir(folder)
            print(f"> folder `{folder}` was created.")

        else:
            print(f"> folder `{folder}` exists, continuing.")


######################################################################################################
# A custom function that connects to MySQL server, execute query and returns the results as dataframe.
class mysql_qry(): 
    def __init__(self,
                 username:str = read_json()["sql"]["username"],
                 password:str= read_json()["sql"]["password"],
                 adress:str= read_json()["sql"]["adress"],
                 port:str= read_json()["sql"]["port"],
                 db_name:str = read_json()["database"]["db_name"]):
        """
        username:str -> Username credentials for the MySQL server.
        password:str -> Password credentials for the MySQL server.
        adress:str -> IP adress of the MySQL server.
        port:str -> Port of the MySQL server.
        qry:str -> SQL query to be executed.
        """
    
        # Setting up MySQL connenction
        connection_mysql = f"mysql+pymysql://{username}:{password}@{adress}:{port}/{db_name}"
        self.engine = create_engine(connection_mysql)
        print(f"> Established connecntion to the {db_name} database.")

    # Executing the query
    def run_qry(self,
                qry):
        qry_df = pd.read_sql(qry, self.engine)
        #print(f"> Query executed successfully: \n{qry}")

        # Returing the table as pd.dataframe
        return qry_df

    
    # Closing the connenction
    def close_conn(self):
        self.engine.dispose()
        print("> MySQL connenction terminated.")
    


######################################################
# A small helper to send output to both CMD and a file
class OutputTee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


# ##########################
# Running Pipeline with test
def run_pipeline(test_pipeline,
                 pipeline_name:str = ""
                 ):
    """
    test_pipeline -> the pipeline uninitited unittest pipeline we want to run
    pipeline_name : str -> pipeline name in string format.
    """
    current_time = datetime.now().strftime("%Y-%m-%d-%H-%M")
    db = read_json()["database"]["db_name"]
    reports_path = os.path.join("reports", db)
    report_name = f"{pipeline_name}_[{current_time}]_report_.txt"
    create_folders([reports_path])
    

    # f is the text file -> sys.stdout is the CMD consol
    with open(os.path.join(reports_path, report_name), "w", encoding="utf-8") as f:
        # sys.stdout is the CMD console
        # f is your text file
        dual_stream = OutputTee(sys.stdout, f)
        
        runner = unittest.TextTestRunner(
            stream=dual_stream, 
            verbosity=2, 
            descriptions=True
        )

        # Initialize the runner
        runner = unittest.TextTestRunner(
                stream=dual_stream, 
                verbosity=2, 
                descriptions=True
                )
        
        # unitest  loader object
        loader = unittest.TestLoader()

        # Load tests from the specific class
        suite = loader.loadTestsFromTestCase(test_pipeline) 
        
        # Run with high verbosity for detail
        result = runner.run(suite)
        
        # Custom detailed summary
        print("\n--- PIPELINE EXECUTION SUMMARY ---")
        if result.wasSuccessful():
            print("Final Status: SUCCESS V")
        else:
            print(f"Final Status: FAILED X ({len(result.failures) + len(result.errors)} issues found)")


#####################################################################################
# Function that create the required folders and defines paths for the data processing 
def init_paths(analysis_name : str,
               database_name : str) -> list:
    
    """
    Function that create the required folders and defines paths for the data processing.
    the folders final path will be: `main_folder/analysis_name/database_name`, for example.
    the output will be list of 4 paths:"data_raw/...", "data_temp/...", "results_figures/...", "results_tables/..."
    analysis name : str -> name of the anlysis performed, for example "lpa_analysis".
    database_name : name of the database, imported from the `config.json` file.
    """
    
    # Defining and creating the main folders
    paths = ["data_temp", "results_figures", "results_tables"]
    create_folders(paths + ["data_raw"])
    create_folders([os.path.join("data_raw", database_name)])

    # Defining and creating the sub-folders according to the input arguments.
    for i in [analysis_name, database_name]:
        paths = [os.path.join(j, i) for j in paths]
        create_folders(paths)
        
    return paths


########################################################################
# Function that create report file and update it's values based on input
def reports(report_name : str,
            step : str):
    """
    report_name : str -> name of the report file in the format of `analysis-database_report` (.txt will be added later).
    step_name : str / list -> the line which be added to the report, if several lines the input should be a list.
    """
    # Creating reports folder
    create_folders(["reports"])
    
    # Initilizing paths and names
    report_path = os.path.join("reports", report_name+".txt")
    report_bool = os.path.exists(report_path)
    name_splited = report_name.split("-")
    steps_dict = {True:step, False:[step]}
    steps = steps_dict[isinstance(step, (list, tuple))]
    
    lines_count = 1
    if report_bool:
        with open(report_path, "r") as report:
            lines_count = len(report.readlines())


    # Creating new report file and adding header, if file doesnt exists
    # If report file already exists -> add the required line only
    with open(report_path, "a") as report:
        if report_bool is False:
            try:
                report.write(f"> Report File of {name_splited[0].replace("_"," ")} Analysis ({name_splited[1]} database).\n")
            except:
                report.write(f"> Report File of {report_name}\n")

        for i in steps:
            report.write(f"{lines_count} {i}\n")
            lines_count += 1