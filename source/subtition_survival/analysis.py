# Custom Functions Imports
from source.helpers import init_paths, reports, read_json, mysql_qry


# Python Modules Import
from distutils.util import strtobool
import matplotlib.pyplot as plt
from datetime import datetime
from os.path import exists
from scipy import stats
from tqdm import tqdm
import seaborn as sns
import pandas as pd
import numpy as np
import regex as re
import json
import ast
import os
import re



#####################################
# Subtitution survival analysis class
class sub_mut():
    def __init__(self):  
        # Getting the config information and database name
        self.config_info = read_json("config.json")
        self.db_name = self.config_info["database"]["db_name"]
        self.analysis_type = "substitution_sruvival"


        # Initilizing the analysis paths and folders
        # folder output order: "data_raw", "data_temp", "results_figures", "results_tables"
        self.paths = init_paths(self.analysis_type, self.db_name)
        self.current_time = datetime.now().strftime('[%d-%m-%Y %H;%M]')
        self.report_name = f"{self.analysis_type}-{self.db_name}_report_{self.current_time}"

        # Importing the required raw tables from the MySQL server
        req_tables = self.config_info[self.analysis_type]["req_tables"].split(",")
        sql_engine = mysql_qry()
        imported_tabeles = []

        # Importing each required table from the MySQL server
        for i in tqdm(req_tables, desc="Importing MySQL tables", unit="tables"):
            i_table_path = os.path.join(self.paths[0], f"{i}.csv")
            i_exists = os.path.exists(i_table_path)

            # Cheecking if the raw table already exsits, if so importing
            if i_exists is False:
                temp_qry = f"""
                            SELECT * FROM {self.db_name}.{i};
                            """
                
                # Incase of invalid input error
                try:
                    temp_df = sql_engine.run_qry(temp_qry)
                    temp_df.to_csv(i_table_path)
                    imported_tabeles.append(i)

                except:
                    print(f"> Invalid database or table name. (db={self.db_name}, table={i})")
                    continue
        
        # I table already exists -> moving on
        else:
            print(f"> table `{i}.csv` alread exists.")

        # Verifing the existance of all of the rquired tables
        
        if (np.sort(["clones", "clone_stats", "sample_metadata"]) == np.sort([i.split(".")[0] for i in os.listdir(self.paths[0])])).all():
            # Finishing report
            print(f"> Finished importing raw tabels ({os.listdir(self.paths[0])}).")
            sql_engine.close_conn()

            # Adding to the report -> Successfully imported (or verified) MySQL tables.
            reports(self.report_name,"Successfully imported (or verified) MySQL tables.")
        
        else:
            raise Exception("Verify config file for correct tables values.")


    # Creation of mutation table
    def mutation_table(self):
        """
        * Creating orginized metadata dataframe with the information provided by the config.json file.
        * Saving the metadata_df into folder.
        * If the dataframe already exists, load it without processing.
        """
        path_processed_dir = self.paths[1]
        path_metdadata_df = os.path.join(path_processed_dir, "sample_metadata_df.csv")

        if exists(path_metdadata_df):
            print("> sample_metadata_df.csv already created, continuing...")
            metadata_df = pd.read_csv(path_metdadata_df, index_col=0)
                                                                                        
        else:
            print("> Creating sample_metadata_df.csv.")
        
            metadata_keys_og = self.config_info["substitution_sruvival"]["req_metadata"].split(",")
            metadata_keys_new =  self.config_info["substitution_sruvival"]["new_metadata_labels"].split(",")
            meta_dict = dict(zip(metadata_keys_og, metadata_keys_new))

            metadata_df = pd.read_csv(os.path.join(self.paths[0], "sample_metadata.csv"), index_col=0)
            metadata_og = metadata_df[metadata_df["key"].isin(metadata_keys_og)]
            metadata_og = metadata_og.replace({"key":meta_dict})

            sample_ids = np.sort(metadata_og["sample_id"].unique())
            metadata_df = pd.DataFrame({"sample_id":sample_ids})
            metadata_df[metadata_keys_new] = np.nan

            for i in sample_ids:
                temp_sid = i
                for j in metadata_keys_new:
                    cond_sid = (metadata_og["sample_id"] == i)
                    cond_key = (metadata_og["key"] == j)
                    metadata_df.loc[metadata_df["sample_id"]==i,j] = metadata_og.loc[(metadata_og["sample_id"]==i)&(metadata_og["key"]==j),"value"].values
            metadata_df.to_csv(path_metdadata_df)
            print("> Done.")

        # Creation of filtred metadata table
        """
        * Creating custom function to pull metadata from metadata_df
        * sample_id validation
        """
        def assign_metadata(sample_id, meta_df):
            meta_list = meta_df.columns[1:]
            meta_sample = meta_df.loc[meta_df["sample_id"]==sample_id, meta_list].values.flatten()
            return meta_sample
        
        clone_stats = pd.read_csv(os.path.join(self.paths[0], "clone_stats.csv"), index_col=0)
        metalist_sids = np.sort(metadata_df.sample_id.unique())
        clones_sids = np.sort(clone_stats.dropna().sample_id.unique()).astype("int")

        values_missing = np.setdiff1d(clones_sids, metalist_sids)
        values_common = np.intersect1d(metalist_sids, clones_sids)

        if len(values_missing) > 0:
            print("> missing sample_id from metadata file:",values_missing)
            clone_stats = clone_stats[clone_stats["sample_id"].isin(values_common)]
            raise TypeError("verify metadata sample_id values") 
        
        # Merging clones and clones status > adding the relevent metadata to the dataframe
        """
        * loading clones_merged if exists, if not creating and saving
        * custom function that extract mutations infromation from the "mutation" json in each row
        * Adding the germline infromation to the clone_stats df
        * Dropping null sample_id rows (cannot assign metadata for those rows)
        * converting "sample_id" values to int instead of floats
        * assiging the metadata into the merged table (applying assign_metadata)
        * renaming id_x to id after merging (left had "id" colum while right had "id"=="clone_id")
        * reseting the index
        """
        
        path_clones_merged = os.path.join(self.paths[1], "clones_merged.csv")

        if exists(path_clones_merged):
            clones_merged = pd.read_csv(path_clones_merged)
            print("> clones_merged dataframe exists, loading and continuing....")

        else: 
            clones = pd.read_csv(os.path.join(self.paths[0], "clones.csv"))
            clones_merged = clone_stats.merge(right=clones[["id","germline"]],
                                                how="left",
                                                left_on="clone_id",
                                                right_on="id")    
            
            clones_merged = clones_merged[clones_merged["sample_id"].notnull()]        
            clones_merged[list(metadata_df.columns)[1:]] = list(clones_merged["sample_id"].apply(assign_metadata, args=(metadata_df,)))
            clones_merged.rename({"id_x":"id"},axis="columns",inplace=True)
            clones_merged.reset_index(drop=True, inplace=True)
            clones_merged.to_csv(path_clones_merged)
            print("> clones_merged dataframe created and saved, continuing....")

        # Creating the mutation dataframe
        """
        * Creating df with the relevent mutations infromation for each clone (mut_df)
        * Cleaning the mut_df and adding relevent data
        * Saving the mut_df
        """

        path_mut_df = os.path.join(self.paths[-1], "mut_df.csv")

        if exists(path_mut_df):
            mut_df = pd.read_csv(path_mut_df,index_col=0)
            print("> mut_df dataframe exists, loading and continuing....")

        else: 
            def mut_regall(string):
                pattern = r"'pos': (?P<position>\d+), 'from_nt': '(?P<from_nt>[\w]+)', 'from_aa': '(?P<from_aa>[\w\*]+)', 'to_nt': '(?P<to_nt>['\w\*]+)', 'to_aas': \[(?P<to_aas>['\w,\s\*]+)], 'unique': (?P<unique>\d+), 'total': (?P<total>\d+), 'intermediate_aa': '(?P<intermediate_aa>[\w\d\*])'"
            
                tjson = json.loads(string)
                
                if "ALL" in str(tjson["regions"].keys()):
                    all_value = str(tjson["regions"]["ALL"])
                    find = re.findall(pattern,all_value)
                    return find
                
                else:
                    else_value = str(tjson["regions"])
                    return else_value
                    
            clones_merged["regions_all"] = clones_merged["mutations"].apply(mut_regall)
            clones_raval = clones_merged.copy()
            ra_val = []
            
            for i in range(0,len(clones_raval)):
                length = len(clones_raval.loc[i,"regions_all"]) # length of the list, number of mutations is the colum
                value = clones_raval.loc[i,"regions_all"] # the value mutations themselves list of lists/ list / np.nan
                id_value = clones_raval.loc[i,"id"] # id value of the row
                clone_id = clones_raval.loc[i,"clone_id"] # clone_id value of the row
                subject_id = clones_raval.loc[i,"subject_id"]# subject_id value of the row
                sample_id = clones_raval.loc[i,"sample_id"] # sample_id value of the row
                funct = clones_raval.loc[i,"functional"] # functional value of the clone
                total_cnt = clones_raval.loc[i,"total_cnt"] # target of the antibody
                unique_cnt = clones_raval.loc[i,"unique_cnt"] # unique sequence is clone
                germline = clones_raval.loc[i,"germline"] #germline sequence
                top_seq = clones_raval.loc[i,"top_copy_seq_sequence"] #top copy of sequence
                metadata = clones_raval.loc[i,metadata_df.columns[1:]].values.flatten().tolist() #metadata list value
                ins_val = [id_value, clone_id, subject_id, sample_id, funct, total_cnt,unique_cnt, germline, top_seq] + metadata
                
                # if single row of mutation
                if length == 1:
                    to_aas = value[0][4].replace(" ","").replace("''","").split(",")
                    
                    if (len(to_aas) == 1):
                        temp_list = list(value[0])
                        ra_val.append(ins_val + temp_list) 
                        
                    else:
                        for i in range(0,len(to_aas)):
                            temp_list = list(value[0])
                            temp_list[4] = to_aas[i]
                            ra_val.append(ins_val + temp_list)
                
                # if multiple rows of mutations
                if length > 1:
                    for j in range(0,length):
                            sub_value = list(value[j]) #each row
                            temp_list = sub_value
                            
                            # making sure that the length of the list is correct, in some rows there is missing values
                            if len(sub_value) == 8:
                                to_aas = sub_value[4].replace(" ","").replace("'","").split(",")
                                
                                if len(to_aas) == 1:
                                    ra_val.append(ins_val + temp_list)
                                elif len(to_aas) > 1:
                                    for aa in set(to_aas): # set() removes duplicate values
                                        temp_list[4] = aa
                                        ra_val.append(ins_val + temp_list)
                                                
                elif length == 0:
                    ra_val.append(ins_val + np.full(shape=len(value), fill_value=np.nan).tolist())
            
            mut_df_cols = ["id", "clone_id", "subject_id", "sample_id", "functional", "total_cnt","unique_cnt", "germline", "top_seq"]
            mut_info_cols = ["pos","from_nt","from_aa","to_nt","to_aas","unique","total","intermidiate_aa"]
            
            mut_df = pd.DataFrame(data=ra_val, columns = mut_df_cols + metadata_df.columns[1:].tolist() + mut_info_cols)
            mut_df["to_aas"] = mut_df["to_aas"].str.replace("'","") #cleaning to_aas string
            mut_df.replace({"to_aas":{"None":np.nan}}, inplace=True) #turining none values to np.nan
            mut_df.dropna(axis=0,subset=["pos","to_aas"], ignore_index=True, inplace=True) #dropping null rows of "pos" and "to_aas"

            # custom function to round numbers upward
            def round_up(number):
                num_dec = number
                num_round = round(number)
                
                if num_round < num_dec:
                    value = num_round + 1
                else:
                    value = num_round
                return value
            
            mut_df.insert(6,"pos_aa",np.nan) #inserting amino acid position column
            mut_df.insert(6,"pos_nt",np.nan) #inserting nucleotide position column
            mut_df.loc[:,"pos_nt"] = mut_df.loc[:,"pos"].apply(int)+1 #filling the pos_nt column
            mut_df.loc[:,"pos_aa"] = ((mut_df.loc[:,"pos_nt"])/3).apply(round_up) #fillint the pos_aa column 
            mut_df.astype({"pos_nt":"int","pos_aa":"int"})
            mut_df.drop(axis=1,columns="pos",inplace=True) #dropping the og column (it was -1 in position...)
            mut_df["syn"] = (mut_df["from_aa"] == mut_df["to_aas"]).apply(int) #creating syn column

            mut_df.to_csv(path_mut_df)
            print("> mut_df dataframe created and saved, continuing....")
            
        reports(self.report_name, "Successfully created mutation dataframe (mut_df.csv)")


    def filter_mutdf(self,
                        aa_range : tuple = (1,104),
                        unique_seq_filt : int = 1,
                        save_filtdf : bool = False
                        ):
        
        ###
        """
        aa_range : tuple / list -> range of the amino acids position which will be included in the analysis.
        unique_seq_filt : int -> Threshold (>) of required unique sequences per clone.
        save_filtdf : bool -> To save the filtred dataframe into the temp_data folder.
        * only functional clones
        * Include only non-synonymous mutations. 
        """
        try:
            mut_df = pd.read_csv(os.path.join(self.paths[-1], "mut_df.csv"))

        except:
            raise Exception(f"> No `mut_df.csv` table found in {self.paths[-1]} folder, please run the sub_mut.mutation_table() method.")

        filt_pos_aa = (mut_df["pos_aa"].between(aa_range[0], aa_range[1], inclusive='both')) # from aa positions 1->104
        filt_unique_cnt = (mut_df["unique_cnt"] > unique_seq_filt) # only clones with more than 1 unique sequence
        filt_functional = (mut_df["functional"] == 1) # only functional clones
        filt_syn = (mut_df["syn"] == 0) # only non-syn mutations (substitutions)

        self.filt_mut_df = mut_df[filt_pos_aa & filt_functional &  filt_syn & filt_unique_cnt]

        if save_filtdf:
            self.filt_mut_df.to_csv(os.path.join(self.paths[1], "filt_mut_df.csv"))

        reports(self.report_name, "Successfully filtred mutation dataframe (filt_mut_df.csv)")

        return self.filt_mut_df

    # Creation of correlation dataframe and correlationsmap
    def correlate(self,
                    remap_list : list = None,
                    aa_range : tuple = (1,104)):
        """
        remap_list = if we want to remap the metadata columns values, we can input list for dicts with values to change (same order as metdata).
        aa_range : tuple / list -> range of the amino acids position which will be included in the analysis.
        """

        ###
        """
        * Custom function that creates matrix of all possible conditions across the metadata
        """
        def metadata_cond(dic_input):
            mt_dic = dic_input
            mt_keys = mt_dic.keys()
            mt_keys_len = np.array([len(mt_dic[i]) for i in mt_keys])

            n_columns = len(mt_keys)
            n_rows = np.prod(mt_keys_len)
            meta_df = pd.DataFrame(np.zeros((n_rows,n_columns)), columns = mt_keys)

            for i,v,l in zip(range(0,len(mt_keys)), mt_keys, mt_keys_len):
                array_length = n_rows
                unique_val = v
                
                if i == 0:          
                    reps_numbers = np.prod(mt_keys_len[i+1:])

                    temp_array = []
                    for val_i in mt_dic[unique_val]:
                        temp_array+=[val_i for i in range(1,reps_numbers+1)]

                    meta_df[v] = temp_array

                else:
                    reps_numbers = np.prod(mt_keys_len[i+1:])
                    temp_array = []
                    
                    for val_i in mt_dic[unique_val]:
                        temp_array += [val_i for k in range(1,reps_numbers+1)]

                    final_array = temp_array.copy()
                    for num in range(1,int(array_length/len(temp_array))):
                        final_array += temp_array.copy()

                    meta_df[v] = final_array   
                    
            return meta_df
        
        ###
        """
        * Creating the metadata dic
        * Using the custom function metadata_cond to create the conditions matrix
        """
        catagories_dict = {}
        change_order = bool(strtobool(self.config_info["substitution_sruvival"]["fig_order_change"]))
        if change_order:
            val_list = self.config_info["substitution_sruvival"]["fig_order"].split(",")
            
            for i in val_list:
                temp_val_list = list(self.filt_mut_df[i].unique())
                temp_val_list.sort()
                catagories_dict[i] = temp_val_list


        else:
            unique_subj = list(self.filt_mut_df["subject_id"].unique())
            
            catagories_dict["subject_id"] = unique_subj
            val_list = list(metadata_df.columns[1:])
            val_list.insert(0,"subject_id")
            
            for i in val_list:
                catagories_dict[i] = np.sort(self.filt_mut_df[i].unique())   
        
        ###
        """
        * option to remap the metadata values if remap == True
        * creation of remaping dictionary (if remap == True)
        * Creating metadata matrix
        """
        remap = bool(strtobool(self.config_info["substitution_sruvival"]["fig_remap"]))
        if remap:
            remap_name = self.db_name  + ";" + ".".join(val_list)
            remap_path = "remapping\\"+remap_name+"-remap"+".txt"
            catag_path = "remapping\\"+remap_name+"-catag"+".txt"


            custom_list = remap_list
            remap_vals = {old_val:new_val for old_dic,new_dic in zip(catagories_dict.values(),custom_list) for old_val,new_val in zip(old_dic,new_dic)}
            catagories_dict = {i:j for i,j in zip(list(catagories_dict.keys()),custom_list)}
                    
            self.filt_mut_df.loc[:,val_list] = self.filt_mut_df.loc[:,val_list].replace(remap_vals)
            catagories_dict = {i:j for i,j in zip(list(catagories_dict.keys()),[i for i in catagories_dict.values()])}
            
            sep = "."
        
        else:
            sep = "|"

        cond_matrix = metadata_cond(catagories_dict).astype("str")
        cond_matrix.to_csv(os.path.join(self.paths[1], "metadata_matrix.csv"))
        
        ###
        """
        * Defining custom function to get frequences of surviving non-syn clones.
        """
        def get_output(df, cname, prange, method):
            df_output = pd.DataFrame({"pos_aa":pos_range})
            total_unique = len(df["clone_id"].unique())
            temp_num = df.groupby("pos_aa").agg({"clone_id":"nunique"}).reset_index()
            
            if method == "freq":
                temp_freq = pd.DataFrame({"pos_aa":temp_num["pos_aa"],cname:temp_num["clone_id"]/total_unique})
            elif method == "count":
                temp_freq = pd.DataFrame({"pos_aa":temp_num["pos_aa"],cname:temp_num["clone_id"]}) 
        
            df_output = df_output.merge(temp_freq, on="pos_aa", how="left")
            df_output.fillna(0,inplace=True)               
            return df_output
        
        ###
        """
        """
        # https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
        range_all = range(aa_range[0],aa_range[1]+1) #FR1, CDR1, FR2, CDR2, FR3

        df_input = self.filt_mut_df.copy()
        pos_range = range_all
        output_summerylis = []
        df_output = pd.DataFrame({"pos_aa":pos_range})

        for i in range(0,len(cond_matrix)):
            df_temp = df_input.copy()
            temp_vals = cond_matrix.iloc[i,:].values

            for name, val in zip(val_list,temp_vals): 
                df_temp = df_temp[df_temp[name] == val]

            col_name = sep.join(temp_vals)
            output_summerylis.append([col_name,len(df_temp["clone_id"].unique())])
            
            df_getoutput = get_output(df_temp,col_name,pos_range,method="freq")
            df_output = df_output.merge(df_getoutput, on="pos_aa", how="left")

        df_output.to_csv(os.path.join(self.paths[-1], "substitution_fractions.csv"))

        ###
        """
        * Creates dataframe with the filtred dataset values and their unique clones.
        * Giving report on small datasets.
        """
        # creates dataframe with the filtred dataset values and their unique clones.
        output_summery = pd.DataFrame(output_summerylis, columns=["dataset","nunique"])
        small_datasets = list(output_summery.loc[output_summery["nunique"]<int(self.config_info["substitution_sruvival"]["clone_filter"]),"dataset"].values)

        if len(small_datasets) > 0:
            for i in small_datasets:
                print("Dropped "+i+" (unique clones:"+ str(output_summery.loc[output_summery["dataset"]==i,"nunique"].values[0])+")")
                df_output.drop(columns=[i],inplace=True)
                output_summery = output_summery[output_summery["dataset"]!= i]

        output_summery.sort_values("nunique").to_csv(os.path.join(self.paths[-1], "summery_nclones.csv"))

        ###
        """
        * Importing statistics and re modules
        * Create empty dataframe with the names of the output_df datasets
        * Getting spearman_r values between each of the fequencies dataset and putting them into the spearman output df.
        """

        names = list(df_output.columns)[1:] #order by name of the filtered loop
        spearmanr_df = pd.DataFrame(data=np.nan,index=names,columns=names)
        spearman_list = []
        
        for i in names:
            for j in names:
                re_pattern = r"statistic=np\.float64\(([\d\.\-e]+)\), pvalue=np.float64\(([\d\.\-e]+)\)\)" #updated
                temp_string = str(stats.spearmanr(df_output[i],df_output[j]))
                temp_spearmanr = re.findall(re_pattern,temp_string)[0][0]
                temp_pvalue = re.findall(re_pattern,temp_string)[0][1]

                spearmanr_df.loc[i,j] = float(temp_spearmanr)
                spearman_list.append([i+ " X " +j, temp_spearmanr, temp_pvalue])

                # if some value isnt statisticly sifnificant print the value.
                if float(temp_pvalue) > 0.05:
                    print(i+ " x " +j +" Not Significant")

        # spearman_r matrix output
        spearman_summery = pd.DataFrame(spearman_list,columns=["dataset", "spearman_r", "p_value"])
        spearman_summery.to_csv(os.path.join(self.paths[-2], f"spearman_heatmap_{self.current_time}.csv"))
        reports(self.report_name, "Successfully saved correlation dataframe.")

        ###
        """
        * importing visualision module
        * creating heatmap graph
        """

        if len(spearmanr_df)/4 < 5:
            size_graph = 5
        else:
            size_graph = len(spearmanr_df)/4

        fig, ax = plt.subplots(1,1, figsize=(size_graph ,size_graph ))

        ax = sns.heatmap(spearmanr_df,
                        cbar_kws={'label': 'Spearman r'},
                        xticklabels = True, 
                        yticklabels = True,
                        cmap = sns.color_palette("coolwarm", as_cmap=True, n_colors=1),
                        vmin = round(spearmanr_df.min().min(),1)-0.05, 
                        vmax = 1)
        ax.invert_yaxis()
        ax.axis('scaled')

        labels_xy = labels=spearmanr_df.columns

        hide_ticklabels = False
        if hide_ticklabels:
            ax.tick_params(
                axis='both',       # changes apply to both x and y axes
                which='both',       # both major and minor ticks are affected
                bottom=False,       # ticks along the bottom edge are off
                top=False,          # ticks along the top edge are off
                left=False,         # ticks along the left edge are off
                right=False,        # ticks along the right edge are off
                labelbottom=False,  # labels along the bottom edge are off
                labelleft=False)    # labels along the left edge are off
            
        save_path = os.path.join(self.paths[-2], f"correlation_heatmap_{self.current_time}.png")
        plt.savefig(save_path, bbox_inches='tight')
        reports(self.report_name, "Successfully created correlation heatmap figure.")
