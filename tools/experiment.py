import requests
import os
import scanpy as sc
import pandas as pd

class Spatial_CxG_Experiment():
    '''
    Class representing a single spatial transcriptomics experiment , 
    as found on the cellxgene datasets website and defined by its dataset id
    '''

    def __init__(
            self,
            dataset_id: str,
            h5ad_folder: str,
    ):
        '''
        Args:
            dataset_id (str): cellxgene dataset_id of the experiment
            h5ad_folder (str): Folder path for storing the downloaded h5ad file
        '''
        super().__init__
        self.dataset_id = dataset_id
        self.dataset_url = self.get_dataset_url()
        self.h5ad_folder = h5ad_folder
        self.h5ad_disk_path = os.path.join(self.h5ad_folder, f'{self.dataset_id}.h5ad')
        self.adata = self.get_adata()
        self.obs_columns = list(self.adata.obs.columns)

    def get_dataset_url(self):
        '''
        Generate download link for the experiments h5ad file
        
        Returns:
            str: Download link for the h5ad file based on the datasets id
        '''
        prefix = 'https://datasets.cellxgene.cziscience.com/'
        suffix = '.h5ad'
        return prefix + self.dataset_id + suffix

    def download_experiment(self):
        '''Downloads the experiments h5ad file from cellxgene'''
        # create folder if necessary
        if not os.path.isdir(self.h5ad_folder):
            os.makedirs(self.h5ad_folder, exist_ok=True)
        #download the respective file
        response = requests.get(self.dataset_url, stream=True)
        if response.status_code == 200:
            with open(self.h5ad_disk_path, 'wb') as file:
                file.write(response.content)
            print('File downloaded')
        else:
            print(f'Download Failed. Status code: {response.status_code}')

    def get_adata(self):
        '''
        Downloads, if not yet downloaded, and reads the experiments h5ad file into memory
        
        Returns:
            AnnData object(anndata._core.anndata.AnnData)
        '''
        if not os.path.isfile(self.h5ad_disk_path):
            self.download_experiment()
        else:
            pass
        return sc.read_h5ad(self.h5ad_disk_path)

    def get_experiment_obs_dict_rev(self, obs_dict: dict[str, list[str]]):
        '''
        For the purpose of standardizing the column names of an AnnData object
        
        Args:
            obs_dict (dict): Dictionary of the shape {'new_column_1': [old_column_1, old column_2, ..., old column_n]}
        
        Returns:
            dict: Dictionary of the shape {'old_column_1': 'new_column_1'}
        '''
        obs_dict_rev = {}
        for column in self.adata.obs.columns:
            for key, value in obs_dict.items():
                if column in value:
                    obs_dict_rev[value] = key
        return obs_dict_rev



class Experiment_Collection():
    '''
    Class representing a collection of spatial transcriptomics experiments,
    as defined by the class Spatial_CxG_Experiment()
    '''

    def __init__(
            self,
            dataset_ids_csv: str,
            h5ad_folder: str,
    ):
        '''
        Args:
            datasets_ids_csv (str): Filepath to a csv, containing all the dataset ids in a column 'dataset_id'
            h5ad_folder (str): Folder path for storing the downloaded h5ad files
        '''
        super().__init__()
        self.dataset_ids = pd.read_csv(dataset_ids_csv)['dataset_id'].tolist()
        self.h5ad_folder = h5ad_folder
        self.experiments = self.get_datasets()
        self.all_obs_columns = self.get_all_obs_columns()
        self.obs_intersection = self.get_obs_intersection()
        self.obs_difference = self.get_obs_difference()
    
    def get_datasets(self):
        '''
        Loops over the dataset ids found in dataset_ids_csv and constructs a class object Spatial_CxG_Experiment() based on that
        Returns:
            list: Contains Spatial_CxG_Experiment() for every dataset_id in dataset_ids
        '''
        experiments = []
        for dataset_id in self.dataset_ids:
            spatial_experiment = Spatial_CxG_Experiment(
                dataset_id=dataset_id,
                h5ad_folder=self.h5ad_folder,
            )
            experiments.append(spatial_experiment)
        return experiments

    def print_experiment_obs_columns(self):
        '''
        Prints:
            pandas.core.indexes.base.Index: adata.obs.columns for every experiment
        '''
        for experiment in self.experiments:
            print(experiment.obs_columns)

    def get_all_obs_columns(self):
        '''
        Returns:
            list: Contains the column names (str) of every experiments adata.obs dataframe 
        '''
        nested_list = [experiment.obs_columns for experiment in self.experiments]
        return sorted(list(set([item for sublist in nested_list for item in sublist])))

    def get_obs_intersection(self):
        '''
        Generate a list of column names (str) shared across every experiments obs_columns
        
        Returns:
            list: Intersection of the set Spatial_CxG_Experiment().obs_columns of all experiments
        '''
        nested_list = []
        for experiment in self.experiments:
            nested_list.append(experiment.obs_columns)

        intersection = []
        for i in self.all_obs_columns:
            if all(i in sublist for sublist in nested_list):
                intersection.append(i)
        return sorted(intersection)
    
    def get_obs_difference(self):
        '''
        Generate a list of column names (str) not shared across every experiments obs_columns
        
        Returns:
            list: Difference of the set Spatial_CxG_Experiment().obs_columns of all experiments
        '''
        return sorted(list(set(self.all_obs_columns) - set(self.obs_intersection)))

    def rename_obs(self, obs_dict: dict[str, list[str]]):
        '''
        Not yet implemented!!!

        Renames every Spatial_CxG_Experiment().adata.obs.columns according to a standardized format defined in a dictionary,
        The key is the new standard name for the column and the values is a list of possible old column names to standardize
        
        Args:
            obs_dict (dict): Dictionary of the shape {'new_column_1': [old_column_1, old column_2, ..., old column_n]},
        '''
        for experiment in self.experiments:
            obs_dict_rev = experiment.get_experiment_obs_dict_rev(obs_dict)
            print(obs_dict_rev)

    def get_obs_unique_values(self, obs_column: str):
        '''
        Generate a list of all the unique values contained in every Spatial_CxG_Experiment().adata.obs[obs_column]

        Args:
            obs_column (str): Name of a column present in one or more Spatial_CxG_Experiment().obs_columns

        Returns:
            list: All the possible values of a column in this Experiment_Collection.datasets instance
        '''
        all_values = []
        for experiment in self.experiments:
            if obs_column in experiment.obs_columns:
                values = list(experiment.adata.obs[obs_column].unique())
                for value in values:
                    all_values.append(value)
            else:
                print(f'Experiment with id {experiment.dataset_id} does not contain obs column {obs_column}')
        return list(set(all_values))