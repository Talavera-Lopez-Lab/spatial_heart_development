import requests
import os
import scanpy as sc
import pandas as pd

class Spatial_CxG_Experiment():

    def __init__(
            self,
            dataset_id,
            h5ad_folder,
    ):
        super().__init__
        self.dataset_id = dataset_id
        self.dataset_url = self.get_dataset_url()
        self.h5ad_folder = h5ad_folder
        self.h5ad_disk_path = os.path.join(self.h5ad_folder, f'{self.dataset_id}.h5ad')
        self.adata = self.get_adata()
        self.obs_columns = list(self.adata.obs.columns)

    def get_dataset_url(self):
        prefix = 'https://datasets.cellxgene.cziscience.com/'
        suffix = '.h5ad'
        return prefix + self.dataset_id + suffix

    def download_experiment(self):
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
        if not os.path.isfile(self.h5ad_disk_path):
            self.download_experiment()
        else:
            pass
        return sc.read_h5ad(self.h5ad_disk_path)

    def get_experiment_obs_dict_rev(self, obs_dict):
        obs_dict_rev = {}
        for column in self.adata.obs.columns:
            for key, value in obs_dict.items():
                if column in value:
                    obs_dict_rev[value] = key
        return obs_dict_rev



class Experiment_Collection():

    def __init__(
            self,
            dataset_ids_csv,
            h5ad_folder,
    ):
        super().__init__()
        self.dataset_ids = pd.read_csv(dataset_ids_csv)['dataset_id'].tolist()
        self.h5ad_folder = h5ad_folder
        self.experiments = self.get_datasets()
        self.all_obs_columns = self.get_all_obs_columns()
        self.obs_intersection = self.get_obs_intersection()
        self.obs_difference = self.get_obs_difference()
    
    def get_datasets(self):
        experiments = []
        for dataset_id in self.dataset_ids:
            spatial_experiment = Spatial_CxG_Experiment(
                dataset_id=dataset_id,
                h5ad_folder=self.h5ad_folder,
            )
            experiments.append(spatial_experiment)
        return experiments

    def print_experiment_obs_columns(self):
        for experiment in self.experiments:
            print(experiment.obs_columns)

    def get_all_obs_columns(self):
        nested_list = [experiment.obs_columns for experiment in self.experiments]
        return sorted(list(set([item for sublist in nested_list for item in sublist])))

    def get_obs_intersection(self):
        nested_list = []
        for experiment in self.experiments:
            nested_list.append(experiment.obs_columns)

        intersection = []
        for i in self.all_obs_columns:
            if all(i in sublist for sublist in nested_list):
                intersection.append(i)
        return sorted(intersection)
    
    def get_obs_difference(self):
        return sorted(list(set(self.all_obs_columns) - set(self.obs_intersection)))

    def rename_obs(self, obs_dict):
        for experiment in self.experiments:
            obs_dict_rev = experiment.get_experiment_obs_dict_rev(obs_dict)
            print(obs_dict_rev)

    def get_obs_unique_values(self, obs_column):
        all_values = []
        for experiment in self.experiments:
            if obs_column in experiment.obs_columns:
                values = list(experiment.adata.obs[obs_column].unique())
                for value in values:
                    all_values.append(value)
            else:
                print(f'Experiment with id {experiment.dataset_id} does not contain obs column {obs_column}')
        return list(set(all_values))