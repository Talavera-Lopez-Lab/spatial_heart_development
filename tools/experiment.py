import requests
import os
import scanpy as sc

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