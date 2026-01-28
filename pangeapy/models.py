import os
import pathlib
import requests
import time
import pandas as pd
from functools import reduce


## models information
_anno_models_url ='https://figshare.com/ndownloader/files/58652170'
_meta_models_url ='https://figshare.com/ndownloader/files/58652173'
_anno_models_info_file = 'anno_models.csv'
_meta_models_info_file = 'meta_models.csv'
_anno_model_path = os.getenv('ANNO_MODEL_PATH', default = os.path.join(str(pathlib.Path.home()), '.pangea'))
_meta_model_path = os.getenv('META_MODEL_PATH', default = os.path.join(str(pathlib.Path.home()), '.pangea'))



def _get_url(filename, url):
    max_retries = 5
    wait_time = 3  # seconds
    
    for attempt in range(max_retries):
        try:
            with requests.get(url, stream=True, timeout=30) as r:
                # Case 1: Success (200)
                if r.status_code == 200:
                    with open(filename, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            f.write(chunk)
                    return True
                
                # Case 2: Preparing (202) -> Wait and Retry
                elif r.status_code == 202:
                    print(f"Server is preparing the file (Status 202). Retrying in {wait_time}s... ({attempt+1}/{max_retries})")
                    time.sleep(wait_time)
                    continue  # Jump to next iteration of the loop
                
                # Case 3: Actual Failure (404, 500, etc.)
                else:
                    print(f"Failed to download {url}. Status code: {r.status_code}")
                    return False
                    
        except Exception as e:
            print(f"Error downloading {url}: {e}")
            return False
            
    print("Max retries exceeded.")
    return False

class CellModels():
    def __init__(self, 
                 reset_ref_path = False,
                 model_savedir=None, 
                 automatic_download = True,
                 modelinfo_verbose = False):
        
        if model_savedir is None:
            model_savedir = _anno_model_path

        if reset_ref_path:
            if os.path.isdir(model_savedir+'/anno_models'):
                os.removedirs(model_savedir+'/anno_models')
        os.makedirs(model_savedir, exist_ok=True)
        
        self.model_savedir = model_savedir
        self.models_info_file_path = model_savedir+'/'+_anno_models_info_file
        
        # Importing model metadata
        file_exists = os.path.isfile(self.models_info_file_path) and os.path.getsize(self.models_info_file_path) > 0
        
        # Only download if file is missing, empty, or reset is requested
        if not file_exists or reset_ref_path:
            print(f"Downloading metadata to {self.models_info_file_path}...")
            success = _get_url(self.models_info_file_path, _anno_models_url)
            
            if not success:
                # If download failed and we don't have a backup, we must stop
                if not os.path.isfile(self.models_info_file_path):
                    raise RuntimeError("Could not download anno_models.csv and no local copy exists.")
                else:
                    print("Download failed. Using existing local copy.")

        # Configuring model metadata file
        _modeldf = pd.read_csv(self.models_info_file_path, header=None)
        _modeldf.columns = ['models', 'source']
        _modeldf['Level'] = [i.split("/")[-2] for i in _modeldf['models']]
        _modeldf['cell'] = [i.split("/")[-1].split("_v")[0] for i in _modeldf['models']]
        _modeldf['version'] = ['v'+i.split("/")[-1].split("_v")[1].split("_")[0] for i in _modeldf['models']]
        _modeldf['location'] = [model_savedir+'/'+i for i in _modeldf['models']]

        self.modelinfo = _modeldf.copy()

        # Check-up download status
        _downstat = self.check_downloads()
        if _downstat:
            if modelinfo_verbose:
                print(self.modelinfo.value_counts(['version','Level']))

        # Automatically download anno_models
        if automatic_download:
            if not _downstat:
                self.download_anno_models()

    # Checking-up model download status
    def check_downloads(self):
        print('Checking-up download status of anno_models')
        
        _model_savedir = self.model_savedir
        _modeldf = self.modelinfo

        statusls = [os.path.isfile(_model_savedir+'/'+i) for i in _modeldf['models']]
        return bool(reduce((lambda x, y : x*y), statusls))
        
    # Downloading annotation models from repository
    def download_anno_models(self):
        print('Downloading annotation models...')

        _model_savedir = self.model_savedir
        _modeldf = self.modelinfo

        dummies = [os.makedirs('/'.join([_model_savedir] + i.split("/")[:-1]), exist_ok=True) for i in _modeldf['models']]
        dummies = [_get_url(_model_savedir+'/'+i, j) for i,j in zip(_modeldf['models'],_modeldf['source'])]


class MetaModels():
    def __init__(self, 
                 reset_ref_path = False,
                 model_savedir= None, 
                 automatic_download = True,
                 modelinfo_verbose = False):
        
        if model_savedir is None:
            model_savedir = _meta_model_path
        if reset_ref_path:
            if os.path.isdir(model_savedir+'/meta_models'):
                os.removedirs(model_savedir+'/meta_models')
        os.makedirs(model_savedir, exist_ok=True)
        
        self.model_savedir = model_savedir
        self.models_info_file_path = model_savedir+'/'+_meta_models_info_file
        
        # Importing model metadata
        if os.path.isfile(self.models_info_file_path):
            os.remove(self.models_info_file_path)
        _get_url(self.models_info_file_path, _meta_models_url)

        # Configuring model metadata file
        _modeldf = pd.read_csv(self.models_info_file_path, header=None)
        _modeldf.columns = ['models', 'source']
        _modeldf['type'] = [i.split("/")[-1].split("_v")[0] for i in _modeldf['models']]
        _modeldf['version'] = ['v'+i.split("/")[-1].split("_v")[1].split(".")[0] for i in _modeldf['models']]
        _modeldf['location'] = [model_savedir+'/'+i for i in _modeldf['models']]

        self.modelinfo = _modeldf.copy()

        # Check-up download status
        _downstat = self.check_downloads()
        if _downstat:
            if modelinfo_verbose:
                print(self.modelinfo.value_counts(['version','Level']))

        # Automatically download anno_models
        if automatic_download:
            if not _downstat:
                self.download_meta_models()

    # Checking-up model download status
    def check_downloads(self):
        print('Checking-up download status of meta_models')
        
        _model_savedir = self.model_savedir
        _modeldf = self.modelinfo

        statusls = [os.path.isfile(_model_savedir+'/'+i) for i in _modeldf['models']]
        return bool(reduce((lambda x, y : x*y), statusls))
        
    # Downloading annotation models from repository
    def download_meta_models(self):
        print('Downloading metadata models...')

        _model_savedir = self.model_savedir
        _modeldf = self.modelinfo

        dummies = [os.makedirs('/'.join([_model_savedir] + i.split("/")[:-1]), exist_ok=True) for i in _modeldf['models']]
        dummies = [_get_url(_model_savedir+'/'+i, j) for i,j in zip(_modeldf['models'],_modeldf['source'])]

