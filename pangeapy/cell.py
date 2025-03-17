import numpy as np
import pandas as pd
import celltypist
import pickle as pkl
import scanpy as sc
from scipy.stats import gmean
from tqdm import tqdm



from .models import CellModels


celltypist.logger.set_level(40)

def read_pickle(filename):
    with open(filename, "rb") as file:
        data = pkl.load(file)
        return data

def process(token):
    return token['text']


class CellAnnotator(CellModels):
    def __init__(self, **kwargs):
        # model_savedir=None, 
        # automatic_download = True,
        # modelinfo_verbose = False

        super().__init__(**kwargs)

    def _predict(self, _adata, _model, majority_voting, celltypist_kwargs):
        """
        Predict celltype labels using celltypist model.
        """

        _pred = celltypist.annotate(_adata, model = _model, majority_voting = majority_voting, **celltypist_kwargs)
        _predres = _pred.predicted_labels.copy()
        _predres['conf_score'] = _pred.probability_matrix.max(axis = 1)
        _predres['cert_score'] = _pred.probability_matrix.max(axis = 1) / _pred.probability_matrix.sum(axis = 1)

        return _predres

    def _annotate_by_sample(self, 
                            adata, 
                            sample_name,
                            target_level, 
                            n_cutoff, 
                            majority_voting, 
                            anno_key,
                            score_key,
                            score_cut,
                            force_l2,
                            celltypist_kwargs):
        modelinfo = self.modelinfo.copy()
        totalidx = adata.obs.index

        _resdic = {}
        for i_level in range(target_level):
            level = 'Level'+str(i_level+1)
            level_prev = 'Level'+str(i_level)
            

            idxs = modelinfo[modelinfo['Level'] == level].index

            results = []
            for idx in idxs:
                _cell = modelinfo.loc[idx, 'cell']
                _model_loc = modelinfo.loc[idx, 'location']
                _model = read_pickle(_model_loc)

                if level == 'Level1':
                    cell_idx = totalidx.copy()
                else:
                    res_prev = _resdic[level_prev].copy()
                    anno_key_prev = level_prev+'|'+anno_key
                    score_key_prev = level_prev+'|'+score_key

                    if anno_key_prev not in res_prev.columns: 
                        continue
                    if force_l2:
                        cell_idx = res_prev[(res_prev[anno_key_prev] == _cell)].index
                    else: 
                        cell_idx = res_prev[(res_prev[anno_key_prev] == _cell) &
                                            (res_prev[score_key_prev] > score_cut)].index
                        if len(cell_idx) < n_cutoff: continue
                
                tdata = adata[cell_idx].copy()
                # print(f'conducting {level} prediction..')
                _predres = self._predict(tdata, _model, majority_voting, celltypist_kwargs)
                results.append(_predres)
            
            if len(results) > 0:
                res = pd.concat(results).copy()
                res.columns = ['|'.join([level, i]) for i in res.columns]
                _resdic[level] = res.copy()
        
        _prediction_results = pd.concat([_resdic[i] for i in _resdic], axis = 1).copy()
        _anno_df = _prediction_results[[i for i in _prediction_results.columns if i.endswith(anno_key)]].copy()
        _score_df = _prediction_results[[i for i in _prediction_results.columns if i.endswith(score_key)]].copy() 
        _prediction_results['PG_annotations'] = ['|'.join([i for i in _anno_df.loc[j] if str(i) != 'nan']) for j in _anno_df.index]
        _prediction_results['PG_combined_score'] = [gmean([i for i in _score_df.loc[j] if str(i) != 'nan']) for j in _score_df.index]
        _prediction_results['Sample'] = sample_name
        
        return _prediction_results
    
    def annotate(self, adata, target_level = 2, n_cutoff = 50, majority_voting = True, 
              anno_key = 'majority_voting', score_key = 'conf_score', force_l2 = False,
              score_cut = 0.5, sample_key = None, celltypist_kwargs = {}):
        
        """
        Predict cell type labels with pre-trained Pangeapy cell type reference models.
        Using celltypist engine (1e4 - logtransformed count required either in adata.X or adata.raw.X)

        By default, hierarchical cell type prediction is conducted (Level1 -> Level2).
        Sublevel annotation is conducted only for the cell types with >= 50 cell counts.

        Parameters
        ----------
        adata : anndata, input single-cell data
            The single-cell data for which we want to get the label predictions.
            Make sure they have 1e4-normalized log1p transformed expression matrix.
        
        target_level: int, 1 or 2 (default)
            Levels of cell type annotation. 
            If it is set to 1, sub-type level (Level2) annotation would not be conducted.

        

        Returns
        -------
        y_pred : ndarray of shape (n_samples,)
            Vector containing the class labels for each sample.
        """
        # 
        resdic = {}
        
        modelinfo = self.modelinfo.copy()
        assert int(modelinfo['Level'].max().split("Level")[1]) >= target_level
        assert anno_key in ['majority_voting', 'predicted_labels']
        n_cutoff = np.max([n_cutoff, 50])

        

        if sample_key is not None:
            
            sample_list = adata.obs[sample_key].unique()
            sample_list = [i for i in sample_list if adata.obs.value_counts(sample_key).loc[i] > n_cutoff]

            def _annotation_process(index):
                _sample = sample_list[index]
                print(f'start cell type annotation: {_sample}')
                adata_sample = adata[adata.obs[sample_key] == _sample].copy()
                
                # Initialize and re-run PCA sample-wise
                if adata_sample.X.min() < 0:
                    assert adata.raw is not None
                    adata_sample = adata_sample.raw.to_adata().copy()
                    adata_sample.raw = adata_sample.copy()

                if 'X_pca' in adata_sample.obsm:
                    del adata.obsm['X_pca']

                _return_df = self._annotate_by_sample( adata = adata_sample, 
                                                            sample_name = _sample,
                                                            target_level=target_level,
                                                            n_cutoff=n_cutoff,
                                                            majority_voting=majority_voting,
                                                            anno_key=anno_key,
                                                            score_key=score_key,
                                                            score_cut = score_cut,
                                                            force_l2=force_l2,
                                                            celltypist_kwargs=celltypist_kwargs)
                return _sample, _return_df
            
            resdic = dict([_annotation_process(i) for i in range(len(sample_list))])
                
        else:
            adata_sample = adata.copy()
            _sample = 'single-run'
            resdic[_sample] = self._annotate_by_sample( adata = adata_sample, 
                                                        sample_name = _sample,
                                                        target_level=target_level,
                                                        n_cutoff=n_cutoff,
                                                        majority_voting=majority_voting,
                                                        anno_key=anno_key,
                                                        score_key=score_key,
                                                        score_cut=score_cut,
                                                        force_l2=force_l2,
                                                        celltypist_kwargs=celltypist_kwargs)

        _prediction_results = pd.concat([resdic[i] for i in resdic], axis = 0, join='outer')             
        
        

        return _prediction_results
                

        

    


