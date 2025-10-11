import numpy as np
import pandas as pd
import celltypist
import pickle as pkl
import scanpy as sc
import anndata as ad

from scipy import sparse
from scipy.stats import gmean, zscore

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
                            compute_uncertainty,
                            uncertainty_kwargs,
                            celltypist_kwargs):
        modelinfo = self.modelinfo.copy()
        totalidx = adata.obs.index

        if sample_name is not None:
            print(f'{sample_name}: processing...')

        _resdic = {}
        for i_level in range(target_level):
            level = 'Level'+str(i_level+1)
            level_prev = 'Level'+str(i_level)
            
            print(f'conducting {level} annotation...')

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
                    
                    if (len(cell_idx) < n_cutoff) and not force_l2: continue
                    if len(cell_idx) < 1: continue 
                
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
        _prediction_results = _prediction_results.loc[adata.obs.index]

        if compute_uncertainty:
            
            print('computing uncertainty scores...')

            _uncert_kwargs = dict( 
                input_score_key = 'Level1|conf_score',
                input_graph_key = "distances",      
                n_neighbors = 10,
                ascending = True,            
                negate = True,
                uncert_cutoff = 1,
                cert_cutoff = 0.8,
                self_vote = .5,            
                output_score_key = "score_uncert",
                output_pred_key = "pred_uncert",
                output_consen_key = "consensus_uncert")

            if uncertainty_kwargs:
                _uncert_kwargs.update(uncertainty_kwargs)


            try:
                _uncert_results = self._compute_uncertainty(adata, input_pred_res=_prediction_results, **_uncert_kwargs)
            except:
                _uncert_results = pd.DataFrame()
            
            _prediction_results = _prediction_results.join(_uncert_results, how = 'left')

        return _prediction_results
    
    def _compute_uncertainty(
        adata: ad.AnnData,
        input_pred_res: pd.DataFrame,
        input_score_key: str,
        input_graph_key: str = "distances",      # use "connectivities" if you prefer
        n_neighbors: int = 10,
        ascending: bool = True,            # distances: True keeps smallest; connectivities: set False to keep largest
        negate: bool = True,
        uncert_cutoff: float = 1,
        cert_cutoff: float = 0.8,
        self_vote: float = .5,            # decisive self-vote weight
        output_score_key: str = "score_uncert",
        output_pred_key: str = "pred_uncert",
        output_consen_key: str = "consensus_uncert",
    ):
        """
        1) Z-score (optionally negate) the score_key.
        2) Use neighbors graph to compute neighbor-averaged score (row-stochastic).
        3) Threshold to raw boolean, then majority-vote over n_neighbors neighbors (+ self_vote).
        """
        # 1) z-score the score column
        z = zscore(np.asarray(input_pred_res[input_score_key], dtype=float))
        if negate: z = -z

        # 2) load sparse adjacency (binary)
        if input_graph_key not in adata.obsp:
            print(f"⚠️ Graph '{graph_key}' not found; computing neighbors...")
            sc.pp.neighbors(adata)

        A = adata.obsp[input_graph_key].tocsr().astype(np.float32)
        A.data[:] = 1.0  # binarize

        # 3) Keep only top-n_neighbors per row
        if input_graph_key == 'distances': ascending = True
        if input_graph_key == 'connectivities': ascending = False

        for i in range(A.shape[0]):
            s, e = A.indptr[i], A.indptr[i+1]
            if e - s > n_neighbors:
                if ascending: 
                    top_idx = np.argsort(A.data[s:e])[:n_neighbors]
                else:
                    top_idx = np.argsort(A.data[s:e])[-n_neighbors:]

                mask = np.zeros(e - s, dtype=bool)
                mask[top_idx] = True
                A.data[s:e] = np.where(mask, A.data[s:e], 0)
 
        A.eliminate_zeros()

        # 3) row-normalize adjacency
        deg = np.asarray(A.sum(axis=1)).ravel()
        deg[deg == 0] = 1.0  # avoid division by zero
        P = sparse.diags(1.0 / deg) @ A

        # 4) neighbor average and classification
        score_uncert = P @ z
        pred_raw = (score_uncert > uncert_cutoff) & (input_pred_res[input_score_key] < cert_cutoff)

        # 5) One-hot encoding
        cats = pd.Categorical(pred_raw)
        y = cats.codes
        onehot = sparse.csr_matrix(
            (np.ones(len(y), dtype=np.float32), (np.arange(len(y)), y)),
            shape=(len(y), len(cats.categories))
        )

        # 6) Aggregate neighbor votes (decisive self-vote)
        if self_vote > 0: 
            V = A @ onehot + onehot*self_vote
        else:
            V = A @ onehot
        V = V.tocsr()

        ## 7) majority voting
        mv_codes = np.asarray(V.argmax(axis=1)).ravel()

        top_raw = V.max(axis=1)
        top = np.array(top_raw.todense() if hasattr(top_raw, "todense") else top_raw).ravel().astype(float)

        total_raw = V.sum(axis=1)
        total = np.array(total_raw.todense() if hasattr(total_raw, "todense") else total_raw).ravel().astype(float)

        conf = np.where(total > 0, top / total, np.nan)

        # 8) labeling

        out_df = pd.DataFrame(
                    {
                        output_score_key: np.asarray(score_uncert).ravel(),
                        output_pred_key: pd.Categorical.from_codes(mv_codes, categories=cats.categories),
                        output_consen_key: conf,
                    },
                    index = adata.obs.index,
                )

        return out_df

    def annotate(self, adata, target_level = 2, n_cutoff = 50, force_l2 = False, score_cut = 0.5,  
                 majority_voting = True, anno_key = 'majority_voting', score_key = 'conf_score', sample_key = None, 
                 compute_uncertainty = False, uncertainty_kwargs = {}, celltypist_kwargs = {}):
        
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
                                                            target_level = target_level,
                                                            n_cutoff = n_cutoff,
                                                            majority_voting = majority_voting,
                                                            anno_key = anno_key,
                                                            score_key = score_key,
                                                            score_cut = score_cut,
                                                            force_l2 = force_l2,
                                                            compute_uncertainty=compute_uncertainty,
                                                            uncertainty_kwargs = uncertainty_kwargs,
                                                            celltypist_kwargs=celltypist_kwargs)
                return _sample, _return_df
            
            resdic = dict([_annotation_process(i) for i in range(len(sample_list))])
                
        else:
            adata_sample = adata.copy()
            _sample = 'single-run'
            resdic[_sample] = self._annotate_by_sample( adata = adata_sample, 
                                                        sample_name = _sample,
                                                        target_level = target_level,
                                                        n_cutoff = n_cutoff,
                                                        majority_voting = majority_voting,
                                                        anno_key = anno_key,
                                                        score_key = score_key,
                                                        score_cut = score_cut,
                                                        force_l2 = force_l2,
                                                        compute_uncertainty=compute_uncertainty,
                                                        uncertainty_kwargs = uncertainty_kwargs,
                                                        celltypist_kwargs = celltypist_kwargs)

        _prediction_results = pd.concat([resdic[i] for i in resdic], axis = 0, join='outer')             
        
        

        return _prediction_results
                

        

    


