import pandas as pd
import numpy as np
import anndata as ad

from .models import MetaModels

import pickle as pkl
def read_pickle(filename):
    with open(filename, "rb") as file:
        data = pkl.load(file)
        return data


class MetaAnnotator(MetaModels):
    def __init__(self, **kwargs):
        # model_savedir=None, 
        # automatic_download = True,
        # modelinfo_verbose = False
        super().__init__()

    def annotate(self, 
                 pred_res,
                 anno_key = 'majority_voting',
                 sample_key = 'Sample',
                 sample_n_cut = 50,
                 organ_prediction = True,
                 organ_prediction_cutoff = .5,
                 phenotype_class = 'auto-detect'):
        

        _sample_count = pred_res[sample_key].value_counts()
        _samples = _sample_count[_sample_count >= sample_n_cut].index

        if len(_samples) < 1:
            print(f'insufficient cell counts... all of the sample counts are below {sample_n_cut}')
            return None
        
        meta_results = MetaResults()
        meta_results.conduct_organ_prediction = organ_prediction
        resdic = {}
        for _sample in _samples:
            _pred_res_by_sample = pred_res[pred_res[sample_key] == _sample].copy()
            _metares = self._annotate_by_sample(pred_res=_pred_res_by_sample,
                                                anno_key=anno_key,
                                                organ_prediction=organ_prediction,
                                                organ_prediction_cutoff=organ_prediction_cutoff,
                                                phenotype_class=phenotype_class)
            if _metares == None:
                print(f'insufficient cell counts for {_sample}')
            resdic[_sample] = _metares
        
        meta_results.resultsdic = resdic

        return meta_results
    
    def _annotate_by_sample(self, 
                        pred_res,
                        anno_key,
                        organ_prediction,
                        organ_prediction_cutoff,
                        phenotype_class):

        _cell_compositions = self._get_cell_compositions(pred_res=pred_res, 
                                                         anno_key=anno_key)
        
        _metares = MetaAnnotation(modelinfo=self.modelinfo[['version','type']])
        _metares.cellcomposition = _cell_compositions

        if _cell_compositions is None:
            _metares.log = 'unable to calculate cell compositions'
            return _metares

        # 1) Organ prediction
        if organ_prediction:
            _organ_predres = self._predict_organ_identity(_cell_compositions)
            _organ_pred, _organ_pred_prob = _organ_predres.idxmax(), _organ_predres.max()

            _metares.organ_prediction = _organ_pred, _organ_pred_prob
            _metares.organ_prediction_df = _organ_predres


        # 2) Meta phenotyping 
        if phenotype_class == 'auto-detect':
            assert organ_prediction == True
            if (_organ_pred_prob >= organ_prediction_cutoff) & (_organ_pred == 'Blood'):
                phenotype_class = 'Blood'
            else:
                phenotype_class = 'Tissue'
        
        if phenotype_class == 'Blood':
            _phenotype_predres, _input_feature_ref = self._predict_blood_phenotypes(_cell_compositions)

        if phenotype_class == 'Tissue':
            _phenotype_predres, _input_feature_ref = self._predict_tissue_phenotypes(_cell_compositions)
        _pheno_pred, _pheno_pred_prob = _phenotype_predres.idxmax(), _phenotype_predres.max()
                
        _metares.phenotype_prediction_df = _phenotype_predres
        _metares.phenotype_prediction = _pheno_pred, _pheno_pred_prob        
        _metares.ref_input_feature = _input_feature_ref
        _metares.no_input_feature = len(set([i.split("|")[0] for i in _input_feature_ref if i in _cell_compositions]))
        _metares.no_input_feature_ref = len(set([i.split("|")[0] for i in _input_feature_ref]))

        _metares.log = 'done'
        return _metares


    def _get_cell_compositions(self,
                               pred_res, 
                               anno_key,
                               l1_cutoff= 0.5,
                               l2_cutoff= 0.5,
                               n_cut = 50,
                               N_cut = 500):
    
        l1key = 'Level1|'+anno_key
        l2key = 'Level2|'+anno_key
        l1score = 'Level1|conf_score'
        l2score = 'Level2|conf_score'

        if l1key not in pred_res.columns:
            return None

        pred_res = pred_res[(pred_res[l1score] > l1_cutoff) & 
                            (pred_res[l2score] > l2_cutoff)].copy()
        if pred_res.shape[0] < N_cut:
            return None
        
        # Level1 proportions
        p1 = pred_res.value_counts(l1key, normalize=True)

        # Level2 proportions
        if l2key in pred_res.columns:
            p2ls = []
            pred_res_l2 = pred_res[~pred_res[l2key].isna()].copy()
            l2_celltypes = pred_res_l2[l1key].unique()
            
            ## calculating Level2 composition in each Level1 cell types
            for cell in l2_celltypes:
                _predres_level2 = pred_res_l2[pred_res_l2[l1key] == cell].copy()
                if _predres_level2.shape[0] < n_cut:
                    continue

                p2 = _predres_level2.value_counts(l2key, normalize=True)
                p2.index = [cell+'|'+i for i in p2.index]
                p2ls.append(p2)  
            prop = pd.concat([p1]+p2ls)
        else:
            prop = p1.copy()
        
        return prop

    def _predict_organ_identity(self, cell_composition):
        _modelinfo = self.modelinfo.copy()
        _model_loc = _modelinfo[_modelinfo['type'] == 'Organ_predictor']['location'].values[0]
        _meta_model = read_pickle(_model_loc)
        _proportion = cell_composition

        X = np.array([_proportion.loc[i] if i in _proportion.index else 0 for i in _meta_model.feature_names_in_]).reshape(1,-1)
        _organ_predres = pd.Series(_meta_model.predict_proba(X).flatten(), index =  _meta_model.classes_)

        return _organ_predres
        
    def _predict_blood_phenotypes(self, cell_composition):
        _modelinfo = self.modelinfo.copy()
        _model_loc = _modelinfo[_modelinfo['type'] == 'Blood_predictor']['location'].values[0]
        _meta_model = read_pickle(_model_loc)
        _proportion = cell_composition

        X = np.array([_proportion.loc[i] if i in _proportion.index else 0 for i in _meta_model.feature_names_in_]).reshape(1,-1)
        _phenotype_predres = pd.Series(_meta_model.predict_proba(X).flatten(), index =  _meta_model.classes_)
        

        return _phenotype_predres, _meta_model.feature_names_in_
        
    def _predict_tissue_phenotypes(self, cell_composition):
        _modelinfo = self.modelinfo.copy()
        _model_loc = _modelinfo[_modelinfo['type'] == 'Tissue_predictor']['location'].values[0]
        _meta_model = read_pickle(_model_loc)
        _proportion = cell_composition

        X = np.array([_proportion.loc[i] if i in _proportion.index else 0 for i in _meta_model.feature_names_in_]).reshape(1,-1)
        _phenotype_predres = pd.Series(_meta_model.predict_proba(X).flatten(), index =  _meta_model.classes_)
        

        return _phenotype_predres, _meta_model.feature_names_in_
    

class MetaAnnotation():
    def __init__(self, modelinfo):
        self.modelinfo = modelinfo

class MetaResults():
    def integrate(self):
        resdic = self.resultsdic

        _indexlist = [i for i in resdic if resdic[i].log == 'done']
        integrated = pd.DataFrame(index = _indexlist)

    

        if self.conduct_organ_prediction:
            integrated['Organ_pred'] = [resdic[i].organ_prediction[0] for i in _indexlist]
            integrated['Organ_prob'] = [resdic[i].organ_prediction[1] for i in _indexlist] 
    
        integrated['Pheno_pred'] = [resdic[i].phenotype_prediction[0] 
                                    if (resdic[i].no_input_feature / resdic[i].no_input_feature_ref) > .5
                                    else np.nan for i in _indexlist] 
        integrated['Pheno_prob'] = [resdic[i].phenotype_prediction[1] 
                                    if (resdic[i].no_input_feature / resdic[i].no_input_feature_ref) > .5
                                    else np.nan for i in _indexlist]

        return integrated
    
    