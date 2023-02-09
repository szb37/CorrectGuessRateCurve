"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

Model family analyis naming convention:
analysis_name: {model_family}_cgrc{cgrc_param_set}
subanalysis_name: {model_family}_cgrc{cgrc_param_set}_{model_name}_trial{model_sim_id}


NEW:
name of run: {analysis_name}_{cgrc_param_set}
subanalysis_name: {analysis_name}_{cgrc_param_set}_{model_name}_trial{model_sim_id}
"""

import src.toy_models.model_defs as model_defs
import src.dataframe_classes as df_class
import src.config as config
import src.folders as folders
import src.cgrc.stats as stats
import src.cgrc.core as cgrc
import src.figures as figures
import src.miscs as miscs
from tqdm.contrib.itertools import product as tqdmproduct
from statistics import mean, median
from scipy.stats import bernoulli
import pandas as pd
import itertools
import random
import os


class Controllers():

    @staticmethod
    def run_cgrc_model_family(models, analysis_name, cgrc_param_set, n_patients, n_trials):
        ''' Runs CGRC analysis on a family of models
            Args:
                models(list): list of model model_defs
                analysis_name (str): prefix string of outputfiles
                cgrc_param_set (str): CGRC parameter set, see config.py
                n_patients (int): number of patients in each simulated trail
                n_trials (int): number of trials to simulate for each member of the model family
        '''

        assert isinstance(models, list)
        assert all([model.is_valid() for model in models])
        assert isinstance(analysis_name, str)
        assert isinstance(cgrc_param_set, str)
        assert isinstance(n_patients, int)
        assert isinstance(n_trials, int)

        cgrc_parameters = config.cgrc_parameters[cgrc_param_set]

        trial_data_dir, trial_stats_dir, cgrc_data_dir, cgrc_stats_dir, cgrc_plots_dir = miscs.create_analysis_dirs(
            analysis_name,
            trial_data_subdir=True,
            incl_cgrc_plots_dir=config.save_figs,
        )

        for model, model_sim_id in tqdmproduct(models, range(n_trials), desc='Toy model(s) analysis'):

            model_name = model['name']
            subanalysis_name = f"{analysis_name}_{cgrc_param_set}_{model_name}_trial{model_sim_id}"
            add_columns = {'model_name': model_name, 'model_sim_id': model_sim_id}

            # Generate pseudodata according to model specifications
            ToyModelsDataGenerator.get_toymodel_data(
                output_dir=trial_data_dir,
                output_prefix=subanalysis_name,
                model=model,
                n=n_patients,
                model_sim_id=model_sim_id
            )

            # Get trial stats for psuedodata
            stats.Controllers.get_trial_stats(
                input_dir=trial_data_dir,
                input_fname=subanalysis_name + '__trial_data.csv',
                output_dir=trial_stats_dir,
                output_prefix=subanalysis_name,
                add_columns=add_columns,
            )

            # Get CGRC from psuedodata
            cgrc.CorrectGuessRateCurve.get_cgrc_data(
                input_dir=trial_data_dir,
                input_fname=subanalysis_name + '__trial_data.csv',
                output_dir=cgrc_data_dir,
                output_prefix=subanalysis_name,
                cgrc_parameters=cgrc_parameters,
                add_columns=add_columns,
            )

            # Get psudeodata CGRC stats
            stats.Controllers.get_cgrc_stats(
                input_dir=cgrc_data_dir,
                input_fname=subanalysis_name + '__cgrc_data.csv',
                output_dir=cgrc_stats_dir,
                output_prefix=subanalysis_name,
                add_columns=add_columns,
                do_stratas = False,
            )

            # Make CGRC figues
            if config.save_figs:
                figures.Controllers.plot_VScgr_twinx(
                    input_dir=cgrc_stats_dir,
                    input_fname=subanalysis_name + '__cgrc_model_comps.csv',
                    output_dir=cgrc_plots_dir,
                    output_prefix=subanalysis_name,
                )

        # Get summary table
        #summary_df = ToyModelsAnalyis.get_model_family_summary(
        #    analysis_name=analysis_name,
        #    analysis_name=analysis_name,
        #    cgrc_param_set=cgrc_param_set,
        #)
        #print('\n', summary_df.to_string(index=False))



class ToyModelsDataGenerator():

    @staticmethod
    def get_toymodel_data(output_dir, output_prefix, model, n, model_sim_id=0, enforce_all_strata=True):
        ''' Generate pseudo-experimental data according to toy model definition
            Args:
                output_dir (str): subfolder name where results will be saved
                output_prefix (str): appended to the begining of output files
                model (dict): model specifiction; see model_defs.py
                n (int): number of datapoints in each simulation
                model_sim_id (int, optional): id of current simulations
                enforce_all_strata(bool, optional): if True, then all 4 strata is guaranteed to have
                    at least 1 point in the sample; otherwise errors may occur
        '''

        assert isinstance(output_dir, str)
        assert isinstance(output_prefix, str)
        assert model.is_valid()
        assert isinstance(n, int)
        assert isinstance(model_sim_id, int)

        # Generate conditions
        conditions = ['AC' if condition==1 else 'PL' for condition in bernoulli.rvs(model['p_act'], size=n)]

        # Generate guesses
        guesses=[]
        for condition in conditions:
            if condition=='AC':
                guesses.append(bernoulli.rvs(model['p_sea'], size=1)[0])
            elif condition=='PL':
                guesses.append(bernoulli.rvs(model['p_sep'], size=1)[0])
            else:
                assert False
        guesses = ['AC' if guess==1 else 'PL' for guess in guesses]

        # Ensure all strata represented if needed
        if enforce_all_strata:
            conditions, guesses = ToyModelsDataGenerator.enforce_all_strata(conditions, guesses)

        # Generate scores
        scores=[]
        for idx, condition in enumerate(conditions):
            score = random.gauss(model['nhist'][0], model['nhist'][1])

            if (conditions[idx]=='PL' and guesses[idx]=='PL'):
                pass
            elif (conditions[idx]=='AC' and guesses[idx]=='PL'):
                score += random.gauss(model['dte'][0], model['dte'][1])
            elif (conditions[idx]=='PL' and guesses[idx]=='AC'):
                score += random.gauss(model['aeb'][0], model['aeb'][1])
            elif (conditions[idx]=='AC' and guesses[idx]=='AC'):
                score += random.gauss(model['dte'][0], model['dte'][1]) + random.gauss(model['aeb'][0], model['aeb'][1])
            else:
                assert False

            scores.append(round(score))

        # Format & save output
        toy_model_data_df = ToyModelsDataGenerator.get_TrialDataDf(
            model['name'],
            model_sim_id,
            conditions,
            guesses,
            scores)

        if config.save_csvs:
            toy_model_data_df.to_csv(os.path.join(output_dir, output_prefix+'__trial_data.csv'), index=False)

        return toy_model_data_df

    @staticmethod
    def get_TrialDataDf(model_name, model_sim_id, conditions, guesses, scores):
        ''' Returns TrialDataDf with aux data filled
            Args:
                model_name (str): simulation id
                model_sim_id (int): simulation id
                conditions (list):
                guesses(list):
                scores(list):
        '''

        assert (isinstance(model_name, str) or model_name is None)
        assert isinstance(model_sim_id, int)
        assert isinstance(conditions, list)
        assert isinstance(guesses, list)
        assert isinstance(scores, list)

        n = len(scores)

        toy_model_data_df = df_class.TrialDataDf()

        toy_model_data_df['trial'] = [model_name]*n
        toy_model_data_df['subject_id'] = list(range(n))
        toy_model_data_df['scale'] = ['scale1']*n
        toy_model_data_df['tp'] = ['wk8']*n
        toy_model_data_df['condition'] = conditions
        toy_model_data_df['guess'] = guesses
        toy_model_data_df['baseline'] = [0]*n
        toy_model_data_df['score'] = scores
        toy_model_data_df['delta_score'] = scores

        # Add new column and rearrange order
        toy_model_data_df.add_columns({'model_sim_id':model_sim_id})
        column_order = [
            'trial',
            'model_sim_id',
            'subject_id',
            'scale',
            'tp',
            'condition',
            'guess',
            'baseline',
            'score',
            'delta_score',]
        toy_model_data_df = toy_model_data_df.reindex(columns=column_order)

        toy_model_data_df.__class__ = df_class.TrialDataDf
        toy_model_data_df.set_column_types()

        assert toy_model_data_df.is_valid()
        return toy_model_data_df

    @staticmethod
    def enforce_all_strata(conditions, guesses):
        ''' Overwrites first 4 elements of conditions/guesses to ensure that there is at least
            on row in each treatment/guess strata

            Args:
                conditions(list)
                conditions(guesses)
        '''

        assert isinstance(conditions, list)
        assert len(conditions) >= 4
        assert isinstance(guesses, list)
        assert len(guesses) >= 4

        conditions[0] = 'PL'
        guesses[0] = 'PL'

        conditions[1] = 'AC'
        guesses[1] = 'PL'

        conditions[2] = 'PL'
        guesses[2] = 'AC'

        conditions[3] = 'AC'
        guesses[3] = 'AC'

        return conditions, guesses


class ToyModelsAnalyis():

    @staticmethod
    def get_model_family_summary(analysis_name, cgrc_param_set):
        ''' Construct summary table for model family
            Args:
                analysis_name (str): analysis_name+postfix
                analysis_name (str): name of model family (i.e. list of model definitions, see model_defs)
                cgrc_param_set (int): number defning which CGRC parameter set to use, see config
        '''

        assert isinstance(analysis_name, str)
        assert isinstance(analysis_name, str)
        assert isinstance(cgrc_param_set, str)

        df = df_class.ModelFamilyResultsDf()

        trial_data = Helpers.get_concateneted_df_type(
            target_dir=os.path.join(folders.trial_data_dir, analysis_name),
            df_type='__trial_data')

        # Get n_patients and n_trials assuming each member of the model family has same n_patient and n_trial
        tmp_mid = trial_data.model_name.to_list()[0]
        tmp_tid = trial_data.model_sim_id.to_list()[0]
        n_patients = trial_data.loc[(trial_data.model_name == tmp_mid) & (
            trial_data.model_sim_id == tmp_tid)].shape[0]
        n_trials = len(trial_data.model_sim_id.unique().tolist())

        # Get unadjusted model components
        unadj_model_comps = Helpers.get_concateneted_df_type(
            target_dir=os.path.join(folders.trial_stats_dir, analysis_name),
            df_type='__model_comps')

        # Get adjusted model components
        cgradj_model_comps = Helpers.get_concateneted_df_type(
            target_dir=os.path.join(folders.cgrc_stats_dir, analysis_name),
            df_type='__cgrc_model_comps')

        cgradj_model_comps = cgradj_model_comps.loc[(
            cgradj_model_comps.cgr == 0.5)]
        assert cgradj_model_comps.shape[0] > 0
        cgradj_n_trials = len(
            cgradj_model_comps.model_sim_id.unique().tolist())
        assert n_trials == cgradj_n_trials
        cgradj_model_names = cgradj_model_comps.model_name.unique().tolist()
        unadj_model_names = unadj_model_comps.model_name.unique().tolist()
        assert cgradj_model_names == unadj_model_names

        for model_name in unadj_model_comps.model_name.unique().tolist():

            filtered_unadj_model_comps = unadj_model_comps.loc[
                (unadj_model_comps.model_name == model_name) &
                (unadj_model_comps.model_type == 'without_guess') &
                (unadj_model_comps.component == 'conditionAC')
            ]

            trial_model_data = trial_data.loc[(
                trial_data.model_name == model_name)]
            cgr = round(trial_model_data.loc[(
                trial_model_data.condition == trial_model_data.guess)].shape[0]/(n_trials*n_patients), 3)

            row = {}
            row['model'] = model_name
            row['n_trials'] = n_trials
            row['n_patients'] = n_patients
            row['cgrc_param_set'] = cgrc_param_set
            row['cgr'] = cgr
            row['avg_trt_p'] = round(miscs.get_estimate(
                filtered_unadj_model_comps.p.tolist()), 3)
            row['sig_trt_rate'] = round(
                sum([el <= 0.05 for el in filtered_unadj_model_comps.p.tolist()])/n_trials, 3)
            row['avg_trt_es'] = round(miscs.get_estimate(
                filtered_unadj_model_comps.est.tolist()), 3)

            filtered_cgradj_model_comps = cgradj_model_comps.loc[
                (cgradj_model_comps.model_name == model_name) &
                (cgradj_model_comps.model_type == 'without_guess') &
                (cgradj_model_comps.component == 'conditionAC')
            ]

            # Calculate average p/es aross n_cgr_trials (avg corresponds to p/es of single trial)
            model_sim_ids = filtered_cgradj_model_comps.model_sim_id.unique().tolist()
            trial_ps = []
            trial_efs = []
            for model_sim_id in model_sim_ids:
                tmp = filtered_cgradj_model_comps.loc[filtered_cgradj_model_comps.model_sim_id == model_sim_id]
                assert tmp.shape[0] == config.cgrc_parameters[cgrc_param_set]['n_cgrc_trials']
                trial_ps.append(mean(tmp.p.tolist()))
                trial_efs.append(mean(tmp.est.tolist()))

            row['cgradj_avg_trt_p'] = round(miscs.get_estimate(trial_ps), 3)
            row['cgradj_sig_trt_rate'] = round(
                sum([p <= 0.05 for p in trial_ps])/n_trials, 3)
            row['cgradj_avg_trt_es'] = round(miscs.get_estimate(trial_efs), 3)

            df = df.append(row, ignore_index=True)

        df.__class__ = df_class.ModelFamilyResultsDf
        df.set_column_types()
        df.to_csv(os.path.join(folders.summary_dir, analysis_name +
                  '_{}__summary_table.csv'.format(config.estimator)), index=False)

        return df

    @staticmethod
    def get_concateneted_df_type(target_dir, df_type):
        ''' Concatenets and returns all dataframes of a certain kind from target_dir
            Args:
                target_dir(str): path to folder
                df_type(str): dataframe type; see assertions for valid types
        '''

        assert df_type in [
            '__trial_data',
            '__model_comps',
            '__strata_contrast',
            '__cgrc_model_comps',
            '__cgrc_strata_contrast']

        # Get all trial stats
        if df_type == '__model_comps':
            master_df = df_class.ModelComponentsDf()
        elif df_type == '__strata_contrast':
            master_df = df_class.StrataContrastDf()
        elif df_type == '__cgrc_model_comps':
            master_df = df_class.ModelComponentsDf()
            master_df.add_columns({'cgr': None, 'cgr_sim_id': None})
        elif df_type == '__cgrc_strata_contrast':
            master_df = df_class.StrataContrastDf()
            master_df.add_columns({'cgr': None, 'cgr_sim_id': None})
        elif df_type == '__trial_data':
            master_df = df_class.TrialDataDf()
        else:
            assert False

        master_df.add_columns({'model_name': None, 'model_sim_id': None})

        for fpath in [fpath for fpath in os.listdir(target_dir) if (df_type in fpath)]:
            df = pd.read_csv(os.path.join(target_dir, fpath))
            master_df = pd.concat([master_df, df], sort=False)

        if df_type == '__model_comps':
            master_df.__class__ = df_class.ModelComponentsDf
        elif df_type == '__strata_contrast':
            master_df.__class__ = df_class.StrataContrastDf
        elif df_type == '__cgrc_model_comps':
            master_df.__class__ = df_class.ModelComponentsDf
        elif df_type == '__cgrc_strata_contrast':
            master_df.__class__ = df_class.StrataContrastDf
        elif df_type == '__trial_data':
            master_df.__class__ = df_class.TrialDataDf
        else:
            assert False

        master_df.set_column_types()
        return master_df
