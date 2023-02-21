"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

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
import pingouin
import random
import os


class Controllers():

    @staticmethod
    def run_toymodels_cgrc(models, analysis_name, cgrc_param_set, n_patients, n_trials):
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

        conditions = [condition for condition in bernoulli.rvs(model['p_act'], size=n)]
        guesses = [guess for guess in bernoulli.rvs(model['p_unb'], size=n)]

        for idx in range(n):
            if guesses[idx]==1:
                guesses[idx] = conditions[idx]
            elif guesses[idx]==0 and conditions[idx]==0:
                guesses[idx] = 1
            elif guesses[idx]==0 and conditions[idx]==1:
                guesses[idx] = 0
            else:
                assert False

        # convert numerics to labels
        conditions = ['AC' if condition==1 else 'PL' for condition in conditions]
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
        toymodel_data_df = ToyModelsDataGenerator.get_TrialDataDf_w_metadata(
            model['name'],
            model_sim_id,
            conditions,
            guesses,
            scores)

        toymodel_data_df.add_columns({'model_name':model['name']})

        if config.save_csvs:
            toymodel_data_df.to_csv(os.path.join(output_dir, output_prefix+'__trial_data.csv'), index=False)

        return toymodel_data_df

    @staticmethod
    def get_TrialDataDf_w_metadata(model_name, model_sim_id, conditions, guesses, scores):
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

        toymodel_data_df = df_class.TrialDataDf()

        toymodel_data_df['trial'] = [model_name]*n
        toymodel_data_df['subject_id'] = list(range(n))
        toymodel_data_df['scale'] = ['scale1']*n
        toymodel_data_df['tp'] = ['wk8']*n
        toymodel_data_df['condition'] = conditions
        toymodel_data_df['guess'] = guesses
        toymodel_data_df['baseline'] = [0]*n
        toymodel_data_df['score'] = scores
        toymodel_data_df['delta_score'] = scores

        # Add new column and rearrange order
        toymodel_data_df.add_columns({'model_sim_id':model_sim_id})
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
        toymodel_data_df = toymodel_data_df.reindex(columns=column_order)

        toymodel_data_df.__class__ = df_class.TrialDataDf
        toymodel_data_df.set_column_types()

        assert toymodel_data_df.is_valid()
        return toymodel_data_df

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

        """
        conditions[4] = 'PL'
        guesses[4] = 'PL'
        conditions[5] = 'AC'
        guesses[5] = 'PL'
        conditions[6] = 'PL'
        guesses[6] = 'AC'
        conditions[7] = 'AC'
        guesses[7] = 'AC'
        """

        return conditions, guesses


class ToyModelsAnalyis():

    @staticmethod
    def get_toymodels_summary(analysis_name):
        ''' Construct summary table for model family
          Args:
              analysis_name (str):
        '''

        assert isinstance(analysis_name, str)

        df = df_class.ToymodelsAnalysisDf()

        # Get trial _data
        trial_data = ToyModelsAnalyis.get_concateneted_df_type(
          target_dir=os.path.join(folders.trial_data_dir, analysis_name),
          df_type='__trial_data')
        cgr_trial_data = ToyModelsAnalyis.get_concateneted_df_type(
          target_dir=os.path.join(folders.cgrc_data_dir, analysis_name),
          df_type='__cgrc_data')

        # Get unadjusted and adjusted model components
        model_comps = ToyModelsAnalyis.get_concateneted_df_type(
          target_dir=os.path.join(folders.trial_stats_dir, analysis_name),
          df_type='__model_comps')
        cgr_model_comps = ToyModelsAnalyis.get_concateneted_df_type(
          target_dir=os.path.join(folders.cgrc_stats_dir, analysis_name),
          df_type='__cgrc_model_comps')
        cgr_model_comps = cgr_model_comps.loc[(cgr_model_comps.cgr == 0.5)]

        # Sanity checks
        assert (
            cgr_model_comps.shape[0] > 0)
        assert (
            len(trial_data.model_sim_id.unique().tolist()) == len(cgr_model_comps.model_sim_id.unique().tolist()))
        assert (
            cgr_model_comps.model_name.unique().tolist() == model_comps.model_name.unique().tolist())

        # Calculate stats for each model
        for model_name in model_comps.model_name.unique().tolist():
            row = {}

            model_data = trial_data.loc[(trial_data.model_name==model_name)]
            cgr_model_data = cgr_trial_data.loc[(cgr_trial_data.model_name==model_name)]

            conditions, intercepts, cgr_conditions, cgr_intercepts = ToyModelsAnalyis.get_intercepts_conditions(
                model_name, model_comps, cgr_model_comps)

            condition_gs = ToyModelsAnalyis.get_effectsize(model_data)

            n_trials = len((model_data.model_sim_id.unique()))
            intercept_ests = intercepts.est.tolist()
            condition_ests = conditions.est.tolist()
            condition_ps = conditions.p.tolist()

            row['model'] = model_name
            row['cgr'] = round(
                model_data.loc[(model_data.condition == model_data.guess)].shape[0]/(model_data.shape[0]), 3)

            row['int'] = round(mean(intercept_ests), 2)
            row['trt'] = round(mean(condition_ests), 2)
            row['trt_p'] = round(mean(condition_ps), 3)
            row['sig_prop'] = round(sum([el <= 0.05 for el in condition_ps])/n_trials, 2)
            row['trt_g'] = round(mean(condition_gs), 2)

            # Calculate average effects across across the n_cgr_sims
            cgr_sim_ids = cgr_model_data.cgr_sim_id.unique().tolist()
            model_sim_ids = cgr_model_data.model_sim_id.unique().tolist()

            cgr_condition_ps = [mean(cgr_conditions.loc[cgr_conditions.model_sim_id==model_sim_id].p.tolist()) for model_sim_id in model_sim_ids]
            cgr_condition_ests = [mean(cgr_conditions.loc[cgr_conditions.model_sim_id==model_sim_id].est.tolist()) for model_sim_id in model_sim_ids]
            cgr_intercept_ests = [mean(cgr_intercepts.loc[cgr_intercepts.model_sim_id==model_sim_id].est.tolist()) for model_sim_id in model_sim_ids]
            cgr_condition_gs = [round(mean(ToyModelsAnalyis.get_cgr_effectsize(
                cgr_model_data.loc[cgr_model_data.model_sim_id==model_sim_id])),2) for model_sim_id in model_sim_ids]

            row['cgr_int'] = round(mean(cgr_intercept_ests), 2)
            row['cgr_trt'] = round(mean(cgr_condition_ests), 2)
            row['cgr_trt_p'] = round(mean(cgr_condition_ps), 3)
            row['cgr_sig_prop'] = round(sum([el <= 0.05 for el in cgr_condition_ps])/n_trials, 2)
            row['cgr_trt_g'] = round(mean(cgr_condition_gs), 2)

            df = df.append(row, ignore_index=True)

        df.__class__ = df_class.ToymodelsAnalysisDf
        df.set_column_types()
        df.to_csv(os.path.join(folders.summary_dir, f'{analysis_name}__summary_table.csv'), index=False)

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
            '__cgrc_data',
            '__cgrc_model_comps',
            '__cgrc_strata_contrast']

        master_df_list=[]


        for fpath in [fpath for fpath in os.listdir(target_dir) if (df_type in fpath)]:
            master_df_list.append(pd.read_csv(os.path.join(target_dir, fpath)))

        master_df = pd.concat(master_df_list, axis=0)

        if df_type == '__model_comps':
            master_df.__class__ = df_class.ModelComponentsDf
        elif df_type == '__strata_contrast':
            master_df.__class__ = df_class.StrataContrastDf
        elif df_type == '__cgrc_data':
            master_df.__class__ = df_class.CGRCurveDf
        elif df_type == '__cgrc_model_comps':
            master_df.__class__ = df_class.ModelComponentsDf
        elif df_type == '__cgrc_strata_contrast':
            master_df.__class__ = df_class.StrataContrastDf
        elif df_type == '__trial_data':
            master_df.__class__ = df_class.TrialDataDf
        else:
            assert False

        assert master_df.is_valid()
        return master_df

    @staticmethod
    def get_intercepts_conditions(model_name, model_comps, cgr_model_comps):

        conditions = model_comps.loc[
          (model_comps.model_name == model_name) &
          (model_comps.model_type == 'without_guess') &
          (model_comps.component == 'conditionAC')
        ]

        intercepts = model_comps.loc[
          (model_comps.model_name == model_name) &
          (model_comps.model_type == 'without_guess') &
          (model_comps.component == 'intercept')
        ]

        cgr_conditions= cgr_model_comps.loc[
            (cgr_model_comps.model_name == model_name) &
            (cgr_model_comps.model_type == 'without_guess') &
            (cgr_model_comps.component == 'conditionAC')
        ]

        cgr_intercepts = cgr_model_comps.loc[
            (cgr_model_comps.model_name == model_name) &
            (cgr_model_comps.model_type == 'without_guess') &
            (cgr_model_comps.component == 'intercept')
        ]

        return conditions, intercepts, cgr_conditions, cgr_intercepts

    @staticmethod
    def get_effectsize(model_data):

        model_sim_ids = model_data.model_sim_id.unique()

        gs = [round(pingouin.compute_effsize(
            model_data.loc[(model_data.condition=='AC') & (model_data.model_sim_id==model_sim_id)].delta_score,
            model_data.loc[(model_data.condition=='PL') & (model_data.model_sim_id==model_sim_id)].delta_score,
        eftype='cohen'), 3) for model_sim_id in model_sim_ids]

        return gs

    @staticmethod
    def get_cgr_effectsize(model_data):

        cgr_sim_ids = model_data.cgr_sim_id.unique()

        gs = [round(pingouin.compute_effsize(
            model_data.loc[(model_data.condition=='AC') & (model_data.cgr_sim_id==cgr_sim_id)].delta_score,
            model_data.loc[(model_data.condition=='PL') & (model_data.cgr_sim_id==cgr_sim_id)].delta_score,
        eftype='cohen'), 3) for cgr_sim_id in cgr_sim_ids]

        return gs
