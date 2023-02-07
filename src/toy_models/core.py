"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

Model family analyis naming convention:
analysis_name: {model_family}_cgrc{cgrc_param_set}
subanalysis_name: {model_family}_cgrc{cgrc_param_set}_{model_name}_trial{model_sim_id}
"""

from tqdm.contrib.itertools import product as tqdmproduct
import src.toy_models.model_defs as model_defs
from statistics import mean, median
import src.dataframe_classes as df_class
import src.constants as constants
import src.folders as folders
import src.figures as figures
import src.cgrc.core as cgrc
import src.cgrc.stats as stats
import src.miscs as miscs
import pandas as pd
import itertools
import random
import os


class Controllers():

    @staticmethod
    def run_cgrc_model_family(model_family_name, cgrc_param_set, postfix, n_patients, n_trials, save_figs=True):
        ''' Runs CGRC analysis on a family of models
            Args:
                model_family_name (str): name of model family (i.e. list of model definitions, see model_defs)
                cgrc_param_set (int): number defning which CGRC parameter set to use, see constants
                postfix (str): str appended to the end of the output files
                n_patients (int): number of patients in each simulated trail
                n_trials (int): number of trials to simulate for each member of the model family
                save_figs (bool, optional): save output figures
        '''

        assert isinstance(model_family_name, str)
        assert isinstance(cgrc_param_set, int)
        assert isinstance(postfix, str)
        assert isinstance(n_patients, int)
        assert isinstance(n_trials, int)
        assert isinstance(save_figs, bool)

        cgrc_parameters = constants.cgrc_parameters[cgrc_param_set]
        models = eval('model_defs.{}'.format(model_family_name))
        analysis_name = model_family_name + '_{}'.format(postfix)

        trial_data_dir, trial_stats_dir, cgrc_data_dir, cgrc_stats_dir, cgrc_plots_dir = miscs.create_analysis_dirs(
            analysis_name,
            trial_data_subdir=True,
            incl_cgrc_plots_dir=save_figs,
        )

        for model, model_sim_id in tqdmproduct(models, range(n_trials), desc='Model family CGRC analysis'):

            model_name = list(model.keys())[0]
            model_specs = model[model_name]
            add_columns = {'model_name': model_name, 'model_sim_id': model_sim_id}
            subanalysis_name = Helpers.get_subanalysis_name(
                analysis_name=analysis_name,
                model_name=model_name,
                model_sim_id=model_sim_id
            )

            # Generate pseudodata according to model specifications
            ToyModelsDataGenerator.get_pseudodata(
                output_dir=trial_data_dir,
                output_prefix=subanalysis_name,
                model_specs=model_specs,
                n_datapoints=n_patients,
                model_name=model_name,
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
            )

            # Make CGRC figues
            if save_figs:
                figures.Controllers.plot_VScgr_twinx(
                    input_dir=cgrc_stats_dir,
                    input_fname=subanalysis_name + '__cgrc_model_components.csv',
                    output_dir=cgrc_plots_dir,
                    output_prefix=subanalysis_name,
                )

        # Get summary table
        #summary_df = ToyModelsAnalyis.get_model_family_summary(
        #    model_family_name=model_family_name,
        #    analysis_name=analysis_name,
        #    cgrc_param_set=cgrc_param_set,
        #)
        #print('\n', summary_df.to_string(index=False))


class ToyModelsAnalyis():
    ''' Functions for model family CGRC analysis '''

    @staticmethod
    def get_model_family_summary(analysis_name, model_family_name, cgrc_param_set):
        ''' Construct summary table for model family
            Args:
                analysis_name (str): model_family_name+postfix
                model_family_name (str): name of model family (i.e. list of model definitions, see model_defs)
                cgrc_param_set (int): number defning which CGRC parameter set to use, see constants
        '''

        assert isinstance(analysis_name, str)
        assert isinstance(model_family_name, str)
        assert isinstance(cgrc_param_set, int)

        df = df_class.ModelFamilyResultsDf()

        trial_data = Helpers.get_concateneted_df_type(
            target_dir=os.path.join(folders.data_trial_data, analysis_name),
            df_type='__trial_data')

        # Get n_patients and n_trials assuming each member of the model family has same n_patient and n_trial
        tmp_mid = trial_data.model_name.to_list()[0]
        tmp_tid = trial_data.model_sim_id.to_list()[0]
        n_patients = trial_data.loc[(trial_data.model_name == tmp_mid) & (
            trial_data.model_sim_id == tmp_tid)].shape[0]
        n_trials = len(trial_data.model_sim_id.unique().tolist())

        # Get unadjusted model components
        unadj_model_components = Helpers.get_concateneted_df_type(
            target_dir=os.path.join(folders.data_trial_stats, analysis_name),
            df_type='__model_components')

        # Get adjusted model components
        cgradj_model_components = Helpers.get_concateneted_df_type(
            target_dir=os.path.join(folders.data_cgrc_stats, analysis_name),
            df_type='__cgrc_model_components')

        cgradj_model_components = cgradj_model_components.loc[(
            cgradj_model_components.cgr == 0.5)]
        assert cgradj_model_components.shape[0] > 0
        cgradj_n_trials = len(
            cgradj_model_components.model_sim_id.unique().tolist())
        assert n_trials == cgradj_n_trials
        cgradj_model_names = cgradj_model_components.model_name.unique().tolist()
        unadj_model_names = unadj_model_components.model_name.unique().tolist()
        assert cgradj_model_names == unadj_model_names

        for model_name in unadj_model_components.model_name.unique().tolist():

            filtered_unadj_model_comps = unadj_model_components.loc[
                (unadj_model_components.model_name == model_name) &
                (unadj_model_components.model_type == 'without_guess') &
                (unadj_model_components.component == 'conditionAC')
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

            filtered_cgradj_model_comps = cgradj_model_components.loc[
                (cgradj_model_components.model_name == model_name) &
                (cgradj_model_components.model_type == 'without_guess') &
                (cgradj_model_components.component == 'conditionAC')
            ]

            # Calculate average p/es aross n_cgr_trials (avg corresponds to p/es of single trial)
            model_sim_ids = filtered_cgradj_model_comps.model_sim_id.unique().tolist()
            trial_ps = []
            trial_efs = []
            for model_sim_id in model_sim_ids:
                tmp = filtered_cgradj_model_comps.loc[filtered_cgradj_model_comps.model_sim_id == model_sim_id]
                assert tmp.shape[0] == constants.cgrc_parameters[cgrc_param_set]['n_cgrc_trials']
                trial_ps.append(mean(tmp.p.tolist()))
                trial_efs.append(mean(tmp.est.tolist()))

            row['cgradj_avg_trt_p'] = round(miscs.get_estimate(trial_ps), 3)
            row['cgradj_sig_trt_rate'] = round(
                sum([p <= 0.05 for p in trial_ps])/n_trials, 3)
            row['cgradj_avg_trt_es'] = round(miscs.get_estimate(trial_efs), 3)

            df = df.append(row, ignore_index=True)

        df.__class__ = df_class.ModelFamilyResultsDf
        df.set_column_types()
        df.to_csv(os.path.join(folders.data_summary_tables, analysis_name +
                  '_{}__summary_table.csv'.format(constants.estimator)), index=False)

        return df


class ToyModelsDataGenerator():
    ''' Functions to get simulated toy model data '''

    @staticmethod
    def get_pseudodata(output_dir, output_prefix, model_specs, n_datapoints, model_name=None, model_sim_id=None, min_strata_size=4, round_digits=0):
        ''' High level function to generate pseudo-experimental toy model data
            Args:
                output_dir (str): subfolder name where results will be saved
                output_prefix (str): appended to the begining of output files
                model_specs (dict): model specifiction; see model_defs.py
                n_datapoints (int): number of datapoints in each trial
                model_name (str): name of model
                model_sim_id (int): id of current simulations
                min_strata_size (int, optional): the minimum sample size of each strata
                round_digits (int, optional): generated scores rounded to
        '''

        assert isinstance(output_dir, str)
        assert isinstance(output_prefix, str)
        assert isinstance(output_prefix, str)
        assert Helpers.is_valid_model_specs(model_specs)
        assert isinstance(n_datapoints, int)
        assert isinstance(model_sim_id, int)
        assert (min_strata_size is None) or (isinstance(min_strata_size, int))
        assert isinstance(round_digits, int)

        df = Helpers.get_TrialDatadDf(n_datapoints)

        #import pdb; pdb.set_trace()

        if min_strata_size is not None:
            df = Helpers.enforce_min_strata_size(df, min_strata_size)

        #import pdb; pdb.set_trace()

        conditions = []
        guesses = []
        scores = []
        delta_scores = []

        for idx, row in df.iterrows():

            condition, guess, score = ToyModelsDataGenerator.get_pseudodatapoint(
                oc_nh=model_specs['oc_nh'],
                gs_nh=model_specs['gs_nh'],
                se=model_specs['se'],
                dte=model_specs['dte'],
                pte=model_specs['pte'],
                ate=model_specs['ate'],
                oc2gs=model_specs['oc2gs'],
                forced_condition=row.condition,
                forced_guess=row.guess,
            )

            conditions.append(condition)
            guesses.append(guess)
            scores.append(round(score, round_digits))
            delta_scores.append(round(score, round_digits))

        df['condition'] = conditions
        df['guess'] = guesses
        df['score'] = scores
        df['delta_score'] = delta_scores


        #row.condition = condition
        #row.guess = guess
        #row.score = round(score, round_digits)
        #row.delta_score = round(score, round_digits)

        if model_name is not None:
            df['model_name'] = model_name

        if model_sim_id is not None:
            df['model_sim_id'] = model_sim_id

        df.__class__ = df_class.TrialDataDf
        df.set_column_types()
        df.to_csv(os.path.join(output_dir, output_prefix +
                  '__trial_data.csv'), index=False)

        return df

    @staticmethod
    def get_pseudodatapoint(oc_nh, gs_nh, se, dte, pte, ate, oc2gs, forced_condition=None, forced_guess=None):
        ''' Generate toy model datapoint
            Args:
                oc_nh (tuple): (mean, std) of the outcomes's natural history (in terms of change score)
                gs_nh (tuple): (mean, std) of treatment guess's natural history
                se (tuple): (mean, std) of treatment allocation's contribution to treatment guess probability
                dte (tuple): (mean, std) of treatment's contribution to outcomes
                pte (tuple): (mean, std) of placebo guess's contribution to guess
                ate (tuple): (mean, std) of active guess's contribution to guess
                oc2gs (tuple): (mean, std) of outcome's contribution to guess
                forced_condition (bool, optional): force treatment to be active/placebo
                forced_guess (bool, optional): force guess to be active/placebo
        '''

        assert isinstance(oc_nh, tuple)
        assert isinstance(gs_nh, tuple)
        assert isinstance(se, tuple)
        assert isinstance(dte, tuple)
        assert isinstance(pte, tuple)
        assert isinstance(ate, tuple)
        assert isinstance(oc2gs, tuple)
        assert forced_guess in [None, 'PL', 'AC']
        assert forced_condition in [None, 'PL', 'AC']

        # get conditon
        if forced_condition == 'PL':
            condition = 'PL'
        elif forced_condition == 'AC':
            condition = 'AC'
        elif forced_condition is None:
            condition = random.choice(['PL', 'AC'])
        else:
            assert False

        # get guess; if roll>=0.5, then guess is active
        if condition == 'PL':
            roll = random.gauss(gs_nh[0], gs_nh[1])
        elif condition == 'AC':
            roll = random.gauss(gs_nh[0], gs_nh[1])
            roll = min(1, roll)
            roll = max(0, roll)
            roll += random.gauss(se[0], se[1])
        else:
            assert False

        roll = min(1, roll)
        roll = max(0, roll)

        if (forced_guess is None) and (roll >= 0.5):
            guess = 'AC'
        elif (forced_guess is None) and (roll < 0.5):
            guess = 'PL'
        elif forced_guess is not None:
            guess = forced_guess
        else:
            assert False

        # get score
        score = random.gauss(oc_nh[0], oc_nh[1])
        if condition == 'AC':
            score += random.gauss(dte[0], dte[1])
        else:
            assert condition == 'PL'

        if guess == 'AC':
            score += random.gauss(ate[0], ate[1])
        elif guess == 'PL':
            score += random.gauss(pte[0], pte[1])
        else:
            assert False

        return condition, guess, score


class Helpers():
    ''' Various helper functions '''

    @staticmethod
    def get_TrialDatadDf(n, equal_sample_per_condition=True):
        ''' Returns TrialDataDf with aux data filled
            Args:
                n (int): number of patients, i.e. n of rows
                equal_sample_per_condition (bool, optional): force to have equal n in active / placebo arms
        '''

        assert isinstance(n, int)
        assert isinstance(equal_sample_per_condition, bool)

        df = df_class.TrialDataDf()

        if equal_sample_per_condition:
            n_pl = round(n/2)
            n_ac = n-n_pl
            df.condition = ['PL' for i in range(
                n_pl)] + ['AC' for i in range(n_ac)]
        else:
            df.condition = [None for i in range(n)]

        df.subject_id = [idx for idx in range(n)]

        df.trial = 'mock'
        df.scale = 'scale1'
        df.tp = 'wk8'
        df.baseline = 0
        df.score = None
        df.delta_score = None
        df.guess = None

        df['baseline'] = df['baseline'].astype('object')
        df['subject_id'] = df['subject_id'].astype('object')

        return df

    @staticmethod
    def enforce_min_strata_size(df, min_strata_size):
        ''' Adds rows to df to ensure that there is at least rows in each treatment/guess strata '''

        assert isinstance(min_strata_size, int)

        condition_idx = df.columns.get_loc('condition')
        guess_idx = df.columns.get_loc('guess')

        stratas = [
            ('PL', 'PL'),
            ('AC', 'PL'),
            ('PL', 'AC'),
            ('AC', 'AC'),
        ]

        idx = 0
        for strata, sample_idx in itertools.product(stratas, range(min_strata_size)):
            df.iloc[idx, condition_idx] = strata[0]
            df.iloc[idx, guess_idx] = strata[1]
            idx += 1

        return df

    @staticmethod
    def get_subanalysis_name(analysis_name, model_name, model_sim_id):
        ''' Return name of subanalysis '''

        assert isinstance(analysis_name, str)
        assert isinstance(model_name, str)
        assert isinstance(model_sim_id, int)

        return '{}_{}_trial{}'.format(analysis_name, model_name, model_sim_id)

    @staticmethod
    def get_concateneted_df_type(target_dir, df_type):
        ''' Concatenets and returns all dataframes of a certain kind from target_dir
            Args:
                target_dir(str): path to folder
                df_type(str): dataframe type; see assertions for valid types
        '''

        assert df_type in [
            '__trial_data',
            '__model_components',
            '__strata_contrast',
            '__cgrc_model_components',
            '__cgrc_strata_contrast']

        # Get all trial stats
        if df_type == '__model_components':
            master_df = df_class.ModelComponentsDf()
        elif df_type == '__strata_contrast':
            master_df = df_class.StrataContrastDf()
        elif df_type == '__cgrc_model_components':
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

        if df_type == '__model_components':
            master_df.__class__ = df_class.ModelComponentsDf
        elif df_type == '__strata_contrast':
            master_df.__class__ = df_class.StrataContrastDf
        elif df_type == '__cgrc_model_components':
            master_df.__class__ = df_class.ModelComponentsDf
        elif df_type == '__cgrc_strata_contrast':
            master_df.__class__ = df_class.StrataContrastDf
        elif df_type == '__trial_data':
            master_df.__class__ = df_class.TrialDataDf
        else:
            assert False

        master_df.set_column_types()
        return master_df

    @staticmethod
    def is_valid_model_specs(model):
        ''' Check if dict is a valid model parameters dictionary
            Returns True if it is; otherwise will throw False Assertion
            Args:
                model(dict): dictionary containing model parameters
        '''

        assert isinstance(model, dict)
        assert len(model.keys()) == 7
        assert all([p in model.keys()
                   for p in ['oc_nh', 'gs_nh', 'se', 'dte', 'pte', 'ate', 'oc2gs']])

        for model_parameter in model.values():
            assert isinstance(model_parameter, tuple)
            assert len(model_parameter) == 2
            assert all([(isinstance(el, float) or isinstance(el, int))
                       for el in model_parameter])

        return True
