"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import src.dataframe_classes as df_class
from statistics import mean, median
import src.config as config
import pandas as pd
import src.folders as folders
import pingouin
import os


def create_dir(target_dir):
    ''' ensures target_dir exists '''

    if os.path.isdir(target_dir):
        pass
    else:
        os.makedirs(target_dir)


def create_empty_dir(target_dir):
    ''' ensures target_dir exists and that it is empty '''

    if os.path.isdir(target_dir):
        [os.remove(os.path.join(target_dir, filename)) for filename in os.listdir(target_dir)]
    else:
        os.mkdir(target_dir)


def get_trial_scales(input_df, trial_scales):
    ''' return info about trial '''

    assert (isinstance(trial_scales, dict) or (trial_scales is None))
    assert (isinstance(input_df, pd.DataFrame) or (input_df is None))

    if trial_scales is None:
        trials = input_df.trial.unique().tolist()
        scales = input_df.scale.unique().tolist()
    else:
        trials = list(trial_scales.keys())
        scales = list(trial_scales.values())[0]

    return trials, scales


def create_analysis_dirs(analysis_name, trial_data_subdir=False, incl_cgrc_plots_dir=True):
    ''' creates all empty directory for analysis '''

    assert isinstance(trial_data_subdir, bool)
    assert isinstance(incl_cgrc_plots_dir, bool)

    if trial_data_subdir:
        trial_data_dir = os.path.abspath(os.path.join(folders.trial_data_dir, analysis_name))
        create_empty_dir(trial_data_dir)
    else:
        trial_data_dir = os.path.abspath(folders.trial_data_dir)
        create_dir(trial_data_dir)

    trial_stats_dir = os.path.abspath(os.path.join(
        folders.trial_stats_dir, analysis_name))
    create_empty_dir(trial_stats_dir)

    cgrc_data_dir = os.path.abspath(os.path.join(
        folders.cgrc_data_dir, analysis_name))
    create_empty_dir(cgrc_data_dir)

    cgrc_stats_dir = os.path.abspath(os.path.join(
        folders.cgrc_stats_dir, analysis_name))
    create_empty_dir(cgrc_stats_dir)

    if incl_cgrc_plots_dir:
        cgrc_plots_dir = os.path.abspath(
            os.path.join(folders.cgrc_plots, analysis_name))
        create_empty_dir(cgrc_plots_dir)
    else:
        cgrc_plots_dir = None

    return trial_data_dir, trial_stats_dir, cgrc_data_dir, cgrc_stats_dir, cgrc_plots_dir


def get_estimate(data, metric=config.estimator):
    ''' choose which estimator '''

    assert metric in ['mean', 'median']

    if metric == 'mean':
        return mean(data)

    if metric == 'median':
        return median(data)

    assert False


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

    return master_df


def get_aeb_summary_table(analysis_name):
    ''' Construct summary table for model family
        Args:
            analysis_name (str):
    '''

    assert isinstance(analysis_name, str)

    # Get unadjusted and adjusted model components
    model_comps = get_concateneted_df_type(
      target_dir=os.path.join(folders.trial_stats_dir, analysis_name),
      df_type='__model_comps')
    cgr_model_comps = get_concateneted_df_type(
      target_dir=os.path.join(folders.cgrc_stats_dir, analysis_name),
      df_type='__cgrc_model_comps')
    cgr_model_comps = cgr_model_comps.loc[(cgr_model_comps.cgr == 0.5)]

    # Sanity checks
    assert (
        cgr_model_comps.shape[0] > 0)
    assert (
        cgr_model_comps.model_name.unique().tolist() == model_comps.model_name.unique().tolist())

    df = pd.DataFrame(columns=[
                'model',
                'trt',
                'trt_p',
                'sig_prop',
                'cgr_trt',
                'cgr_trt_p',
                'cgr_sig_prop',
            ])

    # Calculate stats for each model
    for model_name in model_comps.model_name.unique().tolist():

        row={}

        conditions = model_comps.loc[
          (model_comps.model_name == model_name) &
          (model_comps.model_type == 'without_guess') &
          (model_comps.component == 'conditionAC')
        ]

        n_trials =  max(model_comps.model_sim_id.tolist())+1

        condition_ests = conditions.est.tolist()
        condition_ps = conditions.p.tolist()

        row['model'] = model_name
        row['trt'] = round(mean(condition_ests), 3)
        row['trt_p'] = round(mean(condition_ps), 3)
        row['sig_prop'] = round(sum([el <= 0.05 for el in condition_ps])/n_trials, 2)

        # Calculate average effects across across the n_cgr_sims
        cgr_conditions= cgr_model_comps.loc[
            (cgr_model_comps.model_name == model_name) &
            (cgr_model_comps.model_type == 'without_guess') &
            (cgr_model_comps.component == 'conditionAC')
        ]

        cgr_sim_ids = cgr_model_comps.cgr_sim_id.unique().tolist()
        model_sim_ids = cgr_model_comps.model_sim_id.unique().tolist()

        cgr_condition_ests = [mean(cgr_conditions.loc[cgr_conditions.model_sim_id==model_sim_id].est.tolist()) for model_sim_id in model_sim_ids]
        cgr_condition_ps = [mean(cgr_conditions.loc[cgr_conditions.model_sim_id==model_sim_id].p.tolist()) for model_sim_id in model_sim_ids]

        row['cgr_trt'] = round(mean(cgr_condition_ests), 2)
        row['cgr_trt_p'] = round(mean(cgr_condition_ps), 3)
        row['cgr_sig_prop'] = round(sum([el <= 0.05 for el in cgr_condition_ps])/n_trials, 2)

        df = df.append(row, ignore_index=True) # TODO remove .append from loop to optimize performance

    return df
