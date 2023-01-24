"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2020, DrugNerdsLab
:License: MIT
"""

from statistics import mean, median
import src.constants as constants
import src.folders as folders
import pandas as pd
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
        [os.remove(os.path.join(target_dir, filename))
         for filename in os.listdir(target_dir)]
    else:
        os.mkdir(target_dir)


def get_study_scales(input_df, study_scales):
    ''' return info about study '''

    assert (isinstance(study_scales, dict) or (study_scales is None))
    assert (isinstance(input_df, pd.DataFrame) or (input_df is None))

    respondents = constants.respondents
    guessers = constants.guessers

    if study_scales is None:
        studies = input_df.study.unique().tolist()
        scales = input_df.scale.unique().tolist()
    else:
        studies = list(study_scales.keys())
        scales = list(study_scales.values())[0]

    return studies, scales, respondents, guessers


def create_analysis_dirs(analysis_name, trial_data_subdir=False, incl_cgrc_plots_dir=True):
    ''' creates all empty directory for analysis '''

    assert isinstance(trial_data_subdir, bool)
    assert isinstance(incl_cgrc_plots_dir, bool)

    if trial_data_subdir:
        trial_data_dir = os.path.abspath(os.path.join(
            folders.data_trial_data, analysis_name))
        create_dir(trial_data_dir)
    else:
        trial_data_dir = os.path.abspath(folders.data_trial_data)

    trial_stats_dir = os.path.abspath(os.path.join(
        folders.data_trial_stats, analysis_name))
    create_dir(trial_stats_dir)

    cgrc_data_dir = os.path.abspath(os.path.join(
        folders.data_cgrc_data, analysis_name))
    create_dir(cgrc_data_dir)

    cgrc_stats_dir = os.path.abspath(os.path.join(
        folders.data_cgrc_stats, analysis_name))
    create_dir(cgrc_stats_dir)

    if incl_cgrc_plots_dir:
        cgrc_plots_dir = os.path.abspath(
            os.path.join(folders.cgrc_plots, analysis_name))
        create_dir(cgrc_plots_dir)
    else:
        cgrc_plots_dir = None

    return trial_data_dir, trial_stats_dir, cgrc_data_dir, cgrc_stats_dir, cgrc_plots_dir


def get_estimate(data, metric=constants.estimator):
    ''' choose which estimator '''

    assert metric in ['mean', 'median']

    if metric == 'mean':
        return mean(data)

    if metric == 'median':
        return median(data)

    assert False
