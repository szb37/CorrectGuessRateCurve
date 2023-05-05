"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

from sklearn.neighbors import KernelDensity
from itertools import product as product
import src.dataframe_classes as df_class
import src.config as config
import src.folders as folders
import src.figures as figures
import src.cgrc.stats as stats
import src.miscs as miscs
import pandas as pd
import numpy as np
import pingouin
import random
import os


class Controllers():

    @staticmethod
    def run_cgrc_trial(trial_name, postfix, cgrc_param_set, do_stratas=False, trial_scales=None, save_figs=False):
        ''' Run CGRC pipeline on a TrialDataDf
            Args:
                trial_name (str): name of trial
                cgrc_param_set (int): which paramter set used, see config
                do_stratas (bool, optional): calculate starta compariosn stats or not.
                    By default these are not needed for downstream computation, but maybe useful in other contexts
                trial_scales (dict, optional): key-values of trials-scales for which stats will be calculated,
                    e.g. trial_scales = {'trial2':['scale1', 'scale3']}, will process scales 1 and 3 of trial2.
                    If left None, will automatically calculate results for all scales across all trials
                save_figs (bool, optional): save figures or not
        '''

        assert isinstance(trial_name, str)
        assert isinstance(postfix, str)
        assert isinstance(cgrc_param_set, str)
        assert (isinstance(trial_scales, dict) or (trial_scales is None))

        cgrc_parameters = config.cgrc_parameters[cgrc_param_set]
        analysis_name = trial_name + '_{}'.format(postfix)

        trial_data_dir, trial_stats_dir, cgrc_data_dir, cgrc_stats_dir, cgrc_plots_dir = miscs.create_analysis_dirs(
            analysis_name,
            incl_cgrc_plots_dir=save_figs,
        )

        # Get trial stats
        stats.Controllers.get_trial_stats(
            input_dir=trial_data_dir,
            input_fname=trial_name + '__trial_data.csv',
            output_dir=trial_stats_dir,
            output_prefix=analysis_name,
            trial_scales=trial_scales,
        )

        # Get CGRC from trial data
        CorrectGuessRateCurve.get_cgrc_data(
            input_dir=trial_data_dir,
            input_fname=trial_name + '__trial_data.csv',
            output_dir=cgrc_data_dir,
            output_prefix=analysis_name,
            cgrc_parameters=cgrc_parameters,
            trial_scales=trial_scales,
        )

        # Get CGRC stats
        stats.Controllers.get_cgrc_stats(
            input_dir=cgrc_data_dir,
            input_fname=analysis_name + '__cgrc_data.csv',
            output_dir=cgrc_stats_dir,
            output_prefix=analysis_name,
            trial_scales=trial_scales,
            do_stratas=do_stratas,
        )

        # Make CGRC figues
        if save_figs:
            figures.Controllers.plot_VScgr_twinx(
                input_dir=cgrc_stats_dir,
                input_fname=analysis_name + '__cgrc_model_comps.csv',
                output_dir=cgrc_plots_dir,
                output_prefix=analysis_name,
                trial_scales=trial_scales,
            )

    @staticmethod
    def get_cgrc_comparison_table_v1(trial_name, analysis_name, trial_data_dir, trial_stats_dir, cgrc_data_dir, cgrc_stats_dir, trial_scales=None,):
        ''' Compute and save summary table
            Args:
            trial_name (str): name of trial
            analysis_name (str): trial name_postfix - cna be used distinct analysis of the same datset
            trial_scales (dict, optional): key-values of trials-scales for which stats will be calculated,
                e.g. trial_scales = {'trial2':['scale1', 'scale3']}, will process scales 1 and 3 of trial2.
                If left None, will automatically calculate results for all scales across all trials
         '''

        assert isinstance(trial_name, str)
        assert isinstance(analysis_name, str)
        assert isinstance(trial_data_dir, str)
        assert isinstance(trial_stats_dir, str)
        assert isinstance(cgrc_data_dir, str)
        assert isinstance(cgrc_stats_dir, str)
        assert (isinstance(trial_scales, dict) or (trial_scales is None))

        df = pd.DataFrame(
            columns=[
                'scale',
                'trial_cgr',
                'trial_est',
                'trial_p',
                'trial_g',
                'cgr_cgr',
                'cgrc_est',
                'cgrc_p',
                'cgrc_g',
            ])

        trial_data = pd.read_csv(os.path.join(
            trial_data_dir, trial_name+'__trial_data.csv'))

        trial_stats = pd.read_csv(os.path.join(
            trial_stats_dir, analysis_name+'__model_comps.csv'))
        trial_stats = trial_stats.loc[
            (trial_stats.model_type == 'without_guess') &
            (trial_stats.component == 'conditionAC')
        ]

        cgrc_data = pd.read_csv(os.path.join(
            cgrc_data_dir, analysis_name+'__cgrc_data.csv'))
        cgrc_data = cgrc_data.loc[
            (cgrc_data.cgr == 0.5)
        ]

        cgrc_stats = pd.read_csv(os.path.join(
            cgrc_stats_dir, analysis_name+'__cgrc_model_comps.csv'))
        cgrc_stats = cgrc_stats.loc[
            (cgrc_stats.cgr == 0.5) &
            (cgrc_stats.model_type == 'without_guess') &
            (cgrc_stats.component == 'conditionAC')
        ]

        assert cgrc_stats.scale.unique().tolist() == trial_stats.scale.unique().tolist()

        for scale in trial_stats.scale.unique().tolist():

            tmp_trial_data = trial_data.loc[(trial_data.scale == scale)]
            tmp_trial_stats = trial_stats.loc[(trial_stats.scale == scale)]
            tmp_cgrc_data = cgrc_data.loc[(cgrc_data.scale == scale)]
            tmp_cgrc_stats = cgrc_stats.loc[(cgrc_stats.scale == scale)]

            # Compute trial effect size
            trial_hedges_g = pingouin.compute_effsize(
                tmp_trial_data.loc[tmp_trial_data.condition ==
                                   'AC'].delta_score,
                tmp_trial_data.loc[tmp_trial_data.condition ==
                                   'PL'].delta_score,
                eftype='hedges'
            )

            # Compute CGRC effect size
            cgrc_hedges_g = []
            for cgr_sim_id in cgrc_data.cgr_sim_id.unique().tolist():
                cgrc_hedges_g.append(
                    pingouin.compute_effsize(
                        tmp_cgrc_data.loc[(tmp_cgrc_data.condition == 'AC') & (
                            tmp_cgrc_data.cgr_sim_id == cgr_sim_id)].delta_score,
                        tmp_cgrc_data.loc[(tmp_cgrc_data.condition == 'PL') & (
                            tmp_cgrc_data.cgr_sim_id == cgr_sim_id)].delta_score,
                        eftype='hedges'
                    ))

            # Compute trial CGR
            n_plpl = tmp_trial_data.loc[(tmp_trial_data.condition == 'PL') & (
                tmp_trial_data.guess == 'PL')].shape[0]
            n_acac = tmp_trial_data.loc[(tmp_trial_data.condition == 'AC') & (
                tmp_trial_data.guess == 'AC')].shape[0]
            trial_cgr = (n_plpl+n_acac)/tmp_trial_data.shape[0]

            # Compute CGRC CGR
            cgrc_cgrs = []
            for cgr_sim_id in cgrc_data.cgr_sim_id.unique().tolist():
                n_plpl = tmp_cgrc_data.loc[(tmp_cgrc_data.condition == 'PL') & (
                    tmp_cgrc_data.guess == 'PL') & (tmp_cgrc_data.cgr_sim_id == cgr_sim_id)].shape[0]
                n_acac = tmp_cgrc_data.loc[(tmp_cgrc_data.condition == 'AC') & (
                    tmp_cgrc_data.guess == 'AC') & (tmp_cgrc_data.cgr_sim_id == cgr_sim_id)].shape[0]
                cgrc_cgrs.append(
                    (n_plpl+n_acac)/tmp_cgrc_data.loc[(tmp_cgrc_data.cgr_sim_id == cgr_sim_id)].shape[0])

            # Add results to output table
            row = {}
            row['scale'] = scale
            row['trial_cgr'] = round(trial_cgr, 2)
            row['trial_est'] = '{}±{}'.format(
                round(miscs.get_estimate(tmp_trial_stats.est.tolist()), 1),
                round(miscs.get_estimate(tmp_trial_stats.se.tolist()), 1)
            )
            row['trial_p'] = round(miscs.get_estimate(
                tmp_trial_stats.p.tolist()), 4)
            row['trial_g'] = round(trial_hedges_g, 2)
            row['cgr_cgr'] = round(miscs.get_estimate(cgrc_cgrs), 2)
            row['cgrc_est'] = '{}±{}'.format(
                round(miscs.get_estimate(tmp_cgrc_stats.est.tolist()), 1),
                round(miscs.get_estimate(tmp_cgrc_stats.se.tolist()), 1)
            )
            row['cgrc_p'] = round(miscs.get_estimate(
                tmp_cgrc_stats.p.tolist()), 4)
            row['cgrc_g'] = round(miscs.get_estimate(cgrc_hedges_g), 1)
            df = df.append(row, ignore_index=True)

        df.to_csv(os.path.join(folders.summary_dir, f'{analysis_name}__summary_table_v1.csv'), index=False)
        return df


class CorrectGuessRateCurve():
    ''' Calculate correct guess rate curve '''

    @staticmethod
    def get_cgrc_data(input_dir, input_fname, output_dir, output_prefix, cgrc_parameters, strata_sampling='all_prop', trial_scales=None, add_columns=None, n_sample=None):
        ''' Get CGRC data
            Args:
                input_dir (str): folder of input file
                input_fname (str): filename of input file
                output_dir (str): where output is written
                output_prefix (str): prefix of output files
                cgrc_parameters (dict): dictionary of CGRC parameters; see config
                trial_scales (dict, optional): key-values of trials-scales for which stats will be calculated,
                    e.g. trial_scales = {'trial2':['scale1', 'scale3']}, will process scales 1 and 3 of trial2.
                    If left None, will automatically calculate results for all scales across all trial
                strata_sampling (str, optional): method how to assign sample size between strata;
                    must be in ['all_prop', 'all_equal', 'active_equal', 'active_prop']
        '''

        assert isinstance(input_dir, str)
        assert isinstance(input_fname, str)
        assert isinstance(output_dir, str)
        assert isinstance(output_prefix, str)
        assert isinstance(cgrc_parameters, dict)
        assert (isinstance(trial_scales, dict) or (trial_scales is None))
        assert (isinstance(add_columns, dict)) or (add_columns is None)
        assert isinstance(n_sample, int) or (n_sample is None)

        trial_data_df = pd.read_csv(os.path.join(input_dir, input_fname))
        trials, scales = miscs.get_trial_scales(input_df=trial_data_df, trial_scales=trial_scales)
        n_cgrc_trials = cgrc_parameters['n_cgrc_trials']

        # CGR values need to be rounded, otherwise downstream R df filtering may not work
        cgr_values = [round(el, 5) for el in cgrc_parameters['cgr_values']]

        cgrc_dfs=[]

        for trial, scale in product(trials, scales):

            df_filtered = trial_data_df.loc[(trial_data_df.scale == scale)]
            df_filtered = df_filtered.loc[(df_filtered.trial == trial)]
            if df_filtered.shape[0] == 0:
                continue

            if n_sample is None:
                total_sample_size = df_filtered.shape[0]
            else:
                total_sample_size = n_sample

            kdes = CorrectGuessRateCurve.get_kdes(df_filtered=df_filtered)

            for cgr, cgr_sim_id in product(cgr_values, range(n_cgrc_trials)):

                sample_sizes = CorrectGuessRateCurve.get_strata_sample_sizes(
                    total_sample_size=total_sample_size,
                    df_filtered=df_filtered,
                    correct_guess_rate=cgr,
                    strata_sampling=strata_sampling,)

                cgrc_datapoint_df = CorrectGuessRateCurve.get_cgrc_datapoint_KDE(
                    df_filtered=df_filtered,
                    sample_sizes=sample_sizes,
                    kdes=kdes,
                    cgr=cgr,)

                cgrc_datapoint_df.cgr_sim_id = cgr_sim_id
                cgrc_datapoint_df.trial = trial
                cgrc_datapoint_df.scale = scale
                cgrc_dfs.append(cgrc_datapoint_df)

                del cgrc_datapoint_df, sample_sizes

        master_cgrc_df = pd.concat(cgrc_dfs, axis=0)
        master_cgrc_df.__class__ = df_class.CGRCurveDf
        master_cgrc_df.add_columns(add_columns)
        master_cgrc_df.set_column_types()

        master_cgrc_df.to_csv(os.path.join(
            output_dir, output_prefix+'__cgrc_data.csv'), index=False)

    @staticmethod
    def get_cgrc_datapoint_KDE(df_filtered, sample_sizes, cgr, kdes):
        ''' Generate CGRC datapoints from KDEs
            Args:
                df_filtered (pd.DataFrame):
                sample_sizes (dict): sample size of each strata
                kdes (dict): dictionary of KDE objects, where keys are strata labels
                cgr(float): correct guess rate
        '''

        assert isinstance(sample_sizes, dict)
        assert isinstance(kdes, dict)

        plpl_df = df_class.CGRCurveDf()
        acpl_df = df_class.CGRCurveDf()
        plac_df = df_class.CGRCurveDf()
        acac_df = df_class.CGRCurveDf()

        plpl_scores = [round(score) for score in kdes['PLPL'].sample(
            sample_sizes['PLPL']).reshape(1, -1).tolist()[0]]
        acpl_scores = [round(score) for score in kdes['ACPL'].sample(
            sample_sizes['ACPL']).reshape(1, -1).tolist()[0]]
        plac_scores = [round(score) for score in kdes['PLAC'].sample(
            sample_sizes['PLAC']).reshape(1, -1).tolist()[0]]
        acac_scores = [round(score) for score in kdes['ACAC'].sample(
            sample_sizes['ACAC']).reshape(1, -1).tolist()[0]]

        plpl_df.delta_score = plpl_scores
        acpl_df.delta_score = acpl_scores
        plac_df.delta_score = plac_scores
        acac_df.delta_score = acac_scores

        plpl_df.condition = 'PL'
        acpl_df.condition = 'AC'
        plac_df.condition = 'PL'
        acac_df.condition = 'AC'

        plpl_df.guess = 'PL'
        acpl_df.guess = 'PL'
        plac_df.guess = 'AC'
        acac_df.guess = 'AC'

        df = pd.concat([plpl_df, acpl_df, plac_df, acac_df,], sort=False)
        df.cgr = cgr
        return df

    @staticmethod
    def get_kdes(df_filtered):
        ''' Get strata KDEs '''

        kdes = {}

        plpl_scores = df_filtered.loc[(df_filtered.condition == 'PL') & (
            df_filtered.guess == 'PL'), 'delta_score']
        kdes['PLPL'] = KernelDensity(kernel='gaussian').fit(
            np.array(plpl_scores).reshape(-1, 1))

        acpl_scores = df_filtered.loc[(df_filtered.condition == 'AC') & (
            df_filtered.guess == 'PL'), 'delta_score']
        kdes['ACPL'] = KernelDensity(kernel='gaussian').fit(
            np.array(acpl_scores).reshape(-1, 1))

        plac_scores = df_filtered.loc[(df_filtered.condition == 'PL') & (
            df_filtered.guess == 'AC'), 'delta_score']
        kdes['PLAC'] = KernelDensity(kernel='gaussian').fit(
            np.array(plac_scores).reshape(-1, 1))

        acac_scores = df_filtered.loc[(df_filtered.condition == 'AC') & (
            df_filtered.guess == 'AC'), 'delta_score']
        kdes['ACAC'] = KernelDensity(kernel='gaussian').fit(
            np.array(acac_scores).reshape(-1, 1))

        return kdes

    @staticmethod
    def get_strata_ratio(df_filtered):
        """ Returns the ratio of each strata in df """

        strata_ratio = {'PLPL': 0, 'ACPL': 0, 'PLAC': 0, 'ACAC': 0, }
        n_all = df_filtered.shape[0]

        ss_plpl = df_filtered.loc[(df_filtered.condition == 'PL') & (
            df_filtered.guess == 'PL')].shape[0]
        ss_acpl = df_filtered.loc[(df_filtered.condition == 'AC') & (
            df_filtered.guess == 'PL')].shape[0]
        ss_plac = df_filtered.loc[(df_filtered.condition == 'PL') & (
            df_filtered.guess == 'AC')].shape[0]
        ss_acac = df_filtered.loc[(df_filtered.condition == 'AC') & (
            df_filtered.guess == 'AC')].shape[0]

        strata_ratio['PLPL'] = round(ss_plpl/n_all, 2)
        strata_ratio['ACPL'] = round(ss_acpl/n_all, 2)
        strata_ratio['PLAC'] = round(ss_plac/n_all, 2)
        strata_ratio['ACAC'] = round(ss_acac/n_all, 2)

        return strata_ratio

    @staticmethod
    def get_strata_sample_sizes(total_sample_size, df_filtered, correct_guess_rate, strata_sampling='all_prop'):
        ''' Get number of datapoints in each strata
            Args:
                total_sample_size (int):
                df_filtered (pd.DataFrame):
                correct_guess_rate (float):
                strata_sampling (str, optional): strata sampling method
        '''

        assert strata_sampling in ['all_prop',
                                   'all_equal', 'active_equal', 'active_prop']
        sample_sizes = {}

        if strata_sampling == 'all_prop':
            # sample_size_bb: sample size for blind breaking
            sample_size_bb = round(correct_guess_rate*total_sample_size)
            # sample_size_nbb: sample size for NON blind breaking
            sample_size_nbb = total_sample_size-sample_size_bb
            strata_ratio = CorrectGuessRateCurve.get_strata_ratio(df_filtered)

            # Ratio of PLPL within BB (break blind) cases
            plpl_correct_guess_rate = strata_ratio['PLPL'] / \
                (strata_ratio['PLPL']+strata_ratio['ACAC'])
            sample_sizes['PLPL'] = round(
                sample_size_bb*(plpl_correct_guess_rate))
            sample_sizes['ACAC'] = sample_size_bb-sample_sizes['PLPL']

            # Ratio of ACPL within NBB (non break blind) cases
            acpl_ncorrect_guess_rate = strata_ratio['ACPL'] / \
                (strata_ratio['ACPL']+strata_ratio['PLAC'])
            sample_sizes['ACPL'] = round(
                sample_size_nbb*(acpl_ncorrect_guess_rate))
            sample_sizes['PLAC'] = sample_size_nbb-sample_sizes['ACPL']

        elif strata_sampling == 'all_equal':
            # sample_size_bb: sample size for blind breaking
            sample_size_bb = round(correct_guess_rate*total_sample_size)
            # sample_size_nbb: sample size for NON blind breaking
            sample_size_nbb = total_sample_size-sample_size_bb
            sample_sizes['PLPL'] = round(sample_size_bb/2)
            sample_sizes['ACAC'] = sample_size_bb-sample_sizes['PLPL']

            sample_sizes['ACPL'] = round(sample_size_nbb/2)
            sample_sizes['PLAC'] = sample_size_nbb-sample_sizes['ACPL']

        elif strata_sampling == 'active_equal':

            # sample_size of active stratas
            acs_sample_size = round(total_sample_size/2)
            pls_sample_size = total_sample_size - \
                acs_sample_size  # sample_size of placebo stratas

            sample_sizes['PLPL'] = round(pls_sample_size/2)
            sample_sizes['PLAC'] = pls_sample_size-sample_sizes['PLPL']

            sample_sizes['ACAC'] = round(correct_guess_rate*acs_sample_size)
            sample_sizes['ACPL'] = acs_sample_size-sample_sizes['ACAC']

        elif strata_sampling == 'active_prop':

            strata_ratio = CorrectGuessRateCurve.get_strata_ratio(df_filtered)

            acs_sample_size = round(
                total_sample_size*(strata_ratio['ACPL']+strata_ratio['ACAC']))
            pls_sample_size = total_sample_size-acs_sample_size

            sample_sizes['PLPL'] = round(
                pls_sample_size*(strata_ratio['PLPL']/(strata_ratio['PLPL']+strata_ratio['PLAC'])))
            sample_sizes['PLAC'] = pls_sample_size - sample_sizes['PLPL']

            sample_sizes['ACAC'] = round(correct_guess_rate*acs_sample_size)
            sample_sizes['ACPL'] = acs_sample_size-sample_sizes['ACAC']

        else:
            assert False

        return sample_sizes
