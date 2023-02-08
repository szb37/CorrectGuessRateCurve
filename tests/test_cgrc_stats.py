"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import src.folders as folders
import src.dataframe_classes as df_class
import src.cgrc.stats as stats
from rpy2.robjects import r
from unittest import mock
import pandas as pd
import unittest
import os


class StatsUnitTests(unittest.TestCase):

    # Define mock dataframes
    model_summary = df_class.ModelSummaryDf()
    model_summary.add_columns({'cgr': None, 'cgr_sim_id': None})
    model_summary = model_summary.append({
        'trial': 'test',
        'scale': 'foo',
        'cgr': 50,
        'cgr_sim_id': 50,
        'model_type': 'test_type',
        'df1': 1,
        'df2': 2,
        'f': 3,
        'adjr2': 4,
        'p': 5},
        ignore_index=True)

    model_components = df_class.ModelComponentsDf()
    model_components.add_columns({'cgr': None, 'cgr_sim_id': None})
    model_components = model_components.append({
        'trial': 'test',
        'scale': 'foo',
        'cgr': 50,
        'cgr_sim_id': 50,
        'model_type': 'test_type',
        'component': 'tadaa',
        'est': 1,
        'se': 2,
        't': 3,
        'p': 4},
        ignore_index=True)

    strata_summary = df_class.StrataSummaryDf()
    strata_summary.add_columns({'cgr': None, 'cgr_sim_id': None})
    strata_summary = strata_summary.append({
        'trial': 'test',
        'scale': 'foo',
        'cgr': 40,
        'cgr_sim_id': 50,
        'strata': 'tadaa',
        'est': 4,
        'se': 5,
        'df': 6,
        'lower_CI': 7,
        'upper_CI': 8},
        ignore_index=True)

    strata_contrast = df_class.StrataContrastDf()
    strata_contrast.add_columns({'cgr': None, 'cgr_sim_id': None})
    strata_contrast = strata_contrast.append({
        'trial': 'test',
        'scale': 'foo',
        'cgr': 60,
        'cgr_sim_id': 50,
        'contrast': 'tadaa',
        'type': 'tadaaa',
        'est': 10,
        'se': 20,
        'df': 30,
        't': 40,
        'p': 50,
        'p_adj': 60},
        ignore_index=True)

    all_dfs = {
        'model_summary': model_summary,
        'model_components': model_components,
        'strata_summary': strata_summary,
        'strata_contrast': strata_contrast
    }

    model_summary2 = df_class.ModelSummaryDf()
    model_summary2.add_columns({'cgr': None, 'cgr_sim_id': None})
    model_summary2 = model_summary2.append({
        'trial': 'test',
        'scale': 'foo',
        'cgr': 50,
        'cgr_sim_id': 50,
        'model_type': 'test_type',
        'df1': 5,
        'df2': 6,
        'f': 7,
        'adjr2': 8,
        'p': 9},
        ignore_index=True)

    model_components2 = df_class.ModelComponentsDf()
    model_components2.add_columns({'cgr': None, 'cgr_sim_id': None})
    model_components2 = model_components2.append({
        'trial': 'test',
        'scale': 'foo',
        'cgr': 50,
        'cgr_sim_id': 50,
        'model_type': 'test_type',
        'component': 'tadaa',
        'est': 4,
        'se': 5,
        't': 6,
        'p': 7},
        ignore_index=True)

    strata_summary2 = df_class.StrataSummaryDf()
    strata_summary2.add_columns({'cgr': None, 'cgr_sim_id': None})
    strata_summary2 = strata_summary2.append({
        'trial': 'test',
        'scale': 'foo',
        'cgr': 40,
        'cgr_sim_id': 50,
        'strata': 'tadaa',
        'est': 1,
        'se': 2,
        'df': 3,
        'lower_CI': 4,
        'upper_CI': 5},
        ignore_index=True)

    strata_contrast2 = df_class.StrataContrastDf()
    strata_contrast2.add_columns({'cgr': None, 'cgr_sim_id': None})
    strata_contrast2 = strata_contrast2.append({
        'trial': 'test',
        'scale': 'foo',
        'cgr': 60,
        'cgr_sim_id': 50,
        'contrast': 'tadaa',
        'type': 'tadaaa',
        'est': 1,
        'se': 2,
        'df': 3,
        't': 4,
        'p': 5,
        'p_adj': 6},
        ignore_index=True)

    all_dfs2 = {
        'model_summary': model_summary2,
        'model_components': model_components2,
        'strata_summary': strata_summary2,
        'strata_contrast': strata_contrast2
    }

    def test_get_df_filtered(self):

        input_fpath = os.path.join(folders.fixtures, 'get_df_filtered_input.csv').replace('\\', '/')
        stats.Helpers.load_df_into_R_space(input_fpath)

        stats.Helpers.get_df_filtered(trial='a', scale='foo')
        temp = stats.Helpers.r2pyjson('df_filtered')
        self.assertEqual(temp['delta_score'], [0, 3, 4, 5, 6, 7])

        stats.Helpers.get_df_filtered(trial='b', scale='foo')
        temp = stats.Helpers.r2pyjson('df_filtered')
        self.assertEqual(temp['delta_score'], [1, 8])

        stats.Helpers.get_df_filtered(trial='a', scale='tadaa')
        temp = stats.Helpers.r2pyjson('df_filtered')
        self.assertEqual(temp['delta_score'], [2])

        stats.Helpers.get_df_filtered(trial='a', scale='foo', cgr=0, cgr_sim_id=0)
        temp = stats.Helpers.r2pyjson('df_filtered')
        self.assertEqual(temp['delta_score'], [0, 3, 4])

        stats.Helpers.get_df_filtered(trial='a', scale='foo', cgr=1, cgr_sim_id=0)
        temp = stats.Helpers.r2pyjson('df_filtered')
        self.assertEqual(temp['delta_score'], [5])

        stats.Helpers.get_df_filtered(trial='a', scale='foo', cgr=0, cgr_sim_id=1)
        temp = stats.Helpers.r2pyjson('df_filtered')
        self.assertEqual(temp['delta_score'], [6])

        stats.Helpers.get_df_filtered(trial='a', scale='foo', cgr=1, cgr_sim_id=1)
        temp = stats.Helpers.r2pyjson('df_filtered')
        self.assertEqual(temp['delta_score'], [7])

        stats.Helpers.get_df_filtered(trial='all', scale='foo')
        temp = stats.Helpers.r2pyjson('df_filtered')
        self.assertEqual(temp['delta_score'], [0, 1, 3, 4, 5, 6, 7, 8, 9])

        stats.Helpers.get_df_filtered(trial='all', scale='tadaa')
        temp = stats.Helpers.r2pyjson('df_filtered')
        self.assertEqual(temp['delta_score'], [2])

        stats.Helpers.get_df_filtered(trial='all', scale='nonexistent')
        temp = stats.Helpers.r2pyjson('df_filtered')
        self.assertEqual(temp['delta_score'], [])

    @mock.patch.object(stats.StatsCore, 'get_stats')
    def test_get_trial_stats_SingleDf(self, mock_get_stats):
        """ Single df is 'concatanated' """

        mock_get_stats.side_effect = [
            StatsUnitTests.all_dfs,
            stats.NoSelectedRows,
        ]

        stats.Controllers.get_trial_stats(
            input_dir=folders.fixtures,
            input_fname='get_trial_data_stats_input1.csv',
            output_dir=folders.tmp_dir,
            output_prefix='get_trial_data_stats1',
            trial_scales={'tads': ['bdi', 'rads']},
        )

        output_model_summary = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_trial_data_stats1__model_summary.csv'))
        output_model_components = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_trial_data_stats1__model_components.csv'))
        output_strata_summary = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_trial_data_stats1__strata_summary.csv'))
        output_strata_contrast = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_trial_data_stats1__strata_contrast.csv'))

        # The +2s correspond to the GRC columns
        self.assertEqual(output_model_summary.shape, (1, 10))
        self.assertEqual(output_model_summary.f[0], 3)
        self.assertEqual(output_model_summary.adjr2[0], 4)
        self.assertEqual(output_model_summary.p[0], 5)

        self.assertEqual(output_model_components.shape, (1, 10))
        self.assertEqual(output_model_components.se[0], 2)
        self.assertEqual(output_model_components.t[0], 3)
        self.assertEqual(output_model_components.p[0], 4)

        self.assertEqual(output_strata_summary.shape, (1, 10))
        self.assertEqual(output_strata_summary.est[0], 4)
        self.assertEqual(output_strata_summary.se[0], 5)
        self.assertEqual(output_strata_summary.df[0], 6)

        self.assertEqual(output_strata_contrast.shape, (1, 12))
        self.assertEqual(output_strata_contrast.est[0], 10)
        self.assertEqual(output_strata_contrast.se[0], 20)
        self.assertEqual(output_strata_contrast.df[0], 30)
        self.assertEqual(output_strata_contrast.p_adj[0], 60)

    @mock.patch.object(stats.StatsCore, 'get_stats')
    def test_get_trial_stats_MultiDf(self, mock_get_stats):

        mock_get_stats.side_effect = [
            StatsUnitTests.all_dfs,
            StatsUnitTests.all_dfs2,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
        ]

        stats.Controllers.get_trial_stats(
            input_dir=folders.fixtures,
            input_fname='get_trial_data_stats_input1.csv',
            output_dir=folders.tmp_dir,
            output_prefix='get_trial_data_stats2',
            trial_scales={'tads': ['bdi', 'rads']},
        )

        output_model_summary = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_trial_data_stats2__model_summary.csv'))
        output_model_components = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_trial_data_stats2__model_components.csv'))
        output_strata_summary = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_trial_data_stats2__strata_summary.csv'))
        output_strata_contrast = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_trial_data_stats2__strata_contrast.csv'))

        # The +2s correspond to the GRC columns
        self.assertEqual(output_model_summary.shape, (2, 10))
        self.assertEqual(output_model_summary.f[0], 3)
        self.assertEqual(output_model_summary.adjr2[0], 4)
        self.assertEqual(output_model_summary.p[0], 5)
        self.assertEqual(output_model_summary.f[1], 7)
        self.assertEqual(output_model_summary.adjr2[1], 8)
        self.assertEqual(output_model_summary.p[1], 9)

        self.assertEqual(output_model_components.shape, (2, 10))
        self.assertEqual(output_model_components.se[0], 2)
        self.assertEqual(output_model_components.t[0], 3)
        self.assertEqual(output_model_components.p[0], 4)
        self.assertEqual(output_model_components.se[1], 5)
        self.assertEqual(output_model_components.t[1], 6)
        self.assertEqual(output_model_components.p[1], 7)

        self.assertEqual(output_strata_summary.shape, (2, 10))
        self.assertEqual(output_strata_summary.est[0], 4)
        self.assertEqual(output_strata_summary.se[0], 5)
        self.assertEqual(output_strata_summary.df[0], 6)
        self.assertEqual(output_strata_summary.est[1], 1)
        self.assertEqual(output_strata_summary.se[1], 2)
        self.assertEqual(output_strata_summary.df[1], 3)

        self.assertEqual(output_strata_contrast.shape, (2, 12))
        self.assertEqual(output_strata_contrast.est[0], 10)
        self.assertEqual(output_strata_contrast.se[0], 20)
        self.assertEqual(output_strata_contrast.df[0], 30)
        self.assertEqual(output_strata_contrast.p_adj[0], 60)
        self.assertEqual(output_strata_contrast.est[1], 1)
        self.assertEqual(output_strata_contrast.se[1], 2)
        self.assertEqual(output_strata_contrast.df[1], 3)
        self.assertEqual(output_strata_contrast.p_adj[1], 6)

    @mock.patch.object(stats.StatsCore, 'get_stats')
    def test_get_CGRC_stats_SingleDf(self, mock_get_stats):

        mock_get_stats.side_effect = [
            StatsUnitTests.all_dfs,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
        ]

        stats.Controllers.get_cgrc_stats(
            input_dir=folders.fixtures,
            input_fname='get_CGRC_stats1__cgrc.csv',
            output_dir=folders.tmp_dir,
            output_prefix='get_CGRC_stats1',
            trial_scales={'tads': ['bdi', 'rads']},
        )

        output_model_summary = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_CGRC_stats1__cgrc_model_summary.csv'))
        output_model_components = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_CGRC_stats1__cgrc_model_components.csv'))
        output_strata_summary = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_CGRC_stats1__cgrc_strata_summary.csv'))
        output_strata_contrast = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_CGRC_stats1__cgrc_strata_contrast.csv'))

        self.assertEqual(output_model_summary.shape, (1, 10))
        self.assertEqual(output_model_summary.f[0], 3)
        self.assertEqual(output_model_summary.adjr2[0], 4)
        self.assertEqual(output_model_summary.p[0], 5)

        self.assertEqual(output_model_components.shape, (1, 10))
        self.assertEqual(output_model_components.se[0], 2)
        self.assertEqual(output_model_components.t[0], 3)
        self.assertEqual(output_model_components.p[0], 4)

        self.assertEqual(output_strata_summary.shape, (1, 10))
        self.assertEqual(output_strata_summary.est[0], 4)
        self.assertEqual(output_strata_summary.se[0], 5)
        self.assertEqual(output_strata_summary.df[0], 6)

        self.assertEqual(output_strata_contrast.shape, (1, 12))
        self.assertEqual(output_strata_contrast.est[0], 10)
        self.assertEqual(output_strata_contrast.se[0], 20)
        self.assertEqual(output_strata_contrast.df[0], 30)
        self.assertEqual(output_strata_contrast.p_adj[0], 60)

    @mock.patch.object(stats.StatsCore, 'get_stats')
    def test_get_CGRC_stats_MultiDf(self, mock_get_stats):

        mock_get_stats.side_effect = [
            StatsUnitTests.all_dfs,
            StatsUnitTests.all_dfs2,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
            stats.NoSelectedRows,
        ]

        stats.Controllers.get_cgrc_stats(
            input_dir=folders.fixtures,
            input_fname='get_CGRC_stats2__cgrc.csv',
            output_dir=folders.tmp_dir,
            output_prefix='get_CGRC_stats2',
            trial_scales={'tads': ['bdi', 'rads']},
        )

        output_model_summary = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_CGRC_stats2__cgrc_model_summary.csv'))
        output_model_components = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_CGRC_stats2__cgrc_model_components.csv'))
        output_strata_summary = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_CGRC_stats2__cgrc_strata_summary.csv'))
        output_strata_contrast = pd.read_csv(os.path.join(
            folders.tmp_dir, 'get_CGRC_stats2__cgrc_strata_contrast.csv'))

        self.assertEqual(output_model_summary.shape, (2, 10))
        self.assertEqual(output_model_summary.f[0], 3)
        self.assertEqual(output_model_summary.adjr2[0], 4)
        self.assertEqual(output_model_summary.p[0], 5)
        self.assertEqual(output_model_summary.f[1], 7)
        self.assertEqual(output_model_summary.adjr2[1], 8)
        self.assertEqual(output_model_summary.p[1], 9)

        self.assertEqual(output_model_components.shape, (2, 10))
        self.assertEqual(output_model_components.se[0], 2)
        self.assertEqual(output_model_components.t[0], 3)
        self.assertEqual(output_model_components.p[0], 4)
        self.assertEqual(output_model_components.se[1], 5)
        self.assertEqual(output_model_components.t[1], 6)
        self.assertEqual(output_model_components.p[1], 7)

        self.assertEqual(output_strata_summary.shape, (2, 10))
        self.assertEqual(output_strata_summary.est[0], 4)
        self.assertEqual(output_strata_summary.se[0], 5)
        self.assertEqual(output_strata_summary.df[0], 6)
        self.assertEqual(output_strata_summary.est[1], 1)
        self.assertEqual(output_strata_summary.se[1], 2)
        self.assertEqual(output_strata_summary.df[1], 3)

        self.assertEqual(output_strata_contrast.shape, (2, 12))
        self.assertEqual(output_strata_contrast.est[0], 10)
        self.assertEqual(output_strata_contrast.se[0], 20)
        self.assertEqual(output_strata_contrast.df[0], 30)
        self.assertEqual(output_strata_contrast.p_adj[0], 60)
        self.assertEqual(output_strata_contrast.est[1], 1)
        self.assertEqual(output_strata_contrast.se[1], 2)
        self.assertEqual(output_strata_contrast.df[1], 3)
        self.assertEqual(output_strata_contrast.p_adj[1], 6)


class StatsIntegrationTests(unittest.TestCase):

    def test_get_model_stats1(self):

        r('df_filtered=read.csv("'+folders.fixtures.replace('\\', '/') +
          '//get_stats_data_input1.csv")')
        model_summary, model_components = stats.StatsCore.get_model_stats(add_cgrc_columns=False)

        # Reference outputs are manually checked against R output (codebase/tests/stats_calc_reference.r)
        # If new pseudodata is generated, check references again
        ref_components = pd.read_csv(os.path.join(
            folders.fixtures, 'get_stats_data_output1_model_components.csv'))
        ref_summary = pd.read_csv(os.path.join(
            folders.fixtures, 'get_stats_data_output1_model_summary.csv'))

        # Drop empty columns
        model_components.drop(['trial', 'scale'], axis=1, inplace=True)
        model_summary.drop(['trial', 'scale'], axis=1, inplace=True)
        ref_components.drop(['trial', 'scale', 'guesser', 'respondent',
                            'cgr', 'cgr_sim_id'], axis=1, inplace=True)
        ref_summary.drop(['trial', 'scale', 'guesser', 'respondent',
                         'cgr', 'cgr_sim_id'], axis=1, inplace=True)

        pd.testing.assert_frame_equal(ref_components, model_components)
        pd.testing.assert_frame_equal(ref_summary, model_summary)

    def test_get_strata_stats1(self):

        r('df_filtered=read.csv("'+folders.fixtures.replace('\\', '/') + '//get_stats_data_input1.csv")')
        strata_summary, strata_contrast = stats.StatsCore.get_strata_stats(add_cgrc_columns=False)

        # Reference outputs are manually checked against R output (codebase/tests/stats_calc_reference.r)
        # If new pseudodata is generated, check references again
        ref_summary = pd.read_csv(os.path.join(folders.fixtures, 'get_stats_data_output1_strata_summary.csv'))
        ref_contrast = pd.read_csv(os.path.join(folders.fixtures, 'get_stats_data_output1_strata_contrasts.csv'))

        # Drop empty columns
        strata_summary.drop(['trial', 'scale'], axis=1, inplace=True)
        strata_contrast.drop(['trial', 'scale'], axis=1, inplace=True)
        ref_summary.drop(['trial', 'scale', 'cgr', 'cgr_sim_id'], axis=1, inplace=True)
        ref_contrast.drop(['trial', 'scale', 'cgr', 'cgr_sim_id'], axis=1, inplace=True)

        strata_summary.df = strata_summary.df.astype('int64')
        strata_contrast.df = strata_contrast.df.astype('int64')

        pd.testing.assert_frame_equal(ref_contrast, strata_contrast)
        pd.testing.assert_frame_equal(ref_summary, strata_summary)
