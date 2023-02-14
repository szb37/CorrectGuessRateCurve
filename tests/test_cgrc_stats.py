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


class StatsIntegrationTests(unittest.TestCase):

    def test_get_model_stats1(self):

        r('df_filtered=read.csv("'+folders.fixtures.replace('\\', '/') +
          '//get_stats_data_input1.csv")')
        model_summary, model_comps = stats.StatsCore.get_model_stats()

        # Reference outputs are manually checked against R output (codebase/tests/stats_calc_reference.r)
        # If new pseudodata is generated, check references again
        ref_components = pd.read_csv(os.path.join(
            folders.fixtures, 'get_stats_data_output1_model_components.csv'))
        ref_summary = pd.read_csv(os.path.join(
            folders.fixtures, 'get_stats_data_output1_model_summary.csv'))

        # Drop empty columns
        model_comps.drop(['trial', 'scale'], axis=1, inplace=True)
        model_summary.drop(['trial', 'scale'], axis=1, inplace=True)
        ref_components.drop(['trial', 'scale', 'guesser', 'respondent',
                            'cgr', 'cgr_sim_id'], axis=1, inplace=True)
        ref_summary.drop(['trial', 'scale', 'guesser', 'respondent',
                         'cgr', 'cgr_sim_id'], axis=1, inplace=True)

        # convert frame type
        model_comps.__class__=pd.core.frame.DataFrame
        model_summary.__class__=pd.core.frame.DataFrame

        pd.testing.assert_frame_equal(ref_components, model_comps)
        pd.testing.assert_frame_equal(ref_summary, model_summary)

    def test_get_strata_stats1(self):

        r('df_filtered=read.csv("'+folders.fixtures.replace('\\', '/') + '//get_stats_data_input1.csv")')
        strata_summary, strata_contrast = stats.StatsCore.get_strata_stats()

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

        # convert frame type
        strata_contrast.__class__=pd.core.frame.DataFrame
        strata_summary.__class__=pd.core.frame.DataFrame

        pd.testing.assert_frame_equal(ref_contrast, strata_contrast)
        pd.testing.assert_frame_equal(ref_summary, strata_summary)
