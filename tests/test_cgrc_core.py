"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import src.cgrc.core as cgrc
import src.cgrc.stats as stats
import src.dataframe_classes as df_class
import src.config as config
import src.folders as folders
from numpy import linspace as linspace
from rpy2.robjects import r
from unittest import mock
import pandas as pd
import unittest
import os


class CorrectGuessRateCurveUnitTests(unittest.TestCase):

    def test_get_strata_sample_sizes_method_all_equal(self):

        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.5,
            strata_sampling='all_equal')
        self.assertEqual(sample_sizes['PLPL'], 25)
        self.assertEqual(sample_sizes['ACAC'], 25)
        self.assertEqual(sample_sizes['ACPL'], 25)
        self.assertEqual(sample_sizes['PLAC'], 25)

        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.8,
            strata_sampling='all_equal')
        self.assertEqual(sample_sizes['PLPL'], 40)
        self.assertEqual(sample_sizes['ACAC'], 40)
        self.assertEqual(sample_sizes['ACPL'], 10)
        self.assertEqual(sample_sizes['PLAC'], 10)

        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.25,
            strata_sampling='all_equal')
        self.assertEqual(sample_sizes['PLPL'], 12)
        self.assertEqual(sample_sizes['ACAC'], 13)
        self.assertEqual(sample_sizes['ACPL'], 38)
        self.assertEqual(sample_sizes['PLAC'], 37)

    def test_get_strata_sample_sizes_method_active_equal(self):

        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.5,
            strata_sampling='active_equal')
        self.assertEqual(sample_sizes['PLPL'], 25)
        self.assertEqual(sample_sizes['PLAC'], 25)
        self.assertEqual(sample_sizes['ACAC'], 25)
        self.assertEqual(sample_sizes['ACPL'], 25)

        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.8,
            strata_sampling='active_equal')
        self.assertEqual(sample_sizes['PLPL'], 25)
        self.assertEqual(sample_sizes['PLAC'], 25)
        self.assertEqual(sample_sizes['ACAC'], 40)
        self.assertEqual(sample_sizes['ACPL'], 10)

        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.35,
            strata_sampling='active_equal')
        self.assertEqual(sample_sizes['PLPL'], 25)
        self.assertEqual(sample_sizes['PLAC'], 25)
        self.assertEqual(sample_sizes['ACAC'], 18)
        self.assertEqual(sample_sizes['ACPL'], 32)

    @mock.patch.object(cgrc.CorrectGuessRateCurve, 'get_strata_ratio')
    def test_get_strata_sample_sizes_method_all_prop(self, mock_get_strata_ratio):

        mock_get_strata_ratio.return_value = {
            'PLPL': 0.25,
            'ACPL': 0.25,
            'PLAC': 0.25,
            'ACAC': 0.25,
        }
        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.5,
            strata_sampling='all_prop')
        self.assertEqual(sample_sizes['PLPL'], 25)
        self.assertEqual(sample_sizes['ACAC'], 25)
        self.assertEqual(sample_sizes['ACPL'], 25)
        self.assertEqual(sample_sizes['PLAC'], 25)

        mock_get_strata_ratio.return_value = {
            'PLPL': 0.2,
            'ACPL': 0.4,
            'PLAC': 0.1,
            'ACAC': 0.3,
        }
        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.5,
            strata_sampling='all_prop')
        self.assertEqual(sample_sizes['PLPL'], 20)
        self.assertEqual(sample_sizes['ACPL'], 40)
        self.assertEqual(sample_sizes['PLAC'], 10)
        self.assertEqual(sample_sizes['ACAC'], 30)

        mock_get_strata_ratio.return_value = {
            'PLPL': 0.3,
            'ACPL': 0.15,
            'PLAC': 0.15,
            'ACAC': 0.4,
        }
        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.7,
            strata_sampling='all_prop')
        self.assertEqual(sample_sizes['PLPL'], 30)
        self.assertEqual(sample_sizes['ACPL'], 15)
        self.assertEqual(sample_sizes['PLAC'], 15)
        self.assertEqual(sample_sizes['ACAC'], 40)

    @mock.patch.object(cgrc.CorrectGuessRateCurve, 'get_strata_ratio')
    def test_get_strata_sample_sizes_method_active_prop(self, mock_get_strata_ratio):

        mock_get_strata_ratio.return_value = {
            'PLPL': 0.25,
            'ACPL': 0.25,
            'PLAC': 0.25,
            'ACAC': 0.25,
        }
        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.5,
            strata_sampling='active_prop')
        self.assertEqual(sample_sizes['PLPL'], 25)
        self.assertEqual(sample_sizes['PLAC'], 25)
        self.assertEqual(sample_sizes['ACAC'], 25)
        self.assertEqual(sample_sizes['ACPL'], 25)

        mock_get_strata_ratio.return_value = {
            'PLPL': 0.2,
            'ACPL': 0.4,
            'PLAC': 0.1,
            'ACAC': 0.3,
        }
        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.5,
            strata_sampling='active_prop')
        self.assertEqual(sample_sizes['PLPL'], 20)
        self.assertEqual(sample_sizes['PLAC'], 10)
        self.assertEqual(sample_sizes['ACAC'], 35)
        self.assertEqual(sample_sizes['ACPL'], 35)

        mock_get_strata_ratio.return_value = {
            'PLPL': 0.3,
            'ACPL': 0.15,
            'PLAC': 0.15,
            'ACAC': 0.4,
        }
        sample_sizes = cgrc.CorrectGuessRateCurve.get_strata_sample_sizes(
            total_sample_size=100,
            df_filtered=None,
            correct_guess_rate=0.75,
            strata_sampling='active_prop')

        self.assertEqual(sample_sizes['PLPL'], 30)
        self.assertEqual(sample_sizes['PLAC'], 15)
        self.assertEqual(sample_sizes['ACAC'], 41)
        self.assertEqual(sample_sizes['ACPL'], 14)

    def test_get_strata_ratio(self):

        df_filtered = pd.read_csv(os.path.join(
            folders.fixtures, 'get_strata_ratio_input1.csv'))
        strata_ratio = cgrc.CorrectGuessRateCurve.get_strata_ratio(
            df_filtered=df_filtered)
        self.assertEqual(strata_ratio['PLPL'], 0.25)
        self.assertEqual(strata_ratio['ACPL'], 0.25)
        self.assertEqual(strata_ratio['PLAC'], 0.25)
        self.assertEqual(strata_ratio['ACAC'], 0.25)

        df_filtered = pd.read_csv(os.path.join(
            folders.fixtures, 'get_strata_ratio_input2.csv'))
        strata_ratio = cgrc.CorrectGuessRateCurve.get_strata_ratio(
            df_filtered=df_filtered)
        self.assertEqual(strata_ratio['PLPL'], round(1/20, 2))
        self.assertEqual(strata_ratio['ACPL'], round(6/20, 2))
        self.assertEqual(strata_ratio['PLAC'], round(3/20, 2))
        self.assertEqual(strata_ratio['ACAC'], round(10/20, 2))

    def test_get_cgrc_datapoint_KDE(self):

        df_filtered = pd.read_csv(os.path.join(
            folders.fixtures, 'get_kdes_input1.csv'))
        kdes = cgrc.CorrectGuessRateCurve.get_kdes(df_filtered=df_filtered)
        df = cgrc.CorrectGuessRateCurve.get_cgrc_datapoint_KDE(
            df_filtered=df_filtered,
            sample_sizes={'PLPL': 200, 'ACPL': 210, 'PLAC': 220, 'ACAC': 230},
            kdes=kdes,
            cgr=666)

        self.assertEqual(df.cgr.unique().tolist(), [666])

        self.assertEqual(df.loc[(df.condition == 'PL')
                         & (df.guess == 'PL')].shape[0], 200)
        self.assertEqual(df.loc[(df.condition == 'AC')
                         & (df.guess == 'PL')].shape[0], 210)
        self.assertEqual(df.loc[(df.condition == 'PL')
                         & (df.guess == 'AC')].shape[0], 220)
        self.assertEqual(df.loc[(df.condition == 'AC')
                         & (df.guess == 'AC')].shape[0], 230)

        self.assertTrue(
            abs(df.loc[(df.condition == 'PL') & (df.guess == 'PL')].delta_score.mean()-10) < 0.75)
        self.assertTrue(
            abs(df.loc[(df.condition == 'AC') & (df.guess == 'PL')].delta_score.mean()-6) < 0.75)
        self.assertTrue(
            abs(df.loc[(df.condition == 'PL') & (df.guess == 'AC')].delta_score.mean()-100) < 0.75)
        self.assertTrue(
            abs(df.loc[(df.condition == 'AC') & (df.guess == 'AC')].delta_score.mean()-20) < 0.75)


class CorrectGuessRateCurveIntegrationTests(unittest.TestCase):

    def test_get_cgrc1_prop(self):

        analysis_name = 'get_cgrc1'
        cgrc.CorrectGuessRateCurve.get_cgrc_data(
            input_dir=folders.fixtures,
            input_fname='get_cgrc_input1__trial_data.csv',
            output_dir=folders.tmp_dir,
            output_prefix=analysis_name,
            trial_scales={'test': ['tadaa']},
            cgrc_parameters={
                'cgr_values': linspace(0, 1, 5).tolist(),
                'n_cgrc_trials': 1, },
            strata_sampling='all_prop',
            n_sample=300
        )

        df = pd.read_csv(os.path.join(
            folders.tmp_dir, analysis_name+'__cgrc_data.csv'))

        temp = df.loc[df.cgr == 0]  # ~60
        self.assertTrue(abs(temp.delta_score.mean()-60) < 0.65)

        temp = df.loc[df.cgr == 0.25]  # ~50
        self.assertTrue(abs(temp.delta_score.mean()-50) < 0.65)

        temp = df.loc[df.cgr == 0.5]  # ~40
        self.assertTrue(abs(temp.delta_score.mean()-40) < 0.65)

        temp = df.loc[df.cgr == 0.75]  # ~30
        self.assertTrue(abs(temp.delta_score.mean()-30) < 0.65)

        temp = df.loc[df.cgr == 1]  # ~20
        self.assertTrue(abs(temp.delta_score.mean()-20) < 0.65)

    def test_get_cgrc2_prop(self):

        analysis_name = 'get_cgrc2'
        cgrc.CorrectGuessRateCurve.get_cgrc_data(
            input_dir=folders.fixtures,
            input_fname='get_cgrc_input2__trial_data.csv',
            output_dir=folders.tmp_dir,
            output_prefix=analysis_name,
            trial_scales={'test': ['tadaa']},
            cgrc_parameters={
                'cgr_values': linspace(0, 1, 5).tolist(),
                'n_cgrc_trials': 1, },
            strata_sampling='all_prop',
            n_sample=800
        )

        df = pd.read_csv(os.path.join(
            folders.tmp_dir, analysis_name+'__cgrc_data.csv'))

        temp = df.loc[df.cgr == 0]    # ((3*50)+(5*70))/8=62.5
        self.assertTrue(abs(temp.delta_score.mean()-62.5) < 0.65)

        # (1/4)*(((3*10)+(5*30))/8) + (3/4)*(((3*50)+(5*70))/8)=52.5
        temp = df.loc[df.cgr == 0.25]
        self.assertTrue(abs(temp.delta_score.mean()-52.5) < 0.65)

        # (1/2)*(((3*10)+(5*30))/8) + (1/2)*(((3*50)+(5*70))/8)=42.5
        temp = df.loc[df.cgr == 0.5]
        self.assertTrue(abs(temp.delta_score.mean()-42.5) < 0.65)

        # (3/4)*(((3*10)+(5*30))/8) + (1/4)*(((3*50)+(5*70))/8)=32.5
        temp = df.loc[df.cgr == 0.75]
        self.assertTrue(abs(temp.delta_score.mean()-32.5) < 0.65)

        temp = df.loc[df.cgr == 1]    # ((3*10)+(5*30))/8=22.5
        self.assertTrue(abs(temp.delta_score.mean()-22.5) < 0.65)

    def test_get_cgrc2_equal(self):

        analysis_name = 'get_cgrc3'
        cgrc.CorrectGuessRateCurve.get_cgrc_data(
            input_dir=folders.fixtures,
            input_fname='get_cgrc_input2__trial_data.csv',
            output_dir=folders.tmp_dir,
            output_prefix=analysis_name,
            trial_scales={'test': ['tadaa']},
            cgrc_parameters={
                'cgr_values': linspace(0, 1, 5).tolist(),
                'n_cgrc_trials': 1, },
            strata_sampling='all_equal',
            n_sample=800,
        )

        df = pd.read_csv(os.path.join(
            folders.tmp_dir, analysis_name+'__cgrc_data.csv'))

        temp = df.loc[df.cgr == 0]    # ((50)+(70))/2=60
        self.assertTrue(abs(temp.delta_score.mean()-60) < 0.65)

        # (1/4)*(((10)+(30))/2) + (3/4)*(((50)+(70))/2)=50
        temp = df.loc[df.cgr == 0.25]
        self.assertTrue(abs(temp.delta_score.mean()-50) < 0.65)

        # (1/2)*(((10)+(30))/2) + (1/2)*(((50)+(70))/2)=40
        temp = df.loc[df.cgr == 0.5]
        self.assertTrue(abs(temp.delta_score.mean()-40) < 0.65)

        # (3/4)*(((10)+(30))/2) + (1/4)*(((50)+(70))/2)=30
        temp = df.loc[df.cgr == 0.75]
        self.assertTrue(abs(temp.delta_score.mean()-30) < 0.65)

        temp = df.loc[df.cgr == 1]    # ((10)+(30))/2=20
        self.assertTrue(abs(temp.delta_score.mean()-20) < 0.65)

    def test_get_cgrc_sample_size(self):

        analysis_name = 'get_cgrc4'

        cgrc.CorrectGuessRateCurve.get_cgrc_data(
            input_dir=folders.fixtures,
            input_fname='get_cgrc_input4__trial_data.csv',
            output_dir=folders.tmp_dir,
            output_prefix=analysis_name,
            trial_scales={'test': ['tadaa']},
            cgrc_parameters={
                'cgr_values': [0.25, 0.75],
                'n_cgrc_trials': 1, },
        )
        df = pd.read_csv(os.path.join(
            folders.tmp_dir, analysis_name+'__cgrc_data.csv'))
        self.assertEqual(df.shape[0], 2*1*119)

        cgrc.CorrectGuessRateCurve.get_cgrc_data(
            input_dir=folders.fixtures,
            input_fname='get_cgrc_input4__trial_data.csv',
            output_dir=folders.tmp_dir,
            output_prefix=analysis_name,
            trial_scales={'test': ['tadaa']},
            cgrc_parameters={
                'cgr_values': [0.25, 0.75],
                'n_cgrc_trials': 1, },
            n_sample=100,
        )
        df = pd.read_csv(os.path.join(
            folders.tmp_dir, analysis_name+'__cgrc_data.csv'))
        self.assertEqual(df.shape[0], 2*100)

        cgrc.CorrectGuessRateCurve.get_cgrc_data(
            input_dir=folders.fixtures,
            input_fname='get_cgrc_input4__trial_data.csv',
            output_dir=folders.tmp_dir,
            output_prefix=analysis_name,
            trial_scales={'test': ['tadaa']},
            cgrc_parameters={
                'cgr_values': [0.25, 0.75],
                'n_cgrc_trials': 2, },
        )
        df = pd.read_csv(os.path.join(
            folders.tmp_dir, analysis_name+'__cgrc_data.csv'))
        self.assertEqual(df.shape[0], 2*2*119)

    def test_get_cgrc_stats(self):

        analysis_name = 'get_cgrc3'
        stats.Controllers.get_cgrc_stats(
            input_dir=folders.fixtures,
            input_fname='get_CGRC_stats3__cgrc.csv',
            output_dir=folders.tmp_dir,
            output_prefix=analysis_name,
            trial_scales={'test': ['tadaa']},
            do_stratas=True,
        )

        # Check model components
        df = pd.read_csv(os.path.join(
            folders.tmp_dir, analysis_name+'__cgrc_model_components.csv'))

        temp = df.loc[(df.cgr == 0) & (df.model_type == 'without_guess') & (
            df.component == 'intercept')]
        self.assertTrue(abs(temp.est.tolist()[0]-50) < 0.85)
        temp = df.loc[(df.cgr == 0) & (df.model_type == 'without_guess') & (
            df.component == 'conditionAC')]
        self.assertTrue(abs(temp.est.tolist()[0]-20) < 0.85)
        temp = df.loc[(df.cgr == 0) & (df.model_type == 'with_guess') & (
            df.component == 'intercept')]
        self.assertTrue(abs(temp.est.tolist()[0]-50) < 0.85)
        temp = df.loc[(df.cgr == 0) & (df.model_type == 'with_guess') & (
            df.component == 'conditionAC')]
        self.assertTrue(abs(temp.est.tolist()[0]-20) < 0.85)

        temp = df.loc[(df.cgr == 0.5) & (df.model_type ==
                                         'without_guess') & (df.component == 'intercept')]
        self.assertTrue(abs(temp.est.tolist()[0]-30) < 0.85)
        temp = df.loc[(df.cgr == 0.5) & (df.model_type == 'without_guess') & (
            df.component == 'conditionAC')]
        self.assertTrue(abs(temp.est.tolist()[0]-20) < 0.85)
        temp = df.loc[(df.cgr == 0.5) & (df.model_type ==
                                         'with_guess') & (df.component == 'intercept')]
        self.assertTrue(abs(temp.est.tolist()[0]-10) < 0.85)
        temp = df.loc[(df.cgr == 0.5) & (df.model_type == 'with_guess') & (
            df.component == 'conditionAC')]
        self.assertTrue(abs(temp.est.tolist()[0]-60) < 0.85)
        temp = df.loc[(df.cgr == 0.5) & (df.model_type ==
                                         'with_guess') & (df.component == 'guessAC')]
        self.assertTrue(abs(temp.est.tolist()[0]-40) < 0.85)
        temp = df.loc[(df.cgr == 0.5) & (df.model_type == 'with_guess') & (
            df.component == 'conditionAC:guessAC')]
        self.assertTrue(abs(abs(temp.est.tolist()[0])-80) < 0.85)

        temp = df.loc[(df.cgr == 1) & (df.model_type == 'without_guess') & (
            df.component == 'intercept')]
        self.assertTrue(abs(temp.est.tolist()[0]-10) < 0.85)
        temp = df.loc[(df.cgr == 1) & (df.model_type == 'without_guess') & (
            df.component == 'conditionAC')]
        self.assertTrue(abs(temp.est.tolist()[0]-20) < 0.85)
        temp = df.loc[(df.cgr == 1) & (df.model_type == 'with_guess') & (
            df.component == 'intercept')]
        self.assertTrue(abs(temp.est.tolist()[0]-10) < 0.85)
        temp = df.loc[(df.cgr == 1) & (df.model_type == 'with_guess') & (
            df.component == 'conditionAC')]
        self.assertTrue(abs(temp.est.tolist()[0]-20) < 0.85)

        # Check model summary
        df = pd.read_csv(os.path.join(
            folders.tmp_dir, analysis_name+'__cgrc_model_summary.csv'))

        self.assertTrue(all([df2 in [296, 298] for df2 in df.df2.tolist()]))
        self.assertTrue(df.shape == (10, 10))
        self.assertTrue(all([abs(p) < 0.01 for p in df.p.tolist()]))
        temp = df.loc[(df.cgr.isin([0.25, 0.5, 0.75])) &
                      (df.model_type == 'with_guess')]
        self.assertTrue(all([abs(r2-1) < 0.05 for r2 in temp.adjr2.tolist()]))
        self.assertTrue(all([df1 == 3 for df1 in temp.df1.tolist()]))

        # Check strata summary
        df = pd.read_csv(os.path.join(
            folders.tmp_dir, analysis_name+'__cgrc_strata_summary.csv'))
        self.assertTrue(df.shape == (12, 10))
        temp = df.loc[(df.strata == 'PLPL')]
        self.assertTrue(all([abs(est-10) < 1 for est in temp.est.tolist()]))
        temp = df.loc[(df.strata == 'ACPL')]
        self.assertTrue(all([abs(est-70) < 1 for est in temp.est.tolist()]))
        temp = df.loc[(df.strata == 'PLAC')]
        self.assertTrue(all([abs(est-50) < 1 for est in temp.est.tolist()]))
        temp = df.loc[(df.strata == 'ACAC')]
        self.assertTrue(all([abs(est-30) < 1 for est in temp.est.tolist()]))

        # Check strata contrast
        df = pd.read_csv(os.path.join(
            folders.tmp_dir, analysis_name+'__cgrc_strata_contrast.csv'))

        self.assertTrue(df.shape == (12, 12))
        self.assertTrue(all([abs(p) < 0.001 for p in df.p.tolist()]))
        self.assertTrue(all([df == 296 for df in df.df.tolist()]))

        self.assertTrue(all(
            [abs(abs(est)-40) < 1 for est in df.loc[df.contrast == 'PLPLvsPLAC'].est.tolist()]))
        self.assertTrue(all(
            [abs(abs(est)-40) < 1 for est in df.loc[df.contrast == 'ACPLvsACAC'].est.tolist()]))
        self.assertTrue(all(
            [abs(abs(est)-60) < 1 for est in df.loc[df.contrast == 'PLPLvsACPL'].est.tolist()]))
        self.assertTrue(all(
            [abs(abs(est)-20) < 1 for est in df.loc[df.contrast == 'PLACvsACAC'].est.tolist()]))

    def test_main_analyss(self):

        cgrc.Controllers.run_cgrc_trial(
            trial_name='sbmd',
            postfix='tmp',
            cgrc_param_set='test',
            trial_scales=config.sbmd_test,
            save_figs=False,
            do_stratas=True,
        )
