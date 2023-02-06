"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2021, DrugNerdsLab
:License: MIT
"""

import src.toy_models.core as toy_models
import src.my_dataframes as mydfs
import src.folders as folders
from statistics import mean
import unittest


class IntegrationTests(unittest.TestCase):

    def test_get_pseudodata_case1(self):
        toy_models.Controllers.run_cgrc_model_family(
            model_family_name='default_models',
            n_trials=1,
            n_patients=50,
            postfix='test',
            cgrc_param_set=1,
            save_figs=False,
        )


class DataGeneratorUnitTests(unittest.TestCase):
    """ Tests are probabilistic, ~1% chance of failure; thrs are set based on num experiments """

    def test_get_pseudodata_case1(self):

        score_plpl = []
        score_acpl = []
        score_plac = []
        score_acac = []
        n_plpl = []
        n_acpl = []
        n_plac = []
        n_acac = []

        model_specs = {
            'oc_nh': (20, 2),
            'gs_nh': (0.5, 0.25),
            'se': (0.15, 0.05),
            'dte': (2, 1),
            'pte': (0, 0),
            'ate': (4, 1),
            'oc2gs': (0, 0),
        }
        analysis_name = 'temp'
        n_datapoints = 256

        for idx in range(32):

            df = toy_models.ToyModelsDataGenerator.get_pseudodata(
                output_dir=folders.tmp_dir,
                output_prefix='test_get_pseudodata_case1',
                model_specs=model_specs,
                n_datapoints=n_datapoints,
                model_name='test',
                trial_id=idx,
                min_strata_size=4,
                round_digits=5,
            )

            n_plpl.append(df.loc[(df.condition == 'PL')
                          & (df.guess == 'PL')].shape[0])
            n_acpl.append(df.loc[(df.condition == 'AC')
                          & (df.guess == 'PL')].shape[0])
            n_plac.append(df.loc[(df.condition == 'PL')
                          & (df.guess == 'AC')].shape[0])
            n_acac.append(df.loc[(df.condition == 'AC')
                          & (df.guess == 'AC')].shape[0])

            score_plpl.append(
                mean(df.loc[(df.condition == 'PL') & (df.guess == 'PL'), 'score']))
            score_acpl.append(
                mean(df.loc[(df.condition == 'AC') & (df.guess == 'PL'), 'score']))
            score_plac.append(
                mean(df.loc[(df.condition == 'PL') & (df.guess == 'AC'), 'score']))
            score_acac.append(
                mean(df.loc[(df.condition == 'AC') & (df.guess == 'AC'), 'score']))

        self.assertTrue(abs(n_datapoints/4 - mean(n_plpl)) <= 14)
        self.assertTrue(abs(n_datapoints/2*0.2782 - mean(n_acpl)) <= 11)
        self.assertTrue(abs(n_datapoints/4 - mean(n_plac)) <= 14)
        self.assertTrue(abs(n_datapoints/2*0.7218 - mean(n_acac)) <= 15)

        self.assertTrue(abs(20 - mean(score_plpl)) <= 0.5)
        self.assertTrue(abs(22 - mean(score_acpl)) <= 0.75)
        self.assertTrue(abs(24 - mean(score_plac)) <= 0.5)
        self.assertTrue(abs(26 - mean(score_acac)) <= 0.5)

    def test_get_pseudodata_case2(self):

        score_plpl = []
        score_acpl = []
        score_plac = []
        score_acac = []
        n_plpl = []
        n_acpl = []
        n_plac = []
        n_acac = []

        model_specs = {
            'oc_nh': (20, 2),
            'gs_nh': (0.5, 0.125),
            'se': (0.25, 0.08),
            'dte': (3, 1),
            'pte': (0, 0),
            'ate': (5, 1),
            'oc2gs': (0, 0),
        }
        analysis_name = 'temp'
        n_datapoints = 256

        for idx in range(32):

            df = toy_models.ToyModelsDataGenerator.get_pseudodata(
                output_dir=folders.tmp_dir,
                output_prefix='test_get_pseudodata_case2',
                model_specs=model_specs,
                n_datapoints=n_datapoints,
                model_name='test',
                trial_id=idx,
                min_strata_size=None,
                round_digits=5,
            )

            n_plpl.append(df.loc[(df.condition == 'PL')
                          & (df.guess == 'PL')].shape[0])
            n_acpl.append(df.loc[(df.condition == 'AC')
                          & (df.guess == 'PL')].shape[0])
            n_plac.append(df.loc[(df.condition == 'PL')
                          & (df.guess == 'AC')].shape[0])
            n_acac.append(df.loc[(df.condition == 'AC')
                          & (df.guess == 'AC')].shape[0])

            score_plpl.append(
                mean(df.loc[(df.condition == 'PL') & (df.guess == 'PL'), 'score']))
            score_acpl.append(
                mean(df.loc[(df.condition == 'AC') & (df.guess == 'PL'), 'score']))
            score_plac.append(
                mean(df.loc[(df.condition == 'PL') & (df.guess == 'AC'), 'score']))
            score_acac.append(
                mean(df.loc[(df.condition == 'AC') & (df.guess == 'AC'), 'score']))

        self.assertTrue(abs(n_datapoints/4 - mean(n_plpl)) <= 14)
        self.assertTrue(abs(n_datapoints/2*0.046 - mean(n_acpl)) <= 4)
        self.assertTrue(abs(n_datapoints/4 - mean(n_plac)) <= 14)
        self.assertTrue(abs(n_datapoints/2*0.954 - mean(n_acac)) <= 16)

        self.assertTrue(abs(20 - mean(score_plpl)) <= 0.5)
        self.assertTrue(abs(23 - mean(score_acpl)) <= 2.1)
        self.assertTrue(abs(25 - mean(score_plac)) <= 0.6)
        self.assertTrue(abs(28 - mean(score_acac)) <= 0.5)

    def test_min_strata_size(self):

        model_specs = {
            'oc_nh': (20, 3),
            'gs_nh': (0, 0),
            'se': (0, 0),
            'dte': (3, 1),
            'pte': (0, 0),
            'ate': (0, 0),
            'oc2gs': (0, 0),
        }
        analysis_name = 'temp'
        n_datapoints = 64

        df = toy_models.ToyModelsDataGenerator.get_pseudodata(
            output_dir=folders.tmp_dir,
            output_prefix='test_min_strata_size',
            model_specs=model_specs,
            n_datapoints=n_datapoints,
            model_name='test',
            trial_id=1,
            min_strata_size=None,
            round_digits=5,
        )
        self.assertEqual(df.loc[(df.guess == 'PL')].shape[0], 64)
        self.assertEqual(df.loc[(df.guess == 'AC')].shape[0], 0)

        df = toy_models.ToyModelsDataGenerator.get_pseudodata(
            output_dir=folders.tmp_dir,
            output_prefix='test_min_strata_size',
            model_specs=model_specs,
            n_datapoints=n_datapoints,
            model_name='test',
            trial_id=1,
            min_strata_size=5,
            round_digits=5,
        )
        self.assertEqual(df.loc[(df.condition == 'PL')
                         & (df.guess == 'AC')].shape[0], 5)
        self.assertEqual(df.loc[(df.condition == 'AC')
                         & (df.guess == 'AC')].shape[0], 5)
        self.assertEqual(df.loc[(df.guess == 'PL')].shape[0], 54)

        df = toy_models.ToyModelsDataGenerator.get_pseudodata(
            output_dir=folders.tmp_dir,
            output_prefix='test_min_strata_size',
            model_specs=model_specs,
            n_datapoints=n_datapoints,
            model_name='test',
            trial_id=1,
            min_strata_size=15,
            round_digits=5,
        )
        self.assertEqual(df.loc[(df.condition == 'PL')
                         & (df.guess == 'AC')].shape[0], 15)
        self.assertEqual(df.loc[(df.condition == 'AC')
                         & (df.guess == 'AC')].shape[0], 15)
        self.assertEqual(df.loc[(df.guess == 'PL')].shape[0], 34)


class HelpersUnitTests(unittest.TestCase):
    """ Tests are probabilistic, ~0.003% chance of failure """

    def test_get_TrialDatadDf(self):

        n = 50
        df = toy_models.Helpers.get_TrialDatadDf(n)
        self.assertEqual(df.shape, (50, 11))

        n = 30
        df = toy_models.Helpers.get_TrialDatadDf(n)
        self.assertEqual(df.shape, (30, 11))
