"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import src.toy_models.core as toy_models
import src.toy_models.model_defs as model_defs
import src.config as config
import src.folders as folders
from statistics import mean
from math import sqrt
import unittest
import mock


class ToyModelsDataGeneratorUnitTests(unittest.TestCase):

    def test_enforce_all_strata(self):

        conditions = ['AC', 'AC', 'AC', 'AC']
        guesses    = ['AC', 'AC', 'AC', 'AC']
        conditions, guesses = toy_models.ToyModelsDataGenerator.enforce_all_strata(conditions, guesses)
        assert conditions == ['PL', 'AC', 'PL', 'AC']
        assert guesses    == ['PL', 'PL', 'AC', 'AC']

        conditions = ['PL', 'AC', 'PL', 'PL', 'whatever', 0]
        guesses    = ['AC', 'PL', 'PL', 'AC', 0, 'whatever']
        conditions, guesses = toy_models.ToyModelsDataGenerator.enforce_all_strata(conditions, guesses)
        assert conditions == ['PL', 'AC', 'PL', 'AC', 'whatever', 0]
        assert guesses    == ['PL', 'PL', 'AC', 'AC', 0, 'whatever']

        with self.assertRaises(AssertionError):
            toy_models.ToyModelsDataGenerator.enforce_all_strata([0, 'whatever'], ['AC', 'AC', 'AC', 'AC'])

        with self.assertRaises(AssertionError):
            toy_models.ToyModelsDataGenerator.enforce_all_strata(['AC', 'AC', 'AC', 'AC'], [0, 'whatever'])

    def test_get_toymodel_data_extremes(self):

        model = model_defs.ModelDefinition(n_patients=100, p_act=0, p_unb=0)
        df = toy_models.ToyModelsDataGenerator.get_toymodel_data(
            output_dir=folders.tmp_dir,
            output_prefix='test_get_toymodel_data_extremes',
            model=model,
            model_sim_id=0,
            enforce_all_strata=False
        )
        self.assertEqual(df.condition.to_list().count('AC'), 0)
        self.assertEqual(df.condition.to_list().count('PL'), 100)
        self.assertEqual(df.guess.to_list().count('AC'), 100)
        self.assertEqual(df.guess.to_list().count('PL'), 0)


        model = model_defs.ModelDefinition(n_patients=100, p_act=0, p_unb=1)
        df = toy_models.ToyModelsDataGenerator.get_toymodel_data(
            output_dir=folders.tmp_dir,
            output_prefix='test_get_toymodel_data_extremes',
            model=model,
            model_sim_id=0,
            enforce_all_strata=False
        )
        self.assertEqual(df.condition.to_list().count('AC'), 0)
        self.assertEqual(df.condition.to_list().count('PL'), 100)
        self.assertEqual(df.guess.to_list().count('AC'), 0)
        self.assertEqual(df.guess.to_list().count('PL'), 100)

        model = model_defs.ModelDefinition(n_patients=100, p_act=1, p_unb=0)
        df = toy_models.ToyModelsDataGenerator.get_toymodel_data(
            output_dir=folders.tmp_dir,
            output_prefix='test_get_toymodel_data_extremes',
            model=model,
            model_sim_id=0,
            enforce_all_strata=False
        )
        self.assertEqual(df.condition.to_list().count('AC'), 100)
        self.assertEqual(df.condition.to_list().count('PL'), 0)
        self.assertEqual(df.guess.to_list().count('AC'), 0)
        self.assertEqual(df.guess.to_list().count('PL'), 100)

        model = model_defs.ModelDefinition(n_patients=100, p_act=1, p_unb=1)
        df = toy_models.ToyModelsDataGenerator.get_toymodel_data(
            output_dir=folders.tmp_dir,
            output_prefix='test_get_toymodel_data_extremes',
            model=model,
            model_sim_id=0,
            enforce_all_strata=False
        )
        self.assertEqual(df.condition.to_list().count('AC'), 100)
        self.assertEqual(df.condition.to_list().count('PL'), 0)
        self.assertEqual(df.guess.to_list().count('AC'), 100)
        self.assertEqual(df.guess.to_list().count('PL'), 0)

    # Test is probabilistic; try running again if fail
    def test_get_toymodel_data_normals(self, n=1000, thr=2.5):

        model_params = [
            {'n_patients':n, 'p_act':0.5, 'p_unb':0.5},  # mid / mid
            {'n_patients':n, 'p_act':0.2, 'p_unb':0.5},  # low / mid
            {'n_patients':n, 'p_act':0.5, 'p_unb':0.2},  # mid / low
            {'n_patients':n, 'p_act':0.85, 'p_unb':0.5}, # high / mid
            {'n_patients':n, 'p_act':0.5, 'p_unb':0.85}, # mid / high
        ]

        for model_param in model_params:

            n_patients = model_param['n_patients']
            p_act = model_param['p_act']
            p_unb = model_param['p_unb']

            df = toy_models.ToyModelsDataGenerator.get_toymodel_data(
                output_dir=folders.tmp_dir,
                output_prefix='test_get_toymodel_data_normals',
                model=model_defs.ModelDefinition(
                    n_patients=n_patients,
                    p_act=p_act,
                    p_unb=p_unb,),
                model_sim_id=0,
                enforce_all_strata=False
            )

            se = sqrt((1-p_act)*(p_act)/n)
            self.assertTrue(abs(p_act - df.condition.to_list().count('AC')/n) <= thr*se)
            self.assertTrue(abs((1-p_act) - df.condition.to_list().count('PL')/n) <= thr*se)

            se = sqrt((1-p_unb)*(p_unb)/n)
            n_unblind = sum(
                [guess==condition for condition, guess in zip(df.condition.to_list(), df.guess.to_list())])
            self.assertTrue(abs(p_unb - n_unblind/n) <= thr*se)


class IntegrationTests(unittest.TestCase):

    def test_get_toymodel_data_case1(self):
        toy_models.Controllers.run_toymodels_cgrc(
            analysis_name = 'test',
            models = model_defs.test,
            cgrc_param_set = 'cgrA_low',
            n_trials = 2
        )
