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
    """ Tests are probabilistic, ~1% chance of failure; thrs are set based on num experiments """

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

        model = model_defs.ModelDefinition(p_act=0, p_sep=0)
        df = toy_models.ToyModelsDataGenerator.get_toymodel_data(
            output_dir=folders.tmp_dir,
            output_prefix='test_get_toymodel_data_extremes',
            model=model,
            n=100,
            model_sim_id=0,
            enforce_all_strata=False
        )
        self.assertEqual(df.condition.to_list().count('AC'), 0)
        self.assertEqual(df.condition.to_list().count('PL'), 100)
        self.assertEqual(df.guess.to_list().count('AC'), 0)
        self.assertEqual(df.guess.to_list().count('PL'), 100)

        model = model_defs.ModelDefinition(p_act=0, p_sep=1)
        df = toy_models.ToyModelsDataGenerator.get_toymodel_data(
            output_dir=folders.tmp_dir,
            output_prefix='test_get_toymodel_data_extremes',
            model=model,
            n=100,
            model_sim_id=0,
            enforce_all_strata=False
        )
        self.assertEqual(df.condition.to_list().count('AC'), 0)
        self.assertEqual(df.condition.to_list().count('PL'), 100)
        self.assertEqual(df.guess.to_list().count('AC'), 100)
        self.assertEqual(df.guess.to_list().count('PL'), 0)

        model = model_defs.ModelDefinition(p_act=1, p_sea=0)
        df = toy_models.ToyModelsDataGenerator.get_toymodel_data(
            output_dir=folders.tmp_dir,
            output_prefix='test_get_toymodel_data_extremes',
            model=model,
            n=100,
            model_sim_id=0,
            enforce_all_strata=False
        )
        self.assertEqual(df.condition.to_list().count('AC'), 100)
        self.assertEqual(df.condition.to_list().count('PL'), 0)
        self.assertEqual(df.guess.to_list().count('AC'), 0)
        self.assertEqual(df.guess.to_list().count('PL'), 100)

        model = model_defs.ModelDefinition(p_act=1, p_sea=1)
        df = toy_models.ToyModelsDataGenerator.get_toymodel_data(
            output_dir=folders.tmp_dir,
            output_prefix='test_get_toymodel_data_extremes',
            model=model,
            n=100,
            model_sim_id=0,
            enforce_all_strata=False
        )
        self.assertEqual(df.condition.to_list().count('AC'), 100)
        self.assertEqual(df.condition.to_list().count('PL'), 0)
        self.assertEqual(df.guess.to_list().count('AC'), 100)
        self.assertEqual(df.guess.to_list().count('PL'), 0)

    # Test is probabilistic; try running again if fail
    def test_get_toymodel_data_normals(self, n=1000, thr=3): # Case 1 - normal p_act; normal p_se(a/p)

        model_params = [
            {'p_act':0.5, 'p_se':0.5},  # mid p_act; mid p_se(a/p)
            {'p_act':0.2, 'p_se':0.5},  # low p_act; mid p_se(a/p)
            {'p_act':0.5, 'p_se':0.2},  # mid p_act; low p_se(a/p)
            {'p_act':0.85, 'p_se':0.5}, # high p_act; mid p_se(a/p)
            {'p_act':0.5, 'p_se':0.85}, # mid p_act; high p_se(a/p)
        ]

        for model_param in model_params:

            p_act = model_param['p_act']
            p_se = model_param['p_se']

            df = toy_models.ToyModelsDataGenerator.get_toymodel_data(
                output_dir=folders.tmp_dir,
                output_prefix='test_get_toymodel_data_normals',
                model=model_defs.ModelDefinition(
                    p_act=p_act,
                    p_sep=p_se,
                    p_sea=p_se,),
                n=n,
                model_sim_id=0,
                enforce_all_strata=False
            )

            se = sqrt((1-p_act)*(p_act)/n)
            self.assertTrue(abs(p_act - df.condition.to_list().count('AC')/n) <= thr*se)
            self.assertTrue(abs((1-p_act) - df.condition.to_list().count('PL')/n) <= thr*se)

            se = sqrt((1-p_se)*(p_se)/n)
            self.assertTrue(abs(p_se - df.guess.to_list().count('AC')/n) <= thr*se)
            self.assertTrue(abs((1-p_se) - df.guess.to_list().count('PL')/n) <= thr*se)


""" Heritage code below """
class IntegrationTests(unittest.TestCase):

    @unittest.skip('wip')
    def test_get_toymodel_data_case1(self):

        toy_models.Controllers.run_toymodels_cgrc(
            model_family_name='default_models',
            n_trials=1,
            n_patients=50,
            postfix='test',
            cgrc_param_set=1,
            save_figs=False,
        )
