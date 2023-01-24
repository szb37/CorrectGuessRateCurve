"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2021, DrugNerdsLab
:License: MIT
"""

import src.cgrc.covariate_selector as covariate_selector
import src.folders as folders
from rpy2.robjects import r
import unittest
import os


class CovariateSelectorCoreUnitTests(unittest.TestCase):

    def test_is_all_covs_valid(self):

        # case1; [is any cov non-sig?]: True; [is all non-sig cov enforced?]: True
        model_summary = {
            'age': {'p': 0.5},
            'sex': {'p': 0.6},
            'eth': {'p': 0.7},
            'intercept': {'p': 0.8}, }
        isAllCovariateValid = covariate_selector.CovariateSelectorCore.is_all_covs_valid(
            model_summary,
            enforce_covs=['age', 'sex', 'eth'],
            enforce_int=True)
        self.assertTrue(isAllCovariateValid)

        # case2; [is any cov non-sig?]: False; [is all non-sig cov enforced?]: True
        model_summary = {
            'age': {'p': 0.005},
            'sex': {'p': 0.006},
            'eth': {'p': 0.007},
            'intercept': {'p': 0.8}, }
        isAllCovariateValid = covariate_selector.CovariateSelectorCore.is_all_covs_valid(
            model_summary,
            enforce_covs=['tadaa'],
            enforce_int=True)
        self.assertTrue(isAllCovariateValid)

        # case3; [is any cov non-sig?]: True; [is all non-sig cov enforced?]: False
        model_summary = {
            'age': {'p': 0.005},
            'sex': {'p': 0.6},
            'eth': {'p': 0.007},
            'intercept': {'p': 0.8}, }
        isAllCovariateValid = covariate_selector.CovariateSelectorCore.is_all_covs_valid(
            model_summary,
            enforce_covs=[],
            enforce_int=True)
        self.assertFalse(isAllCovariateValid)

        # case4; [is any cov non-sig?]: False; [is all non-sig cov enforced?]: False
        model_summary = {
            'age': {'p': 0.005},
            'sex': {'p': 0.6},
            'eth': {'p': 0.007},
            'intercept': {'p': 0.8}, }
        isAllCovariateValid = covariate_selector.CovariateSelectorCore.is_all_covs_valid(
            model_summary,
            enforce_covs=['age'],
            enforce_int=True)
        self.assertFalse(isAllCovariateValid)

    def test_rm_least_sig(self):

        # case1; [is any cov non-sig?]: True; [is least sig enforced?]: True
        model_summary = {
            'age': {'p': 0.5},
            'sex': {'p': 0.6},
            'eth': {'p': 0.7},
            'intercept': {'p': 0.8}, }
        try_covs = covariate_selector.CovariateSelectorCore.rm_least_sig(
            model_summary=model_summary,
            try_covs=['age', 'sex', 'eth'],
            enforce_covs=['eth'],
            enforce_int=True)
        self.assertTrue(try_covs == ['age', 'eth'])

        # case2; [is any cov non-sig?]: False; [is least sig enforced?]: True
        model_summary = {
            'age': {'p': 0.005},
            'sex': {'p': 0.006},
            'eth': {'p': 0.007},
            'intercept': {'p': 0.8}, }
        with self.assertRaises(AssertionError):
            try_covs = covariate_selector.CovariateSelectorCore.rm_least_sig(
                model_summary=model_summary,
                try_covs=['age', 'sex', 'eth'],
                enforce_covs=['eth'],
                enforce_int=True)

        # case3; [is any cov non-sig?]: True; [is least sig enforced?]: False
        model_summary = {
            'age': {'p': 0.005},
            'sex': {'p': 0.6},
            'eth': {'p': 0.7},
            'intercept': {'p': 0.8}, }
        try_covs = covariate_selector.CovariateSelectorCore.rm_least_sig(
            model_summary=model_summary,
            try_covs=['age', 'sex', 'eth'],
            enforce_covs=['age'],
            enforce_int=True)
        self.assertTrue(try_covs == ['age', 'sex'])

        # case4; [is any cov non-sig?]: False; [is least sig enforced?]: False
        model_summary = {
            'age': {'p': 0.005},
            'sex': {'p': 0.6},
            'eth': {'p': 0.7},
            'intercept': {'p': 0.8}, }
        try_covs = covariate_selector.CovariateSelectorCore.rm_least_sig(
            model_summary=model_summary,
            try_covs=['age', 'sex', 'eth'],
            enforce_covs=['age'],
            enforce_int=True)
        self.assertTrue(try_covs == ['age', 'sex'])


class CovariateSelectorCoreIntegrationTests(unittest.TestCase):

    def test_get_adjusted_model1(self):

        # variable A is sig, B is not
        input_fpath = os.path.join(
            folders.fixtures, 'get_adjusted_model_input1.csv').replace('\\', '/')
        r('df_filtered = read.csv("'+input_fpath+'")')

        model_name = 'test_model1'

        covariate_selector.Controllers.get_adjusted_model(
            model_name=model_name,
            try_covs=[],
            enforce_covs=['A', 'B'])
        summary = covariate_selector.CovariateSelectorCore.get_pystandrd_model_summary(
            model_name)
        self.assertTrue('A' in summary.keys())
        self.assertTrue('B' in summary.keys())

        covariate_selector.Controllers.get_adjusted_model(
            model_name=model_name,
            try_covs=['A', 'B'],
            enforce_covs=[])
        summary = covariate_selector.CovariateSelectorCore.get_pystandrd_model_summary(
            model_name)
        self.assertTrue('A' in summary.keys())
        self.assertTrue('B' not in summary.keys())

        covariate_selector.Controllers.get_adjusted_model(
            model_name=model_name,
            try_covs=['A'],
            enforce_covs=['B'])
        summary = covariate_selector.CovariateSelectorCore.get_pystandrd_model_summary(
            model_name)
        self.assertTrue('A' in summary.keys())
        self.assertTrue('B' in summary.keys())

        covariate_selector.Controllers.get_adjusted_model(
            model_name=model_name,
            try_covs=['B'],
            enforce_covs=['A'])
        summary = covariate_selector.CovariateSelectorCore.get_pystandrd_model_summary(
            model_name)
        self.assertTrue('A' in summary.keys())
        self.assertTrue('B' not in summary.keys())
