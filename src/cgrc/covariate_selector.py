"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2020, DrugNerdsLab
:License: MIT
"""

import src.cgrc.stats as stats
from rpy2.robjects import r


class Controllers():

    @staticmethod
    def get_adjusted_model(model_name, try_covs, enforce_covs, enforce_int=True, dependent_variable='delta_score'):
        """ Creates regression model in R where all covariates in enforce_covs are included and all
            try_covs are included if they are significant. The code does backward eliminiation of
            covs, i.e. construct model, removes the cov with highest p and repeat untill all covs are
            significant (or in enforce_covs). Resulting model is in R global space

            Code assumes that df_filtered exists in R global namespace and is the source data

            Args:
                model_name (str): name of the model
                try_covs (list of strs): list of covariates that will be tested for the model. All
                    elements must be columns in df_filtered
                enforce_covs (list of strs): list of covariates that will be added for the model
                    regardless of their p-value. All elements must be columns in df_filtered
                enforce_int (bool, optional): if True model intercept is kept even if its unsignificant
                    TODO: test code for enforce_int=False case
                dependent_variable (str, optional): dep vairable of model; must be columns in df_filtered
        """

        assert isinstance(enforce_covs, list)
        assert isinstance(try_covs, list)
        assert all([isinstance(el, str)
                   for el in enforce_covs]) or enforce_covs == []
        assert all([isinstance(el, str) for el in try_covs]) or try_covs == []

        enforce_covs_str = CovariateSelectorCore.get_formatted_covs_str(
            enforce_covs)

        while True:

            try_covs_str = CovariateSelectorCore.get_formatted_covs_str(
                try_covs)
            covs_str = enforce_covs_str+try_covs_str
            covs_str = covs_str[1:]  # rm leading + sign

            if covs_str == '':
                assert False  # there should be at least one indep variable

            r('{}=lm(formula=delta_score~{}, df_filtered)'.format(model_name, covs_str))
            model_summary = CovariateSelectorCore.get_pystandrd_model_summary(
                model_name)

            isAllCovariateValid = CovariateSelectorCore.is_all_covs_valid(
                model_summary,
                enforce_covs,
                enforce_int=enforce_int
            )

            if isAllCovariateValid:
                return

            try_covs = CovariateSelectorCore.rm_least_sig(
                model_summary,
                try_covs,
                enforce_covs,
                enforce_int=enforce_int
            )


class CovariateSelectorCore():
    """ Methods to select covariates with progressive backwards elminiation (=eliminate least sig) """

    @staticmethod
    def is_all_covs_valid(model_summary, enforce_covs, enforce_int=True):
        """ Returns True if all non-enforced model components are significant; False otherwise """

        isAllCovariateValid = True
        for comp, ests in model_summary.items():

            if comp in enforce_covs:
                continue

            if enforce_int and comp == 'intercept':
                continue

            if ests['p'] > 0.05:
                isAllCovariateValid = False

        return isAllCovariateValid

    @staticmethod
    def rm_least_sig(model_summary, try_covs, enforce_covs, enforce_int=True):
        """ Return try_covs list with the least significant component removed """

        least_sig_comp = 'pseudo'
        least_sig_p = 0

        for comp, ests in model_summary.items():

            if comp in enforce_covs:
                continue

            if enforce_int and comp == 'intercept':
                continue

            if ests['p'] < 0.05:
                continue

            if ests['p'] > least_sig_p:
                least_sig_comp = comp
                least_sig_p = ests['p']

        assert least_sig_comp != 'pseudo'
        try_covs.remove(least_sig_comp)

        return try_covs

    @staticmethod
    def get_formatted_covs_str(covs_list=[]):
        """ Returns string of covariates separated by '+' """

        covs_str = ''
        for cov in covs_list:
            covs_str += '+{}'.format(cov)

        return covs_str

    @staticmethod
    def get_pystandrd_model_summary(model_name):
        """ Convert R model summary to python dict """

        r('model_summary_R = summary({})'.format(model_name))
        model_summary_temp = stats.Helpers.r2pyjson('model_summary_R')
        r('rm(model_summary_R)')

        comps = [el for el in model_summary_temp['aliased'].keys()]
        model_summary = {}

        for idx, comp in enumerate(comps):

            try:
                model_summary_temp['coefficients'][idx]['Estimate']
            except IndexError:
                continue

            model_summary[comp] = {
                'est': model_summary_temp['coefficients'][idx]['Estimate'],
                'se': model_summary_temp['coefficients'][idx]['Std. Error'],
                't': model_summary_temp['coefficients'][idx]['t value'],
                'p': model_summary_temp['coefficients'][idx]['Pr(>|t|)']}

        model_summary = CovariateSelectorCore.conv_covs_names(model_summary)
        return model_summary

    @staticmethod
    def conv_covs_names(model_summary):
        """ Convert R components name to standard """

        if '(Intercept)' in model_summary.keys():
            model_summary['intercept'] = model_summary.pop('(Intercept)')

        if 'sexF' in model_summary.keys():
            model_summary['sex'] = model_summary.pop('sexF')
        elif 'sexM' in model_summary.keys():
            model_summary['sex'] = model_summary.pop('sexM')

        if 'conditionPL' in model_summary.keys():
            model_summary['condition'] = model_summary.pop('conditionPL')
        elif 'conditionAC' in model_summary.keys():
            model_summary['condition'] = model_summary.pop('conditionAC')

        if 'guessPL' in model_summary.keys():
            model_summary['guess'] = model_summary.pop('guessPL')
        elif 'guessAC' in model_summary.keys():
            model_summary['guess'] = model_summary.pop('guessAC')

        if 'conditionAC:guessAC' in model_summary.keys():
            model_summary['condition*guess'] = model_summary.pop(
                'conditionAC:guessAC')
        elif 'conditionPL:guessPL' in model_summary.keys():
            model_summary['condition*guess'] = model_summary.pop(
                'conditionPL:guessPL')
        if 'guessAC:conditionAC' in model_summary.keys():
            model_summary['condition*guess'] = model_summary.pop(
                'guessAC:conditionAC')
        elif 'guessPL:conditionPL' in model_summary.keys():
            model_summary['condition*guess'] = model_summary.pop(
                'guessPL:conditionPL')

        return model_summary
