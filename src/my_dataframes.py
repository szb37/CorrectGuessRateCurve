"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2020, DrugNerdsLab
:License: MIT
"""

import pandas as pd
import math
import copy
import abc


class MyDataframes(metaclass=abc.ABCMeta):
    ''' ABC of CGRC related DFs '''

    @abc.abstractmethod
    def set_column_types(self):
        pass

    @abc.abstractmethod
    def check_assumptions(self):
        pass

    def add_univalue_columns(self, add_columns):
        """ Add columns with a fixed value """

        if add_columns is None:
            return

        assert isinstance(add_columns, dict)

        for col, value in add_columns.items():
            assert isinstance(col, str)
            self[col] = value


class TrialDataDf(pd.DataFrame, MyDataframes):
    ''' Minimal set of columns for CGRC analysis '''

    def __init__(self):
        super(TrialDataDf, self).__init__(columns=[
            'study', 'subject_id', 'scale', 'tp', 'condition', 'guess','baseline', 'score', 'delta_score',
        ])

    def set_column_types(self):
        self['study'] = self['study'].astype('str')
        self['subject_id'] = self['subject_id'].astype('int')
        self['scale'] = self['scale'].astype('str')
        self['tp'] = self['tp'].astype('str')
        self['baseline'] = self['baseline'].astype('float')
        self['score'] = self['score'].astype('float')
        self['delta_score'] = self['delta_score'].astype('float')
        self['condition'] = self['condition'].astype('str')
        self['guess'] = self['guess'].astype('str')

    def check_assumptions(self):

        # All conditions / guesses are either 'PL' or 'AC'
        assert all([condition in ['PL', 'AC'] for condition in self.condition.to_list()])
        assert all([guess in ['PL', 'AC'] for guess in self.guess.to_list()])

        # Assert all baseline / scores / delta_scores are numbers
        assert all([not (math.isnan(score) or score is None) for score in self.score.to_list()])
        assert all([not (math.isnan(delta_score) or delta_score is None) for delta_score in self.delta_score.to_list()])
        assert all([not (math.isnan(baseline) or baseline is None) for baseline in self.baseline.to_list()])




        # Check whether all combinations of subject_id/study/scale/tp are unique by checking if
        # n of rows change after duplicate removal
        temp = copy.deepcopy(self)
        cols_to_keep = ['study', 'subject_id', 'scale', 'tp']
        [cols_to_keep.append(col) for col in ['guesser',
                                              'respondent'] if col in self.columns]
        temp = temp.loc[:, cols_to_keep]

        n_before = temp.shape[0]
        temp.drop_duplicates(inplace=True)
        n_after = temp.shape[0]
        assert n_before == n_after


class CGRCurveDf(pd.DataFrame, MyDataframes):
    ''' DataFrame for holding CGRC data '''

    def __init__(self):
        super(CGRCurveDf, self).__init__(columns=[
            'study', 'scale', 'cgr', 'cgr_trial_id', 'condition', 'guess', 'baseline', 'score', 'delta_score'])

    def set_column_types(self):
        """ Set column types """

        self['study'] = self['study'].astype('str')
        self['scale'] = self['scale'].astype('str')
        self['cgr'] = self['cgr'].astype('float')
        self['cgr_trial_id'] = self['cgr_trial_id'].astype('int')
        self['condition'] = self['condition'].astype('str')
        self['guess'] = self['guess'].astype('str')
        self['baseline'] = self['baseline'].astype('float')
        self['score'] = self['score'].astype('float')
        self['delta_score'] = self['delta_score'].astype('float')

        if 'guesser' in self.columns:
            self['guesser'] = self['guesser'].astype('str')

        if 'respondent' in self.columns:
            self['respondent'] = self['respondent'].astype('str')

        if 'model_id' in self.columns:
            self['model_id'] = self['model_id'].astype('int')

        if 'trial_id' in self.columns:
            self['trial_id'] = self['trial_id'].astype('int')


class ModelSummaryDf(pd.DataFrame, MyDataframes):
    ''' DataFrame for holding model summary output '''

    def __init__(self):
        super(ModelSummaryDf, self).__init__(columns=[
            'study', 'scale', 'guesser', 'respondent', 'model_type', 'df1', 'df2', 'f', 'adjr2', 'p'])

    def set_column_types(self):
        """ Set column types """

        self['study'] = self['study'].astype('str')
        self['scale'] = self['scale'].astype('str')
        self['model_type'] = self['model_type'].astype('str')
        self['df1'] = self['df1'].astype('float')
        self['df2'] = self['df2'].astype('float')
        self['f'] = self['f'].astype('float')
        self['adjr2'] = self['adjr2'].astype('float')
        self['p'] = self['p'].astype('float')

        if 'cgr' in self.columns:
            self['cgr'] = self['cgr'].astype('float')

        if 'cgr_trial_id' in self.columns:
            self['cgr_trial_id'] = self['cgr_trial_id'].astype('int')

        if 'guesser' in self.columns:
            self['guesser'] = self['guesser'].astype('str')

        if 'respondent' in self.columns:
            self['respondent'] = self['respondent'].astype('str')

        if 'model_id' in self.columns:
            self['model_id'] = self['model_id'].astype('int')

        if 'trial_id' in self.columns:
            self['trial_id'] = self['trial_id'].astype('int')


class ModelComponentsDf(pd.DataFrame, MyDataframes):
    ''' DataFrame for holding model components output '''

    def __init__(self):
        super(ModelComponentsDf, self).__init__(columns=[
            'study', 'scale', 'guesser', 'respondent', 'model_type', 'component', 'est', 'se', 't', 'p'])

    def set_column_types(self):
        """ Set column types """

        self['study'] = self['study'].astype('str')
        self['scale'] = self['scale'].astype('str')
        self['model_type'] = self['model_type'].astype('str')
        self['component'] = self['component'].astype('str')
        self['est'] = self['est'].astype('float')
        self['se'] = self['se'].astype('float')
        self['t'] = self['t'].astype('float')
        self['p'] = self['p'].astype('float')

        if 'cgr' in self.columns:
            self['cgr'] = self['cgr'].astype('float')

        if 'cgr_trial_id' in self.columns:
            self['cgr_trial_id'] = self['cgr_trial_id'].astype('int')

        if 'guesser' in self.columns:
            self['guesser'] = self['guesser'].astype('str')

        if 'respondent' in self.columns:
            self['respondent'] = self['respondent'].astype('str')

        if 'model_id' in self.columns:
            self['model_id'] = self['model_id'].astype('int')

        if 'trial_id' in self.columns:
            self['trial_id'] = self['trial_id'].astype('int')


class StrataSummaryDf(pd.DataFrame, MyDataframes):
    ''' DataFrame for holding model output by strata '''

    def __init__(self):
        super(StrataSummaryDf, self).__init__(columns=[
            'study', 'scale', 'guesser', 'respondent', 'strata', 'est', 'se', 'df', 'lower_CI', 'upper_CI'])

    def set_column_types(self):
        """ Set column types """

        self['study'] = self['study'].astype('str')
        self['scale'] = self['scale'].astype('str')
        self['strata'] = self['strata'].astype('str')
        self['est'] = self['est'].astype('float')
        self['se'] = self['se'].astype('float')
        self['df'] = self['df'].astype('float')
        self['lower_CI'] = self['lower_CI'].astype('float')
        self['upper_CI'] = self['upper_CI'].astype('float')

        if 'cgr' in self.columns:
            self['cgr'] = self['cgr'].astype('float')

        if 'cgr_trial_id' in self.columns:
            self['cgr_trial_id'] = self['cgr_trial_id'].astype('int')

        if 'guesser' in self.columns:
            self['guesser'] = self['guesser'].astype('str')

        if 'respondent' in self.columns:
            self['respondent'] = self['respondent'].astype('str')

        if 'model_id' in self.columns:
            self['model_id'] = self['model_id'].astype('int')

        if 'trial_id' in self.columns:
            self['trial_id'] = self['trial_id'].astype('int')


class StrataContrastDf(pd.DataFrame, MyDataframes):
    ''' DataFrame for holding model output by strata '''

    def __init__(self):
        super(StrataContrastDf, self).__init__(columns=[
            'study', 'scale', 'guesser', 'respondent', 'contrast', 'type', 'est', 'se', 'df', 't', 'p', 'p_adj'])

    def set_column_types(self):
        """ Set column types """

        self['study'] = self['study'].astype('str')
        self['scale'] = self['scale'].astype('str')
        self['guesser'] = self['guesser'].astype('str')
        self['type'] = self['type'].astype('str')
        self['est'] = self['est'].astype('float')
        self['se'] = self['se'].astype('float')
        self['df'] = self['df'].astype('float')
        self['t'] = self['t'].astype('float')
        self['p'] = self['p'].astype('float')
        self['p_adj'] = self['p_adj'].astype('float')

        if 'cgr' in self.columns:
            self['cgr'] = self['cgr'].astype('float')

        if 'cgr_trial_id' in self.columns:
            self['cgr_trial_id'] = self['cgr_trial_id'].astype('int')

        if 'guesser' in self.columns:
            self['guesser'] = self['guesser'].astype('str')

        if 'respondent' in self.columns:
            self['respondent'] = self['respondent'].astype('str')

        if 'model_id' in self.columns:
            self['model_id'] = self['model_id'].astype('int')

        if 'trial_id' in self.columns:
            self['trial_id'] = self['trial_id'].astype('int')


class ModelFamilyResultsDf(pd.DataFrame, MyDataframes):
    ''' DataFrame for holding toy models data '''

    def __init__(self):
        super(ModelFamilyResultsDf, self).__init__(columns=[
            'model',
            'n_trials',
            'n_patients',
            'cgrc_param_set',
            'cgr',
            'avg_trt_p',
            'sig_trt_rate',
            'avg_trt_es',
            'cgradj_avg_trt_p',
            'cgradj_sig_trt_rate',
            'cgradj_avg_trt_es',
        ])

    def set_column_types(self):
        """ Set column types """

        self['model'] = self['model'].astype('str')
        self['n_trials'] = self['n_trials'].astype('int')
        self['cgr'] = self['cgr'].astype('float')  # correct guess rate
        self['avg_trt_p'] = self['avg_trt_p'].astype(
            'float')  # average treatment p
        self['sig_trt_rate'] = self['sig_trt_rate'].astype(
            'float')  # % of trials with significant treatment p
        self['avg_trt_es'] = self['avg_trt_es'].astype(
            'float')  # average treatment effect size

        if 'cgrc_param_set' in self.columns:
            self['cgrc_param_set'] = self['cgrc_param_set'].astype('int')

        if 'cgradj_avg_trt_p' in self.columns:
            self['cgradj_avg_trt_p'] = self['cgradj_avg_trt_p'].astype('float')

        if 'cgradj_sig_trt_rate' in self.columns:
            self['cgradj_sig_trt_rate'] = self['cgradj_sig_trt_rate'].astype(
                'float')

        if 'cgradj_avg_trt_es' in self.columns:
            self['cgradj_avg_trt_es'] = self['cgradj_avg_trt_es'].astype(
                'float')
