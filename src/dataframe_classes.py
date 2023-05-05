"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

import src.dataframe_classes as df_class
a =  df_class.TrialDataDf()
"""

import pandas as pd
import math
import copy
import abc


class CGRCDataFrame(metaclass=abc.ABCMeta):
    ''' ABC of CGRC related DFs '''

    @abc.abstractmethod
    def set_column_types(self):
        pass

    @abc.abstractmethod
    def is_valid(self):
        pass

    def add_columns(self, add_columns):
        """ Add columns with a fixed value
            Args:
                add_columns (dict): keys are the new column names, values are the value of the column
        """

        assert (isinstance(add_columns, dict) or (add_columns is None))

        if add_columns is None:
            return

        for col, value in add_columns.items():
            self[col] = value


class TrialDataDf(pd.DataFrame, CGRCDataFrame):
    ''' Minimal set of columns for CGRC analysis '''

    def __init__(self):
        super(TrialDataDf, self).__init__(columns=[
            'trial',
            'subject_id',
            'scale',
            'tp',
            'condition', 'guess',
            'baseline', 'score','delta_score',])
        self.set_column_types()

    def set_column_types(self):
        self['trial'] = self['trial'].astype('str')
        self['subject_id'] = self['subject_id'].astype('int')
        self['scale'] = self['scale'].astype('str')
        self['tp'] = self['tp'].astype('str')
        self['condition'] = self['condition'].astype('str')
        self['guess'] = self['guess'].astype('str')
        self['baseline'] = self['baseline'].astype('float')
        self['score'] = self['score'].astype('float')
        self['delta_score'] = self['delta_score'].astype('float')

    def is_valid(self):

        # All conditions / guesses are either 'PL' or 'AC'
        assert all([condition in ['PL', 'AC'] for condition in self.condition.to_list()])
        assert all([guess in ['PL', 'AC'] for guess in self.guess.to_list()])

        # Assert all baseline / scores / delta_scores are numbers
        assert all([
            (isinstance(baseline, float) and not math.isnan(baseline)) for baseline in self.baseline.to_list()])
        assert all([
            (isinstance(score, float) and not math.isnan(score)) for score in self.score.to_list()])
        assert all([
            (isinstance(delta_score, float) and not math.isnan(delta_score)) for delta_score in self.delta_score.to_list()])

        self._def_check_entry_uniqueness()
        return True

    def _def_check_entry_uniqueness(self):
        # All combinations of trial/subject_id/scale/tp should be unique (add model_sim_id when )

        if 'model_sim_id' in self.columns:
            temp = copy.deepcopy(self).loc[:, ['trial', 'model_sim_id', 'subject_id', 'scale', 'tp']]
            assert temp.shape[0] == temp.drop_duplicates().shape[0]
        else:
            temp = copy.deepcopy(self).loc[:, ['trial', 'subject_id', 'scale', 'tp']]
            assert temp.shape[0] == temp.drop_duplicates().shape[0]


class CGRCurveDf(TrialDataDf):
    ''' DataFrame for holding CGRC data '''

    def __init__(self):
        super(TrialDataDf, self).__init__(columns=[
            'trial',
            'cgr',
            'cgr_sim_id',
            'scale',
            'tp',
            'condition',
            'guess',
            'baseline',
            'score',
            'delta_score'])
        self.set_column_types()

    def set_column_types(self):
        self['trial'] = self['trial'].astype('str')
        self['cgr'] = self['cgr'].astype('float')
        self['cgr_sim_id'] = self['cgr_sim_id'].astype('int')
        #self['subject_id'] = self['subject_id'].astype('int')
        self['scale'] = self['scale'].astype('str')
        self['tp'] = self['tp'].astype('str')
        self['condition'] = self['condition'].astype('str')
        self['guess'] = self['guess'].astype('str')
        self['baseline'] = self['baseline'].astype('float')
        self['score'] = self['score'].astype('float')
        self['delta_score'] = self['delta_score'].astype('float')
        if 'model_sim_id' in self.columns:
            self['model_sim_id'] = self['model_sim_id'].astype('int')

    def is_valid(self):
        #super(CGRCurveDf, self).is_valid()
        return True


class ModelSummaryDf(pd.DataFrame, CGRCDataFrame):
    ''' DataFrame for holding model summary output '''

    def __init__(self):
        super(ModelSummaryDf, self).__init__(columns=[
            'trial', 'scale', 'model_type', 'df1', 'df2', 'f', 'adjr2', 'p'])

    def set_column_types(self):

        self['trial'] = self['trial'].astype('str')
        self['scale'] = self['scale'].astype('str')
        self['model_type'] = self['model_type'].astype('str')
        self['df1'] = self['df1'].astype('float')
        self['df2'] = self['df2'].astype('float')
        self['f'] = self['f'].astype('float')
        self['adjr2'] = self['adjr2'].astype('float')
        self['p'] = self['p'].astype('float')

        if 'cgr' in self.columns:
            self['cgr'] = self['cgr'].astype('float')

        if 'cgr_sim_id' in self.columns:
            self['cgr_sim_id'] = self['cgr_sim_id'].astype('int')

    def is_valid(self):
        self.set_column_types()
        return True


class ModelComponentsDf(pd.DataFrame, CGRCDataFrame):
    ''' DataFrame for holding model components output '''

    def __init__(self):
        super(ModelComponentsDf, self).__init__(columns=[
            'trial', 'scale', 'model_type', 'component', 'est', 'se', 't', 'p'])

    def set_column_types(self):

        self['trial'] = self['trial'].astype('str')
        self['scale'] = self['scale'].astype('str')
        self['model_type'] = self['model_type'].astype('str')
        self['component'] = self['component'].astype('str')
        self['est'] = self['est'].astype('float')
        self['se'] = self['se'].astype('float')
        self['t'] = self['t'].astype('float')
        self['p'] = self['p'].astype('float')

        if 'cgr' in self.columns:
            self['cgr'] = self['cgr'].astype('float')

        if 'cgr_sim_id' in self.columns:
            self['cgr_sim_id'] = self['cgr_sim_id'].astype('int')

        if 'model_sim_id' in self.columns:
            self['model_sim_id'] = self['model_sim_id'].astype('int')

    def is_valid(self):
        self.set_column_types()
        return True


class StrataSummaryDf(pd.DataFrame, CGRCDataFrame):
    ''' DataFrame for holding model output by strata '''

    def __init__(self):
        super(StrataSummaryDf, self).__init__(columns=[
            'trial', 'scale', 'strata', 'est', 'se', 'df', 'lower_CI', 'upper_CI'])

    def set_column_types(self):

        self['trial'] = self['trial'].astype('str')
        self['scale'] = self['scale'].astype('str')
        self['strata'] = self['strata'].astype('str')
        self['est'] = self['est'].astype('float')
        self['se'] = self['se'].astype('float')
        self['df'] = self['df'].astype('float')
        self['lower_CI'] = self['lower_CI'].astype('float')
        self['upper_CI'] = self['upper_CI'].astype('float')

        if 'cgr' in self.columns:
            self['cgr'] = self['cgr'].astype('float')

        if 'cgr_sim_id' in self.columns:
            self['cgr_sim_id'] = self['cgr_sim_id'].astype('int')

        if 'model_sim_id' in self.columns:
            self['model_sim_id'] = self['model_sim_id'].astype('int')

    def is_valid(self):
        self.set_column_types()
        return True


class StrataContrastDf(pd.DataFrame, CGRCDataFrame):
    ''' DataFrame for holding model output by strata '''

    def __init__(self):
        super(StrataContrastDf, self).__init__(columns=[
            'trial', 'scale', 'contrast', 'type', 'est', 'se', 'df', 't', 'p', 'p_adj'])

    def set_column_types(self):

        self['trial'] = self['trial'].astype('str')
        self['scale'] = self['scale'].astype('str')
        self['type'] = self['type'].astype('str')
        self['est'] = self['est'].astype('float')
        self['se'] = self['se'].astype('float')
        self['df'] = self['df'].astype('float')
        self['t'] = self['t'].astype('float')
        self['p'] = self['p'].astype('float')
        self['p_adj'] = self['p_adj'].astype('float')

        if 'cgr' in self.columns:
            self['cgr'] = self['cgr'].astype('float')

        if 'cgr_sim_id' in self.columns:
            self['cgr_sim_id'] = self['cgr_sim_id'].astype('int')

        if 'model_sim_id' in self.columns:
            self['model_sim_id'] = self['model_sim_id'].astype('int')

    def is_valid(self):
        self.set_column_types()
        return True
