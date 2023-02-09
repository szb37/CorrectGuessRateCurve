"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

There is inconsistency in how the model and strata stats are extracted from R:
- for strata, R output is converted to JSON which is then parsed
- for model, R output is directly extracted form the R output

TODO:
- sync data extract methods, see above
- move all R related code to r files and call them when appropiate from python
"""

import src.dataframe_classes as df_class
import src.config as config
import src.miscs as miscs
from itertools import product as product
from rpy2.robjects import r
import pandas as pd
import json
import copy
import os


r('library(lmerTest)')
r('library("emmeans")')
r('library("dplyr")')
r('library(RJSONIO)')


class Controllers():

    @staticmethod
    def get_trial_stats(input_dir, input_fname, output_dir, output_prefix, do_stratas=False, trial_scales=None, add_columns=None):
        ''' Gets stats from processedDf and save results as CSV
            Args:
                input_dir (str): input file's directory
                input_fname (str): input filename
                output_dir (str): where to save outputs
                output_prefix (str): prefix to label outputs
                trial_scales (dict): key-values of trials-scales for which stats will be calculated
                    E.g. trial_scales = {'trial2':['scale1', 'scale3']}, will process scales 1 and 3 of trial2
                add_columns(dict, optional): any key-value pair to add as column to output df
        '''

        assert isinstance(input_dir, str)
        assert isinstance(input_fname, str)
        assert isinstance(output_dir, str)
        assert isinstance(output_prefix, str)
        assert isinstance(do_stratas, bool)
        assert (isinstance(trial_scales, dict) or (trial_scales is None))
        assert (isinstance(add_columns, dict)) or (add_columns is None)

        # Initalize output
        model_summary_list = []
        model_comps_list = []
        strata_summary_list = []
        strata_contrast_list = []

        # Read dataframe into python & R (and set baseline to be PL)
        input_trial_data_fpath = os.path.join(input_dir, input_fname).replace('\\', '/')
        trial_data_df = pd.read_csv(input_trial_data_fpath)
        Helpers.load_df_into_R_space(input_trial_data_fpath)

        # Get trial_scales from trial_data input
        trials, scales = miscs.get_trial_scales(input_df=trial_data_df, trial_scales=trial_scales)

        for trial, scale in product(trials, scales):

            Helpers.get_df_filtered(trial, scale)
            if r('nrow(df_filtered)')[0]==0:
                continue

            model_summary, model_comps = StatsCore.get_model_stats(add_columns)

            model_summary.trial = trial
            model_summary.scale = scale
            model_summary_list.append(model_summary)

            model_comps.trial = trial
            model_comps.scale = scale
            model_comps_list.append(model_comps)

            if do_stratas:
                strata_summary, strata_contrast = StatsCore.get_strata_stats(add_columns)

                strata_summary.trial = trial
                strata_summary.scale = scale
                strata_summary_list.append(strata_summary)

                strata_contrast.trial = trial
                strata_contrast.scale = scale
                strata_contrast_list.append(strata_contrast)

            r('rm(df_filtered)')

        # Clear R namespace
        r('rm(list = ls())')

        model_summary, model_comps, strata_summary, strata_contrast = Helpers.normalize_trial_stats_dfs(
            model_summary_list,
            model_comps_list,
            strata_summary_list,
            strata_contrast_list,
            do_stratas)
        Helpers.save_trial_stats_dfs(
            input_dir,
            input_fname,
            output_dir,
            output_prefix,
            model_summary,
            model_comps,
            strata_summary,
            strata_contrast,
            do_stratas)

    @staticmethod
    def get_cgrc_stats(input_dir, input_fname, output_dir, output_prefix, do_stratas=False, trial_scales=None, add_columns=None):
        """ Gets stats from processedDf (with BBC specific columns, cgr and cgr_sim_id) and save results as CSV
            Args:
                input_dir (str): input file's directory
                input_fname (str): input filename
                output_dir (str): where to save outputs
                output_prefix (str): prefix to label outputs
                trial_scales (dict): key-values of trials-scales for which stats will be calculated
                    E.g. trial_scales = {'trial2':['scale1', 'scale3']}, will process scales 1 and 3 of trial2
                add_columns(dict, optional): any key-value pair to add as column to output df
        """

        assert isinstance(input_dir, str)
        assert isinstance(input_fname, str)
        assert isinstance(output_dir, str)
        assert isinstance(output_prefix, str)
        assert (isinstance(trial_scales, dict) or (trial_scales is None))
        assert (isinstance(add_columns, dict)) or (add_columns is None)

        # Read dataframe into R and set baseline to be PL
        input_fpath = os.path.join(input_dir, input_fname).replace('\\', '/')
        Helpers.load_df_into_R_space(input_fpath)

        # Extracts cgrs and cgr_sim_ids from df
        cgrc_data_df = pd.read_csv(input_fpath)
        cgr_sim_ids = cgrc_data_df.cgr_sim_id.unique().tolist()
        cgrs = cgrc_data_df.cgr.unique().tolist()

        # Initalize output
        model_summary_list = []
        model_comps_list = []
        strata_summary_list = []
        strata_contrast_list = []

        trials, scales = miscs.get_trial_scales(input_df=cgrc_data_df, trial_scales=trial_scales)

        for trial, scale, cgr, cgr_sim_id, in product(trials, scales, cgrs, cgr_sim_ids):

            Helpers.get_df_filtered(trial, scale)
            if r('nrow(df_filtered)')[0]==0:
                continue

            model_summary, model_comps = StatsCore.get_model_stats(add_columns)

            model_summary.trial = trial
            model_summary.scale = scale
            model_summary_list.append(model_summary)

            model_comps.trial = trial
            model_comps.scale = scale
            model_comps_list.append(model_comps)

            if do_stratas:
                strata_summary, strata_contrast = StatsCore.get_strata_stats(add_columns)

                strata_summary.trial = trial
                strata_summary.scale = scale
                strata_summary_list.append(strata_summary)

                strata_contrast.trial = trial
                strata_contrast.scale = scale
                strata_contrast_list.append(strata_contrast)

            r('rm(df_filtered)')

        # Clear R namespace
        r('rm(list = ls())')

        model_summary, model_comps, strata_summary, strata_contrast = Helpers.normalize_trial_stats_dfs(
            model_summary_list,
            model_comps_list,
            strata_summary_list,
            strata_contrast_list,
            do_stratas)
        Helpers.save_trial_stats_dfs(
            input_dir,
            input_fname,
            output_dir,
            output_prefix,
            model_summary,
            model_comps,
            strata_summary,
            strata_contrast,
            do_stratas)


class StatsCore():
    ''' Functions to extract starta/model stats using R '''

    @staticmethod
    def get_model_stats(add_columns):
        ''' Returns model stats dataframe '''

        assert (isinstance(add_columns, dict) or (add_columns is None))

        # Initalize output
        model_summary = df_class.ModelSummaryDf()
        model_summary.add_columns(add_columns)
        model_comps = df_class.ModelComponentsDf()
        model_comps.add_columns(add_columns)

        r("without_guess=lm(formula='delta_score~condition', df_filtered)")
        r("with_guess=lm(formula='delta_score~condition+guess+condition*guess', df_filtered)")

        for model_type in ['without_guess', 'with_guess']:

            # extract model summary df
            r('model_sum = summary({})'.format(model_type))
            model_summary_fromR = Helpers.get_model_summary_stats()
            model_summary_fromR['trial'] = [None]
            model_summary_fromR['scale'] = [None]
            model_summary_fromR['model_type'] = [model_type]
            model_summary = pd.concat(
                [model_summary, pd.DataFrame.from_dict(model_summary_fromR)], sort=False)

            # extract model components df
            comps = r('rownames(model_sum$coefficients)')
            for comp in comps:

                model_comps_fromR = Helpers.get_model_component_stats(
                    comp=comp)

                if comp == '(Intercept)':
                    comp = 'intercept'

                model_comps_fromR['trial'] = [None]
                model_comps_fromR['scale'] = [None]
                model_comps_fromR['model_type'] = [model_type]
                model_comps_fromR['component'] = [comp]
                model_comps = pd.concat(
                    [model_comps, pd.DataFrame.from_dict(model_comps_fromR)], sort=False)

            r('rm(model_sum)')

        r('rm(without_guess)')
        r('rm(with_guess)')

        model_summary.index = range(model_summary.index.shape[0])
        model_comps.index = range(model_comps.index.shape[0])

        return model_summary, model_comps

    @staticmethod
    def get_strata_stats(add_columns):
        ''' Returns strata contrast and summary dataframes '''

        assert (isinstance(add_columns, dict) or (add_columns is None))

        strata_summary = df_class.StrataSummaryDf()
        strata_contrast = df_class.StrataContrastDf()

        py_df_filtered = Helpers.r2pyjson('df_filtered')
        if 'cgr' in py_df_filtered.keys():
            assert len(set(py_df_filtered['cgr'])) == 1
            if (0 in set(py_df_filtered['cgr'])) or (1 in set(py_df_filtered['cgr'])):
                # By definition cannot calc strata stats if CGR==0 or 1
                return strata_summary, strata_contrast

        # Calc strata outputs
        r("with_guess=lm(formula='delta_score~condition+guess+condition*guess', df_filtered)")

        # stratas with fixed guess
        r('emmFixGuess = emmeans(with_guess, specs = pairwise ~ condition|guess)')
        # stratas with fixed condition
        r('emmFixCond  = emmeans(with_guess, specs = pairwise ~ guess|condition)')

        # Calc strata outputs - w Tukey p-value adj for multiple comparisons
        r('emm_Tukey = emmeans(with_guess, specs= c("condition", "guess"))')
        tukey_contrasts = Helpers.r2pyjson('summary(pairs(emm_Tukey))')
        tukey_contrasts = Helpers.format_Tukey_contrast(tukey_contrasts)
        r('rm(with_guess)')
        r('rm(emm_Tukey)')

        for emmName in ['emmFixCond', 'emmFixGuess']:

            # Get strata contrast numbers
            r("contrastSummary = summary({}$contrasts)".format(emmName))
            starata_contrast_fromR = Helpers.r2pyjson('contrastSummary')
            starata_contrast_fromR = Helpers.format_comparisons(
                starata_contrast_fromR, tukey_contrasts)
            r('rm(contrastSummary)')

            # Format strata contrast to be compatible with df
            starata_contrast_fromR['est'] = starata_contrast_fromR.pop(
                'estimate')
            starata_contrast_fromR['se'] = starata_contrast_fromR.pop('SE')
            starata_contrast_fromR['t'] = starata_contrast_fromR.pop('t.ratio')
            starata_contrast_fromR['p'] = starata_contrast_fromR.pop('p.value')
            starata_contrast_fromR['contrast'] = starata_contrast_fromR.pop(
                'comparison')
            starata_contrast_fromR['p_adj'] = starata_contrast_fromR.pop(
                'adj_p')
            starata_contrast_fromR['trial'] = [None, None]
            starata_contrast_fromR['scale'] = [None, None]
            strata_contrast = pd.concat(
                [strata_contrast, pd.DataFrame.from_dict(starata_contrast_fromR)], sort=False)

            # Get strata summary
            if emmName == 'emmFixCond':
                continue  # No need to calculate strata stats

            r("strataSummary = summary({}$emmeans)".format(emmName))
            starata_summary_fromR = Helpers.r2pyjson('strataSummary')
            r('rm(strataSummary)')

            # Format strata summary to be compatible with df
            starata_summary_fromR['est'] = starata_summary_fromR.pop('emmean')
            starata_summary_fromR['se'] = starata_summary_fromR.pop('SE')
            starata_summary_fromR['lower_CI'] = starata_summary_fromR.pop(
                'lower.CL')
            starata_summary_fromR['upper_CI'] = starata_summary_fromR.pop(
                'upper.CL')
            starata_summary_fromR['trial'] = [None, None, None, None]
            starata_summary_fromR['scale'] = [None, None, None, None]
            starata_summary_fromR['strata'] = [
                starata_summary_fromR['condition'][0] +
                starata_summary_fromR['guess'][0],
                starata_summary_fromR['condition'][1] +
                starata_summary_fromR['guess'][1],
                starata_summary_fromR['condition'][2] +
                starata_summary_fromR['guess'][2],
                starata_summary_fromR['condition'][3] +
                starata_summary_fromR['guess'][3],
            ]
            del starata_summary_fromR['condition']
            del starata_summary_fromR['guess']

            strata_summary = pd.concat(
                [strata_summary, pd.DataFrame.from_dict(starata_summary_fromR)], sort=False)

        r('rm(emmFixGuess)')
        r('rm(emmFixCond)')
        del strata_contrast['condition']
        del strata_contrast['guess']

        strata_contrast.index = range(strata_contrast.index.shape[0])
        strata_summary.index = range(strata_summary.index.shape[0])

        return strata_summary, strata_contrast


class Helpers():
    ''' Various helper functions '''

    @staticmethod
    def normalize_trial_stats_dfs(model_summary_list, model_comps_list, strata_summary_list, strata_contrast_list, do_stratas=False):

        # Concatanate dataframes
        model_summary = pd.concat(model_summary_list, axis=0)
        model_comps = pd.concat(model_comps_list, axis=0)
        if do_stratas:
            strata_summary = pd.concat(strata_summary_list, axis=0)
            strata_contrast = pd.concat(strata_contrast_list, axis=0)

        # Set class
        model_summary.__class__ = df_class.ModelSummaryDf
        model_comps.__class__ = df_class.ModelComponentsDf
        if do_stratas:
            strata_summary.__class__ = df_class.StrataSummaryDf
            strata_contrast.__class__ = df_class.StrataContrastDf

        # Set columns
        model_summary.set_column_types()
        model_comps.set_column_types()
        if do_stratas:
            strata_summary.set_column_types()
            strata_contrast.set_column_types()

        # Check validty
        #assert master_model_summary.is_valid()
        #assert master_model_comps.is_valid()
        #if do_stratas:
        #    assert master_strata_summary.is_valid()
        #    assert master_strata_contrast.is_valid()

        if do_stratas:
            return model_summary, model_comps, strata_summary, strata_contrast

        return model_summary, model_comps, None, None

    @staticmethod
    def save_trial_stats_dfs(input_dir, input_fname, output_dir, output_prefix, model_summary, model_comps, strata_summary, strata_contrast, do_stratas=False):

        model_summary.to_csv(os.path.join(output_dir, output_prefix+'__model_summary.csv'), index=False)
        model_comps.to_csv(os.path.join(output_dir, output_prefix+'__model_comps.csv'), index=False)

        if do_stratas:
            strata_summary.to_csv(os.path.join(output_dir, output_prefix+'__strata_summary.csv'), index=False)
            strata_contrast.to_csv(os.path.join(output_dir, output_prefix+'__strata_contrast.csv'), index=False)

    @staticmethod
    def load_df_into_R_space(input_fpath, relevel_sex=False):
        """ Loads dataframe into R global space

            Args:
                input_fpath (str): filepath to dataframe saved as CSV
                relevel_sex (bool, optional): if True, then 'sex' is converted to factorial with male as baseline
        """

        r('df = read.csv("'+input_fpath+'")')

        r('df$condition = as.factor(df$condition)')
        r('df <- within(df, condition <- relevel(condition, ref="PL"))')

        r('df$guess = as.factor(df$guess)')
        r('df <- within(df, guess <- relevel(guess, ref="PL"))')

        if relevel_sex:
            r('df$sex = as.factor(df$sex)')
            r('df <- within(df, sex <- relevel(sex, ref="M"))')

    @staticmethod
    def r2pyjson(r_var):
        """ Converts an R object to Python json

            Args:
                r_var(str): name of the variable in R global space

            Returns:
                JSON of r_var (in Python namespace)
        """

        rjson = r('toJSON({})'.format(r_var))
        return json.loads(rjson[0])

    @staticmethod
    def get_df_filtered(trial, scale, cgr=None, cgr_sim_id=None):
        """ Selects subset of R dataframe.
            It is assumed that 'df' exists in R global space and it is an instance of the XYZ
            dataframe types. This functions filters df by trial, scale, cgr and cgr_sim_id
            Nothing is returned, but df_filtered is created in R global space.

            Args:
                trial (str): name of the trial. If 'all', then
                scale (str): name of the scale.
                cgr (float, optional): break blind ratio, ignored if None
                cgr_sim_id (int, optional): trial index if bbc_engine DF is the input, ignored if None
        """

        if (trial in ['all', 'sbmd']) and (cgr is None) and (cgr_sim_id is None):
            filter_string = 'df_filtered = filter(df, scale=="{}")'.format(scale)

        elif (trial in ['all', 'sbmd']) and (cgr is not None) and (cgr_sim_id is not None):
            filter_string = 'df_filtered = filter(df, scale=="{}" & cgr=={} & cgr_sim_id=={})'.format(
                scale, cgr, cgr_sim_id)

        elif (trial not in ['all', 'sbmd']) and (cgr is None) and (cgr_sim_id is None):
            filter_string = 'df_filtered = filter(df, trial=="{}" & scale=="{}")'.format(
                trial, scale)

        elif (trial not in ['all', 'sbmd']) and (cgr is not None) and (cgr_sim_id is not None):
            filter_string = 'df_filtered = filter(df, trial=="{}" & scale=="{}" & cgr=={} & cgr_sim_id=={})'.format(
                trial, scale, cgr, cgr_sim_id)

        else:
            assert False  # Invalid input

        r('{}'.format(filter_string))

    @staticmethod
    def get_model_summary_stats():
        """ Returns stats for model summary """

        f = r('model_sum$fstatistic[1]')
        df1 = r('model_sum$fstatistic[2]')
        df2 = r('model_sum$fstatistic[3]')
        adjr2 = r('model_sum$adj.r.squared')
        r('fstats = model_sum$fstatistic')
        p = r('pf(fstats[1], fstats[2], fstats[3], lower=FALSE)')

        summaryStatsDict = {'f': f[0], 'df1': df1[0],
                            'df2': df2[0], 'adjr2': adjr2[0], 'p': p[0]}
        return summaryStatsDict

    @staticmethod
    def get_model_component_stats(comp):
        """ Returns model component stats """

        n_row = r('which(rownames(model_sum$coefficients)=="{}")'.format(comp))

        n_col = r(
            'which(colnames(model_sum$coefficients)=="{}")'.format("Estimate"))
        est = r('model_sum$coefficients[{},{}]'.format(n_row[0], n_col[0]))

        n_col = r(
            'which(colnames(model_sum$coefficients)=="{}")'.format("Std. Error"))
        se = r('model_sum$coefficients[{},{}]'.format(n_row[0], n_col[0]))

        n_col = r('which(colnames(model_sum$coefficients)=="{}")'.format("t value"))
        t = r('model_sum$coefficients[{},{}]'.format(n_row[0], n_col[0]))

        n_col = r(
            'which(colnames(model_sum$coefficients)=="{}")'.format("Pr(>|t|)"))
        p = r('model_sum$coefficients[{},{}]'.format(n_row[0], n_col[0]))

        componentStatsDict = {'est': est[0], 'se': se[0], 't': t[0], 'p': p[0]}
        return componentStatsDict

    @staticmethod
    def format_comparisons(py_contrast_dict, tukey_contrasts):
        """ Helper function to format compariosn strings """

        contrast = py_contrast_dict
        contrast['comparison'] = []
        contrast['type'] = []
        contrast['adj_p'] = []

        # Fixed guess case
        if 'guess' in contrast.keys():

            for idx, guess in enumerate(contrast['guess']):
                condition1 = contrast['contrast'][idx][:2]
                condition2 = contrast['contrast'][idx][5:]
                assert condition1 in ['PL', 'AC']
                assert condition2 in ['PL', 'AC']

                comparison = condition1 + guess + 'vs' + condition2 + guess
                contrast['type'].append('fixGuess')
                contrast['comparison'].append(comparison)

                tukey_idx = Helpers.find_contrast_idx(
                    tukey_contrasts, comparison)
                contrast['adj_p'].append(tukey_contrasts['p.value'][tukey_idx])

        # Fixed condition case
        if 'condition' in contrast.keys():

            for idx, condition in enumerate(contrast['condition']):
                guess1 = contrast['contrast'][idx][:2]
                guess2 = contrast['contrast'][idx][5:]
                assert guess1 in ['PL', 'AC']
                assert guess2 in ['PL', 'AC']

                comparison = condition + guess1 + 'vs' + condition + guess2
                contrast['type'].append('fixCondition')
                contrast['comparison'].append(comparison)

                tukey_idx = Helpers.find_contrast_idx(
                    tukey_contrasts, comparison)
                contrast['adj_p'].append(tukey_contrasts['p.value'][tukey_idx])

        return contrast

    @staticmethod
    def find_contrast_idx(tukey_contrasts, comparison):
        """ Returns ID of comparison (str) within tukey_contrasts """

        strata1 = comparison[:4]
        strata2 = comparison[6:]
        # same comparison wit strata order switched
        alt_comparison = strata2+'vs'+strata1

        id = None
        for idx, contrast in enumerate(tukey_contrasts['contrast']):
            if (contrast == comparison) or (contrast == alt_comparison):
                id = idx

        if id is None:
            assert False

        return id

    @staticmethod
    def format_Tukey_contrast(tukey_contrasts):
        """ Format Tukey contrast output """

        formatted_contrasts = copy.deepcopy(tukey_contrasts['contrast'])

        for idx, original_contrast in enumerate(tukey_contrasts['contrast']):

            if original_contrast in ['PL,PL - AC,PL', 'AC,PL - PL,PL']:
                formatted_contrasts[idx] = 'PLPLvsACPL'

            elif original_contrast in ['PL,PL - PL,AC', 'PL,AC - PL,PL']:
                formatted_contrasts[idx] = 'PLPLvsPLAC'

            elif original_contrast in ['PL,PL - AC,AC', 'AC,AC - PL,PL']:
                formatted_contrasts[idx] = 'PLPLvsACAC'

            elif original_contrast in ['AC,PL - PL,AC', 'PL,AC - AC,PL']:
                formatted_contrasts[idx] = 'ACPLvsPLAC'

            elif original_contrast in ['AC,PL - AC,AC', 'AC,AC - AC,PL']:
                formatted_contrasts[idx] = 'ACPLvsACAC'

            elif original_contrast in ['PL,AC - AC,AC', 'AC,AC - PL,AC']:
                formatted_contrasts[idx] = 'PLACvsACAC'

            else:
                assert False

        tukey_contrasts['contrast'] = formatted_contrasts
        return tukey_contrasts

    @staticmethod
    def add_metadata_process_cgrc(all_dfs, do_stratas, trial, scale, cgr, cgr_sim_id):
        """ Add metadata to model_summary, model_comps, strata_summary, strata_contrast """

        all_dfs['model_summary'].trial = trial
        all_dfs['model_summary'].scale = scale
        all_dfs['model_summary'].cgr = cgr
        all_dfs['model_summary'].cgr_sim_id = cgr_sim_id

        all_dfs['model_comps'].trial = trial
        all_dfs['model_comps'].scale = scale
        all_dfs['model_comps'].cgr = cgr
        all_dfs['model_comps'].cgr_sim_id = cgr_sim_id

        if do_stratas:
            all_dfs['strata_summary'].trial = trial
            all_dfs['strata_summary'].scale = scale
            all_dfs['strata_summary'].cgr = cgr
            all_dfs['strata_summary'].cgr_sim_id = cgr_sim_id

            all_dfs['strata_contrast'].trial = trial
            all_dfs['strata_contrast'].scale = scale
            all_dfs['strata_contrast'].cgr = cgr
            all_dfs['strata_contrast'].cgr_sim_id = cgr_sim_id

        return all_dfs


class NoSelectedRows(Exception):
    pass
