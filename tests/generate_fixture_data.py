"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

Functions to generate faux data for testing purposes
"""

import src.config as config
import src.folders as folders
import pandas as pd
import numpy as np
import itertools
import random
import os


def gen_faux_processed_bbc1():

    df = pd.DataFrame(columns=[
        'trial', 'subject_id', 'tp', 'scale', 'score', 'delta_score', 'condition', 'guess', 'guess_conf', 'guesser', 'respondent'])

    conditions = ['PL', 'AC']
    guesses = ['PL', 'AC']

    for condition, guess in itertools.product(conditions, guesses):
        for idx in range(500):

            if (condition == 'PL') and (guess == 'PL'):
                delta_score = random.gauss(mu=10, sigma=1)
            elif (condition == 'AC') and (guess == 'PL'):
                delta_score = random.gauss(mu=70, sigma=1)
            elif (condition == 'PL') and (guess == 'AC'):
                delta_score = random.gauss(mu=50, sigma=1)
            elif (condition == 'AC') and (guess == 'AC'):
                delta_score = random.gauss(mu=30, sigma=1)
            else:
                assert False

            df = df.append({
                'trial': 'test',
                'subject_id': idx,
                'tp': 'wk8',
                'scale': 'tadaa',
                'score': 0,
                'delta_score': delta_score,
                'condition': condition,
                'guess': guess,
                'guesser': 'self',
                'respondent': 'self'},
                ignore_index=True)

    df.to_csv(os.path.join(folders.fixtures,
              'pseudodata_get_bbc_from_processed1.csv'), index=False)


def gen_faux_processed_bbc2():

    df = pd.DataFrame(columns=[
        'trial', 'subject_id', 'tp', 'scale', 'score', 'delta_score', 'condition', 'guess', 'guess_conf', 'guesser', 'respondent'])

    conditions = ['PL', 'AC']
    guesses = ['PL', 'AC']

    for condition, guess in itertools.product(conditions, guesses):

        for idx in range(300):

            if (condition == 'PL') and (guess == 'PL'):
                delta_score = random.gauss(mu=10, sigma=1)
            elif (condition == 'PL') and (guess == 'AC'):
                delta_score = random.gauss(mu=50, sigma=1)
            elif ((condition == 'AC') and (guess == 'PL')) or ((condition == 'AC') and (guess == 'AC')):
                continue
            else:
                assert False

            df = df.append({
                'trial': 'test',
                'subject_id': idx,
                'tp': 'wk8',
                'scale': 'tadaa',
                'score': 0,
                'delta_score': delta_score,
                'condition': condition,
                'guess': guess,
                'guesser': 'self',
                'respondent': 'self'},
                ignore_index=True)

        for idx in range(500):

            if ((condition == 'PL') and (guess == 'PL')) or ((condition == 'PL') and (guess == 'AC')):
                continue
            elif (condition == 'AC') and (guess == 'PL'):
                delta_score = random.gauss(mu=70, sigma=1)
            elif (condition == 'AC') and (guess == 'AC'):
                delta_score = random.gauss(mu=30, sigma=1)
            else:
                assert False

            df = df.append({
                'trial': 'test',
                'subject_id': idx,
                'tp': 'wk8',
                'scale': 'tadaa',
                'score': 0,
                'delta_score': delta_score,
                'condition': condition,
                'guess': guess,
                'guesser': 'self',
                'respondent': 'self'},
                ignore_index=True)

    df.to_csv(os.path.join(folders.fixtures,
              'pseudodata_get_bbc_from_processed222.csv'), index=False)


def gen_faux_get_model_stats1():


def gen_faux_bbc_1():

    df = pd.DataFrame(columns=[
        'trial', 'scale', 'respondent', 'guesser', 'cgr', 'cgr_sim_id', 'condition', 'guess', 'delta_score'])

    cgrs = np.linspace(0, 1, 15)
    guesser = 'self'
    respondent = 'self'

    condition = 'PL'
    guess = 'PL'

    for cgr, idx in itertools.product(cgrs, range(100)):

        delta_score = random.gauss(mu=10, sigma=1)

        df = df.append({
            'trial':,
            'scale':,
            'respondent':,
            'guesser':,
            'cgr':,
            'cgr_sim_id':,
            'condition':,
            'guess':,
            'delta_score':},
            ignore_index=True)

    df.to_csv(os.path.join(folders.fixtures,
              'pseudodata_get_model_stats1.csv'), index=False)


def gen_faux_cov_selector data():

    df1 = pd.DataFrame(columns=['delta_score', 'A', 'B'])
    # df2 = pd.DataFrame(columns=['delta_score', 'A', 'B'])
    # df3 = pd.DataFrame(columns=['delta_score', 'A', 'B'])
    # df4 = pd.DataFrame(columns=['delta_score', 'A', 'B'])

    for idx in range(2000):

        delta_score = random.gauss(50, 10)
        A = random.gauss(10,  5)

        df = df.append({
            'delta_score': delta_score+A,
            'A': A,
            'B': random.gauss(10,  5), },
            ignore_index=True)

    df.to_csv(os.path.join(folders.fixtures,
              'pseudodata_cov_selector1.csv'), index=False)
