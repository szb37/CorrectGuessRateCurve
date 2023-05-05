"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import numpy as np

cgrc_parameters = {
    'test': {  # Quick testing
        'cgr_values': [0.5],
        'n_cgrc_trials': 2, },
    # Use cgrA (=correct guess rate ADJUSTMENT) if only perfect blinding estimate is needed
    'cgrA_low': {
        'cgr_values': [0.5],
        'n_cgrc_trials': 32, },
    'cgrA_mid': {
        'cgr_values': [0.5],
        'n_cgrc_trials': 64, },
    'cgrA_high': {
        'cgr_values': [0.5],
        'n_cgrc_trials': 96, },
    # Use cgrC (=correct guess rate CURVE) when plotting the whole CGR curve
    'cgrC_low': {
        'cgr_values': np.linspace(0, 1, 13).tolist(),
        'n_cgrc_trials': 32, },
    'cgrC_mid': {
        'cgr_values': np.linspace(0, 1, 13).tolist(),
        'n_cgrc_trials': 64, },
    'cgrC_high': {
        'cgr_values': np.linspace(0, 1, 13).tolist(),
        'n_cgrc_trials': 96, },
}

sbmd_all = {'sbmd': ['PANAS', 'mood', 'creativity', 'WEMWB', 'SCS', 'CPS', 'energy', 'focus', 'temper', 'QIDS', 'STAIT']}
sbmd_acutes = {'sbmd': ['CPS', 'PANAS', 'mood', 'creativity', 'energy', 'focus', 'temper']}
sbmd_postacutes = {'sbmd': ['QIDS', 'WEMWB', 'STAIT']}
sbmd_plots = {'sbmd': ['PANAS', 'mood', 'creativity', 'energy', 'CPS', 'WEMWB']}
sbmd_test = {'sbmd': ['PANAS']}

save_csvs = True
save_figs = False
estimator = 'mean'

trial_cgrs = {
    'sbmd': 0.72,
    'mock': 0.7,
}
