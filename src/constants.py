"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import numpy as np


cgrc_parameters = {
    0: {  # use for testing
        'cgr_values': [0.5],
        'n_cgrc_trials': 4, },
    1: {  # use when only perfect blinding adjustment is needed; corresponds to cgrA (CGR adjustment; calc only CGR=0.5)
        'cgr_values': [0.5],
        'n_cgrc_trials': 4, },
    2: {  # use for full analysis; corresponds to cgrC (CGR curve calc range of CGR values)
        'cgr_values': np.linspace(0, 1, 13).tolist(),
        'n_cgrc_trials': 4, },
}

sbmd_all = {'sbmd': ['PANAS', 'mood', 'creativity', 'WEMWB', 'SCS', 'CPS', 'energy', 'focus', 'temper', 'QIDS', 'STAIT']}
sbmd_acutes = {'sbmd': ['CPS', 'PANAS', 'mood', 'creativity', 'energy', 'focus', 'temper']}
sbmd_postacutes = {'sbmd': ['QIDS', 'WEMWB', 'STAIT']}
sbmd_plots = {'sbmd': ['PANAS', 'mood', 'creativity', 'energy', 'CPS', 'WEMWB']}
sbmd_test = {'sbmd': ['PANAS']}
sbmd_tmp = {'sbmd': ['creativity']}


trial_cgrs = {
    'sbmd': 0.72,
    'mock': 0.7,
}

estimator = 'mean'
