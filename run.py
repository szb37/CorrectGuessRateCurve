"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import src.toy_models.model_defs as model_defs
import src.toy_models.core as toy_models
import src.config as config
import src.cgrc.core as cgrc
import time

start = time.time()

if False: # speed test, ~3.8s
    toy_models.Controllers.run_toymodels_cgrc(
        analysis_name = 'test',
        models = model_defs.test,
        cgrc_param_set = 'cgrA_low',
        n_patients = 50,
        n_trials = 2
    )

if True: # run default models
    toy_models.Controllers.run_toymodels_cgrc(
        analysis_name = 'default_models',
        models = model_defs.default_models,
        cgrc_param_set = 'cgrA_mid',
        n_patients = 120,
        n_trials = 100
    )


end = time.time()
print("Execution time was:", round(end-start, 3),'s')
