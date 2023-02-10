"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

Functions that generate CGRC are not optimized, may take a while to run

GitHub repo only contains source data, all simulation data (and results) used in the paper can be accessed from:
https://drive.google.com/drive/folders/1sQe1GQ97DbkIiw3YeOcXwA1KZIsp79zT?usp=share_link
(replace content of codebase/data/ with content of shared drive folder)
"""

import src.toy_models.model_defs as model_defs
import src.toy_models.core as toy_models
import src.config as config
import src.cgrc.core as cgrc
import time

start = time.time()

if False: # speed test, ~3.5s
    toy_models.Controllers.run_toymodels_cgrc(
        analysis_name = 'tmp',
        models = model_defs.speed_test,
        cgrc_param_set = 'cgrA_low',
        n_patients = 50,
        n_trials = 2
    )


if True: # santy check
    toy_models.Controllers.run_toymodels_cgrc(
        analysis_name = 'tmp',
        models = model_defs.models5,
        cgrc_param_set = 'cgrA_low',
        n_patients = 50,
        n_trials = 10
    )



end = time.time()
print("Execution time was:", round(end-start, 3),'s')
