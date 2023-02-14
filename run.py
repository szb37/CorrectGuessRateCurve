"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

Functions that generate CGRC are not optimized, may take a while to run

GitHub repo only contains source data, all simulation data (and results) used in the paper can be accessed from:
https://drive.google.com/drive/folders/1sQe1GQ97DbkIiw3YeOcXwA1KZIsp79zT?usp=share_link
(copy drive folder to codebase/data/)
"""

import src.toy_models.model_defs as model_defs
import src.toy_models.core as toy_models
import src.config as config
import src.cgrc.core as cgrc
import time

start = time.time()

if True: # speed test, ~3.8s
    toy_models.Controllers.run_toymodels_cgrc(
        analysis_name = 'test',
        models = model_defs.test,
        cgrc_param_set = 'cgrA_low',
        n_patients = 50,
        n_trials = 2
    )

if False: # santy check
    toy_models.Controllers.run_toymodels_cgrc(
        analysis_name = 'tmp',
        models = model_defs.tmp,
        cgrc_param_set = 'cgrA_mid',
        n_patients = 200,
        n_trials = 100
    )


end = time.time()
print("Execution time was:", round(end-start, 3),'s')
