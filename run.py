"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

Functions that generate CGRC are not optimized, may take a while to run

GitHub repo only contains source data, all simulation data (and results) used in the paper can be accessed from:
https://drive.google.com/drive/folders/1sQe1GQ97DbkIiw3YeOcXwA1KZIsp79zT?usp=share_link
(replace content of codebase/data/ with content of shared drive folder)
"""

import src.toy_models.core as toy_models
import src.constants as constants
import src.cgrc.core as cgrc
import time

start = time.time()

# CGR adjusted analysis of all toy models (including robustness analysis) - reproduces Table 1 and Supp Table 1
if True:
    toy_models.Controllers.run_cgrc_model_family(
        model_family_name='test_model',
        postfix='tmp',
        cgrc_param_set=1,
        save_figs=False,
        n_trials=8,
        n_patients=80,
    )

if False:  # CGR curves of default toy models - reproduces Figure 3
    toy_models.Controllers.run_cgrc_model_family(
        model_family_name='default_models',
        postfix='cgrC',
        cgrc_param_set=2,
        save_figs=True,
        n_trials=1,
        n_patients=100,
    )

if False:  # CGR adjusted outcomes of self-blinding microdose trial - reproduces Table 2
    cgrc.Controllers.run_cgrc_trial(
        trial_name='sbmd',
        postfix='cgrA_tmp',
        trial_scales=constants.sbmd_all,
        cgrc_param_set=1,
        save_figs=False,
    )

if False:  # CGR curves of self-blinding microdose trial outcomes - reproduces Figure 4
    cgrc.Controllers.run_cgrc_trial(
        trial_name='sbmd',
        postfix='cgrC_tmp',
        trial_scales=constants.sbmd_plots,
        cgrc_param_set=2,
        save_figs=True,
    )

end = time.time()

print("Execution time was:", (end-start)*10**3, "ms")
