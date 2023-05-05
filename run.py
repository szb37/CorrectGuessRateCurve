"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

The self-blinding MD (SBMD) trial data is included in this repo (data/01_trial_data/).
AEB model data, CGR curve data for either AEB or SBMD is not stored, only their relevant statistics.
These can be found here for both the AEB models and SBMD:
https://drive.google.com/drive/folders/1sQe1GQ97DbkIiw3YeOcXwA1KZIsp79zT?usp=share_link
"""

import src.cgrc.core as cgrc
import src.config as config
import src.folders as folders
import src.figures as figures
import src.miscs as miscs
import os


if False: # CGR adjustment for self-blinding MD trial data (Table 2)

    trial_name = 'sbmd'
    postfix = 'tmp'
    analysis_name = f'{trial_name}_{postfix}'

    cgrc.Controllers.run_cgrc_trial(
        trial_name=trial_name,
        postfix=postfix,
        trial_scales=config.sbmd_all,
        cgrc_param_set = 'cgrA_low',
    )

    df = cgrc.Controllers.get_cgrc_comparison_table_v1(
        trial_name = trial_name,
        analysis_name = analysis_name,
        trial_data_dir = folders.trial_data_dir,
        trial_stats_dir = os.path.join(folders.trial_stats_dir, analysis_name),
        cgrc_data_dir = os.path.join(folders.cgrc_data_dir, analysis_name),
        cgrc_stats_dir = os.path.join(folders.cgrc_stats_dir, analysis_name),
        trial_scales = config.sbmd_all
    )

if False: # CGR curves for self-blinding MD trial data (Figure 4)

    trial_name = 'sbmd'
    postfix = 'tmp'
    analysis_name = f'{trial_name}_{postfix}'

    cgrc.Controllers.run_cgrc_trial(
        trial_name=trial_name,
        postfix=postfix,
        trial_scales=config.sbmd_plots,
        cgrc_param_set = 'cgrC_low',
        save_figs = True,
    )

if False: # CGR adjustment for the AEB models (Table 1)

    trial_name = 'defaults'
    postfix = 'cgrA'
    analysis_name = f'{trial_name}_{postfix}'

    summary_df = miscs.get_aeb_summary_table(
        analysis_name = analysis_name,
    )
    print('\n', summary_df.to_string(index=False))
    summary_df.to_csv('agyfasz', index=False)
