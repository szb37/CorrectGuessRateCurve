"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import os


src = os.path.dirname(os.path.abspath(__file__))
codebase = os.path.abspath(os.path.join(src, os.pardir))

data = os.path.abspath(os.path.join(codebase, 'data'))
data_trial_data = os.path.abspath(os.path.join(data, '01_trial_data'))
data_trial_stats = os.path.abspath(os.path.join(data, '02_trial_stats'))
data_cgrc_data = os.path.abspath(os.path.join(data, '03_cgrc_data'))
data_cgrc_stats = os.path.abspath(os.path.join(data, '04_cgrc_stats'))
data_summary_tables = os.path.abspath(os.path.join(data, '05_summary_tables'))

export_plot = os.path.abspath(os.path.join(codebase, 'export figs'))
strata_plots = os.path.abspath(os.path.join(export_plot, 'stratification'))
cgrc_plots = os.path.abspath(os.path.join(export_plot, 'CGR_curves'))

test = os.path.abspath(os.path.join(codebase, 'tests'))
fixtures = os.path.abspath(os.path.join(test, 'fixtures'))
tmp_dir = os.path.join(fixtures, 'tmp')
