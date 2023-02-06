"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2020, DrugNerdsLab
:License: MIT

The following parameters are all random normal variables defining the toy models (mean, std):
    oc_nh: Outcomes's natural history; score space
    gs_nh: Guess's natural history; probability space
    se: Treatment's contribution to guess; probability space
    dte: Treatment's contribution to outcomes; score space
    pte: Placebo guess's contribution to guess; score space
    ate: Active guess's contribution to guess; score space
    dte: Treatment's contribution to outcomes; score space
    gs2oc: Guess's contribution to outcome; score space
    oc2gs: Outcome's contribution to guess; probability space

Model families: off_off, on_off, off_on, on_on




mu_TRT, sigma_TRT, mu_SE, sigma_SE

    oc_nh: Outcomes's natural history; score space
    gs_nh: Guess's natural history; probability space
    se: Treatment's contribution to guess; probability space
    dte: Treatment's contribution to outcomes; score space
    pte: Placebo guess's contribution to guess; score space
    ate: Active guess's contribution to guess; score space
    dte: Treatment's contribution to outcomes; score space
    gs2oc: Guess's contribution to outcome; score space
    oc2gs: Outcome's contribution to guess; probability space





"""


""" Models with direct drug effect OFF and placebo leak OFF """
off_off_0 = {
    'off_off_0': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0, 0),
        'dte': (0, 0),
        'pte': (0, 0),
        'ate': (0, 0),
        'oc2gs': (0, 0),
    },
}

""" Models with direct drug effect ON and placebo leak OFF """
on_off_m3 = {
    'on_off_m3': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0, 0),
        'dte': (7, 2),
        'pte': (0, 0),
        'ate': (0, 0),
        'oc2gs': (0, 0),
    },
}
on_off_m2 = {
    'on_off_m2': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0, 0),
        'dte': (8, 2),
        'pte': (0, 0),
        'ate': (0, 0),
        'oc2gs': (0, 0),
    },
}
on_off_m1 = {
    'on_off_m1': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0, 0),
        'dte': (9, 2),
        'pte': (0, 0),
        'ate': (0, 0),
        'oc2gs': (0, 0),
    },
}
on_off_0 = {
    'on_off_0': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0, 0),
        'dte': (10, 2),
        'pte': (0, 0),
        'ate': (0, 0),
        'oc2gs': (0, 0),
    },
}
on_off_p1 = {
    'on_off_p1': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0, 0),
        'dte': (11, 2),
        'pte': (0, 0),
        'ate': (0, 0),
        'oc2gs': (0, 0),
    },
}
on_off_p2 = {
    'on_off_p2': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0, 0),
        'dte': (12, 2),
        'pte': (0, 0),
        'ate': (0, 0),
        'oc2gs': (0, 0),
    },
}
on_off_p3 = {
    'on_off_p3': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0, 0),
        'dte': (13, 2),
        'pte': (0, 0),
        'ate': (0, 0),
        'oc2gs': (0, 0),
    },
}

""" Models with direct drug effect OFF and placebo leak ON """
off_on_m3 = {
    'off_on_m3': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (0, 0),
        'pte': (0, 0),
        'ate': (5, 2),
        'oc2gs': (0, 0),
    },
}
off_on_m2 = {
    'off_on_m2': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (0, 0),
        'pte': (0, 0),
        'ate': (6, 2),
        'oc2gs': (0, 0),
    },
}
off_on_m1 = {
    'off_on_m1': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (0, 0),
        'pte': (0, 0),
        'ate': (7, 2),
        'oc2gs': (0, 0),
    },
}
off_on_0 = {
    'off_on_0': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (0, 0),
        'pte': (0, 0),
        'ate': (8, 2),
        'oc2gs': (0, 0),
    },
}
off_on_p1 = {
    'off_on_p1': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (0, 0),
        'pte': (0, 0),
        'ate': (9, 2),
        'oc2gs': (0, 0),
    },
}
off_on_p2 = {
    'off_on_p2': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (0, 0),
        'pte': (0, 0),
        'ate': (10, 2),
        'oc2gs': (0, 0),
    },
}
off_on_p3 = {
    'off_on_p3': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (0, 0),
        'pte': (0, 0),
        'ate': (11, 2),
        'oc2gs': (0, 0),
    },
}

""" Models with direct drug effect ON and placebo leak ON """
on_on_m3 = {
    'on_on_m3': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (7, 2),
        'pte': (0, 0),
        'ate': (5, 2),
        'oc2gs': (0, 0),
    },
}
on_on_m2 = {
    'on_on_m2': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (8, 2),
        'pte': (0, 0),
        'ate': (6, 2),
        'oc2gs': (0, 0),
    },
}
on_on_m1 = {
    'on_on_m1': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (9, 2),
        'pte': (0, 0),
        'ate': (7, 2),
        'oc2gs': (0, 0),
    },
}
on_on_0 = {
    'on_on_0': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (10, 2),
        'pte': (0, 0),
        'ate': (8, 2),
        'oc2gs': (0, 0),
    },
}
on_on_p1 = {
    'on_on_p1': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (11, 2),
        'pte': (0, 0),
        'ate': (9, 2),
        'oc2gs': (0, 0),
    },
}
on_on_p2 = {
    'on_on_p2': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (12, 2),
        'pte': (0, 0),
        'ate': (10, 2),
        'oc2gs': (0, 0),
    },
}
on_on_p3 = {
    'on_on_p3': {
        'oc_nh': (20, 4),
        'gs_nh': (0.5, 0.1),
        'se': (0.135, 0.03),
        'dte': (13, 2),
        'pte': (0, 0),
        'ate': (11, 2),
        'oc2gs': (0, 0),
    },
}

""" Model families (=lists of model specifications) """
test_model = [
    off_off_0,
]

default_models = [
    off_off_0,
    on_off_0,
    off_on_0,
    on_on_0,
]

all_models = [
    on_off_m3,
    off_on_m3,
    on_on_m3,

    on_off_m2,
    off_on_m2,
    on_on_m2,

    on_off_m1,
    off_on_m1,
    on_on_m1,

    off_off_0,
    on_off_0,
    off_on_0,
    on_on_0,

    on_off_p1,
    off_on_p1,
    on_on_p1,

    on_off_p2,
    off_on_p2,
    on_on_p2,

    on_off_p3,
    off_on_p3,
    on_on_p3,
]
