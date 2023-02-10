"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

name: model's name
p_act: probability of AC treatment
p_sea: probability of side effects in AC group; i.e. guessing AC when treatment is AC
p_sep: probability of side effects in PL group; i.e. guessing AC when treatment is PL
nhist: mean and std of the outcome's natural history
dte: mean and std of the direct treatment effect
aeb: mean and std of the activated expectancy bias
"""


class ModelDefinition(dict):
    ''' A dictionary defining toy models '''

    def __init__(self, name=None, p_act=0.5, p_sea=0.9, p_sep=0.5, nhist=(20,4), dte=(10,2), aeb=(8,2)):
        super(ModelDefinition, self).__init__()
        self['name']=name
        self['p_act']=float(p_act)
        self['p_sea']=float(p_sea)
        self['p_sep']=float(p_sep)
        self['nhist']=(float(nhist[0]), float(nhist[1]))
        self['dte']=(float(dte[0]), float(dte[1]))
        self['aeb']=(float(aeb[0]), float(aeb[1]))
        assert self.is_valid()

    def is_valid(self):
        assert (isinstance(self['name'], str) or self['name'] is None)

        assert isinstance(self['p_act'], float)
        assert (0 <= self['p_act'] <=1)

        assert isinstance(self['p_sea'], float)
        assert (0 <= self['p_sea'] <=1)

        assert isinstance(self['p_sep'], float)
        assert (0 <= self['p_sep'] <=1)

        assert isinstance(self['nhist'], tuple)
        assert isinstance(self['nhist'][0], float)
        assert isinstance(self['nhist'][1], float)

        assert isinstance(self['dte'], tuple)
        assert isinstance(self['dte'][0], float)
        assert isinstance(self['dte'][1], float)

        assert isinstance(self['aeb'], tuple)
        assert isinstance(self['aeb'][0], float)
        assert isinstance(self['aeb'][1], float)

        return True


# Define models to be simulated
off_off_0 = ModelDefinition(
    name='off_off_0',
    dte=(0, 0),
    aeb=(0, 0),
    )
on_on_0 = ModelDefinition(
    name='on_on_0',
    )

on_off_10 = ModelDefinition(
    name='on_off_10',
    dte=(10, 2),
    aeb=(0, 0),
    )
on_off_7 = ModelDefinition(
    name='on_off_7',
    dte=(7, 2),
    aeb=(0, 0),
    )
on_off_5 = ModelDefinition(
    name='on_off_5',
    dte=(5, 2),
    aeb=(0, 0),
    )




# Define usefull lists of models
speed_test = [off_off_0]
models710 = [on_off_7, on_off_10]
models5 = [on_off_5]



test_models = [off_off_0]
default_models = []
all_models = []
