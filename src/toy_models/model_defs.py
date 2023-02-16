"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT

name: model's name
p_act: probability of AC treatment
p_unb: probability of unblinding; i.e. prob of correct treatment guess
nhist: mean and std of the outcome's natural history
dte: mean and std of the direct treatment effect
aeb: mean and std of the activated expectancy bias
"""

class ModelDefinition(dict):
    ''' A dictionary defining toy models '''

    def __init__(self, name=None, p_act=0.5, p_unb=0.7, nhist=(10,2), dte=(10,2), aeb=(10,2)):
        super(ModelDefinition, self).__init__()
        self['name']=name
        self['p_act']=float(p_act)
        self['p_unb']=float(p_unb)
        self['nhist']=(float(nhist[0]), float(nhist[1]))
        self['dte']=(float(dte[0]), float(dte[1]))
        self['aeb']=(float(aeb[0]), float(aeb[1]))
        assert self.is_valid()

    def is_valid(self):
        assert (isinstance(self['name'], str) or self['name'] is None)

        assert isinstance(self['p_act'], float)
        assert (0 <= self['p_act'] <=1)

        assert isinstance(self['p_unb'], float)
        assert (0 <= self['p_unb'] <=1)

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
dte0aeb0 = ModelDefinition(
    name='dte0aeb0',
    dte=(0, 0),
    aeb=(0, 0),)
dte10aeb0 = ModelDefinition(
    name='dte10aeb0',
    dte=(10, 2),
    aeb=(0, 0),)
dte0aeb8 = ModelDefinition(
    name='dte0aeb8',
    dte=(0, 0),
    aeb=(8, 2),)
dte10aeb8 = ModelDefinition(
    name='dte10aeb8',
    dte=(10, 2),
    aeb=(8, 2),)


# Define usefull lists of models
test = [dte0aeb0]
test_models = [dte0aeb0, dte10aeb8]
default_models = [dte0aeb0, dte10aeb0, dte0aeb8, dte10aeb8]
