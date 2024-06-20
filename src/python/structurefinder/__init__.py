import pluginplay as pp
from simde import TotalEnergy

def __init__(self):
    pp.ModuleBase.__init__(self)
    self.description("Geometry Optimizer Plugin")
    self.satisfies_property_type(TotalEnergy())
    self.add_input('molecule')
