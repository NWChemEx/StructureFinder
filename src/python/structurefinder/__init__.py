# Copyright 2024 NWChemEx Community
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from .pyberny import load_pyberny_modules
from .lj_potential.lennard_jones_potential_module import load_lennard_jones_potential
from .fire.backward_euler.backward_euler import load_backwardeulerfire_modules

def load_modules(mm):
    """
    Loads the collection of all modules provided by StructureFinder.
    """
    load_pyberny_modules(mm)
    load_lennard_jones_potential(mm)
    load_backwardeulerfire_modules(mm)
