from nwchemex import chemist
from berny import Berny, geomlib, optimize
from berny.solvers import MopacSolver

# Molecule for testing
#h2 = chemist.Molecule()
#a1       = chemist.Atom("H", 1, 1.0079, 0.0, 0.0, 0.0)
#a2       = chemist.Atom("H", 1, 1.0079, 0.0, 0.0, 1.0)
#h2.push_back(a1)
#h2.push_back(a2)

def optimize_pyberny(molecule):
    """
    This function takes a chemical system and uses PyBerny
    to optimize the provided geometry (molecule) and 
    returns the energy of the optimized system.
    """
    # Convert a Chemical System to an XYZ coordinate string
    xyz = str(molecule.size()) + "\n\n" + str(molecule.nuclei)
    
    # Loads the geometry string into the Berny optimizer
    # object.
    optimizer = Berny(geomlib.loads(xyz, fmt='xyz'))
    
    # Uses MOPAC/OpenMOPAC to perform geometry optimization.
    # The program "mopac" must be in system PATH
    solver = MopacSolver()
    next(solver)
    for geom in optimizer:
        energy, gradients = solver.send((list(geom), geom.lattice))
        optimizer.send((energy, gradients))

    relaxed = geom
    xyz_opt = relaxed.dumps(fmt='xyz')
    # Optimized energy is of type "float"
    return energy, xyz_opt
