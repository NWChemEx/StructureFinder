from berny import Berny, geomlib


def optimize_pyberny(molecule):
    """
    This function takes a chemical system and uses PyBerny
    to optimize the provided geometry (molecule) and
    returns the energy of the optimized system.
    """
    # Convert a Chemical System to an XYZ coordinate string
    xyz = ""
    xyz += (str(molecule.size()) + "\n\n")
    for i in range(molecule.size()):
        xyz += (molecule.at(i).name + " " + str(molecule.at(i).x) + " " +
                str(molecule.at(i).y) + " " + str(molecule.at(i).z) + "\n")

    # Loads the geometry string into the Berny optimizer
    # object.
    optimizer = Berny(geomlib.loads(xyz, fmt='xyz'))

    for geom in optimizer:
        energy = calculate_energy(geom)
        gradients = calculate_gradient(geom)
        optimizer.send((energy, gradients))

    relaxed = geom
    xyz_opt = relaxed.dumps(fmt='xyz')

    # Optimized energy is of type "float"
    return energy, xyz_opt
