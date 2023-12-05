import chemist
import itertools

element_symbols = [
    'X', 'H', 'He', 
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe'
]

atomic_masses = [
    0.0,     # Dummy
    1.0079,  # Hydrogen
    4.0026,  # Helium
    6.94,    # Lithium
    9.0122,  # Beryllium
    10.81,   # Boron
    12.011,  # Carbon
    14.007,  # Nitrogen
    15.999,  # Oxygen
    18.998,  # Fluorine
    20.180,  # Neon
    22.990,  # Sodium
    24.305,  # Magnesium
    26.982,  # Aluminium
    28.085,  # Silicon
    30.974,  # Phosphorus
    32.06,   # Sulfur
    35.45,   # Chlorine
    39.948,  # Argon
    39.098,  # Potassium
    40.078,  # Calcium
    44.956,  # Scandium
    47.867,  # Titanium
    50.942,  # Vanadium
    51.996,  # Chromium
    54.938,  # Manganese
    55.845,  # Iron
    58.933,  # Cobalt
    58.693,  # Nickel
    63.546,  # Copper
    65.38,   # Zinc
    69.723,  # Gallium
    72.630,  # Germanium
    74.922,  # Arsenic
    78.971,  # Selenium
    79.904,  # Bromine
    83.798,  # Krypton
    85.468,  # Rubidium
    87.62,   # Strontium
    88.906,  # Yttrium
    91.224,  # Zirconium
    92.906,  # Niobium
    95.95,   # Molybdenum
    97.907,  # Technetium
    101.07,  # Ruthenium
    102.91,  # Rhodium
    106.42,  # Palladium
    107.87,  # Silver
    112.41,  # Cadmium
    114.82,  # Indium
    118.71,  # Tin
    121.76,  # Antimony
    127.60,  # Tellurium
    126.90,  # Iodine
    131.29   # Xenon
]

def get_atomic_mass(x):
    """Return the atomic mass of an element or atom.

    Parameters
    ----------
    x : str or int
        Element symbol or atomic number

    Returns
    -------
    mass : float
        Atomic mass of the atom
    """
    mass = 0.0
    if type(x) == str:
        atomic_number = get_atomic_number(x)
        mass = atomic_masses[atomic_number]
    elif type(x) == int:
        mass = atomic_masses[x]
    else:
        print("Error: Invalid input type")
    return mass


def get_atomic_number(symbol):
    """Return the atomic number for a given element symbol.

    Parameters
    ----------
    symbol : str
        Element symbol (case insensitive)
    
    Returns
    -------
    atomic_number : int
        Atomic number of the element
    >>> print(get_atomic_number('H'))
    1
    """
    symbol = symbol.capitalize()
    return element_symbols.index(symbol)


def get_molecule(symbols, coordinates):
    """Return a chemist.Molecule object from a list of symbols and coordinates.

    Parameters
    ----------
    symbols : list of str
        List of element symbols
    coordinates : list of lists of floats
        List of lists containing the coordinates of the atoms in the molecule.
        The inner lists contain the x, y, and z coordinates of each atom.

    Returns
    -------
    mol : chemist.Molecule
        chemist.Molecule object
    """
    natom = len(symbols)
    mol = chemist.Molecule()
    if type(coordinates[0]) == list or str(type(coordinates[0])).startswith('numpy'):
        coordinates = list(itertools.chain(*coordinates))
    for i in range(natom):
        x, y, z = [float(x) for x in coordinates[3 * i:3 * i + 3]]
        symbol = symbols[i]
        atomic_number = get_atomic_number(symbol)
        mass = get_atomic_mass(atomic_number)
        atom = chemist.Atom(symbol, atomic_number, mass, x, y, z)
        mol.push_back(atom)
    return mol


def get_molecule_coordinates(mol):
    """Return the coordinates of a molecule as a list of lists.

    Parameters
    ----------
    mol : chemist.Molecule
        chemist.Molecule object

    Returns
    -------
    coordinates : list of lists
        List of lists containing the coordinates of the atoms in the molecule.
        The inner lists contain the x, y, and z coordinates of each atom.
    """
    return [[mol.at(i).x, mol.at(i).y, mol.at(i).z] for i in range(mol.size())]


def get_molecule_symbols(mol):
    """Return the symbols of a molecule as a list.

    Parameters
    ----------
    mol : chemist.Molecule
        chemist.Molecule object
    
    Returns
    -------
    symbols : list of str
        List of element symbols of the atoms in the molecule.
    """
    natom = mol.size()
    symbols = []
    for i in range(natom):
        symbols.append(get_symbol(mol.at(i).Z))
    return symbols


def get_periodic_table():
    """Return the periodic table as a list.
    Returns
    -------
    element_symbols : list of str
    
    Notes
    -----
    Up to elements with atomic number 54 (Xenon).
    """
    return element_symbols


def get_symbol(atomic_number):
    """Returns the element symbol for a given atomic number.

    Parameters
    ----------
    atomic_number : int
        Atomic number of the element
    
    Returns
    -------
    symbol : str
        Element symbol
    
    Notes
    -----
    Returns 'X' for atomic_number=0
    >>> print(get_symbol(1))
    H
    """
    return element_symbols[atomic_number]

