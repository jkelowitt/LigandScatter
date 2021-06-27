"""
@Author: Jackson Elowitt
@Start Date: 6/13/2021
@Contact: jkelowitt@protonmail.com
@Site: https://github.com/jkelowitt/LigandScatter

Generates a set of files where ligands are placed at random positions around a molecule
within a sphere centered at a user defined moiety.
"""
from copy import deepcopy
from glob import glob

from numpy import *
from tqdm import tqdm

from classes import *
from functions import *
from parsing import *

w = [Atom("O", (0, 0, 0)),
     Atom("H", (1, 0, 0)),
     Atom("H", (0, 1, 0))]

water = Molecule("Water", w)

e = [Atom("C", (-1.50610, 0.76158, -0.10063)),
     Atom("C", (-0.37068, 1.77049, -0.18112)),
     Atom("H", (-1.89961, 0.71316, 0.93696)),
     Atom("H", (-1.13618, -0.24392, -0.39267)),
     Atom("H", (-2.32724, 1.05681, -0.78791)),
     Atom("H", (0.00934, 1.79810, -1.22659)),
     Atom("H", (0.44016, 1.45189, 0.51114)),
     Atom("O", (-0.84407, 3.03504, 0.18818)),
     Atom("H", (-0.06799, 3.64958, 0.11823))]

ethanol = Molecule("Ethanol", e)


# show_structure(water)


def center_of_mass(molecule: Molecule, atoms: list[int], weighted=True):
    """
    Determine the center of mass for a collection of atoms within a molecule

    Parameters
    ----------
    molecule: The molecule for analysis
    atoms: The list of atom indices to have the center found
    weighted: If true, the covalent radii of the atoms will be taken into account
        when determining the

    Example
    -------

    w = [Atom("O", (0, 0, 0)),
        Atom("H", (1, 0, 0)),
        Atom("H", (0, 1, 0))]

    water = Molecule("Water", w)
    """
    assert len(atoms) > 0, "You must have at least one atom in the collection in order to determine the center of mass"

    com = array([0, 0, 0])

    total_mass = 0 if weighted else 1
    for a in atoms:
        atom = molecule.atoms[a]
        com = add(com, atom.pos)

        # Perform covalent radius bias
        if weighted:
            com = multiply(com, atom.cov_radius)  # Required or unsafe casting error occurs
            total_mass += atom.cov_radius

    com /= total_mass

    return com


def select_file(prompt):
    """Allows the user to select a file from a selected directory"""
    files = []
    while not files:
        # Get the directory from the user
        print(f"\n{prompt}")
        print("Valid file types:")
        for ext in parsing_dict:
            print(f"\t.{ext}")
        directory = input("Directory: ")

        # Get all the files
        # which have a parsing function.
        for ext in parsing_dict:
            files += glob(f"{directory}/*.{ext}")

        files.sort(key=lambda x: len(x))

        # Check that there are valid files to be found.
        if not files:
            print("\nNo valid files found in the selected directory.")

    choice = make_choice_list(files)
    return choice


def make_molecule_from_file(file):
    """Given a file, return a Molecule object with the molecule contained in the file"""
    xyz = parsing_dict[file[file.index(".") + 1:]](file)
    name = file.split("\\")[-1]
    name = name[:name.index(".")]
    atoms = [Atom(a[0], (a[1], a[2], a[3])) for a in xyz]
    molecule = Molecule(name, atoms)
    molecule = center_on_atom(molecule, 0)
    return molecule


def main():
    """Main Function"""

    base_prompt = "Enter the directory which contains the base file."
    ligand_prompt = "Enter the directory which contains the ligand file."

    base_file = select_file(base_prompt)
    ligand_file = select_file(ligand_prompt)

    # Make the base and ligand compounds
    base = make_molecule_from_file(base_file)
    ligand = make_molecule_from_file(ligand_file)

    if yes_no("Show the base molecule"):
        show_structure(base)

    # Get user input
    selection = input("What atoms make up the moiety (ex. '1,4,2'): ").replace(" ", "")
    center = center_of_mass(base, list(map(int, selection.split(","))))

    rad = verified_input("What is the radius of the moiety (Å): ", float)
    ligand_count = verified_input("How many ligands would you like to add: ", int)
    file_count = verified_input("How many files do you want to end up with: ", int)
    output_dir = input("What would you like to name the destination folder: ")

    moiety = Sphere(pos=center, radius=rad)

    bond_check_mo = deepcopy(base)

    # Make a molecule which contains the same ligands as our final compound, but
    # guarantee that none of them overlap.
    for i in range(ligand_count):
        far_away = ((i + 1) * 100_000, (i + 1) * 100_000, (i + 1) * 100_000)
        bond_check_mo.add_molecule(ligand, global_pos=far_away)

    # Perform ligand additions and writing
    for i in tqdm(range(file_count), desc="Adding Ligands and Saving File"):
        error = True
        while error:
            # Deep copy so we get no mutability issues
            mo = deepcopy(base)

            # Add randomly oriented ligands
            for _ in range(ligand_count):
                spun_ligand = randomly_orient(ligand)
                ligand_pos = moiety.sample()
                mo.add_molecule(spun_ligand, ligand_pos)

            # Check that no ligands overlap
            error = check_bonds(bond_check_mo, mo) # 0 -> No errors -> break

        write_job_to_com(mo, title=f"{mo.name}_{ligand.name}_{i + 1}", output=output_dir)


if __name__ == "__main__":
    print("Ligand Scatter".center(50, "~"))
    print("Author: Jackson Elowitt")
    print("Repo: https://github.com/jkelowitt/LigandScatter")
    print("Version: v2")
    print("".center(50, "~"))

    main()

    input("\nFiles saved to folder in executable's location. Press enter to close. ")
