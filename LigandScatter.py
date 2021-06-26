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


def main():
    """Main Function"""

    base_prompt = "Enter the directory which contains the base file."
    ligand_prompt = "Enter the directory which contains the ligand file."

    base_file = select_file(base_prompt)
    ligand_file = select_file(ligand_prompt)

    # Make the base and ligand compounds
    base_xyz = parsing_dict[base_file[base_file.index(".") + 1:]](base_file)
    base_atoms = [Atom(a[0], (a[1], a[2], a[3])) for a in base_xyz]
    molecule_name = input("\nWhat is the name of the base compound: ")
    base_compound = Molecule(molecule_name, base_atoms)
    base = center_on_atom(base_compound, 0)

    ligand_xyz = parsing_dict[ligand_file[ligand_file.index(".") + 1:]](ligand_file)
    ligand_atoms = [Atom(a[0], (a[1], a[2], a[3])) for a in ligand_xyz]
    molecule_name = input("\nWhat is the name of the ligand compound: ")
    ligand_compound = Molecule(molecule_name, ligand_atoms)
    ligand = center_on_atom(ligand_compound, 0)

    if yes_no("Show the base molecule"):
        show_structure(base)

    # Get user input
    selection = input("What atoms make up the moiety (ex. '1,4,2'): ").replace(" ", "")
    center = center_of_mass(base, list(map(int, selection.split(","))))

    rad = verified_input("What is the radius of the moiety (Å): ", float)
    ligand_count = verified_input("How many ligands would you like to add: ", int)
    file_count = verified_input("How many files do you want to end up with: ", int)

    moiety = Sphere(pos=center, radius=rad)

    # Perform ligand additions and writing
    for i in tqdm(range(file_count), desc="Adding Ligands and Saving File"):
        mo = deepcopy(base)
        for _ in range(ligand_count):
            spun_ligand = randomly_orient(ligand)
            ligand_pos = moiety.sample()
            mo.add_molecule(spun_ligand, ligand_pos)
        write_job_to_com(mo, title=f"{mo.name}_{ligand_count}{ligand.name}_{i + 1}")


if __name__ == "__main__":
    print("Ligand Scatter".center(50, "~"))
    print("Author: Jackson Elowitt")
    print("Repo: https://github.com/jkelowitt/LigandScatter")
    print("Version: v1")
    print("".center(50, "~"))

    main()

    input("\nCalculating and saving complete. Press enter to close. ")
