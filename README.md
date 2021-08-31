# LigandScatter

## Purpose

The purpose of ligand scatter is to generate gaussian input files where a specified ligand is randomly placed in a
sphere around a certain moiety of a desired molecule. Because the position *and* orientation of the ligand is randomly
chosen on each placement, if a large number of files are generated, we can assume that we have statistically all
possible ligand approaches to the moiety.

## Example

Suppose we want to optimize the interactions between the alcohol group on cyclopentanol, and water.

I have two files, one named cyclopentanol.xyz, and the other named water.xyz

![cyclopentanol](https://puu.sh/I7rsW/b40f019ced.png)

![water](https://puu.sh/I7rrW/3296773133.png)

I run the executable and follow the instructions. The first step is to input the directory which contains the base
molecule. In this case, we're looking for the folder which holds the cyclopentanol file.

![step1.1](https://puu.sh/I7rU6/f4e1fe1516.png)

![step1.2](https://puu.sh/I7rVk/38ad5729b9.png)

When prompted, give the number which represents the file we want to use. In this case, 2.

![step1.3](https://puu.sh/I7rVO/5072cabffc.png)

Repeat this process for the ligand. In this case, the water molecule.

![step2](https://puu.sh/I7rX9/4a6b404de8.png)

We now have the option to view the base molecule. We will need to know the indices of the alcohol moiety, so we will
select yes.

![step3.1](https://puu.sh/I7rYT/650232e50d.png)

Once you press enter, the base molecule will be shown:

![step3.2](https://puu.sh/I7rZf/3d51624eb4.png)

The screen that pops up is interactable, so you can drag, zoom, etc.

![step3.3](https://puu.sh/I7rZF/39fd44dcf2.png)

We can see here that the alcohol consists of atoms 13 and 14. We will need to remember these numbers.

Once we have the numbers, we need to close the figure. When prompted, we type in the numbers representing the alcohol
moiety as indicated:

![step4.1](https://puu.sh/I7s2J/492c770333.png)

The center of the moiety is defined as the weighted center of 'mass' of the atoms in the moiety. The weight for each
atom is equal to that atom's covalent/ionic radius. The ligand will be scattered in a sphere around this center of mass.
The next step is to define the radius of that sphere in angstroms (Technically this number is unitless, but most file
formats give position data in angstroms). In this case, I will use 3 angstroms.

![step5.1](https://puu.sh/I7s7M/66084fd799.png)

Next we give the number of ligands we want to add. This is not the number of files we want to make, this is the number
of ligands we will be adding to the molecule for each file. In this case, lets go with 2.

![step6.1](https://puu.sh/I7s8J/94771369bb.png)

Next we give the number of files we want to end up with. This number should be large, so I'm going with 10,000.

![step7.1](https://puu.sh/I7s9Y/f706b04fdb.png)

Next we name the folder we want the files to be written to. This folder will be made in the same location as the
executable / script.

![step8.1](https://puu.sh/I7saS/532a7f0e95.png)

Now we are prompted to change any of the output file settings. I'm going to change the job to Opt, instead of Opt Freq.
I type out 'n' to not use the default settings, then 'job' to edit the job.

![step9.1](https://puu.sh/I7sHE/0baf5d89b4.png)

![step9.2](https://puu.sh/I7sHW/3bd1e4810e.png)

Since I don't want to change anything else, I type exit to save the settings.

![step9.3](https://puu.sh/I7sIh/54ed823b00.png)

Once we press enter, the program will:

1. Sample a random location from the moiety.
2. Randomly orient a ligand
3. Add the ligand to the base molecule
4. Repeat steps 1-3 for the desired number of ligands
5. Check that no molecules are too close together (determined by covalent/ionic radius)
    1. If this check fails, the program will start from scratch.
6. Write the molecule and settings to a .com file in the specified output folder
7. Repeat steps 1-7 till the specified file count is reached.

Once the program finishes, it will tell you so, and you are safe to close the window and view your files.


![video_of_output](https://youtu.be/VtDmOqaP7BQ)
