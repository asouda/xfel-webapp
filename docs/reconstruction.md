# Reconstruct the model from the omit mFo−DFc maps

This workflow rebuilds an omit model one heavy atom at a time. Each iteration refines the current omit model, calculates the weighted difference map \(mF_o-DF_c\), and finds the highest peak in the clipped map. The reference-model atom closest to that peak is inserted into the refined coordinates, at the peak rather than at the reference coordinate. That refined model, with the new atom, is the omit model of the next iteration.

The coefficient is the omit-map coefficient. The [omit-map page](omit-map.md) is the program sequence that produces `DELFWT` / `PHDELWT`, and the [theoretical background](theory.md) is why a missing atom is a positive peak, how the voxel is chosen, and how the model is updated between iterations. This page is the files, the loop, and the worked example.

## Worked example

The starting model is the single conformation B of the oxygen-evolving complex in the Ibrahim *et al.* two-flash dataset `6w1v`, with O6 removed. The reference model is the deposited `6w1v` structure, whose major site `OEY` A416 includes O6. The observations are the `6w1v` structure factors. This is the same comparison as Fig. 2 of Wang, Armstrong, and Batista, [Proc. Natl. Acad. Sci. U.S.A. 118, e2023982118 (2021)](https://doi.org/10.1073/pnas.2023982118), run repeatedly so that each omitted atom can be put back where the difference density is highest.

| Role | What the form needs | Download |
| --- | --- | --- |
| Starting omit coordinates | `6w1v_omitO6_single_confB.pdb`. Ligand `OEY`, chain `A`, residue **416**, with O6 absent. Both coordinate files must be PDB. | [6w1v_omitO6_single_confB.pdb](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v_omitO6_single_confB.pdb) |
| Reference model | The deposited `6w1v` coordinates. Bundle 1 contains the cluster, including O6 on `OEY` A416. | [bundle 1](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-pdb-bundle1.pdb), [bundle 2](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-pdb-bundle2.pdb) |
| Observed structure factors | `6w1v-sf.cif`, or an MTZ of the same reflections. | [6w1v-sf.cif.gz](https://github.com/asouda/xfel-webapp/raw/main/data/6w1v/6w1v-sf.cif.gz) |

The structure-factor file in the repository is gzip-compressed. Point the form at the `.gz` file. A name ending in `.gz` is decompressed with `gunzip`.

Set the site from the omit file: chain `A`, residue `416`, label `OEC`, and list D1-Glu189 (`A 189`) as a neighbour. Confirm the residue with **List hetero residues**. The deposited entry uses residue **415** for the minor conformation `OEX` and residue **416** for the major conformation `OEY`. The omit file on this form is only the major conformation, and O6 is the atom removed from it.

The default number of atoms is 1, the REFMAC cycle count is 0, and the distance and angle lists are empty. The default PNG contour is a single mesh at \(3.5\sigma\).

With the default choice, every heavy atom that is in the reference and absent from the omit model is a candidate, not only O6. The deposited `6w1v` model also contains the minor conformation `OEX` A415, which is absent from this omit file. The atom that is added is the missing heavy atom closest to the highest peak. A peak on the omitted O6 is expected to lie nearest O6 when that oxygen is closer to the peak than any atom of the other conformation.

## What one iteration does

The programs inside each iteration are the omit-map programs: `cif2mtz` when the reflections are mmCIF, `mtzdump`, `ctruncate` when the file holds intensities, `cad`, `refmac5`, `fft`, then the map padding, Coot, and PyMOL. `libcheck` is added only when distance restraints are widened in a monomer library. The commands are the ones on the [omit-map page](omit-map.md). The reconstruction adds the steps below, and it repeats them.

The input omit file is copied to `{omit model name}_working.pdb` in the output directory and is not modified. Each iteration copies that working file to `r-iteration-01`, `r-iteration-02`, and so on, as `{omit model name}_omit.pdb`, and runs the omit map in that subdirectory. Coordinate, map, CCP4, PNG, and density-section files in the iteration directory then receive `_step01`, `_step02`, and so on, before the extension. Logs, MTZ files, and the PyMOL script keep their names. The paths inside those scripts and logs are updated to the new file names.

1. **Calculate the current difference.** REFMAC runs for the number of cycles on the form. Zero cycles, the default, calculate \(mF_o-DF_c\) without moving atoms. A custom dictionary is read only when the cycle count is greater than zero. When external harmonic restraints are on, each atom added in an earlier iteration is held at the coordinates already written into the working model. The default sigma is 0.1 Å. Those restraints are written only when the cycle count is greater than zero and at least one atom has already been added.
2. **Pause, if a check was requested.** A map-check interval of \(n > 0\) pauses when the iteration number is a multiple of \(n\), after the map exists and before the peak is chosen. The viewer shows the current residual map and the current omit model. **Continue iterations** resumes the search. An interval of 0 does not pause. The final residual map is shown when the run finishes, whether or not checks were requested.
3. **Find the highest peak in the clip.** The search reads the Coot fragment, not the whole cell. A voxel is a candidate when its density equals the strongest finite value in a cube of the neighborhood size centered on that voxel. The neighborhood size is an odd integer; the default is 3. The selected peak is the candidate with the largest density. The reported density is the value on that voxel. The reported position is the Cartesian centre of the voxel.
4. **Optionally move the position onto a spline.** When the spline maximum is on, the position moves to the maximum of a spline interpolant of the chosen order inside that neighborhood. Orders run from 0 through 5. The default order is 3, a cubic spline, and the spline is off by default. Orders below 2 stay on the voxel. The density in the log remains the voxel value.
5. **Write the cross sections.** `density_sections_stepNN.txt` in the iteration directory samples the clip through the peak, from −2 Å to +2 Å in steps of 0.05 Å. The offset is the distance in ångströms along each line. The columns are the Cartesian \(x\), \(y\), and \(z\) axes, the six face diagonals of the cube, and the four body diagonals. A face diagonal is divided by \(\sqrt{2}\) and a body diagonal by \(\sqrt{3}\), so a step of 1 Å is 1 Å along the diagonal.
6. **Match the peak to a reference atom.** Candidates are reference atoms that are not hydrogen or deuterium and have occupancy greater than zero. By default the atom name must also be absent from the same chain, residue number, and insertion code in the working omit model. An atom name already present at any altloc is not a candidate in that mode. The form can instead make every heavy reference atom eligible, including atoms already in the omit model. A selected atom is still inserted as an additional record at the peak. Atoms added earlier in the run are excluded in either mode. The match is the candidate with the smallest Euclidean distance to the peak. A distance beyond 1 Å is a warning, and the atom is still added. A distance beyond 3 Å stops the run and the atom is not added. The run also stops when the clip has no finite peak, or when no eligible heavy atom remains.
7. **Insert the atom into the refined model.** The reference PDB record is copied into `{omit model name}_refined_stepNN.pdb`, the coordinates REFMAC wrote for this iteration. The atom name, residue, chain, altloc, occupancy, element, and displacement parameters stay those of the reference atom. The orthogonal \(x\), \(y\), and \(z\) are the peak coordinates. The result is `{omit model name}_refined_added_{atom name}_stepNN.pdb`. That file is the omit model of the next iteration, and it replaces `{omit model name}_working.pdb`. Every atom REFMAC moved stays where this refinement left it. The loop returns to step 1 until the requested number of atoms has been added.

When that count is reached, one further iteration refines the finished model and calculates its residual map. It does not search for another peak. A run that asks for one atom therefore writes two iteration directories: the search that adds the atom, and the refinement of the model that contains it. The same last iteration is what runs when the search stops early, because the clip has no finite peak, no eligible heavy atom remains, or the nearest atom is farther than 3 Å. The last iteration is the one that inserts no atom. Its REFMAC coordinate file, `{omit model name}_refined_stepNN.pdb`, is copied to `{omit model name}_reconstructed_final.pdb` in the output directory. `{omit model name}_working.pdb` remains the model that entered that refinement.

The figure for an iteration shows the refined coordinates of that iteration. The atom added in that iteration is not in the figure, because it is inserted after the map is drawn. It is in the omit model of the next iteration.

## The viewer

The title is **Reconstruction Step:** followed by the iteration number. The status line reads **Running step** and that number, or **Paused for inspection after step** and that number.

**Reference site** draws the reference coordinates in ball and stick, including atoms that are still absent from the omit model. Bonds are drawn from covalent radii: two heavy atoms are bonded when the distance is at least 0.40 Å, no longer than the sum of the covalent radii plus 0.45 Å, and no longer than 2.60 Å. A long metal–metal contact is not drawn. Hydrogens are not bonded. Different altlocs are not bonded to each other. The Mn/Ca–O name list used by the other difference workflows is not used here.

The four switches for Mn/Ca–O bonds, ray tracing, map recomputation, and fragment recomputation are not on this form. Every iteration recomputes the map and the clip and ray-traces the figure. The Mn1–O5 and Glu OE2 distance switch remains.

## Files

| Role | Form field | Required |
| --- | --- | --- |
| Starting omit coordinates | Starting omit-model coordinates | Yes. PDB or `.ent`. This is the model the first iteration reads. |
| Reference model | Reference model | Yes. PDB or `.ent`. It supplies the atoms that may be added, and it can be drawn in the viewer. It is not refined. |
| Observed structure factors | Observed structure factors, mmCIF or MTZ | Yes. \(F_o\) comes from this file. |

The two name fields are labels in the output files. They do not have to be PDB codes. Choosing a coordinate file fills the matching name from the file name, and fills the output directory so that the last path segment is the omit-model name. A directory you edit by hand is left as you edited it.

## Parameters that change this reconstruction

| Choice | Effect |
| --- | --- |
| Number of atoms to reconstruct | How many heavy atoms to insert. At least 1. The default is 1 |
| Atoms eligible for placement | Only heavy atoms absent from the omit model, the default, or every heavy atom in the reference model |
| Map check interval | `0` shows the viewer when the run finishes. A positive integer pauses after that many search iterations |
| Neighborhood size | Odd integer. Side of the cube in which a voxel must be the maximum. Default 3 |
| Spline maximum | Off by default. On moves the inserted coordinates to the spline maximum inside the neighborhood |
| Spline order | Integer from 0 through 5. Default 3. Used only when the spline maximum is on |
| REFMAC cycles | `0` leaves heavy atoms fixed. A larger integer refines the current omit model before each map. Either way, the next iteration starts from the refined coordinate file with the new atom inserted |
| External harmonic restraints | Off by default. On holds each atom added in an earlier iteration at the coordinates in the omit model of this iteration |
| Harmonic sigma (Å) | Width of that restraint. Default 0.1. Applied only when the cycle count is greater than zero |
| Distance and angle restraints | Same fields as the omit map. Empty by default on this form |
| Custom dictionary | Extra monomer library, under Parameters. Read only when the cycle count is greater than zero |
| Resolution source and the two cutoffs | Shell of the `DELFWT` coefficients. The explicit pair stored on the form, used only in explicit mode, is 30.0–1.89 Å |
| FFT grid resolution for a final map | Same meaning as on the omit map |
| PNG sigma levels | The worked example uses a single \(3.5\sigma\) mesh |

The 1 Å warning and the 3 Å stop are not form fields.

## What the output directory contains

Each iteration has its own subdirectory.

| Path | Contents |
| --- | --- |
| `r-iteration-01/`, `r-iteration-02/`, … | One omit-map calculation. PDB, CCP4, map, and PNG names in the directory end in `_stepNN`. Logs, MTZ files, and the PyMOL script do not |
| `r-iteration-NN/{omit model name}_omit_stepNN.pdb` | The omit model copied in at the start of the iteration, before REFMAC. The first iteration is the starting model. Each later iteration is the previous iteration’s refined coordinates with the atom added there |
| `r-iteration-NN/{omit model name}_refined_stepNN.pdb` | Coordinates after REFMAC for that iteration, before any atom is added |
| `r-iteration-NN/{omit model name}_refined_added_{atom name}_stepNN.pdb` | The refined coordinates with the atom inserted at the peak. Written only when an atom is added. This is the omit model of the next iteration |
| `r-iteration-NN/reconstruction.log` | The peak, the matched atom, the distance, and whether the atom was added |
| `r-iteration-NN/density_sections_stepNN.txt` | Cross sections through the peak. Written when a peak was found. The header records the voxel, the density, the peak coordinates, the atom, and the action |
| `{omit model name}_working.pdb` | The omit model that will enter the next iteration. After the last iteration, this is still the model that entered that refinement |
| `{omit model name}_reconstructed_final.pdb` | A copy of the last iteration’s refined coordinates, the iteration that inserts no atom |

The resolution tag in the map names replaces the decimal point with `p`, as on the omit map. A finer FFT grid adds `_grid` and that limit.

## Reading the result

Positive density is present in the observations and absent from the omit model of that iteration. After an atom is inserted at a peak, the next difference should be weaker there, because \(F_c\) now includes that atom. Density that remains is either another missing atom, a coordinate error, or noise.

A warning above 1 Å means the nearest missing reference atom does not sit on the peak. The atom is still added, at the peak. A stop above 3 Å means the highest peak is far from every remaining reference atom, so the run does not invent an identity for it.

The model in an iteration figure is the REFMAC coordinate set for that iteration, before the new atom is inserted. With zero cycles the heavy atoms are unmoved and REFMAC may still add riding hydrogens. The coordinates that enter the next iteration are `{omit model name}_refined_added_{atom name}_stepNN.pdb`: those refined coordinates, with the new atom at the peak. The file `{omit model name}_reconstructed_final.pdb` is the REFMAC coordinate set of the last iteration, the one that inserts no atom.
