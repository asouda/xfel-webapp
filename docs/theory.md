# Theoretical background

The three workflows turn deposited observations and atomic coordinates into a map on a grid. This page is the crystallography underneath that calculation: what a structure factor is, where its amplitude and phase come from, how the Fourier transform rebuilds the density, what a refinement cycle changes, and which numbers are already stored in a validation coefficient file.

The formulas for the three maps, and the files and programs each workflow uses, are on the [Fo−Fo](fo-fo.md), [omit-map](omit-map.md), and [2mFo−DFc](two-fo-fc.md) pages. Contours, clipping, and the browser viewer are on the [background](concepts.md) page.

## How structure factors are measured

A crystal repeats the same electron density \(\rho(\mathbf{r})\) on a lattice. An X-ray wave scatters from those electrons. In the far field, the amplitude scattered in one direction is the Fourier component of \(\rho\) at the matching scattering vector. The pattern on the detector is therefore the density expressed in reciprocal space, not a lens image of the molecule. At X-ray wavelengths the crystal itself is the diffraction grating, and the recorded spots are the samples of that transform.

Because the density is periodic, waves scattered from successive unit cells cancel except in the directions that satisfy Bragg’s law,

\[
2\,d\sin\theta = \lambda,
\]

where \(d\) is the spacing of the planes, \(\theta\) is the glancing angle, and \(\lambda\) is the wavelength. Those directions are the points of the reciprocal lattice. Each point is labelled by integer Miller indices \(hkl\), and the scattering vector that lands on it is

\[
\mathbf{s} = h\mathbf{a}^{*} + k\mathbf{b}^{*} + l\mathbf{c}^{*}, \qquad |\mathbf{s}| = \frac{2\sin\theta}{\lambda} = \frac{1}{d}.
\]

The structure factor \(F(hkl)\) is the complex Fourier coefficient at that point:

\[
F(hkl) = \int_{\mathrm{cell}} \rho(\mathbf{r})\, e^{2\pi i\, \mathbf{h}\cdot\mathbf{r}}\, dV = |F(hkl)|\, e^{i\varphi(hkl)}.
\]

\(|F|\) sets the strength of the spot and \(\varphi\) sets the relative timing of the scattered wave. A pixel on the detector can record the strength. It cannot record \(\varphi\).

An XFEL measurement is a still from one crystal orientation, so each frame contains partial reflections. Indexing, integration, scaling, and merging across many crystals produce one merged intensity, and an uncertainty, for each unique \(hkl\). The deposited structure-factor mmCIF is that list, together with the unit cell and the space group. For `6w1o` and `6w1v` the intensity columns are `_refln.intensity_meas` and `_refln.intensity_sigma`.

## How amplitudes are obtained

The measured intensity is proportional to the square of the amplitude, after Lorentz, polarisation, and scale corrections:

\[
I(hkl) \propto |F(hkl)|^{2}.
\]

Squaring removes the phase, so the experiment yields \(|F|\) and not \(F\). A plain square root is a poor estimator when \(I\) is weak or negative, which is common for a reflection whose true amplitude is near the noise. `ctruncate` applies the French–Wilson procedure: given \(I\) and \(\sigma(I)\), and Wilson’s distribution of intensities, it returns the posterior expectation of \(|F|\) and a standard deviation \(\sigma(|F|)\). Those are the amplitudes the later programs call `F` and `SIGF`, or `FP` and `SIGFP` after REFMAC has renamed them.

`cif2mtz` is the step that makes this possible inside CCP4. It copies the mmCIF reflection loop into an MTZ, mapping `intensity_meas` / `intensity_sigma` onto `I` / `SIGI`. The omit-map workflow then runs `ctruncate` on those two columns. The Fo−Fo workflow reaches the same amplitudes by a different route: a zero-cycle `refmac5` job reads the MTZ and writes `FP` and `SIGFP` without moving atoms.

The `9z35` structure-factor file is stored as anomalous intensities, `_refln.pdbx_I_plus` and `_refln.pdbx_I_minus`, with a free-R flag. The sharpened-map workflow does not convert those intensities. It uses the file so `mtzdump` can read the resolution range. The amplitudes that enter the Fourier transform are already present in the validation coefficient file, described below.

## How phases are calculated from the model

The missing phase is computed from the atomic coordinates. Atom \(j\), at position \(\mathbf{r}_j\), with occupancy \(\mathrm{occ}_j\) and isotropic displacement parameter \(B_j\), contributes

\[
\mathrm{occ}_j\, f_j(|\mathbf{s}|)\, \exp\!\big(-B_j\,|\mathbf{s}|^{2}/4\big)\, \exp(2\pi i\, \mathbf{h}\cdot\mathbf{r}_j).
\]

\(f_j(|\mathbf{s}|)\) is the atomic scattering factor, the Fourier transform of that atom’s electron cloud. It decreases as the scattering angle increases, so high-resolution terms are weaker for every atom. The exponential in \(B_j\) is the Debye–Waller factor: a larger displacement parameter damps the high-resolution contribution of that atom still further.

The calculated structure factor is the sum over the model, plus a bulk-solvent term for the disordered solvent in the crystal:

\[
F_c(hkl) = \sum_j \mathrm{occ}_j\, f_j(|\mathbf{s}|)\, e^{-B_j |\mathbf{s}|^{2}/4}\, e^{2\pi i\, \mathbf{h}\cdot\mathbf{r}_j} + F_{\mathrm{solvent}}.
\]

Its modulus is \(|F_c|\) and its phase is \(\varphi_c = \arg(F_c)\). `refmac5` evaluates this sum from the coordinate file and the cell. With `NCYC 0` the coordinates stay where they were deposited, and the MTZ still receives \(F_c\) and \(\varphi_c\).

Which model supplies \(\varphi_c\) depends on the workflow.

- Fo−Fo borrows one phase for the whole difference. The default is `PHWT` from the light validation file, the phase stored with the weighted \(2mF_o-DF_c\) coefficient of the activated structure. The form can switch that to the dark validation phase, or to the phase written by either zero-cycle REFMAC job.
- The omit map calculates \(\varphi_c\) and \(|F_c|\) from the omit coordinates. An atom deleted from that file is absent from \(F_c\), so it can show up in the difference.
- The sharpened \(2mF_o-DF_c\) map uses the phase `PHWT` already stored in the validation file. That phase is the model phase from the refinement that produced the deposition.

## How the density is obtained

The electron density on a real-space grid is the inverse Fourier sum. With \(V\) the unit-cell volume,

\[
\rho(\mathbf{r}) = \frac{1}{V} \sum_{hkl} F(hkl)\, e^{-2\pi i\, \mathbf{h}\cdot\mathbf{r}}.
\]

`fft` evaluates that sum over the resolution shell selected on the form and writes a CCP4 map of the asymmetric unit. The complex number \(F(hkl)\) is the coefficient the workflow has assembled:

| Map | Coefficient inside the sum |
| --- | --- |
| Fo−Fo | \(\big(\|F_{\mathrm{light}}\| - \|F_{\mathrm{dark}}\|\big)\, e^{i\varphi_{\mathrm{ref}}}\), after `scaleit` has matched the overall scale and \(B\) factor of the two amplitude sets |
| mFo−DFc | `DELFWT` with phase `PHDELWT`, where the amplitude is \(m\|F_o\| - D\|F_c\|\) |
| 2mFo−DFc | `FWT` with phase `PHWT`, where the amplitude is \(2m\|F_o\| - D\|F_c\|\), multiplied in `fft` by the sharpening factor \(\exp\!\big(-\Delta B\,(\sin\theta/\lambda)^2\big)\) |

The weight \(m\) is the figure of merit, the expected cosine of the phase error. The factor \(D\) comes from the \(\sigma_A\) treatment of coordinate error: as the model becomes a poorer description of the data, \(D\) falls and \(|F_c|\) is down-weighted. REFMAC stores the already weighted difference as `DELFWT` / `PHDELWT`. The validation file stores the already weighted \(2mF_o-DF_c\) coefficient as `FWT` / `PHWT`.

`fft` writes one asymmetric unit. `extends` expands that map by the space-group symmetry and adds a border, 5 Å unless the form says otherwise, so the grid covers the coordinate file that will be drawn. Coot then cuts a sphere around the chosen residue. PyMOL, and the browser viewer, contour that fragment. A level of \(n\) in sigma mode is \(n\) times the root-mean-square of the map that was written.

## How refinement cycles are performed

REFMAC adjusts the model so that \(|F_c|\) agrees with \(|F_o|\), while bond lengths, angles, planes, chirality, and non-bonded contacts stay near their dictionary values. One cycle is:

1. Compute \(F_c\) and \(\varphi_c\) from the current coordinates, occupancies, and \(B\) factors, including bulk solvent.
2. Scale the observed amplitudes to the calculated ones and build the maximum-likelihood target in \(|F_o|\), plus the geometry restraints.
3. Shift coordinates, \(B\) factors, and occupancies along the gradient of that target.

`NCYC` is the number of times this loop runs. `NCYC 0`, the omit-map default, stops after the first calculation. Atoms stay put, and the MTZ still contains \(F_o\), \(F_c\), `DELFWT`, and `PHDELWT`. A positive cycle count repeats the shift. The figure is then drawn on the shifted coordinates, written as `{omit}_refined.pdb`.

The residual reported for a cycle is the crystallographic \(R\) factor,

\[
R = \frac{\sum \big|\,|F_o| - |F_c|\,\big|}{\sum |F_o|}.
\]

Reflections tagged with a free-R flag are left out of the target function. Their residual is \(R_{\mathrm{free}}\). In the `6w1o` and `6w1v` intensity files that flag is the reflection status: `_refln.status` is `o` for an observation used in refinement and `f` for a free-set reflection. `cif2mtz` stores the free set as the MTZ column `FREE`. The `9z35` structure-factor file stores it explicitly as `_refln.pdbx_r_free_flag`. Before REFMAC, the omit workflow runs `cad` to copy `FREE` from the intensity MTZ onto the amplitude MTZ, so the free set stays free.

`MAKE HYDR Y` and `MAKE HOUT Y` tell REFMAC to build riding hydrogens and to write them into the coordinate file. External distance and angle restraints, or a monomer library written by `libcheck`, are placed in front of those keywords when the form lists them. They matter once `NCYC` is greater than zero: refinement will try to explain leftover density by moving nearby atoms, and the restraints limit how far selected cluster distances and angles can go.

The Fo−Fo workflow also calls `refmac5`, but only with `NCYC 0`, and only to place each amplitude set in a REFMAC-style MTZ. Those coordinates are not refined against each other.

## What the validation coefficient files store

A validation coefficient file is the output of a finished refinement, written so that the standard maps can be redrawn without repeating it. It is a reflection mmCIF with the cell and space group, and one row per \(hkl\). It does not contain the measured intensities.

For the entries used here, `{id}_validation_2fo-fc_map_coef.cif` has three coefficient columns:

| mmCIF column | MTZ label after `cif2mtz` | What the number is |
| --- | --- | --- |
| `_refln.pdbx_FWT` | `FWT` | \(2m\|F_o\| - D\|F_c\|\), the weighted amplitude |
| `_refln.pdbx_PHWT` | `PHWT` | the phase that goes with `FWT`, the model phase \(\varphi_c\) |
| `_refln.fom` | `FOM` | the figure of merit \(m\) already included in `FWT` |

The sharpened-map workflow transforms `FWT` and `PHWT`. It leaves `FOM` unread, so the figure of merit is applied once. The Fo−Fo workflow reads the same file only when the phase source is a validation file, and then only to copy `PHWT` into the combined MTZ.

Some depositions also include `{id}_validation_fo-fc_map_coef.cif`. For `9z35` its columns are:

| mmCIF column | What the number is |
| --- | --- |
| `_refln.pdbx_DELFWT` | \(m\|F_o\| - D\|F_c\|\) |
| `_refln.pdbx_DELPHWT` | the phase of that weighted difference |
| `_refln.fom` | the figure of merit |

The omit workflow does not open this companion file. It asks REFMAC to recompute `DELFWT` and `PHDELWT` from the omit coordinates and the observed amplitudes of the dataset on the form. That recalculation is what makes the result an omit map rather than a redraw of the deposited difference.
