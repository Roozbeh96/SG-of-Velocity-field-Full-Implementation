# Stochastic Generation of 2D-TBL Velocity Field

Tools to stochastically generate turbulent boundary-layer (TBL) velocity fields—from synthetic profile generation through reorganization and vortex addition. Available for MATLAB and Python. The full framework implemented here, and its validation against the high-Reynolds-number rough-wall boundary layer of Saddoughi & Veeravalli, is described in:

**R. Ehsani and M. Guala, A stochastic modeling framework to generate 2-D rough-wall high-Reynolds-number turbulent boundary layers, arXiv:2609.10236 (2026).** https://arxiv.org/abs/2609.10236

"The paper" below (Figures 5, 6 and 7) refers to this preprint. The individual parts of the model are described in the following papers:

1. https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/stochastic-modelling-of-the-instantaneous-velocity-profile-in-roughwall-turbulent-boundary-layers/492F3CD03C8C3E7ED306E9117B848B5E

2. https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/stochastic-modal-velocity-field-in-roughwall-turbulence/800D97E60609F4FE4D92847CCDB1C7A0

3. https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/stochastic-generation-of-velocity-fields-to-reproduce-the-energy-spectrum-in-wall-turbulence/1D7BDC29D1D577DAE7672E86E21641B7

4. https://journals.ametsoc.org/view/journals/bams/105/1/BAMS-D-23-0066.1.xml

The generated field is compared with the measurements of Saddoughi & Veeravalli (U<sub>∞</sub> = 10 and 50 m/s), digitized in `Matlab/src/Data/Saddoughi` and used by both packages:

S. G. Saddoughi and S. V. Veeravalli, Local isotropy in turbulent boundary layers at high Reynolds number, J. Fluid Mech. 268, 333–372 (1994). https://doi.org/10.1017/S0022112094001370

## Contents
- [MATLAB package](#matlab-package)
- [Python package](#python-package)
- [Sample of stochatically generated velocity field](#Stochatically-generated-velocity-field)
- [Citation](#citation)
- [License](#license)

## MATLAB package
The MATLAB package (R2025b) is in the `Matlab/src` folder.

<details>
<summary>Click to see the documentation for MATLAB package</summary>

### MATLAB files
Image processing toolbox should be installed.<br>
All the code files are in `Matlab/src` folder:

`Matlab/src`
- main.m: defining initial conditions. The documentation for the code is provided in the files as comments.
- stochastic_generation.m: defining an object and transfers to different functions
- read_data.m: reads the Saddoughi data (`Matlab/src/Data/Saddoughi`) used in the plots. Set the paths of the JSON files to the location of the repository on your computer.

`Matlab/src/utils/SG_VelProf.m`
- SGVP: code to generate velocity profiles. `genGaussCop.m` generates correlated random numbers from a Gaussian copula.

`Matlab/src/utils/SG_VelField.m`, `Matlab/src/utils/adding_shear_layer.m`
- adding_shear_layer: smoothing the velocity jump at neighboring UMZs.
- SGVF: code to reorganize the generated profiles into correlated field and generate low-resolution velocity field.

`Matlab/src/utils/increasing_res.m`
- Increasing the resolution of the generated velocity field. 

`Matlab/src/utils/Lambda_ci.m`
- Lambda_ci: In this function, swirling strength of the generated velocity field are calculated. 

`Matlab/src/utils/cdf_omega.m`
- The distribution of the vorticity of the generated velocity field for the percentile matching algorithm is computed.

`Matlab/src/utils/VorX_generator.m`
- In this function, the vortex features (radius, azimuthal velocity, correlation coefficient) are generated. The vorticity of the generated vortices are employed in the percentile matching algorithm.

`Matlab/src/utils/SG_VorX.m`
- SGVorX: code to add near wall vortex cores.
- Lambdaci_VorX: code to insert vortices into shear layers

`Matlab/src/utils/MVG_kernel.m`
- Convolving the small-scale 2-D Gaussian filter with the generated velocity field.

`Matlab/src/utils/spec_analysis_E11.m`
- spec_analysis_E11: In this function, Power spectrum density ($E_{11}$) of the generated velocity field is computed.
- The same function for the vertical velocity ($E_{22}$), and the co-spectra ($E_{12}$) are provided.

`Matlab/src/utils/str_func_D11.m`
- str_func_D11: In this function, structure function ($D_{11}$) of the generated velocity field is computed.
- The same function for the vertical velocity ($D_{22}$) is provided.

`Matlab/src/utils/production_comp.m`
- Production of turbulent kinetic energy.

`Matlab/src/utils/plots.m`
- After generating the object (and running `read_data.m`), plot Figures 5, 6 and 7 of the paper: mean velocity and second moments, spectra and structure function, and co-spectrum, compared with the Saddoughi data (U<sub>∞</sub> = 50 m/s). Run each section to draw its figure. The figures use LaTeX text (Computer Modern font).

</details>

## Python package
The Python package (`sgvf`, Python 3.11–3.13) is a port of the MATLAB package. It uses the same steps, parameters and seeds, and is in the `Python` folder.

<details>
<summary>Installation</summary>

1. 🛠️ Installing Poetry

To install [Poetry](https://python-poetry.org/) (Python dependency management and packaging tool), run the following command in your terminal:

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

After installation, make sure Poetry is in your `PATH`
- macOS/Linux
```bash
export PATH="$HOME/.local/bin:$PATH"
```

Verify installation:

```bash
poetry --version
```

Optional: keep the virtual environment inside the project (works great with VS Code):

```bash
poetry config virtualenvs.in-project true
```

2. Installing the environment

Go to the `Python` folder and install the dependencies (numpy, scipy, matplotlib, tqdm, ipykernel) together with the `sgvf` package:

```bash
cd Python
poetry install
```

Any Python version from 3.11 to 3.13 works. To pick a specific interpreter, run `poetry env use python3.12` (for example) before `poetry install`.

</details>

<details>
<summary>Usage</summary>

Run the whole pipeline (the same steps as `main.m`) and save the object to `Gen_sample.pkl`:

```bash
cd Python
poetry run python main.py
```

The physical parameters (u<sub>τ</sub>, z<sub>0</sub>, k<sub>s</sub>, δ, λ<sub>T</sub>, ν) are set in the `if __name__ == "__main__":` block at the end of `main.py`; the Saddoughi (U<sub>∞</sub> = 50 m/s) values are used by default.

Options: `--n-prof` (number of profiles, default `1.6e5`), `--output` (output file, `''` to skip saving), `--seed` (seed of the field reorganization and vortex insertion; by default every run is different, like `rng('shuffle')` in MATLAB), `--parallel` (see below).

**Parallel or series computing.** The steps that use `parfor` in the MATLAB version (`increasing_res`, `lambda_ci` and `cdf_omega`) can run in series (the default) or in parallel. Both give identical results; parallel is faster on a multi-core computer.

- Command line: add `--parallel` to run these three steps in parallel, e.g. `poetry run python main.py --parallel`. Without it, they run in series.
- Python: pass `parallel=True` (and optionally `n_workers`, the number of workers; `None` uses all CPU cores) to each function:
  ```python
  increasing_res(Gen_sample, L, Kr, parallel=True)
  lambda_ci(Gen_sample, parallel=True)
  cdf_omega(Gen_sample, parallel=True)
  ```
  `increasing_res` and `lambda_ci` use threads; `cdf_omega` uses processes, so a script that calls it with `parallel=True` needs an `if __name__ == "__main__":` guard (as in `main.py`).
- VS Code: `Python/.vscode/launch.json` has the run configurations "main (series)" and "main (parallel)" (and "plots", see below). Choose one in the Run and Debug panel.

**Experimenting with the steps on a saved field (`main_jupyter.py`).** `main_jupyter.py` runs the same pipeline as `main.py`, but every step is a function with its own `# %%` cell, so it can be used in two ways:

- No `.pkl` file yet: run it as a Python file, e.g. `poetry run python main_jupyter.py` (same options as `main.py`) or "Run Python File" in VS Code. It goes through the whole pipeline from zero to the end and saves `Gen_sample.pkl`.
- With a saved `.pkl` file: run its sections interactively in VS Code ("Run Cell" above a `# %%` line). Run the first two cells once: the first defines the functions, the second loads `Gen_sample` from `PKL_FILE` (skipped if it is already in memory, e.g. loaded by `run_plots.py` in the same interactive window). Then run any step cell on its own, without regenerating or reorganizing the field:

  | Cell | What it does |
  |---|---|
  | Increase the resolution | rebuilds `Gen_u_HRVF`/`Gen_w_HRVF` from the saved LRVF (one step back: the HRVF without vortices) |
  | Add vortices | computes `Lambda_ci` and the c.d.f. of omega of the current HRVF and inserts the vortices |
  | Gaussian filter | applies the small-scale 2-D Gaussian filter |
  | Analysis | spectra (E11, E22, E12), structure functions (D11, D22) and production of the current HRVF |
  | Save | saves the object to `Gen_sample_new.pkl` (the original file is kept) |

  For example, "Increase the resolution" followed by "Analysis" gives the spectra and structure functions of the field without vortices, and "Add vortices", "Gaussian filter" and "Analysis" give them for the final field again. `PARALLEL` and `SEED` in the loading cell set the options of the steps.

  The steps overwrite the object in memory; results of later steps (e.g. `E11` after going back) keep their old values until their step is rerun. Adding the vortices again reproduces the saved field only if the same seed was used; to get the saved field back, run `del Gen_sample` and the loading cell again.

**Plotting the results.** `run_plots.py` loads the saved object and draws Figures 5, 6 and 7 of the paper (the figures of `plots.m`), comparing the generated field with the Saddoughi data, and the additional figures listed below. It can be used in two ways:

- Interactively, like MATLAB sections (recommended): in VS Code, click "Run Cell" above a `# %%` line (or press Shift+Enter). Run the first cell once to read `Gen_sample.pkl`; the object then stays in memory, and each figure has its own cell that you can rerun as often as you like. Changes saved in `sgvf/plots.py` are reloaded automatically, so you can edit a figure and rerun its cell without reading the file again. The last cell saves all figures to `Python/figures`. This needs VS Code's Jupyter extension; `ipykernel` is installed with the environment.
- As a script: `python run_plots.py`, or "plots" in the Run and Debug panel of VS Code, draws all figures, saves them to `Python/figures` and opens them.

The file to plot and the output folder are set by `PKL_FILE` and `SAVE_DIR` at the top of `run_plots.py`. Each figure cell saves its figure as a PNG at 600 dpi (`fig.savefig(SAVE_DIR / "name.png", dpi=600)`); `plot_all(..., dpi=600)` does the same for all figures.

Besides Figures 5, 6 and 7, `run_plots.py` has cells for:

| Cell | What it shows |
|---|---|
| Compensated co-spectrum | k<sub>1</sub><sup>7/3</sup>\|E<sub>12</sub>\| of the smoothed generated co-spectrum and of Saddoughi's data; a −7/3 range is a plateau. The shaded band is the range in `fit_range` and the dotted line the constant `C` (default 0.15). |
| Co-spectrum with a power-law fit | Figure 7 with a straight line fitted in log–log to the smoothed generated co-spectrum (`fit_range`) and to Saddoughi's points (`sadd_fit_range`). The slopes ± their standard errors are shown in the legend and printed with R<sup>2</sup>, to compare with −7/3. |
| Distribution of u′w′ | `plot_uw_distribution`, see Python files. |
| In-plane divergence | `plot_divergence`, see Python files. |

The fit ranges are in k<sub>1</sub>ε<sup>1/2</sup>S<sup>−3/2</sup>, e.g. `plots.plot_E12(Gen_sample, sadd, fit_range=(5.5, 22), sadd_fit_range=(1, float("inf")))`. `sgvf.plots.fit_power_law(x, y, fit_range)` returns the slope, its standard error and R<sup>2</sup> of such a fit.

The same can be done from your own code:

```python
from sgvf import StochasticGeneration
from sgvf.plots import plot_all

Gen_sample = StochasticGeneration.load("Gen_sample.pkl")
plot_all(Gen_sample, save_dir="figures", show=True)
```

Each step can also be run on its own, exactly like in `main.m`:

```python
from sgvf import StochasticGeneration, increasing_res, lambda_ci, vorx_generator, cdf_omega

Gen_sample = StochasticGeneration(u_tau, z_0, ks, delta, lambda_T, nu, z_min, z_max, res_z)
Gen_sample.SGVP(ro_uw, N_prof)
Gen_sample.SGVF(Delta_x, M)
increasing_res(Gen_sample, L, Kr)
lambda_ci(Gen_sample)
vorx = vorx_generator(Gen_sample)
cdf_omega(Gen_sample)
Gen_sample.SGVorX(**vorx)
```

</details>

<details>
  <summary>Click to see the documentation for Python files</summary>

### Python files
All the code files are in the `Python` folder. Each MATLAB file has a Python counterpart with the same name in snake case:

`Python/main.py`
- Defines the initial conditions and runs the whole pipeline (`main.m`). The documentation for the code is provided in the files as comments. `main(u_tau, z_0, ks, delta, lambda_T, nu, ...)` can also be called from other scripts.

`Python/main_jupyter.py`
- Same pipeline as `main.py`, split into step functions (`setup`, `sgvp`, `sgvf`, `increase_resolution`, `add_vortices`, `gaussian_filter`, `analysis`, `save`). Run as a script for the whole pipeline, or run its `# %%` cells on a saved object to repeat single steps (see Usage).

`Python/run_plots.py`
- Plots Figures 5, 6 and 7, the compensated co-spectrum, the co-spectrum with a power-law fit, the distribution of u′w′ and the divergence profiles from a saved object, interactively (`# %%` cells) or as a script (see Usage).

`Python/.vscode/launch.json`
- VS Code run configurations: "main (series)", "main (parallel)" and "plots".

`Python/sgvf/stochastic_generation.py`
- `StochasticGeneration`: the object that holds the parameters and results and passes them to the different functions (`stochastic_generation.m`). `save`/`load` write and read the object with pickle.

`Python/sgvf/sg_vel_prof.py`
- `sg_vel_prof`: generates the velocity profiles (`SG_VelProf.m`). `gen_gauss_cop.py` generates correlated random numbers from a Gaussian copula (`genGaussCop.m`).

`Python/sgvf/sg_vel_field.py`, `Python/sgvf/adding_shear_layer.py`
- `adding_shear_layer`: smooths the velocity jump at neighboring UMZs.
- `sg_vel_field`: reorganizes the generated profiles into a correlated low-resolution velocity field.

`Python/sgvf/increasing_res.py`
- Increases the resolution of the generated velocity field (makima interpolation and band-stop filter). It can run in series or in parallel (`parallel=True`).

`Python/sgvf/lambda_ci.py`
- `lambda_ci`: computes the swirling strength of the generated velocity field. It can run in series or in parallel (`parallel=True`).

`Python/sgvf/cdf_omega.py`
- Computes the distribution of the vorticity of the generated velocity field for the percentile matching algorithm. It can run in series or in parallel (`parallel=True`).

`Python/sgvf/vorx_generator.py`
- Generates the vortex features (radius, azimuthal velocity, correlation coefficient) and splits the vortices into deciles of their vorticity for the percentile matching algorithm.

`Python/sgvf/sg_vorx.py`
- `sg_vorx`: adds near-wall vortex cores (`SG_VorX.m`).
- `lambdaci_vorx`: inserts vortices into shear layers (`Lambdaci_VorX.m`).

`Python/sgvf/mvg_kernel.py`
- `small_scale_sigma`: standard deviation of the small-scale 2-D Gaussian filter (the `fzero` section of `main.m`).
- `mvg_kernel`: convolves the 2-D Gaussian filter with the field (`MVG_kernel.m`).

`Python/sgvf/spectral.py`
- `spec_analysis_E11`, `spec_analysis_E22`, `spec_analysis_E12`: power spectral densities ($E_{11}$, $E_{22}$) and co-spectrum ($E_{12}$) of the generated velocity field.

`Python/sgvf/structure_function.py`
- `str_func_D11`, `str_func_D22`: structure functions ($D_{11}$, $D_{22}$) of the generated velocity field.

`Python/sgvf/production.py`
- `production_comp`: production of turbulent kinetic energy.

`Python/sgvf/read_data.py`
- `read_saddoughi`: reads the Saddoughi data (`read_data.m`). The JSON files are shared with the MATLAB package (`Matlab/src/Data`).

`Python/sgvf/plots.py`
- The figures of `plots.m`, used by `run_plots.py`: `plot_moments` (Figure 5: mean velocity and second moments), `plot_spectra` (Figure 6: E11 spectra and D11 with the premultiplied E11 inset) and `plot_E12` (Figure 7: co-spectrum; with `fit_range`/`sadd_fit_range` it adds the power-law fits). `plot_E12_compensated` shows the co-spectrum compensated by k<sub>1</sub><sup>7/3</sup>, and `fit_power_law` fits a power law in log–log and returns the slope with its standard error. `plot_uw_distribution` shows the distribution of the instantaneous u′w′/u<sub>τ</sub><sup>2</sup> with elevation, from the wall upward: (a) the PDFs at selected elevations stacked by z/δ, (b) the log<sub>10</sub> PDF at every elevation with the mean profile. `plot_divergence` shows profiles of the normalized in-plane divergence ∂u′/∂x + ∂w′/∂z versus z/δ: (a) the correlation coefficient of ∂u′/∂x and ∂w′/∂z, (b) R = ⟨div<sup>2</sup>⟩/(⟨(∂u′/∂x)<sup>2</sup>⟩ + ⟨(∂w′/∂z)<sup>2</sup>⟩), (c) div<sub>rms</sub>/ω′<sub>rms</sub>, with reference lines for a strictly 2D divergence-free field (−1, 0, 0), an x–z plane of 3D locally isotropic turbulence (−1/2, 1/2, 1/√5) and gradients with no continuity constraint (0, 1). `plot_all` draws all six. The text uses the Computer Modern font (`cmr10` with mathtext), like MATLAB's LaTeX interpreter.

### Differences from the MATLAB package
- Indices are 0-based: MATLAB `Gen_u_HRVF(50,:,1)` is `Gen_u_HRVF[49, :, 0]` in Python. Array shapes are the same as in MATLAB.
- The random number generators are different, so the generated samples are not identical to MATLAB's even with the same seeds. The statistics are the same.
- The logs (`log_data_SGVP`, `near_wall_VorX_data`, `prograde_VorX_data`, `retrograde_VorX_data`, `hist_corr`) are numpy arrays (one row per entry) instead of cells of dictionaries. The column names are in `sgvf.sg_vel_prof.LOG_COLUMNS`, `sgvf.sg_vorx.NEAR_WALL_COLUMNS` and `sgvf.sg_vorx.LAMBDACI_COLUMNS`.
- Memory: with the default parameters each high-resolution field (`Gen_u_HRVF`, `Gen_w_HRVF`, `Lambda_ci`, `omega`) is about 3.2 GB, the same as in MATLAB. Use a machine with enough memory, or reduce `N_prof`. You can leave large fields out of the saved file with `Gen_sample.save(path, exclude=("Lambda_ci", "omega"))`.
</details>

## Stochatically Generated Velocity Field
<img src = "Gen.gif">

## Citation
If you use this code, please cite:

```bibtex
@misc{ehsani2026stochasticmodelingframeworkgenerate,
      title={A stochastic modeling framework to generate 2-D rough-wall high-Reynolds-number turbulent boundary layers},
      author={Roozbeh Ehsani and Michele Guala},
      year={2026},
      eprint={2609.10236},
      archivePrefix={arXiv},
      primaryClass={physics.flu-dyn},
      url={https://arxiv.org/abs/2609.10236},
}
```

## License
This work is licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/) (Attribution-NonCommercial-ShareAlike 4.0 International). See [LICENSE](LICENSE).
