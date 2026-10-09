# Activation Background Analysis Pipeline — Summary and User Guide

## Overview

This pipeline estimates the activation-induced detector background in a space instrument after irradiation by trapped protons. It takes as input the raw output of Geant4 Monte Carlo simulations (activation score files for a set of monochromatic proton energies) and a SPENVIS orbital proton flux file, and produces time-resolved background spectra and count rates.

---

## Pipeline architecture

```
0. Geant4 Activation run	                                                 SPENVIS flux file
    0_run.py                                                                         │
    Geant4 score CSVs                                                                │
    (one per proton energy)                                                          │
        │                                                                            │
        ▼                                                                            │
1. Activation parser                                                                 │
    activation_parser.py                                                             │
    → results.pkl                                                                    │
        │                                                                            │
        ▼                                                                            │
2. build_decay_chains.py                                                             │
    (uses decay_chain_builder)                                                       │
    → DecayChains/*.dat                                                              │
        │                                                                            │
        ▼                                                                            │
3. compute_activities.py                                                             │
    (uses bateman)                                                                   │
    → activities.pkl                                                                 │
    → active_isotopes.pkl —————————————————————————-┐                                │
        |                                           |                                │
        |                                           ▼                                │
        |                              4. Geant4 PostActivation run                  │
        |                                  (scripts not yet in the repository)       │
        |                                           |                                │
        ▼                                           ▼                                │
 6. accumulate_spectra.py  ←── result_spectra/{vol}_{iso}_S-mode.dat                 │
                           ←── 5. SPENVIS parser  ←──--------------------------------┘
						          spenvis_parser.py
    → spectra.pkl
    → count_rate.dat
        │
        ├──▶ 7. activation_history.py
        │        → time history of count rate (single orbit or long-term)
        │
        └──▶ 8. average_spectrum.py
                 → steady-state out-of-belt spectrum with line IDs
```

---

## Script reference

### 0. `0_run.py`

Runs the Geant4 activation simulation once per monochromatic energy. Each energy runs in its own work directory, `output/work/<E>_MeV/`, which links in every file of the mass-model directory. The per-thread ntuple files are merged into `output/scorefile_<E>_MeV_<nprim>.csv`, a plain CSV with the columns `EventID,Isotope,Lifetime,Volume`. The source geometry, the beam area and the primaries of each run are recorded in `output/run_info.json`, which the later steps read.

**Inputs:** Geant4 executable, mass-model directory\
**Outputs:** `output/scorefile_*_MeV_*.csv`, `output/run_info.json`, `output/work/`

```bash
python 0_run.py                                   # CUSP defaults
python 0_run.py --config config_cusp.toml         # same, from the config file
python 0_run.py --energies 10 30 100 --nprim 1e6 -t 8 --outdir test/
python 0_run.py --g4-args=-n --outdir output_nact # extra executable options
python 0_run.py --dry-run                         # write macros, print commands
```

Every setting has a CUSP default. It can be set in a TOML file (`--config`; see `config_cusp.toml`, which documents all the keys) and overridden on the command line. Energies already done are skipped on rerun unless `--overwrite` is given. The script refuses to add runs with a different source geometry to an existing output directory.

**Source geometry.** Primaries start on a sphere of radius R and point inwards with a cosine law limited to θ ≤ θmax. For an omnidirectional integral flux F, the number of primaries for an exposure Δt is N = Δt · F · πR² sin²θmax. This holds only if:
- the whole mass model lies inside the sphere of radius R sin θmax around the centre;
- the source sphere lies inside the Geant4 world. Primaries that start outside the world are not tracked through the geometry.

The CUSP defaults are R = 13 cm, centre (0, −1.5, 4.8) cm and θmax = 90°. The world is a 50 cm box, and the payload's bounding box is x −49…47, y −70…40, z −46…142 mm.

| Key / argument | CUSP default | Description |
| --- | --- | --- |
| `executable` | `../Debug/cusp-activation` | Geant4 application |
| `geometry_dir` | `../gdml-mass-model` | Mass-model files |
| `threads` / `-t` | `0` (all cores) | Worker threads |
| `g4_args` / `--g4-args=` | none | Extra executable arguments, e.g. `-n` |
| `pre_commands` / `--pre-command` | `/cusp/stepping/verbose 0` | Macro commands before `/run/beamOn` |
| `energies` | 5 … 700 MeV | Primary energies |
| `nprim` | `1e8` | Primaries per energy |
| `radius_cm`, `centre_cm`, `thetamax_deg` | 13, (0, −1.5, 4.8), 90 | Source sphere |
| `outdir` | `output` | Output directory |

---

### 1. `activation_parser.py`

Parses the score CSV files into a pandas DataFrame indexed by `(energy_MeV, volume, isotope)`. It has these columns:
- `count`
- `efficiency` \[nuclides/primary\]
- `efficiency_err` (Poisson)
- `half_life_s` \[s\]

Isotope names are converted to the canonical form defined in `nuclides.py`, e.g. `Al26[228.305]` or `Ta180[77.200X]`. Every step uses this form.

The number of primaries of each file is read from `run_info.json`. `df.attrs` carries `nprim` (`{energy: primaries}`, including energies with no recorded nuclide) and `source` (the source geometry). Both are saved in the pickle.

**Inputs:** score CSV files (+ `run_info.json`)\
**Output:** `results.pkl`

```bash
# Whole output directory of 0_run.py
python activation_parser.py --dir output/ --save-pkl results.pkl

# Single file (energy inferred from filename)
python activation_parser.py --file output/scorefile_10_MeV_1e+08.csv --save-pkl results.pkl

# Score files without run_info.json
python activation_parser.py --dir old_results/ --nprim 1e8 --save-pkl results.pkl
```

**Key arguments:**

| Argument | Default | Description |
| --- | --- | --- |
| `--file` / `--dir` | — | Single CSV or directory of CSVs |
| `--nprim` | from `run_info.json` | Primaries per file. Only needed without `run_info.json`, and an error if it contradicts it |
| `--save-pkl` | — | Output pickle path |
| `--save-fig` | — | Summary figure path (otherwise shown on screen) |
| `--top` | `3` | Entries in statistics report |

**Library usage:**

```python
from activation_parser import parse_score_files, load, all_isotopes, active_volumes

df = parse_score_files([("scorefile_10_MeV.csv", 10.0, 1e8),
                         ("scorefile_20_MeV.csv", 20.0, 1e8)])
# df.loc[(10.0, "mask_phys")]          → all isotopes in mask_phys at 10 MeV
# df.xs("Re182", level="isotope")      → Re182 across all energies and volumes
```

---

### Helper library: `decay_chain_builder.py`

Builds the radioactive decay chains of a nuclide, linearised into paths, each with its branching ratio. It writes one `.dat` file per nuclide.

**Data sources.** Decay modes and branching ratios come from the Geant4 RadioactiveDecay library (`z<Z>.a<A>` files). Half-lives come from `ENSDFSTATE.dat` (G4ENSDFSTATE), as in Geant4 itself. Both are found through `$GEANT4_DATA_DIR`, which `geant4.sh` sets, unless `--db`/`--ensdf` or `$G4RADIOACTIVEDATA`/`$G4ENSDFSTATEDATA` are given.

**Level matching** follows `G4VRadioactiveDecay::LoadDecayTable`, so the chains describe what the Geant4 decay simulation will do:
- A level with a floating-level letter (`Ta178[0.000X]`) needs both energy and letter to match.
- A level without a letter takes the first library level with the same energy. It keeps its own ENSDF half-life.
- An excited level missing from the library decays by IT to the ground state.
- A floating E = 0 level missing from the library has no decay data. Geant4 kills it without emitting anything, so it is terminal.

A warning is printed for each level whose data was borrowed or assumed.

**Corrections to RadioactiveDecay6.1.2:**
- **Half-lives:** the `P`-line half-life is never used. 80 of those values are wrong, e.g. swapped between ground state and isomer for Sc44, Eu152 and Cu68.
- **Branching:** each mode's total comes from its summary line (a fraction) and is shared among daughter levels in proportion to the detail lines (in percent). Previously the two were summed, so an isomer's IT fraction could drop from 98.8% to 44%.
- **Self-transitions:** spurious IT branches of ten ground states onto themselves (e.g. Y92 50%, Tb154 18%) are dropped.

**Chain structure.** Nuclides with T½ < `--min-halflife` are not recorded: a step goes directly to the first long-lived (or stable) nuclides that such prompt decays lead to. The daughter of step *i* is therefore always the parent of step *i+1*. Branches below `--prune` are dropped, and the rest are renormalised to 1.

Each file starts with a `# params:` line (format version, `prune`, `min_halflife`, data libraries). An existing file is rebuilt automatically when that line differs.

**Inputs:** Geant4 RadioactiveDecay and G4ENSDFSTATE libraries\
**Output:** `DecayChains/<isotope>.dat`

```bash
# Single isotope
python decay_chain_builder.py Os172 --outdir DecayChains/

# Multiple isotopes, explicit data paths
python decay_chain_builder.py Re182 Cs134 Ba133 \
    --db /path/to/RadioactiveDecay6.1.2 --ensdf /path/to/G4ENSDFSTATE3.0
```

**Key arguments:**

| Argument | Default | Description |
| --- | --- | --- |
| `--db` | `$G4RADIOACTIVEDATA`, then `$GEANT4_DATA_DIR/RadioactiveDecay*` | Decay library |
| `--ensdf` | `$G4ENSDFSTATEDATA`, then `$GEANT4_DATA_DIR/G4ENSDFSTATE*` | Half-lives |
| `--outdir` | `DecayChains` | Output directory |
| `--prune` | `0.001` | Discard branches with cumulative BR below this (0.1%, as in the paper) |
| `--min-halflife` | `1e-6 s` | Nuclides shorter-lived than this are passed through |
| `--overwrite` | `False` | Rebuild files even if up to date |

---

### 2. `build_decay_chains.py`

Convenience wrapper. It reads `results.pkl`, extracts all unique isotopes, and builds their chains with `decay_chain_builder`. Files that are up to date are skipped; a file for a nuclide that has no chain any more is deleted. The summary reports how many library half-lives were replaced by ENSDF values.

**Inputs:** `results.pkl`\
**Output:** `DecayChains/*.dat`

```bash
source ~/geant4-v11.4.1/bin/geant4.sh       # sets GEANT4_DATA_DIR
python build_decay_chains.py results.pkl --outdir DecayChains/
```

---

### Helper library: `bateman.py`

Solves the Bateman equations of a linear chain (paper Eq. A7) for N₀(0) = 1. `chain_populations(lambdas, times)` returns Nₖ(t), and `chain_activities` returns λₖNₖ(t).
- **Method:** the analytic solution (Eq. A8–A9) in double precision. Wherever its terms cancel by more than a factor of 10³, the value is recomputed with mpmath, raising the precision until two successive evaluations agree to 10⁻¹³. The cancellation can exceed a hundred orders of magnitude: for close decay constants, for half-lives spanning many decades, or at times short compared with the half-lives.
- **Equal decay constants** (different nuclides with the same tabulated half-life) are separated by a relative 10⁻¹².
- **Accuracy:** better than 10⁻¹² against the exact formula evaluated with 1500 digits (random chains of 2–8 members, half-lives from 0.1 s to 10²⁰ s) and against a stiff ODE solver (Radau) on the longest chains of a CUSP run.

---

### 3. `compute_activities.py`

Computes the activities after an instantaneous irradiation, in Bq/primary, for every `(energy_MeV, volume, isotope)`. The chains depend only on the nuclide produced. So the unit activities a_{r,i}(t) of every nuclide *i* in the chains of each produced (root) nuclide *r* are solved once, then combined with the production yields:

A_{E,v,i}(t) = Σ_r eff_{E,v,r} · a_{r,i}(t)

The activity of a nuclide is assigned to the volume where its root was produced. A produced nuclide without a chain file only contributes its own activity λe^{−λt}, and is listed in a warning. This is the case for nuclides with no decay data, such as Ta178[0.000X]. With `--errors`, the Poisson uncertainty of the yields is propagated.

For a CUSP run (17 energies, 37,624 yield entries, 758 produced nuclides, 201 times) this takes about 11 s.

**Inputs:** `results.pkl`, `DecayChains/*.dat`\
**Outputs:** `activities.pkl`, `activities_err.pkl` (with `--errors`), `unit_activities.pkl`, `active_isotopes.pkl`

```bash
python compute_activities.py results.pkl --chains DecayChains/ --outdir output/ --errors
```

**Key arguments:**

| Argument | Default | Description |
| --- | --- | --- |
| `--chains` | `DecayChains` | Directory with chain `.dat` files |
| `--tmin`, `--tmax`, `--per-decade` | `0.1`, `1e9`, `20` | Log-spaced time grid \[s\] (201 points) |
| `--times` | — | Explicit evaluation times \[s\], overriding the grid |
| `--threshold` | `1e-15` Bq/primary | Min peak activity for `active_isotopes` |
| `--errors` | off | Also save the 1σ statistical uncertainties |
| `--outdir` | `.` | Output directory |
| `--fmt` | `pkl` | `pkl` or `parquet` for the activity tables |

**Output structure:**
- **`activities.pkl`:** MultiIndex `(energy_MeV, volume, isotope)` × time columns `t_<seconds>`, in Bq/primary. `attrs` carries `times`, and `nprim` and `source` from `results.pkl`.
- **`activities_err.pkl`:** same layout, 1σ uncertainty.
- **`unit_activities.pkl`:** MultiIndex `(root, isotope)` × time columns, in Bq per root nucleus.

```python
from compute_activities import load_outputs, times_of
df, active_isos = load_outputs("output/")
t = times_of(df)                                  # evaluation times [s]
df.loc[(100.0, "PV-Absorber_001")]                # all nuclides in one volume at 100 MeV
df.groupby(level="isotope").sum()                 # summed over energies and volumes
```

---

### 4. Geant4 post-activation run

Execute the post-activation Geant4 run, using the list of active isotopes obtained in the previous step `active_isotopes.pkl`. The scripts for this step are not yet in the repository.

**Inputs:** `active_isotopes.pkl`

**Outputs:** `result_spectra/` directory

---

### 5. `spenvis_parser.py`

Parses a SPENVIS AP9/AE9 integral proton flux output file and returns:

- Mean integral flux in each energy band \[protons/cm²/s\]

- Number of simulation primaries per band (given simulation geometry)

- Plots: differential flux spectrum, total flux vs time, primaries per band

**Input:** SPENVIS AP9/AE9 output text file\
**Outputs:** flux values, plots

```bash
python spenvis_parser.py AP8MIN.AP8.output_mean_flux_550km_SSO.txt \
    --R 13 --thetamax 90 --timestep 60 \
    --in-belt-only \
    --save-diff diff_flux.pdf \
    --save-time flux_time.pdf \
    --save-prim primaries.pdf
```

**Key arguments:**

| Argument | Default | Description |
| --- | --- | --- |
| `--energies` | `[7,10,...,400]` MeV | Band edge energies |
| `--R` | `5000` cm | Source sphere radius (simulation geometry) |
| `--thetamax` | `0.8021` deg | Cone half-angle |
| `--timestep` | `60` s | Irradiation timestep duration |
| `--in-belt-only` | `False` | Average over in-belt steps only |

**Primaries formula:**

```
N = timestep × F_band [p/cm²/s] × π × R² × sin²(θ_max)
```

**Library usage:**

```python
from spenvis_parser import parse_spenvis, band_fluxes_to_primaries

_, band_fluxes = parse_spenvis("AP8MIN.AP8.output_mean_flux_550km_SSO.txt", energies=[7,10,...,400])
n_prim = band_fluxes_to_primaries(energies, band_fluxes,
                                   timestep=60, R=13, thetamax=90)
```

---

### 6. `accumulate_spectra.py`

Computes the background spectrum and count rate at each decay time after **one SPENVIS time step** (Δt, 60 s) of irradiation at the orbit-mean proton spectrum. This is the reference curve R(τ) used by steps 7 and 8.

- **Primaries per simulation energy** (paper Eq. 4, s_j ΔE_j): each simulation energy E_j stands for the band between the geometric midpoints to its neighbours; the first and last energy are the outer edges. N_j = Δt · (F(>e_j) − F(>e_{j+1})) · A_beam, with F the orbit-mean integral flux from SPENVIS.
- **Beam area** A_beam = πR² sin²θmax comes from the source geometry recorded by `0_run.py` (`run_info.json` → `activities.pkl` attrs). `--R`/`--thetamax` override it.
- **Energies** are all those in `activities.pkl`; `--energies` selects a subset.
- **Δt** is the SPENVIS time step, as in `activation_history.py`.
- Activities are summed over energies for each (volume, isotope) and multiplied by its spectrum per decay (paper Eq. 6).
- `--threshold` (default 0) drops (volume, isotope) activities after the single step that are not above it. It is only a speed-up: the dropped contributions add up over many steps in steps 7 and 8.

**Inputs:** `activities.pkl`, SPENVIS file, `result_spectra/*.dat`\
**Outputs:** `spectra.pkl` (with the normalisation in its attrs), `count_rate.dat`, plots

```bash
python accumulate_spectra.py output/activities.pkl AP8MIN.AP8.output_mean_flux_550km_SSO.txt \
    --spectra-dir result_spectra/ \
    --outdir output/ \
    --save-plot count_rate.pdf \
    --save-spectra spectra.pdf
```

**Spectrum file naming convention:**

```
result_spectra/{volume}_{isotope}_S-mode.dat
```

Each file contains 2047 values of \[counts/s/keV/decay\], one per S-mode energy channel defined by:

```python
binning_S = np.arange(18., 4114., 2.)   # edges [keV]
en_s = binning_S[:-1] + np.diff(binning_S / 2.)   # 2047 centres
```

**Output structure:**

`spectra.pkl` — DataFrame indexed by decay time \[s\], columns = energy channels. Values are total background spectrum \[counts/s/keV\] at each time. `attrs`: `timestep_s`, `duty` (1, or the in-belt fraction with `--in-belt-only`), `beam_area_cm2`, `nprim_per_step` (`{E: N_j}`), `in_belt_only`, `activity_threshold`.

`count_rate.dat` — two-column ASCII: `time_s count_rate_cps`.

```python
from accumulate_spectra import load_outputs
spectra_df, count_rate = load_outputs("output/")
spectra_df.loc[100.0]      # spectrum at t = 100 s after irradiation
```

---

### 7. `activation_history.py`

Computes the time history of the activation-induced count rate by convolving the orbital proton flux time series with the reference curve R(τ) of step 6.

**Physical model:** the SPENVIS flux F is constant during each Δt step. The rate at the end of step j is

```
C[j] = Σ_{k=0}^{j}  K_k · F[j-k] / F_mean,     K_k = (1/Δt) ∫_{kΔt}^{(k+1)Δt} R(τ) dτ
```

evaluated via FFT. The step-averaged kernel K_k counts the decays within one step of the irradiation. Sampling R at the lags (k+1)Δt instead would lose them: a nuclide with mean life 9 s (Al26m) would get ~0 instead of its equilibrium rate, and one with 10 min about 95%. In long-term mode the SPENVIS flux is tiled over the requested duration, assuming periodic repetition.

**Out of belt:** a step is in the belt when the total flux (first SPENVIS column) is above `--belt-threshold` (default 0). The running average and the summary rates use the out-of-belt steps only; the all-step mean is printed for comparison.

**Inputs:** `count_rate.dat`, SPENVIS file\
**Outputs:** time-history ASCII (`t_s`, `count_rate_cps`, `out_of_belt_avg_count_rate_cps`), plot

```bash
# Single SPENVIS period
python activation_history.py output/count_rate.dat AP8MIN.AP8.output_mean_flux_550km_SSO.txt \
    --save-plot history.pdf --save-dat history.dat

# Long-term: 3 years, 1-week out-of-belt running average
python activation_history.py output/count_rate.dat AP8MIN.AP8.output_mean_flux_550km_SSO.txt \
    --duration 3y --avg-window 1w \
    --save-plot history_3yr.pdf
```

**Duration/window strings:** `s`, `m`, `h`, `d`, `w`, `mo`, `y`(e.g. `3y`, `18mo`, `2.5w`, `90m`).

---

### 8. `average_spectrum.py`

Computes the steady-state **out-of-belt** background spectrum after a mission of duration T, and identifies the isotopes responsible for the most prominent spectral lines.

**Physical basis:** with the same step-averaged kernel as step 7, K_k(E) = (1/Δt) ∫_{kΔt}^{(k+1)Δt} S(τ, E) dτ, the spectrum averaged over the out-of-belt steps j is

```
⟨S(E)⟩ = Σ_{k < T/Δt} w_k K_k(E),     w_k = mean_{j out of belt} F[j-k] / F_mean
```

The weights w_k are exact (circular correlation of the out-of-belt mask with the flux) for lags within one SPENVIS period, and equal to their period mean for longer lags. `--all-orbit` averages over all steps, which reduces to (⟨F⟩/F_mean/Δt) ∫₀^T S dτ. Nuclides with half-lives much longer than T count with the activity built up to T. Without `--duration`, T is the end of the time grid of step 3 (1e9 s). Δt, the primaries and the flux normalisation come from the attrs of `spectra.pkl`; the SPENVIS file must be the one used in step 6. The total rate agrees with the last-period mean of step 7 to ~0.1%.

Line identification pipeline for each peak:

1. Find the N tallest well-separated peaks (log-scale separation, greedy)
2. Query IAEA LiveChart for nuclides with a known gamma line at that energy
3. Intersect with isotopes present in the simulation
4. Pick the largest contributor among the database-confirmed candidates

**Inputs:** `spectra.pkl`, SPENVIS file\
**Optional:** `activities.pkl`, `result_spectra/` (for line identification)\
**Outputs:** average spectrum ASCII, publication-quality plot

```bash
python average_spectrum.py output/spectra.pkl AP8MIN.AP8.output_mean_flux_550km_SSO.txt \
    --duration 3y \
    --activities output/activities.pkl \
    --spectra-dir result_spectra/ \
    --n-label 6 \
    --window 5 \
    --save-plot avg_spectrum.pdf \
    --save-dat avg_spectrum.dat
```

**Key arguments:**

| Argument | Default | Description |
| --- | --- | --- |
| `--duration` | end of time grid | Mission duration T (`3y`, `18mo`, …) |
| `--all-orbit` | off | Average over all steps instead of the out-of-belt ones |
| `--belt-threshold` | `0` p/cm²/s | Total flux above which a step is in the belt |
| `--activities` | — | For per-isotope line identification |
| `--spectra-dir` | — | For per-isotope line identification |
| `--n-label` | `5` | Number of peaks to label |
| `--prominence` | `0.5` | Min log₁₀-prominence for peak detection |
| `--min-sep` | `0.10` dex | Min peak separation in log₁₀(E) decades |
| `--window` | `5` keV | ±energy window for DB query |
| `--cache` | `.gamma_line_cache.pkl` | Cache for DB responses |

---

## Shared geometry parameters

The source geometry and the simulation energies are recorded once, by `0_run.py` in `run_info.json`, and travel with the data: `results.pkl` and `activities.pkl` carry them in their attrs, and `spectra.pkl` carries the resulting normalisation (Δt, beam area, primaries per step). `accumulate_spectra.py` and `average_spectrum.py` read them from there; only the stand-alone `spenvis_parser.py` CLI still takes `--R`, `--thetamax`, `--timestep` and `--energies` (its defaults are those of an older 50 m source: pass `--R 13 --thetamax 90` for CUSP).

`--in-belt-only` must be the same in `accumulate_spectra.py` and `activation_history.py` (`average_spectrum.py` reads it from `spectra.pkl`).

---

## Common file formats

| File | Format | Producer | Consumer(s) |
| --- | --- | --- | --- |
| `results.pkl` | pandas pickle | `activation_parser` | `build_decay_chains`, `compute_activities` |
| `DecayChains/*.dat` | ASCII | `decay_chain_builder` / `build_decay_chains` | `compute_activities` |
| `activities.pkl` | pandas pickle | `compute_activities` | `accumulate_spectra`, `average_spectrum` |
| `active_isotopes.pkl` | Python pickle (dict) | `compute_activities` | (reference) |
| `spectra.pkl` | pandas pickle | `accumulate_spectra` | `average_spectrum` |
| `count_rate.dat` | 2-col ASCII | `accumulate_spectra` | `activation_history`, `average_spectrum` |
| `result_spectra/*.dat` | ASCII (2047 values) | Geant4 simulation | `accumulate_spectra`, `average_spectrum` |
| `.gamma_line_cache.pkl` | Python pickle | `average_spectrum` | `average_spectrum` (cache) |

---

## Complete example run

```bash
SPENVIS=AP8MIN.AP8.output_mean_flux_550km_SSO.txt

# 0. Run Geant4
python 0_run.py

# 1. Parse Geant4 score files (primaries from output/run_info.json)
python activation_parser.py --dir output/ --save-pkl results.pkl

# 2. Build decay chains for all activated isotopes
python build_decay_chains.py results.pkl --outdir DecayChains/

# 3. Solve radioactive decay equations
python compute_activities.py results.pkl --chains DecayChains/ --outdir output/ --errors

# 4. Geant4 post-activation run (scripts not yet in the repository)
#    -> result_spectra/{vol}_{iso}_S-mode.dat

# 5. Inspect the SPENVIS orbital flux (optional)
python spenvis_parser.py $SPENVIS --R 13 --thetamax 90

# 6. Spectra and count rate after one step at the mean flux
python accumulate_spectra.py output/activities.pkl $SPENVIS \
    --spectra-dir result_spectra/ --outdir output/

# 7a. Rate history over one SPENVIS period
python activation_history.py output/count_rate.dat $SPENVIS \
    --save-plot history_orbit.pdf

# 7b. Long-term (3-year) history, out-of-belt running average
python activation_history.py output/count_rate.dat $SPENVIS \
    --duration 3y --avg-window 1w --save-plot history_3yr.pdf

# 8. Steady-state out-of-belt spectrum with line identification
python average_spectrum.py output/spectra.pkl $SPENVIS --duration 3y \
    --activities output/activities.pkl \
    --spectra-dir result_spectra/ \
    --save-plot avg_spectrum.pdf \
    --save-dat avg_spectrum.dat
```

If the SPENVIS file tabulates the flux only up to an energy below the highest simulation energy (AP8: 400 MeV), step 6 clips the bands there and warns: protons above it are ignored.

---

## Dependencies

```
numpy, scipy, pandas, matplotlib, mpmath
```

Python ≥ 3.11 (`0_run.py` uses `tomllib`). `average_spectrum.py` queries IAEA LiveChart with `urllib` (requires internet access; results are cached locally after the first run). `myUtilities.prettifyPlot` is used for the plot style when it is importable, and skipped otherwise.