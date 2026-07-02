# Smart eDNA sampling strategies — analysis code

Analysis code for *Smart eDNA sampling strategies*, evaluating how sampling
effort affects the accuracy of Gaussian Process (GP) predictions and spatial
parameter recovery for environmental DNA (eDNA) surveys. The analysis combines
empirical eDNA data from 12 fish species collected along the U.S. West Coast
(Guri et al. 2024, *ICES Journal of Marine Science*; Guri et al. 2025) with
simulated spatial fields covering a range of spatial autocorrelation
strengths, and progressively thins both to quantify how many samples are
needed to reliably estimate spatial structure and predict eDNA concentration.

## Repository structure

```
Code/            # analysis scripts (this README)
Data/            # input data — pulled from Zenodo, see "Data & Output" below
Output/          # simulation and model output — pulled from Zenodo, see below
Plots/           # manuscript figures produced by 03_Plots.R
```

`Data/` and `Output/` are not part of this repository — they are hosted on
Zenodo (see below) and are expected as siblings of `Code/` once downloaded.
`Plots/` must exist before running `03_Plots.R`, which saves the manuscript's
main and supplementary figures there (`ggsave()` does not create missing
folders).

## Pipeline

Scripts are numbered in the order they should be run:

| Script | Purpose |
|---|---|
| `00_Download_data.R` | Downloads the `Data/` and `Output/` folders from the Zenodo archive (see below) and unpacks them as siblings of `Code/`. Run this once, before anything else. |
| `01_GP_sim.R` | **Simulations.** Simulates raw spatial GP fields at the empirical prediction grid over a range of length-scales (`Output/Raw_GP_fileds_simulated/`), then progressively thins each field down to a series of sample sizes *N* and fits the Stan GP model (`GP.stan`) to each thinned sample, saving the thinned samples, estimated GP parameters (`μ`, `α`, `ρ`, `σ`), and grid-wide predictions to `Output/GP_<N>/`. The same thin-and-fit procedure is repeated on the empirical multi-species eDNA data (`Output/GP_multifish/`). |
| `02_Data_manip.R` | **Data manipulation.** Collects the many per-iteration/per-species files produced by `01_GP_sim.R` from `Output/GP_<N>/`, `Output/Raw_GP_fileds_simulated/`, `Output/simulated_parameters/`, and `Output/GP_multifish/`, and compiles them into a small number of compact combined data frames (`.rds`/`.csv`) saved to `Output/Complete_db/`. |
| `03_Plots.R` | Reads the compiled data frames in `Output/Complete_db/` and produces the manuscript's main and supplementary figures, saved to `Plots/`. |

`GP.stan` is the Stan model fitted by `01_GP_sim.R`. `GP.rds` is an
auto-generated compiled-model cache created by `rstan` (via
`rstan_options(auto_write = TRUE)`) the first time `GP.stan` is compiled — it
does not need to be run or edited directly, and will be regenerated
automatically if deleted.

## Data & Output (Zenodo)

`Data/` and `Output/` are archived on Zenodo:

> **Zenodo record:** `[10.5281/zenodo.21141367]`

Run `00_Download_data.R` to fetch and unpack them automatically, or download
the archive manually from the Zenodo record above and place the `Data/` and
`Output/` folders as siblings of `Code/`. See the `README.md` included in
that archive for a full description of the folder/file structure.

## Requirements

- R (≥ 4.3)
- R packages: `dplyr`, `tidyr`, `tibble`, `purrr`, `readr`, `stringr`,
  `ggplot2`, `PNWColors`, `broom`, `here`, `fields`, `MASS`,
  `spatstat.geom`, `rstan`
- A working Stan toolchain (via `rstan`) to compile and sample `GP.stan`

`01_GP_sim.R` fits several thousand Stan models (one per depth stratum ×
length-scale × iteration × sample size / species) and was run on an HPC
cluster; expect it to be computationally intensive on a laptop.

## Citation

If you use this code or data, please cite the associated manuscript and the
source of the empirical eDNA data:

> Guri, G. et al. Smart eDNA sampling strategies. *In review.*
>
> Guri, G., Shelton, A.O., Kelly, R.P., et al. (2024). Predicting trawl
> catches using environmental DNA. *ICES Journal of Marine Science*, 81(8),
> 1536–1548. https://doi.org/10.1093/icesjms/fsae097
