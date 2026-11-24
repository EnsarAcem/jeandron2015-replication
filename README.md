# jeandron2015-replication

# Computational reproduction and robustness checks of Jeandron et al. (2015)

This repository contains the code and reproduced figures for a computational reproduction and robustness reproduction of:
Jeandron, A., Saidi, J. M., Kapama, A., Burhole, M., Birembano, F., Vandevelde, T., Gasparrini, A., Armstrong, B., Cairncross, S., & Ensink, J. H. J. (2015). Water supply interruptions and suspected cholera incidence: A time-series regression in the Democratic Republic of the Congo. *PLoS Medicine, 12*(10), e1001893. https://doi.org/10.1371/journal.pmed.1001893

The work was done for a Replication Game of the Institute for Replication (I4R) by Ensar Acem (Kadir Has University).

## Files
- `pmed.1001893.s004.r` – the authors' R script, with our corrections (missing `attrdl()` function, Models A and E)
- `robustness_checks.R` – robustness checks and Figure 6
- `figure1.tif` to `figure6.tif` – reproduced figures
- Data: download S1 Table from the paper’s Supporting Information

## How to run
1. Put the data file in the R working directory and install the packages `splines`, `tsModel`, `dlnm`, `gplots`, `mvmeta`, `MASS` and `FluMoDL` (archived; install from source and set the path in the script).
2. Run `pmed.1001893.s004.r`.
3. In the same R session, run `robustness_checks.R`.

Analyses were run in R 4.5.1.
