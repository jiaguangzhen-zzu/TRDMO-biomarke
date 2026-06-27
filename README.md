# TRDMO: Temporal Knowledge Transfer Dynamic Multi-objective Optimization
Source MATLAB codes and datasets for critical disease biomarker identification based on dynamic gene networks.

## File Overview
### MATLAB Code Files
- `TRDMO.m`: Main execution script of the proposed algorithm
- Baseline comparison methods: `DNSGA2.m`, `MOEAPSL.m`, `RBM.m`
- Evolutionary operators: `creatpop.m`, `Cross1.m`, `Muation.m`, `TournamentSelection.m`, `NDSort.m`, `CrowdingDistance.m`
- Performance metrics: `HV.m`, `IGD.m`
- Network construction & training: `calculate_second_order_connectivity.m`, `DAE.m`, `ModelTraining.m`
- Result processing & plotting: `extract.m`, `Fin_POP.m`, `FinPOP_fig.m`

### /data Folder
TCGA multi-cancer gene expression profiles and dynamic gene interaction networks for experiments.

## Running Environment
- MATLAB R2021b or newer
- Required Toolboxes: Global Optimization Toolbox, Bioinformatics Toolbox

## Quick Start
1. Set MATLAB working directory to repository root
2. Run `TRDMO.m` to reproduce the proposed method
3. Execute baseline scripts to compare existing algorithms

## Citation
If this code helps your research, cite our paper:
```bibtex
@article{trdmo2026biomarker,
  title={Temporal Knowledge Transfer Dynamic Multi-objective Optimization for Critical Biomarker Identification},
  author={XXX},
  journal={XXX},
  year={2026}
}
