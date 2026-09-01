# Process-Based Modeling in Ocean and Climate Science

**Daniela Osorio Rodriguez, PhD**  
Postdoctoral Research Fellow, University of Southern California  
[danieosro@gmail.com](mailto:danieosro@gmail.com) | [github.com/daniosro](https://github.com/daniosro)

---

## About this book

This Jupyter Book documents three process-based computational models at the intersection of ocean science, climate, and biotechnology. Each model encodes known physical, chemical, or biological mechanisms as differential equations or optimization constraints, and uses those equations to predict system behavior under conditions that have not been directly measured.

All models are implemented in open, reproducible formats with documented code, described parameters, and executable outputs.

---

## Models

### 1. AWL Ship Reactor — Dong et al. (2025)

A box model of the **Accelerated Weathering of Limestone (AWL)** reactor developed by Dong et al. (*Science Advances*, 2025, eadr7250). The reactor dissolves limestone in CO2-rich seawater aboard cargo ships, converting flue gas CO2 to stable bicarbonate ions dissolved in the ocean — permanent carbon sequestration without pipelines or injection wells.

**Implements:** Eqs. 4-11 from the paper — calcite dissolution kinetics, Henry's law CO2 saturation, and mass balance ODEs for a 4-unit packed column stirred tank reactor in parallel and counterflow configurations.

**Language:** Python (NumPy, SciPy, Matplotlib)

---

### 2. Coupled Wastewater Treatment and Carbon Mineralization — CREW

A coupled model of biological wastewater treatment (Activated Sludge Model No. 1, Henze et al. 1987) and inorganic carbon mineralization, inspired by [CREW's](https://www.crewcarbon.com) process intensification technology. CREW makes wastewater treatment cheaper and carbon-negative by reacting biogenic CO2 with alkaline minerals, permanently sequestering it as dissolved bicarbonate.

**Implements:** ASM1 mass balance ODEs for COD, nitrogen, oxygen, and biomass dynamics coupled to a first-order CO2 mineralization rate law. Demonstrates the synergistic coupling between biological CO2 production and mineral dissolution.

**Language:** Python (NumPy, SciPy, Matplotlib, Plotly)

---

### 3. Colonial *Trichodesmium* Optimization Model

A proteome-allocation optimization model for colonial *Trichodesmium* — a nitrogen-fixing, filamentous marine cyanobacterium critical to ocean biogeochemical cycling. The model predicts how the cell optimally allocates its protein investment across all metabolic functions to maximize growth rate under a range of temperatures, with phosphorus supplied by the colonial microbiome as dissolved organic phosphorus (DOP).

**Scientific context:** From Osorio-Rodriguez et al. (2026, *ISME Journal*) — "Microbiome-mediated nutrient acquisition and cell morphology jointly shape the niche and warming responses of N2-fixing cyanobacteria." The model examines how the *Trichodesmium* microbiome impacts growth under oligotrophic conditions and how thermal optima and cell morphology determine biogeographic distributions.

**Implements:** Nonlinear programming (NLP) optimization via JuMP/Ipopt in Julia. Maximizes log growth rate subject to mass balance, energetic, stoichiometric, membrane, and density constraints across a 18-35°C temperature sweep.

**Language:** Julia (JuMP, Ipopt, DataFrames, Plots)

---

## Mathematical framework

All three models share the same fundamental approach: **process-based modeling with empirical kinetic rate laws**.

Each model encodes known mechanisms — enzyme kinetics, dissolution rate laws, mass balance constraints — as mathematical equations, then solves those equations to predict system behavior. This approach is more interpretable than black-box machine learning models: every parameter has a physical meaning, every equation represents a known process, and the model's predictions can be traced back to specific assumptions.

The key modeling elements shared across all three notebooks:

| Element | AWL Reactor | CREW Model | *Trichodesmium* |
|---------|-------------|------------|-----------------|
| Kinetics | Calcite dissolution rate law | Monod + first-order | Michaelis-Menten |
| Mass balance | ODE per reactor unit | ODE per state variable | Algebraic constraints |
| Solver | SciPy LSODA | SciPy LSODA | JuMP/Ipopt NLP |
| Key parameter | Solid holdup, grain size | HRT, SRT, mineral loading | Temperature, nutrient availability |
| Output | Alkalinity, pH, efficiency | COD removal, CO2 sequestration | Growth rate, proteome allocation |

---

## Reproducibility

All code is available on [GitHub](https://github.com/daniosro). Each notebook documents:
- What the model does in plain language
- What each parameter represents and its units
- What assumptions are made and where they may fail
- How to run the code and interpret the output

Python notebooks require: `numpy`, `scipy`, `matplotlib`, `plotly`, `pandas`  
Julia notebooks require: `JuMP`, `Ipopt`, `DataFrames`, `Plots`, `CSV`
