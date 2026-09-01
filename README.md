# Process-Based Modeling in Ocean and Climate Science

**Daniela Osorio Rodriguez, PhD**  
Postdoctoral Research Fellow, University of Southern California  
[danieosro@gmail.com](mailto:danieosro@gmail.com)

---

## Interactive Book

**[View the interactive Jupyter Book](https://daniosro.github.io/Modeling_portfolio/)**

The book includes fully interactive Plotly figures — hover for exact values, zoom, and toggle series — for all Python notebooks.

---

## Models

### 1. AWL Ship Reactor — Dong et al. (2025)

**[View notebook](https://daniosro.github.io/Modeling_portfolio/notebooks/awl_reactor_model.html)**

A box model of the Accelerated Weathering of Limestone (AWL) reactor (Dong et al., *Science Advances*, 2025, eadr7250). The reactor dissolves limestone in CO₂-rich ship flue gas and seawater, converting CO₂ to stable bicarbonate ions — permanent ocean carbon sequestration without pipelines or injection wells.

**Implements:** Eqs. 4-11 from the paper — two-segment calcite dissolution kinetics, Henry's law CO₂ saturation, and mass balance ODEs for a 4-unit packed column stirred tank reactor in parallel and counterflow configurations.

**Language:** Python (NumPy, SciPy, Matplotlib, Plotly)

---

### 2. Coupled Wastewater Treatment and Carbon Mineralization

**[View notebook](https://daniosro.github.io/Modeling_portfolio/notebooks/crew_process_model.html)**

A coupled model of biological wastewater treatment (Activated Sludge Model No. 1, Henze et al. 1987) and inorganic carbon mineralization, inspired by [CREW's](https://www.crewcarbon.com) process intensification technology. Demonstrates the synergistic coupling between biological CO₂ production and mineral dissolution for permanent carbon sequestration.

**Implements:** ASM1 mass balance ODEs for COD, nitrogen, oxygen, and biomass dynamics coupled to a first-order CO₂ mineralization rate law.

**Language:** Python (NumPy, SciPy, Matplotlib, Plotly)

---

### 3. Colonial *Trichodesmium* Optimization Model

**[View notebook](https://daniosro.github.io/Modeling_portfolio/notebooks/trichodesmium_temperature_limitedPO4.html)**

A proteome-allocation optimization model for colonial *Trichodesmium* — a nitrogen-fixing marine cyanobacterium central to ocean biogeochemical cycling. Predicts how cells optimally allocate protein investment across all metabolic functions to maximize growth rate under a range of temperatures, with phosphorus supplied by the colonial microbiome as dissolved organic phosphorus (DOP).

From Osorio-Rodriguez et al. (in review) — "Microbiome-mediated nutrient acquisition and cell morphology jointly shape the niche and warming responses of N₂-fixing cyanobacteria."

**Implements:** Nonlinear programming (NLP) optimization via JuMP/Ipopt. Maximizes log growth rate subject to mass balance, energetic, stoichiometric, membrane, and density constraints across an 18–35°C temperature sweep.

**Language:** Julia (JuMP, Ipopt, DataFrames, Plots)

---

## Mathematical framework

All three models share the same core approach: **process-based modeling with empirical kinetic rate laws** — encoding known physical, chemical, or biological mechanisms as differential equations or optimization constraints, then solving them to predict system behavior.

| | AWL Reactor | CREW Model | *Trichodesmium* |
|--|--|--|--|
| Kinetics | Calcite dissolution rate law | Monod + first-order | Michaelis-Menten |
| Math | ODE per reactor unit | ODE per state variable | NLP constraints |
| Solver | SciPy LSODA | SciPy LSODA | JuMP/Ipopt |
| Key parameter | Solid holdup, grain size | HRT, SRT, mineral loading | Temperature, nutrient supply |
| Output | Alkalinity, pH, efficiency | COD removal, CO₂ sequestration | Growth rate, proteome allocation |

---

## Running locally

**Python notebooks (AWL and CREW):**
```bash
pip install numpy scipy matplotlib plotly pandas jupyterlab
jupyter lab
```

**Julia notebook (Trichodesmium):**
Requires Julia with JuMP, Ipopt, DataFrames, Plots, and CSV packages.

---

## References

- Dong et al. (2025). Accelerated weathering of limestone on cargo ships for ocean carbon dioxide removal. *Science Advances*, eadr7250.
- Henze et al. (1987). Activated Sludge Model No. 1. IWA Publishing.
- Naviaux et al. (2019). Calcite dissolution rates in seawater. *Geochimica et Cosmochimica Acta*, 246:363–384.
- Osorio-Rodriguez et al. (in review). Microbiome-mediated nutrient acquisition and cell morphology jointly shape the niche and warming responses of N₂-fixing cyanobacteria.
