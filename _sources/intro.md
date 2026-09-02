# Process-Based Modeling in Ocean and Climate Science
---

## About this book

This Jupyter Book documents three process-based computational models at the intersection of ocean science, climate, and biotechnology. Each model encodes known physical, chemical, or biological mechanisms as differential equations or optimization constraints, and uses those equations to predict system behavior under conditions that have not been directly measured.

All models are implemented in open, reproducible formats with documented code, described parameters, and executable outputs.

---

## Models

### 1. Accelerated Weathering of Limestone (AWL) Ship Reactor — Dong et al. (2025)

A box model of the **Accelerated Weathering of Limestone (AWL)** reactor developed by Dong et al. (*Science Advances*, 2025). The reactor dissolves limestone in CO2-rich seawater aboard cargo ships, converting flue gas CO2 to stable bicarbonate ions dissolved in the ocean. The model assesses the feasibility of the process for permanent carbon sequestration without pipelines or injection wells.

**Implements:** Eqs. 4-11 from the paper (calcite dissolution kinetics, Henry's law CO2 saturation, and mass balance ODEs for a 4-unit packed column stirred tank reactor in parallel and counterflow configurations).

**Language:** Python (NumPy, SciPy, Matplotlib)

---

### 2. Coupled Wastewater Treatment and Carbon Mineralization — CREW

A coupled model of biological wastewater treatment (Activated Sludge Model No. 1, Henze et al. 1987) and inorganic carbon mineralization, inspired by [CREW's](https://www.crewcarbon.com) process intensification technology. CREW treats wastewater by reacting biogenic CO2 with alkaline minerals for permanent carbon sequestration as dissolved bicarbonate.

**Implements:** ASM1 mass balance ODEs for COD, nitrogen, oxygen, and biomass dynamics coupled to a first-order CO2 mineralization rate law. Demonstrates the synergistic coupling between biological CO2 production and mineral dissolution.

**Language:** Python (NumPy, SciPy, Matplotlib, Plotly)

---

### 3. Single-cell proteome allocation model for Colonial *Trichodesmium* 

This notebook implements the proteome-allocation optimization model for colonial Trichodesmium described in Osorio-Rodriguez et al. (2026, in prep.). The model predicts how a colonial Trichodesmium cell optimally allocates its proteome (the total protein investment across all metabolic functions) to maximize growth rate under a range of temperatures, with phosphorus supplied exclusively as inorganic phosphate (PO₄) and phycosphere-related phosphorus (PRP, here called DOP).

**Scientific context:** From Osorio-Rodriguez et al. (2026, *in prep.*) — "Microbiome-mediated nutrient acquisition and cell morphology jointly shape the niche and warming responses of N2-fixing cyanobacteria." The model examines how the *Trichodesmium* microbiome impacts growth under oligotrophic conditions and how thermal optima and cell morphology determine biogeographic distributions.

**Implements:** Nonlinear programming (NLP) optimization via JuMP/Ipopt in Julia. Maximizes log growth rate subject to mass balance, energetic, stoichiometric, membrane, and density constraints across a 18-35°C temperature sweep.

**Language:** Julia (JuMP, Ipopt, DataFrames, Plots)

---

### 4. Physics-Informed Machine Learning *Trichodesmium* Surrogate Model

A physics-informed neural network (PIML) surrogate trained on outputs from the *Trichodesmium* optimization model. Standard neural networks trained on biological data can predict impossible states, such as negative growth rates, or proteome fractions that exceed the cell's total protein budget. Physics-informed machine learning fixes this by embedding known biological constraints directly into the training loss, penalizing violations of the proteome budget, Arrhenius temperature scaling, and growth rate bounds during network training.

The surrogate predicts growth rate, proteome allocation fractions, and elemental ratios from temperature and DOP inputs in milliseconds (roughly 10,000x faster than the full Julia optimization) while remaining physically consistent. 

**Implements:** Two-network comparison  (standard neural network vs physics-informed network) with interactive Plotly figures showing training curves, prediction accuracy, and constraint satisfaction. Includes a response surface heatmap of the full temperature x DOP parameter space.

**Connection to real-time monitoring:** The same PIML approach applies directly to wastewater treatment and carbon sequestration platforms: train a surrogate on process model outputs, enforce mass balance and kinetic constraints in the loss, deploy for real-time prediction from sensor streams.

**Language:** Python (PyTorch, NumPy, SciPy, Plotly)

---

## Mathematical framework

All three models share the same fundamental approach: **process-based modeling with empirical kinetic rate laws**.

Each model encodes known mechanisms, such as enzyme kinetics, dissolution rate laws, and mass balance constraints as mathematical equations, then solves those equations to predict system behavior. This approach is more interpretable than black-box machine learning models: every parameter has a physical meaning, every equation represents a known process, and the model's predictions can be traced back to specific assumptions.

The key modeling elements shared across all notebooks:

| Element | AWL Reactor | CREW Model | *Trichodesmium* | *Trichodesmium* ML |
|---------|-------------|------------|-----------------|
| Kinetics | Calcite dissolution rate law | Monod + first-order | Michaelis-Menten | Neural network + physics penalties |
| Mass balance | ODE per reactor unit | ODE per state variable | Algebraic constraints | Backpropagation + constraint loss |
| Solver | SciPy LSODA | SciPy LSODA | JuMP/Ipopt NLP | PyTorch Adam |
| Key parameter | Solid holdup, grain size | HRT, SRT, mineral loading | Temperature, nutrient availability | Penalty weight lambda |
| Output | Alkalinity, pH, efficiency | COD removal, CO2 sequestration | Growth rate, proteome allocation | Growth rate, proteome fractions (real-time) |

---

## Reproducibility

All code is available on [GitHub](https://github.com/daniosro). Each notebook documents:
- What the model does in plain language
- What each parameter represents and its units
- What assumptions are made and where they may fail
- How to run the code and interpret the output

Python notebooks require: `numpy`, `scipy`, `matplotlib`, `plotly`, `pandas`  
Julia notebooks require: `JuMP`, `Ipopt`, `DataFrames`, `Plots`, `CSV`
