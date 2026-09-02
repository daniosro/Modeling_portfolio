# Process-Based Modeling in Ocean and Climate Science

## Interactive Book

**[View the interactive Jupyter Book](https://daniosro.github.io/Modeling_portfolio/)**

The book includes fully interactive Plotly figures (hover for exact values, zoom, and toggle series) for all Python notebooks.

---

## Models

### 1. Accelerated Weathering of Limestone Ship Reactor — Dong et al. (2025)

**[View notebook](https://daniosro.github.io/Modeling_portfolio/notebooks/awl_reactor_model.html)**

A box model of the Accelerated Weathering of Limestone (AWL) reactor (Dong et al., *Science Advances*, 2025). The reactor dissolves limestone in CO₂-rich ship flue gas and seawater, converting CO₂ to stable bicarbonate ions, evaluating the potential for permanent ocean carbon sequestration without pipelines or injection wells.

**Implements:** Eqs. 4-11 from the paper (two-segment calcite dissolution kinetics, Henry's law CO₂ saturation, and mass balance ODEs for a 4-unit packed column stirred tank reactor in parallel and counterflow configurations).

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

A proteome-allocation optimization model for colonial *Trichodesmium*, one of the main nitrogen-fixing marine cyanobacterium. Predicts how cells optimally allocate protein investment across all metabolic functions to maximize growth rate under a range of temperatures, with phosphorus supplied by the colonial microbiome as dissolved organic phosphorus (DOP).

From Osorio-Rodriguez et al. (in review) — "Microbiome-mediated nutrient acquisition and cell morphology jointly shape the niche and warming responses of N₂-fixing cyanobacteria."

**Implements:** Nonlinear programming (NLP) optimization via JuMP/Ipopt. Maximizes log growth rate subject to mass balance, energetic, stoichiometric, membrane, and density constraints across an 18–35°C temperature sweep.

**Language:** Julia (JuMP, Ipopt, DataFrames, Plots)

### 4. Physics-Informed Machine Learning Trichodesmium Surrogate Model

A physics-informed neural network (PIML) surrogate trained on outputs from the Trichodesmium optimization model. Standard neural networks trained on biological data can predict impossible states, such as negative growth rates, or proteome fractions that exceed the cell’s total protein budget. Physics-informed machine learning fixes this by embedding known biological constraints directly into the training loss, penalizing violations of the proteome budget, Arrhenius temperature scaling, and growth rate bounds during network training.

The surrogate predicts growth rate, proteome allocation fractions, and elemental ratios from temperature and DOP inputs in milliseconds (roughly 10,000x faster than the full Julia optimization) while remaining physically consistent.

Implements: Two-network comparison (standard neural network vs physics-informed network) with interactive Plotly figures showing training curves, prediction accuracy, and constraint satisfaction. Includes a response surface heatmap of the full temperature x DOP parameter space.

Connection to real-time monitoring: The same PIML approach applies directly to wastewater treatment and carbon sequestration platforms: train a surrogate on process model outputs, enforce mass balance and kinetic constraints in the loss, deploy for real-time prediction from sensor streams.

Language: Python (PyTorch, NumPy, SciPy, Plotly)

---

## Mathematical framework

All models share the same fundamental approach: **process-based modeling with empirical kinetic rate laws**.

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
