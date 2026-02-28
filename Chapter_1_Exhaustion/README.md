# Mathematical Model of Killer Cell Exhaustion and Target Cell Response

> **Question to be answered**:
> - Q1. How the exhaustion of killer immune cell influence the killing output? This can be answered by the analytical solution of the ODEs (see: [Analytical Jupyter Notebook](./Analytical_Solution.ipynb))
> - Q2. Within the context of cell migration and cell adhesion, what is the influence of exhaustion to the killing output? Here are two 
scenario. 
>   - the same adhesion for all cell-cell interactions: we can know what is the pure effects from cell migration
>   - stronger ligands-receptor binding has stronger adhesion parameter: combination effects.
> - Q3. Any optimised combination of this system? like two kinds of ligand-receptor binding, sort of stuff. What is the antibody's role in this story?

## 1. The Model

### Cytotoxic interactions
1. **Killer Cell Cytotoxicity ($ \Omega $)**: The capacity of a killer cell (e.g., a Cytotoxic T Lymphocyte or a Natural Killer cell) to induce death factor accumulation in a target cell.
    $$
    \frac{d\Omega_{t}}{dt} = - \tilde{\rho} \Omega_{t}
    $$
    The exhaustion of the killer cell is modelled as an exponential decay of its cytotoxic ability, $\Omega_{t}$, at a rate $\tilde{\rho}$.

2.  **Target Cell Death Factor ($ \chi $)**
    $$
    \frac{d\chi_{t}}{dt} = \lambda \Omega_{t} - \phi \chi_{t}
    $$
    This equation describes the change in the target cell's death factor, $\chi_{t}$. The term $\Delta \Omega_{t}$ represents the accumulation of the death factor, which is directly proportional to the killer cell's current cytotoxicity. The term $-\phi \chi_{t}$ represents the natural clearance or repair of this factor by the target cell at a rate $\phi$.

3. **The ligand-receptor binding strength ($\lambda$)**
    $$
     \frac{\tilde{\rho} - \rho_{min}}{\rho_{max} - \rho_{min}} = \lambda\Rightarrow \tilde{\rho} = \rho_{min} + \lambda \cdot (\rho_{max} + \rho_{min})
    $$
---
- Variables: $\Omega$ & $\chi$
- Parameters: $\lambda$, $\phi$, $\rho_{min}$, $\rho_{max}$, $\Omega_{0}$, & $\chi_{0}$
---


### Cell Migration and Adhesion

We model two-dimensional (2D) cell motility using the Ornstein-Uhlenbeck (OU) process. Each cell is characterized by a time-dependent polarity vector, $\vec{\mu}(t)$, which dictates its intrinsic direction and persistence of motion. The dynamics of this vector are governed by the stochastic differential equation:
$$
\frac{d\vec{\mu}(t)}{dt} = -\gamma \vec{\mu}(t) + g \vec{\eta}_{\mu}
$$
where:
- The term $-\gamma \vec{\mu}(t)$ represents the relaxation of cell polarity towards zero at a rate $\gamma$.
- $g \vec{\eta}_{\mu}$ is a white Gaussian noise term with intensity $g$, introducing stochasticity.

The cell's velocity is determined by its internal polarity, mechanical forces from neighboring cells, and a random stochastic component. The equation for the cell's position, $\vec{r}(t)$, is:
$$
\frac{d\vec{r}(t)}{dt} = \beta \vec{\mu}(t) + \frac{1}{\lambda} \sum_{j} \vec{F}_{LJ}(\vec{r}_{ij}) + b \vec{\eta}_r
$$
where:
- $\beta \vec{\mu}(t)$ is the velocity component driven by cell polarity.
- $\frac{1}{\lambda} \sum_{j} \vec{F}_{LJ}(\vec{r}_{ij})$ is the velocity resulting from the sum of Lennard-Jones (LJ) forces from neighboring cells.
- $b \vec{\eta}_r$ is a white Gaussian noise term representing random fluctuations.

Mechanical interactions, such as adhesion and repulsion, are modeled using a pairwise Lennard-Jones (LJ) force. This force is a function of the distance between cell centers and is crucial for simulating cell-cell contact, synapse formation, and preventing cellular overlap. The LJ force between two cells is given by:
$$
    \vec{F}_{LJ}(\vec{r}) = 24\epsilon \left[ 2\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right] \frac{\hat{r}}{r}
$$
where:
- $r$ is the center-to-center distance.
- $\hat{r}$ is the unit vector pointing between the cells.
- $\epsilon$ is the interaction strength (depth of the potential well), representing adhesion.
- $\sigma$ is the interaction distance (finite distance at which the inter-particle potential is zero).

The parameters $\epsilon$ and $\sigma$ are specific to the types of interacting cells (e.g., killer-killer, target-target, or killer-target pairs).


## 2. Simulation results

### 2-1. Analytical solutions of the cytotoxic interaction ODEs (w/o ABM).
