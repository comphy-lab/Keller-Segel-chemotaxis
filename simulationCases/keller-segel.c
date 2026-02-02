/**
# Keller-Segel Chemotaxis Model

The Keller-Segel model describes chemotaxis: the directed movement of cells
in response to a chemical gradient. This is a fundamental process in biology,
including embryonic development, wound healing, and immune response.

## Physical Setup

The model couples cell density $\rho$ and chemoattractant concentration $c$:
$$
\partial_t \rho = \nabla^2 \rho - \chi \nabla \cdot (\rho \nabla c)
$$
$$
\partial_t c = D \nabla^2 c + \alpha \rho - \beta c
$$

where:
- $\chi$: chemotactic sensitivity
- $D$: chemoattractant diffusion coefficient
- $\alpha$: production rate
- $\beta$: degradation rate

## Implementation

We use a Cartesian (multi)grid, the generic time loop, and the
time-implicit diffusion solver from Basilisk.

## Author

Vatsal Sanjay  
Email: vatsalsy@comphy-lab.org  
CoMPhy Lab  
Last updated: Jan 30, 2026

## Note

This file is currently configured as a Brusselator example and requires
adaptation to implement the actual Keller-Segel chemotaxis equations.
*/

#include "grid/multigrid.h"
#include "run.h"
#include "diffusion.h"

/**
## Variables

We define scalar fields for the concentrations `C1` and `C2`.
In the Keller-Segel model these would represent cell density $\rho$ and
chemoattractant concentration $c$ respectively. */

scalar C1[], C2[];

/**
## Parameters

Placeholder parameters from [Pena and Perez-Garcia, 2001](/src/references.bib#pena2001).
These need to be adapted for the Keller-Segel model:

- `k`: Reaction/chemotaxis rate (placeholder: 1.0)
- `ka`: Production parameter (placeholder: 4.5)
- `D`: Diffusion coefficient (placeholder: 8.0)
- `mu`: Control parameter
- `kb`: Derived parameter
*/

double k = 1., ka = 4.5, D = 8.;
double mu, kb;

/**
The generic time loop requires a timestep `dt`. We store the statistics
of the diffusion solvers in `mgd1` and `mgd2` for monitoring convergence. */

double dt;
mgstats mgd1, mgd2;

/**
### main()

Main simulation driver (currently configured as Brusselator example).

We configure:
- Grid resolution: 128 × 128
- Domain size: 64 × 64
- Diffusion solver tolerance: 1e-4

TODO: Adapt parameter sweep for Keller-Segel chemotaxis regimes. */

int main()
{
  init_grid (128);
  size (64);
  TOLERANCE = 1e-4;

  /**
  Run three cases (currently Brusselator parameters):
  - $\mu = 0.04$
  - $\mu = 0.1$
  - $\mu = 0.98$
  */

  mu = 0.04; run();
  mu = 0.1;  run();
  mu = 0.98; run();
}

/**
## Initial Conditions

### event init()

Initialize concentration fields (currently using Brusselator setup).

For Brusselator physics:
- Calculate $\nu = \sqrt{1/D}$: characteristic wavenumber
- Compute $k_b^{crit} = (1 + ka \cdot \nu)^2$: critical bifurcation parameter
- Set $kb = k_b^{crit}(1 + \mu)$ based on control parameter $\mu$

TODO: Replace with Keller-Segel initial conditions (e.g., localized cell density,
uniform chemoattractant). */

event init (i = 0)
{
  double nu = sqrt(1./D);
  double kbcrit = sq(1. + ka*nu);
  kb = kbcrit*(1. + mu);

  /**
  Initialize near stationary solution $C_1 = ka$, $C_2 = kb/ka$ with
  random perturbation in $[-0.01, 0.01]$. */

  foreach() {
    C1[] = ka ; 
    C2[] = kb/ka + 0.01*noise();
  }
}

/**
## Outputs

### event movie()

Generate animation frames showing the evolution of $C_1$ concentration.

We output PPM images every 10 iterations for video generation. The
`spread` parameter sets the color scale to $\pm$ twice the standard
deviation. Progress information (iteration, time, timestep, and solver
iterations) is printed to stderr for monitoring.

TODO: Adapt visualization for Keller-Segel (cell density and chemotactic field). */

event movie (i = 1; i += 10)
{
  output_ppm (C1, linear = true, spread = 2, file = "f.mp4", n = 200);
  fprintf (stderr, "%d %g %g %d %d\n", i, t, dt, mgd1.i, mgd2.i);
}

/**
### event final()

Save final steady-state or aggregation pattern as a PNG image.

The filename encodes the $\mu$ parameter value for identification. */

event final (t = 3000)
{
  char name[80];
  sprintf (name, "mu-%g.png", mu);
  output_ppm (C1, file = name, n = 200, linear = true, spread = 2);
}

/**
## Time Integration

### event integration()

Advance the system by one timestep (currently Brusselator implementation).

For chemotaxis, this would use operator splitting or coupled implicit solvers
to handle the chemotactic flux term $\chi \nabla \cdot (\rho \nabla c)$.

#### Current Algorithm (Brusselator)

1. Set adaptive timestep (max 1.0 for stability)
2. Solve $C_1$ with implicit diffusion:
   $$
   \partial_t C_1 = \nabla^2 C_1 + k k_a + k (C_1 C_2 - k_b - 1) C_1
   $$
3. Solve $C_2$ with implicit diffusion (coefficient $D$):
   $$
   \partial_t C_2 = D \nabla^2 C_2  + k k_b C_1 - k C_1^2 C_2
   $$

TODO: Implement Keller-Segel chemotactic coupling. */

event integration (i++)
{
  dt = dtnext (1.);

  /**
  Solve for $C_1$ with source term $r = k \cdot ka$ and coefficient
  $\beta = k(C_1 C_2 - k_b - 1)$. */

  scalar r[], beta[];
  
  foreach() {
    r[] = k*ka;
    beta[] = k*(C1[]*C2[] - kb - 1.);
  }
  mgd1 = diffusion (C1, dt, r = r, beta = beta);

  /**
  Solve for $C_2$ with diffusion coefficient $D$ and source/sink terms. */

  foreach() {
    r[] = k*kb*C1[];
    beta[] = - k*sq(C1[]);
  }
  const face vector c[] = {D, D};
  mgd2 = diffusion (C2, dt, c, r, beta);
}

/**
## Results

Currently showing Brusselator [Turing
patterns](http://en.wikipedia.org/wiki/The_Chemical_Basis_of_Morphogenesis).

TODO: Replace with Keller-Segel chemotaxis results showing cell aggregation
patterns, blow-up phenomena, or stable traveling waves depending on parameters.

<center>
<table>
<tr>
<td>![](brusselator/mu-0.04.png)</td>
<td>![](brusselator/mu-0.1.png)</td>
 <td>![](brusselator/mu-0.98.png)</td>
</tr>
<tr>
<td>$\mu=0.04$ (weak instability)</td> 
<td>$\mu=0.1$ (stripe patterns)</td> 
<td>$\mu=0.98$ (hexagonal patterns)</td>
</tr>
</table>
</center>

![Animation of pattern formation](brusselator/f.mp4)
*/
