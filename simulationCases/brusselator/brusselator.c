/**
# Coupled Reaction-Diffusion Equations: Brusselator Model

The [Brusselator](http://en.wikipedia.org/wiki/Brusselator) is a
theoretical model for a type of autocatalytic reaction. The
Brusselator model was proposed by Ilya Prigogine and his collaborators
at the Free University of Brussels.

## Physical Setup

Two chemical compounds with concentrations $C_1$ and $C_2$ interact
according to the coupled reaction--diffusion equations:
$$
\partial_t C_1 = \nabla^2 C_1 + k(ka - (kb + 1)C_1 + C_1^2 C_2)
$$
$$
\partial_t C_2 = D \nabla^2 C_2  + k(kb C_1 - C_1^2 C_2)
$$

We use the same parameters as [Pena and Perez-Garcia, 2001](/src/references.bib#pena2001).

## Implementation

We use a Cartesian (multi)grid, the generic time loop, and the
time-implicit diffusion solver from Basilisk.

## Author

Vatsal Sanjay  
Email: vatsalsy@comphy-lab.org  
CoMPhy Lab  
Last updated: Jan 30, 2026
*/

#include "grid/multigrid.h"
#include "run.h"
#include "diffusion.h"

/**
## Variables

We define scalar fields for the chemical concentrations `C1` and `C2`. */

scalar C1[], C2[];

/**
## Parameters

Model parameters from [Pena and Perez-Garcia, 2001](/src/references.bib#pena2001):

- `k`: Reaction rate constant (default: 1.0)
- `ka`: Parameter controlling $C_1$ production (default: 4.5)
- `D`: Diffusion coefficient ratio for $C_2$ (default: 8.0)
- `mu`: Control parameter for bifurcation analysis
- `kb`: Derived parameter based on `mu`
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

Main simulation driver. Initializes the grid and solver parameters,
then runs simulations for multiple control parameter values.

We configure:
- Grid resolution: 128 × 128
- Domain size: 64 × 64
- Diffusion solver tolerance: 1e-4

Here $\mu$ is the control parameter. For $\mu > 0$ the system is
supercritical (Hopf bifurcation). We test several values of $\mu$ to
observe different pattern formation regimes. */

int main()
{
  init_grid (128);
  size (64);
  TOLERANCE = 1e-4;

  /**
  Run three cases covering different bifurcation regimes:
  - $\mu = 0.04$: Weak instability
  - $\mu = 0.1$: Stripe patterns
  - $\mu = 0.98$: Hexagonal patterns
  */

  mu = 0.04; run();
  mu = 0.1;  run();
  mu = 0.98; run();
}

/**
## Initial Conditions

### event init()

Initialize concentration fields near the unstable stationary solution.

The marginal stability is obtained for `kb = kbcrit`. We calculate:
- $\nu = \sqrt{1/D}$: characteristic wavenumber
- $k_b^{crit} = (1 + ka \cdot \nu)^2$: critical bifurcation parameter
- $kb = k_b^{crit}(1 + \mu)$: actual parameter based on control value $\mu$
*/

event init (i = 0)
{
  double nu = sqrt(1./D);
  double kbcrit = sq(1. + ka*nu);
  kb = kbcrit*(1. + mu);

  /**
  The (unstable) stationary solution is $C_1 = ka$ and $C_2 = kb/ka$. We
  perturb it with random noise in $[-0.01, 0.01]$ to trigger pattern formation. */

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
deviation, highlighting pattern structures. Progress information (iteration,
time, timestep, and solver iterations) is printed to stderr for monitoring. */

event movie (i = 1; i += 10)
{
  output_ppm (C1, linear = true, spread = 2, file = "f.mp4", n = 200);
  fprintf (stderr, "%d %g %g %d %d\n", i, t, dt, mgd1.i, mgd2.i);
}

/**
### event final()

Save final steady-state pattern as a PNG image.

The image filename encodes the $\mu$ value for easy identification of
different bifurcation regimes. */

event final (t = 3000)
{
  char name[80];
  sprintf (name, "mu-%g.png", mu);
  output_ppm (C1, file = name, n = 200, linear = true, spread = 2);
}

/**
## Time Integration

### event integration()

Advance the system by one timestep using operator splitting.

We separate each reaction-diffusion equation into:
$$
\partial_t C = D \nabla^2 C + r + \beta C
$$
where `r` and `beta` are the source term and linear coefficient respectively.

#### Algorithm

1. Set adaptive timestep (max 1.0 for stability of reactive terms)
2. Solve $C_1$ with implicit diffusion:
   $$
   \partial_t C_1 = \nabla^2 C_1 + k k_a + k (C_1 C_2 - k_b - 1) C_1
   $$
3. Solve $C_2$ with implicit diffusion (diffusion coefficient $D$):
   $$
   \partial_t C_2 = D \nabla^2 C_2  + k k_b C_1 - k C_1^2 C_2
   $$
*/

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
  Solve for $C_2$ with anisotropic diffusion coefficient $D$ and
  source/sink terms depending on current $C_1$ field. */

  foreach() {
    r[] = k*kb*C1[];
    beta[] = - k*sq(C1[]);
  }
  const face vector c[] = {D, D};
  mgd2 = diffusion (C2, dt, c, r, beta);
}
