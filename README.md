# Keller-Segel chemotaxis

Basilisk reaction-diffusion simulations for Keller-Segel chemotaxis and related models.

## Requirements
- Basilisk source tree (installed locally)
- qcc available in PATH
- Build tools (Xcode CLI on macOS or build-essential on Linux)

## Quick start
1. Install Basilisk and generate `.project_config`:
   - macOS: `./reset_install_requirements.sh`
   - Linux/no-darcs: `./reset_install_requirements-no-darcs.sh`
2. Run a case:
   - `./simulationCases/runCases.sh keller-segel`
   - `./simulationCases/runCases.sh brusselator`
3. Clean outputs:
   - `./simulationCases/cleanup.sh keller-segel`

Outputs are written to `simulationCases/<case>/`.

## Structure
- `simulationCases/` case entry points and run scripts
- `simulationCases/<case>/` output directories (including legacy outputs)
- `src-local/` project-specific headers and helpers
- `postProcess/` analysis and plotting utilities
- `basilisk/` local Basilisk checkout (not tracked)
