# GillespieClockModel.jl
Stochastic simulation of cyanobacterial circadian clock dynamics in random light environments

The code implements a mass-action adaptation of the model in Chew et al. Nature Communications 9: 3004 (2018). Light inputs are assumed to modulate the phosphorylation rate. A mass action-reaction model with piece-wise constant light inputs is read from a file and simulated using Gillespie-type stochastic simulation.

## Light inputs

Several input functions are provided:

- constant
- square wave
- sine wave
- noisy day start
- noisy day end
- Caribbean 1
- Caribbean 2

## System requirements

Tested using Julia 1.10+. Needs Catalyst.jl, ModelingToolkit.jl, JumpProcesses.jl, Plots.jl, CSV.jl, DataFrames.jl

## Installation

No installation required.

## Usage

Run `julia randenv.jl`

## Demo

Expected output shown in `output.pdf`
