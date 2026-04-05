# Getting started with \`mparse\`

## Overview

The `mparse` function is the core engine for defining custom stochastic
disease models in SimInf. Instead of writing complex C code manually,
`mparse` allows you to describe your model’s transitions using a simple,
human-readable string syntax in R. The function then parses this
description, generates optimized C code, compiles it, and returns a
`SimInf_model` object ready for simulation.

This approach offers the best of both worlds: the ease of defining
models in R and the computational speed of compiled C code. It is
particularly powerful for models with complex propensity functions,
multiple compartments, or node-specific parameters.

In this vignette, we will explore: - The basic syntax for defining
transitions. - How to define variables and population sizes. - How to
incorporate global and local data. - How to run and visualize the
resulting model.
