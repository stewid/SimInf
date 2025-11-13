# Definition of the SIS model

Class to handle the SIS
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model.md).

## Details

The SIS model contains two compartments; number of susceptible (S), and
number of infectious (I). Moreover, it has two state transitions, \$\$S
\stackrel{\beta S I / N}{\longrightarrow} I\$\$ \$\$I \stackrel{\gamma
I}{\longrightarrow} S\$\$ where \\\beta\\ is the transmission rate,
\\\gamma\\ is the recovery rate, and \\N=S+I\\.

## Examples

``` r
## Create an SIS model object.
model <- SIS(u0 = data.frame(S = 99, I = 1),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the SIS model and plot the result.
set.seed(22)
result <- run(model)
plot(result)
```
