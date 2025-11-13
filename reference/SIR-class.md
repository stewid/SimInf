# Definition of the SIR model

Class to handle the SIR
[`SimInf_model`](http://stewid.github.io/SimInf/reference/SimInf_model.md).

## Details

The SIR model contains three compartments; number of susceptible (S),
number of infectious (I), and number of recovered (R). Moreover, it has
two state transitions, \$\$S \stackrel{\beta S I / N}{\longrightarrow}
I\$\$ \$\$I \stackrel{\gamma I}{\longrightarrow} R\$\$ where \\\beta\\
is the transmission rate, \\\gamma\\ is the recovery rate, and
\\N=S+I+R\\.

## Examples

``` r
## Create an SIR model object.
model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
             tspan = 1:100,
             beta = 0.16,
             gamma = 0.077)

## Run the SIR model and plot the result.
set.seed(22)
result <- run(model)
plot(result)
```
