# Extract the events from a `SimInf_model` object

Extract the scheduled events from a `SimInf_model` object.

## Usage

``` r
events(object, ...)

# S4 method for class 'SimInf_model'
events(object, ...)
```

## Arguments

- object:

  The `model` to extract the events from.

- ...:

  Additional arguments affecting the generated events.

## Value

[`SimInf_events`](http://stewid.github.io/SimInf/reference/SimInf_events-class.md)
object.

## Examples

``` r
## Create an SIR model that includes scheduled events.
model <- SIR(u0     = u0_SIR(),
             tspan  = 1:(4 * 365),
             events = events_SIR(),
             beta   = 0.16,
             gamma  = 0.077)

## Extract the scheduled events from the model and display summary
summary(events(model))
#> Number of scheduled events: 466692
#>  - Exit: 182535 (n: min = 1 max = 1 avg = 1.0)
#>  - Enter: 182685 (n: min = 1 max = 1 avg = 1.0)
#>  - Internal transfer: 0
#>  - External transfer: 101472 (n: min = 1 max = 1 avg = 1.0)

## Extract the scheduled events from the model and plot them
plot(events(model))
```
