# Display the distribution of individual events over time

Plot the number of individuals that are affected by the scheduled events
over time. The events are aggregated per time-point and event type, and
four panels are produced in a 2x2 grid, one for each event type: 'Exit',
'Enter', 'Internal transfer', and 'External transfer'. All panels share
the same y-axis scale to ease comparison between the event types, and an
event type with no events is displayed as a line at zero. If the
time-points in the events are dates, the x-axis is displayed as dates.

## Usage

``` r
# S4 method for class 'SimInf_individual_events'
plot(x, frame.plot = FALSE, ...)
```

## Arguments

- x:

  The `SimInf_individual_events` object to plot.

- frame.plot:

  Draw a frame around each plot. Default is FALSE.

- ...:

  Additional arguments affecting the plot, passed on to the plot
  function for each panel.
