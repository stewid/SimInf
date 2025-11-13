# Scheduled events

NB: This vignette is work-in-progress and not yet complete.

## Overview

This vignette describes how births, deaths and movements can be
incorporated into a model as scheduled events at predefined time-points.
Events can, for example, be used to simulate disese spread among
multiple subpopulations (e.g., farms) when individuals can move between
the subpopulations and thus transfer infection, see Figure 1. In SimInf,
we use `node` to denote a subpopulation.

$\ $

![\*\*Figure 1.\*\* Illustration of movements between nodes. Each time
step depicts movements during one time unit, for example, a day. The
network has \*N=4\* nodes where node \*1\* is infected and nodes
\*2\*--\*4\* are non-infected. Arrows indicate movements of individuals
from a source node to a destination node and labels denote the size of
the shipment. Here, infection may spread from node \*1\* to node \*3\*
at \*t=2\* and then from node \*3\* to node \*2\* at
\*t=3\*.](img/temporal-network.svg)

**Figure 1.** Illustration of movements between nodes. Each time step
depicts movements during one time unit, for example, a day. The network
has *N=4* nodes where node *1* is infected and nodes *2*–*4* are
non-infected. Arrows indicate movements of individuals from a source
node to a destination node and labels denote the size of the shipment.
Here, infection may spread from node *1* to node *3* at *t=2* and then
from node *3* to node *2* at *t=3*.

## A first example

Let us define the **6** movement events in Figure 1 to include them in
an SIR model. Below is a `data.frame`, that contains the movements.
Interpret it as follows:

1.  In time step **1** we move **9** individuals from node **3** to node
    **2**
2.  In time step **1** we move **2** individuals from node **3** to node
    **4**
3.  In time step **2** we move **8** individuals from node **1** to node
    **3**
4.  In time step **2** we move **3** individuals from node **4** to node
    **3**
5.  In time step **3** we move **5** individuals from node **3** to node
    **2**
6.  In time step **3** we move **4** individuals from node **4** to node
    **2**

``` r
events <- data.frame(
  event      = rep("extTrans", 6),  ## Event "extTrans" is
                                    ##  a movement between nodes
  time       = c(1, 1, 2, 2, 3, 3), ## The time that the event happens
  node       = c(3, 3, 1, 4, 3, 4), ## In which node does the event occur
  dest       = c(4, 2, 3, 3, 2, 2), ## Which node is the destination node
  n          = c(9, 2, 8, 3, 5, 4), ## How many individuals are moved
  proportion = c(0, 0, 0, 0, 0, 0), ## This is not used when n > 0
  select     = c(4, 4, 4, 4, 4, 4), ## Use the 4th column in
                                    ## the model select matrix
  shift      = c(0, 0, 0, 0, 0, 0)) ## Not used in this example
```

and have a look at the `data.frame`

``` r
events
```

    ##      event time node dest n proportion select shift
    ## 1 extTrans    1    3    4 9          0      4     0
    ## 2 extTrans    1    3    2 2          0      4     0
    ## 3 extTrans    2    1    3 8          0      4     0
    ## 4 extTrans    2    4    3 3          0      4     0
    ## 5 extTrans    3    3    2 5          0      4     0
    ## 6 extTrans    3    4    2 4          0      4     0

Now, create an SIR model where we turn off the disease dynamics (beta=0,
gamma=0) to focus on the scheduled events. Let us start with different
number of individuals in each node.

``` r
library(SimInf)

model <- SIR(u0 = data.frame(S = c(10, 15, 20, 25),
                             I = c(5,  0,  0,  0),
                             R = c(0,  0,  0,  0)),
             tspan = 0:3,
             beta = 0,
             gamma = 0,
             events = events)
```

The compartments that an event operates on, is controlled by the select
value specified for each event together with the model select matrix
(E). Each row in E corresponds to one compartment in the model, and the
non-zero entries in a column indicate which compartments to sample
individuals from when processing an event. Which column to use in E for
an event is determined by the event select value. In this example, we
use the 4^(th) column which means that all compartments can be sampled
in each movement event (see below).

``` r
model@events@E
```

    ## 3 x 4 sparse Matrix of class "dgCMatrix"
    ##   1 2 3 4
    ## S 1 . . 1
    ## I . 1 . 1
    ## R . . 1 1

In another case you might be interested in only targeting the
susceptibles, which means for this model that we select the first
column. Now, let us run the model and generate data from it. For
reproducibility, we first call the
[`set.seed()`](https://rdrr.io/r/base/Random.html) function and specify
the number of threads to use since there is random sampling involved
when picking inviduals from the compartments.

``` r
set.seed(1)
set_num_threads(1)
result <- run(model)
```

And plot (Figure 2) the number of individuals in each node.

``` r
plot(result, range = FALSE)
```

![\*\*Figure 2.\*\* Number of susceptible, infected and recovered
individuals in each
node.](scheduled-events_files/figure-html/unnamed-chunk-6-1.png)

**Figure 2.** Number of susceptible, infected and recovered individuals
in each node.

$\ $

Or use the
[`trajectory()`](http://stewid.github.io/SimInf/reference/trajectory.md)
function to more easily inspect the outcome in each node in detail.

``` r
trajectory(result)
```

    ##    node time  S I R
    ## 1     1    0 10 5 0
    ## 2     2    0 15 0 0
    ## 3     3    0 20 0 0
    ## 4     4    0 25 0 0
    ## 5     1    1 10 5 0
    ## 6     2    1 17 0 0
    ## 7     3    1  9 0 0
    ## 8     4    1 34 0 0
    ## 9     1    2  6 1 0
    ## 10    2    2 17 0 0
    ## 11    3    2 16 4 0
    ## 12    4    2 31 0 0
    ## 13    1    3  6 1 0
    ## 14    2    3 25 1 0
    ## 15    3    3 12 3 0
    ## 16    4    3 27 0 0

## Varying probability of picking individuals

It is possible to assign different probabillities for the compartments
that an event sample individuals from. If the weights in the select
matrix $E$ are non-identical, individuals are sampled from a biased urn.
To illustrate this, let us create movement events between two nodes for
the built-in SIR model, where we start with 300 individuals ($S = 100$,
$I = 100$, $R = 100$) in the first node and then move them, one by one,
to the second node.

``` r
u0 <- data.frame(S = c(100, 0),
                 I = c(100, 0),
                 R = c(100, 0))
```

``` r
events <- data.frame(
  event = rep("extTrans", 300), ## Event "extTrans" is a movement between nodes
  time = 1:300,                 ## The time that the event happens
  node = 1,                     ## In which node does the event occur
  dest = 2,                     ## Which node is the destination node
  n = 1,                        ## How many individuals are moved
  proportion = 0,               ## This is not used when n > 0
  select = 4,                   ## Use the 4th column in the model select matrix
  shift = 0)                    ## Not used in this example
```

Now, create the model. Then run it, and plot the number of individuals
in the second node.

``` r
model <- SIR(u0 = u0,
             tspan = 1:300,
             events = events,
             beta = 0,
             gamma = 0)
```

``` r
plot(run(model), index = 2)
```

![\*\*Figure 3.\*\* The individuals have an equal probability of being
selected regardless of
compartment.](scheduled-events_files/figure-html/unnamed-chunk-11-1.png)

**Figure 3.** The individuals have an equal probability of being
selected regardless of compartment.

$\ $

The probability to sample an individual from each compartment is

$$p_{S} = \frac{w_{S}*S}{w_{S}*S + w_{I}*I + w_{R}*R}$$

$$p_{I} = \frac{w_{I}*I}{w_{S}*S + w_{I}*I + w_{R}*R}$$

$$p_{R} = \frac{w_{R}*R}{w_{S}*S + w_{I}*I + w_{R}*R}$$

Where $w_{S}$, $w_{I}$ and $w_{R}$ are the weights in E. These
probabilities are applied sequentially, that is the probability of
choosing the next item is proportional to the weights amongst the
remaining items. Let us now double the weight to sample individuals from
the $I$ compartment and then run the model again.

``` r
model@events@E[2, 4] <- 2
plot(run(model), index = 2)
```

![\*\*Figure 4.\*\* The individuals in the \$I\$ compartment are more
likely of being selected for a movement
event.](scheduled-events_files/figure-html/unnamed-chunk-12-1.png)

**Figure 4.** The individuals in the $I$ compartment are more likely of
being selected for a movement event.

$\ $

And a much larger weight to sample individuals from the $I$ compartment.

``` r
model@events@E[2, 4] <- 10
plot(run(model), index = 2)
```

![\*\*Figure 5.\*\* The individuals in the \$I\$ compartment are even
more likely of being selected for a movement event compared to the
previous
example.](scheduled-events_files/figure-html/unnamed-chunk-13-1.png)

**Figure 5.** The individuals in the $I$ compartment are even more
likely of being selected for a movement event compared to the previous
example.

$\ $

Increase the weight for the $R$ compartment and run the model again.

``` r
model@events@E[3, 4] <- 4
plot(run(model), index = 2)
```

![\*\*Figure 6.\*\* The individuals in the \$I\$ and \$R\$ compartments
are more likely of being selected for a movement event compared to
individuals in the \$S\$
compartment.](scheduled-events_files/figure-html/unnamed-chunk-14-1.png)

**Figure 6.** The individuals in the $I$ and $R$ compartments are more
likely of being selected for a movement event compared to individuals in
the $S$ compartment.
