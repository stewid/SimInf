##' @param replicates Number of model replicates to simulate (default
##'     \code{NULL}, treated as \code{1L}). When \code{replicates >
##'     1L}, each replicate is simulated independently using its own
##'     initial state (from \code{u0}), but shares the same parameters
##'     (\code{gdata}, \code{ldata}), scheduled events, and structure
##'     (transitions, compartments).
##'
##'     The \code{u0} argument must contain initial states for all
##'     replicates. Each replicate requires \emph{n} nodes, where
##'     \emph{n} is the number of nodes in the model. Thus, \code{u0}
##'     must have \code{replicates * n} nodes total. Nodes are grouped
##'     by replicate: the first \emph{n} nodes belong to replicate 1,
##'     the next \emph{n} to replicate 2, and so on. This allows
##'     different starting conditions per replicate if desired.
##'
##'     \code{ldata} remains unchanged: its data per node is shared
##'     across all replicates. Scheduled events are also shared—the
##'     same event schedule applies to each replicate.
##'
##'     Use this when you need multiple independent stochastic
##'     trajectories from the same model in a single simulation run.
##'     For identical starting conditions across replicates, simply
##'     repeat the same node pattern in \code{u0}.
