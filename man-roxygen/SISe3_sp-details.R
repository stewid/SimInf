##' @details
##' The \code{SISe3_sp} model contains two compartments in three age
##' categories: \strong{S}usceptible (\eqn{S_1, S_2, S_3}) and
##' \strong{I}nfected (\eqn{I_1, I_2, I_3}). Additionally, it includes
##' a continuous \strong{environmental} compartment (\eqn{\varphi}) to
##' model shedding of a pathogen to the environment. Moreover, it
##' includes a spatial coupling of the environmental contamination
##' among proximal nodes to capture between-node spread unrelated to
##' moving infected individuals.
##'
##' The model is defined by six state transitions:
##'
##' \deqn{S_1 \stackrel{\upsilon_1 \varphi S_1}{\longrightarrow} I_1}{
##' S_1 -- upsilon_1 phi S_1 --> I_1}
##' \deqn{I_1 \stackrel{\gamma_1 I_1}{\longrightarrow} S_1}{ I_1 --
##' gamma_1 I_1 --> S_1}
##' \deqn{S_2 \stackrel{\upsilon_2 \varphi S_2}{\longrightarrow} I_2}{
##' S_2 -- upsilon_2 phi S_2 --> I_2}
##' \deqn{I_2 \stackrel{\gamma_2 I_2}{\longrightarrow} S_2}{ I_2 --
##' gamma_2 I_2 --> S_2}
##' \deqn{S_3 \stackrel{\upsilon_3 \varphi S_3}{\longrightarrow} I_3}{
##' S_3 -- upsilon_3 phi S_3 --> I_3}
##' \deqn{I_3 \stackrel{\gamma_3 I_3}{\longrightarrow} S_3}{ I_3 --
##' gamma_3 I_3 --> S_3}
##'
##' where the transition rate from susceptible to infected in age
##' category \eqn{k} is proportional to the environmental
##' contamination \eqn{\varphi} and the transmission rate
##' \eqn{\upsilon_k}. The recovery rate \eqn{\gamma_k} moves
##' individuals from infected back to susceptible.
##'
##' The environmental infectious pressure \eqn{\varphi(t)} in each
##' node evolves according to:
##'
##' \deqn{\frac{d \varphi_i(t)}{dt} = \frac{\alpha \left(I_{i,1}(t) +
##' I_{i,2}(t) + I_{i,3}(t)\right)}{N_i(t)} +
##' \sum_k{\frac{\varphi_k(t) N_k(t) - \varphi_i(t) N_i(t)}{N_i(t)}
##' \cdot \frac{D}{d_{ik}}} - \beta(t) \varphi_i(t)}{
##' dphi(t)/dt=
##' alpha (I_1+I_2+I_3)/N+
##' D*sum_k(phi_k*N_k-phi_i*N_i)/(d_ik*N_i)-beta*phi_i}
##'
##' where \eqn{\alpha} is the average shedding rate of the pathogen to
##' the environment per infected individual and \eqn{N = S_1 + S_2 +
##' S_3 + I_1 + I_2 + I_3} the size of the node. Next comes the
##' spatial coupling among proximal nodes, where \eqn{D} is the rate
##' of the local spread and \eqn{d_{ik}} the distance between holdings
##' \eqn{i} and \eqn{k}. The seasonal decay and removal of the
##' pathogen is captured by \eqn{\beta(t)}. The environmental
##' infectious pressure \eqn{\varphi(t)}{phi(t)} in each node is
##' evolved each time unit by the Euler forward method. The value of
##' \eqn{\varphi(t)}{phi(t)} is saved at the time-points specified in
##' \code{tspan}.
##'
##' \strong{Seasonal Decay (\eqn{\beta(t)}):}
##' The decay rate \eqn{\beta(t)} is piecewise constant, defined by four
##' intervals determined by the parameters \code{end_t1}, \code{end_t2},
##' \code{end_t3}, and \code{end_t4} (days of the year, where
##' \code{0 <= day < 365}). The year is divided into four intervals based
##' on the sorted order of these endpoints. The interval that wraps around
##' the year boundary (from the last endpoint to day 365, then from day 0
##' to the first endpoint) receives the same rate as the interval
##' preceding the first endpoint. Three orderings are supported:
##'
##' \strong{Case 1:} \code{end_t1 < end_t2 < end_t3 < end_t4}
##' \itemize{
##'   \item Interval 1: \code{[0, end_t1)} with rate \code{beta_t1}
##'   \item Interval 2: \code{[end_t1, end_t2)} with rate \code{beta_t2}
##'   \item Interval 3: \code{[end_t2, end_t3)} with rate \code{beta_t3}
##'   \item Interval 4: \code{[end_t3, end_t4)} with rate \code{beta_t4}
##'   \item Interval 1 (wrap-around): \code{[end_t4, 365)} with rate \code{beta_t1}
##' }
##'
##' \strong{Case 2:} \code{end_t3 < end_t4 < end_t1 < end_t2}
##' \itemize{
##'   \item Interval 3: \code{[0, end_t3)} with rate \code{beta_t3}
##'   \item Interval 4: \code{[end_t3, end_t4)} with rate \code{beta_t4}
##'   \item Interval 1: \code{[end_t4, end_t1)} with rate \code{beta_t1}
##'   \item Interval 2: \code{[end_t1, end_t2)} with rate \code{beta_t2}
##'   \item Interval 3 (wrap-around): \code{[end_t2, 365)} with rate \code{beta_t3}
##' }
##'
##' \strong{Case 3:} \code{end_t4 < end_t1 < end_t2 < end_t3}
##' \itemize{
##'   \item Interval 4: \code{[0, end_t4)} with rate \code{beta_t4}
##'   \item Interval 1: \code{[end_t4, end_t1)} with rate \code{beta_t1}
##'   \item Interval 2: \code{[end_t1, end_t2)} with rate \code{beta_t2}
##'   \item Interval 3: \code{[end_t2, end_t3)} with rate \code{beta_t3}
##'   \item Interval 4 (wrap-around): \code{[end_t3, 365)} with rate \code{beta_t4}
##' }
##'
##' These different orderings allow the model to handle seasonal patterns
##' where, for example, a winter peak crosses the year boundary.
