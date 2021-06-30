## This file is part of SimInf, a framework for stochastic
## disease spread simulations.
##
## Copyright (C) 2015 Pavol Bauer
## Copyright (C) 2017 -- 2019 Robin Eriksson
## Copyright (C) 2015 -- 2019 Stefan Engblom
## Copyright (C) 2015 -- 2021 Stefan Widgren
##
## SimInf is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SimInf is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(SimInf)

## Specify the number of threads to use.
max_threads <- set_num_threads(1)

## For debugging
sessionInfo()

## Define a tolerance
tol <- 1e-8

## Expected phi
phi_exp <- c(0.810011, 0.185349185610407, 0.0424465987871056,
             0.00975507058676853, 0.00227629753002237,
             0.000565394139653014, 0.000173994321933323,
             8.44545979644941e-05, 6.39707811474e-05,
             5.9284740887299e-05, 5.82127251826398e-05,
             5.79674823748003e-05, 5.79113786866402e-05,
             6.95407791782812e-05, 9.7263016625088e-05,
             0.000112148945427821, 0.000120142199101708,
             0.000124434313160683, 0.0001267390371,
             0.000127976597976407, 0.00012864112742251,
             0.000128997957856597, 0.000129189564051683,
             0.000129292450270145, 0.000129347696782111,
             0.000129377362340141, 0.000131884062252197,
             0.000138101356693643, 0.000141703770219126,
             0.000143791074163116, 0.000145000496047577,
             0.000145701257093339, 0.000146107290793875,
             0.000146342554107908, 0.000146478869952767,
             0.000146557853833598, 0.000146603618531295,
             0.000146630135429847, 0.000146645499803353,
             0.000117375145714429, 7.32919208901144e-05,
             6.27631756303809e-05, 6.02485121232044e-05,
             5.964791514813e-05, 5.95044698221693e-05,
             5.94702096403935e-05, 5.94620270102002e-05,
             5.94600726879108e-05, 5.94596059216608e-05,
             5.94594944401854e-05, 5.94594678141829e-05,
             5.94594614548841e-05, 5.91621626358214e-05,
             5.81846832104543e-05, 5.79610672710232e-05,
             5.79099111166633e-05, 5.78982082294025e-05,
             5.78955309841726e-05, 5.78949185163378e-05,
             5.78947784033424e-05, 5.78947463499832e-05,
             5.78947390172028e-05, 5.78947373396978e-05,
             5.78947369559385e-05, 5.78947368681467e-05,
             6.3973684217256e-05, 9.42736699178863e-05,
             0.000110543764665068, 0.00011928026987988,
             0.000123971485546074, 0.000126490513908745,
             0.000127843149204895, 0.000128569469825775,
             0.000128959480086351, 0.000129168902755875,
             0.000129281355834533, 0.000129341739435452,
             0.000129374163441121, 0.000130685471048975,
             0.000137406871141659, 0.00014130137265387,
             0.000143557917692025, 0.000144865400945854,
             0.000145622980533292, 0.000146061935930711,
             0.000146316274674941, 0.000146463643171141,
             0.000146549031159267, 0.000146598506513255,
             0.000146627173433579, 0.000146643783568821,
             0.000130521651182138, 7.64318047728573e-05,
             6.35130989329136e-05, 6.04276222370442e-05,
             5.96906934332844e-05, 5.95146868983662e-05,
             5.94726498655133e-05, 5.94626098284664e-05,
             5.94602118870071e-05)

## Check phi from the SISe model
sis_e <- SISe(u0      = data.frame(S = 100, I = 0),
              tspan   = seq(from = 1, to = 700, by = 7),
              events  = NULL,
              phi     = 1,
              upsilon = 0,
              gamma   = 0,
              alpha   = 0,
              beta_t1 = 0.19,
              beta_t2 = 0.085,
              beta_t3 = 0.075,
              beta_t4 = 0.185,
              end_t1  = 91,
              end_t2  = 182,
              end_t3  = 273,
              end_t4  = 365,
              epsilon = 0.000011)

sis_e <- run(sis_e)
sis_e_phi_obs <- trajectory(sis_e, "phi", format = "matrix")[1, ]
stopifnot(all(abs(sis_e_phi_obs - phi_exp) < tol))
sis_e_phi_obs <- trajectory(sis_e, "phi")$phi
stopifnot(all(abs(sis_e_phi_obs - as.numeric(phi_exp)) < tol))

## Run with sparse V
punchcard(sis_e) <- data.frame(time = seq(from = 1, to = 700, by = 7),
                               node = 1,
                               phi = TRUE)
sis_e <- run(sis_e)
sis_e_phi_obs <- trajectory(sis_e, "phi", format = "matrix")[1, ]
stopifnot(all(abs(sis_e_phi_obs - phi_exp) < tol))
sis_e_phi_obs <- trajectory(sis_e, "phi")$phi
stopifnot(all(abs(sis_e_phi_obs - as.numeric(phi_exp)) < tol))

if (SimInf:::have_openmp() && max_threads > 1) {
    set_num_threads(2)
    sis_e <- run(sis_e)
    set_num_threads(1)
    sis_e_phi_obs <- trajectory(sis_e, "phi", format = "matrix")[1, ]
    stopifnot(all(abs(sis_e_phi_obs - phi_exp) < tol))
    sis_e_phi_obs <- trajectory(sis_e, "phi")$phi
    stopifnot(all(abs(sis_e_phi_obs - as.numeric(phi_exp)) < tol))
}

## Check phi from the SISe3 model
sis_e3 <- SISe3(u0      = data.frame(S_1 = 10, I_1 = 0,
                                     S_2 = 20, I_2 = 0,
                                     S_3 = 70, I_3 = 0),
                tspan   = seq(from = 1, to = 700, by = 7),
                events    = NULL,
                phi       = 1,
                upsilon_1 = 0,
                upsilon_2 = 0,
                upsilon_3 = 0,
                gamma_1   = 0,
                gamma_2   = 0,
                gamma_3   = 0,
                alpha     = 0,
                beta_t1   = 0.19,
                beta_t2   = 0.085,
                beta_t3   = 0.075,
                beta_t4   = 0.185,
                end_t1    = 91,
                end_t2    = 182,
                end_t3    = 273,
                end_t4    = 365,
                epsilon   = 0.000011)

sis_e3 <- run(sis_e3)
sis_e3_phi_obs <- trajectory(sis_e3, "phi", format = "matrix")[1, ]
stopifnot(all(abs(sis_e3_phi_obs - phi_exp) < tol))
sis_e3_phi_obs <- trajectory(sis_e3, "phi")$phi
stopifnot(all(abs(sis_e3_phi_obs - as.numeric(phi_exp)) < tol))

## Run with sparse V
punchcard(sis_e3) <- data.frame(time = seq(from = 1, to = 700, by = 7),
                                node = 1,
                                phi = TRUE)
sis_e3 <- run(sis_e3)
sis_e3_phi_obs <- trajectory(sis_e3, "phi", format = "matrix")[1, ]
stopifnot(all(abs(sis_e3_phi_obs - phi_exp) < tol))
sis_e3_phi_obs <- trajectory(sis_e3, "phi")$phi
stopifnot(all(abs(sis_e3_phi_obs - as.numeric(phi_exp)) < tol))

if (SimInf:::have_openmp() && max_threads > 1) {
    set_num_threads(2)
    sis_e3 <- run(sis_e3)
    set_num_threads(1)
    sis_e3_phi_obs <- trajectory(sis_e3, "phi", format = "matrix")[1, ]
    stopifnot(all(abs(sis_e3_phi_obs - phi_exp) < tol))
    sis_e3_phi_obs <- trajectory(sis_e3, "phi")$phi
    stopifnot(all(abs(sis_e3_phi_obs - as.numeric(phi_exp)) < tol))
}

## Check decay of phi for various configurations of the intervals.
## [1, 2)   1*0.97
## [2, 3)   1*0.97*0.97
## [3, 4)   1*0.97*0.97*0.95
## [4, 5)   1*0.97*0.97*0.95*0.95
## [5, 6)   1*0.97*0.97*0.95*0.95*0.93
## [6, 7)   1*0.97*0.97*0.95*0.95*0.93*0.93
## [7, 8)   1*0.97*0.97*0.95*0.95*0.93*0.93*0.91
## [8, 9)   1*0.97*0.97*0.95*0.95*0.93*0.93*0.91*0.91
## [9, 10)  1*0.97*0.97*0.95*0.95*0.93*0.93*0.91*0.91*0.97
## [10, 11) 1*0.97*0.97*0.95*0.95*0.93*0.93*0.91*0.91*0.97*0.97
phi <- c(0.97, 0.9409, 0.893855, 0.84916225, 0.7897208925,
         0.734440430025, 0.66834079132275, 0.608190120103702,
         0.589944416500591, 0.572246084005574)

model <- SISe(u0 = data.frame(S = 1, I = 0),
              tspan = 1:10,
              phi = 1,
              upsilon = 0,
              gamma = 0,
              alpha = 0,
              beta_t1 = 0.03,
              beta_t2 = 0.05,
              beta_t3 = 0.07,
              beta_t4 = 0.09,
              end_t1 = 3,
              end_t2 = 5,
              end_t3 = 7,
              end_t4 = 9,
              epsilon = 0)
stopifnot(all(abs(phi - trajectory(run(model), ~phi)$phi) < tol))

model <- SISe(u0 = data.frame(S = 1, I = 0),
              tspan = 1:10,
              phi = 1,
              upsilon = 0,
              gamma = 0,
              alpha = 0,
              beta_t1 = 0.05,
              beta_t2 = 0.07,
              beta_t3 = 0.09,
              beta_t4 = 0.03,
              end_t1 = 5,
              end_t2 = 7,
              end_t3 = 9,
              end_t4 = 3,
              epsilon = 0)
stopifnot(all(abs(phi - trajectory(run(model), ~phi)$phi) < tol))

model <- SISe(u0 = data.frame(S = 1, I = 0),
              tspan = 1:10,
              phi = 1,
              upsilon = 0,
              gamma = 0,
              alpha = 0,
              beta_t1 = 0.07,
              beta_t2 = 0.09,
              beta_t3 = 0.03,
              beta_t4 = 0.05,
              end_t1 = 7,
              end_t2 = 9,
              end_t3 = 3,
              end_t4 = 5,
              epsilon = 0)
stopifnot(all(abs(phi - trajectory(run(model), ~phi)$phi) < tol))
