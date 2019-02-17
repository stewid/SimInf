## SimInf, a framework for stochastic disease spread simulations
## Copyright (C) 2015  Pavol Bauer
## Copyright (C) 2015 - 2019  Stefan Engblom
## Copyright (C) 2015 - 2019  Stefan Widgren
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

library("SimInf")

## For debugging
sessionInfo()

## Define a tolerance
tol = 1e-8

## Expected phi
phi_exp <- c(1, 0.228812574827662, 0.0523896281322291, 0.0120297167737883,
0.00279666361731157, 0.000684437209448165, 0.000201227557942374,
9.06846888450544e-05, 6.53960261079013e-05, 5.96107912188877e-05,
5.82873150402961e-05, 5.79845461417288e-05, 5.79152823291854e-05,
6.39789936374658e-05, 9.42765209017355e-05, 0.000110545295549532,
0.000119281091914435, 0.000123971926951566, 0.000126490750928962,
0.000127843276476948, 0.000128569538166677, 0.000128959516783166,
0.000129168922460856, 0.00012928136641546, 0.000129341745117061,
0.000129374166491958, 0.000130685472705078, 0.000137406872101236,
0.000141301373209865, 0.000143557918014179, 0.000144865401132516,
0.000145622980641447, 0.000146061935993378, 0.000146316274711252,
0.00014646364319218, 0.000146549031171458, 0.000146598506520319,
0.000146627173437672, 0.000146643783571192, 0.000130521651183348,
7.64318047731465e-05, 6.35130989329827e-05, 6.04276222370606e-05,
5.96906934332883e-05, 5.95146868983672e-05, 5.94726498655136e-05,
5.94626098284665e-05, 5.94602118870071e-05, 5.9459639167682e-05,
5.94595023805956e-05, 5.94594697106539e-05, 5.94594619078333e-05,
5.94594600442239e-05, 5.82526953215486e-05, 5.79766262605224e-05,
5.79134705143992e-05, 5.78990225054351e-05, 5.78957172644106e-05,
5.78949611312813e-05, 5.78947881522746e-05, 5.78947485802262e-05,
5.78947395274109e-05, 5.78947374564171e-05, 5.78947369826401e-05,
5.78947368742551e-05, 5.78947368494601e-05, 9.10094753200943e-05,
0.000108790999633954, 0.000118339092764896, 0.000123466104421938,
0.000126219140883875, 0.000127697430825022, 0.000128491223853306,
0.000128917464575247, 0.000129146341809699, 0.000129269241349217,
0.000129335234355685, 0.000129370670427455, 0.000129389698431324,
0.000136656076909902, 0.000140866348814994, 0.000143305856964351,
0.000144719352373896, 0.000145538357333289, 0.000146012903708877,
0.00014628786451345, 0.000146447181806639, 0.000146539493145154,
0.00014659298001433, 0.000146623971279545, 0.000146641928182509,
0.000146652332738819, 8.02844230341808e-05, 6.44332502244339e-05,
6.06473892479069e-05, 5.97431821267293e-05, 5.95272231881794e-05,
5.94756440067648e-05, 5.94633249429036e-05, 5.94603826834443e-05)

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

sis_e <- run(sis_e, threads = 1)
sis_e_phi_obs <- trajectory(sis_e, "phi", as.is = TRUE)[1,]
stopifnot(all(abs(sis_e_phi_obs - phi_exp) < tol))
sis_e_phi_obs <- trajectory(sis_e, "phi")$phi
stopifnot(all(abs(sis_e_phi_obs - as.numeric(phi_exp)) < tol))

## Run with sparse V
punchcard(sis_e) <- data.frame(time = seq(from = 1, to = 700, by = 7),
                               node = 1,
                               phi = TRUE)
sis_e <- run(sis_e, threads = 1)
sis_e_phi_obs <- trajectory(sis_e, "phi", as.is = TRUE)[1,]
stopifnot(all(abs(sis_e_phi_obs - phi_exp) < tol))
sis_e_phi_obs <- trajectory(sis_e, "phi")$phi
stopifnot(all(abs(sis_e_phi_obs - as.numeric(phi_exp)) < tol))

if (SimInf:::have_openmp()) {
    sis_e <- run(sis_e, threads = 2)
    sis_e_phi_obs <- trajectory(sis_e, "phi", as.is = TRUE)[1,]
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

sis_e3 <- run(sis_e3, threads = 1)
sis_e3_phi_obs <- trajectory(sis_e3, "phi", as.is = TRUE)[1,]
stopifnot(all(abs(sis_e3_phi_obs - phi_exp) < tol))
sis_e3_phi_obs <- trajectory(sis_e3, "phi")$phi
stopifnot(all(abs(sis_e3_phi_obs - as.numeric(phi_exp)) < tol))

## Run with sparse V
punchcard(sis_e3) <- data.frame(time = seq(from = 1, to = 700, by = 7),
                                node = 1,
                                phi = TRUE)
sis_e3 <- run(sis_e3, threads = 1)
sis_e3_phi_obs <- trajectory(sis_e3, "phi", as.is = TRUE)[1,]
stopifnot(all(abs(sis_e3_phi_obs - phi_exp) < tol))
sis_e3_phi_obs <- trajectory(sis_e3, "phi")$phi
stopifnot(all(abs(sis_e3_phi_obs - as.numeric(phi_exp)) < tol))

if (SimInf:::have_openmp()) {
    sis_e3 <- run(sis_e3, threads = 2)
    sis_e3_phi_obs <- trajectory(sis_e3, "phi", as.is = TRUE)[1,]
    stopifnot(all(abs(sis_e3_phi_obs - phi_exp) < tol))
    sis_e3_phi_obs <- trajectory(sis_e3, "phi")$phi
    stopifnot(all(abs(sis_e3_phi_obs - as.numeric(phi_exp)) < tol))
}

## Check decay of phi for various configurations of the intervals.
## [1, 2)  1*0.97
## [2, 3)  1*0.97*0.97
## [3, 4)  1*0.97*0.97*0.95
## [4, 5)  1*0.97*0.97*0.95*0.95
## [5, 6)  1*0.97*0.97*0.95*0.95*0.93
## [6, 7)  1*0.97*0.97*0.95*0.95*0.93*0.93
## [7, 8)  1*0.97*0.97*0.95*0.95*0.93*0.93*0.91
## [8, 9)  1*0.97*0.97*0.95*0.95*0.93*0.93*0.91*0.91
## [9, 10) 1*0.97*0.97*0.95*0.95*0.93*0.93*0.91*0.91*0.97
phi <- c(1, 0.97, 0.9409, 0.893855, 0.84916225,
         0.7897208925, 0.734440430025, 0.66834079132275,
         0.608190120103702, 0.589944416500591)

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
