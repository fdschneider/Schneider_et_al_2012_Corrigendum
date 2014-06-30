# =======================
#
# Analyses and simulations performed for Corrigendum to Schneider, Scheu & Brose 2012  *Ecology Letters*, 15:436-443, Mai 2012. doi: 10.1111/j.1461-0248.2012.01750.x 
#
# The MIT License (MIT)
# 
# Copyright (c) 2014 the authors
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# =====================

rm(list=ls())


#install.packages('deSolve')    #run if deSolve-package is not installed
library(deSolve)

#setwd("")  # set your working directory
raw_data <- read.delim("data\\data_schneider.txt", dec = ".", header = TRUE) 

# removing replicate 43, because the springtail and predator populations were profoundly affected by very low water content.
# removing replicates of controls and single predator treatments for clarity
exp_data <- raw_data[raw_data$ID != 43 & raw_data$treat_name %in% c("full","ko_lith","ko_pard","ko_hypo"),] 



#	calculating experimental interaction strengths
# 	----------------------------------------------

# function for calculation of interaction strengths with 
#    x = biomass or vector of biomasses of target species in absence of the influencing species, 
#    y = biomass or vector of biomasses of target species in presence of the influencing species, and
#    n = individual density of influencing species (for per capita effects)

IS <- function(x, y, n = NULL) {
lx <- length(x)   # length of vector x
ly <- length(y)   # length of vector y

# repeating vectors for systematic combination of all values in x with all values in y
X <- rep(x, each = ly)  
Y <- rep(y, times = lx)
N <- rep(n, each = ly)

# calculating interaction strengths and per capita Interaction strengths
IS.result <- log10((X+1)/(Y+1))
pc.IS.result <- IS.result/N

if(lx > 1 | ly > 1) {
# creating return object with mean and standard deviation of population level and per capita interaction strengths, as well as p of a student's t-test 
out <- c(mean.IS = mean(IS.result), 
	sd.IS = sd(IS.result), 
	p = round(t.test(log10(x+1),log10(y+1))$p.value, digits = 3), 
	mean.pc.IS = mean(pc.IS.result),
	sd.pc.IS = sd(pc.IS.result) )
} else {
# in case of only one value in x and y, simplified output
out <- c(IS = IS.result, pc.IS = pc.IS.result )
}
return(out)
}


# applying the function to all replicate combinations of interest
#"lith>pard","lith>hypo","lith>het","pard>lith","pard>hypo","pard>het","hypo>lith", "hypo>pard", "hypo>het"

exp_IS <- with(exp_data, data.frame(rbind( 
	IS(B1_pard[treat_name == "full"], B1_pard[treat_name == "ko_lith"], (N0_lith[treat_name == "full"] + N1_lith[treat_name == "full"])/2 ),
	IS(B1_hypo[treat_name == "full"], B1_hypo[treat_name == "ko_lith"], (N0_lith[treat_name == "full"] + N1_lith[treat_name == "full"])/2  ),
	IS(B1_het[treat_name == "full"], B1_het[treat_name == "ko_lith"], (N0_lith[treat_name == "full"] + N1_lith[treat_name == "full"])/2  ),
	IS(B1_lith[treat_name == "full"], B1_lith[treat_name == "ko_pard"], (N0_pard[treat_name == "full"] + N1_pard[treat_name == "full"])/2),
	IS(B1_hypo[treat_name == "full"], B1_hypo[treat_name == "ko_pard"], (N0_pard[treat_name == "full"] + N1_pard[treat_name == "full"])/2),
	IS(B1_het[treat_name == "full"], B1_het[treat_name == "ko_pard"], (N0_pard[treat_name == "full"] + N1_pard[treat_name == "full"])/2),
	IS(B1_lith[treat_name == "full"], B1_lith[treat_name == "ko_hypo"], (N0_hypo[treat_name == "full"] + N1_hypo[treat_name == "full"])/2),
	IS(B1_pard[treat_name == "full"], B1_pard[treat_name == "ko_hypo"], (N0_hypo[treat_name == "full"] + N1_hypo[treat_name == "full"])/2),
	IS(B1_het[treat_name == "full"], B1_het[treat_name == "ko_hypo"], (N0_hypo[treat_name == "full"] + N1_hypo[treat_name == "full"])/2) 
	)))

	
	

#	functional response equations
# 	-----------------------------

# attack rate (Eq. 2 in Schneider et al. 2012)  < encounter rate >     <                 success probability                > 
a <- function(mj, mi, parms.fr) with(parms.fr,  a0 * mi^ai * mj^aj  *  ( (mj / (mi * Ropt)) * exp( 1 - (mj/(mi*Ropt)) ) )^gam  )

# handling time (Eq. 3 in Schneider et al. 2012)
h <- function(mj, mi, parms.fr) with(parms.fr,  h0 * mj^hj * mi^hi)

# predator interference (Eq. 4, in Schneider et al. 2012)
cj <- function(mj, parms.fr) with(parms.fr, c0 * mj^cj)

# metabolic demands (Eq. 8 in Schneider et al. 2012 , after Ehnes et al 2011 Ecology Letters)
x <- function(mj, parms.fr) with(parms.fr, exp(i0) * mj^ax *  exp( ea / (boltz * (273.15+Temp)) ) * 3 / 7 / mj )

# per capita assimilation efficiency
e <- function(mj, mi, parms.fr)  with(parms.fr, e0 * mi/mj )



#	Parameter definitions
# ---------------------

# setting timesteps for simulation
t0 <- as.POSIXct(0, origin = "2008-06-10")		#June 10th 2008, springtails initialised
t1 <- as.POSIXct(0, origin = "2008-07-01")		#July 07th 2008, predators initialised
t2 <- as.POSIXct(0, origin = "2008-08-18")		#August 08th 2008, experiment terminated

exptime <- as.numeric(t2-t1)*24       		# 48 days of experiment = 1152 hours
times <- seq(0, exptime, length = 100)					#  population densities are calculated for each step in 'times'	


# estimation of growth parameters
# raw data:
growth <- data.frame(
			t = rep( as.numeric( c(0,		t1-t0,					t2-t0) ), each = 5)*24 , 
			Bb = c(100,100,100,100,100, 	1805,668,970,513,603, 	12090,5491,6583,5947,6892))

growthfit <- nls( Bb ~ K / ( 1 + exp( a - r*t ) ) ,  # logistic growth function fitted by nonlinear least squares
  nls.control(maxiter = 1000),
  data = growth,
  start = list(r = 0.001, K = 10000, a = 1),
  trace = FALSE )
  
plot(growth)
lines(seq(0,2000, length = 100), predict(growthfit, newdata = list(t = seq(0,2000, length = 100))), col = "red" )

######## parameters for simulation of population dynamics ########

parms <- data.frame(
r = coef(growthfit)[1], K = coef(growthfit)[2], 	# intrinsic growth parameters of springtails (in Eq. 5, Schneider et al. 2012)
ai = 0.25, aj = 0.25, a0 = .15, 					# attack rate parameters, encounter rate (Eq. 2, Schneider et al. 2012 )
Ropt = 200, gam = 1, 								# attack rate parameters, success rate (Eq. 2, Schneider et al. 2012 )
h0= 8, hj = -0.25, hi = 0.25, 						# handling time (Eq. 3, Schneider et al. 2012 )
c0 = 1, cj = 0.5, 									# predator interference (Eq. 4, Schneider et al. 2012 )
e0 = 0.85,											# conversion efficiency (Eq. 7, Schneider et al. 2012 )
boltz = 0.00008617343, Temp = 16.2, i0 = 23.055335, ax = 0.695071 , ea = -0.68642) # metabolic rate (Eq. 8, Schneider et al. 2012 )


#times <- seq(0, maxtime, length = tsteps/20)

N <- c(N1 = 912, N2 = 350, N3 = 4, N4 = 2)			# initial densities of springtails, mites, spiders, centipedes
m <- c(0.104, 0.157, 25.69, 125.3)  # mean body mass of springtails, mites, spiders, centipedes mg


#  defining differential equations for population dynamics
#  ---------------------------------------------------------

# extinction threshold
extinct = 1 

# dynamic equations (function handed to ode() )
dynamics <- function(t, n, parms) {

n[which( n < extinct)] <- 0

	# Heteromurus nitidus (Springtails)
	dN1 <-  with(as.list(c(n, parms)), { N1*r*(1-N1/K) - 	(N2*a12*N1)/(1+c2*(N2-1)+(a12*h12*N1)) - (N3*a13*N1)/(1+c3*(N3-1)+(a13*h13*N1)+(a23*h23*N2)) - (N4*a14*N1)/(1+c4*(N4-1)+(a14*h14*N1)+(a24*h24*N2)+(a34*h34*N3)) })

	# Hypoaspis miles (Mites)
	dN2 <- with(as.list(c(n, parms)), {(N2*a12*N1*e12)/(1+c2*(N2-1)+(a12*h12*N1)) - (N3*a23*N2)/(1+c3*(N3-1)+(a13*h13*N1)+(a23*h23*N2)) - (N4*a24*N2)/(1+c4*(N4-1)+(a14*h14*N1)+(a24*h24*N2)+(a34*h34*N3)) - x2*N2})

	# Pardosa lugubris (Spiders)
	dN3 <- with(as.list(c(n, parms)), { (N3*a13*N1*e13)/(1+c3*(N3-1)+(a13*h13*N1)+(a23*h23*N2)) + (N3*a23*N2*e23)/(1+c3*(N3-1)+(a13*h13*N1)+(a23*h23*N2)) - (N4*a34*N3)/(1+c4*(N4-1)+(a14*h14*N1)+(a24*h24*N2)+(a34*h34*N3)) - x3*N3})

	# Lithobius forficatus (Centipedes)
	dN4 <- with(as.list(c(n, parms)), {(N4*a14*N1*e14)/(1+c4*(N4-1)+(a14*h14*N1)+(a24*h24*N2)+(a34*h34*N3)) + (N4*a24*N2*e24)/(1+c4*(N4-1)+(a14*h14*N1)+(a24*h24*N2)+(a34*h34*N3)) + (N4*a34*N3*e34)/(1+c4*(N4-1)+(a14*h14*N1)+(a24*h24*N2)+(a34*h34*N3)) - x4*N4})

dn = c(dN1 = dN1, dN2 = dN2, dN3 = dN3, dN4 = dN4)

return(list(dn))

}

# body-mass dependent parameters for all predator-prey pairs
# functions take following arguments:  "predator mass", "prey mass", "parameter set"
a12 <- a(m[2], m[1], parms)
a13 <- a(m[3], m[1], parms)
a14 <- a(m[4], m[1], parms)
a23 <- a(m[3], m[2], parms)
a24 <- a(m[4], m[2], parms)
a34 <- a(m[4], m[3], parms)

h12 <- h(m[2], m[1], parms)
h13 <- h(m[3], m[1], parms)
h14 <- h(m[4], m[1], parms)
h23 <- h(m[3], m[2], parms)
h24 <- h(m[4], m[2], parms)
h34 <- h(m[4], m[3], parms)

c2 <- cj(m[2], parms)
c3 <- cj(m[3], parms)
c4 <- cj(m[4], parms)

x2 <- x(m[2], parms)
x3 <- x(m[3], parms)
x4 <- x(m[4], parms)

e12 <- e(m[2], m[1], parms)
e13 <- e(m[3], m[1], parms)
e14 <- e(m[4], m[1], parms)
e23 <- e(m[3], m[2], parms)
e24 <- e(m[4], m[2], parms)
e34 <- e(m[4], m[3], parms)


#  running population dynamics ( using ode() from the deSolve package)
#  ---------------------------------------------------------

# full module with all three predators
n = N 
full <- ode(n, times, dynamics, parms, method = rkMethod("rk45dp7") )
full <- as.data.frame(full)
full[full < extinct] <- 0

# knockout module without Lithobius
n = N * c(1,1,1,0)  # here, initial population density of lithobius is set 0
ko_lith <- ode(n, times, dynamics, parms, method = rkMethod("rk45dp7") )
ko_lith <- as.data.frame(ko_lith)
ko_lith[ko_lith < extinct] <- 0

# knockout module without Pardosa
n = N * c(1,1,0,1) 
ko_pard <- ode(n, times, dynamics, parms, method = rkMethod("rk45dp7") )
ko_pard <- as.data.frame(ko_pard)
ko_pard[ko_pard < extinct] <- 0

# knockout module without Hypoaspis
n = N * c(1,0,1,1) 
ko_hypo <- ode(n, times, dynamics, parms, method = rkMethod("rk45dp7") )
ko_hypo <- as.data.frame(ko_hypo)
ko_hypo[ko_hypo < extinct] <- 0



#  calculating simulated Interaction strengths
#  -----------------------------

# collecting final population densities in one table
sim_out <-  data.frame(
					full = as.numeric(full[times == exptime, -1 ] * m/1000),
					ko_lith = as.numeric(ko_lith[times == exptime, -1 ]* m/1000),
					ko_pard = as.numeric(ko_pard[times == exptime, -1 ]* m/1000),
					ko_hypo = as.numeric(ko_hypo[times == exptime, -1 ]* m/1000),
					row.names = c("het","hypo","pard","lith")
	)

# interaction strength for all species pairs (returns population and per capita IS)
sim_IS <- data.frame(rbind(
IS(sim_out[3,1], sim_out[3,2], (full[times == exptime,5]+N[4]/2)),
IS(sim_out[2,1], sim_out[2,2], (full[times == exptime,5]+N[4]/2)),
IS(sim_out[1,1], sim_out[1,2], (full[times == exptime,5]+N[4]/2)), 
IS(sim_out[4,1], sim_out[4,3], (full[times == exptime,4]+N[3]/2)), 
IS(sim_out[2,1], sim_out[2,3], (full[times == exptime,4]+N[3]/2)),
IS(sim_out[1,1], sim_out[1,3], (full[times == exptime,4]+N[3]/2)), 
IS(sim_out[4,1], sim_out[4,4], (full[times == exptime,3]+N[2]/2)),
IS(sim_out[3,1], sim_out[3,4], (full[times == exptime,3]+N[2]/2)),
IS(sim_out[1,1], sim_out[1,4], (full[times == exptime,3]+N[2]/2)) 
 ))

sim_IS

#  visualisations of time series
#  -----------------------------

layout(matrix(c(1,2,5, 3,4,0), nrow = 2, byrow = T), width = c(1,1,0.5))

## plotting function for timeseries
timeseries <- function(output, ymax = "auto", ...) {
	palts <- c("#EE801D","#C5CF56", "#A6A11B", "#6E7320", "#45501E")
	if(ymax == "auto") ymax <- 1.2*max(output[,-1]) else ymax <- ymax
	
	S <- length(output[1,-1])
	plot(output[,2]~output$time, col = palts[1], type = "l", bty = "n", xaxs = "i", ylim = c(.1, ymax), ylab = "Pop. density", xlab = "time [h]", ...)
		for(i in 2:S) lines(output[,i+1]~output$time, col = palts[i])
}


timeseries(ko_lith, main = "mites & spiders", log = "y", ymax = 7000)

timeseries(ko_pard, main = "mites & centipedes",log = "y", ymax = 7000)

timeseries(ko_hypo, main = "spiders & centipedes", log = "y", ymax = 7000)

timeseries(full, main = "full module", log = "y", ymax = 7000)

# predicted vs experimental interaction strengths
pal <- c("#EE801D","#C5CF56", "#A6A11B", "#6E7320", "#45501E")

plot(NA, NA, xlim = c(-0.25,0.25), ylim = c(-0.25, 0.25), type = "n", xaxs = "i", yaxs = "i", tck = 0.03, las = 1, cex.axis= 0.75, cex.lab = 0.7)
abline(h = 0, col = "grey75"); abline(v = 0, col = "grey75")
abline(a = 0, b = 1,  col = "grey75", lty = 2)
points(exp_IS$mean.IS ~ sim_IS$IS, pch = 21, cex = 1.5, col = "white", lwd = 1.5, bg = pal[c(4,4,4,3,3,3,2,2,2)])
#text(0.2, 0.23, labels = round(rho, digits = 4))
box()


# file output of predicted vs experimental interaction strengths 
#  -------------------------------------------------------------

pdf("manuscript/schneider_fig4b.pdf",  width = 4, height = 4, paper = "special")
#png("manuscript/schneider_fig4b.png",  width = 900, height = 900, res = 230)

par(mfrow = c(1,1), mar = c(4.5,4.5,3.5,3.5)+0.1, bty = "o")

pal <- c("#EE801D", "#C5CF56", "#8A911E", "#45501E")

plot(IS_table$mean.IS ~ IS_table$IS, 
     xlim = c(-0.2,0.2),  xaxs = "i", xlab = expression(paste("simulated ", italic(IS[ij]))),
     ylim = c(-0.2, 0.2), yaxs = "i", ylab = expression(paste("experimental ", italic(IS[ij]))),
     type = "n", tck = 0.03, las = 1, cex.axis= 0.75, cex.lab = 1)
abline(h = 0, col = "grey75"); abline(v = 0, col = "grey75")
abline(a = 0, b = 1,  col = "grey75", lty = 2)
points(IS_table$mean.IS ~ IS_table$IS, pch = 21, cex = 1.5, col = "white", lwd = 1.5, bg = pal[c(4,4,4,3,3,3,2,2,2)])
#text(-.2,0.18, label = paste("corCoef: ", round(cor(normal(new$IS_table$mean.IS), normal(new$IS_table$IS ), method = "pearson"), digits = 4), sep = ""), pos = 4)

#abline(lm(new$IS_table$mean.IS ~ -1 + new$IS_table$IS))

dev.off()


# correlation coefficients
# ------------------------

 r = cor.test(exp_IS$mean.IS, sim_IS$IS, alternative = c("greater"), method = "pearson")
 pc.r = cor.test(exp_IS$mean.pc.IS, sim_IS$pc.IS.N4, alternative = c("greater"), method = "pearson") 

 rho = cor.test(exp_IS$mean.IS, sim_IS$IS, alternative = c("greater"), method = "spearman")$estimate 
text(0.2, 0.2, labels = round(r$estimate, digits = 3))


# result table
# ------------

IS.name = c("lith>pard","lith>hypo","lith>het","pard>lith","pard>hypo","pard>het","hypo>lith", "hypo>pard", "hypo>het")
(IS_table = cbind(IS.name, round(exp_IS, digits = 4), round(sim_IS, digits = 4) ))

 
 

