#
#  ********************************
#  **  Written by Yoan Diekmann  **
#  **    y.diekmann@ucl.ac.uk    **
#  ********************************
#
#
#
library(LaplacesDemon)
#
add.alpha <- function(alpha, col=1) {
    if(missing(alpha)) stop("Please provide a vector of alpha values.")
    col_t = col2rgb(col)/255
    apply( as.matrix(alpha), 2, function(x) rgb(col_t[1], col_t[2], col_t[3], alpha=x) )
}
#
#
#
# *****************************
# *         FUNCTIONS         *
# *****************************
#
E_of_T_given_D <- function(x, N_e, G, n, mu) {
    #
    lambda_is = vector(mode='numeric', length=(n-1))
    for (i in 1:(n-1)) {
        lambda_is[i] = ((i+1)*(i+(2*N_e*mu))) / (2*N_e*G)
    }
    #
    total = 0
    for (i in 1:(n-1)) {
        f_i = lambda_is[i] * exp(-lambda_is[i]*x)
        tmp = rep(1, (n-1))
        #
        for (j in 1:(n-1)) {
            if (i==j) { next }
            tmp[j] = (lambda_is[j] / (lambda_is[j] - lambda_is[i]))
        }
        total = total + f_i * prod(tmp)
    }
    return(total)
}
#
calculate_area <- function(xs, ys) {
    stopifnot(length(xs)==length(ys))
    #
    area = vector(mode='numeric', length=(length(xs)-1))
    area[1] = (xs[2]-xs[1]) * ((ys[1]+ys[2])/2)
    for (i in 2:(length(xs)-1)) {
        area[i] = area[i-1] + (xs[i+1]-xs[i]) * ((ys[i]+ys[i+1])/2)
    }
    return(area)
}
#
#
#
# *****************************
# *       Demographics        *
# *****************************
#
# https://en.wikipedia.org/wiki/Medieval_demography
# https://en.wikipedia.org/wiki/Demographics_of_Italy
years_italy = c(1000, 1100, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1861, 1871, 1881, 1901, 1911, 1921, 1931, 1936,
          1951, 1961, 1971, 1981, 1991, 2001, 2011, 2015)
italy = c(7000000, 7500000, 8000000, 9000000, 10000000, 8000000, 10000000, 10500000, 11000000, 22182377, 27303509,
          28953480, 32965504, 35845048, 39943528, 41651000, 42943602, 47515537, 50623569, 54136547, 56556911, 56778031,
          56995744, 59433744, 60656000)
#
#
#
# LS-fitting of base
bases = seq(from=4, to=18, length.out=1000)
rs = exp(seq(-3, -8, length.out=1001))
#
pf = matrix(nrow=1001, ncol=1000)
for (i in 1:1001) {
    for (j in 1:1000) {
        pf[i, j] = (sum(((italy[length(italy)]*
                    bases[j]^(bases[j]^(rs[i]*(years_italy-years_italy[length(years_italy)]))-1)) - italy)^2))^(1/2) } }
# ROW == rs (but on x-axis in image!)
r_which = which(pf==min(pf, na.rm=TRUE), arr.ind = TRUE)[1]
r = rs[r_which]
# COL == bases (but on y-axis in image!)
b_which = which(pf==min(pf, na.rm=TRUE), arr.ind = TRUE)[2]
b = bases[b_which]
#
#pdf('population_size_param_fitting.pdf') #, width=5, height=4
png('population_size_param_fitting.png', width = 480, height = 480)
pal <- colorRampPalette(rev(c('white', 'grey', 'blue', 'cyan', 'yellow', 'red'))) #, bias=2 
image(pf, xaxt='n', yaxt='n', xlab='r', ylab='b', col=pal(100))
axis(1, at=seq(0, 1, length.out = 11), labels = round(exp(seq(-3, -8, length.out=11)), 3), cex=0.8)
axis(2, at=seq(0, 1, length.out = 11), labels = seq(from=4, to=18, length.out=11), cex=0.8) # , las=2
contour(pf, levels=c(min(pf, na.rm=TRUE)*1.05, min(pf, na.rm=TRUE)*1.01), labels=c('5%', '1%'), add=TRUE)
lines(c(r_which/1001, r_which/1001), par('usr')[3:4])
lines(par('usr')[1:2], c(b_which/1000, b_which/1000))
dev.off()
#
#
#
pdf('population_size.pdf', width=7, height=6)
plot(years_italy, italy, type='l', ylim=c(0, 63000000), xlim=c(500, 2014),
     xlab='year AD', ylab='census population size Italy', yaxs='i')
#
xs = years_italy[1]:years_italy[length(years_italy)]
ys = sapply((years_italy[1]-years_italy[length(years_italy)]):0, function(x) italy[length(italy)]*b^(b^(r*x)-1) )
lines(xs, ys, col='blue', lty='solid')
#
xs = 500:years_italy[1]
ys = sapply((years_italy[1]-years_italy[length(years_italy)]-500):(years_italy[1]-years_italy[length(years_italy)]),
                                                                        function(x) italy[length(italy)]*b^(b^(r*x)-1) )
lines(xs, ys, col='green', lty='solid')
#
legend('bottomright', c('demographic data', 'hyper-exponential fit', 'extrapolation'),
       lty=c('solid', 'solid', 'solid'), lwd=c(1, 1, 1), col=c('black', 'blue', 'green'), bty = 'n', cex=0.9)
dev.off()
#
#
#
# from https://en.wikipedia.org/wiki/Sicily#Demographics
# http://www.tacitus.nu/historical-atlas/population/italy.htm
years_sicily = c(1770, 1800, 1816, 1825, 1838, 1848, 1858, 1861, 1871, 1881, 1901, 1911, 1921, 1931, 1936, 1951, 1961,
                1971, 1981, 1991, 2001, 2011, 2014)
sicily = c(1300000, 1700000, 1600000, 1700000, 1900000, 2100000, 2300000, 2409000, 2590000, 2933000, 3568000, 3812000,
           4223000, 3906000, 4000000, 4487000, 4721000, 4681000, 4907000, 4966000, 4969000, 5002904, 5092080)
conv_factor_s = sicily[length(sicily)] / cal_mes[length(cal_mes)]
#
#
#
# http://www.tuttitalia.it/calabria/provincia-di-reggio-calabria/statistiche/censimenti-popolazione/
# http://www.tuttitalia.it/calabria/provincia-di-reggio-calabria/statistiche/popolazione-andamento-demografico/
# http://www.tuttitalia.it/sicilia/provincia-di-messina/statistiche/censimenti-popolazione/
# http://www.tuttitalia.it/sicilia/provincia-di-messina/statistiche/popolazione-andamento-demografico/
years_cal_mes = c(1861, 1871, 1881, 1901, 1911, 1921, 1931, 1936, 1951, 1961, 1971, 1981, 1991, 2001, 2011, 2014)	
calabria = c(323862, 354332, 375528, 437209, 470400, 525304, 565518, 578262, 639471, 609140, 578323, 573093, 576693,
             564223, 550967, 557993)
messina = c(399688, 426402, 467233, 550895, 545974, 613028, 605456, 627093, 667963, 685260, 654703, 669323, 646871,
            662450, 649824, 645296)
cal_mes = calabria + messina
conv_factor_i = italy[length(italy)] / cal_mes[length(cal_mes)]
#
#
#
census_pop_size = cal_mes
census_pop_today = census_pop_size[length(census_pop_size)]
census_pop_past = census_pop_size[1]
#
#
#
# *****************************
# *     tMRCA computation     *
# *****************************
#
N_e = census_pop_today / 100
#G = 25
n = 3
mu = 0.062633
#
xs = seq(0, 1000, length.out = 10001)
ys = sapply(xs, E_of_T_given_D, N_e, 1, n, mu)
mean_tMRCA = sum(xs*ys) / 10 # cause 1000 / 10001
areas = calculate_area(xs, ys)
lower = max(xs[c(areas<0.025, FALSE)])
upper = min(xs[areas>0.975])
select = lower <= xs & xs <= upper
#
pdf('tMRCA_generations.pdf', width=6, height=5)
plot(xs, ys, type='l', xlim=c(0, 100), xlab='generations', ylab='density')
polygon(c(lower, xs[select], upper), c(0, ys[select], 0), col=add.alpha(0.15), border=NA)
text(60, 0.05, pos=4, paste('mode =', round(xs[which(ys==max(ys))], 2)))
text(60, 0.046, pos=4, paste('mean =', round(sum(xs*ys)/10, 2)))
text(60, 0.042, pos=4, paste('95% CI = [', lower, ', ', upper, ']', sep=''))
dev.off()
#
#
#
# N_e doesn't really matter!
# N_es = seq(100, 2000, length.out = 101)
# plot(N_es, sapply(N_es, function(i) sum(xs*sapply(xs, E_of_T_given_D, i, 1, n, mu)) / 10), type='l')
#
#
#
# *****************************
# *       Forward sims        *
# *****************************
#
pop_size_at_time_t <- function(t, b, r) {
    (italy[length(italy)]*b^(b^(r*(t-years_italy[length(years_italy)]))-1)) / conv_factor_i
}
#
if (!file.exists(paste('sims.Robj'))) {
    nsim = 5*10^5
    xs = seq(1, 1000, length.out = 10000)
    ys = sapply(xs, E_of_T_given_D, N_e, 1, n, mu)
    ngen = round(replicate(nsim, sample(xs, 1, prob = ys))) # prob is normalised to 1 internally
    Gs = runif(nsim, min = 25, max = 31)
    t_at_gen = 2015 - ngen*Gs
    N = pop_size_at_time_t(t_at_gen, b, r)
    #
    # simulate allele counts in WF populations initially 1 mutant allele
    init = 1
    #
    nbgen = round(max(ngen))
    freq <- matrix(init, nsim, nbgen+1)
    #
    for (i in 1:nbgen) {
        Nnew = pop_size_at_time_t(t_at_gen+i*Gs, b, r)
        freq[,i+1] <- rbinom(nsim, round(Nnew), freq[,i]/N)
        N = Nnew
    }
    res = rep(0, nsim)
    for(i in 1:nsim) { res[i] = freq[i, ngen[i]+1] }
    res_ = res[res >= 3]
    hpd = p.interval(res_, HPD=TRUE, MM=FALSE, prob=0.95, plot=FALSE, PDF=FALSE)
    lower = as.numeric(hpd[1, 1])
    upper = as.numeric(hpd[1, 2])
    #
    dens = density(res_, from=0, to=250)
    #
    save(dens, res_, lower, upper, file = paste('sims.Robj'))
    #
} else {
    load(file = 'sim.Robj')
}
#
#
#
pdf('carriers.pdf', width=6, height=5)
plot(dens$x, dens$y, type='l', xlim=c(0, 120), xlab='nb. of carriers', ylab='density', main='')
select = lower <= dens$x & dens$x <= upper
polygon(c(lower, dens$x[select], upper), c(0, dens$y[select], 0), col=add.alpha(0.15), border=NA)
text(70, 0.08, pos=4, paste('mode =', round(dens$x[which(dens$y==max(dens$y))])))
text(70, 0.074, pos=4, paste('mean =', round(mean(res_))))
text(70, 0.068, pos=4, paste('95% HPD = [', lower, ', ', upper, ']', sep=''))
dev.off()
#







#
#  **************************  #
#  **************************  #
#  *** REDO WTH EXP MODEL ***  #
#  **************************  #
#  **************************  #
#
# LS-fitting of exponent
ers = exp(seq(-6.5, -12, length.out=1001))
#
lss = vector(mode = 'numeric', length = 1001)
for (i in 1:1001) {
    lss[i] = (sum(((italy[length(italy)]^(ers[i]*(years_italy-years_italy[length(years_italy)])+1)) - italy)^2))^(1/2) }
er = ers[which(lss==min(lss))]
#
pdf('population_size_param_fitting_EXP.pdf', width=5, height=4)
plot(ers, lss, type='l', xlab='r', ylab='RSSD')
lines(c(er, er), par('usr')[3:4])
dev.off()
#
#
#
#pdf('population_size_EXP.pdf', width=7, height=6)
pdf('population_size_EXP_only.pdf', width=7, height=6)
plot(years_italy, italy, type='l', ylim=c(0, 63000000), xlim=c(500, 2014),
     xlab='year AD', ylab='census population size Italy', yaxs='i')
#
xs = years_italy[1]:years_italy[length(years_italy)]
ys = sapply((years_italy[1]-years_italy[length(years_italy)]):0, function(x) italy[length(italy)]^(er*x+1))
lines(xs, ys, col='grey', lty='solid')
#
xs = 500:years_italy[1]
ys = sapply((years_italy[1]-years_italy[length(years_italy)]-500):(years_italy[1]-years_italy[length(years_italy)]),
            function(x) italy[length(italy)]^(er*x+1))
lines(xs, ys, col='red', lty='solid')
#
legend('topleft', c('demographic data', 'exponential fit', 'extrapolation'),
                    lty=c('solid', 'solid', 'solid'), lwd=c(1, 1, 1),
                    col=c('black', 'grey', 'red')) # , bty = 'n', cex=0.9
dev.off()
#
#
#
pop_size_at_time_t <- function(t, er) {
    (italy[length(italy)]^(er*(t-years_italy[length(years_italy)])+1)) / conv_factor_i
}
#
if (!file.exists(paste('sims_EXP.Robj'))) {
    nsim = 10^5 # 5*
    xs = seq(1, 1000, length.out = 10000)
    ys = sapply(xs, E_of_T_given_D, N_e, 1, n, mu)
    ngen = round(replicate(nsim, sample(xs, 1, prob = ys))) # prob is normalised to 1 internally
    Gs = runif(nsim, min = 25, max = 31)
    t_at_gen = 2015 - ngen*Gs
    N = pop_size_at_time_t(t_at_gen, er)
    #
    # simulate allele counts in WF populations initially 1 mutant allele
    init = 1
    #
    nbgen = round(max(ngen))
    freq <- matrix(init, nsim, nbgen+1)
    #
    for (i in 1:nbgen) {
        Nnew = pop_size_at_time_t(t_at_gen+i*Gs, er)
        freq[,i+1] <- rbinom(nsim, round(Nnew), freq[,i]/N)
        N = Nnew
    }
    res = rep(0, nsim)
    for(i in 1:nsim) { res[i] = freq[i, ngen[i]+1] }
    res_ = res[res >= 3]
    hpd = p.interval(res_, HPD=TRUE, MM=FALSE, prob=0.95, plot=FALSE, PDF=FALSE)
    lower = as.numeric(hpd[1, 1])
    upper = as.numeric(hpd[1, 2])
    #
    dens = density(res_, from=0, to=250)
    #
    save(dens, res_, lower, upper, file = paste('sims_EXP.Robj'))
    #
} else {
    load(file = 'sim_EXP.Robj')
}
#
#
#
pdf('carriers_EXP.pdf', width=6, height=5)
plot(dens$x, dens$y, type='l', xlim=c(0, 200), xlab='nb. of carriers', ylab='density', main='')
select = lower <= dens$x & dens$x <= upper
polygon(c(lower, dens$x[select], upper), c(0, dens$y[select], 0), col=add.alpha(0.15), border=NA)
text(120, 0.045, pos=4, paste('mode =', round(dens$x[which(dens$y==max(dens$y))])))
text(120, 0.042, pos=4, paste('mean =', round(mean(res_))))
text(120, 0.039, pos=4, paste('95% HPD = [', lower, ', ', upper, ']', sep=''))
dev.off()
#









