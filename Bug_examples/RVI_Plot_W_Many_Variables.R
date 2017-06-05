library(MASS)
data(birthwt)
x <- rep(1:14, 14)
y <- as.integer(gl(14,14))
coords <- cbind(x[-(190:196)], y[-(190:196)])
formula <- formula(low ~ age + lwt + race + smoke + ftv +  bwt)
rvi.plot(formula, "gaussian", birthwt, coords, maxlevel = 4,
         detail = TRUE, wavelet = "d4")
rvi.plot(carlina.horrida ~ aridity + land.use, "poisson",
         carlinadata, coords, maxlevel = 4,
         detail = TRUE, wavelet = "d4")
