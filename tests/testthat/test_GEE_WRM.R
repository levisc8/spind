context('Test GEE and WRM outputs')

test_that('Coefficient estimates + SEs are as expected', {
  skip_on_cran()
  data(musdata)
  coords <- musdata [ ,4:5]

  formula <- as.formula(musculus ~ pollution + exposure)

  mgee <- GEE(formula, family = "poisson", data = musdata,
              coord = coords, corstr = "fixed", plot = F, scale.fix = FALSE)


  mwrm <- WRM(musculus ~ pollution + exposure, family = "poisson",
              data = musdata, coord = coords, level = 1, plot = F)
  expect_equal(as.vector(mgee$b), c(-1.904747, 3.362160, -1.463483),
               tolerance = 1e-7)
  expect_equal(as.vector(mwrm$b), c(-1.935992, 3.184136, -1.228579),
               tolerance = 1e-7)
  expect_equal(as.vector(mgee$s.e.), c(1.3109113, 0.9141623, 0.8801014),
               tolerance = 1e-7)
  expect_equal(as.vector(mwrm$s.e.), c(1.917704, 1.225106, 1.506331),
               tolerance = 1e-7)

  expect_true(is.list(mgee))
  expect_is(mgee, 'GEE')
  expect_true(is.list(mwrm))
  expect_is(mwrm, 'WRM')
})

test_that('predict.GEE and .WRM interact properly with parent functions',{
  skip_on_cran()
  data(musdata)
  coords <- musdata [ ,4:5]

  formula <- as.formula(musculus ~ pollution + exposure)

  mgee <- GEE(formula, family = "poisson", data = musdata,
              coord = coords, corstr = "fixed", plot = F, scale.fix = FALSE)


  mwrm <- WRM(musculus ~ pollution + exposure, family = "poisson",
              data = musdata, coord = coords, level = 1, plot = F)
  expect_match(class(predict.GEE(mgee, newdata = musdata)), 'numeric')
  expect_match(class(predict.WRM(mwrm, newdata = musdata)), 'numeric')
})

test_that('fails correctly', {
  skip_on_cran()

  expect_error(GEE(formula, family = "poisson", data = musdata,
               coord = coord, corstr = "fixed", plot = F, scale.fix = FALSE),
               "formula: specified notation is missing")

  formula <- as.formula(musculus ~ pollution + exposure)

  coord <- cbind(rnorm(400), runif(400))

  expect_error(GEE(formula, family = "poisson", data = musdata,
                   coord = coord,plot = F),
               "coordinates not integer")

  coord <- cbind(round(runif(300)), round(runif(300)))

  expect_error(GEE(musculus ~ pollution + exposure, family = "poisson",
               data = musdata, coord = coord, plot = F),
               "length of data does not match length of coordinates")


  coord <- cbind(rnorm(400), runif(400))

  expect_error(WRM(formula, family = "poisson", data = musdata,
                   coord = coord,plot = F),
               "coordinates not integer")

  coord <- cbind(round(runif(300)), round(runif(300)))

  expect_error(WRM(musculus ~ pollution + exposure, family = "poisson",
                   data = musdata, coord = coord, plot = F),
               "error in dimension")

})



