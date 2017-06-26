context('adjusted.actuals')

test_that('adjusted.actuals returns correct values',{

  data(hook)
  data<- hook[,1:2]
  coord<- hook[,3:4]
  aa<-adjusted.actuals(data,coord,plot.maps=F)

  expect_match(class(aa),'numeric')

  expect_equal(length(aa), dim(data)[1])
  expect_true(max(aa) <= 1)
  expect_true(min(aa) >= 0)
})

test_that('fails correctly',{

  data(hook)
  data <- hook[,1:2]
  coord <- cbind(seq(0,50, length.out = dim(hook)[1]),
                 seq(0,50, length.out = dim(hook)[1]))

  expect_error(adjusted.actuals(data,coord,plot.maps=F),
               'coordinates not integer')

  coord <- cbind(1:200,1:200)
  expect_error(adjusted.actuals(data,coord,plot.maps=F),
               'coordinates and data have different dimensions')

})