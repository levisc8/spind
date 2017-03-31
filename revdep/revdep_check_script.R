library(devtools)

revdep_check(pkg='test2',recursive=T,bioconductor=T)
revdep_check_print_problems(pkg='test2')

readRDS('test2/revdep/checks.rds')
