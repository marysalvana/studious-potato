
source("./Functions/load_packages.R")
source("./Functions/cov_func.R")
source("./Functions/auxiliary_functions.R")

sourceCpp("./Functions/spatially_varying_parameters2-IBEX.cpp")
sourceCpp("./Functions/distR.cpp")


root <- gsub('.{7}$', '', getwd())


