.First.lib <- function(libname, package) {
	library.dynam("msR", package)
	library.dynam("stats_popR", package)
   library.dynam("output_cmdR", package)
}
