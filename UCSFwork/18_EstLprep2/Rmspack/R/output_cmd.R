### Function to call the program output_cmd lines for estlike program ###



`output_cmd` <-
  function(datafile, param)
  {
    return(.Call("outputcmd",datafile,param, PACKAGE="output_cmdR"))
  
  }
