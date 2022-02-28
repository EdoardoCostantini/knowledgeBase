# Project:   knowledgeBase
# Objective: Using the tryCatch function
# Author:    Edoardo Costantini
# Notion:
# Other ref: https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
# Created:   2022-02-28
# Modified:  2022-02-28

# Define a function that reads what is in a url with tryCatch statements
readUrl <- function(url) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one 
      # R expression in the "try" part then you'll have to 
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      
      message("This is the 'try' part")
      
      objct <- readLines(con = url, warn=FALSE)
      length(objct)
      # The return value of `readLines()` is the actual value 
      # that will be returned in case there is no condition 
      # (e.g. warning or error). 
      # You don't need to state the return value via `return()` as code 
      # in the "try" part is not wrapped insided a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error = function(x) {
      # message(paste("URL does not seem to exist:", url))
      # message("Here's the original error message:")
      # message(cond)
      # Choose a return value in case of error
      err <- paste0("Error: ", x)
      return(err)
    },
    warning=function(x) {
      # message(paste("URL caused a warning:", url))
      # message("Here's the original warning message:")
      # message(cond)
      # Choose a return value in case of warning
      warn <- paste0("Warning: ", x)
      return(warn)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      message(paste("Processed URL:", url))
      message("Some other message at the end")
    }
  )    
  return(out)
}

# Define some URL inputs
urls <- c(
  "http://stat.ethz.ch/R-manual/R-devel/library/base/html/connections.html",
  "http://en.wikipedia.org/wiki/Xz",
  "xxxxx" # will give an error
)

# Use the reading function
y <- lapply(urls, readUrl)
y
