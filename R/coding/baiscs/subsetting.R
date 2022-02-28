# Project:   knowledgeBase
# Objective: Understanding subsetting
# Author:    Edoardo Costantini
# Notion:
# Other ref: http://adv-r.had.co.nz/Subsetting.html
# Created:   2022-02-28
# Modified:  2022-02-28

# Data Types --------------------------------------------------------------
# Atomic vector
x <- c(2.1, 4.2, 3.3, 5.4)

x[c(3, 1)]
x[order(x)]
x[c(1, 1)] # Duplicated indices yield duplicated values
x[c(2.1, 2.9)] # trucated to integer

x[c(-1, -3)] # omits elements
x[c(-1, 2)]  # can't mix positive and negative

x[c(TRUE, TRUE, NA, FALSE)]
x[c(TRUE, TRUE, NA)] #recycles

# Character vectors
(y <- setNames(x, letters[1:4]))

# Matrices and arrays -----------------------------------------------------

(a <- matrix(1:9, nrow = 3))
  colnames(a) <- c("A", "B", "C")
  
a[c(TRUE, FALSE, TRUE), c("B", "A")]
a[0, 1:2]

(vals <- outer(1:5, 1:5, FUN = "paste", sep = ","))
vals[c(4, 15)] # can use a vector to subset a matrix based on "count" by column

select <- matrix(ncol = 2, byrow = TRUE, c(
  1, 1,
  3, 1,
  2, 4
))
vals[select] # uses first column to define the row to select, and second to define the column to select

# Data Frames -------------------------------------------------------------

df <- data.frame(x = 1:3, y = 3:1, z = letters[1:3])

df[df$x == 2, ]

df["x"]
  str(df["x"])
df[, "x"]
  str(df[, "x"])
df[["x"]]
  str(df[["x"]])
  
# Subsetting Operators ----------------------------------------------------

a <- list(a = 1, b = 2)

a[1]
a[[1]]

b <- list(a = list(b = list(c = list(d = 1))))
b[["a"]][["b"]][["c"]][["d"]]
b[[c("a", "b", "c", "d")]] # indexes recursively  <- !!!!!

# Simplifying vs. preserving subsetting -----------------------------------
# Look at table under the homonymous section in the source link

z <- factor(c("a", "b"))
z[1]
z[1, drop = TRUE] # drops unused levels

a <- matrix(1:4, nrow = 2)
a[1, , drop = FALSE] # keeps
a[1, ] # simplifies



