<!--
% Introduction to R
% Feb 2021
% Thu Nov 25 14:35:25 2021
-->

<!-- vim: tw=66 nosmartindent fmr=<<<,>>> cms=<!--%s-->

<!-- <<< readme.md -->
# readme.md

# Advanced R 2 and 3 Feb 2022

# [https://bit.ly/3pPBwwC](https://bit.ly/3pPBwwC)

## Requirements

* A reasonably capable computer.
* Installation of R version 4. (4.1.2 is the most recent).
* Installation of a recent version of RStudio (2021.09.2-382 is
  the most recent.

## Timings for both Wed 2 Feb and Thu 3 Feb

---------------------        ---------------
Course                       0930  to   1100
Coffee break                 1100  to   1130
Course                       1130  to   1300
Lunch (not provided)         1300  to   1400
Course                       1400  to   1530
Coffee break                 1530  to   1600
Course                       1600  to   1730
---------------------        ---------------

## As soon as you join the MS Teams meeting

* Point your web browser to the url shown at the top
  which is the same as

- [https://bit.ly/3pPBwwC](https://bit.ly/3pPBwwC) which is the same as 
  [https://streptomyces.github.io/advancedR/](https://streptomyces.github.io/advancedR/).



## How will this course work?

- I will talk about commands and syntax.
- You will run some commands along with me to see them
  in action.
- I will explain the commands and the syntax you have
  just seen in action.

### Sometimes there will be things for you to do on your own.

+ These tasks will be described in plain english.
+ You will have to think how to achieve them in R using
   as many steps as you need.
+ Finally, I will show you how I would have done the
  tasks but keep in mind that, if you got the right results,
  your way of doing them may also be correct.

## Limitations

* I use R on a daily basis but not RStudio.
* The focus of this course is on R syntax and
  techniques rather than statistics.
* The aim is to clarify the workings of R so that you
  can build your own data solutions and workflows.
* This is my understanding and my view of R. It is
  not complete but I do believe that it is not wrong.
* As you progress, your view will not be identical to
  mine. Do try to keep it consistent though.

<!-- >>> -->

<!-- <<< rstudio.md -->
# rstudio.md

# RStudio

* Start RStudio.

## Getting the scripts and data we will be using

In Rstudio, use the drop down menu to do

    File -> New File -> R Script

In the blank R script we will write and run
the indented lines shown below. Please do this
with me. Resist the temptation to charge ahead.

~~~ {.r}
    setwd("u:/")
    unlink("Rtrain", recursive = TRUE)
    dir.create("Rtrain")

    setwd("Rtrain")

    getwd()

    unlink("*")

    list.files()
    
    download.file("https://github.com/streptomyces/advancedR/raw/master/adr.zip",
    "adr.zip")

    list.files()

    unzip("adr.zip")
    
    setwd("adr")

    list.files()
~~~

### Four panes

* Source editor
* Console
* Environment and History
* File, Plots, Help etc.

### Getting started

* Type commands in the _source_ editor.

* To run a command from the _source_ editor place the
  cursor anywhere on that line and press `Control-Enter`.

* You can also select multiple lines and then press
  `Control-Enter` to run all the selected lines.

* Finally, you can use `Control-Shift-Enter` to _source_
  (i.e. run) the entire script.

* The action happens in the _Console_ pane. i.e. any
  output from R is shown in the console frame. You can
  also type commands directly in the console.

* R keeps a history of your commands which you can see in
  the history tab of RStudio.

* You can select and run commands from the history tab as
  well. Just press enter on the current line or select some
  lines and press enter.

* However you run a command, it is like typing it into
  the console and pressing Enter.

* Matching parentheses and quotes are automatically
  inserted. You can disable this in options. I do.
  (`Tools -> Global Options -> Code -> Editing`)

#### My preferred arrangement of RStudio panes

Use the drop down menus to do the following.

    Tools -> Global Options -> Pane layout

Arrange to have *Source* on the top left and *Console*
on the top right. With this arrangement we can minimise
two bottom panes most of the time and have more of
screen space.

<!-- >>> -->

<!-- <<< r.md -->
# r.md

# R

* Comments begin with #. Everything after a # until the end
  of the line is ignored by R.

### Some syntax identifying features

~~~ 
<- and =        Assignment                  x <- 42            
                                            x = 42

->              Rightward assignment        42 -> x

Function        Unquoted word followed      mean()
call            by parentheses

Object name     Unquoted word               gene.lens

String          Quoted alphanumeric         "whiA"
literals        characters

Numeric         Unquoted digits and         2345, 42L, 3.14
literals        scientific notation         1e6, 1e-6
~~~


* Parentheses, (), are required in function calls even if
you are not passing any arguments to the function being
called. `ls()` works but `ls` shows you the definition of
the function `ls`. Functions are objects too.

* The parentheses makes it easy to identify individual
function calls in long and complex R statements where
function calls are embedded within other function
calls.

* Unquoted words which are not reserved words are assumed
to be object names.

* Almost no restrictions on names but see `help("Quotes")` and
  `help("make.names")`. Use meaningful names.

~~~ 
    mean <- c(2,3,4,5); # example of bad name.
    mean(mean);
~~~

* Commands can continue over multiple lines.

* Semicolons are only needed if you put two commands on one line.

~~~ 
    x <- 3; y <- x * 3; x; y
~~~

* If you get inside a complicated command which you
  cannot finish, try Escape (Control-C in Linux) to bail out. This
  usually happens because of unmatched parentheses or quotes.

## Names are free!

- While doing the tasks during this course you will not always be
  explicitly asked name your objects.

- Name your objects  as needed. e.g. if a task says make a vector
  of the numbers 1 to 10 do

~~~ 
x <- seq(1,10);
~~~

instead of

~~~ 
seq(1,10);
~~~

- Then you will have a named you can use in the
  following steps of the task.

- Use as many intermediate objects as you need to complete the
  task.

- For example, if asked to generate 20 normally distributed
  numbers, take their log2 and, find the median.

~~~ 
x <- rnorm(20, mean = 30, sd = 1.2)
xl2 <- log2(x);
median(xl2);
~~~

The same in one step.

~~~ 
median(log2(rnorm(20, mean = 30, sd = 1.2)))
~~~

*Work at a level of complexity that you are comfortable
with.*

## Data analysis cycle

#### Experiment

1. Experiment
2. Data
3. Data analysis
4. Go to 1.

#### Data analysis

1. Data file(s) from experiment
2. Cleaning up and formatting
3. Reading into objects
4. Apply functions to manipulate the objects
5. Repeat above as required
6. Make plots
7. Write out to files for further work or communication

<!-- >>> -->

<!-- <<< datatypes.r -->
# Datatypes

~~~ {.r}

########################
### Basic data types ###
########################

# Character (strings of characters)
# Numbers
  ## Integers
  ## Real (double precision floating point numbers)
# Logical. TRUE or FALSE only.
# Complex
# Raw

# We will not bother with the Complex and Raw types.

# Mostly we don't have to worry whether our numbers
# are integers or doubles but there are some situations
# where we should. There are functions to coerce one
# into another if needed.

# Explain x <- 42L

x <- 42;

y <- 42L;

typeof(x);

typeof(y);

### More data types ###

# Lists
# Functions (closures)
# Environments
# Builtin

~~~

<!-- >>> -->

<!-- <<< objects.r -->
# Objects

~~~ {.r}

#######################
### Data structures ###
#######################

# Types which are combinations of the basic data types.

# Can be arbitrarily complex.

# Several provided by R as classes.

###############
### Classes ###
###############

# Data structures with associated functions (methods).

# e.g. the class data.frame has the following methods.

    methods(class = "data.frame");

# [             [[            [[<-          [<-           $<-          
# aggregate     anyDuplicated as.data.frame as.list       as.matrix    
# by            cbind         coerce        dim           dimnames     
# dimnames<-    droplevels    duplicated    edit          format       
# formula       head          initialize    is.na         Math         
# merge         na.exclude    na.omit       Ops           plot         
# print         prompt        rbind         row.names     row.names<-  
# rowsum        show          slotsFromS3   split         split<-      
# stack         str           subset        summary       Summary      
# t             tail          transform     type.convert  unique       
# unstack       within 

# A class is a recipe for making objects and "constructor"
# methods are usually provided.

###############
### Objects ###
###############

# Objects are "instances" of classes.

# Objects of a class have all the methods of the class
# available to them.

# Objects can be of more than one class. Then they have
# methods of all the classes available to them.

# There are no simple variables in R.

# Even if you need to store just one value you use a data
# structure that is capable of storing multiple values.

~~~

<!-- >>> -->

<!-- <<< variables.r -->
# Variables

~~~ {.r}

####################################
### Names, symbols and, bindings ###
####################################

# Names (or symbols) are handles or labels for objects
# through which we can access the objects.

# After an assignment like: x <- 42

#     the name "x" has the value 42.

#     the value 42 has a label attached to it with "x"
#     written on the label.

#     The name "x" is bound to the value 42.

#     Explain unbound (non-existent) names.

# It is not accurate to imagine x to be a container
# containing the value 42.

# Nothing prevents two different names from refering to
# the same object. In fact this is what happens after

y <- x;

# What happens when rm is invoked?

rm(x);

~~~

<!-- >>> -->

<!-- <<< vectors.r -->
# Vectors

~~~ {.r}

###############
### Vectors ###
###############

# Ordered collection of values.

# Contiguous cells containing data.

# All values have to be of the same type.

# They can hold values of any of the basic types.

# One way of making a vector is by using the c() function.
# (concatenate).

x <- c(10, 20, 30)

# Cells are accessible through subscripting / indexing.
# e.g. if vector x contains 10, 20 and 30 then

x[1] # refers to 10,
x[2] # refers to 20 and,
x[3] # refers to 30

# Indexing begins from one, not zero.


x <- c(2.9, 4.1, 3.9, 4.5, 3.7, 45.3, 21.6);

x;

x[1];

x[7];

x[4:6]   # ":" is the range operator. Step by 1 only.


# The simplest data type in R is a vector.

y <- 5.4; # is the same as y <- c(5.4);

y;

y[1];

z <- c("ftsZ", "sigE", "bldN", "whiA", "whiB", "rdlA", "chpA");
z;
z[3];

length(x)
length(z)


### Named Vectors ###

# Elements of a vector can be named

names(x) <- z;
x;

x["whiA"]; # is the same as x[4];

v <- c("whiA", "rdlA");

x[v];

#################################
### Do the following yourself ###
#################################

# 1. Make the vectors x and z as discussed above.

x <- c(2.9, 4.1, 3.9, 4.5, 3.7, 45.3, 21.6);
z <- c("ftsZ", "sigE", "bldN", "whiA", "whiB", "rdlA", "chpA");

# 2. Name the elements of x as the strings in z.

# 3. Examine the output of

unname(x)

# 4. Make a copy of x in vecx. We will use it in the next
# script.


~~~

<!-- >>> -->

<!-- <<< functions.r -->
# Functions

~~~ {.r}

#################
### Functions ###
#################

# Functions do things for us.

# Need arguments to work on.

# May assume default arguments if none are provided.

# May return something which you may wish to store as a
# new object.

# If you are not assigning the returned value to an object
# then the returned value may be printed on the console.
# invisible()

# Some functions used primarily for their side-effects
# rather than their return value.

# User defined functions.

bmi <- function(kilograms, metres) {
  return(kilograms / metres ** 2);
}

bmi(m = 1.65, 65);

### Components of functions ###

#  1. Formal argument list.
#  2. Body
#  3. Environment


### Arguments

# Arguments are values or objects a function acts on. The
# process of getting arguments into functions is called
# argument passing.

# Supplied arguments are matched to formal arguments in a
# three pass process.

#  1. Exact matching on tags
#  2. Partial matching on tags
#  3. Positional matching.

# We will use the seq() function of R as an example.

# positional arguments
seq(2, 10);

seq(2, 10, 2);

# Positional and named arguments
seq(2, 10, by = 2);

# Named arguments only
seq(from = 2, to = 10, by = 2);

# Positional and named arguments
seq(2, to = 10, by = 2);

# Positional and named arguments
seq(2, by = 2, 10);

# Partial matching of argument names
seq(2, fr = 4, to = 20);


#######################################
### Calling bmi() is different ways ###
#######################################
bmi(65, 1.68);
bmi(kilograms = 65, metres = 1.68)
bmi(metres = 1.68, kilograms = 65)
bmi(kilo = 65, met = 1.68)
bmi(m = 1.68, k = 65)
bmi(m = 1.68, 65)


####################################
### Return values from functions ###
####################################

# If a function returns something and you want to keep it
# you have to assign it to a (new) object.

x <- seq(2, 10);
x;
y = seq(3, 30, by = 3);
y;


# Some functions do not return anything (NULL).
x <- seq(10, 100, by = 5);
cat(x, "\n");
catret <- cat(x, "\n");
catret;

# Demonstrate invisible here.

bmii <- function(kilograms, metres) {
  invisible(kilograms / metres ** 2);
}

bmi(m = 1.65, 65);
bmii(m = 1.65, 65);
bmx <- bmii(m = 1.65, 65);
bmx

##############################################
### A note about assignment operators in R ###
##############################################

x <- seq(10, 100, by = 5);
x
# Works, but for the wrong reason.
x <- seq(10, 100, by <- 5);
x

# Fails. Wrong sign in by argument.
rm(x);
x <- seq(10, by <- 5, 100);
x

# Works. But unintended result.
rm(x);
x <- seq(10, by <- 5, -1);
x

# Alt-minus key combination gives you the
# assignment operator in RStudio
 
# Use <- or = on the command line.
# Always use = in function calls and function definitions.

# Function calls as arguments

rep(4, times = 4);

rep(seq(1,4), times = 4);

rep(seq(1,4), each = 4);

#################################
### Do the following yourself ###
#################################

# 1. Run the code below.
#
x <- y <- z <- 42
#
# 2. Check the values of x, y and z.


# 1. Generate the sequence -20 to 50 increasing in steps
# of 5 and name it ds.

# 2. Check the value of ds.

######################################################

# 1. Write a function named cel2fah which takes one
# argument, the temperature in degrees Celsius, and returns
# the same temperature in Fahrenheit.

# 2. Hints:
# - Multiply Celsius by 1.8 then add 32 to get Fahrenheit.

cel2fah <- function(cel) {
fah <- (cel * 1.8) + 32;
return(fah)
}


### Write your function here ###

# 3. Test your function by calling it as below.

#    a. cel2fah(100)
#    b. cel2fah(0)
#    c. cel2fas(ds)

# 4. Try the below.

curve(cel2fah, -20, 50, xlab = "Celsius",
      ylab = "Fahrenheit", lwd = 3, col = "royalblue")

# Functions are objects. We passed a function as an argument
# to another function.

# Another example of passing a function as an argument to
# another function.

deg2rad <- function(deg) {
return(deg * 2 * pi / 360);
}
curve(sin, deg2rad(0), deg2rad(360), col = "red", lwd = 3,
      ylab = "f(x)", xlab = "x in radians");
curve(cos, deg2rad(0), deg2rad(360), col = "blue", lwd = 3, add = T);


########################
### ... (three dots) ###
########################

dotdemo <- function(x, ...) {
  cat("x is: ", x, "\n");
  dots <- list(...);
  cat("Then: ", dots[[1]], "\n");
}
dotdemo("stuff", "morestuff");


dotdemo <- function(x, ...) {
  cat("x is: ", x, "\n");
  dots <- list(...);
  cat("y: ", dots$y, "\n");
  cat("Then: ", dots[[2]], "\n");
}
dotdemo(x = "stuff", y = "morestuff", 23);


##############################
### Operators as functions ###
##############################

x <- 2;
y <- 3;

x + y

'+'(x,y)
'*'(x,y)


'['(vecx, 4);
'['(vecx, "whiA");

~~~

<!-- >>> -->

<!-- <<< scopeenv.r -->
# Scope and Environment

~~~ {.r}

##############################
### Scope and Environments ###
##############################

# Scope of a variable is where, in a running program, it is
# accessible.

x <- 42;

# Functions and environments define the scope of variables
# in them. In an R session the default environment is the global
# environment named .GlobalEnv

# Lexical scope: Functions are evaluated in the environment
# in which they are defined.
#
# Closures
#
# The environment provides values for any unbound symbols in
# in the body of a function when the function is evaluated.

rm(list = ls());

rate <- 0.2;
withvat <- function(exvat) {
  return(exvat + (exvat * rate))
}


vatfuncgen <- function(vatrate) {
  vatfun <- function(exvat) {
    return(exvat + (exvat * vatrate))
  }
  return(vatfun)
}

standard <- vatfuncgen(0.2);
reduced <- vatfuncgen(0.1);

standard(100)
reduced(100)

environment(withvat);
environment(standard);
environment(reduced);

ls(envir = environment(standard));
ls(envir = environment(reduced));
ls(envir = .GlobalEnv);

parent.env(environment(standard));

superlow <- standard;
sl.env <- new.env();
sl.env$vatrate <- 0.05;
environment(superlow) <- sl.env;

superlow(100);
parent.env(environment(superlow));
parent.env(parent.env(environment(superlow)));
parent.env(parent.env(parent.env(environment(superlow))));

### assign() and get() ###

# Standard VAT goes up to 0.25.

assign("vatrate", 0.25, envir = environment(standard));

standard(100);

get("vatrate", envir = environment(superlow));

# The effective environment is almost always a nesting of
# environments.

#################################
### The search path. search() ###
#################################

search();

attach(environment(standard));

search();

ls(); # Notice that there is no vatfun in the listing.

vatfun(100); # But you can still call vatfun.

~~~

<!-- >>> -->

<!-- <<< queryingobjects.r -->
# Querying Objects

~~~ {.r}

######################################
### Querying the nature of objects ###
######################################

# Mode

# The mode of an object signifies the type of data in it.
# numeric, character, logical, complex and raw.
x <- seq(1, 5, by = 0.5);
mode(x);
y <- as.integer(x);
mode(y);

# Type and Storage mode

# These two are approximately the same.
storage.mode(x);
typeof(x);

storage.mode(y);
typeof(y);

# Class

# The class of an object is the kind of data structure it is. For
# vectors it is the same as the mode. list, data.frame, matrix,
# function. Class of an object determines how functions treat it.

class(x);
class(y);

# str() is useful to display the internal structure of R objects.
# Especially useful for data frames and more complex objects.

hx <- hist(rnorm(300), plot = F);
str(hx)


~~~

<!-- >>> -->

<!-- <<< packages.r -->
# Packages

~~~ {.r}

################
### Packages ###
################

# R classes and functions written by others for you to
# use.

# Provide namespaces which roughly translate into
# environments on the search() hierarchy thus avoiding
# symbol clashes.

# Packaged to be loaded on demand to minimise resource
# use.

# When you load a package additional functions defined
# in the package become available for you to use.
# Sometimes existing functions might be over-ridden.

# CRAN has a repository of R packages.

# Bioconductor is a repository of packages for
# bioinformatics. Bioconductor has its own installer.

### Installing edgeR ###
# We don't need to do this because edgeR is already installed
# on these machines by computing.

#   if (!requireNamespace("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#   BiocManager::install("edgeR", version = "3.8")

#################################
### Using packages. library() ###
#################################

library("edgeR");

# Load the namespace of the package with name
# ‘package’ and attach it on the search list.

# The reverse is

detach("package:edgeR", unload = TRUE)

# Sometimes you will see a function being called as

#        packagename::funcname()
# e.g.   readxl::read_excel()

# In such cases the package is loaded (into memory) but it
# is not attached to the search path. This is good for one
# time use of a function and to avoid namespace conflicts
# resulting in masking.

# Also good for accessing masked functions.

# e.g. stats::filter() to use the filter function in the
# stats package while it is masked by the one in dplyr
# package.




~~~

<!-- >>> -->

<!-- <<< documentation.r -->
# Documentation

~~~ {.r}

#############################
### Reading documentation ###
#############################


# Go to the "Help" tab in the bottom right
# pane then search for "seq". You can also
# do the below to see the help in the "Help"
# tab.

help("seq")
help("plot")


# help.search() shows all the packages which
# have anything to do with "sequence".

help.search("sequence")


### Sometimes there is runnable example code. ###

example("hist")


### There are lots of bundled example data sets ###

data() # lists available data sets.

data("cars") # loads the named data set.
ls()
head(cars)

### Finding methods ###

# All definitions for a generic method in the search
# path.

methods("plot"); # plot is a generic function.
?plot.histogram

# All methods available in a class.
methods(class = "histogram");
methods(class = "data.frame");

# Packages are not the same as classes. Classes are
# defined in packages and more than one may be defined
# in one package.

# Before calling library("edgeR");

library(help = "edgeR")

# Above should tell you what classes and functions are
# available in the edgeR package. Then you can see
# the help for those as below.

help("DGEList", package = "edgeR")
help("DGEList-class", package = "edgeR")


# After calling library("edgeR");
library("edgeR");
methods(class = "DGEList");
help("DGEList-class")
help("DGEList")
help("subsetting", package = "edgeR");

methods(class = "DGEGLM");
help("DGEGLM-class")
help("DGEGLM")

### Operator Precedence

help("Syntax");

### Reserved words

help("Reserved");

~~~

<!-- >>> -->

<!-- <<< recycling.r -->
# Recycling

~~~ {.r}

######################################
### Recycling in Vector operations ###
######################################


#################################
### Do the following yourself ###
#################################

# 1. Store the sequence from 1 to 5 in vector x.
# 2. Store the sequence from 21 to 25 in vector y.
# 3. Examine the result of y * x.

#############################################################
  
# 4. Store the sequence from 21 to 30 in vector z.
# 5. Confirm the lengths of x and z using the function length().
# 6. Examine the result of z * x.

#############################################################
  
# 7. Store the sequence from 1 to 7 in vector u.
# 8. Examine the result of z * u.



#############################################################
#############################################################
# In operations involving two vectors of unequal length, elements of the
# shorter vector get recycled.
# If the longer vector is not an integer multiple of the shorter vector
# you get a warning but the operation is valid and successful.


~~~

<!-- >>> -->

<!-- <<< commonfunc.r -->
# Common Functions

~~~ {.r}

#############################
### Some Common Functions ###
#############################

# c(): concatenate.

x <- c(2.17, 3.14, 13, 29, 48.356);
y <- c(200,300);
z <- c(x, y);
z;

# length(): Return the length of the named object.
length(x);
length(y);

# min() and max(): Minimum / Maximum element of an object.
min(x);
min(y);

# mean() and median():
mean(x);
median(y);

# summary(): Some key statistics about an object.
summary(y);


#####################
### Rounding etc. ###
#####################

# round()
pi
round(pi)
round(pi, 3);

# signif()
signif(pi)
signif(pi, 6);
signif(pi, 3);

# floor()
floor(pi)

# ceiling()
ceiling(pi)

############################
### paste() and paste0() ###
############################

# Always return a character vector.

x <- LETTERS[1:10]
y <- seq(1,10);

paste(x);
paste(y);   # A character vector.
typeof(paste(y));

paste(x, y);
paste0(x, y);  # No separator.

y <- seq(1,5);
paste0(x, y);  # y is shorter than x, hence recycling.

z <- c("-rep1", "-rep2");
paste0(x, y, z); # z is also recycled as needed.


# You can "collapse" the resulting vector to a single
# string.

paste0(x, y, collapse = " ");

#################################
### Do the following yourself ###
#################################

mtl <- LETTERS[1:8]
mtn <- seq(1,12)

# Using rep() and paste0() and the two vectors made
# above, make a vector of all the addresses on a 96
# well microtitre plate.

# A1 to A12 then B1 to B12 ... H1 to H12.

# Hint: call to rep() will go inside the call to paste0()
# and will also use the named argument "each".




~~~

<!-- >>> -->

<!-- <<< lists.r -->
# Lists

~~~ {.r}

#############
### Lists ###
#############

# Lists are generic vectors.

# Generally, when we say "vector" we mean "atomic
# vectors".

# Individual elements of a list can refer to any type of
# R object of any complexity, including other lists.

# This allows for the creation of objects of arbitrary
# complexity.

# Functions returning a lot of related information
# often return their results as lists.

x <- seq(1,20, by = 5);
y <- c("one", "two", "three");
i <- list(x, y);

i;
i[1];   # list
i[[1]]; # vector
class(i[1]);
class(i[[1]]);

names(i) <- c("numbers", "strings");

i["numbers"];
i[["numbers"]];
names(i);


j <- list(numbers = x, strings = y, serial = seq(1,20));

j[["serial"]][3:10];
j$serial;
j$serial[15:20];

################
### unlist() ###
################

unlist(j); # all atomic components of j in a vector
# Notice coercion above.

k <- list(x = x, y = seq(200, 400, by = 10));
k
unlist(k)
class(unlist(k))

#################################
### Do the following yourself ###
#################################

# 1. Make a vector of characters A to N by appropriately
# subsetting LETTERS.

# 2. Make a vector of integers from 10 to 20.

# 3. Make a list li, with two named elements, alpha
# containing the vector made in step 1 above and, beta
# containing the vector made in step 2 above.

# 4. Display the names of the elements in the list made
# above.

# 5. In one step, extract the fourth element from the
# alpha member of li.



~~~

<!-- >>> -->

<!-- <<< attributes.r -->
# Attributes

~~~ {.r}

##################
### Attributes ###
##################

# Objects can have attributes attached to them.

# Think of them as arbitrary name = value pairs.

# attributes()

# attr()

# The following are used by R.

# 1. Names. names().

# 2. Dimensions. dim(). Matrices and Arrays which we
# will see later on are vectors with the dim attribute
# attached to them.

# 3. Dimnames. dimnames(). If dimensions are named
# the names are stored in the dimnames attribute of the
# array or matrix.

# 4. Class. class(). The name of the class an object
# belongs to is stored in and attribute named class.

###################################################
### Step and think through the statements below ###
###################################################

x <- seq(1,10);
x
class(x);
class(x) <- c(class(x), "adr");
class(x);

attr(x, "purpose")
attr(x, "purpose") <- "testing";
attr(x, "purpose")
attr(x, "purpose") <- c(attr(x, "purpose"), "training")
attr(x, "purpose")

attributes(x);

~~~

<!-- >>> -->

<!-- <<< dataframes.r -->
# Data Frames and Tibbles

~~~ {.r}

###############################
### Data Frames and tibbles ###
###############################

# Like two dimensional tables where columns are vectors.

# Columns can be of different types.

# All values in any one column have to be of the same type.

# All columns have to be of the same length.

# Commonly used for getting data into R from text files.

# Tibbles don't have row names.

##############################
### Reading MS Excel files ###
##############################

# Package named readxl provides functions for reading
# Microsoft Excel files.

# If the package readxl is not installed then you can
# install it by doing the following.

# install.packages("readxl");

# Determine for format of the excel file. .xlsx or .xls
excel_format("data/file.xlsx")

# See names of sheets in the excel file.
excel_sheets("data/file.xlsx")


# Read the sheet named "hyphal_width".
hwt <- read_excel("data/file.xlsx", sheet = "hyphal_width")

# Read the first sheet in the file.
hwt <- read_xlsx("data/file.xlsx", sheet = 1)

hwt

class(hwt)

# Demonstrate consequences on not having a header in the
# excel spreadsheet.

nht <- read_excel("data/noheader.xlsx", sheet = "hyphal_width")
nht

nht <- read_excel("data/noheader.xlsx", sheet = "hyphal_width",
                 col_names = FALSE);
nht # Notice the column names.

# You can use colnames() to change the column names to more
# meaningful ones.

colnames(nht) <- c("hw", "strain", "microscope");

#########################
### Reading CSV files ###
#########################

# In the past I have had problems reading MS Excel files in
# R. In such cases you can write a csv from from Excel and
# then use read.csv() or read_csv() to read the csv file
# into a data frame or tibble.

# In MS Excel, after making sure that there is a line of
# header at the top, export your worksheet of interest
# as a csv file.

# Use read.csv() in R to read the csv file into a data
# frame. read.csv() assumes the presence of a header.

# Use read_csv() in R to read the csv file into a tibble.
# read_csv() assumes the presence of a header.

# Some other named arguments to read.csv() which might
# be useful are stringsAsFactors and row.names. We will
# come to these later.

# data frame
hwf <- read.csv("data/hw.csv");
head(hwf);
nrow(hwf);

# tibble
hwt <- read_csv("data/hw.csv");
hwt;

# If your csv file or spreadsheet has no column names make
# sure you use "col_names = FALSE" or "header = FALSE".

nhf <- read.csv("data/noheader.csv");
head(nhf);
nrow(nhf);

nht <- read_csv("data/noheader.csv", col_names = FALSE);
nht;

colnames(nht) <- c("hw", "strain", "microscope");

#################################
### Do the following yourself ###
#################################

# The sheet named "expression" in the file file.xlsx
# contains the columns named "gene", "control" and
# "treatment".

# 1. Read this sheet into a tibble named "expt";

# 2. Try the commands colnames(), nrow(), ncol()
# with expt as the only argument.

# 3. Try the command class(expt).



~~~

<!-- >>> -->

<!-- <<< factors.r -->
# Factors

~~~ {.r}

###############
### Factors ###
###############


# Think of them as categories.


numbers <- c(1.200261, 1.479827, 1.482392, 1.716793, 1.518791, 1.000030,
             1.933209, 1.841415, 1.315890, 1.849663);

category <- c("A", "A", "B", "B", "B", "A",
              "C", "C", "A", "B");

factr <- factor(category);

typeof(category);

class(category);

typeof(factr);

class(factr);

#####################################################
  
tapply(numbers, factr, mean);

tapply(numbers, category, mean);

levels(factr);

levels(category)


# Allow you to assign individual elements of a vector to
# groups or categories.

# Using factors you can do calculations such as mean(),
# sd() etc. on a group-by-group basis.

#################################
### Do the following yourself ###
#################################

# 1. Find out the class of hwf$strain.

# 2. Determine the number of distinct strains in hwf.

# 3. Use tapply() to calculate the mean of the hyphal
# widths for each of the strains in hwf.

# 4. Use tapply() to calculate the standard deviation of
# the hyphal widths for each of the strains in hwf.


### Tidyverse style ###

# Tidyverse functions read_excel and read_csv etc. do not
# make factors when they read data in. If you need to, you
# can convert the columns you want to factors as below.
# read_csv() has an option (col_types) to let you specify
# column types when you are reading data in.

mutate(hwt, strain = factor(strain),
       microscope = factor(microscope))

# group_by() and summarise() by piping
group_by(hwt, strain) %>%
summarise(grmean = mean(hw));

group_by(hwt, microscope) %>%
summarise(grmean = mean(hw));

# Via an intermediate object.
bygr <- group_by(hwt, strain);
summarise(bygr, grmean = mean(hw));


~~~

<!-- >>> -->

<!-- <<< models.r -->
# Statistical Models

~~~ {.r}

###############################
### Statistical Models in R ###
###############################

# Normal
# rnorm(): random generation function
# pnorm(): probability density function
# dnorm(): density function
# qnorm(): quantile function
# Similarly for uniform, binomial, gamma, t, etc.

### Example: Normal distribution.

x <- rnorm(500, mean = 15, sd = 5);

pnorm(35, mean = 15, sd = 5);  # P(X <= 35).

pnorm(35, mean = 15, sd = 5, lower.tail = FALSE);  # P(X > 35)

dnorm(seq(1,10), mean = 5, sd = 1)


###### Begin ignore ######

# Below is to generate plots for explaining. We will see
# some of these functions later we cover plotting. For now,
# simply ignore.

par(mfrow = c(2,1), bg = "lightgrey");
bar <- 15; xsd <- 3;
ddn <- function(x) {
  return(dnorm(x, mean = bar, sd = xsd))
}
cfrom <- bar - 4 * xsd;
cto <- bar + 4 * xsd;
curve(ddn, from = cfrom , to = cto, add = FALSE, lwd = 3,
      main = paste("Mean:", bar,  "and SD:", xsd),
      ylab = "Density",
      xlab = "x"); grid(col = "black");
ppn <- function(x) {
  return(pnorm(x, mean = bar, sd = xsd, lower.tail = FALSE))
}
ppnlt <- function(x) {
  return(pnorm(x, mean = bar, sd = xsd, lower.tail = TRUE))
}
curve(ppn, from = cfrom , to = cto, add = FALSE, lwd = 3,
      col = "darkred", ylab = "Probability", xlab = "x")
curve(ppnlt, from = cfrom , to = cto, add = TRUE, lwd = 3,
      col = "darkgreen")
grid(col = "black")
legend(cfrom, 0.9, "lower.tail FALSE", fill = "darkred", bty = "n")
legend(cfrom, 0.8, "P(X > x)", bty = "n")
legend(cto, 0.9, "lower.tail TRUE", fill = "darkgreen", bty = "n", xjust = 1)
legend(cto, 0.8, "P(X <= x)", bty = "n", xjust = 1)

###### End ignore ######

#################################
### Do the following yourself ###
#################################

# 1. Generate 2000 normally distributed random numbers
# from a population where the mean is 18 and standard
# deviation is 3 and keep them in a vector named x.

# 2. Confirm that x has 2000 elements without listing
# them all.

# 3. Find the minimum, maximum, mean and median of x.

# 4. The function for standard deviation is sd. Use this
# function to find the standard deviation of x. Should
# be approximately 3.

# 5. The function for square root is sqrt. Standard
# error of the mean (SEM) is calculated as the SD
# divided by the square root of the sample size. Without
# using any object other than x, find the SEM. Should
# be approximately 0.06.



~~~

<!-- >>> -->

<!-- <<< logicals.r -->
# Truth in R

~~~ {.r}

##################
### Truth in R ###
##################

x <- rnorm(10, mean = 5, sd = 1);
x;

x <= 5;
x > 5;

g <- x <= 5;
g

# There is a data type called "logical".

typeof(g);

x[g];

#############################################################

### any() and all() ###


any(g);
all(g);

fl <- as.integer(c(1, 2, 3, 4, 5, -6));
any(fl);
all(fl);

fl <- as.integer(c(1, 2, 3, 4, 5, 0));
any(fl);
all(fl);

fl <- as.integer(c(-1, 0, 0, 0));
any(fl);
all(fl);

### Strings cannot be used as logicals.

s <- c("str1", "", "str3");

as.logical(s);         # NA

####################################
### Testing for numeric equality ###
####################################

x <- seq(1,10, by = 0.25);
y <- x;

x == y;
identical(x,y);
all.equal(x,y);


# Divide x by 39.473 then multiply the result by 39.473
z <- x / 39.473
x <- z * 39.473

### identical() and all.equal()
x == y;
identical(x,y);
all.equal(x,y);

# all.equal() never returns FALSE. It either returns TRUE or a string.

x <- rnorm(10, mean = 5, sd = 0.0001);
y <- rnorm(10, mean = 5.0001, sd = 0.0001);

# Below returns a string or TRUE
all.equal(x,y);
all.equal(x,y, tolerance = 1e-4);


# Below return FALSE or TRUE
isTRUE(all.equal(x,y));
isTRUE(all.equal(x,y, tolerance = 1e-4))


# It might be best to decide how much difference you are
# willing to ignore and then do something like below.

any(abs(x-y) > 1e-3)
any(abs(x-y) > 1e-4)
x
y
abs(x - y)

#################################
### Do the following yourself ###
#################################

# 0. rm(m,n);

# 1. Store the value of 0.5 - 0.3 into m

# 2. Store the value of 0.3 - 0.1 into n

# 3. Use == to test the equality of m and n.

# 4. Use identical() to test the equality of m and n.

# 5. Use all.equal() to test the equality of m and n.

# 6. What is the output of isTRUE(all.equal(m,n))?

# 7. What is the magnitude of the difference between
# m and n?

# Refer to the documentation of all.equal() to find out
# the default value of the tolerance argument.

#############################################################
#############################################################

# Use isTRUE(all.equal()), not ==, when comparing floating
# point numbers.

# mean(abs(x-y))/mean(abs(x))



~~~

<!-- >>> -->

<!-- <<< subsetting.r -->
# Subsetting

~~~ {.r}

##############################
### Subsetting data frames ###
##############################

### Subscripts are used to get subsets of data frames
#   (and other types of objects) in R.

df <- data.frame(x = runif(26, min = 3, max = 12),
                 y = runif(26, min = 30, max = 120));
rownames(df) <- LETTERS;
head(df);

# dataframe[rownum, colnum]

df[3, 2]

# dataframe["rowname", "colname"]

df["C", "y"]

### You get all columns if no columns are specified.
# dataframe[rownum, ]
# dataframe["rowname", ]

df[5,]
df["E",]

### You get all rows if no rows are specified.
# dataframe[, colnum]
# dataframe[, "colname"]

df[, "x"]
df[, 1]

### Rows and columns may be specified as vectors.

rows <- seq(17,24);
coln <- c(1);
df[rows, coln];

colhead <- c("x");
df[rows, colhead];

### The dollar notation can only be used with single column names.
#   This is actually very common.
# dataframe$colname

df$x 


### If you ask for multiple columns you always get a data frame in return.
#   Even if you have asked for a single row.

df[c(1,2,3), c("x", "y")]
df[2, c("x", "y")]

### Single columns are returned as vectors.

#################################
### Do the following yourself ###
#################################

# 1. Make a vector, ronum, containing the numbers 21 to 30
# and 51 to 60 by concatenating (remember c()?) two
# calls to seq().

# 2. Use the vector ronum, made above to get a subset data
# frame from hwf consisting of only row numbers in ronum

# 3. From hwf, extract the strain associated with the
# hyphal width in row number 58 and store it in a vector
# named strain58. Check the class and mode of strain58.


### Tidyverse style ###

slice(hwt, ronum);

# hwt[ronum,] # also works.

hwt[58, "strain"]


~~~

<!-- >>> -->

<!-- <<< nainf.r -->
# NA, NaN, Inf, NULL

~~~ {.r}

# Include na.omit, na.fail etc.

##########################
### NA, NaN, Inf, NULL ###
##########################

x <- seq(from = 2, to = 10, by = 2);
x <- c(x, NA);
x;
mean(x);
sum(x);
sd(x);

####################################################

books <- c("Animal Rights", "A Theory of Justice",
           "The Wealth of Nations", "What Money Can't Buy");

authors <- c("Peter Singer", "John Rawls", NA,
             "Michael Sandel");

sold <- c(0.6, 0.2, NA, 0.7);

ba <- data.frame(books = books, authors = authors, sold = sold);

ba

# What is the type of NA?

#################################
### Do the following yourself ###
#################################

# 1. Call the mean(), sum() and sd() functions as above
# but with the additional named argument "na.rm = TRUE"
# e.g. sum(x, na.rm = TRUE)

# 2. Examine the output of summary(x).

# 3. Examine the output of is.na(x).

 
### Inf and -Inf are reserved words in R.

#################################
### Do the following yourself ###
#################################

# 1. Examine the output of 2/0

# 2. Examine the output of -2/0

# 3. Examine the output of is.finite(2/0)

# 4. Examine the output of is.finite(0/2)

# 5. Store the output of 2/0 in i.

# 6. Test that i is infinite.


### NaN (Not a Number) is a reserved word in R.

0/0
x <- log2(-10)
x

# 2. Examine the output of is.nan(x).

# 3. Examine the output of is.na(x).

### NULL

# Complete absence of any value or object.

######################################################
### Think about the output of the following blocks ###
######################################################

is.null(0)

retval <- cat("Just some text\n");
retval

is.null(0/0)

is.null(authors[3]);




~~~

<!-- >>> -->

<!-- <<< matrices.r -->
# Matrices

~~~ {.r}

################
### Matrices ###
################

# Like data frames, these are two dimensional tables with rows and
# columns.

# All elements in a matrix have to be of the same type.

### Making a matrix by using matrix().

x <- seq(1,50);
x;

m <- matrix(x, nrow = 5, ncol = 10);

rownames(m) <- paste("R", seq(1, 5), sep = "");
colnames(m) <- paste("C", seq(1, 10), sep = "");

rownames(m);
colnames(m);


### Making a matrix by adding the dim attribute to a vector.

m <- seq(51,100);
dim(m) <- c(5,10);
m;

attributes(m);
class(m);

rn <- paste("R", seq(1, 5), sep = "");
cn <- paste("C", seq(1, 10), sep = "");
dimnames(m) <- list(rn, cn);

attributes(m);
str(m);

###########################
### Subsetting Matrices ###
###########################

# matrix[rownum, colnum]

# If you do not specify any rows to get, you get all rows.

# If you do not specify any columns to get, you get all columns.

# Just like subsetting data frames, except that if you ask for a single
# row you get a vector.

m[1,2];

m["R1", "C2"];

row3 <- m[3,]; row3;

class(row3);

col1 <- m[, 1]; col1;

class(col1);

### Using more than one rows and columns for subsetting a matrix

# To extract any arbitrary part of your matrix as a matrix you
# may use vectors of rownums and colnums as below.

m[c(2,3), c(2,4)];

s <- m[c("R2","R3"), c("C2", "C4")];

s;

class(s);

attributes(s);

#################################
### Do the following yourself ###
#################################

# 1. Store the sequence 1 to 28 to a vector d.

# 2. Use matrix() to make a matrix x from d having 7 columns and
#    4 rows such that the first row reads 1,2,3,4,5,6,7 as shown
#    below. Refer to the help of matrix() to find out the named
#    argument you need to use.

#   1    2    3    4    5    6    7
#   8    9   10   11   12   13   14
#  15   16   17   18   19   20   21
#  22   23   24   25   26   27   28

# 3. Set the column names to Mon, Tue, Wed...

# 4. Set the row names to Week1, Week2, Week3...

# 5. Examine the matrix x.

#       Mon Tue Wed Thu Fri Sat Sun
# Week1   1   2   3   4   5   6   7
# Week2   8   9  10  11  12  13  14
# Week3  15  16  17  18  19  20  21
# Week4  22  23  24  25  26  27  28

# 6. Extract the date on the Friday of the third week.
# 7. Store all the dates falling on Wednesday to wed.
# 8. Examine wed.



~~~

<!-- >>> -->

<!-- <<< arrays.r -->
# Arrays

~~~ {.r}

##############
### Arrays ###
##############

# Like matrices but can have more than just two dimensions.
# You can think of a 3 dimensional array as a collection
# of matrices.
# Function named array().

a <- array(rnorm(60), c(3,2,10),
           dimnames=list (
             paste("R", seq(1, 3), sep = ""),
             paste("C", seq(1, 2), sep = ""),
             paste("M", seq(1, 10), sep = "")
           )
);




a["R1","C2","M5"];   # vector with one value

a["R2","C2",];    # vector

a[,,"M6"];     # matrix


# Another way to make the same array.

rm(a);
a <- rnorm(60);
dim(a) <- c(2,3,10);
a;

class(a);

attributes(a);

str(a);

a[1,2,5];   # vector with just one value
a[2,2,];    # vector
a[,,6];     # matrix

class(a[1,2,5]);
class(a[2,2,]);
class(a[,,6]);    # matrix



~~~

<!-- >>> -->

<!-- <<< newcol.r -->
# Calculated column

~~~ {.r}

################################
### Calculating a new column ###
################################

# The file data/expression.csv has the same contents as the
# "expression" sheet of file.xlsx.

data <- read.csv("data/expression.csv", header = TRUE, row.names = 1)
head(data)
nf <- cbind(data, logFC = data$treatment - data$control)
head(nf)


#################################
### Do the following yourself ###
#################################

# 1. Instead of getting a new data frame as above, you
# can assign a new column to your original data frame
# (named "data" in this case). This is done by assigning
# a vector to a new column name. Delete data then read
# the file into the data again. Now add a column named
# logFC to data by directly assigning to data$logFC.


# 2. Confirm that "logFC" column has been added to data.


# 3. The products of all the genes in "data" are written
# in the file "products.csv".


# 4. Read the file "data/products.csv" into a data frame
# named temp. Use the named arguments row.names = 1,
# stringsAsFactors = FALSE.


# 5. Add the product column from temp to data but keep
# in mind that the order of genes in temp is not the
# same as in data. Hint: Use rownames(data) to subset
# temp.

# The data frame "data" is needed in the later
# exercises. Keep it.

###################
### dplyr style ###
###################

tib <- read_csv("data/expression.csv", col_names = TRUE);

tib1 <- mutate(tib, logFC = treatment - control);
tib1

# Or, we could simply overwrite tib.

tib <- mutate(tib, logFC = treatment - control);
tib

# dplyr joins #

temptib <- read_csv("data/products.csv")
temptib
left_join(tib, temptib, by = "gene")

# Unwanted column in temptib

temptib <- mutate(temptib, unwanted = rep("unwanted",
                  nrow(temptib)));
temptib
left_join(tib, temptib, by = "gene")
left_join(tib, select(temptib, gene, product),
                  by = "gene")

### Piping ###
# Useful for discarding intermediate objects.
# Easier to apply on several objects.

select(tib, gene, logFC) %>% 
  left_join(select(temptib, gene, product), by = "gene")

tib <- tib %>% left_join(select(temptib, gene, product),
                         by = "gene")
tib;



~~~

<!-- >>> -->

<!-- <<< sorting.r -->
# Sorting

~~~ {.r}

###########################
### Sorting a data frame ###
############################

# Three related functions: sort(), order(), rank()

### sort() ###

set.seed(56)
x <- as.integer(runif(10, 1, 50) * 7)
x
sx <- sort(x)
sx

### The order of numbers in x itself is unchanged by the sort
### operation. sx is a new vector.

### rank() ###

rx <- rank(x)

rx

### order() ###

ox <- order(x);

ox;

# For the same vector x, compare the output of rank() and order()
# x, ox and rx.


# To sort a data frame by some column, we cannot simply use
# sort().

# 1. First we get the order of rows by calling order() on the
#    column we wish to sort by.

# 2. Then we use the order to sort the entire data frame.

ordvec <- order(data$logFC);
sdata <- data[ordvec, ];
head(sdata);
tail(sdata);

#################################
### Do the following yourself ###
#################################

# Run the following block if you have lost or broken the
# data frame named "data" made earlier.
data <- read.csv("data/expression.csv", header = TRUE, row.names = 1)
data$logFC <- data$treatment - data$control;
temp <- read.csv("data/products.csv", row.names = 1,
                    stringsAsFactors = FALSE);
data$product <- temp[rownames(data), "product"]


# 1. It would be more meaningful to sort data on absolute log
#    fold change and to have the highest change at the top.

# 2. The function for getting absolute values is abs().

# 3. Consult the help for order() to find the named argument you
#    need to use to order in decreasing order.

# 4. Now sort data by decreasing value of absolute logFC.



# If you have ties in the column you are sorting by then you
# might wish to use a second column to break the ties.
 
# Below, I am copying data to data2 and setting some control
# values to 50.

data2 <- data;
data2[c("SCO0500", "SCO0501", "SCO0502", "SCO0503") , "control"] <- 50;

# Suppose we wish to sort this data frame in decreasing order by
# the control value but in increasing order by the absolute
# logFC.

ordvec <- order(-data2$control, abs(data2$logFC))
sdata2 <- data2[ordvec, ]
head(sdata2)

###################
### dplyr style ###
###################

stib <- arrange(tib, logFC);
stib

stib <- arrange(tib, desc(abs(logFC)));
stib



~~~

<!-- >>> -->

<!-- <<< condsubset.r -->
# Conditional subsetting

~~~ {.r}

##############################
### Conditional subsetting ###
##############################

# Run the following block if you have lost or broken the
# data frame named "data" made earlier.
data <- read.csv("data/expression.csv", header = TRUE,
                 row.names = 1)
data$logFC <- data$treatment - data$control;
temp <- read.csv("data/products.csv", row.names = 1,
                    stringsAsFactors = FALSE);
data$product <- temp[rownames(data), "product"]


head(data);
nrow(data)


up1.tf <- data$logFC >= 1;
up1 <- data[up1.tf, ];
nrow(up1);
head(up1)

up1c7.tf <- data$logFC >= 1 & data$control >= 7 ;
up1c7 <- data[up1c7.tf, ];
nrow(up1c7)
head(up1c7)

############################
### which() and subset() ###
############################

which17 <- which(data$logFC >= 1 & data$control >= 7);
class(which17)
which17
up1c7 <- data[which17, ];
nrow(up1c7)
head(up1c7)

################
### subset() ###
################

fc1c7 <- subset(data, logFC >= 1 & control >= 7);
fc1c7;
fc1c7 <- subset(data, data$logFC >= 1 & data$control >= 7);
fc1c7;


#################################
### Do the following yourself ###
#################################

# 1. Make a new data frame containing just the control and
#    treatment columns from data.

# 2. Use rowMeans() on this new data frame to get the average of
#    control and treatment for each row.

# 3. Add a column named "avexp" to data containing the vector
#    obtained in the step above.

# 4. Select rows from data where avexp is more than or equal to 6
#    and abs(logFC) is more than or equal to 4.

# 5. How many rows do you get in the step above?

####################################
### Vectorised or non-vectorised ###
####################################

set.seed(3141593); # So that we have the same random numbers.
x <- runif(10, min = 3, max = 6);
y <- runif(10, min = 12, max = 20);
df <- data.frame(x = x, y = y);
df
subset(df, x > 4 & y > 15);
subset(df, x > 4 && y > 15);

subset(df, x < 4 | y < 15);
subset(df, x < 4 || y < 15);

# What is the output of 
# 1. df$x > 4 & df$y > 15
# 2. df$x > 4 && df$y > 15
# 3. df$x < 4 | df$y < 15
# 4. df$x < 4 || df$y < 15

###################
### dplyr style ###
###################

### group_by() and summarise() ###

group_by(stib, lfc3 = abs(logFC) >= 3) %>%
summarise(count = n());

stib %>% group_by(lfc3 = abs(logFC) >= 3) %>%
          summarise(count = n());

# Take care not to put the pipe at the start of a line.

### filter() ###

# Note dplyr::filter() masks functions of the same
# name in base and stats packages.

# At least two-fold upregulated.

stib %>% filter(logFC >= 1)
filter(stib, logFC >= 1)

stib %>% filter(logFC >= 1) %>% summarise(count = n())

# Reminder to demonstrate the consequence of putting
# the pipe at the start of a line.

stib %>% filter(logFC >= 1 & control >= 7) %>%
  summarise(count = n())

stib %>% filter(logFC >= 1 & control >= 7
   & (grepl("regulator|sigma|transcription", product,
      ignore.case = T))
    ); 

# Below is the same as above.
filter(stib, logFC >= 1 & control >= 7
   & (grepl("regulator|sigma|transcription", product,
      ignore.case = T))
    ); 

# The outputs of the above two are assignable.

hiup2 <- filter(stib, logFC >= 1 & control >= 7
   & (grepl("regulator|sigma|transcription", product,
      ignore.case = T))
    ); 
hiup2



~~~

<!-- >>> -->

<!-- <<< regex.r -->
# Regular Expressions

~~~ {.r}

###########################
### Regular Expressions ###
###########################


# Suppose we want all genes with logFC more than 2.5 AND which
# are not "hypothetical proteins".

notHypoth2.5 <- data[
  abs(data$logFC) >= 2.5 &
    data$product != "hypothetical protein"
  , ];

head(notHypoth2.5);
nrow(notHypoth2.5);


# In reality the situation is usually not so simple.

# There might be more than one space between "hypothetical" and
# "protein"?

# Inconsistent capitalisation, e.g. "hypothetical" starting with
# a capital "H".

# Alternative spellings.

# Most characters in a RE match themselves but some have special
# meanings.

# To use the special ones literally you need to escape them with
# "\". e.g. \$ for the dollar symbol.

### Quantifiers ###

# ?     Zero or once.                              
# *     Zero or more.                              
# +     Once or more.                              
# {n}   Exactly n times.                           
# {m,n} At least m times but no more than n times. 
# {m,}  At least m times.                          


### Metacharacters ###
 
# .   Any character                               
# ^   Beginning of the string.                    
# $   End of the string.                          


### Character classes ###

# []  Group of characters e.g. [A-Z].             
# [^] Not in the group of characters e.g. [^0-9]. 

### Abbreviations ###

# \s   Any white space [ \t\n\r\f]
# \d   Any digit [0-9]
# \w   Any word character [0-9a-zA-Z_]
# There are others.

f <- c(
"favourite colour",
"Favourite Colour",
"favorite color",
"favorite colored dress",
"favourite  coloured dress",
"nice colored dress"
);

# The indexes which match.
grep("colour", f);

# The values which match.
grep("colour", f, value = TRUE);

# One or zero 'u'.
grep("colou?r", f, value = TRUE);

# $ to anchor match at end of string.
grep("colou?r$", f, value = TRUE);

grep("favou?rite colou?r$", f, value = TRUE);

# + after the space mean one or more spaces.
# And, case-insensitive matching.
grep("favou?rite +colou?r", f, value = TRUE, ignore.case = TRUE);

# using \s for space. Otherwise same as above.
grep("favou?rite\\s+colou?r", f, value = TRUE, ignore.case = TRUE);

# logical (TRUE or FALSE) output.
grepl("favou?rite +colou?r", f, ignore.case = TRUE);

# Substitution
spacefixed <- gsub(" {2,}", " ", f);
spellfixed <- gsub("vori", "vouri", spacefixed);


tof.not.hypo <- !grepl(
"hypothetical +protein",
data$product, ignore.case = TRUE);

notHypoth2.5 <- data[abs(data$logFC) >= 2.5
& tof.not.hypo , ]

nrow(notHypoth2.5)

head(notHypoth2.5)

# There is a lot more to regular expressions than we
# have demonstrated above. They are a truly general
# purpose high value skill to acquire. All programming
# languages worth learning / using have them in some
# form or another.

### Some useful regular expression related functions in
### R.

# grep()
# regexpr()
# gregexpr()
# sub()
# gsub()



~~~

<!-- >>> -->

<!-- <<< apply.r -->
# Apply

~~~ {.r}

######################
### apply() et al. ###
######################

### apply() ###
x <- runif(60, min = 20, max = 30)
dim(x) <- c(15, 4);
x

apfun <- function(arg1) {
return(arg1 - mean(arg1));
}
apply(x, 1, apfun);
t(apply(x, 1, apfun));

y <- cbind(x, t(apply(x, 1, apfun)));
y


### lapply() and sapply() ###

strepgenes <- read.csv("data/strepGenes.txt", header = F,
                       stringsAsFactors = F);

canoname <- function(x) {
spl <- strsplit(x, "\\s+", perl = TRUE);
cano <- spl[[1]][!grepl("SVEN_|SVEN15_|vnz_|^-$",
                        spl[[1]], perl = TRUE)];
return(cano);
}

lapply(strepgenes[[1]], canoname);

unlist(lapply(strepgenes[[1]], canoname));

sapply(strepgenes[[1]], canoname);

unname(sapply(strepgenes[[1]], canoname));

temp <- cbind(strepgenes,
cano = unname(sapply(strepgenes[[1]], canoname)));

### tapply() ###

# We saw this when doing Factors.



~~~

<!-- >>> -->

<!-- <<< attach.r -->
# Attach

~~~ {.r}

##############################
### Attaching a data frame ###
##############################

# It is a convenience feature allowing you to refer to df$column
# as just column.


rm(list = ls());

control <- seq(1, 100);

df <- read.csv("data/expression.csv", row.names = "gene");

head(df);

attach(df);    # Note the warning following this command.

head(control);   # This control is not coming from the data frame df.

rm(control);

head(control);

search()

detach(df)

search()



~~~

<!-- >>> -->

<!-- <<< clear.r -->
# Clearing your workspace

~~~ {.r}

################################
### Clearing your work space ###
################################

### Scenario 1 ###

x <- rnorm(20, mean = 20, sd = 3);
y <- rnorm(20, mean = 10, sd = 3);

# The above x and y are two objects sitting around from older work.

# Below you create x and y again for some new analysis.
# But the assignment to y fails because of a syntax error.

x <- rnorm(20, mean = 212, sd = 3);
y <- rnorm(20, mean = 210, , sd = 3);   # Assignment fails.

z <- x - y;  # This is being evaluated using the new x and the old y!

z;


### Scenario 2 ###

rm(list=ls());

x <- rnorm(20, mean = 212, sd = 3);
y <- rnorm(20, mean = 210, , sd = 3);  # Assignment fails.

z <- x - y; # Fails.
z;

# You make the correction needed and run again.

x <- rnorm(20, mean = 210, sd = 3);
y <- rnorm(20, mean = 212, sd = 3);

z <- x - y;

z;



~~~

<!-- >>> -->

<!-- <<< rnaseq.r -->
# RNA-seq

~~~ {.r}

###############
### RNA-Seq ###
###############

# Use of the package edgeR to analyse RNASeq data

# edgeR is a Bioconductor package you can read about it
# here
# https://www.bioconductor.org/packages/release/bioc/html/edgeR.html
# Start with the User Guide. You may have to read the
# Reference Manual to do somethings.

# Using these two it is possible to develop a sequence
# of R commands which will do the analysis you desire on
# the data you have. This is a general pattern in
# Bioconductor packages.

### The experiment ###

# Hfq is a small RNA binding protein originally
# identified as the host factor essential for the
# replication of the bacteriophage Qβ. Two biological
# replicates each of the wild type and a Hfq deletion
# mutant of a Pseudomonas strain STR25 were subjected to
# RNA-Seq analysis.

# Reads from the sequencer (in fastq files) were aligned
# to the Pseudomonas genomic sequence using bowtie2 to
# get SAM files.

# bowtie2 --phred33 -x bwt2ndx/STR25 \
# -I 100 -X 400 -p 12 --no-unal \
# -1 rs1-R1.fastq -2 rs1-R2.fastq -S sam/rs1.sam \
# 2> "bowtie2Reports/rs1.bwt2repo"
 
# Then samtools was used to convert these SAM files to BAM files.
 
# samtools view -b -o bam/rs1.bam --threads 8 sam/rs1.sam
 
# Then, in R, the function featureCounts() of the
# Rsubread package was used to get the counts of reads
# overlapping each of the 6003 genes. These counts have
# been saved in the files
 
# featureCounts/rs1
# featureCounts/rs2
# featureCounts/rs3
# featureCounts/rs4

# It is also possible to do the above using samtools
# functions bedcov and depth. Another possibility is to
# use tools in the Rsubread package to carry out the
# alignment as well as feature counting. Anyhow, we have
# to arrive at counts of reads covering each gene.
 
# This is where we will pick up this analysis now and
# take it to completion ending in a table of log fold
# changes for all the genes in the genome. We will use
# functions provided in BioConductor package edgeR for
# doing this.
 
### Analysis in R ###

rm(list = ls());
library("edgeR");

fcpath = "data/featureCounts";
list.files(path = fcpath);
files <- list.files(path = fcpath);
files


# use of rep() to make groups.
temp <- c("wt", "hfq");
groups <- rep(temp, each = 2);

fg <- data.frame(files = files, group = groups);
fg;

# Both columns of fg are factors. So we can call levels() on them.
levels(fg$files);
levels(fg$group);

# labels for each sample.
lbs <- paste(groups, c("A", "B"), sep = "");
lbs
# Examine the values of fg and lbs at this stage
# to make sure they are as you want them to be.

# Now we are ready to read the data from the files into R.

d <- readDGE(fg, path = fcpath, labels = lbs, header = FALSE);

### Examine d here. What class does d belong to? ###
d
class(d)
names(d);

# Saving individual objects
saveRDS(d, file = "after_readDGE.rds");
# d <- readRDS("after_readDGE.rds");

# Saving a list of objects.
save(files, groups, fg, lbs, d, file = "after_readDGE")
# load("after_readDGE");

# Saving entire workspace.
save.image(file = "after_readDGE.image")
# load("after_readDGE.image");
# In a new session of R, you will need to load the
# required libraries again.

### Add gene length information.

gene.lens <- read.csv(file = "data/str25.genelengths",
                      stringsAsFactors = F);

d$genes <- gene.lens;

all.equal(rownames(d$counts), d$genes$gene)


### Normalisation factors.

d <- calcNormFactors(d);

# Drawing of the design matrix
#            wt        hfq
# wtA         1         0
# wtB         1         0
# hfqA        0         1
# hfqB        0         1

# Make the design matrix. Alternative 1.
dm = matrix(c(1,1,0,0,0,0,1,1), nrow = 4);
rownames(dm) = c("wtA", "wtB", "hfqA", "hfqB");
colnames(dm) = c("wt", "hfq");
dm

# Make the design matrix. Alternative 2.
des <- model.matrix(~0+fg$group)
des
rownames(des) <- lbs;
colnames(des) <- levels(fg$group);
des


names(d);
d <- estimateDisp(d, design = des);
names(d)


#################################
### Do the following yourself ###
#################################

# 1. Examine the output of methods(class = class(d))

# 2. Amongst a lot of other things, the above will show
# you a method named rpkm().

# 3. Have a look at the help for rpkm().

# 4. Now use rpkm() to get the RPKMs and store them in
# wt.hfq.rpkm.

# 5.  Examine the top of wt.hfq.rpkm using head().

# 6. write.csv(wt.hfq.rpkm, file =
# "../wt_hfq_RPKM.csv")


###############################
### Differential Expression ###
###############################

# The object d made above has all the data required to
# determine which genes are differentially expressed
# between the WT and Hfq, to what extent and in which
# direction.

# Use of exactTest() and topTags()

et <- exactTest(d, pair = c("wt", "hfq"));

byFC <- topTags(et, sort.by = "logFC");
byFC;

unsorted_lfc <- topTags(et, sort.by = "none", n = nrow(d));
saveRDS(unsorted_lfc$table, "unsorted_lfc.rds")

byPV <- topTags(et, n = nrow(d), sort.by = "PValue");
write.csv(byPV, file = "../wt_hfq_DE.csv", row.names = F);

# head(byPV); # Fails.

class(byPV);
names(byPV);

class(byPV$table);

head(byPV$table);

de.table <- byPV$table

# Hfq is SS_0520 and it has been deleted in the strain
# we are referring to as hfq. So it should be
# significantly down-regulated.

de.table["SS_0520",];


# Collecting the RPKMs of the most highly changed genes
 
# "wt.hfq.rpkm" is a matrix. "de.table" is a data frame.
# "de.table" is ordered such that the most significantly
# changed genes are at the top.
 
# How will you get the RPKMs of the 30 most
# significantly changed genes?
 
#################################
### Do the following yourself ###
#################################

# 1. From de.table get a vector v, of the first 30 row
# names.

# 2. Create top30.rpkm by subscripting wt.hfq.rpkm to
# get the rows in vector v made above and all columns.

##################
### Why logFC? ###
##################

# Linear scale is asymmetric.

# Two fold up-regulation is 2

# Two fold down-regulation is 0.5

# So the entire down-regulation side is squeezed between
# 0 and 1.

plot(de.table$logFC, pch = 20)
plot(2**(de.table$logFC), pch = 20)

lfc <- sample(de.table$logFC, length(de.table$logFC));
plot(lfc, pch = 20)
plot(2**(lfc), pch = 20)

# Table of linear to log2 values.
x <- c(seq(2,16), seq(20, 200, by = 20));
lx <- log2(x);
lin2log <- data.frame(x = x, log2.x = lx)
lin2log

# Table of log2 to linear values
log2lin <- data.frame(log2.x = seq(0, 10, by = 0.5),
                      x = 2**(seq(0,10, by = 0.5))
                      );
log2lin


# Keep the following in mind when working with logFC and
# linearFC.

# 1. You cannot take the log of a negative number.

# 2. Log of numbers less than 1 is negative.

# 3. Log of 1 is zero.

# 4. Log of zero is -Inf.

# 5. Raise the base to the logFC value to get the linear
# fold change.



~~~

<!-- >>> -->

<!-- <<< devices.r -->
# Plotting devices

~~~ {.r}

########################
### Plotting devices ###
########################

### Plotting Devices ###

# Plotting functions such as plot() draw to a plotting
# device.

# If no plotting device(s) exists then one of the
# default type is created.

# Devices are identified by numbers.

# At any time you have only one active plotting device
# although more may be in existence.

# dev.* functions control plotting devices.

# e.g. dev.new() will give you a new plotting device of
# the default type.

# e.g. dev.set() is called with an argument which is the
# device number you wish to make active.

# You can use windows() to get a plotting device for
# your current display. On Linux this command is x11()
# and on OSX it is quartz().

# dev.off() can be used to close any plotting device. If
# no argument is given, the active plotting device is
# closed.

# graphics.off() is used to close all plotting devices.

# It is important to close devices if they are connected
# to a file such as a png or a jpeg file.

#################################
### Do the following yourself ###
#################################

# 1. x <- seq(1, 1e4, by = 5)

# 2. y <- log10(x)

# 3. Get a new plotting device by using dev.new()

# 4. Get a new plotting device by using windows()

# 5. Examine the output of dev.list()

# 6. plot(x,y)

# 7. Use dev.set() to make the other plotting device
# active.

# 8. Again, plot(x,y)

# 9. Use dev.off() to close the active device.

# 10. Use graphics.off() to close all plotting devices.



~~~

<!-- >>> -->

<!-- <<< graphics.r -->
# Base graphics

~~~ {.r}

#####################
### Base graphics ###
#####################

# Several commands, plot(), boxplot(), points(), lines(),
# abline() and several others.

# plot() is usually the starting function call.

# points(), lines(), abline() etc. can add to the plot
# initially made by plot.

# plot() can behave differently depending upon what
# arguments are passed to it.

graphics.off()
x <- rnorm(1000, mean  = 20, sd = 2);
plot(x);
boxplot(x);

m <- matrix(x, nrow = 200, ncol = 5);
plot(m[,1], m[,2]);

plot(m[,1], m[,3]);

boxplot(m)


~~~

<!-- >>> -->

<!-- <<< par.r -->
# Plotting parameters

~~~ {.r}

#################################
### Plotting parameters par() ###
#################################

# e.g. colour, title, margins, plots per device, font
# sizes, etc.

# The function par() is used to set the plotting
# parameters.

# Calls to par() affect the active plotting device.

# Calls to par() may be made before, in-between and
# along with calls to plot(), points, lines() etc.

# The changed parameters affect all subsequent commands
# acting on the plot.

# help("par") lists a lot of parameters that can be set
# via par().

graphics.off();
par(mfrow = c(1,2), bty = "n");
x <- rnorm(50);
par(pch = c(19));
plot(x);
# Examine the plot here and then proceed.
par(pch = c(20), col = "red");
plot(x, main = "This time in red")


#################################
### Do the following yourself ###
#################################

# 1. Get a new plotting device so that all parameters
#    are at their default values.

# 2. x <- seq(2,100, by = 2).

# 3. Use the pch argument to par to set the plotting
#    character to 19.

# 4. Plot x.

# 5. Now set the plotting character to 23 and the
#    plotting colour to "darkred" (the argument to par()
#    for setting the colour is col).

# 6. points() is used to add points to an existing
#    plot. Use points to add points for x/2 to the plot.

# 7. Use points to add points for x/1.75 to the plot
#    but this time use the argument col = "darkgreen" as
#    well so that they get rendered in dark green.



~~~

<!-- >>> -->

<!-- <<< colours.r -->
# Colours

~~~ {.r}

##########################
### Specifying colours ###
##########################

# Using the 10 digits in base 10 (0 to 9), the highest
# two digit number we can express is 99.

10^2

### Counting in base 16 (hexadecimal or hex) ###

# Using 16 digits in base 16 (0 to 9 then A to F), the
# highest two digit number we can express is 255.

16^2

# 0 decimal is 00 in hex and 255 decimal is FF in hex.

# In binary, we need 8 bits (one byte) to count upto
# 255.

2^8

# So two hex digits can be used to express all the
# values possible using a single byte (0 to 255)

#  Colours are made by mixing red, green and,
#  blue channels.

#  Each channel can be anything from 0 to 255
#  (00 to FF in hex).

# An alpha channel can be added for controlling
# transparency of the colour.

#  Be careful with the spelling of colo(u)r.

# The R function rgb() lets you generate colours without
# the need to think in hex or even on a scale of 0 to
# 255.

# Fully opaque red.
rgb(255,0,0, maxColorValue = 255);
rgb(255,0,0,255, maxColorValue = 255);

# Fully opaque green. maxColorValue defaults to 1.
rgb(0,1,0, maxColorValue = 1);
rgb(0,1,0,1, maxColorValue = 1);

# Find it easier to think in percent?
rgb(0,100,0, maxColorValue = 100);
rgb(0,100,0,100, maxColorValue = 100);

barplot(c(1,2), col = c(
  rgb(100,100,0, maxColorValue = 100),
  rgb(0,100,100,100, maxColorValue = 100)
));  


reds <- seq(0.1, 1, 0.05);
greens <- seq(1, 0.1, -0.05);
blues <- c(0);

colmix <- rgb(reds, greens, blues);
names(colmix) <- paste("CO", seq(1, length(colmix)), sep = "");
colmix
ht <- runif(19, min = 5, max = 20)
barplot(ht, col = colmix)


### Function generator ###

crpfun <- colorRampPalette(c(rgb(0,0,1,1), rgb(1,0,0,1)),
                           alpha = TRUE);

crpfun <- colorRampPalette(c("blue", "red"),
                           alpha = TRUE);


crpfun(19)
ht <- runif(19, min = 5, max = 20)
barplot(ht, col = crpfun(19))



# colorRampPalette(c(rgb(0,0,1,1), rgb(1,0,0,1)), alpha = TRUE)(19)



~~~

<!-- >>> -->

<!-- <<< ggplot.r -->
# ggplot2

~~~ {.r}

###############
### ggplot2 ###
###############

# Part of a collection of packages called tidyverse.

# The "gg" in ggplot2 stands for Grammar of Graphics.

# Much more structured and formal than the collection of
# graphics commands in base R.

### Three essential components of a ggplot ###

# 1. Data. Usually a dataframe or a tibble.

# 2. Aesthetic mappings. Think of them as mapping of
# data to axes and colour and shape of points. Function:
# aes().

# 3. Layers. Actual rendering on the plotting device.
# Functions: geom_*(). The kind of plot you want.
# Multiple layers are allowed and indeed, common.

### Function ggplot() ###

# Usually called with two arguments, data and aesthetic
# mapping.

# Data, of course, is your dataframe and you get an
# aesthetic mapping by calling aes(). Calls to aes() can
# be used as argument to ggplot() or to geom_*(). If
# included in ggplot() the values are inherited by all
# layers.

### Layers ###

# To the result of the call to ggplot() we can add
# layers.

### Theme ###

# Function theme().

# Control the overall appearance of a plot.


rm(list = ls());
df <- read.csv("data/tfaG.csv");
head(df);

ggplot(df) +
  geom_point(aes(hour, wt), colour = "red") +
  geom_point(aes(hour, tfaG), colour = "blue")

p1 <- ggplot(df)
p1 + geom_point(aes(hour, wt), colour = "red")
p1 + geom_point(aes(hour, tfaG), colour = "blue")

p1 + geom_point(aes(hour, wt), colour = "red") + 
  geom_point(aes(hour, tfaG), colour = "blue")

p1 + geom_path(aes(hour, wt), colour = "red") + 
  geom_path(aes(hour, tfaG), colour = "blue")



p2 <- p1 +
geom_path(aes(hour, wt), colour = "darkred", size = 1.2) +
geom_path(aes(hour, tfaG), colour = "blue", size = 1.2) +
geom_point(aes(hour, wt), colour = "darkred", size = 3) +
geom_point(aes(hour, tfaG), colour = "blue", size = 3) +
xlab("Hours of growth") +
ylab("Log2 expression") +
ggtitle("Expression of spoF in wt and tfaG deletion strains") +
theme_bw() +
theme(plot.title = element_text(size = 24, face = "bold",
                                hjust = 0.5)) +
theme(axis.title = element_text(size = 15, face = "bold",
                                vjust = 0.5)) +
theme(axis.text = element_text(size = 12, face = "plain",
                                vjust = 0.5))

graphics.off()
p2


# The above is not the best way to organise your data for
# ggplot2. If you have a dataframe like.
# 
#   hour       wt     tfaG
# 1    1 2.597158 2.830137
# 2    2 2.636187 2.800135
# 3    3 2.837971 2.917881
# 4    4 2.711013 2.703749
# 5    5 2.882070 2.838809
# 6    6 2.724614 2.028673


# Reform it to something like this.

#     hour strain  logexpr
# 1      1     wt 2.597158
# 2      2     wt 2.636187
# 3      3     wt 2.837971
# ...
# 118   58   tfaG 2.246144
# 119   59   tfaG 2.755333
# 120   60   tfaG 2.775126

# Think like this: everything on the x-axis goes in one
# column no matter how many categories those values
# belong to. The categories can be specified in other
# columns. Similarly for the values intended for the
# y-axis.


# Below is one way of getting from df to gdf.
#
# There may be functions in available packages to do
# this.

gdf <- data.frame(hour = rep(df$hour, 2),
      strain = c(rep(colnames(df)[2], 60),
      rep(colnames(df)[3], 60)),
      logexpr = c(df$wt, df$tfaG)
)
head(gdf)
class(gdf$strain)

### Plotting begins ###

ggplot(gdf, aes(hour, logexpr, colour = strain)) +
  geom_point() +
  geom_path() +
  theme_bw()


ggplot(gdf, aes(hour, logexpr, colour = strain)) +
  geom_point(size = 3) +
  geom_path(size = 1.2) +
  theme_bw()


### We will save the plot objects as we build them.
# gdfcols <- c("#129628", "#961254");
p1 <- ggplot(gdf, aes(hour, logexpr, colour = strain,
                      shape = strain))
p2 <- p1 + geom_point(size = 3) +
  geom_path(size = 1.2)
p2
# scale_colour_manual(values = c(gdfcols))


# Axis labels and grid lines.

p2 +
  scale_x_continuous(name = "Hours of growth") +
  scale_y_continuous(name = "Log2 expression")

xbr <- seq(0, 60, by = 5)
p3 <- p2 +
  scale_y_continuous(name = "Log2 expression") +
  scale_x_continuous(name = "Duration of culture",
                     breaks = xbr,
                     labels = paste0(xbr, "h"))

p3


# Legend

# p4 <- p3 + scale_colour_discrete(name = "spoF expression in",
#                            labels = c("tfaG deletion strain",
#                                      "Wild type strain"))

# p4


# Our colours

# p3 + labs(colour = "Strain", shape = "Strain")

gdfcols <- c("#129628", "#961254");
# gdfcols <- c("darkblue", "red");
p5 <- p3 + scale_colour_manual(name = "Strain id",
   values = gdfcols,
   labels = c("tfaG deletion strain",
   "M600")) +
scale_shape_discrete(name = "Strain id",
   labels = c("tfaG deletion strain",
   "M600"))

# p3 + labs(colour = "Strain", shape = "Strain")


p5

# Main title

p6 <- p5 + 
ggtitle("Expression of SpoF in M600 and tfaG deletion strains\nover 60 hours of growth in shaken flask");
p6

# Theme
p7 <- p6 +
  theme_bw() +
  theme(plot.title = element_text(size = 24, face = "bold",
                                  hjust = 0.5)) +
  theme(axis.title = element_text(size = 18, face = "bold",
                                  vjust = 0.5)) +
  theme(axis.text = element_text(size = 15, face = "plain",
                                 vjust = 0.5)) +
  theme(legend.text = element_text(size = 15, face = "bold",
                                   hjust = 0.5)) +
  theme(legend.key.size = unit(15, "mm")) +
  theme(legend.title = element_text(size = 18, face = "bold",
                                    hjust = 0.5)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(panel.grid = element_blank())

p7
ggsave("../expression.pdf", p7)
ggsave("../expression.png", p7)

#################################
### Do the following yourself ###
#################################

# A while back we saved a RDS file named
# unsorted_lfc.rds.

# Use readRDS() to read this file into an object named
# ulfc. Determine its class and head() it.

# Use ggplot() to plot the logFC column on the y-axis
# against serial numbers on the x-axis.

# Convert logFC to linear fold change and plot as above.



~~~

<!-- >>> -->

<!-- <<< sporelens.r -->
# Spore lengths

~~~ {.r}

#####################
### Spore Lengths ###
#####################

# Spore lengths were measured in 5 mutants (L, S, T, U,
# V) and the WT of a species of Streptomyces grown in 2
# or 3 different growth media. Now we wish to know what
# effect do the mutations and the growth media have on
# spore lengths.

# In the directory named spores there are files named
# L1, L2, S1, S2, S3, T1, T2, U1, U2, U3, V1, V2, V3,
# WT1, WT2. The beginning letters of the filename
# represents the strain and the following digit
# represents the growth medium.
#
# All of these files contain just one column containing
# spore lengths.

# Please do the following with me.

list.files();
list.dirs();
list.files(path = "data/spores");
list.files(path = "data/spores", full.names = TRUE);
basename("data/spores/WT1");

alist <- list(); # An empty list;

# for loop below.
for(file in list.files(path = "data/spores", full.names = T)) {
  tf <- read.delim(file, col.names = c("spore.length"), header = F);
  id.lst <- basename(file);
  alist[[id.lst]] <- tf$spore.length;
}

names(alist);
head(alist$L1)

# Use of stack() to get data frame of spore lengths and a factor.
sldf <- stack(alist);
head(sldf);
# stack() used the column names "values" and "ind" by default.

# Change column names.
colnames(sldf) <- c("spore.len", "strain");
head(sldf);
class(sldf$strain);
levels(sldf$strain);

# At this time the "strain" column contains both the strain as
# well as the growth medium information. Below we separate the
# two into two different factors.

# Use of regular expressions below.
str <- sub("\\d+$", "", sldf$strain, perl = TRUE);
medium <- sub("^[A-Z]+", "M", sldf$strain, perl = TRUE);

sldf$strain <- factor(str);
sldf$medium <- factor(medium);
head(sldf);

# Now we have the data in a form that can be used in
# function calls. We are going to use this data for plotting
# but this form of data is also the form you need if you
# wish to do ANOVA.

# Run aov() and anova()
# medium.aov <- aov(spore.len ~ medium, data = sldf);
# anova(medium.aov);
#
# strain.aov <- aov(spore.len ~ strain, data = sldf);
# anova(strain.aov)


# Scatter plot spore.len against medium
ggplot(sldf, aes(x = medium, y = spore.len,
                 colour = medium)) +
  geom_point()

# Scatter plot log2(spore.len) against medium
ggplot(sldf, aes(x = medium, y = log2(spore.len),
                 colour = medium)) +
  geom_point()

# Scatter plot log2(spore.len) against strain
ggplot(sldf, aes(x = strain, y = log2(spore.len),
                 colour = strain)) +
  geom_point()

# Separate boxplots of spore.len against growth medium
# using facet_wrap().
ggplot(sldf, aes(x = medium, y = log2(spore.len),
                 colour = medium, fill = medium)) +
  geom_boxplot() + facet_wrap(~strain)


#################################
### Do the following yourself ###
#################################

# 1. Make a boxplot keeping the x axis as strain but
# changing colour and fill to medium. (Don't use
# facet_wrap()).

# 2. Make another boxplot keeping the x axis as medium
# and changing colour and fill to strain.

# 3. How can you keep both the above plots in view so
# that you can compare them?

# 4. Try geom_violin() instead of geom_boxplot().


crpfun <- colorRampPalette(
  c("red", "yellow", "brown", "violet", "blue", "green"),
    alpha = FALSE, space = "Lab");
crpfun(6)

ggplot(sldf, aes(x = medium, y = log2(spore.len),
                 colour = strain, fill = strain)) +
  geom_boxplot() +
  scale_colour_manual(values = (crpfun(6))) +
  scale_fill_manual(values = (crpfun(6)))


# Do we have enough data? Plotting counts.
# geom_bar() counts.

ggplot(sldf, aes(x = medium)) +
  geom_bar(colour = "darkblue", fill = "darkblue") +
  facet_wrap(~strain)

ggplot(sldf, aes(x = strain)) +
  geom_bar(colour = "brown", fill = "brown") +
  facet_wrap(~medium)


# geom_col() does not count. Use this if you have
# counts.

countvec <- integer();
for(s in levels(sldf$strain)) {
s.count <- nrow(subset(sldf, strain == s))
countvec[s] = s.count;
}
countvec

cdf <- data.frame(count = countvec, strain = names(countvec))
cdf

# Or
cdf <- stack(countvec);
colnames(cdf) <- c("count", "strain");
cdf


ggplot(cdf, aes(x = strain, y = count)) +
  geom_col()



~~~

<!-- >>> -->

<!-- <<< histograms.r -->
# Histograms

~~~ {.r}

##################
### Histograms ###
##################

# The data in septaldist.csv contains the inter-septal
# distances measured in the wild type and a mutant
# strain of Streptomyces. The first column contains the
# distances in μm and the second column contains either
# "wt" or "mut" (factor).

df <- read.csv("data/septaldist.csv");
head(df)
tail(df)
wt <- df[df$strain == "wt", 1]
mut <- df[df$strain == "mut", 1]

# We have the two vectors made above and we wish to plot
# their histograms next to each other so that we can see
# the overlap between them. hist() has an argument named
# add which adds to an existing histogram rather than
# plot a new one. So we decide to use this.

hist(wt, breaks = 20, col = "#00009955",
     main = "Histograms of WT and Mutant");

# Notice the add argument below.
hist(mut, breaks = 20, col = "#00990055", add = TRUE);


# Most of the histogram for the mutant is beyond the
# limit of the x axis. So we decide to extend the limits
# of the x-axis.

hist(wt, breaks = 20, col = "#00009955",
     xlim = c( min(c(wt,mut)), max(c(wt,mut)) ),
     main = "Histograms of WT and Mutant");
hist(mut, breaks = 30, col = "#00990055", add = TRUE);

# Now we are losing the tops of the central bars of the
# mutant histogram. So we need to extend the y-axis as
# well. The problem is that the y-axis range gets
# decided in the call to hist(). It cannot be determined
# by the looking at data.

# We need to find out the height of the tallest bar in
# the histogram and adjust the upper limit of the y-axis
# before the histogram is actually drawn. For this we
# need to save the return values of the calls to hist()
# and also suppress actual plotting when hist() is
# called.

hwt <- hist(wt, breaks = 20, plot = FALSE);
hmut <- hist(mut, breaks = 30, plot = FALSE);

# Examine hwt and hmut here.

xlm <- c( min(c(wt,mut)), max(c(wt,mut)) );
ylm <- c( min(c(hwt$counts, hmut$counts)), max(c(hwt$counts, hmut$counts)) );

plot(hwt, col = "#00009955",
     xlim = xlm, ylim = ylm, xlab = "Septal distance",
     main = "Histograms of WT and Mutant"
);

plot(hmut, col = "#00990055", add = TRUE);

### Putting the final plot in a png file.

png(filename = "../intersept.png", width = 1200, height = 800);
plot(hwt, col = "#00009955",
     xlim = xlm, ylim = ylm, xlab = "Septal distance",
     main = "Histograms of WT and Mutant"
);

plot(hmut, col = "#00990055", add = TRUE);
dev.off()

# There are similar functions for jpeg, tiff, pdf,
# postscript (ps) etc.
#
# Journals often want Encapsulated Postscript (eps) files.
# setEPS() calls ps.options() with reasonable defaults.
# Then you can call postscript() just like you called
# png() above to get an eps file of your plot.


### Using the ggplot2 library ###

ggplot(df, aes(sep.dist, fill = strain)) +
  geom_histogram(alpha = 0.3, bins = 60, colour = "grey80")

ggplot(df, aes(sep.dist, fill = strain)) +
  geom_histogram(alpha = 1, bins = 40, colour = "grey80") +
  facet_wrap(~strain)



h1 <- ggplot(df, aes(sep.dist, fill = strain)) +
  geom_histogram(alpha = 0.3, bins = 60, color = "grey70")
h1


# Main title
mt <- "Histogram of interseptal distances in mutant"
mt <- paste(mt, "and wild type strains");

# Axis labels
h1 <- h1 + labs(x = "Interseptal Distance", y = "Count",
                title = mt);
h1

h2 <- ggplot(df, aes(strain, sep.dist, colour = strain)) +
  geom_boxplot()

h3 <- ggplot(df, aes(strain, sep.dist, colour = strain)) +
  geom_violin()


# Writing ggplot2 objects to a file.

pdffn <- c("../hist123.pdf");
pdf(pdffn)
print(h1)
print(h2)
print(h3)
dev.off()

# ggsave()

pdffn <- c("../hist1.pdf");
ggsave(pdffn, h1)

~~~

<!-- >>> -->

<!-- <<< close.md -->
# Closing comments

* R is not the most convenient environment for all kinds
  of data processing.

* Some things might be easier to do in a general purpose
  programming language such as Perl or Python.

* You can find plenty of help and support for Python and
  Perl around the NRP.

* If you are good at learning on your own, or can find
  someone to learn with, then consider Julia very seriously.

* Like any other skill, the best way to maintain and
  advance your R skills is by using it regularly. Try to get
  together with a friend and use R regularly. You will learn
  a lot faster if you practice in pairs. 

* Ability to read documentation quickly is more important
  than you think. Not just in R.

* R manuals webpage.
  https://cran.r-project.org/manuals.html

* PDF of the Introduction to R book.
  https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf

* Tidyverse (dplyr, tibbles, ggplot and more) documentation.
  https://www.tidyverse.org/

* RNA-Seq data analysis using edgeR.
  https://f1000research.com/articles/5-1438

<!-- >>> -->
