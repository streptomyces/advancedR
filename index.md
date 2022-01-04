## Introduction to R

[This training material](http://streptomyces.org.uk/customers/training/aar/)

#### Start

As soon as you arrive

-   Log in on one of the PCs using your institute identity.
-   If you cannot log in with your institute identity then log in as the
    student name mentioned on the PC. (e.g. b26stu16).
-   Then start Firefox and go to the URL shown below.

#### http://streptomyces.org.uk/customers/training/aar/

##### Timings on both Mon 21 May and Tue 22 May

```
Course                        0930  to   1100
Coffee                        1100  to   1130
Course                        1130  to   1300
Lunch                         1300  to   1400
Course                        1400  to   1530
Coffee                        1530  to   1600
Course                        1600  to   1730
```

##### Some requests

-   Put your phones away during the course. There is a break every 90
    minutes. You can use your phones then.
-   Discussing with your neighbours is fine (even encouraged) when you
    are doing the exercises. *But please stop all conversation when the
    instructor begins to speak*.
-   Please ask if something is not clear.

<hr />

#### Before we start

* We are...
* This is not a statistics course. The focus is on R syntax and techniques
rather than statistics. Although real data is used sometimes, we also
use a lot of toy data generated as we go along.
* Accept the syntax, don't fight it.
* How will this course work?
1.  I will introduce methods by talking about them.
2.  You will run some commands along with me to see the methods in
    action.
3.  I will explain the commands and the syntax you have just seen in
    action.
4.  Sometimes there will be little things for you to do on your own.

#### RStudio

Four frames, clockwise starting from top left.

-   Script editor
-   Environment and History
-   File, Plots, Help etc.
-   Console


* Comments begin with [\#]{.code}. Everything after a [\#]{.code} is
ignored by R.
* Type commands in the script editor. To run a command from the script
editor place the cursor on that line and press Control-Enter.
* You can also select multiple lines and then press Control-Enter to run
all the selected lines.
* Finally, you can use Control-Shift-Enter to *source* the entire script.
The action happens in the Console frame. i.e. any output from R is shown
in the console frame. You can type commands directly in the console.
* R keeps a history of your commands which you can see in the history tab
of the top right frame of RStudio.
* You can select and run commands from the history tab as well.
* However you run a command, it is like typing it into the console and
pressing Enter.

#### R

<table>
<tr>
<td>Function</td>
<td>Unquoted word followed by parentheses</td>
<td><span class="code">mean()</span></td>
</tr>

<tr>
Variable
Unquoted word
gene.lens
</tr>

<tr>
<td>String literals</td>
<td>Quoted alphanumeric characters</td>
<td>&quot;whiA&quot;</td>
</tr>

<tr>
<td>Numeric literals</td>
<td>Unquoted digits and scientific notation</td>
<td>2345, 1e6, 1e-6</td>
</tr>
</table>


-   Parentheses, [()]{.code}, are required in function calls even if you
    are not passing any arguments to the function being called.
    [ls()]{.code} works but [ls]{.code} does not, at least not in the
    way you want it to.
-   The parentheses makes it easy to identify individual function calls
    in long and complex R statements where function calls are embedded
    within other function calls.
-   Unquoted words which are not reserved words are assumed to be
    variable names.
-   Numbers are not quoted.
-   Almost no restrictions on variable names. Be sensible.

    ``` {.code}
    mean <- c(2,3,4,5);
    mean(mean);
    ```

    Do yourselves a favour, use sensible variable names.

-   Commands can continue over multiple lines. R is smart enough to know
    that the command is not finished at the end of the line and will let
    you continue typing to finish the command. It even changes the
    command prompt to "+" to indicate this to you. (Demonstrate now).
-   If you get inside a complicated command which you cannot finish, try
    Control-C (Escape in Windows) to bail out. This usually happens
    because of unmatched parentheses or quotes.


<hr />


#### Getting some data to play with



-   Log in on one of the PCs using your institute identity.
-   If you cannot log in with your institute identity then log in as the
    student name mentioned on the PC. (e.g. b26stu16).
-   Start RStudio from the icon on the Windows Desktop. If there is no
    icon then search for RStudio in the windows start menu and run it
    from there.
-   Copy the following set of commands into the script editor and run
    them one by one by pressing Control-Enter.

``` {.code}

setwd("u:/");
unlink("Rtrain", recursive = TRUE);
dir.create("Rtrain");


setwd("Rtrain");

getwd();

unlink("*");

list.files();

download.file("http://streptomyces.org.uk/customers/training/aar/data.zip",
"data.zip", "internal");

list.files();

unzip("data.zip");

list.files();
```

The last command above should show you a listing of the data files.


------------------------------------------------------------------------


#### Objects



-   R is all about manipulating objects.

##### Simple variables

``` {.code}
x = 2
y = "John"
h = 1.73
height = 1.68
```

##### Data Structures

``` {.code}
j <- list(name = "John", height = 1.70, weight = 70);
t <- list(name = "Tom", height = 1.65, weight = 63, size = "medium");

j$name
t$height
```

##### Objects

-   It is possible to associate functions with data structures (we will
    not be doing this).
-   Objects are simply data structures with methods (functions)
    associated with them.

``` {.code}
bmi(j);
bmi(t);
```

##### Classes and their instances

-   Explain the difference and the relationship.

##### Some peculiarities of Object Orientation in R

-   The syntax of R commands is not typical of OO languages.
-   Usually, (Python, C++) [object.function()]{.code}.
-   Or (Perl) [object-&gt;function()]{.code}.
-   In R, [function(object)]{.code}.
-   The function behaves differently depending upon what object it is
    asked to work on. This is achieved by re-defining the function in
    the class of the object.
-   R treats all data types as objects, though some are more like *data
    structures* than *objects*. The distinction between data structures
    and objects is rather fuzzy in R.
-   If you query the class of a vector, R will respond "numeric" or
    "character", but if you specifically ask R whether a vector x is an
    object, R will respond "FALSE".


------------------------------------------------------------------------


#### Functions and operators in R.



##### Functions

-   Do something for us.
-   Need *arguments* to work on.
-   May assume default arguments if none are provided.

<!-- -->

-   May return something which you may wish to keep as a new object.
-   Some functions may be primarily used for their side-effects rather
    than their return value.
-   User defined functions.
-   Programming constructions such as loops and conditionals. e.g. [for,
    while, if]{.code} etc.

##### Operators

-   Act on *operands*.
-   +, -, \*, \^, /, %%, etc.
-   Not just in R but in any programming language, be aware of operator
    precedence and associativity.
-   [2 \* 2 + 3]{.code} is not the same as [2 \* (2 + 3)]{.code}. [\# 7
    and 10]{.hashcol}
-   [10 - 6 - 3]{.code} is evaluated as [(10 - 6) - 3]{.code} *not*
    [10 - (6 - 3)]{.code} because the [-]{.code}operator is left
    associative. i.e. operations are grouped from the left.

##### User Defined Functions

``` {.code}
bmi <- function (weight, height) {
bmi = weight/(height^2);
return(bmi)
}

bmi(62, 1.62);

```


------------------------------------------------------------------------


#### Vectors I



-   Ordered collection of values.
-   All values have to be of the same type.
-   They can hold numeric or character values.
-   If vector x contains 10, 20 and 30 then\
    \
    x\[1\] refers to 10,\
    x\[2\] refers to 20 and,\
    x\[3\] refers to 30\
    \
    and so on.
-   [x &lt;- c(10, 20, 30)]{.code}


------------------------------------------------------------------------


#### Arguments to functions



Arguments are values or objects a function acts on. The process of
getting arguments into functions is called argument passing.

-   Positional arguments.
-   Named arguments.

``` {.code}
### positional arguments
seq(2, 10); # positional arguments

seq(2, 10, 2);

### Positional and named arguments
seq(2, 10, by = 2);

### Named arguments only
seq(from = 2, to = 10, by = 2);

### Positional and named arguments
seq(2, to = 10, by = 2);

### Positional and named arguments
seq(2, by = 2, 10);

# Works, but confusing. Don't do this.
seq(2, from = 4, to = 20);
```


------------------------------------------------------------------------


#### Return values from functions



If a function returns something and you want to keep it you have to
assign it to a (new) object.

``` {.code}
x <- seq(2, 10);
x;
y = seq(3, 30, by = 3);
y;
```

Some functions do not return anything (NULL).

``` {.code}
x <- seq(10, 100, by = 5);
cat(x, "\n");
catret <- cat(x, "\n");
catret;
```

##### A note about assignment operators in R

``` {.code}
x <- seq(10, 100, by = 5);
```

[\#\#\# Works, but for the wrong reason.]{.hashcol}

``` {.code}
x <- seq(10, 100, by <- 5);

x <- seq(10, by = 5, 100);
```

[\#\#\# Fails. Wrong sign in by argument.]{.hashcol}

``` {.code}
x <- seq(10, by <- 5, 100);
```

[\#\#\# Works. But unintended result.]{.hashcol}

``` {.code}
x <- seq(10, by <- 5, -1);
```

-   Use [&lt;-]{.code} or [=]{.code} on the command line.
-   *Always* use [=]{.code} in function calls and function definitions.

##### Function calls as arguments

``` {.code}
rep(4, times = 4);

rep(seq(1,4), times = 4);

rep(seq(1,4), each = 4);
```


##### Do the following yourself.

1.  Generate the sequence 100 to 10 decreasing by 10 and keep it in a
    vector named [ds]{.code}.
2.  Confirm by checking the contents of [ds]{.code}.

------------------------------------------------------------------------

Write a function named [fi2m]{.code} which requires two arguments,
[feet]{.code} and [inches]{.code}, and returns their conversion to
metres.
Hints:
-   Scroll up to slide *Functions and operators in R* and see the
    function [bmi()]{.code} to remind yourself how a function is
    written.
-   Convert feet to inches and add to the inches supplied as argument.
-   Centimetres are inches times 2.54.
-   Divide centimetres by 100 to get metres.
-   Finally, [return(metres)]{.code}.

Call [fi2m]{.code} with the positional arguments 5 and 6.
Call [fi2m]{.code} with the named arguments [inches = 6, feet =
5]{.code}.


##### Scope in R

-   Use [fi2m()]{.code} to discuss scope here.


------------------------------------------------------------------------


#### Some Common Functions



[\#\#\# c(): concatenate.]{.hashcol}

``` {.code}
x <- c(2.17, 3.14, 13, 29, 48.356);
y <- c(200,300);
z <- c(x, y);
z;
```

[\#\#\# rep(): repeat.\
]{.hashcol} [\# You have already seen this above.]{.hashcol}

``` {.code}
```

[\#\#\# seq(): Generate a sequence of numbers.\
]{.hashcol} [\# You have already seen this above.]{.hashcol}

[\#\#\# ls(): list objects in existence.]{.hashcol}

``` {.code}
ls();
```

[\#\#\# objects(): list objects in existence.]{.hashcol}

``` {.code}
objects();
```

[\#\#\# rm(): remove the named objects.]{.hashcol}

``` {.code}
rm(z);
# rm(list = ls());  # Remove all objects in existence.
```

[\#\#\# length(): Return the length of the named object. ]{.hashcol}

``` {.code}

length(x);

length(y);
```

[\#\#\# min() and max(): Minimum / Maximum element of an
object.]{.hashcol}

``` {.code}
min(x);

min(y);
```

[\#\#\# mean() and median():]{.hashcol}

``` {.code}

mean(x);

median(y);
```

[\#\#\# summary(): Some key statistics about an object.]{.hashcol}

``` {.code}
summary(y);
```

[\#\#\# paste(): Join variables to a single string.]{.hashcol}

``` {.code}
m = 3.976;
s = "Here is a number: ";
paste(s, m, sep = " ");
```


------------------------------------------------------------------------


#### Querying the nature of objects



##### Mode

-   The mode of an object signifies the type of data in it.
-   numeric, character, logical, complex and raw.
-   x &lt;- seq(1, 5, by = 0.5);
-   mode(x);

##### Type and Storage mode

-   These two are approximately the same.
-   storage.mode(x);
-   typeof(x);

##### Class

-   The class of an object is the kind of data structure it is.
-   For vectors it is the same as the mode.
-   list, data.frame, matrix, function.
-   Class of an object determines how functions treat it.
-   class(x);

##### *str()*

-   *str()* is useful to display the internal structure of R objects.
    Especially useful for data frames and more complex objects.

##### *methods()*

-   *methods()*. Lists the methods associated with a class. Most classes
    provide methods to access the data inside objects of their type.
    These methods are the safest way of accessing the data embedded
    inside objects.
-   class(x)
-   methods(class = "data.frame")
-   methods(class = class(x))


------------------------------------------------------------------------


#### Reading documentation



``` {.code}
help("seq")
help.search("sequence")
example("hist")
data()
data("cars")
```

-   Documentation related to packages is available only if the package
    is loaded or is explicitly specified. We will see this later when we
    come to the use of packages.


------------------------------------------------------------------------


#### Vectors



-   Ordered collection of values.
-   All values have to be of the same type.
-   They can hold numeric or character values.

``` {.code}
x <- c(2.9, 4.1, 3.9, 4.5, 3.7, 45.3, 21.6);
x;
x[1];
x[7];
x[4:6]       # ":" is the range operator. Step by 1 only.

y <- 5.4; # is the same as y <- c(5.4);
y;
y[1];

z <- c("ftsZ", "sigE", "bldN", "whiA", "whiB", "rdlA", "chpA");
z;
z[3];
```

-   Elements of a vector can be named.

``` {.code}
names(x) <- z;
x;

x["whiA"]; # is the same as x[4];

v <- c("whiA", "rdlA");

x[v];
```


##### Do the following yourself.

1.  Objects x, y and z should still be in existence at this stage. Do a
    listing of objects to confirm this.
2.  Find the type, mode and class of all the three objects. You will
    have to do this one object at a time.
3.  Have a quick look at the help of [names()]{.code} by issuing
    [help("names")]{.code} at the R prompt.



------------------------------------------------------------------------


#### Recycling in Vector operations




##### Do the following yourself.

1.  Store the sequence from 1 to 5 in vector [x]{.code}.
2.  Store the sequence from 21 to 25 in vector [y]{.code}.
3.  Examine the result of [y \* x]{.code}.
4.  Store the sequence from 21 to 30 in vector [z]{.code}.
5.  Confirm the lengths of x and z using the function [length()]{.code}.
6.  Examine the result of [z \* x]{.code}.
7.  Store the sequence from 1 to 7 in vector [u]{.code}.
8.  Examine the result of [z \* u]{.code}.


-   In operations involving two vectors of unequal length, elements of
    the shorter vector get recycled.
-   If the longer vector is not an integer multiple of the shorter
    vector you get a warning but the operation is valid and successful.


------------------------------------------------------------------------


#### Statistical Models in R



Normal
-   rnorm(): random generation function
-   pnorm(): probability density function
-   dnorm(): density function
-   qnorm(): quantile function

Similarly uniform, binomial, gamma, t, etc.
##### Example: Normal distribution.


![](normal.jpg)


``` {.code}
x <- rnorm(500, mean = 15, sd = 5);

pnorm(35, mean = 15, sd = 5);  # P(X <= 35).

pnorm(35, mean = 15, sd = 5, lower.tail = FALSE);  # P(X > 35)

probs <- c(0.95, 0.05);

qnorm(probs, mean = mean(x), sd = sd(x));

qnorm(probs, mean = mean(x), sd = sd(x), lower.tail = FALSE);

# Below does not depend on a distribution.
quantile(x, probs);
```


##### Do the following yourself.

1.  Generate 2000 normally distributed numbers from a population where
    the mean is 18 and standard deviation is 3 and keep them in a vector
    named *x*.
2.  Confirm that *x* has 2000 elements without listing them all.
3.  Find the minimum, maximum, mean and median of *x*.
4.  Find the 85th percentile value in *x*.
5.  The function for standard deviation is [sd]{.code}. Use this
    function to find the standard deviation of *x*. Should be
    approximately 3.
6.  The function for square root is [sqrt]{.code}. Standard error of the
    mean (SEM) is calculated as the SD divided by the square root of the
    sample size. Without using any variable other than *x*, find the
    SEM. Should be approximately 0.06.



------------------------------------------------------------------------


#### Truth in R



``` {.code}
x <- rnorm(10, mean = 5, sd = 1);

x <= 5;

x > 5;

g <- x <= 5;

g

typeof(g);

x[g];
```

[any()]{.code} and [all()]{.code}

``` {.code}

any(g);
all(g);

fl <- c(1, 2, 3, 4, 5, -6);
any(fl);
all(fl);

fl <- c(1, 2, 3, 4, 5, 0);
any(fl);
all(fl);

fl <- as.integer(c(-1, 0, 0, 0));
any(fl);
all(fl);
```

Strings cannot be used as logicals.

``` {.code}
s <- c("str1", "", "str3");

as.logical(s);         # NA
```

##### Testing for numeric equality

``` {.code}
x <- seq(1,10, by = 0.25);
y <- x;

x == y;
identical(x,y);
all.equal(x,y);

z <- x / 2.4e6
x <- z * 2.4e6

x == y;
identical(x,y);
all.equal(x,y);
```

-   [all.equal()]{.code} never returns [FALSE]{.code}. It either returns
    [TRUE]{.code} or a string.

``` {.code}
x <- c(3,3,3);
y <- c(4,4,4);

# Below returns a string.
all.equal(x,y);

# Below returns FALSE.
isTRUE(all.equal(x,y));
```


##### Do the following yourself.

1.  Store the value of [0.5 - 0.3]{.code} into [m]{.code}
2.  Store the value of [0.3 - 0.1]{.code} into [n]{.code}
3.  Use [==]{.code} to test the equality of [m]{.code} and [n]{.code}.
4.  Use [identical()]{.code} to test the equality of [m]{.code} and
    [n]{.code}.
5.  Use [all.equal()]{.code} to test the equality of [m]{.code} and
    [n]{.code}.
6.  What is the output of [isTRUE(all.equal(m,n))]{.code}?


Always use [isTRUE(all.equal())]{.code}, not [==]{.code}, when comparing
floating point numbers.


------------------------------------------------------------------------


#### Data Frames



-   Like two dimensional tables where columns are vectors.
-   Columns can be of different types.
-   All values in any one column have to be of the same type.
-   All columns have to be of the same length.
-   Commonly used for getting data into R from text files.


------------------------------------------------------------------------


#### Reading CSV files



-   There is a package for reading MS Excel file but I have never been
    able to get it to work reliably.
-   In MS Excel, after making sure that there is a line of header at the
    top, export your worksheet of interest as a csv file.
-   Strings should be quoted and columns separated by commas.
-   Use [read.csv()]{.code} in R to read the csv file into a data frame.
-   [read.csv()]{.code} assumes the presence of a header.
-   Some other named arguments to [read.csv()]{.code} which might be
    useful are [stringsAsFactors]{.code} and [row.names]{.code}. We will
    come to these later.


##### Do the following yourself.

1.  Open the file *hw.xlsx* in MS Excel.
2.  Save it as a csv file named *hw.csv*.
3.  In R, use [read.csv()]{.code} to read the contents of the file
    *hw.csv* into a data frame named *hwf*.
4.  Do a listing of objects to confirm that *hwf* exists.
5.  Try the commands\
    [head(), tail(), colnames(), rownames(), nrow(), ncol()]{.code}\
    with *hwf* as the only argument.



------------------------------------------------------------------------


#### Subsetting data frames



-   Subscripts are used to get subsets of data frames (and other types
    of objects) in R.
-   X~row,col~

``` {.code}
dataframe[rownum, colnum]

dataframe["rowname", "colname"]

# You get all columns if no columns are specified.
dataframe[rownum, ]
dataframe["rowname", ]

# You get all rows if no rows are specified.
dataframe[, colnum]
dataframe[, "colname"]

# Rows and columns may be specified as vectors.
rows <- seq(21,30);
coln <- c(1);
dataframe[rows, coln];

colhead <- c("hw");
dataframe[rows, colhead];

# The dollar notation can only be used with single column names.
# This is actually very common.
dataframe$colname
```

-   If you ask for multiple columns you always get a data frame in
    return. Even if you have asked for a single row.
-   Single columns are returned as vectors.


##### Do the following yourself.

1.  Make a vector, *ron*, containing the numbers 21 to 30 and 51 to 60
    by concatenating (remember [c()]{.code}?) two calls to
    [seq()]{.code}.
2.  Use the vector *ron*, made above to get a subset data frame from
    *hwf* consisting of only row numbers in *ron*.
3.  From *hwf*, extract the strain associated with the hyphal width in
    row number 58 and store it in a vector named *strain58*. Check the
    class and mode of *strain58*.



------------------------------------------------------------------------


#### Factors



-   Think of them as categorical variables.

``` {.code}

numbers <- c(1.200261, 1.479827, 1.482392, 1.716793, 1.518791, 1.000030,
           1.933209, 1.841415, 1.315890, 1.849663);

category <- c("A", "A", "B", "B", "B", "A",
          "C", "C", "A", "B");

factr <- factor(category);

typeof(category);

class(category);

typeof(factr);

class(factr);
```

------------------------------------------------------------------------

``` {.code}
tapply(numbers, factr, mean);

tapply(numbers, category, mean);

levels(factr);

levels(category);
```

-   Allow you to assign individual elements of a vector to groups or
    categories.
-   Using factors you can do calculations such as *mean(), sd()* etc. on
    a group-by-group basis.


##### Do the following yourself.

1.  Find out the class of *hwf\$strain*.
2.  Determine the number of distinct strains in *hwf*.
3.  Use [tapply()]{.code} to calculate the mean of the hyphal widths for
    each of the strains in *hwf*.
4.  Use [tapply()]{.code} to calculate the standard deviation of the
    hyphal widths for each of the strains in *hwf*.



------------------------------------------------------------------------


#### NA, NaN, Inf, NULL



``` {.code}
x <- seq(from = 2, to = 10, by = 2);
x <- c(x, NA);
x;
mean(x);
sum(x);
sd(x);
```


##### Do the following yourself.

1.  Run the last three operations above but with the additional named
    argument [na.rm = TRUE]{.code} e.g. [sum(x, na.rm = TRUE)]{.code}
2.  Examine the output of [summary(x)]{.code}.
3.  Examine the output of [is.na(x)]{.code}.


##### [Inf]{.code} and [-Inf]{.code} are reserved words in R.

##### Do the following yourself.


1.  Examine the output of [2/0]{.code}
2.  Examine the output of [-2/0]{.code}
3.  Examine the output of [is.finite(2/0)]{.code}
4.  Examine the output of [is.finite(0/2)]{.code}
5.  Store the output of [2/0]{.code} in [i]{.code}.
6.  Test that [i]{.code} is infinite.


##### [NaN]{.code} (Not a Number) is a reserved word in R.


1.  [x &lt;- c(x, NaN)]{.code}
2.  Examine the output of [is.nan(x)]{.code}.
3.  Examine the output of [is.na(x)]{.code}.


##### You have already seen [NULL]{.code} above.

``` {.code}
retval <- cat("Just some text\n");
retval;
```


------------------------------------------------------------------------


#### Matrices



-   Like data frames, these are two dimensional tables with rows and
    columns.
-   All elements in a matrix have to be of the same type.

Making a matrix by using [matrix()]{.code}.

``` {.code}
x <- seq(1,50);
x;

m <- matrix(x, nrow = 5, ncol = 10);

rownames(m) <- paste("R", seq(1, 5), sep = "");
colnames(m) <- paste("C", seq(1, 10), sep = "");

rownames(m);
colnames(m);
```

Making a matrix by adding the *dim* attribute to a vector.

``` {.code}
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
```

##### Subsetting Matrices

-   matrix\[rownum, colnum\]
-   If you do not specify any rows to get, you get all rows.
-   If you do not specify any columns to get, you get all columns.
-   Just like subsetting data frames, except that if you ask for a
    single row you get a vector.

``` {.code}
m[1,2];

m["R1", "C2"];

row3 <- m[3,]; row3;

class(row3);

col1 <- m[, 1]; col1;

class(col1);
```

##### Using more than one rows and columns for subsetting a matrix

To extract any arbitrary part of your matrix as a matrix you may use
vectors of rownums and colnums as below.

``` {.code}

m[c(2,3), c(2,4)];

s <- m[c("R2","R3"), c("C2", "C4")];

s;

class(s);

attributes(s);
```


##### Do the following yourself.

1.  Store the sequence 1 to 28 to a vector [d]{.code}.
2.  Use [matrix()]{.code} to make a matrix [x]{.code} from [d]{.code}
    having 7 columns and 4 rows such that the first row reads
    1,2,3,4,5,6,7 as shown below. Refer to the help of [matrix()]{.code}
    to find out the named argument you need to use.
3.  Set the column names to Mon, Tue, Wed...
4.  Set the row names to Week1, Week2, Week3...
5.  Examine the matrix [x]{.code}.
6.  Extract the date on the Friday of the third week.
7.  Store all the dates falling on Wednesday to [wed]{.code}.
8.  Examine [wed]{.code}.



------------------------------------------------------------------------


#### Arrays



-   Like matrices but can have more than just two dimensions.
-   You can think of a 3 dimensional array as a collection of matrices.

``` {.code}

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
```

Another way to make the same array.

``` {.code}
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
```


------------------------------------------------------------------------


#### Lists



-   Lists are collections in which you can mix different types of
    values.
-   Lists can contain objects of any complexity, including other lists.
-   This allows you to create objects of arbitrary complexity.
-   Functions which do complicated things often return their results as
    lists.

``` {.code}
x <- seq(1,20, by = 5);
y <- c("one", "two", "three");
i <- list(x, y);

i;
i[1];   # list
i[[1]]; # vector
class(i[1]); class(i[[1]]);

names(i) <- c("numbers", "strings");

i["numbers"];
i[["numbers"]];
names(i);


j <- list(numbers = x, strings = y, serial = seq(1,20));

j[["serial"]][3:10];
j$serial;
j$serial[15:20];
```


##### Do the following yourself.

1.  Make a matrix of numbers 1 to 15 having 3 columns and 5 rows.
2.  Make a vector of numbers 10 to 20.
3.  Make a list with two named elements, m containing the matrix made in
    step 1 above and v containing the vector made in step 2 above.
4.  Display the names of the elements in the list made above.
5.  Extract the third row of the matrix contained in the list.



------------------------------------------------------------------------


#### Headers



##### Noheader. Filename noheader.tsv

  --------- --------- -------
  SCO0001   5.08099   5.217
  SCO0002   5.83      6.23
  SCO0003   6.29      5.13
  --------- --------- -------

[read.table()]{.code} is used to read in data from files into data
frames.

``` {.code}
dfr <- read.table("noheader.tsv");
head(dfr);
class(dfr$V1);
```

-   Columns get default names.
-   Character columns become factors.

------------------------------------------------------------------------

##### Header type 1. Filename header1.tsv

  --------- ----------- -------
  control   treatment   
  SCO0001   5.08099     5.217
  SCO0002   5.83        6.23
  SCO0003   6.29        5.13
  --------- ----------- -------

``` {.code}
dfr1 <- read.table("header1.tsv");
head(dfr1);
dfr1$V1;
rownames(dfr1);
```

-   Assumed presence of header.
-   The first column of the file is assumed to contain row names.
-   No column named "V1".

------------------------------------------------------------------------

##### Header type 2. Filename header2.tsv

  --------- --------- -----------
  gene      control   treatment
  SCO0001   5.08099   5.217
  SCO0002   5.83      6.23
  SCO0003   6.29      5.13
  --------- --------- -----------

``` {.code}
dfr2 <- read.table("header2.tsv");
head(dfr2);
class(dfr2$V1);
class(dfr2$V2);
class(dfr2$V3);
```

-   Absence of header is assumed.
-   Columns types are *coerced* to character because one value is
    character.
-   Remember all values in a column have to be of one type.
-   Character columns become factors.

##### Header type 2 read with option [header = TRUE]{.code}.

  --------- --------- -----------
  gene      control   treatment
  SCO0001   5.08099   5.217
  --------- --------- -----------

``` {.code}
rm(dfr2);
dfr2 <- read.table("header2.tsv", header = TRUE);
head(dfr2);
class(dfr2$gene);
class(dfr2$control);
class(dfr2$treatment);
```

-   Column types are fine.
-   "gene" column has become a factor, which is not useful.

##### Header type 2 read with options [header = TRUE]{.code} and [stringsAsFactors = FALSE]{.code}.

  --------- --------- -----------
  gene      control   treatment
  SCO0001   5.08099   5.217
  --------- --------- -----------

``` {.code}
rm(dfr2);
dfr2 <- read.table("header2.tsv", header = TRUE, stringsAsFactors = FALSE);
head(dfr2);
class(dfr2$gene);
class(dfr2$control);
class(dfr2$treatment);
```

-   We have a column named "gene" and it is not a factor.
-   Other column types are fine.

##### Header type 2 read with options [header = TRUE]{.code} and [row.names = 1]{.code}.

  --------- --------- -----------
  gene      control   treatment
  SCO0001   5.08099   5.217
  --------- --------- -----------

``` {.code}
rm(dfr2);
dfr2 <- read.table("header2.tsv", header = TRUE, row.names = 1);
head(dfr2);
class(dfr2$control);
class(dfr2$treatment);
head(rownames(dfr2));
```

-   No column named gene.
-   We have row names.

There are lots of other important and useful options in *read.table()*.
Please see the help for *read.table()*.


------------------------------------------------------------------------


#### Refresh the web page for the second day.



Refresh your browser window / tab to reload this web page. {#refresh-your-browser-window-tab-to-reload-this-web-page. .red}
----------------------------------------------------------

http://streptomyces.org.uk/customers/training/aar/
--------------------------------------------------


##### Do the following yourself.

1.  Start Rstudio
2.  Find out your working directory.
3.  If required, change it to "u:/Rtrain"
4.  Confirm that you are now in "u:/Rtrain"
5.  List all the files in your working directory.



------------------------------------------------------------------------


#### The for loop



-   [dfr2]{.code} made above has row names as well as column names.
-   We will use a for loop to find the mean of each row.

``` {.code}
m <- c()
n <- c()
for(rn in rownames(dfr2)) {
  rv = as.matrix(dfr2[rn, ]);
  rom = mean(rv);
  m = c(m, rom);
  n = c(n, rn)
}
names(m) <- n;
m
```

-   There is a function named [rowMeans()]{.code} to do what we have
    done above.


##### Do the following yourself.

1.  Use [rowMeans()]{.code} to confirm the results we obtained above.



------------------------------------------------------------------------


#### Calculating a new column



``` {.code}
data <- read.csv("dfops.csv", header = TRUE, row.names = 1);

nf <- cbind(data, logFC = data$treatment - data$control);
head(nf);
```


##### Do the following yourself.

1.  Instead of getting a new data frame as above, you can assign a new
    column to your original data frame (named "data" in your case). This
    is done by assigning a vector to a new column name like below. Do
    it.\
    [data\$logFC &lt;- data\$treatment - data\$control]{.code}.
2.  Confirm that "logFC" column has been added to [data]{.code}.
3.  The products of all the genes in "data" are written in the file
    "products.csv".
4.  Read the file "products.csv" into a data frame named [temp]{.code}.
    Use the named arguments [ row.names = 1, stringsAsFactors =
    FALSE]{.code}.
5.  Add the product column from [temp]{.code} to [data]{.code} but keep
    in mind that the order of genes in [temp]{.code} is not the same as
    in [data.]{.code}\
    Hint: Use rownames(data) to subset temp.



------------------------------------------------------------------------


#### Sorting a data frame



##### Three related functions: *sort(), order(), rank()*

##### [sort()]{.code}

The order of numbers in x itself is unchanged by the sort operation. sx
is a new vector.

``` {.code}
x <- as.integer(runif(10, 1, 50) * 7);
sx <- sort(x);
sx;
```

##### [rank()]{.code}

``` {.code}
rx <- rank(x);

rx
```

##### [order()]{.code}

``` {.code}
ox <- order(x);

ox;
```

For the same vector x, compare the output of [rank()]{.code} and
[order()]{.code}

To sort a data frame by some column, we cannot use sort().

1.  First we get the order of rows by calling [order()]{.code} on the
    column we wish to sort by.
2.  Then we use the order to sort the entire data frame.

``` {.code}
ordvec <- order(data$logFC);
sdata <- data[ordvec, ];
head(sdata);
tail(sdata);
```


##### Do the following yourself.

1.  It would be more meaningful to sort [data]{.code} on absolute log
    fold change and to have the highest change at the top.
2.  The function for getting absolute values is [abs()]{.code}.
3.  Consult the help for [order()]{.code} to find the named argument you
    need to use to order in decreasing order.
4.  Now sort [data]{.code} by decreasing value of absolute logFC.


------------------------------------------------------------------------

-   If you have ties in the column you are sorting by then you might
    wish to use a second column to break the ties.

Below, I am copying [data]{.code} to [data2]{.code} and setting some
control values to 50.

``` {.code}
data2 <- data;
data2[c("SCO0500", "SCO0501", "SCO0502", "SCO0503") , "control"] <- 50; 
```

-   Suppose we wish to sort this data frame in decreasing order by the
    control value but in decreasing order by the absolute logFC.

``` {.code}
ordvec <- order(-data2$control, abs(data2$logFC));
sdata2 <- data2[ordvec, ];
head(sdata2);
```


------------------------------------------------------------------------


#### Conditional subsetting



``` {.code}
data <- read.csv("dfops.csv", header = TRUE, row.names = 1);
data$logFC = data$treatment - data$control;
temp <- read.csv("products.csv", row.names = 1, stringsAsFactors = FALSE);
data$product <- temp[rownames(data), "product"];
head(data);

up1.tf <- data$logFC >= 1;
up1 <- data[up1.tf, ];

up1c7.tf <- data$logFC >= 1 & data$control >= 7 ;
up1c7 <- data[up1c7.tf, ];
```


##### Do the following yourself.

1.  Make a new data frame containing just the control and treatment
    columns from [data]{.code}.
2.  Use [rowMeans()]{.code} on this new data frame to get the average of
    control and treatment for each row.
3.  Add a column named "avexp" to [data]{.code} containing the vector
    obtained in the step above.
4.  Select rows from [data]{.code} where avexp is more than or equal to
    6 and abs(logFC) is more than or equal to 4.
5.  How many rows do you get in the step above.



------------------------------------------------------------------------


#### Regular Expressions



Suppose we want all genes with logFC more than 2.5 AND which are not
"hypothetical proteins".

``` {.code}
# Breaking a command over multiple lines to improve readability. 

notHypoth2.5 <- data[
abs(data$logFC) >= 2.5 &
data$product != "hypothetical protein"
, ];

head(notHypoth2.5);
nrow(notHypoth2.5);
```

-   In reality the situation is usually not so simple.
-   There might be more than one space between "hypothetical" and
    "protein"?
-   Inconsistent capitalisation, e.g. "hypothetical" starting with a
    capital "H".
-   Alternative spellings.

##### Most characters in a RE match themselves but some\
have special meanings.

To use the special ones literally you need to escape them with "\\".
e.g. \\? for the dollar symbol.

##### Quantifiers

  ------- --------------------------------------------
  ?       Zero or once.
  \*      Zero or more.
  +       Once or more.
  {n}     Exactly n times.
  {m,n}   At least m times but no more than n times.
  {m,}    At least m times.
  ------- --------------------------------------------

##### Metacharacters

  -------- ------------------------------------------------
  .        Any character
  \^       Beginning of the string.
  \$       End of the string.
  \[\]     Group of characters e.g. \[A-Z\].
  \[\^\]   Not in the group of characters e.g. \[\^0-9\].
  -------- ------------------------------------------------

``` {.code}
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

# logical (TRUE or FALSE) output.
grepl("favou?rite +colou?r", f, ignore.case = TRUE); 

# Substitution
spacefixed <- gsub(" {2,}", " ", f);
```

------------------------------------------------------------------------

``` {.code}
tof.not.hypo <- !grepl(
"hypothetical +protein",
data$product, ignore.case = TRUE); 

notHypoth2.5 <- data[abs(data$logFC) >= 2.5 
& tof.not.hypo , ];

nrow(notHypoth2.5);

head(notHypoth2.5);
```

There is a lot more to regular expressions than we have demonstrated
above. They are a truly general purpose high value skill to acquire. All
programming languages worth learning / using have them in some form or
another.

##### Some useful *regular expression* related functions in R.

-   grep()
-   regexpr()
-   gregexpr()
-   sub()
-   gsub()


------------------------------------------------------------------------


#### Attaching a data frame



It is a convenience feature allowing you to refer to [df\$column]{.code}
as just [column]{.code}.

``` {.code}
rm(list = ls());

control <- seq(1, 100);

df <- read.csv("dfops.csv");

head(df);

attach(df);    # Note the warning following this command.

head(control);   # This control is not coming from the data frame df.

rm(control);

control;

search();

detach(df);
```


------------------------------------------------------------------------


#### Packages



##### Packages

-   R code written by others for you to use.
-   Packaged to be loaded on demand to minimise resource use.
-   And to avoid name clashes.
-   When you load a package additional functions defined in the package
    become available for you to use. Sometimes existing functions are
    over-ridden.

Installing edgeR

``` {.code}
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
```

Below is how you can get help for a given package. In Rstudio you can
simply search in the "Help" tab.

``` {.code}
help("edgeR");
vignette("edgeR");
help("edgeR", package="edgeR"); # This should work.
```

To "load" a package

``` {.code}
library("edgeR");
```

Installing edgeR

``` {.code}
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
```

Help for functions defined in the loaded package is available as usual.
Once again, in Rstudio you can simply search in the "Help" tab.

``` {.code}
help(estimateCommonDisp); # Should work.
```


------------------------------------------------------------------------


#### RNA-Seq



##### Use of the package *edgeR* to analyse RNASeq data

edgeR is a Bioconductor package you can read about it
[here.](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)\
Start with the *User Guide*. You may have to read the *Reference Manual*
to do somethings.\
Using these two it is possible to develop a sequence of R commands which
will do the analysis you desire on the data you have. This is a general
pattern in Bioconductor packages.

##### The experiment

*Hfq* is a small RNA binding protein originally identified as the host
factor essential for the replication of the bacteriophage Qβ. Two
biological replicates each of the wild type and a *Hfq* deletion mutant
of a *Pseudomonas* strain STR25 were subjected to RNA-Seq analysis.

Reads from the sequencer (in fastq files) were aligned to the
*Pseudomonas* genomic sequence using *bowtie2* to get SAM files.

``` {.code}
bowtie2 --phred33 -x bwt2ndx/STR25 \
-I 100 -X 400 -p 12 --no-unal \
-1 rs1-R1.fastq -2 rs1-R2.fastq -S sam/rs1.sam \
2> "bowtie2Reports/rs1.bwt2repo"
```

Then *samtools* was used to convert these SAM files to BAM files.
``` {.code}
samtools view -b -o bam/rs1.bam --threads 8 sam/rs1.sam
```

Then, in R, the function *featureCounts()* of the *Rsubread* package was
used to get the counts of reads overlapping each of the 6003 genes.
These counts have been saved in the files\

*featureCounts/rs1*\
*featureCounts/rs2*\
*featureCounts/rs3*\
*featureCounts/rs4*\

This is where we will pick up this analysis now and take it to
completion ending in a table of log fold changes for all the genes in
the genome. We will use functions provided in BioConductor package
*edgeR* for doing this.

##### Analysis in R

``` {.code}
rm(list = ls());
library("edgeR");
# help("readDGE");

files <- c(
"featureCounts/rs1",
"featureCounts/rs2",
"featureCounts/rs3",
"featureCounts/rs4"
);
```

[\# use of rep() to make groups.]{.hashcol}

``` {.code}
temp <- c("wt", "hfq");
groups <- rep(temp, each = 2);

fg <- data.frame(files = files, group = groups);
```

[\# Both columns of fg are factors. So we can call levels() on
them.]{.hashcol}

``` {.code}

levels(fg$group);

lbs <- paste(groups, c("A", "B"), sep = "");
```

[\# Examine the contents of fg and lbs at this stage\
]{.hashcol} [\# to make sure they are as you want them to be.]{.hashcol}

``` {.code}

d <- readDGE(fg, path=".", labels = lbs, header = FALSE);
```

[\#\#\# Examine d here. What class does d belong to? \#\#\#]{.hashcol}

``` {.code}

names(d);

d <- calcNormFactors(d);

# Drawing of the design matrix
#            wt        hfq
# wtA         1         0
# wtB         1         0
# hfqA        0         1
# hfqB        0         1

# Make the design matrix
dm = matrix(c(1,1,0,0,0,0,1,1), nrow = 4);
rownames(dm) = c("wtA", "wtB", "hfqA", "hfqB");
colnames(dm) = c("wt", "hfq");

dm

d <- estimateDisp(d, design = dm);

names(d);



```


##### Do the following yourself.

1.  Examine the output of [methods(class = class(d))]{.code}
2.  Amongst a lot of other things, the above will show you a method
    named [rpkm()]{.code}.
3.  Have a quick look at the top of the help for [rpkm()]{.code}.
4.  From the help it is clear that we need to have gene lengths if we
    wish to calculate the RPKM for each gene. These are in the file
    "str25.genelengths"
5.  Use [read.table]{.code} as shown below\

    ``` {.code}
    gene.lens <- read.table(file = "str25.genelengths", 
    sep = "\t", stringsAsFactors = F, header = F, row.names = 1,
    col.names = c("gene", "len"))
    ```

    to read the file into a data frame named [gene.lens]{.code}

6.  Get a vector of gene names in the order in which they are in
    [d]{.code}. Store it in [dgn]{.code}.\
    [dgn &lt;- row.names(as.matrix(d))]{.code}\
7.  Use [dgn]{.code} to subset [gene.lens]{.code} and get a vector of
    gene lengths. Remember to store it in [gl]{.code}. The row names you
    want are in the vector "dgn" and the only column you want is named
    "len".
8.  Now call [rpkm()]{.code}with the arguments [d]{.code} and
    [gene.length = gl]{.code} to get the RPKMs. Remember to store them
    in [wt.hfq.rpkm]{.code}.
9.  Examine the top of [wt.hfq.rpkm]{.code} using [head()]{.code}.
10. [write.csv(wt.hfq.rpkm, file = "wt\_hfq\_RPKM.csv")]{.code}


##### Differential Expression

The object [d]{.code} made above has all the data required to determine
which genes are differentially expressed between the WT and Hfq, to what
extent and in which direction.

[\# Use of exactTest() and topTags()]{.hashcol}

``` {.code}
et <- exactTest(d, pair = c("wt", "hfq"));

byFC <- topTags(et, sort.by = "logFC");
byFC;

byPV <- topTags(et, n = nrow(d), sort.by = "PValue");
write.csv(byPV, file = "wt_hfq_DE.txt");
# write.table(byPV, file = "wt_hfq_DE.tdf", sep = "\t", quote = FALSE);

head(byPV); # Fails.

class(byPV);
names(byPV);

class(byPV$table);

head(byPV$table);

de.table <- byPV$table;
```

##### Collecting the RPKMs of the most highly changed genes

"wt.hfq.rpkm" is a matrix. "de.table" is a data frame. "de.table" is
ordered such that the most significantly changed genes are at the top.

How will you get the RPKMs of the 30 most significantly changed genes?


##### Do the following yourself.

1.  From de.table get a vector v, of the first 30 row names.
2.  Create top30.rpkm by subscripting wt.hfq.rpkm to get the rows in
    vector v made above and all columns.


##### Heatmap

``` {.code}
heatmap(as.matrix(top30.rpkm), col = heat.colors(10));
```


------------------------------------------------------------------------


#### Linear regression



Reload this web page. {#reload-this-web-page. .red}
---------------------

Here we simulate an attempt to determine the ratio of circumference to
the diameter of circles (π). We have a rope which is good for drawing
circles and measuring in straight lines. Measurement of circumferences
will not be accurate.

``` {.code}

rm(list = ls());

rd <- read.table(file = "regrData.tsv", header = TRUE, sep = "\t");

head(rd);

model <- lm(circum ~ dia, data = rd);

plot(circum ~ dia, pch = 20, col = "#0000ddff",
main = expression(pi),
ylim = c(0, max(circum)),
xlim = c(0, max(dia)),
data = rd);

abline(model);

model;
```

Since we know that when the diameter of a circle is zero its
circumference is also zero, we can put this into our model and force the
regression line to go through the origin.

``` {.code}
model <- lm(circum ~ dia - 1, data = rd);

plot(circum ~ dia, pch = 20, col = "#0000ddff", main = expression(pi),
cex.main = 1.5,
ylim = c(0, max(circum)),
xlim = c(0, max(dia)),
data = rd);

abline(model);

model;
```

Once we have a model we can use it to make predictions.

``` {.code}
ndia <- data.frame(dia = seq(20, 200, by = 5));

ndia$circum <- predict(model, ndia);

ndia;
```

It is important to get the modelling formula right. Read
[help("formula")]{.code}.


------------------------------------------------------------------------


#### ANOVA



Spore lengths were measured in 5 mutants (L, S, T, U, V) and the WT of a
species of *Streptomyces* grown in 2 or 3 different growth media. Now we
wish to know what effect do the mutations and the growth media have on
spore lengths.

In the directory named [spores]{.code} there are files named [L1, L2,
S1, S2, S3, T1, T2, U1, U2, U3, V1, V2, V3, WT1, WT2]{.code}. The
beginning letters of the filename represents the strain and the
following digit represents the growth medium. All of these files contain
just one column containing spore lengths.

Please do the following with me.

``` {.code}
list.files();
list.dirs();
list.files(path = "spores");
list.files(path = "spores", full.names = TRUE);
basename("spores/WT1");

alist <- list(); # An empty list;

for(file in list.files(path = "spores", full.names = T)) {
tf <- read.table(file, sep = "\t", col.names = c("spore.length"));
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
```

At this time the "strain" column contains both the strain as well as the
growth medium information. Below we separate the two into two different
factors.

``` {.code}
# Use of regular expressions below.
str <- sub("\\d+$", "", sldf$strain, perl = TRUE);
medium <- sub("^[A-Z]+", "M", sldf$strain, perl = TRUE);

sldf$strain <- factor(str);
sldf$medium <- factor(medium);
head(sldf);
```

Now we have the data in a form that can be used in function calls.

``` {.code}
# Run aov() and plot.
medium.aov <- aov(spore.len ~ medium, data = sldf);
anova(medium.aov);
```

Use [TukeyHSD()]{.code} to get the confidence intervals of the
differences between the means of the levels of a factor (e.g. medium or
strain).

``` {.code}
medium.tsd <- TukeyHSD(medium.aov);
plot(medium.tsd)
```


##### Do the following yourself.

1.  Carry out aov of the dependence of spore length on strain and view
    the result of calling [anova()]{.code}. How does it compare with the
    result we got earlier for the dependence of spore length on medium?
2.  Call [TukeyHSD]{.code} on the result of [aov()]{.code}, store the
    result as [strain.tsd]{.code} and, plot the CIs.


------------------------------------------------------------------------

Plotting the two CI's side by side.

``` {.code}
# strain.tsd will exist only if you have completed the
# task above.
png("anovaCI.png", width = 1920, height = 1080);
par(mfrow = c(1, 2));
plot(strain.tsd);
plot(medium.tsd);
dev.off();
```

Getting the data frame of differences.

``` {.code}
strain.tsd;
class(strain.tsd$strain);
```

From the above we know that [strain.tsd\$strain]{.code} is a matrix. So
below, we convert it to a data frame and sort it by the "diff" column.

``` {.code}
strain.tsdf <- as.data.frame(strain.tsd$strain);
for.ord <- order(abs(strain.tsdf[["diff"]]), decreasing = TRUE);
strain.tsdf[for.ord,];
```

[anova()]{.code} will accept the result of a call to [lm()]{.code} but
[TukeyHSD()]{.code} will not. So use [aov()]{.code} which calls
[lm()]{.code} anyway.

We could have considered both the factors together.

``` {.code}
both.aov <- aov(spore.len ~ strain + medium, data = sldf);
anova(both.aov);
```


------------------------------------------------------------------------


#### Plotting devices and parameters



##### Plotting Devices

-   Plotting functions such as [plot()]{.code} draw to a *plotting
    device*.
-   If no plotting device(s) exists then one of the default type is
    created.
-   Devices are identified by numbers.
-   At any time you have only one *active* plotting device although more
    may be in existence.
-   *dev.\** functions control plotting devices.
-   e.g. [dev.new()]{.code} will give you a new plotting device of the
    default type.
-   e.g. [dev.set()]{.code} is called with an argument which is the
    device number you wish to make active.
-   You can use [windows()]{.code} to get a plotting device for your
    current display. On Linux this command is [x11()]{.code} and on OSX
    it is [quartz()]{.code}.
-   [dev.off()]{.code} can be used to close any plotting device. If no
    argument is given, the active plotting device is closed.
-   [graphics.off()]{.code} is used to close all plotting devices.
-   *It is important to close devices if they are connected to a file
    such as a png or a jpeg file.*


##### Do the following yourself.

1.  [x &lt;- seq(1, 1e4, by = 5)]{.code}
2.  [y &lt;- log10(x)]{.code}
3.  Get a new plotting device by using [dev.new()]{.code}
4.  Get a new plotting device by using [windows()]{.code}
5.  Examine the output of [dev.list()]{.code}
6.  [plot(x,y)]{.code}
7.  Use dev.set() to make the other plotting device active.
8.  Again, [plot(x,y)]{.code}
9.  Use [dev.off()]{.code} to close the active device.
10. Use [graphics.off()]{.code} to close all plotting devices.


##### Plotting parameters

-   e.g. colour, title, margins, plots per device, font sizes, etc.
-   The function [par()]{.code} is used to set the plotting parameters.
-   Calls to [par()]{.code} affect the active plotting device.


##### Do the following yourself.

1.  The plotting parameter [mfrow]{.code} sets the number of plots to
    make on the same device. The arguments to it specify the number of
    rows and the number of columns of plots on the device. e.g. argument
    [c(2,2)]{.code} means 4 plots on the same device.
2.  Use [par(mfrow=c(1,2))]{.code} to set the device to make two plots
    side by side.
3.  Plot x against y.
4.  Use [rnorm()]{.code} to generate 2000 normally distributed numbers
    and store them in a vector named [z]{.code}.
5.  [hist(z)]{.code} to make a histogram of the numbers in [z]{.code}.



------------------------------------------------------------------------


#### Plotting



-   Several commands,[ plot(), boxplot(), points(), lines(),
    abline()]{.code} and several others.
-   [plot()]{.code} is the most basic and most frequently called.
-   Plotting commands can behave differently depending upon what
    arguments they are called to work on.

``` {.code}

par(mfrow = c(2,3));

x <- rnorm(1000, mean  = 20, sd = 2);
plot(x);
boxplot(x);

m <- matrix(x, nrow = 200, ncol = 5);
plot(m[,1], m[,2]);

plot(m[,1], m[,3]);

boxplot(m);
```


##### Do the following yourself.

1.  Use [read.csv()]{.code} to read the file [septaldist.csv]{.code}
    into a data frame [d]{.code}.
2.  Separate the "wt" and "mut" values into two different vectors named
    "wt" and "mut".
3.  Put the two vectors made above in a list named "sepdist".
4.  Switch off all plotting devices.
5.  Use [par(mfrow = ...))]{.code} to set the plotting device to plot 4
    plots in 2 rows and two columns.
6.  [hist(sepdist\$wt)]{.code}
7.  [hist(sepdist\$mut)]{.code}
8.  [boxplot(sepdist)]{.code}



------------------------------------------------------------------------


#### Plotting parameters *par()*



-   In R, you can build up your plot in several steps, each step adding
    a little bit to your plot.
-   You cannot remove something once it is on the plot. So if you get
    something wrong you need to start again.
-   Calls to [par()]{.code} may be made before, in-between and along
    with calls to plot().
-   The changed parameters affect all subsequent commands acting on the
    plot.

``` {.code}
graphics.off();
windows();
par(mfrow = c(1,2));
x <- rnorm(50);
par(pch = c(19));
plot(x);
# Examine the plot here and then proceed.
par(pch = c(23), col = "red");
plot(x, main = "This time in red");
```


##### Do the following yourself.

1.  Get a new plotting device so that all parameters are at their
    default values.
2.  [x &lt;- seq(2,100, by = 2)]{.code}.
3.  Use the [pch]{.code} argument to [par]{.code} to set the plotting
    character to 19.
4.  Plot [x]{.code}.
5.  Now set the plotting character to 23 and the plotting colour to
    "red" (the argument to [par()]{.code} for setting the colour is
    [col]{.code}).
6.  [points()]{.code} is used to add points to an existing plot. Use
    [points]{.code} to add points for [x/2]{.code} to the plot.
7.  Use [points]{.code} to add points for [x/1.75]{.code} to the plot
    but this time use the argument [col = "darkgreen"]{.code} as well so
    that they get rendered in dark green.


##### The basic ideas are

-   Get a plotting device of your choice.
-   Set the plotting parameters how you want them.
-   Make the plot.
-   Add to the plot using functions like [points()]{.code}, changing the
    plotting parameters if you wish to.
-   Close the plotting device (important when the plotting device is a
    file).


------------------------------------------------------------------------


#### Specifying colours more precisely



##### Counting in base 16 (hexadecimal or hex)

-   We need 16 digits!
-   0 to 9 and then A to F
-   0 decimal is 00 in hex and 255 decimal is FF in hex.
-   Try [16\^2]{.code} at the R prompt.

<!-- -->

-   Made by mixing red, green and, blue values.
-   Each value can be anything from 0 to 255.
-   An alpha channel can be added for controlling transparency of the
    colour.
-   Be careful with the spelling of *colo(u)r*.

``` {.code}
# Fully opaque red.
rgb(255,0,0, maxColorValue = 255);
rgb(255,0,0,255, maxColorValue = 255);

# Fully opaque green
rgb(0,1,0, maxColorValue = 1);
rgb(0,1,0,1, maxColorValue = 1);

reds <- seq(0.1, 1, 0.05);
greens <- seq(1, 0.1, -0.05);
blues <- c(0);

colmix <- rgb(reds, greens, blues); # maxColorValue defaults to 1.
names(colmix) <- paste("CO", seq(1, length(colmix)), sep = "");
```


------------------------------------------------------------------------


#### Plotting to a file.



The expression of a gene was measured at 1 hour interval for 60 hours in
the wild type and a strain in which the gene *tfaG* was deleted. This
data is in two columns named "wt" and "tfag". Obviously, there are 60
rows of data in the file. We will plot the expression of this gene
against time for both the wild type and the "tfag" deletion strains.

-   We will make this plot directly to a png file.
-   Using functions such as [png(), jpeg(), tiff()]{.code} or
    [pdf()]{.code} you can set the plotting device to a file of the
    desired type.

[\# Read the data and open a graphics device to a png file.]{.hashcol}

``` {.code}
df <- read.table(file = "tfag.tsv", header = TRUE, sep = "\t");
head(df);
attach(df);

graphics.off();
png(filename = "tfaG.png", width = 1200, height = 800);
```

[\# Set some plotting parameters.]{.hashcol}

``` {.code}
par(cex = 1, lwd = 2);
par(cex.axis = 1.2, cex.lab = 1.4, cex.main = 2, mar = c(5,5,4,2));
```

[\# Store colours in a vector named pcol and give them strain
names.]{.hashcol}

``` {.code}
# pcol <- rainbow(2);
pcol <- c("#cc0000ff", "#00cc00ff");
names(pcol) <- c("wt", "tfag");
```

[\# Main title and axis labels.\
]{.hashcol} [\# Suppress x-axis for later customisation.]{.hashcol}

``` {.code}
plot( smooth(wt), type = "o", col = pcol["wt"],
bg = pcol["wt"], pch = 21, xlab = "Time",
ylab = expression(paste(log[2], " Expression")),
main = expression(paste("Expression of ", italic(impG), " in WT and ",
Delta, italic(tfaG), " backgrounds")), xaxt = "n" );
```

[\# Add points for tfaG.]{.hashcol}

``` {.code}
points(smooth(tfaG), type = "o", col = pcol["tfag"],
bg = pcol["tfag"], pch = 22);
```

[\# Add x-axis.]{.hashcol}

``` {.code}
axis(1, seq(5, 60, by = 5), labels = paste(seq(5, 60, by = 5),
"h", sep = "") );
```

[\# Add a couple of horizontal lines.]{.hashcol}

``` {.code}
th.wt = mean(wt) + 2 * (sd(wt));
abline(h = th.wt, col = pcol["wt"]);
th.tfag = mean(tfaG) + 2 * (sd(tfaG));
abline(h = th.tfag, col = pcol["tfag"]);
```

[\# Add some text labels.]{.hashcol}

``` {.code}
text(5, c(th.wt, th.tfag) - 0.1,
labels = c(paste(c("wt", "tfag"), "mean + 2SD"))
);
```

[\# Add a legend]{.hashcol}

``` {.code}
legend( x = "topright", legend = c("Wild Type",
expression(paste(Delta, italic(tfaG)))),
col = pcol, pt.bg = pcol, cex = 1.2,
pch = c(21, 22) );
```

[\# Close the plotting device.]{.hashcol}

``` {.code}
dev.off();
```

[\# Detach the dataframe.]{.hashcol}

``` {.code}
detach(df);
```


------------------------------------------------------------------------


#### Histograms



The data in [septaldist.csv]{.code} contains the inter-septal distances
measured in the wild type and a mutant strain of *Streptomyces*. The
first column contains the distances in μm and the second column contains
either "wt" or "mut".

``` {.code}
df <- read.csv("septaldist.csv");
head(df)
tail(df)
wt <- df[df$strain == "wt", 1]
mut <- df[df$strain == "mut", 1]
```

We have the two vectors made above and we wish to plot their histograms
next to each other so that we can see the overlap between them.
[hist()]{.code} has an argument named [add]{.code} which adds to an
existing histogram rather than plot a new one. So we decide to use this.

``` {.code}
hist(wt, breaks = 20, col = "#00009955",
main = "Histograms of WT and Mutant");

# Notice the add argument below.
hist(mut, breaks = 20, col = "#00990055", add = TRUE);
```

Most of the histogram for the mutant is beyond the limit of the x axis.
So we decide to extend the limits of the x-axis.

``` {.code}
hist(wt, breaks = 20, col = "#00009955",
xlim = c( min(c(wt,mut)), max(c(wt,mut)) ),
main = "Histograms of WT and Mutant");
hist(mut, breaks = 30, col = "#00990055", add = TRUE);
```

Now we are losing the tops of the central bars of the mutant histogram.
So we need to extend the y-axis as well. The problem is that the y-axis
range gets decided in the call to [hist()]{.code}. It cannot be
determined by the looking at data.

We need to find out the height of the tallest bar in the histogram and
adjust the upper limit of the y-axis before the histogram is actually
drawn. For this we need to save the return values of the calls to
[hist()]{.code} and also suppress actual plotting when [hist()]{.code}
is called.

``` {.code}
hwt <- hist(wt, breaks = 20, plot = FALSE);
hmut <- hist(mut, breaks = 30, plot = FALSE);
```

Examine [hwt]{.code} and [hmut]{.code} here.
``` {.code}
xlm <- c( min(c(wt,mut)), max(c(wt,mut)) );
ylm <- c( min(c(hwt$counts, hmut$counts)), max(c(hwt$counts, hmut$counts)) );

plot(hwt, col = "#00009955",
xlim = xlm, ylim = ylm, xlab = "Septal distance", 
main = "Histograms of WT and Mutant"
);

plot(hmut, col = "#00990055", add = TRUE);
```


------------------------------------------------------------------------


#### Hypothesis testing



###### One sample t-test

There is no data file to be read in. We will generate the data to work
on by calling [rnorm()]{.code}.

``` {.code}

x <- rnorm(50);  # mean defaults to zero and sd defaults to 1;
mean(x);
sd(x);

tt <- t.test(x);

class(tt);
attributes(tt);
names(tt);

tt;
```

###### Two sample t-test

``` {.code}

strainV <- sldf[sldf$strain == "V", "spore.len"]
strainS <- sldf[sldf$strain == "S", "spore.len"]
strainWT <- sldf[sldf$strain == "WT", "spore.len"]

ttv <- t.test(strainV, strainWT);
ttv
names(ttv);
ttv$p.value;

tts <- t.test(strainS, strainWT);
tts
names(tts);
tts$p.value;
```


##### Do the following yourself.

1.  There is a test known as *wilcox.test* which works just like the
    *t.test* used above. It is the non-parametric equivalent of the
    t-test.
2.  Carry out a two sample *wilcox.test* on [strainV]{.code} and
    [strainWT]{.code}.
3.  Carry out a two sample *wilcox.test* on [strainS]{.code} and
    [strainWT]{.code}.



------------------------------------------------------------------------


#### Scripting



Below are the contents of the file script.r. The following command is
issued on the windows command prompt.\
\
[ "C:\\Program Files\\R\\R-3.3.2\\bin\\Rscript" script.r dfops.csv
fromscript.csv ]{.code}

Below are contents of the script file. Do not run them on the R prompt.

``` {.codisp}
args <- commandArgs(trailingOnly = TRUE);
cat(args, "\n");
infile <- args[1];
outfile <- args[2];

cat(infile, outfile, "\n");
idf <- read.csv(file = infile, row.names = 1, header = TRUE);
cat("Read", infile, "\n");

idf$lfc <- idf$treatment - idf$control;
odr <- order(abs(idf$lfc), decreasing = TRUE);
odf <- idf[odr,];
head(odf);
write.csv(odf, file = outfile);
cat("Written", outfile, ".\n");
```


------------------------------------------------------------------------


#### Clearing your work space

[▼▲](Clearing%20your%20work%20space)


``` {.code}
x <- rnorm(20, mean = 20, sd = 3);
y <- rnorm(20, mean = 10, sd = 3);

# The above x and y are two objects sitting around from older work.

# Below you create x and y again for some new analysis.
# But the assignment to y fails because of a syntax error.

x <- rnorm(20, mean = 212, sd = 3);
y <- rnorm(20, mean = 210, , sd = 3);   # Assignment fails.

z <- x - y;  # This is being evaluated using the new x and the old y!

z;
```

------------------------------------------------------------------------

``` {.code}
rm(list=ls());

x <- rnorm(20, mean = 212, sd = 3);
y <- rnorm(20, mean = 210, , sd = 3);  # Assignment fails.

z <- x - y; # Fails.
z;

# You make the correction and run again.

x <- rnorm(20, mean = 210, sd = 3);
y <- rnorm(20, mean = 212, sd = 3);

z <- x - y;

z;
```


------------------------------------------------------------------------


#### Closing comments



-   Ten simple rules for biologists learning to program. [PLOS
    Computational
    Biology](http://dx.plos.org/10.1371/journal.pcbi.1005871)
    ([Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/29300745))
-   Some things might be easier to do in a general purpose programming
    language such as Python or Perl.
-   You can find plenty of help and support for Python and Perl around
    the NRP.
-   Like any other skill, the best way to maintain and advance your R
    skills is by using it regularly.
-   R is not the most convenient environment for all kinds of data
    processing.
-   [R manuals webpage](https://cran.r-project.org/manuals.html).
-   [PDF of the Introduction to R
    book.](https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf)


------------------------------------------------------------------------

