.. -*- mode: rst; noweb-default-code-mode: R-mode; -*-

===========
A Test File
===========
Friedrich Leisch (adapted by Mark Clements)
-------------------------------------------

We load the `ascii` package and set the output type as ``"pandoc"``.

.. code-block:: r

  > library(ascii)
  > options(asciiType="rest")


A simple example: the integers from 1 to 10 are

.. code-block:: r

  > ascii(1:10)

+------+------+------+------+------+------+------+------+------+-------+
| 1.00 | 2.00 | 3.00 | 4.00 | 5.00 | 6.00 | 7.00 | 8.00 | 9.00 | 10.00 |
+------+------+------+------+------+------+------+------+------+-------+

.. the above is just to ensure that 2 code chunks can follow each other

We can also emulate a simple calculator:

.. code-block:: r

  > 1 + 1
  [1] 2
  > 1 + pi
  [1] 4.141593
  > sin(pi/2)
  [1] 1


Now we look at Gaussian data:

.. code-block:: r

  > library(stats)
  > print(x <- rnorm(20))
   [1] -0.7926126 -1.8747680 -2.0707713 -0.7957416 -1.8916423 -0.7869058
   [7]  0.7646434 -0.7948277 -0.1260127 -0.1624380 -0.4139552  0.7975334
  [13]  2.1173623  0.2271348  1.0278714 -0.9185205  0.1550679  1.0805556
  [19]  0.5575418  1.6741352
  > print(t1 <- t.test(x))
  	One Sample t-test
  
  data:  x
  t = -0.42705, df = 19, p-value = 0.6741
  alternative hypothesis: true mean is not equal to 0
  95 percent confidence interval:
   -0.6569038  0.4342688
  sample estimates:
   mean of x 
  -0.1113175 

Note that we can easily integrate some numbers into standard text: The
third element of vector `x` is -2.07077126284279, the
_p_-value of the test is 0.67415.

Now we look at a summary of the famous `iris` data set, and we
want to see the commands in the code chunks:

.. code-block:: r

  > data(iris)
  > ascii(summary(iris),header=TRUE)

+---+---------------+---------------+---------------+---------------+---------------+
|   |  Sepal.Length |  Sepal.Width  |  Petal.Length |  Petal.Width  |       Species |
+===+===============+===============+===============+===============+===============+
| 1 | Min.   :4.300 | Min.   :2.000 | Min.   :1.000 | Min.   :0.100 | setosa    :50 |
+---+---------------+---------------+---------------+---------------+---------------+
| 2 | 1st Qu.:5.100 | 1st Qu.:2.800 | 1st Qu.:1.600 | 1st Qu.:0.300 | versicolor:50 |
+---+---------------+---------------+---------------+---------------+---------------+
| 3 | Median :5.800 | Median :3.000 | Median :4.350 | Median :1.300 | virginica :50 |
+---+---------------+---------------+---------------+---------------+---------------+
| 4 | Mean   :5.843 | Mean   :3.057 | Mean   :3.758 | Mean   :1.199 |               |
+---+---------------+---------------+---------------+---------------+---------------+
| 5 | 3rd Qu.:6.400 | 3rd Qu.:3.300 | 3rd Qu.:5.100 | 3rd Qu.:1.800 |               |
+---+---------------+---------------+---------------+---------------+---------------+
| 6 | Max.   :7.900 | Max.   :4.400 | Max.   :6.900 | Max.   :2.500 |               |
+---+---------------+---------------+---------------+---------------+---------------+

.. code-block:: r

  > library(graphics)
  > pairs(iris)

.. figure:: ReST-test-1-007.jpg

   Pairs plot of the iris data.
.. Note that the caption on the preceding line needs to be indented

.. code-block:: r

  > boxplot(Sepal.Length~Species, data=iris)

.. figure:: ReST-test-1-008.jpg

   Boxplot of sepal length grouped by species
