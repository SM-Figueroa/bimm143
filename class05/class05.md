# Class05: Data Vis with ggplot
Sawyer (PID: A16335277)

# Graphics systems in R

There are many graphic systems in R for making plots and figures.

We have already played a little with **“base R”** graphics and the
`plot()` function.

Today we will start learning about a popular graphics package called
`ggplot2()`.

This is an add on package - i.e. we need to install it. I install it
(like I install any package) with the `install.packages()` function.

## Creating Scatter plots

``` r
plot(cars)
```

![](class05_files/figure-commonmark/unnamed-chunk-1-1.png)

Before I can use the functions from a package I have to load up the
package from my “library”. We use the `library(ggplot2)` command to load
it up.

``` r
#install.packages('ggplot2')
library(ggplot2)
ggplot(cars)
```

![](class05_files/figure-commonmark/unnamed-chunk-2-1.png)

Every ggplot is made up of at least 3 things: - data (the numbers etc.
that will go into your plot) - aes (how the columns of data map to the
plot aesthetics) - geoms (how the plot actually looks, points, bars,
lines etc.)

``` r
ggplot(cars) + 
  aes(x=speed, y=dist) +
  geom_point()
```

![](class05_files/figure-commonmark/unnamed-chunk-3-1.png)

For simple plots ggplot is more verbose - it takes more code - than base
R plot.

Add some more layers to our ggplot:

``` r
ggplot(cars) +
  aes(x=speed, y=dist) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  labs(title="Stopping distance of old cars",
       subtitle = "A silly example plot") +
  theme_bw()
```

    `geom_smooth()` using formula = 'y ~ x'

![](class05_files/figure-commonmark/unnamed-chunk-4-1.png)

``` r
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)
```

            Gene Condition1 Condition2      State
    1      A4GNT -3.6808610 -3.4401355 unchanging
    2       AAAS  4.5479580  4.3864126 unchanging
    3      AASDH  3.7190695  3.4787276 unchanging
    4       AATF  5.0784720  5.0151916 unchanging
    5       AATK  0.4711421  0.5598642 unchanging
    6 AB015752.4 -3.6808610 -3.5921390 unchanging

Q. Use the nrow() function to find out how many genes are in this
dataset. What is your answer?

``` r
nrow(genes)
```

    [1] 5196

Q. Use the colnames() function and the ncol() function on the genes data
frame to find out what the column names are (we will need these later)
and how many columns there are. How many columns did you find?

``` r
colnames(genes)
```

    [1] "Gene"       "Condition1" "Condition2" "State"     

``` r
ncol(genes)
```

    [1] 4

Q. Use the table() function on the State column of this data.frame to
find out how many ‘up’ regulated genes there are. What is your answer?

``` r
table(genes$State)
```


          down unchanging         up 
            72       4997        127 

Q. Using your values above and 2 significant figures. What fraction of
total genes is up-regulated in this dataset?

``` r
table(genes$State)/nrow(genes) * 100
```


          down unchanging         up 
      1.385681  96.170131   2.444188 

``` r
ggplot(genes) + 
    aes(x=Condition1, y=Condition2) +
    geom_point()
```

![](class05_files/figure-commonmark/unnamed-chunk-10-1.png)

``` r
p <- ggplot(genes) + 
    aes(x=Condition1, y=Condition2, col=State) +
    geom_point()
p
```

![](class05_files/figure-commonmark/unnamed-chunk-11-1.png)

``` r
p + scale_colour_manual( values=c("blue","gray","red") )
```

![](class05_files/figure-commonmark/unnamed-chunk-12-1.png)

``` r
p + labs(title="Gene Expression Changes Upon Drug Treatment",
       x = "Control (no drug)", y = "Drug Treatment")
```

![](class05_files/figure-commonmark/unnamed-chunk-13-1.png)

``` r
# File location online
url <- "https://raw.githubusercontent.com/jennybc/gapminder/master/inst/extdata/gapminder.tsv"

gapminder <- read.delim(url)
```

``` r
#install.packages('dplyr')
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
gapminder_2007 <- gapminder %>% filter(year==2007)
```

``` r
ggplot(gapminder_2007) +
  aes(x=gdpPercap, y=lifeExp) +
  geom_point()
```

![](class05_files/figure-commonmark/unnamed-chunk-16-1.png)

``` r
ggplot(gapminder_2007) +
  aes(x=gdpPercap, y=lifeExp) +
  geom_point(alpha=0.5)
```

![](class05_files/figure-commonmark/unnamed-chunk-17-1.png)

``` r
ggplot(gapminder_2007) +
  aes(x=gdpPercap, y=lifeExp, color=continent, size=pop) +
  geom_point(alpha=0.5)
```

![](class05_files/figure-commonmark/unnamed-chunk-18-1.png)

``` r
ggplot(gapminder_2007) + 
  aes(x = gdpPercap, y = lifeExp, color = pop) +
  geom_point(alpha=0.8)
```

![](class05_files/figure-commonmark/unnamed-chunk-19-1.png)

``` r
ggplot(gapminder_2007) + 
  aes(x = gdpPercap, y = lifeExp, size = pop) +
  geom_point(alpha=0.5)
```

![](class05_files/figure-commonmark/unnamed-chunk-20-1.png)

``` r
ggplot(gapminder_2007) + 
  geom_point(aes(x = gdpPercap, y = lifeExp,
                 size = pop), alpha=0.5) + 
  scale_size_area(max_size = 10)
```

![](class05_files/figure-commonmark/unnamed-chunk-21-1.png)

-Q. Can you adapt the code you have learned thus far to reproduce our
gapminder scatter plot for the year 1957? What do you notice about this
plot is it easy to compare with the one for 2007?

``` r
gapminder_1957 <- gapminder %>% filter(year==1957)
ggplot(gapminder_1957) + 
  aes(x=gdpPercap, y=lifeExp, color=continent, size=pop) +
  geom_point(alpha=0.7) +
scale_size_area(max_size = 10)
```

![](class05_files/figure-commonmark/unnamed-chunk-22-1.png)

In this plot there is one clear outlier on the far right, slightly
distorting the visual. It is not easy to compare to this data to that of
2007 because the data are not side by side and do not share axes.

Q. Do the same steps above but include 1957 and 2007 in your input
dataset for ggplot(). You should now include the layer facet_wrap(~year)
to produce the following plot:

``` r
gapminder_both_yrs <- gapminder %>% filter(year==1957 | year==2007)
ggplot(gapminder_both_yrs) + 
  aes(x=gdpPercap, y=lifeExp, color=continent, size=pop) +
  geom_point(alpha=0.7) +
scale_size_area(max_size = 10) +
  facet_wrap(~year)
```

![](class05_files/figure-commonmark/unnamed-chunk-23-1.png)
