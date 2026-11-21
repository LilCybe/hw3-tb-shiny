Download this repository and open R in its folder. Install the packages once:

```r
install.packages(c("shiny", "ggplot2", "uwot", "mclust", "BiocManager"))
BiocManager::install("Biobase", ask = FALSE, update = FALSE)
install.packages("NMF")
```

Then run:

```r
shiny::runApp(".")
```

Keep `app.R` and `TBnanostring.rds` in the same folder. The three dimension
reductions are computed once at startup;
