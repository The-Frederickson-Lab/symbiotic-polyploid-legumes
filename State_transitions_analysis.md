State transitions analysis
================
Tia Harrison
2022-01-20

## Overall setup

### Global setup, load relevant packages

### Ploidy dataset and phylogeny

The datasheet contains information about the symbiotic status (Fixer
column) of the legume species where 1 indicates a legume that forms
nodules with rhizobia (aka a mutualist) and where 0 indicates a plant
that forms no association with rhizobia (nonmutualist). The diPloidyLow
column indicates whether plant species are diploids (coded as 0) or
polyploid (coded as 1). The diPloidyLow values were calculated from
genus-level data. The NewPloidy column also indicates diploidy (0) and
polyploidy (1) but was calculated from a combination of genus and
subfamily level data.

``` r
# Dataset 
ploidy<- read.csv("Legume_ploidy_dataset.csv", row.names=1)

# Phylogeny 
legumes1<-read.tree("Vascular_Plants_rooted.dated.tre")

# Clean up the data 
ploidy1<-ploidy%>%
  rownames_to_column() %>%
  mutate(NewPloidy = ifelse(is.na(diPloidyLow), disfPloidy_corrected, diPloidyLow)) %>%
  mutate(Fixer = ifelse(is.na(Fixer) & numOTUs >= 1, 1, Fixer)) %>%
  mutate(Fixer = ifelse(is.na(Fixer) & numGenera >=1 , 1, Fixer)) %>%
  column_to_rownames()
```

### Pruning for the genus + subfamily level ploidy data

We pruned the Zanne et al (2014) angiosperm phylogeny for legume species
in our dataset.

``` r
# Remove all NAs from the data for the relevant columns 
ploidy_sub <- ploidy1 %>% 
  filter(!is.na(Fixer), !is.na(NewPloidy))

# Make sure the data and tree match up 
TreeOnly_sub <- setdiff(legumes1$tip.label, rownames(ploidy_sub))
DataOnly_sub <- setdiff(rownames(ploidy_sub), legumes1$tip.label)

# Prune the tree 
sub_tree <- drop.tip(legumes1, TreeOnly_sub)
sub_tree$node.label <- NULL

# Get the tip names 
tip_names_sub<-data.frame(sub_tree$tip.label)

# Filter the dataset for tree species 
ploidy_sub2<-ploidy_sub %>%
  rownames_to_column("Species") %>%
  filter(Species %in% tip_names_sub$sub_tree.tip.label) %>%
  dplyr::select(Species, Fixer, NewPloidy)
```

## Estimating state transitions on subfamily and genus level ploidy data

### Simple model

We ran a simple model first with no hidden states. Dual transitions are
not allowed in the model, only one state transition at a time. We
allowed polyploids to revert back to diploids in the model.

Results are interpreted as transitions from ROW to COLUMN Traits: 1 =
nonmutualist diploid 2 = nonmutualist polyploid 3 = mutualist diploid 4
= mutualist polyploid

``` r
# Simple model with two traits and no hidden states 
MK_sub_simple<- corHMM(phy=sub_tree, data=ploidy_sub2, rate.cat=1)
```

    ## State distribution in data:
    ## States:  1   2   3   4   
    ## Counts:  44  24  523 248 
    ## Beginning thorough optimization search -- performing 0 random restarts 
    ## Finished. Inferring ancestral states using marginal reconstruction.

``` r
# Look at the results of the model 
MK_sub_simple
```

    ## 
    ## Fit
    ##       -lnL      AIC    AICc Rate.cat ntax
    ##  -439.7148 895.4295 895.603        1  839
    ## 
    ## Rates
    ##            (1,R1)      (2,R1)      (3,R1)      (4,R1)
    ## (1,R1)         NA 0.018954239 0.005072482          NA
    ## (2,R1) 0.03440833          NA          NA 0.004561596
    ## (3,R1) 0.00236954          NA          NA 0.023474456
    ## (4,R1)         NA 0.000000001 0.077509318          NA
    ## 
    ## Arrived at a reliable solution

``` r
plotMKmodel(MK_sub_simple)
```

![](State_transitions_analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Hidden states model

We included rate.cat=2 to specify hidden states in the model. This model
with hidden states performed better than the model with no hidden states
and so it is reported in the manuscript.

Results are interpreted as transitions from ROW to COLUMN Traits: 1 =
nonmutualist diploid 2 = nonmutualist polyploid 3 = mutualist diploid 4
= mutualist polyploid

States: R1 = rate category 1 R2 = rate category 2

``` r
# Two traits and hidden states model 
MK_sub_hidden<- corHMM(phy=sub_tree, data=ploidy_sub2, rate.cat=2)
```

    ## State distribution in data:
    ## States:  1   2   3   4   
    ## Counts:  44  24  523 248 
    ## Beginning thorough optimization search -- performing 0 random restarts 
    ## Finished. Inferring ancestral states using marginal reconstruction.

``` r
# Look at the results of the model 
MK_sub_hidden
```

    ## 
    ## Fit
    ##       -lnL      AIC     AICc Rate.cat ntax
    ##  -388.9723 813.9446 814.7788        2  839
    ## 
    ## Rates
    ##             (1,R1)       (2,R1)      (3,R1)      (4,R1)      (1,R2)      (2,R2)
    ## (1,R1)          NA 1.4291645132 0.000000001          NA 0.006423362          NA
    ## (2,R1) 7.348182633           NA          NA 0.000000001          NA 0.006423362
    ## (3,R1) 0.001580376           NA          NA 0.000000001          NA          NA
    ## (4,R1)          NA 0.0007919703 0.000000001          NA          NA          NA
    ## (1,R2) 0.024014729           NA          NA          NA          NA 0.032431450
    ## (2,R2)          NA 0.0240147293          NA          NA 0.011605394          NA
    ## (3,R2)          NA           NA 0.024014729          NA 0.004257555          NA
    ## (4,R2)          NA           NA          NA 0.024014729          NA 0.000000001
    ##             (3,R2)      (4,R2)
    ## (1,R1)          NA          NA
    ## (2,R1)          NA          NA
    ## (3,R1) 0.006423362          NA
    ## (4,R1)          NA 0.006423362
    ## (1,R2) 0.156317460          NA
    ## (2,R2)          NA 0.000000001
    ## (3,R2)          NA 0.097033417
    ## (4,R2) 0.288173003          NA
    ## 
    ## Arrived at a reliable solution

``` r
plotMKmodel(MK_sub_hidden, display = "row")
```

![](State_transitions_analysis_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Visualization of transitions

### Prep data for visualization

Pick 8 colours for plotting the different transitions on the tree and
re-label to the corresponding states for ploidy and symbiosis.

``` r
# Load colour palettes 
state_colours1<- inauguration("inauguration_2021")
state_colours2<- inauguration("bernie_mittens")

# Get the colour codes 
state_colours1[c(1,2,3,5)]
```

    ## [1] "#5445b1" "#749dae" "#f3c483" "#cd3341"

``` r
state_colours2[c(1,2,3,4)]
```

    ## [1] "#372421" "#50506D" "#855F4C" "#465952"

``` r
# Combine colours for for plotting 
state_colours8<-setNames(c("#cd3341", "#749dae", "#5445b1", "#372421", "#50506D", "#465952", "#f3c483", "#855F4C"),
                         c("1", "2", "3", "4", "5", "6", "7", "8")) 

# Colours and labels for the legend
legend_colours<-setNames(c("#cd3341", "#749dae", "#5445b1", "#372421", "#50506D", "#465952", "#f3c483", "#855F4C"),
                         c("non-symbiotic, diploid, R1", "non-symbiotic, polyploid, R1", "symbiotic, diploid, R1", "symbiotic, polyploid, R1", "non-symbiotic, diploid, R2", "non-symbiotic, polyploid, R2", "symbiotic, diploid, R2", "symbiotic, polyploid, R2")) 
```

### Genus/Subfamily level ploidy data

We plotted the transition rates determined from the hidden state model
which was the best fitting model, for the largest dataset (genus +
subfamily level ploidy).

``` r
# Extract the data from the model 
phy_sub_hidden = MK_sub_hidden$phy
data_sub_hidden = MK_sub_hidden$data
model_sub_hidden = MK_sub_hidden$solution
model_sub_hidden[is.na(model_sub_hidden)] <- 0
diag(model_sub_hidden) <- -rowSums(model_sub_hidden)

# Run to get simmap 
simmap_sub_hidden <- makeSimmap(tree = phy_sub_hidden, data = data_sub_hidden, model = model_sub_hidden, rate.cat = 2)

# Import phytools plotSimmap for plotting
phytools::plotSimmap(simmap_sub_hidden[[1]], colors=state_colours8, fsize = 0.05, type="fan")
phytools::add.simmap.legend(colors=legend_colours, fsize=0.7, prompt=FALSE, x=-120, y=70, vertical=TRUE)
```

![](State_transitions_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->