# __Analysis__ module

Once the simulation results are available (or a part of them),
they can be analyzed to plot kinetic traces, mass-spectra or
perform sensitivity analysis. 

The analysis of reaction fluxes is done in the `Fluxes` module.

* `Load MC results` has first to be pressed to read the results files.
A short summary is presented below the button when the operation 
is completed.

__News__

* 2020_09_23: Statistics about the mole fractions at final simulation 
time are saved in `Outputs/MoleFrac_tFinal.csv`


## __Kinetics__ tab 

A plot of the time-dependent mole fractions 
of the chemical species and their uncertainty bands.
The plot is zoomable (select the area and double-click; 
double-click again to cancel the zoom).

* `Show error bands` to display the 95% uncertainty bands
resulting from the Monte-Carlo runs.
    
* `Fixed colors` to attribute a fixed color to each species,
independently of the set of species to be displayed.
    
* `Draw PPM scale` to add a series of horizontal lines showing 
the ppm, ppb and ppt mole fraction levels.
    
* The following controls enable to select the classes of
species to be displayed. Classes `Cx` refer to the number
of heavy atoms in a species.
    
## __Pseudo-MS__ tab

A plot of (pseudo) mass spectra. 
The MS for neutrals and ions are presented separately.
__For neutral species, these are not electronic impact (EI)
fragmentation mass spectra, but abundances at the unfragmented 
species mass.__

* `Show error bars` to display the 95% uncertainty bars
resulting from the Monte-Carlo runs.
    
* `Draw PPM scale` to add a series of horizontal lines showing 
the ppm, ppb and ppt mole fraction levels.

* The following controls enable to select the classes of
species to be displayed. Classes `Cx` refer to the number
of heavy atoms in a species.
    
* `Log sampling time` enables to define the time-snapshot
to be used to draw the MS.
    
* `MS Amplitude Range` controls the range of the y-axis 
of the mass-spectra.
    
* `Adjust bar width` enables to tune the width of the MS bars.

* `Size of names` controls the peaks labels size (0 = no labels)
    
## __Sensitivity__ tab

This is where one performs sensitivity analysis (SA)
to identify key reactions. In the Monte-Carlo framework, key reactions
are those which have the largest impact on the prediction uncertainty
of a quantity of interest. 
This will run only if MC runs have been performed.

* `Sensitivity indices` presents a choice of  sensitivity index (SI)
quantifying the dependence between the concentrations and reaction 
rates over the Monte-Carlo sample:

    - `Rank correl.`: Spearman correlation coefficient.
     
    - `dCorr`: distance correlation t-test of multivariate 
    and functional independence (Szekely, G.J. and Rizzo, M.L. (2013). 
    The distance correlation t-test of independence in high dimension. 
    Journal of Multivariate Analysis, Volume 117, pp. 193-213.
    [doi](https://doi.org/10.1016/j.jmva.2013.02.012) ).
    
    - `dHSIC`: d-variable Hilbert Schmidt independence criterion 
    (Gretton, A., K. Fukumizu, C. H. Teo, L. Song, B. Scholkopf 
    and A. J. Smola (2007). A kernel statistical test of independence. 
    In Advances in Neural Information Processing Systems, pp. 585-592.
    [doi](https://dl.acm.org/doi/10.5555/2981562.2981636) ).
    
    __Note__: the `dCorr` and `dHSIC` features are experimental...
    
* `Choose species`: target species for the SA. 
If no species is specified, the sensitivity indices are estimated 
from the sum of the "correlation" matrix over all species. 
For `Rank correl.`, the square of the correlation matrix is used.

    __Warning__: the computation of `dCorr` and `dHSIC` indices might be very
    time consuming, especially if no target species is specified.

* `Run SA`: start computation of sensitivity indices
    
* `Plot type`: 
    
    - `Bar Plot` shows a barplot of the Top 20 reactions 
    with the larges sensitivity index, whether a target species
    is specified or not.
    
    - `Scatterplot` is active only if a target species is specified.
    It draws scatterplots between the MC concentrations of the
    target species and the 20 MC reaction rates with largest SI. 
    The value of the SI is reported in red.
  
## __Sanity checks__ tab


### __Outputs__ tab    

This module checks in MC samples if final concentrations, photorates 
or reaction rates present numerical exceptions (zero, infinity, NaN...). 

The following statistics are reported for each problematic variable:

* `Var=0`: a logical indicaticating if the variance over the MC sample is null

* `Nzero`: fraction of the MC samples with a null value

* `Ninf`: fraction of the MC samples with infinite value

### __Integration__ tab
 
This module plots statistics about the IRKC implicit-explicit integrator
(_cf._ <https://github.com/ppernot/Reactor>):

* `NFE` :  number of evaluations of F_E

* `NFI` :  number of evaluations of F_I

* `NACCPT` :  number of accepted steps

* `NREJCT` :  number of rejected steps

* `NFESIG` : number of evaluations of F_E used 
  in estimating the spectral radius
  
* `MAXM` : maximum number of stages used

### __Spectral radius__ tab
 
This module plots statistics about the IRKC implicit-explicit integrator
(_cf._ <https://github.com/ppernot/Reactor>):

* `SPRAD` :  spectral radius of the jacobian of F_E
