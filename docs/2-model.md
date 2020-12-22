# __Model__ module


The left panel displays the current values of the project's parameters
(initialized from the `control.dat` file present in the project's folder
or the chosen template for new projects).

The right panel contains a series of tabs for the different parts of
the model.


## __Chemistry__ tab

This is where the chemical scheme/network is set-up.

* `Generate`: generate the list of reactions corresponding to the initial
mixture.

    + the `Species` and  `Compo.` boxes enable to define the composition
    of the initial gaz mixture.
    
    + when the composition is typed in, hit the `Generate Reactions` button.
    The code reads the chemistry databases and elaborates iteratively
    a set of reactions specific to the species present in the initial 
    gaz mixture. 
    
    + the results of the Volpert analysis used to generate the reaction
    network are proveded in the right panel, showing the new species 
    added at each iteration. The Volpert index (VlpI) is the number of 
    steps separating a species from the initial mixture.
    
        ```
        # Pseudocode for the 'Volpert' generation of reaction scheme
        # Inputs: 
        # - the initial mixture of the gaz
        # - a database of reactions (including photoprocesses)
        # Outputs:
        # - a reduced and consistent list of reactions
        # - a list of accessible species
    
        1. VlpI = 0 : initialize species list with initial mixture
    
        2. VlPI = 1 :
           a. search and list photoprocess for species in the list
           b. augment the list with the photoprocess products
               (new species are tagged with VlpI = 1)
    
        3. VlpI = VlpI +1 
           a. search and list photoprocess and reactions involving 
              species in the list as reactants
           b. augment the list with the products of these reactions
              (new species are tagged with VlpI)
           c. if (no new species) 
                stop 
              else 
                iterate at 3. 
        ```
    

* `Network` presents a zoomable and active graphical view of the network. 
Two controls are available:

    + `Nodes attraction` controls the strength of the forces in the
    network representation.
    
    + `Max Volpert Index` enables to see the iterative building of the
    reactions network, from the initial mixture, to the final scheme.
    The species appearing at each Volpert iteration are color coded.
    
    + `Coloring/Clustering` offers a choice of options to color the nodes 
    of the graph:
    
        - `volpert` colors by the value of the volpert index
        
        - `charge` colors by the charge of the species
        
        - `edge_betweeneness`, `louvain`, `fast_greedy` and
        `leading_eigen` are community detection algorithms 
        proposed by the [igraph](https://igraph.org/r) package 
        (_experimental_; might be very slow for large networks)
    
    
* `Reactions` lists the reactions. 

    + The list can be copied to the clipboard or downloaded to disk.
    
    + Typing a species name into the `Target sp.` textbox on the left 
    selects the subset of reactions involving this species.

* `Checks` lists the species for which no loss reactions have been 
found. These will act as sinks.

* `Sample` is used to assemble Monte-Carlo samples from the generated network.
The data are gathered from `ChemDBPublic`, a repository of database samples
generated from `MC-ChemDB` (see <https://github.com/ppernot/MC-ChemDB>).

    + `#MC samples` enables to define the desired number 
    (_N_, currently 500 max). 
    
    + `Generate Samples` starts the generation of MC samples. 
    
    _N_+1 samples are generated, with numbers from 0 to _N_, where 0 is the
    nominal sample, containing the unperturbed reference values for
    the chemistry parameters.
    
    The generated files are stored in the `MC_Input` subfolder 
    of the project's folder. 
    __They will overwrite existing files with the same sample numbers.__


## __Irradiation__ tab

The spectrum, intensity and cross-section of the irradiation beam are
defined here.

* `Beam spectrum file` is used to choose a spectrum file on disk.

    __Note__ the spectral resolution has to be conform to the one
    declared for the photolysis cross-sections in the `ChemDB Versions`
    tab.

* `Predefined files` proposes a set of spectra.

* `Spectrum Range` enables to define the irradiation spectral range
(by default the range of the spectrum file, in nanometers)

* `Beam intensity (ph.cm^-2.s^-1)` defines the intensity of the beam.

* `Beam Section (cm^2)` defines the  beam's cross-section.

A summary is provided at the bottom of the controls column,
and a figure of the model spectrum is shown on the right.

## __Reactor__ tab

The reactor is considered as a tube in which a gaz flow
is maintained.

The dimensions of the tube (length, cross-section), 
the working temperatures (gaz and electrons),
the properties of the gaz flow (total pressure,
reactants pressure and reactants flux)
and the duration of the experiment are managed here.


## __ChemDB Versions__ tab

By default, the latest versions of the chemmistry databases
(Photoprocs, Neutrals and Ions) are used, but the user can 
control specific versions through this tab.

The resolution of the photoprocesses cross-sections can be switched
between `1 nm` (default) and `0.1 nm` (high resolution).

