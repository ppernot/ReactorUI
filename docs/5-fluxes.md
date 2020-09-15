# __Fluxes__ module

Performes a flux analysis for a target species 
at the final reaction time.

## Controls

* `Compute fluxes` starts the computation of fluxes. 
__This has to be done once__. 
The results shown in the right tabs depend on the other controls.

* `Choose species` specify the target species

* `Threshold` enables to limit the number of reactions 
kept in the results

* `Curved arrows` replaces straight arrows by curved ones in 
the graph representation (`ViewFlow`)

* `Graph algorithm` enables a choice of representations
for the graph (`ViewFlow`)

## Outputs

* `ViewFlow` tab: a graph representation of the fluxes,
where the size of arrows increases with the flux.
Red arrows correspond to formation reactions for the target product,
and blue arrows to loss reactions.

    The reactions are tagged by their index given in 
    `Model > Chemistry > Reactions`,
    and color coded (orange: photoprocesses; green: reactions).

* `Budget/Target` lists the main production and loss paths for
the target species. 
The main production path from the initial mixture is also estimated.
The relative weight/contribution of each reaction is shown,
ased on the computed fluxes.
