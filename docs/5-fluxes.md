# __Fluxes__ module

Performes a flux analysis for a target species 
at the final reaction time.

## Controls

* `Compute fluxes` starts the computation of fluxes.
This is a time-consuming computation and is not started automatically.
__It has to be done once__. 
The results shown in the right tabs depend only on the other controls.

* `Choose species` specifies the target species

* `Threshold` enables to limit the number of reactions 
kept in the results

* `Curved arrows` replaces straight arrows by curved ones in 
the graph representation (`ViewFlow`)

* `Graph algorithm` enables a choice of representations
for the graph (`ViewFlow`)

* `Nodes attraction` repulsion factor between nodes (`ViewFlow-D3`)


## Outputs

### __Budget/Target__ tab

List of the main production and loss paths for the target species. 
The main production path from the initial mixture is also estimated.
The relative weight/contribution of each reaction is shown;
it is based on the fluxes.

### __ViewFlow-D3__ tab

An interactive graph representation of the fluxes, 
where the size of arrows increases with the flux.

The reactions are tagged by their index given in 
`Model > Chemistry > Reactions`.


### __ViewFlow__ tab

A graph representation of the fluxes, where the size of arrows 
increases with the flux.
Red arrows correspond to formation reactions of the target product,
and blue arrows to its loss reactions.

The reactions are tagged by their index given in 
`Model > Chemistry > Reactions`,
and color coded (orange: photoprocesses; green: reactions).

