# __Run__ module

Once the model is defined, one can run the `reactor` code to
simulate the photochemistry over the defined reaction time.
Two controls are available:

* `# MC runs` defines the number of Monte-Carlo samples
of the chemistry for which the code is to be run.

    + if set to 0, only the nominal run is performed.

    + for positive numbers, an additional control appears,
    asking if one wishes to append thos runs to the exisiting
    ones. In this case, the sample numbers will be chosen 
    accordingly, to avoid double runs with the same sample...

* `# snapshots per run` defines the number of snapshots of the
species concentrations will be saved along the simulation. 
The points are regularly spaced on the logarithmic time scale.

## `Local` tab

The `Start` button in the `Local` tab will launch the simulations
on the local computer. 
The text outputs and error messages should be captured and displayed 
in the right panels.

__Warning__: it might happen that no output is displayed, although
the simulations are perfectly running in the background.
If this happens, check that the `MC_Output` folder is being 
populated according to your expectations...

## `Cloud` tab

There is provision for a `Cloud` mode to launch the
simulations on the Virtual Data cloud. To be done...

To run `reactor` on the Virtual Data openstack cloud:

1. perform a nominal run on the local machine to ensure 
and check that the project is correctly configured;

2. follow the manual procedure explained on 
<https://github.com/ppernot/CloudReactor> in this
[document](https://github.com/ppernot/CloudReactor/blob/master/Doc/CloudReactor.pdf);

3. when the runs are done, proceed to the [Analysis](4-analysis.html) tab.

