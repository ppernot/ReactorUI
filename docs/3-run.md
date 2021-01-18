# __Run__ module

Once the model is defined, one can run the `reactor` code to
simulate the photochemistry over the defined reaction time.
Two controls are available:

* `# MC runs` defines the number of Monte-Carlo samples
of the chemistry for which the code is to be run.

    + if set to 0, only the nominal run is performed.

    + for positive numbers, an additional control appears,
    asking if one wishes to append those runs to the exisiting
    ones. In this case, the sample numbers will be chosen 
    accordingly, to avoid double runs with the same sample...

* `# snapshots per run` defines the number of snapshots of the
species concentrations to be saved along the simulation. 
The points are regularly spaced on the logarithmic time scale.

* `Log relative error` and `Log absolute error` control
the error levels of the  IRKC implicit-explicit integrator
(see <https://github.com/ppernot/Reactor>)

    __Note__: the default values are chosen to provide rapidly 
    reasonable results. Fot a better accuracy on the minor species,
    it might be better to lower these values.

* `Use SR` and `Max. SR`: by default (`Use SR` unchecked), 
the reactor code estimates the spectral radius (SR) of the 
jacobian of the explicit contribution 
(see <https://github.com/ppernot/Reactor>). The estimation can fail,
with a "Convergence failure in estimating the spectral radius."
message. To circumvent this, one might shortcut the automatized 
estimation of SR by checking `Use SR` and providing a value in `Max. SR`.
If an appropriate integration statistics file exists, 
and if the reported value of SPRAD is not constant, SRmax
is estimated as 1.2*max(SPRAD). This choice can be overridden
by providing a value in `Max. SR`. If available, evolution of SPRAD can 
be seen in `Analysis>Sanity checks>Spectral radius`.


## __Local__ tab

* `Log tail only` enables to see the last lines of the reactor code
log when ckecked or the full log when unchecked. It is checked by 
default to enable to visualize the progress of the task.

The `Start` button in the `Local` tab will launch the simulations
on the local computer. 
The text outputs and error messages should be captured and displayed 
in the right panels.

__Warning__: it might happen that no output is displayed, although
the simulations are perfectly running in the background.
If this happens, check that the `MC_Output` folder is being 
populated according to your expectations...

## __Cloud__ tab

There is provision for a `Cloud` mode to launch the
simulations on the Virtual Data cloud. To be done...
For now, one needs to do it externally from the ReactorUI interface.

__Note__: to run `reactor` on the Virtual Data openstack cloud:

1. perform a nominal run on the local machine to ensure 
and check that the project is correctly configured;

2. follow the manual procedure explained on
[CloudReactor](https://github.com/ppernot/CloudReactor) 
in this
[document](https://github.com/ppernot/CloudReactor/blob/master/Doc/CloudReactor.pdf);

3. when the runs are done, proceed to the [Analysis](4-analysis.html) tab.

