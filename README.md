[![DOI](https://zenodo.org/badge/165652991.svg)](https://zenodo.org/badge/latestdoi/165652991)

# ReactorUI
Graphical User Interface for the [reactor](https://github.com/ppernot/Reactor) 
code.

## How to cite

> Z. Peng, N. Carrasco and P. Pernot (2014) 
> "Modeling of synchrotron-based laboratory simulations of 
> Titan's ionospheric photochemistry", _GeoResJ_ __1-2__:33–53.
> <http://dx.doi.org/10.1016/j.grj.2014.03.002>

and
    
> Pernot, P. (2020) ReactorUI: graphical UI for the 
> simulation of photochemical reactors
> (Version 1.3). <https://doi.org/10.5281/zenodo.3946078>
   

## New Release (v_1.3)

* Enables the use of a local chemical DB (for nominal run)

    
## Docker container

The [reactorui](https://hub.docker.com/repository/docker/ppernot1/reactorui)
[Docker](https://www.docker.com/) container has all elements preinstalled.

To run the container:

0. Install [Docker](https://www.docker.com/products/docker-desktop)

1. Type the following commands in a terminal
```
  cd Projects    # This is the home of your `reactor` projects   

  docker run -d -p 3838:3838 --mount type=bind,source=".",target=/Projects --name reactorui ppernot1/reactorui
```
If you want to use a local database (e.g. in repertory `myChemDB`, at the same level as `Projects`), 
it should also be mounted
```
  docker run -d -p 3838:3838 --mount type=bind,source=".",target=/Projects --mount type=bind,source="$(pwd)"/../myChemDB,target=/ChemDBLocal --name reactorui ppernot1/reactorui
```
    
2. Access ReactorUI at http://localhost:3838 in your favorite browser

3. When finished
```
  docker kill reactorui
```

4. For further sessions
```
  docker restart reactorui
```

4. To cleanup
```
  docker remove -v reactorui
```
