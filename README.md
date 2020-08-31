[![DOI](https://zenodo.org/badge/165652991.svg)](https://zenodo.org/badge/latestdoi/165652991)

# ReactorUI
Graphical User Interface for the [reactor](https://github.com/ppernot/Reactor) 
code.

## New Release (v_1.0)

* First release with full functionalities

    + Create new projects / open existing ones
    
    + Modify project parameters (chemistry, irradiation, reactor)
    
    + Run the [reactor](https://github.com/ppernot/Reactor) code 
    (local execution only at the moment)
    
    + Analyse the results

* Requires installation of [ChemDBPublic](http://dx.doi.org/10.5281/zenodo.3946044)  
  from Zenodo. 
  
    + Follow the instruction in ReadMe.txt 

## Docker container

The [reactorui](https://hub.docker.com/repository/docker/ppernot1/reactorui)
[Docker](https://www.docker.com/) container has all elements preinstalled.
It avoids you the trouble of installing `R`, `reactor` and `chemDBPublic`.

To run the container:

0. Install [Docker](https://www.docker.com/products/docker-desktop)

1. Type the following commands in a terminal
```
cd Projects    
docker run -d -p 3838:3838 --mount type=bind,source=".",target=/Projects \
  --name reactorui ppernot1/reactorui
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
