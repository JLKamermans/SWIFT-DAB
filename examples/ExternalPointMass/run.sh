#!/bin/bash

# Generate the initial conditions if they are not present.
if [ ! -e Sphere.hdf5 ]
then
    echo "Generating initial conditions for the point mass potential box example..."
    python makeIC.py 10000
fi

../swift -g -t 2 externalPointMass.yml