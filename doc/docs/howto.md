# How To 

## How to compile the GEANT4 application

The cmake configuration file is CMakeLists.txt

To build the application:

```c
mkdir -p build
cd build
cmake ..
make
```

## How to build the documentation

This documentation is based on [MkDocs](https://www.mkdocs.org/)

Go first into the `doc` directory.

This will build the static site:

`mkdocs build`

This will start the server on http://127.0.0.1:8000/

`mkdocs serve`

## How to deploy the documentation

Use the following command line:

`mkdocs gh-deploy`

This will:

* Build the static site

* Push it to the gh-pages branch of your repository

Wait a coulple of minutes to see the changes on https://giovixo.github.io/g4cusp-gps  