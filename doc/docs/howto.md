# How To 

## How to compile the GEANT4 application

```c
If needed edit the file CMakeLists.txt
mkdir -p build
cmake ..
make
```

## How to develop the documentation

`mkdocs build`

This will build the static site

`mkdocs serve`

This will start the server on http://127.0.0.1:8000/

## How to deploy the documentation

`mkdocs gh-deploy`

This will:

* Build the static site

* Push it to the gh-pages branch of your repository

Wait a coulple of minutes to see the changes on https://giovixo.github.io/g4cusp-gps  