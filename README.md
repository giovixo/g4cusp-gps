# g4cusp-gps

A GEANT4 simulation for CUSP

Source model: General Purpose Source (GPS)

# Git workflow

This repository contain the following branches:

```text
             ┌── gamma (branch mantained by Giovanni De Cesare)
             │
main ────────┼── activation (branch maintained by Riccardo Campana)
             │
             |── post-activation (branh mantained by Riccardo Campana)
             |
             └── gh-pages (documentation based on mkdocs. Branch mantained by Giovanni De Cesare)             
```

<u>VERY IMPORTANT: WORK IN THE BRANCH OF YOUR PROJECT</u>

For example if you project is post-activation, type in the terminal
> git checkout post-activation

Before start, verify that you are in the correct branch, typing:

> git branch

# User manual

The user manual is based on ([mkdocs](https://www.mkdocs.org/))

Localy, type:

> cd doc

> mkdocs serve

On line:

https://giovixo.github.io/g4cusp-gps/


