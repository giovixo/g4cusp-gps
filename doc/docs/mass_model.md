# Mass Model

The mass model is defined by a GDML file, exported from the CUSP freeCAD model. The current mass model is the git branch `20240502_mass_model`. 

**Refactoring needed**
```txt
The GDML mass model must be carefully tested and validated. Currently, the code uses 
different mass models, each associated with a specific Git branch. This is certainly 
not a good code organization and will need to be revised. The main branch should 
contain the mass model in use, while the other Git branches should be used to 
implement various experimental and non-experimental features which, if appropriate, 
may be introduced into the main branch.
```