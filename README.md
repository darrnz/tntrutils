# tntrtools

*Custom R functions for interacting with TNT*

This is a small collection of custom R utility functions for working with TNT. These are mostly complementary to Nick Matzke's TNTR.

The code is clunky and has not been tested extensively. Use with caution!

## Functions available

### writecont.tnt
Export a matrix of continuous data into the TNT matrix format for analysis of continuous characters "as such" (Goloboff et al. 2006).
### writeland.tnt
Export an array or list of arrays with 2D or 3D landmark coordinates into TNT matrix format for analysis with spatial optimisation (Catalano et al. 2010). Multiple arrays in a list will be interpreted as different configurations ("characters").
### write.tree.tnt
Export ape's phylo or multiPhylo objects into the TNT parenthetical tree format. Optionally, it creates a dummy matrix for reading the trees in TNT without previous data in memory.
