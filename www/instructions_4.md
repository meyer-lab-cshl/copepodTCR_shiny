
**Right peptide** and **Right epitope** columns check whether the scheme remains interpretable upon pool drop-out. **Right peptide** should always contain the value “True”; otherwise, recovery was unsuccessful.

Results of the simulation can be assesed with help of generated histplot.
- it shows **# of peptides** for each combination of activated pools; generally this number should correspond to number of peptides sharing the same epitope; two end-position peptides (the first and the last ones) would have lesser number of peptides; all other results mean inconsistent overlap length across peptide list;
- it shows whether there are epitopes erroneously recovered after a drop-out (**Right peptide** = False).

### <span style="color:#015396">STL files</span>

To avoid mixing pools manually, the user might print special punched cards using files with their 3D models.

Each card represents one pool, with holes positioned at the coordinates corresponding to the peptides designated for addition to that pool. 

Produced punched card is placed on the empty tip box, and open holes are filled with tips. This patterned pippette tip array is used to transfer peptides from the plate to the corresponding pool.

The user can adjust parameters to fit their plate: