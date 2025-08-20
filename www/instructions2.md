
- **Peptide** — peptide sequence
- **Address** — pool indices where this peptide should be added
- **Epitope** — checked epitope from this peptide
- **Act pools** — list with pool indices where this epitope is present
- **# of act pools** — number of pools where this epitope is present
- **# of epitopes with these act pools** — number of epitopes that are present in the same pools (= number of possible epitopes upon activation of such pools)
- **# of peptides with these act pools** — number of peptides in which there are epitopes that are present in the same pools (= number of possible peptides upon activation of such pools)

**# of peptides with these act pools** should be equal to number peptides sharing an epitope. For end-position peptides, it would be less. However, if for some epitopes **# of peptides with these act pools** is bigger, than these peptides have bigger overlap than others.

To interpret the results of the experiment, user can find all rows where the **Act Pools** column contains respective combination of activated pools. This way, all possible peptides and epitopes leading to the activation of such a combination of pools are obtained. But we recommend using **Interpretation** tab to interpret results of the experiment.

### <span style="color:#015396">3D masks</span>

To avoid mixing pools manually, the user might print special punched cards using files with their 3D models.

Each card represents one pool, with holes positioned at the coordinates corresponding to the peptides designated for addition to that pool. 

Produced punched card is placed on the empty tip box, and open holes are filled with tips. This patterned pippette tip array is used to transfer peptides from the plate to the corresponding pool.

The user can adjust parameters to fit their plate: