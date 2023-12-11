- **Peptide** — peptide sequence
- **Address** — pool indices where this peptide should be added
- **Epitope** — checked epitope from this peptide
- **Act pools** — list with pool indices where this epitope is present
- **# of pools** — number of pools where this epitope is present
- **# of epitopes** — number of epitopes that are present in the same pools (= number of possible epitopes upon activation of such pools)
- **# of peptides** — number of peptides in which there are epitopes that are present in the same pools (= number of possible peptides upon activation of such pools)
- **Remained** — only upon regime=”with dropouts”, list of pools left after a drop-out
- **# of lost** — only upon regime=”with dropouts”, number of dropped pools
- **Right peptide** — True or False, whether the peptide is present in the list of possible peptides
- **Right epitope** — True or False, whether the peptide is present in the list of possible peptides

To interpret the results of the experiment, user can find all rows where the **Act Pools** column contains respective combination of activated pools. This way, all possible peptides and epitopes leading to the activation of such a combination of pools are obtained.

Without drop-outs, **# of peptides** should be equal to number peptides sharing an epitope. For end-position peptides, it would be less. However, if for some epitopes **# of peptides** is bigger, than these peptides have bigger overlap than others.

If observed combination of activated pools is not present in the table, simulation with drop-outs can be checked.

#### <span style="color:#015396">With drop-outs</span>

It would look as follows: