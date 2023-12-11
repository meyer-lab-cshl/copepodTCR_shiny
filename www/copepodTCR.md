<span style="color:white"> </span>
## <span style="color:#015396">COmbinatorial PEptide POoling Design for TCR specificity</span>
<span style="color:#015396">CopepodTCR</span> is a tool for the design of combinatorial peptide pooling schemes for TCR speficity assays.

<span style="color:#015396">CopepodTCR</span> guides the user through all stages of the experiment design and interpetation:
- selection of parameters for the experiment (**Balance check**)
- examination of peptides (**Overlap check**)
- generation of pooling scheme (**Pooling scheme**)
- generation of punched cards of efficient peptide mixing (**STL files**)
- results interpetation using hierarchical Bayesian model (**Activated pools**)

# <span style="color:#015396">Cite as</span>

Kovaleva V. A., et al. "copepodTCR: Identification of Antigen-Specific T Cell Receptors with combinatorial peptide pooling." bioRxiv (2023): 2023-11.

Or use the following BibTeX entry:

```
@article{
	kovaleva2023copepodtcr,
	title        = {copepodTCR: Identification of Antigen-Specific T Cell Receptors with combinatorial peptide pooling},
	author       = {Kovaleva, Vasilisa A and Pattinson, David J and Barton, Carl and Chapin, Sarah R and Minerva, Anastasia A and Richards, Katherine A and Sant, Andrea J and Thomas, Paul G and Pogorelyy, Mikhail V and Meyer, Hannah V},
	year         = 2023,
	journal      = {bioRxiv},
	publisher    = {Cold Spring Harbor Laboratory},
	pages        = {2023--11}
}
```

# <span style="color:#015396">Description</span>

Identification of a cognate peptide for TCR of interest is crucial for biomedical research. Current computational efforts for TCR specificity did not produce reliable tool, so testing of large peptide libraries against a T cell bearing TCR of interest remains the main approach in the field.

Testing each peptide against a TCR is reagent- and time-consuming. More efficient approach is peptide mixing in pools according to a combinatorial scheme. Each peptide is added to a unique subset of pools ("address"), which leads to matching activation patterns in T cells stimulated by combinatorial pools.

Efficient combinatorial peptide pooling (CPP) scheme must implement:
- use of overlapping peptide in the assay to cover the whole protein space;
- error detection.

Here, we present <span style="color:#015396">CopepodTCR</span> -- a tool for design of CPP schemes. CopepodTCR detects experimental errors and, coupled with a hierarchical Bayesian model for unbiased results interpretation, identifies the response-eliciting peptide for a TCR of interest out of hundreds of peptides tested using a simple experimental set-up.