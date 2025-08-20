<span style="color:white"> </span>
## <span style="color:#015396">COmbinatorial PEptide POoling Design for TCR specificity</span>
<span style="color:#015396">CopepodTCR</span> is a tool for the design of combinatorial peptide pooling schemes for TCR speficity assays.

<span style="color:#015396">CopepodTCR</span> guides the user through all stages of the experiment design and interpetation:
- selection of parameters for the experiment (**Balance check**)
- generation and assessment of peptides (**Peptide check**)
- generation of pooling scheme (**Pooling scheme**)
- generation of punched cards of efficient peptide mixing (**3D masks**)
- results interpetation using hierarchical Bayesian model (**Interpretation**)

# <span style="color:#015396">Cite as</span>

Kovaleva V. A., et al. "copepodTCR: Identification of Antigen-Specific T Cell Receptors with combinatorial peptide pooling." bioRxiv (2023): 2023-11.

Or use the following BibTeX entry:

```
@article{
	kovaleva2023copepodtcr,
	title        = {copepodTCR: Identification of Antigen-Specific T Cell Receptors with combinatorial peptide pooling},
	author       = {Kovaleva, Vasilisa A and Pattinson, David J and Barton, Carl and Chapin, Sarah R and Minervina, Anastasia A and Richards, Katherine A and Sant, Andrea J and Thomas, Paul G and Pogorelyy, Mikhail V and Meyer, Hannah V},
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

The experimental setup starts with defining the protein/proteome of interest and obtaining synthetic peptides tiling its space. Peptide sequences can be generated in silico from a protein of interest and then checked using **Peptide check** tab.
This set of peptides, containing an overlap of a constant length, is entered into <span style="color:#015396">copepodTCR</span>. Parameters for CPP scheme can be selected using **Balance check** tab. To create a pooling scheme with selected parameters, use **Pooling scheme** tab. To make pipetting easier, you can use **3D masks** tab to create mask models which could be further 3D printed and overlay the physical plate or pipette tip box.

Following this scheme, the peptides should be mixed, and the resulting peptide pools should be tested in a T cell activation assay. The activation of T cells can be measured for each peptide pool with the assay of choice, such as flow cytometry- or microscopy-based activation assays detecting transcription and translation of a reporter gene. The experimental measurements for each pool can be entered back into <span style="color:#015396">copepodTCR</span> which employs a Bayesian mixture model to identify activated pools. Based on the activation patterns, it returns the set of overlapping peptides leading to T cell activation (**Interpretation** tab).

For more details, refer to “copepodTCR: Identification of Antigen-Specific T Cell Receptors with combinatorial peptide pooling” ([bioRxiv version](https://www.biorxiv.org/content/10.1101/2023.11.28.569052v1)).

# <span style="color:#015396">CopepodTCR python package</span>

Alternatively, you might want to use <span style="color:#015396">copepodTCR</span> python package as it provides more flexibility. 

It can be installed with pip:
```python
pip install copepodTCR
```

or conda: 
```python
conda install -c vasilisa.kovaleva copepodTCR
```

Documentation: [copepodTCR.readthedocs](https://copepodtcr.readthedocs.io/en/latest/index.html).
