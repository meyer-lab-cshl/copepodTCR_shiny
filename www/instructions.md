The experimental setup starts with defining the protein/proteome of interest and obtaining synthetic peptides tiling its space.

This set of peptides, containing an overlap of a constant length, is entered into copepodTCR. It creates a peptide pooling scheme and, optionally, provides the pipetting scheme to generate the desired pools as either 384-well plate layouts or punch card models which could be further 3D printed and overlay the physical plate or pipette tip box.

Following this scheme, the peptides are mixed, and the resulting peptide pools tested in a T cell activation assay. The activation of T cells is measured for each peptide pool (experimental layout, activation assay, and experimental read out) with the assay of choice, such as flow cytometry- or microscopy-based activation assays detecting transcription and translation of a reporter gene.

The experimental measurements for each pool are entered back into copepodTCR which employs a Bayesian mixture model to identify activated pools.  Based on the activation patterns, it returns the set of overlapping peptides leading to T cell activation (Results interpretation).

# <span style="color:#015396">Branch-and-Bound algorithm</span>

For detailed description of the algorithm and its development refer to Kovaleva et al (2023).

The Branch-and-Bound part of copepodTCR generates a peptide mixing scheme by optimizing the peptide distribution into a predefined number of pools <i>n</i> (in Figure <i>n=6</i>). The distribution of each peptide is encoded into an address (edges in the graph), which connect nodes in the graph (circles) that represent a union between two addresses. The peptide mixing scheme constitutes the path through these unions and connecting addresses that ensure a balanced pool design.

# <span style="color:#015396">Activation model</span>

For detailed description of the model, refer to Kovaleva et al (2023).

To accurately interpret results of T cell activation assay, copepodTCR utilizes a Bayesian mixture model.

The model considers the activation signal to be drawn from two distinct distributions arising from the activated and non-activated pools and provides the probabilities that the value was drawn from either distribution as a criterion for pool classification.

# <span style="color:#015396">CopepodTCR Python package</span>

Can be installed with
```python
pip install copepodTCR
```

or 
```python
conda install -c vasilisa.kovaleva copepodTCR
```

### <span style="color:#015396">Requirements</span>
Required packages should be installed simulataneously with the copepodTCR packages.

But if they were not, here is the list of requirements:
```python
	pip install "pandas>=1.5.3"
	pip install "numpy>=1.23.5"
	pip install "cvxpy>=1.3.2"
	pip install "trimesh>=3.23.5"
	pip install "trimesh>=3.23.5"
	pip install "pymc>=5.9.2"
	pip install "arviz>=0.16.1"
```