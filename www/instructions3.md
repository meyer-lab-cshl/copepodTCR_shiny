- number of rows — number of rows in the plate
- number of columns — number of columns in the plate
- length — length of the plate (in mm)
- width — width of the plate (in mm)
- thickness — thickness of the plate (in mm)
- hole radius — diameter of the well divided by 2
- X offset — margin along the X axis for the A1 well, in mm
- Y offset — margin along the Y axis for the A1 well, in mm
- well spacing — distance between wells, in mm

To better orient tip pattern, the user can add the last hole (with coordinates m-k). It should be used only in absence of peptide in the corresponding well. Also, the user can choose to add marks indicating pool index to the masks.

### <span style="color:#015396">Activated pools</span>

The experiment can be conducted using flow cytometry or microscopy.

After the experiment, copepodTCR can help with data analyzation.
The primitive version of experiment interpetation is decscribed in section Pooling scheme, in the explanation of simulation step.

But also user can analyze the results using Bayesian Mixture model. This model returns the probability of each pool being activated (green) or not (gray).