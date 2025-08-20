# <span style="color:#015396">Usage</span>

The tool consists of five distinctive parts, each of which can be used separately. Each part corresponds to a step in the CPP experimental workflow.

- selection of parameters for the experiment (**Balance check**)
- generation and assessment of peptides (**Peptide check**)
- generation of pooling scheme (**Pooling scheme**)
- generation of punched cards of efficient peptide mixing (**3D masks**)
- results interpetation using hierarchical Bayesian model (**Interpretation**)

### <span style="color:#015396">Balance check</span>

First, the appropriate number of pools and peptide occurrence (number of pools per one peptide) should be selected.

Peptide occurrence affects number of peptides in one pool, and therefore too high peptide occurrence may lead to higher dilution of a single peptide. In Kovaleva et al (2023), we were able to detect signal with the cognate peptide diluted to 1.58μM.

Peptide dilution can be mitigated by increasing number of pools, however, it might increase the complexity of experimental set-up. Consequently, these parameters should be chosen carefully.

To assist with this process, copepodTCR provides the user with possible peptide occurrence values based on given number of pools and number of tested peptides.

Also, we discovered that sensitivity and specificity of the developed Bayesian mixture model for identification of activated pools depend on the negative share — the expected number of non-activated pools.
After setting the desired peptide occurrence, the approximate negative share is calculated as (N - I - 1)/N, where N is the total number of pools and I is the peptide occurrence (it is assumed that one epitope is shared between two peptides).
For exploratory screens, we recommend choosing peptide occurrence resulting in a negative share equal to 0.5–0.6. For focused screens with the prior assumption that a cognate epitope is present in the tested peptide library, peptide occurrence leading to a higher negative share can be chosen.
The app then simulates the pooling and shows how close the peptide distribution is to a perfectly balanced scheme.


### <span style="color:#015396">Peptide check</span>

Here you can generate peptides from an entered amino acid sequence or check the peptides you already have for their overlap length consistency. Inconsistent overlap length in the list of tested peptides can lead to imprecise results interpretation.

- if you have protein sequence(s):  Paste the sequence or upload a FASTA file. Enter the desired peptide length and the overlap length. Select what to do when the last peptide to be generated is too short: either discard this peptide or allow a bigger overlap with the second-to-last peptide. Click "Generate" to see the peptide list, use "Download" or "Send" to export the peptides.

- if you have peptides: Paste or upload a TSV/CSV file with one peptide per line (without a header), then click "Check". Then, copepodTCR returns the distribution of sequence and overlap lengths in your peptide set. It also returns peptide pairs with overlap lengths differing from the most common one.


### <span style="color:#015396">Pooling scheme</span>

Upon parameters selection and peptides check, the user can enter them into copepodTCR and get a peptide pooling scheme.

CopepodTCR returns three tables:
- peptide pooling scheme pool-wise (i.e. the table with peptides in each pool)
- peptide pooling scheme peptide-wise (i.e. the table with pools for each peptide)
- simulation table

During simulation step, copepodTCR simulates results of the experiment for any possible epitope of the provided length and returns a table with every possible epitope and all pools where this epitope is present.

Simulation table looks as follows:
