# <span style="color:#015396">Usage</span>

The tool consists of four distinctive parts, each of which can be used separately. These four parts constitute the stages of the experiment design and implementation.

### <span style="color:#015396">Balance check</span>

First, the appropriate number of pools and peptide occurrence (number of pools per one peptide) should be selected.

Peptide occurrence affects number of peptides in one pool, and therefore too high peptide occurrence may lead to higher dilution of a single peptide. In Kovaleva et al (2023), we were able to detect signal with the cognate peptide diluted to 1.58μM.

Peptide dilution can be mitigated by increasing number of pools, however, it might increase the complexity of experimental set-up. Consequently, these parameters should be chosen carefully.

To assist with this process, copepodTCR provides the user with possible peptide occurrence values based on given number of pools and number of tested peptides. Also, copepodTCR calculates the resulting distribution of peptides in pools and compares this distribution with perfect scenario, where number of peptides per pools is completely balanced.

### <span style="color:#015396">Overlap check</span>

Inconsistent overlap length in the list of tested peptides can lead to imprecise results interpretation. In copepodTCR, the user can check overlap consistency across the entire list of tested peptides.

Upon peptides entry, copepodTCR returns the list of present overlap lengths and corresponding number of peptide pairs. If two overlap lengths are observed, it also returns peptide pairs with all overlap lengths differing from the most common one.

### <span style="color:#015396">Pooling scheme</span>

Upon parameters selection and peptides check, the user can enter them into copepodTCR and get a peptide pooling scheme.

CopepodTCR returns three tables:
- peptide pooling scheme pool-wise (i.e. the table with peptides in each pool)
- peptide pooling scheme peptide-wise (i.e. the table with pools for each peptide)
- simulation table

During simulation step, copepodTCR simulates results of the experiment for any possible epitope of the provided length and returns a table with every possible epitope and all pools where this epitope is present.

The function has two regimes: with and without drop-outs. Without drop-outs, it returns a table as there were no experimental mistakes, and there were zero erroneous non-activated pools. With drop-outs, it returns a table with all possible mistakes (i.e. all possible erroneous non-activated pools). This option needs time to be generated, usually several minutes, although it depends on the number of peptides and peptide occurrence.

#### <span style="color:#015396">Without drop-outs</span>

Simulation table  without drop-outs looks as follows: