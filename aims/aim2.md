# Aim 2

**To evaluate the extent to which low-dimensional representations of single cell data learned by various methods can be used to decompose tissues in the context of autoimmune/rheumatic diseases.**

## Background

Human Cell Atlas' partnership with the Immunological Genome Project (immgenH) will result in single-cell gene expression-based immunocyte phenotyping at unprecedented resolution.
A compendium comprised of bulk gene expression data from autoimmune/rheumatic diseases is exceptionally well-suited for the biologically-motivated evaluation (e.g., benchmarking) of these immunocyte data.
The _objective of this aim_ is to build such a compendium and to evaluate methods, including VAEs and other methods developed through this initiative, for defining cell-type-specific expression signatures from single-cell data by measuring their performance in the context of decomposing bulk, whole-tissue autoimmune/rheumatic disease data.
The quality of the cell-type-specific signatures (and by extension the methods) will be measured by using this information in deconvolution frameworks and in correctly identifying cell-type proportion and phenotype in the context of both simulation and real data experiments.

## Approach

We will generate simulated bulk datasets drawn from purified cell lineages, combining the purified cell expression data at different proportions.
We will use single-cell data-derived cell-type-specific signatures from different methods (e.g., VAEs, other methods developed through this initiative) as input to methods for decomposing bulk data and evaluate performance (e.g., prediction of proportions). 
Bulk transcriptomes of homogenous cell-types perturbed with many stimuli, such as those proposed by Arjun Raj's group, are well-suited for this experiment.
We will use these bulk data to compare cell-type proportion estimates across methods when highly similar cell-states (i.e., the same cell type under different conditions) are simulated at varying proportions.

We estimate the proposed multi-tissue autoimmune/rheumatic disease compendium (real data) will contain over ten thousand samples and include samples from patients with systemic lupus erythematosus (SLE), sarcoidosis, and inflammatory bowel disorders (to name a few). 
This compendium will have desirable properties for evaluating single-cell data-derived signatures, as it will allow us to evaluate cell type proportions and phenotypes based on a body of previous literature. 
For instance, we expect to detect higher proportions of activated macrophages in lupus nephritis samples as compared to controls [[1](https://dx.doi.org/10.4049/jimmunol.1103031)].
In addition, the compendium will include drug studies of highly-targeted therapeutics (e.g., a monoclonal antibody to IFN-gamma [anti-IFNg] in the context of systemic lupus erythematosus [[2](https://dx.doi.org/10.1002/art.39248)]).
In the anti-IFNg example above, we might ask what cell-types change in proportion during the reduction of this cytokine or if the IFN-inducible expression is preferentially altered in one cell type.
Such experiments allow us to not only assess methods for defining cell-type from single-cell data, but also methods that may be useful in decomposing whole-tissue bulk data that are developed through this initiative such as latent space arithmetic (see aim 3, below).
