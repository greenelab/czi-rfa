### Aim 2: Evaluate the extent to which low-dimensional representations of single cell data learned by various methods can be used to decompose bulk tissues in the context of rheumatic diseases.**

The HCA's partnership with the Immunological Genome Project (immgenH) will provide single-cell gene expression-based immunocyte phenotyping at an unprecedented resolution.
Ultimately, we expect that HCA-identified cell populations will translate to an improved understanding of health and disease.
A compendium comprised of bulk gene expression data from autoimmune/rheumatic diseases is exceptionally well-suited for the biologically-motivated evaluation (e.g., benchmarking) of these immunocyte data.
The _objective of this aim_ is to build such a compendium and to evaluate methods, including VAEs and other methods developed through this initiative, for defining cell-type-specific expression signatures from the HCA's single-cell datasets by measuring their performance in the context of decomposing bulk, whole-tissue autoimmune/rheumatic disease data.
We will develop both real and simulated benchmark datasets to evaluate the quality of the cell-type-specific signatures (and by extension the methods) for deconvolution and the identification of correct cell-type proportions and phenotypes.

We will generate simulated bulk datasets drawn from purified cell lineages by combining the purified cell expression data at different proportions.
We will apply methods including existing methods, VAEs (Aim 1), and those developed through this initiative to generate single-cell data-derived cell-type-specific signatures.
We will input these signatures into methods for decomposing bulk data and evaluate performance (e.g., prediction of proportions).
Bulk transcriptomes of homogenous cell-types perturbed with many stimuli, such as those proposed by Arjun Raj's group, are well-suited for this experiment.
We will use these simulated benchmark data to compare cell-type proportion estimates across methods when closely related cell-states (i.e., the same cell type under different conditions) are simulated at varying proportions.

We will also build a multi-tissue autoimmune/rheumatic disease compendium from existing public data of bulk tissues.
We have curated a relevant set of datasets to date of more than 12,000 samples.
This compendium includes samples from patients with systemic lupus erythematosus (SLE), sarcoidosis, and inflammatory bowel disorders among many other diseases.
This compendium will have desirable properties for evaluating single-cell data-derived signatures, as it will allow us to evaluate cell type proportions and phenotypes based on a body of previous literature.
For instance, we expect to detect higher proportions of activated macrophages in lupus nephritis samples as compared to controls [@doi:10.4049/jimmunol.1103031].
In addition, the compendium includes drug studies of highly-targeted therapeutics (e.g., a monoclonal antibody to IFN-gamma [anti-IFNg] in the context of systemic lupus erythematosus [@doi:10.1002/art.39248]).
In the anti-IFNg example above, we cam ask what cell-types change in proportion during the reduction of this cytokine or if the IFN-inducible expression is preferentially altered in one cell type.
Such experiments allow us to not only assess methods for defining cell-type from single-cell data, but also methods that may be useful in decomposing whole-tissue bulk data that are developed through this initiative such as latent space arithmetic (see aim 3, below).

Importantly, a bulk compendium comprised of public data from a disease context will enable HCA participants (methods developers, RNA-seq based assay developers, and others) to directly evaluate approaches in the context where we expect their most immediate impact: application to existing datasets to explain disease-relevant phenomena via a single-cell perspective.
