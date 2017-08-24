## Summary

A latent space is a low-dimensional representation of meaningful sources of variability in a dataset.
Latent spaces can be sampled to provide meaningful interpolations between examples.
Latent spaces of image datasets can be used to transition a picture of a young man to a female representation or to an older version of the subject.
Biologically relevant latent spaces provide an opportunity for hypothesis generation on a genome-wide scale.
For example, we could imagine that we might predict how gene expression would change over time, with drug treatment, in the context of a specific disease, or a combination of these factors.
Certain deep learning methods are capable of constructing meaningful latent spaces in other fields.
Our objective is to determine whether or not these methods can be trained with single cell expression data and the extent to which such methods define biologically-relevant latent spaces.

## Aims

**Aim 1: Develop proof-of-concept unsupervised deep learning methods for single cell transcriptomic data from the HCA.**

**Aim 2: Generate a benchmark dataset of harmonized public data to evaluate the extent to which HCA cell types capture disease states for multiple rheumatic diseases.**

**Aim 3: Evaluate the extent to which models that create reduced representations enable latent space arithmetic in the HCA.**

## Prior Contributions / Preliminary Results (not required)

* ADAGE/eADAGE
* GANs + clinical data
* Variational Autoencoders

## Proposed work and deliverables

* Procedure for training VAEs (backup option: GANs)
    * Bulk version -> single cell
    * Augmented training evaluated.
* Application of deep NN -> HCA data.
* Evaluation of models by application to rheumatic disease compendium
* Evaluation of methods with regards to HCA benchmark datasets (TBD)

## Proposal for evaluation and dissemination.

* Biological grounding in the context of rheumatic disease (Mike Collab Network?)
* Benchmark HCA datasets (ideal: many cell pops +/- some perturbations).
    * Computationally mask one or more cell types from one or more perturbations -> eval extent to which latent space arithmetic produces observations that match real observations.

## Statement of commitment to share

We look forward to sharing our proposals, methods, data, and code with other researchers funded by this RFA, with CZI, and with the broader research community.

* We are already sharing our proposal under a CC-BY license as it is being written.
* We will share our description of methods, written reports, figures, and other such outputs under a CC-BY license.
* We will share any resulting data under a CC0 license to make our intention to share clear, though it is worth noting that in the US much of this data would not seem to be copyrightable.
* We will share source code for our software under the BSD 3-clause license. If we are contributing changes to software with a different open license, we will use that license. To avoid a situation where our work cannot be shard, we will not contribute to software with any license that has not been approved by the OSI as an open source license.

We would consider sharing any of these elements under more permissive licenses if requested by CZI.
