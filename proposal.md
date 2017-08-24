## Summary

Certain types of deep unsupervised neural networks, including Restricted Boltzmann Machines (RBMs) [1], Variational Autoencoders (VAEs) [@arxiv:1312.6114], and Generative Adversarial Networks (GANs) [@arxiv:1406.2661], can generate hypothetical data by learning and decoding a lower dimensional latent space.
This latent space preserves non-linear distances between data points and represents features that drive variation in the input data.
An ideal latent space enables arithmetic operations that use data to produce novel transformations.
Models of image data, such as the "DRAW" network, produce realistic looking fake images through the learning of a lower dimensional representation of different pixel arrangements [4].
GANs trained on fluorescent microscopy images of yeast cells enable interpolation within the learned latent space to model cell-cycle protein localization [5].
A fun and intuitive example of this is FaceApp, which is capable of modifying a picture of an individual to produce an image of the subject at an older age, with a different expression, or of a different gender [6].

Generative deep neural network models that produce meaningful latent spaces can serve two purposes: (1) non-linear dimensionality reduction and (2) the ability to generate hypothetical data.
Biologically relevant latent spaces provide an opportunity for hypothesis generation on a genome-wide scale.
For example, we could imagine that we might predict how gene expression would change for each cell type measured in the human cell atlas with drug treatment, in the context of a specific genetic variant, with a specific disease, or as a combination of these and other factors.

Our _overall objective_ is to determine whether or not unsupervised deep neural network models can be trained with single cell expression data and the extent to which such methods define biological latent spaces that capture disease states and other perturbations.

## Aims

**Aim 1: Develop proof-of-concept deep learning methods to learn a latent space from single cell transcriptomic data from the HCA.**

**Aim 2: Evaluate the extent to which low-dimensional representations of single cell data learned by various methods can be used to decompose bulk tissues in the context of rheumatic diseases.**

**Aim 3: Evaluate the extent to which latent space arithmetic works in HCA benchmark datasets generated under this or other initiatives.**

## Prior Contributions / Preliminary Results (not required)

* ADAGE/eADAGE
* GANs + clinical data
* Variational Autoencoders
* Scleroderma pubs -> cell lineage importance

## Proposed work and deliverables

* Procedure for training VAEs (backup option: GANs)
    * Bulk version -> single cell
    * Augmented training evaluated.
* Application of deep NN -> HCA data.
* Evaluation of models by application to rheumatic disease compendium
* Evaluation of methods with regards to HCA benchmark datasets (TBD)

## Proposal for evaluation and dissemination.

* Biological grounding in the context of rheumatic disease
* Benchmark HCA datasets (ideal: many cell pops +/- some perturbations).
    * Computationally mask one or more cell types from one or more perturbations -> eval extent to which latent space arithmetic produces observations that match real observations.

## Statement of commitment to share

We look forward to sharing our proposals, methods, data, and code with other researchers funded by this RFA, with CZI, and with the broader research community.

* We are already sharing our proposal under a CC-BY license as it is being written.
* We will share our description of methods, written reports, figures, and other such outputs under a CC-BY license.
* We will share any resulting data under a CC0 license to make our intention to share clear, though it is worth noting that in the US much of this data would not seem to be copyrightable.
* We will share source code for our software under the BSD 3-clause license. If we are contributing changes to software with a different open license, we will use that license. To avoid a situation where our work cannot be shard, we will not contribute to software with any license that has not been approved by the OSI as an open source license.

We would consider sharing any of these elements under more permissive licenses if requested by CZI.

# References (non-doi)
1.  Smolensky, P. Information processing in dynamical systems: Foundations of harmony theory. In D. Rumelhart and J. McClelland
(Eds.), Parallel distributed processing, vol. 1, chapter 6, 194–281. Cambridge: MIT Press (1986).
4.  Gregor, K., Danihelka, I., Graves, A., Rezende, D. J., Wierstra, D. DRAW: A Recurrent Neural Network for Image Generation. ArXiv150204623 Cs Stat (2015).
5.  Osokin, A., Chessel, A., Carazo Salas, R. E., & Vaggi, F. GANs for Biological Image Synthesis. ArXiv170804692 Cs Stat (2017).
6.  FaceApp - Free Neural Face Transformation Filters. Accessed 2017-08-24. https://www.faceapp.com
7.  Way, G. P., & Greene, C.S. Extracting a Biologically Relevant Latent Space from Cancer Transcriptomes with Variational Autoencoders. bioRxiv. https://doi.org/10.1101/174474 (2017).
8.  Liao, Q., Leibo, J. Z., Poggio, T. Learning invariant representations and applications to face verification. Advances in Neural Information Processing Systems 26 (2013).
