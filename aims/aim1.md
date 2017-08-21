# Aim 1:

**Develop proof-of-concept deep learning methods to learn a latent space from single cell transcriptomic data from the HCA.**

## Background

Deep generative models can approximate a lower dimensional manifold often represented in a learned latent space.
The latent space of such unsupervised models, including Restricted Boltzmann Machines (RBMs), Variational Autoencoders (VAEs), and Generative Adversarial Networks (GANs), preserves non-linear distances between objects and represents features driving data variation [1,2,3].
Thus, the algorithms serve a dual purpose: (1) dimensionality reduction and (2) ability to generate fake data.
Applied to image data, the algorithms have shown much promise in generating realistic looking fake data from lower dimensional features that shed light on the underlying structure of a complex, high dimensional arrangement of pixels.
Similarly, as applied to single cell gene expression data, such algorithms can learn complex biological structures enabling cell-type segmentation and interpolations between cell states.
We have previously shown the ability to train a VAE on bulk gene expression data that identifies latent features with biological meaning [4].
Additionally, others trained a GAN on images of yeast cells and interpolate the latent space to model protein localization [5].
Therefore, the _objective of this aim_ is to train a generative model on single cell gene expression features from the HCA demonstrating a proof of concept modeling approach.

## Approach

We will train a VAE on single cell transcriptome data as provided by the HCA.
An example of such data would be provided by Arjun Raj.
These data consist of RNAseq of homogenized cell-types under various perturbations.
Much like the augmented training of image data under arbitrary rotations and transformations, training a VAE on cell-types under multiple perturbations will learn invariant features defining cell-type.

In order to select the best model to validate, we will perform an extensive grid search with a held-out validation set to identify optimal VAE architecture and hyperparameters.
Once an optimal model is identified, we will demonstrate proof of concept by validating learned features.
The largest features in single cell gene expression are the abundance of technical and biological zeroes, and we will assess the ability of the VAE to automatically extract this artifact by comparing it to existing batch correction methods [6].

Additionally, we will use the learned VAE features to classify specific perturbations.
By inputing the learned features into a supervised algorithm, we can interrogate how well the algorithm learned features to separate condition.
We will compare performance of the classification task with features derived from raw data and other dimensionality reduction techniques.

## References:

1.  Smolensky, P. Information processing in dynamical systems: Foundations of harmony theory. In D. Rumelhart and J. McClelland
(Eds.), Parallel distributed processing, vol. 1, chapter 6, 194–281. Cambridge: MIT Press (1986).
2.	Kingma, D. P. & Welling, M. Auto-Encoding Variational Bayes. ArXiv13126114 Cs Stat (2013).
3.	Goodfellow, I. J. et al. Generative Adversarial Networks. ArXiv14062661 Cs Stat (2014).
4.	Way, G. P., & Greene, C.S. Extracting a Biologically Relevant Latent Space from Cancer Transcriptomes with Variational Autoencoders. bioRxiv. https://doi.org/10.1101/174474 (2017).
5.  Osokin, A., Chessel, A., Carazo Salas, R. E., & Vaggi, F. GANs for Biological Image Synthesis. ArXiv170804692 Cs Stat (2017).
6.  Bacher, R. & Kendziorski, C. Design and computational analysis of single-cell RNA-sequencing experiments. Genome Biology, 17:63 (2016).

