# Aim 1:

**Develop proof-of-concept deep learning methods to learn a latent space from single cell transcriptomic data.**

## Background

Deep unsupervised models can learn to approximate a data generating manifold often represented in a lower dimensional latent space.
The latent space of such models, including Restricted Boltzmann Machines (RBMs), Variational Autoencoders (VAEs), and Generative Adversarial Networks (GANs), preserves non-linear distances between objects and learns to represent features driving variation in the input data [1,2,3].
Thus, the algorithms serve a dual purpose: (1) non-linear dimensionality reduction and (2) the ability to generate hypothetical data.
Applied to image data, the algorithms have shown much promise in generating realistic looking fake data and learning a lower dimensional representation that captures the underlying structure of high dimensional arrangements of pixels.
Similarly, as applied to single cell gene expression data, such algorithms can learn complex biological features that enable cell-type segmentation and cell state interpolation.
Previously, we have demonstrated the ability to train a VAE on bulk gene expression data to identify latent features with biological meaning [4].
Others trained GANs on images of yeast cells and interpolated the learned latent space to model protein localization [5].
Therefore, the _objective of this aim_ is to train a generative model on single cell gene expression features from the HCA demonstrating a proof of concept modeling approach.

## Approach

We will train a VAE on single cell transcriptome data as provided by the HCA.
An example of such data would be provided by the laboratory of Arjun Raj.
These data are benchmark RNAseq of homogenized cell-types under various exogenous perturbations.
Much like the augmented training of image data under arbitrary rotations and transformations, training a VAE on cell-types under multiple perturbations will learn invariant features defining cell-type.
Aiding in our search for invariant features, we will use a fast bootstrap resampling approach developed in the laboratory of Rob Patro.
Briefly, the approach will resample RNAseq reads from the same sample under different initialized parameters artificially increasing sample size while retaining cell-type features.

In order to select the best model to validate, we will perform an extensive grid search with a held-out validation set to identify the optimal VAE architecture and hyperparameters.
We will select optimal models based on minimizing reconstruction loss plus a KL divergence term which constrains the features to a Gaussian distribution.
Once an optimal model is identified, we will demonstrate proof of concept by validating learned features.
The largest features in single cell gene expression are the abundance of technical and biological zeroes, and we will assess the ability of the VAE to automatically extract this artifact by comparing it to existing batch correction methods [6].

Additionally, we will use the learned VAE features to classify specific perturbations.
By inputing the learned latent features into a supervised algorithm, we can interrogate how well the algorithm's features separate condition.
We will compare performance of the classification task with features derived from raw data and other dimensionality reduction techniques.
Once trained, we will distribute the reproducible algorithm for use by other groups in the CZI collective.

## References:

1.  Smolensky, P. Information processing in dynamical systems: Foundations of harmony theory. In D. Rumelhart and J. McClelland
(Eds.), Parallel distributed processing, vol. 1, chapter 6, 194–281. Cambridge: MIT Press (1986).
2.	Kingma, D. P. & Welling, M. Auto-Encoding Variational Bayes. ArXiv13126114 Cs Stat (2013).
3.	Goodfellow, I. J. et al. Generative Adversarial Networks. ArXiv14062661 Cs Stat (2014).
4.	Way, G. P., & Greene, C.S. Extracting a Biologically Relevant Latent Space from Cancer Transcriptomes with Variational Autoencoders. bioRxiv. https://doi.org/10.1101/174474 (2017).
5.  Osokin, A., Chessel, A., Carazo Salas, R. E., & Vaggi, F. GANs for Biological Image Synthesis. ArXiv170804692 Cs Stat (2017).
6.  Bacher, R. & Kendziorski, C. Design and computational analysis of single-cell RNA-sequencing experiments. Genome Biology, 17:63 (2016).

