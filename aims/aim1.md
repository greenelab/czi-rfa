# Aim 1:

**Develop proof-of-concept deep learning methods to learn a latent space from single cell transcriptomic data.**

## Background

Deep unsupervised models can generate hypothetical data by learning and decoding a lower dimensional latent space.
The latent space of such models, including Restricted Boltzmann Machines (RBMs), Variational Autoencoders (VAEs), and Generative Adversarial Networks (GANs), preserves non-linear distances between data points and represents features that drive variation in the input data [1,2,3].
Thus, the algorithms serve a dual purpose: (1) non-linear dimensionality reduction and (2) the ability to interpolate between data.
Applied to image data, the algorithms, including the recurrent VAE DRAW network, have shown much promise in generating realistic looking fake data through the learning of a lower dimensional representation of different pixel arrangements [4]. 
Additionally, GANs trained on fluorescent microscopy images of yeast cells enable interpolation of the learned latent space to model cell-cycle protein localization [5].
As applied to single cell gene expression data, we hypothesize that such algorithms can learn complex biological features that enable cell-type segmentation and cell state interpolation.
Previously, we have demonstrated the ability to train a VAE on bulk gene expression data to identify latent biological features [6].
Therefore, the _objective of this aim_ is to train a deep neural network generative model built from single cell gene expression from the HCA demonstrating a proof of concept modeling approach.

## Approach

We will train a VAE on single cell transcriptome data from the HCA.
In order to select the best model to validate, we will perform an extensive grid search with a held-out validation set to identify optimal VAE architecture and hyperparameters.
We will select optimal models based on minimizing reconstruction loss plus a KL divergence term which constrains the features to a Gaussian distribution.

In addition, we will compare our optimized VAE to an alternative VAE we will train with an augmentation approach.
A common image processing goal is to disentangle invariant features from a corpus of images.
These invariant features are key to describing the essence of the object in the image, but the algorithms must learn how to extract them under various rotations, zoom, and other augmentations.
To encourage invariant feature learning, images can be resampled with artificial rotation and added back as input to training deep models [7].
Therefore, we will collaborate with Rob Patro's laboratory and pursue a similar approach applied to single cell gene expression: Rapid subsampling and bootstrapping of raw RNAseq reads to overcome transcript depth and rare transcript biases.
We posit that this augmented training approach will improve latent feature generalization not only by isolating specific cell-type invariance, but also by increasing sample size during training.

We will evaluate our models using different transcriptomic data from the HCA.
An example of such data would be provided by the laboratory of Arjun Raj.
These data are bulk RNAseq of homogenized cell-types under various exogenous perturbations.
With this data, we can ask how well our algorithms reconstruct held out, purified cell-type.
Under perturbation, if the models are invariant, they should recover the same purified signal.
Once trained, we will distribute all reproducible algorithms for use by other groups in the CZI collective.

## References:

1.  Smolensky, P. Information processing in dynamical systems: Foundations of harmony theory. In D. Rumelhart and J. McClelland
(Eds.), Parallel distributed processing, vol. 1, chapter 6, 194–281. Cambridge: MIT Press (1986).
2.  Kingma, D. P. & Welling, M. Auto-Encoding Variational Bayes. ArXiv13126114 Cs Stat (2013).
3.  Goodfellow, I. J. et al. Generative Adversarial Networks. ArXiv14062661 Cs Stat (2014).
4.  Gregor, K., Danihelka, I., Graves, A., Rezende, D. J., Wierstra, D. DRAW: A Recurrent Neural Network for Image Generation. ArXiv150204623 Cs Stat (2015).
5.  Osokin, A., Chessel, A., Carazo Salas, R. E., & Vaggi, F. GANs for Biological Image Synthesis. ArXiv170804692 Cs Stat (2017).
6.  Way, G. P., & Greene, C.S. Extracting a Biologically Relevant Latent Space from Cancer Transcriptomes with Variational Autoencoders. bioRxiv. https://doi.org/10.1101/174474 (2017).
7.  Liao, Q., Leibo, J. Z., Poggio, T. Learning invariant representations and applications to face verification. Advances in Neural Information Processing Systems 26 (2013).

