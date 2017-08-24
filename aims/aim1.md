### Aim 1: Develop proof-of-concept unsupervised deep learning methods for single cell transcriptomic data from the HCA.

Previously, we have demonstrated the ability to train a VAE on bulk gene expression data to identify latent biological features [7].
The _objective of this aim_ is to implement and test approaches to adapt deep generative models, such as VAEs or GANs, to HCA-produced single cell gene expression data.

Single cell data pose unique challenges, as well as opportunities, for deep neural network algorithms.
Though each estimate of transcript abundance in each cell may include substantially more error than bulk samples, there are also many observations as numerous cells are often assayed.
Also, the HCA project includes data from multiple cell types, so we expect substantial numbers of observations to be available.
From our experience with generative deep neural networks [@doi:10.1101/159756], we have found that it can be difficult to predict precisely which advances will enable robust training in a specific context, particularly when the specific context has not been tackled before.
We describe our standard approach here, as well as a couple selected techniques that we anticipate will be particularly helpful in this setting such as zero-inflated loss and data augmentation, but we expect to also employ other strategies where warranted.

Our standard approach is to first perform an extensive grid search to identify VAE architectures and hyperparameters that are amenable to this context.
Some parameter settings can be easily ruled out when models fail to train.
Chris Probert noted on our proposal's GitHub repository [@url:https://github.com/greenelab/czi-rfa/issues/11] that numerous manuscripts demonstrated the advantages of zero-inflated loss in the setting of single cell sequencing data [@doi:10.1186/s13059-015-0805-z @arxiv:1610.05857 @doi:10.1186/s13059-017-1188-0].
We will evaluate multiple types of reconstruction loss during training.
These efforts will allow us to identify a subset of parameter combinations that are worth exploring with new approaches designed specifically for this type of genomic data.

We will also develop a data augmentation approach for training deep models on single cell RNA-seq data.
These approaches have not yet been applied to genomic data, but are widely used in image analysis, where the goal is to distinguish relevant from irrelevant features in an image corpus.
To understand data augmentation, imagine pathology slides that have been scanned.
Each slide may be prepared and scanned in a subtly different orientation or at somewhat different magnifications depending on who generated the data.
A deep learning method may learn to identify these subtle differences, which is undesirable.
It's also possible that there are simply too few slides to train a deep learning algorithm.
To address these challenges, deep learning practitioners can apply arbitrary rotations, zooms, and other irrelevant transformations to image data during training.
This process is called data augmentation.

No data augmentation approaches for genomic data have been published to date.
However, we expect that very fast abundance estimates from RNA-seq data [@doi:10.1038/nmeth.4197 @10.1038/nbt.3519] will make a data augmentation approach feasible in this domain as well.
Resampling reads to generate abundance estimates can help to capture the uncertainty in the data, akin to arbitrary rotations.
Subsampling reads to generate abundance estimates can help to capture changes related to sequencing depth but unrelated to the biology, akin to arbitrary zooms.
Therefore, we plan to collaborate with Rob Patro's laboratory via our collaborative network to implement approaches including but not limited to rapid subsampling and bootstrapping to generate augmented training data.
We posit that genomic data augmentation will improve latent feature generalization by separating biological from technical features and increasing the effective sample size during training.

We will evaluate these models on data held out during training.
We will select high-quality models by choosing those that minimize both reconstruction loss and a KL divergence term which constrains the features to a Gaussian distribution [@arxiv:1312.6114].
Models that exhibit desirable properties will be evaluated in the context of a rheumatic disease compendium (aim 2) as well as their suitability for latent space arithmetic (aim 3).

#### References:

1.  Smolensky, P. Information processing in dynamical systems: Foundations of harmony theory. In D. Rumelhart and J. McClelland
(Eds.), Parallel distributed processing, vol. 1, chapter 6, 194–281. Cambridge: MIT Press (1986).
2.  Kingma, D. P. & Welling, M. Auto-Encoding Variational Bayes. ArXiv13126114 Cs Stat (2013).
3.  Goodfellow, I. J. et al. Generative Adversarial Networks. ArXiv14062661 Cs Stat (2014).
4.  Gregor, K., Danihelka, I., Graves, A., Rezende, D. J., Wierstra, D. DRAW: A Recurrent Neural Network for Image Generation. ArXiv150204623 Cs Stat (2015).
5.  Osokin, A., Chessel, A., Carazo Salas, R. E., & Vaggi, F. GANs for Biological Image Synthesis. ArXiv170804692 Cs Stat (2017).
6.  FaceApp - Free Neural Face Transformation Filters. Accessed 2017-08-24. https://www.faceapp.com
7.  Way, G. P., & Greene, C.S. Extracting a Biologically Relevant Latent Space from Cancer Transcriptomes with Variational Autoencoders. bioRxiv. https://doi.org/10.1101/174474 (2017).
8.  Liao, Q., Leibo, J. Z., Poggio, T. Learning invariant representations and applications to face verification. Advances in Neural Information Processing Systems 26 (2013).
