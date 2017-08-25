### Aim 3: Evaluate the extent to which models that create reduced representations enable latent space arithmetic in the HCA.

Latent features describe a lower dimensional representation of input data and there are many dimensionality reduction methods, including principal components analysis (PCA) and non-negative matrix factorization (NMF).
Certain classes of generative deep neural network models, including VAEs and GANs, add a distribution matching constraint that can imbue intuitive mathematical features to the learned latent features [1,2].
For instance, a GAN learned latent features that could be manipulated with arithmetic: Subtracting out the essence of a smile from a woman and adding it to a neutral man resulted in an image vector of a smiling man [3].
We will assess the extent to which single cell gene expression latent spaces aggregated from deep generative neural networks can be manipulated mathematically.
The _objective of this aim_ is to test the hypothesis that arithmetic in the latent space will accurately predict gene expression values of perturbation experiments in specific cell types.


We will construct latent spaces with linear (PCA, NMF) and nonlinear methods (ADAGE, deep VAE models from Aim 1) as well as other relevant methods from groups funded under this initiative.
We will test latent space arithmetic using benchmark data that contains parallel perturbations, which allows us to hold out one or more perturbations for one or more cell types to evaluate each method.
We describe two experiments using data proposed by Arjun Raj's group, but any HCA benchmark datasets with similar properties will be suitable.

The differentiation of cardiomyocytes from fibroblasts is well studied, and the driving transcription factors have been identified [4], and the Raj lab proposes to generate data related to this process.
We will use these data to test the hypothesis that latent space vectors capture key transcription factor networks (Gata4, Mef2c, and Tbx5) known to drive the transition.
First, we will obtain the latent space projections of isolated human fibroblasts and cardiomyocytes.
Subtracting the mean latent space representation vector of fibroblasts from the cardiomyocyte latent vector will provide a differentiation vector.
We will compare the gene composition of this differentiation vector to targets calls from cistrome [@doi:10.1093/nar/gkw983], which are available for each of these TFs.

Latent space arithmetic can also be used to generate new hypothetical data, which we will test via gene expression assays of homogenized cell-types under various perturbations.
We will hold out one or more cell types under a specific perturbation, while mapping other observations into the latent space.
We can calculate the perturbation vector by subtracting the latent space vectors of included cell types from the same cell types after perturbation.
Adding the perturbation vector to one of the cell types for which the perturbation was withheld allows us to generate hypothetical data.
We will compare the results to the held out data to determine the extent to which measurements capture the cell type after perturbation.
We will compare results with each latent space to the baseline method of performing the same subtraction operations in raw gene expression space.

## References:

1.	Kingma, D. P. & Welling, M. Auto-Encoding Variational Bayes. ArXiv13126114 Cs Stat (2013).
2.	Goodfellow, I. J. et al. Generative Adversarial Networks. ArXiv14062661 Cs Stat (2014).
3.	Radford, A., Metz, L. & Chintala, S. Unsupervised Representation Learning with Deep Convolutional Generative Adversarial Networks. ArXiv151106434 Cs (2015).
4.  Ieda, M. et al. Direct Reprogramming of Fibroblasts into Functional Cardiomyocytes by Defined Factors. Cell 142, 375–386 (2010).
