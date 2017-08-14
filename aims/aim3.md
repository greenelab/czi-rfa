# Aim 3:

**To evaluate the extent to which latent space arithmetic works in HCA benchmark datasets generated under this or other initiatives.**

## Background:

Generative models can approximate a lower dimensional manifold.
These manifolds constrain datasets to occupy a small subset of total space and can contain vital information regarding data generating functions.
Remarkably, certain classes of generative models, including variational autoencoders (VAE) and generative adversarial networks (GANs), that add a distribution matching constraint, imbue intuitive mathematical features to the learned latent features [1,2].
For instance, a recent GAN learned latent representations that could be manipulated with arithmetic: Subtracting out the essence of a smile from a woman and adding it to a neutral man resulted in an image vector of a smiling man [3].
While image data have been extensively modeled, less is known about the learning capacity of such models applied to single cell gene expression data.
The _objective of this aim_ is to assess the extent to which single cell gene expression latent spaces can be manipulated and how such manifold preserving properties will enable biological discovery regarding cell type, cellular differentiation, and disease perturbation.
The _working hypothesis of this aim_ is that the latent space will preserve arithmetic operations between state transitions of differentiating hematopoietic cells.

## Preliminary Data:

We have previously modeled bulk gene expression data with denoising autoencoders and non-linear activation functions [4–6].
These models have engineered features that encode information about gene ontology terms and KEGG pathways.
We have shown that these models capture more information than linear models like PCA, ICA, and NMF, and that the learned features are often overlapping and redundant [7].
We have also recently trained a VAE on cancer gene expression datasets and observed that the learned latent features capture interesting biological distributions of features across different tumors [8].
Therefore, we are confident that modeling gene expression with reconstruction cost captures meaningful biology.
Furthermore, we trained VAEs and GANs on electronic health record (EHR) data including a model to generate privacy aware fake datasets to supplement learning or data sharing [9,10].
We will develop similar approaches to extract knowledge from Human Cell Atlas (HCA) single cell data.

## Enabling latent space arithmetic in single cell gene expression with variational autoencoders:

Once a VAE is trained on HCA single cell gene expression data, we will obtain the latent space representation of a heterogeneous population of blood cells.
We will use well characterized gene expression markers and unsupervised clustering algorithms to define core subsets of common blood cell types.
We will define core subsets by positive silhouette widths of well characterized cell types representing a variety of the differentiation trajectories of hematopoiesis.
We will test our hypothesis that the learned latent space vectors can be manipulated mathematically to obtain biological insights and useful properties.
We will calculate the mean latent space vector representing erythroid cells and mean vector representing the population of stem cell progenitors.
We will take the mean of these two vectors, which, if the hypothesis were true, will represent cells in intermediate trajectories.
We will repeat this procedure with mature leukocyte populations.
We will determine how similar the observed mean vectors for observed intermediate cell-types are to the predicted cell types.
As a negative control, we will perform the same arithmetic using raw gene expression features, PCA features, and vanilla autoencoder features
A manifold preserving algorithm, a VAE can capture underlying trajectories between cell states, and can help to define core cell types.

## Potential Pitfalls and Alternative Approaches:

-	Latent space does not capture manifold
    - Compare differentiation trajectory to current state of the art pseudotime algorithms
-	Cell types not defined
    -	Use single cell markers as identified in Villani et al. 2017 and conservative unsupervised clustering algorithms to define core subset [11]
-	Adjusting for Large batch effects
    -	Variance captured in data should be explained by VAE feature(s)
    - Identify highest activated feature separating batch, set activation nodes constant for each cell, and test reconstruction pre and post correction

## References:

1.	Kingma, D. P. & Welling, M. Auto-Encoding Variational Bayes. ArXiv13126114 Cs Stat (2013).
2.	Goodfellow, I. J. et al. Generative Adversarial Networks. ArXiv14062661 Cs Stat (2014).
3.	Radford, A., Metz, L. & Chintala, S. Unsupervised Representation Learning with Deep Convolutional Generative Adversarial Networks. ArXiv151106434 Cs (2015).
4.	Tan, J., Ung, M., Cheng, C. & Greene, C. S. Unsupervised feature construction and knowledge extraction from genome-wide assays of breast cancer with denoising autoencoders. Pac. Symp. Biocomput. Pac. Symp. Biocomput. 132–143 (2015).
5.	Tan, J. et al. Unsupervised Extraction of Stable Expression Signatures from Public Compendia with an Ensemble of Neural Networks. Cell Syst. 5, 63–71.e6 (2017).
6.	Tan, J., Hammond, J. H., Hogan, D. A. & Greene, C. S. ADAGE-Based Integration of Publicly Available Pseudomonas aeruginosa Gene Expression Data with Denoising Autoencoders Illuminates Microbe-Host Interactions. mSystems 1, e00025-15 (2016).
7.	Chen, K. M. et al. PathCORE: Visualizing globally co-occurring pathways in large transcriptomic compendia. (2017). doi:10.1101/147645
8.	Way, G. P. & Greene, C. S. Extracting a Biologically Relevant Latent Space from Cancer Transcriptomes with Variational Autoencoders. bioRxiv 174474 (2017). doi:10.1101/174474
9.	Beaulieu-Jones, B. K., Greene, C. S. & Pooled Resource Open-Access ALS Clinical Trials Consortium. Semi-supervised learning of the electronic health record for phenotype stratification. J. Biomed. Inform. 64, 168–178 (2016).
10.	Beaulieu-Jones, B. K., Wu, Z. S., Williams, C. & Greene, C. S. Privacy-preserving generative deep neural networks support clinical data sharing. (2017). doi:10.1101/159756
11.	Villani, A.-C. et al. Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors. Science 356, eaah4573 (2017).
