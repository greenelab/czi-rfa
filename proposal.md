# Title: Genome-wide hypothesis generation for single-cell expression via latent spaces of deep neural networks

## Primary Focus Area

Computational Biology

## Project Summary

## Keywords

gene expression, deep learning, unsupervised, latent space, hypothesis generating

## Full citations of contributions relevant to proposal

* Tan J, Doing G, Lewis KA, Price CE, Chen KM, Cady KC, Perchuk B, Laub MT, Hogan DA, Greene CS. Unsupervised extraction of stable expression signatures from public compendia with an ensemble of neural networks. _Cell Systems_. 5:63-71. PMID:[28711280](https://www.ncbi.nlm.nih.gov/pubmed/28711280)
* Beaulieu-Jones BK, Wu ZS, Williams C, Greene CS. Privacy-preserving generative deep neural networks support clinical data sharing. _bioRxiv_. DOI: [0.1101/159756](http://doi.org/10.1101/159756)
* Ching T, Himmelstein DS, Beaulieu-Jones BK, Kalinin AA, Do BT, Way GP, Ferro E, Agapow P, Xie W, Rosen GL, Legenrich BJ, Israeli J, Lanchantin J, Woloszynek S, Carpenter AE, Shrikumar A, Xu J, Cofer EM, Harris DJ, DeCaprio D, Qi Y, Kundaje A, Peng Y, Wiley LK, Segler MHS, Gitter A+, Greene CS+. Opportunities and Obstacles for Deep Learning in Biology and Medicine. _bioRxiv_. doi:[10.1101/142760](https://doi.org/10.1101/142760)
* Way GP, Greene CS. Extracting a Biologically Relevant Latent Space from Cancer Transcriptomes with Variational Autoencoders. _bioRxiv_. doi:[10.1101/174474](https://doi.org/10.1101/174474)
* Beaulieu-Jones BK, Greene CS. Reproducibility of computational workflows is automated using continuous analysis. _Nature Biotechnology_. 35:342-346. PMID:[28288103](https://www.ncbi.nlm.nih.gov/pubmed/28288103)

## Collaborative Network

Our collaborative network includes a number of groups responding to this RFA. This group of investigators came together around a set of responses to the RFA that were shared openly on GitHub as they were written. Selected interactions between groups center on topics such as:

* _Deep learning_ methods are expected to be relevant to many proposals. For example, Arjun Raj proposes to apply these methods to analyze images of single cells. We recently completed a large GitHub-authored collaborative review to identify deep learning strategies that address biological problems. We look forward to discussing techniques with Arjun's group and others as they apply these methods.
* _Augmented training_ is a technique used to analyze image data with deep learning. During training, arbitrary patches are selected and arbitrary rotations are applied to training data to avoid learning irrelevant patterns. We propose an innovative equivalent for single cell expression data by sampling reads. This requires fast expression quantification. We plan to work with Rob Patro, an expert in fast gene expression quantification, to implement data augmentation by very rapidly sampling single cell observations.
* _Benchmark data_ need to be diverse and different data will be important for different evaluations. We will contribute a large compendium of harmonized bulk gene expression data from multiple rheumatic diseases. This, combined with HCA single cell assays, will help us to evaluate the extent to which bulk public data measuring disease states can be explained by HCA-identified cell types. Our analysis of latent spaces benefits from datasets that measure parallel perturbations across multiple cell types, like those proposed by Arjun Raj. Such data from the Raj lab and other groups will be needed to effectively evaluate unsupervised deep learning methods for single cell data.
* _Model interpretation_ is a key factor for any unsupervised approach. Elana Fertig and Loyal Goff's applications describe methods for the interpretation of unsupervised models. We are excited to collaborate on methods and infrastructure for model interpretation. Members of the Fertig lab have requested features in other software that we have developed, which are now implemented. We look forward to continued interactions as a collaborative network.
* _Interactive servers_ can empower bench scientists to direct the analysis of single-cell genomics data. Lana Garmire's group developed Granatum, an interactive graphic interface for single cell analysis and is responding to this analyses. We will provide source code for generating VAEs that can be incorporated into Granatum and Granatum-like systems.
* _Reproducible workflows_ will be needed to consistently evaluate methods. For example, Rob Patro's group has previously extended our example configurations to provide a workflow for Salmon. We are happy to work with any participating groups and CZI software developers to provide and configure continuous integration servers for scientific workflows. We have not budgeted this in our application as a centrally provided service makes sense; however, we can modify our budget to cover this if requested.

## List of key personnel including name, organization, role on project

* Casey Greene (PI)
* Qiwen Hu (Postdoc)
* Jaclyn Taroni (Postdoc)
* Gregory Way (Student)

## Proposal

### Summary

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
The _rationale_ is that latent space arithmetic would enable researchers to generate hypothetical cells under simulated perturbation for rapid discovery-oriented analyses probing cell-type, cellular differentiation.

### Aims

**Aim 1: Develop proof-of-concept unsupervised deep learning methods for single cell transcriptomic data from the HCA.**

> Supporting analytical methods and machine learning approaches to solving problems such as multimodal integration, inference of state transitions and developmental trajectories, and representation of spatial relationships at the cellular or molecular level

**Aim 2: Generate a benchmark dataset of harmonized public data to evaluate the extent to which HCA cell types capture disease states for multiple rheumatic diseases.**

> Generating curated benchmark datasets from new or existing data for evaluating computational methods and designing future analysis competitions

**Aim 3: Evaluate the extent to which models that create reduced representations enable latent space arithmetic in the HCA.**

> Supporting analytical methods and machine learning approaches to solving problems such as multimodal integration, inference of state transitions and developmental trajectories, and representation of spatial relationships at the cellular or molecular level

### Prior Contributions / Preliminary Results (not required)

* ADAGE/eADAGE
* GANs + clinical data
* Variational Autoencoders

### Proposed work and deliverables

#### Aim 1: Develop proof-of-concept unsupervised deep learning methods for single cell transcriptomic data from the HCA.

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

##### References:

1.  Smolensky, P. Information processing in dynamical systems: Foundations of harmony theory. In D. Rumelhart and J. McClelland
(Eds.), Parallel distributed processing, vol. 1, chapter 6, 194–281. Cambridge: MIT Press (1986).
2.  Kingma, D. P. & Welling, M. Auto-Encoding Variational Bayes. ArXiv13126114 Cs Stat (2013).
3.  Goodfellow, I. J. et al. Generative Adversarial Networks. ArXiv14062661 Cs Stat (2014).
4.  Gregor, K., Danihelka, I., Graves, A., Rezende, D. J., Wierstra, D. DRAW: A Recurrent Neural Network for Image Generation. ArXiv150204623 Cs Stat (2015).
5.  Osokin, A., Chessel, A., Carazo Salas, R. E., & Vaggi, F. GANs for Biological Image Synthesis. ArXiv170804692 Cs Stat (2017).
6.  FaceApp - Free Neural Face Transformation Filters. Accessed 2017-08-24. https://www.faceapp.com
7.  Way, G. P., & Greene, C.S. Extracting a Biologically Relevant Latent Space from Cancer Transcriptomes with Variational Autoencoders. bioRxiv. https://doi.org/10.1101/174474 (2017).
8.  Liao, Q., Leibo, J. Z., Poggio, T. Learning invariant representations and applications to face verification. Advances in Neural Information Processing Systems 26 (2013).

#### Aim 2: Evaluate the extent to which low-dimensional representations of single cell data learned by various methods can be used to decompose bulk tissues in the context of rheumatic diseases.**

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

#### Aim 3: Evaluate the extent to which models that create reduced representations enable latent space arithmetic in the HCA.

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

##### References:

1.	Kingma, D. P. & Welling, M. Auto-Encoding Variational Bayes. ArXiv13126114 Cs Stat (2013).
2.	Goodfellow, I. J. et al. Generative Adversarial Networks. ArXiv14062661 Cs Stat (2014).
3.	Radford, A., Metz, L. & Chintala, S. Unsupervised Representation Learning with Deep Convolutional Generative Adversarial Networks. ArXiv151106434 Cs (2015).
4.  Ieda, M. et al. Direct Reprogramming of Fibroblasts into Functional Cardiomyocytes by Defined Factors. Cell 142, 375–386 (2010).

### Proposal for evaluation and dissemination.

* Biological grounding in the context of rheumatic disease
* Benchmark HCA datasets (ideal: many cell pops +/- some perturbations).
    * Computationally mask one or more cell types from one or more perturbations -> eval extent to which latent space arithmetic produces observations that match real observations.

### Statement of commitment to share

We look forward to sharing our proposals, methods, data, and code with other researchers funded by this RFA, with CZI, and with the broader research community.

* We are already sharing our proposal under a CC-BY license as it is being written.
* We will share our description of methods, written reports, figures, and other such outputs under a CC-BY license.
* We will share any resulting data under a CC0 license to make our intention to share clear, though it is worth noting that in the US much of this data would not seem to be copyrightable.
* We will share source code for our software under the BSD 3-clause license. If we are contributing changes to software with a different open license, we will use that license. To avoid a situation where our work cannot be shard, we will not contribute to software with any license that has not been approved by the OSI as an open source license.

We would consider sharing any of these elements under more permissive licenses if requested by CZI.

## References (non-doi)
1.  Smolensky, P. Information processing in dynamical systems: Foundations of harmony theory. In D. Rumelhart and J. McClelland
(Eds.), Parallel distributed processing, vol. 1, chapter 6, 194–281. Cambridge: MIT Press (1986).
4.  Gregor, K., Danihelka, I., Graves, A., Rezende, D. J., Wierstra, D. DRAW: A Recurrent Neural Network for Image Generation. ArXiv150204623 Cs Stat (2015).
5.  Osokin, A., Chessel, A., Carazo Salas, R. E., & Vaggi, F. GANs for Biological Image Synthesis. ArXiv170804692 Cs Stat (2017).
6.  FaceApp - Free Neural Face Transformation Filters. Accessed 2017-08-24. https://www.faceapp.com
7.  Way, G. P., & Greene, C.S. Extracting a Biologically Relevant Latent Space from Cancer Transcriptomes with Variational Autoencoders. bioRxiv. https://doi.org/10.1101/174474 (2017).
8.  Liao, Q., Leibo, J. Z., Poggio, T. Learning invariant representations and applications to face verification. Advances in Neural Information Processing Systems 26 (2013).

## Biosketch (Attached)

**Casey Greene, PhD (PI).**
_1.2 Cal Months_.
    * Short bio.
    * cgreene on github.
    * Dr. Greene has primary responsibility for managing and carrying out the proposed research, interpreting findings, disseminating findings, and preparing and submitting progress reports.

## Description of other key personnel

**Gregory Way, MS (GCB Graduate Student).**
_4.2 Cal Months_.
    * Short bio.
    * gwaygenomics on github.
    * Also partially supported by T32 for latent space methods for cancer datasets (thus only 35%)
    * Adaptation of deep learning methods -> single cell (w/ Qiwen)

**Qiwen Hu, PhD (Postdoctoral Researcher).**
_6.0 Cal Months_.
    * Short bio
    * huqiwen0313 on github.
    * Adaptation of deep learning methods -> single cell (w/ Greg)
    * Evaluation of deep learning methods on HCA data.

**Jaclyn Taroni, PhD (Postdoctoral Researcher).**
_2.4 Cal Months_.
    * Short bio.
    * jaclyn-taroni on github.
    * Evaluation for bio relevance to rheumatic disease bulk dataset compendia

## Brief preliminary budget (Attached)
