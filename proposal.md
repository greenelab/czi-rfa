# Title: Genome-wide hypothesis generation for single-cell expression via latent spaces of deep neural networks

## Primary Focus Area

Computational Biology

## Project Summary

The Human Cell Atlas (HCA) aims to provide a comprehensive map of all types of human cells.
Connecting that map to disease states, which will be key to the CZI's mission of curing or managing all diseases in the next eighty years, will require us to see how these cell types change during aging, during disease processes, or in the presence of drugs.
Ideally, we'd be able to apply a transformation to the HCA's reference map to predict and study these states.

Certain types of deep unsupervised neural networks can generate hypothetical data by learning and decoding a lower dimensional latent space.
An ideal latent space enables arithmetic operations that use data to produce realistic output for novel transformations.
For example, FaceApp [@url:https://www.faceapp.com] can modify a picture of an individual to produce an image of the subject at an older age, with a different expression, or of a different gender.

The overall objective of this proposal is to determine how unsupervised deep neural network models can best be trained on single cell expression data from the HCA and the extent to which such models define biological latent spaces that capture disease states and targeted perturbations.
The rationale is that latent space arithmetic for single cell transcriptomes would enable researchers to use predict how the expression of every gene would change in each HCA-identified cell type in numerous conditions including after drug treatment, in the context of a specific genetic variant, with a specific disease, or a combination of these and other factors.

## Keywords

gene expression, deep learning, unsupervised, latent space, hypothesis generating

## Full citations of contributions relevant to proposal

* Tan J, Doing G, Lewis KA, Price CE, Chen KM, Cady KC, Perchuk B, Laub MT, Hogan DA, Greene CS. Unsupervised extraction of stable expression signatures from public compendia with an ensemble of neural networks. _Cell Systems_. 5:63-71. PMID:[28711280](https://www.ncbi.nlm.nih.gov/pubmed/28711280)
* Beaulieu-Jones BK, Wu ZS, Williams C, Greene CS. Privacy-preserving generative deep neural networks support clinical data sharing. _bioRxiv_. doi:[0.1101/159756](http://doi.org/10.1101/159756)
* Ching T, Himmelstein DS, Beaulieu-Jones BK, Kalinin AA, Do BT, Way GP, Ferro E, Agapow P, Xie W, Rosen GL, Legenrich BJ, Israeli J, Lanchantin J, Woloszynek S, Carpenter AE, Shrikumar A, Xu J, Cofer EM, Harris DJ, DeCaprio D, Qi Y, Kundaje A, Peng Y, Wiley LK, Segler MHS, Gitter A+, Greene CS+. Opportunities and Obstacles for Deep Learning in Biology and Medicine. _bioRxiv_. doi:[10.1101/142760](https://doi.org/10.1101/142760)
* Way GP, Greene CS. Extracting a Biologically Relevant Latent Space from Cancer Transcriptomes with Variational Autoencoders. _bioRxiv_. doi:[10.1101/174474](https://doi.org/10.1101/174474)
* Beaulieu-Jones BK, Greene CS. Reproducibility of computational workflows is automated using continuous analysis. _Nature Biotechnology_. 35:342-346. PMID:[28288103](https://www.ncbi.nlm.nih.gov/pubmed/28288103)

## Collaborative Network

Our collaborative network came together around a set of responses to the RFA that were shared openly on GitHub as they were written. Selected interactions between our groups center on topics such as:

* _Deep learning_ methods are expected to be relevant to many proposals. We recently completed a large GitHub-authored collaborative review of identify deep learning strategies for biological problems. Arjun Raj proposes to apply deep learning to images of single cells. We would love to use this RFA as an opportunity to host public discussions of techniques to troubleshoot deep learning challenges as they arise.
* _Data augmentation_ is a technique used in deep learning analysis of images in which arbitrary patches are selected or arbitrary rotations are applied during training to avoid learning irrelevant patterns. We propose an innovative equivalent for single cell expression data: sampling reads. This requires fast expression quantification. We plan to work with Rob Patro, an expert in fast gene expression quantification, to implement data augmentation by very rapidly sampling single cell observations. This parallels the sub-sampling proposed by Elana Fertig, and we will incorporate results from her work into our augmentation strategy.
* _Benchmark data_ need to be diverse and different data will be important for different evaluations. We will contribute a large compendium of harmonized bulk gene expression data from multiple rheumatic diseases. This, combined with HCA single cell assays, will help us to evaluate the extent to which bulk public data measuring disease states can be explained by HCA-identified cell types. Our analysis of latent spaces benefits from datasets that measure parallel perturbations across multiple cell types, like those proposed by Arjun Raj. Such data from the Raj lab and other groups will be needed to effectively evaluate unsupervised deep learning methods for single cell data.
* _Model interpretation_ is key for any unsupervised learning method. Elana Fertig and Loyal Goff propose new interpretation techniques, which we are excited to collaborate on. We have collaborated before via githug: members of the Fertig lab have requested features in software from our group, which we implemented. We look forward to continued interactions as a collaborative network.
* _Interactive servers_ can empower bench scientists to direct the analysis of single-cell genomics data. Lana Garmire's group developed Granatum, an interactive graphic interface for single cell analysis and is responding to this RFA. We will provide source code for generating deep neural network models of single cell data that can be incorporated into Granatum and Granatum-like systems.
* _Reproducible workflows_ will be needed to consistently and fairly evaluate methods. Our continuous analysis approach provides a framework for automated and reproducible analyses. Rob Patro's group has previously extended our example configurations to provide a workflow for Salmon. We are happy to work with any participating groups and CZI software developers to provide and configure continuous integration servers for scientific workflows. We have not budgeted this in our application as a centrally provided service makes sense; however, we can modify our budget to cover this if requested.

## List of key personnel including name, organization, role on project

* Casey Greene (PI)
* Qiwen Hu (Postdoc)
* Jaclyn Taroni (Postdoc)
* Gregory Way (Graduate Student)

## Proposal

### Summary

Certain types of deep unsupervised neural networks can generate hypothetical data by learning and decoding a lower dimensional latent space.
An ideal latent space enables arithmetic operations that use data to produce realistic output for novel transformations.
Models of image data produce realistic looking fake images [@doi:1502.04623] and capture protein localization within the cell-cycle to enable interpolation [@doi:1708.04692].
FaceApp [@url:https://www.faceapp.com], which modifies a picture of an individual to produce an image of the subject at an older age, with a different expression, or of a different gender is an accessible example of latent space transformations.

Our _overall objective_ is to determine how unsupervised deep neural network models can best be trained on single cell expression data from the Human Cell Atlas (HCA) and the extent to which such models define biological latent spaces that capture disease states and targeted perturbations.
The _rationale_ is that latent space arithmetic for genomic data would enable researchers to use predict how the expression of every gene would change in each HCA-identified cell type after drug treatment, in the context of a specific genetic variant, with a specific disease, or a combination of these and other factors.

### Aims

**Aim 1: Develop proof-of-concept unsupervised deep learning methods for single cell transcriptomic data from the HCA.**

> Supporting analytical methods and machine learning approaches for solving the inference of state transitions and developmental trajectories

**Aim 2: Generate a benchmark dataset of harmonized public data to evaluate the extent to which HCA cell types capture disease states for multiple rheumatic diseases.**

> Generating curated benchmark datasets from existing data for evaluating computational methods and designing future analysis competitions

### Prior Contributions / Preliminary Results

, including Restricted Boltzmann Machines (RBMs) [@url:http://dl.acm.org/citation.cfm?id=104290], Variational Autoencoders (VAEs) [@arxiv:1312.6114], and Generative Adversarial Networks (GANs) [@arxiv:1406.2661],

We have developed a number of methods to integrate transcriptomic data compendia with neural networks [@doi:10.1128/mSystems.00025-15 @doi:10.1016/j.cels.2017.06.003 @doi:10.1101/156620].
We've recently moved to generative models, which in the field of image analysis define meaningful latent spaces.
We found that GANs, with some additional steps, can generate realistic shareable samples from individual clinical trials records under a differential privacy framework [@doi:10.1101/159756].
In recent work, we've constructed VAEs over bulk transcriptomic data with the goal of describing a biologically-relevant latent space [@doi:10.1101/174474].
In this work, we will apply these unsupervised deep learning methods to single cell transcriptomic data, while incorporating data augmentation approaches that are novel to genomics.
We also bring expertise automating reproducible workflows to the HCA community [@doi:10.1038/nbt.3780].

### Proposed work and deliverables

#### Aim 1: Develop proof-of-concept unsupervised deep learning methods for single cell transcriptomic data from the HCA.

The _objective of this aim_ is to implement and test approaches to adapt deep generative models, such as VAEs or GANs, to HCA-produced single cell gene expression data.

Single cell data pose unique challenges, as well as opportunities, for deep neural network algorithms.
Though each estimate of transcript abundance in each cell may include substantially more error than bulk samples, there are also many observations as numerous cells are often assayed.
Also, the HCA project includes data from multiple cell types, so we expect substantial numbers of observations to be available.
From our experience with generative deep neural networks [@doi:10.1101/159756], we have found that it can be difficult to predict precisely which advances will enable robust training in a specific context, particularly when the specific context has not been tackled before.
We describe our standard approach here, as well as a couple selected techniques that we anticipate will be particularly helpful in this setting such as zero-inflated loss and data augmentation, but we expect to also employ other strategies where warranted.

Our standard approach is to first perform an extensive grid search to identify VAE architectures and hyperparameters that are amenable to this context.
We will evaluate multiple types of reconstruction loss during training, as Chris Probert noted on our proposal's GitHub repository [@url:https://github.com/greenelab/czi-rfa/issues/11] that numerous manuscripts demonstrated the advantages of zero-inflated loss in the setting of single cell sequencing data [@doi:10.1186/s13059-015-0805-z @arxiv:1610.05857 @doi:10.1186/s13059-017-1188-0].
These efforts will allow us to identify a subset of parameter combinations that are worth exploring with new approaches designed specifically for this type of genomic data.

We will also develop a data augmentation approach for training deep models on single cell RNA-seq data, as no data augmentation approaches for transcriptomic data have been published to date.
These approaches are widely used in image analysis, where the goal is to distinguish relevant from irrelevant features in an image corpus.
To understand data augmentation, imagine pathology slides that have been scanned.
Each slide may be prepared and scanned in a subtly different orientation or at somewhat different magnifications depending on who generated the data.
A deep learning method may learn to identify these subtle differences, which is undesirable, or there may simply be too few slides to train a deep learning algorithm.
To address these challenges, deep learning practitioners can apply arbitrary rotations, zooms, and other irrelevant transformations to image data during training.
This process is called data augmentation.

We expect that very fast abundance estimates from RNA-seq data [@doi:10.1038/nmeth.4197 @10.1038/nbt.3519] will make a data augmentation approach feasible in this domain as well.
Resampling reads to generate abundance estimates can help to capture the uncertainty in the data, akin to arbitrary rotations.
Subsampling reads to generate abundance estimates can help to capture changes related to sequencing depth but unrelated to the biology, akin to arbitrary zooms.
Therefore, we plan to collaborate with Rob Patro's laboratory via our collaborative network to implement approaches including but not limited to rapid subsampling and bootstrapping to generate augmented training data.
We posit that genomic data augmentation will improve latent feature generalization by separating biological from technical features and increasing the effective sample size during training.

We will select high-quality models by choosing those that minimize both reconstruction loss and a KL divergence term which constrains the features to a Gaussian distribution [@arxiv:1312.6114].
We will provide a source code repository with this implementation, as well as a summary of our findings, and we may produce a manuscript on the topic.
Resulting models that exhibit desirable properties will be evaluated in the context of the rheumatic disease compendium (Aim 2) as well as their suitability for latent space arithmetic (see: Evaluation).

#### Aim 2: Generate a benchmark dataset of harmonized public data to evaluate the extent to which HCA cell types capture disease states for multiple rheumatic diseases.

The HCA's partnership with the Immunological Genome Project (immgenH) will provide single-cell gene expression-based immunocyte phenotyping at an unprecedented resolution.
Ultimately, we expect that HCA-identified cell populations will translate to an improved understanding of health and disease.
A compendium comprised of bulk gene expression data from autoimmune/rheumatic diseases is exceptionally well-suited for the biologically-motivated evaluation (e.g., benchmarking) of these immunocyte data.
The _objective of this aim_ is to build such a compendium and to evaluate methods, including VAEs and other methods developed through this initiative, for defining cell-type-specific expression signatures from the HCA's single-cell datasets by measuring their performance in the context of decomposing bulk, whole-tissue autoimmune/rheumatic disease data.
We will develop both real and simulated benchmark datasets to evaluate the quality of the cell-type-specific signatures (and by extension the methods) for deconvolution and the identification of correct cell-type proportions and phenotypes.

We will generate simulated bulk datasets drawn from purified cell lineages by combining the purified cell expression data at different proportions.
We will also build a multi-tissue autoimmune/rheumatic disease compendium from existing public data of bulk tissues.
We have curated a relevant set of datasets to date of more than 12,000 samples.
This compendium includes samples from patients with systemic lupus erythematosus (SLE), sarcoidosis, and inflammatory bowel disorders among many other diseases.
This compendium will have desirable properties for evaluating single-cell data-derived signatures, as it will allow us to evaluate cell type proportions and phenotypes based on a body of previous literature.
For instance, we expect to detect higher proportions of activated macrophages in lupus nephritis samples as compared to controls [@doi:10.4049/jimmunol.1103031].

Importantly, a bulk compendium comprised of public data from a disease context will enable HCA participants (methods developers, RNA-seq based assay developers, and others) to directly evaluate approaches in the context where we expect their most immediate impact: application to existing datasets to explain disease-relevant phenomena via a single-cell perspective.

### Proposal for evaluation and dissemination.

We will apply methods that produce low-dimensional representations including VAEs (Aim 1) and other methods to HCA-produced single cell transcriptomes.
These low-dimensional representations will be used in evaluations to test latent space arithmetic and relevance to disease.
Source code that reproducibly generates these models will be released via GitHub as it is developed.
Models will be disseminated via periodic release on Zenodo or a similar platform.

#### Evaluate the extent to which low-dimensional representations enable latent space arithmetic in the HCA.

Certain classes of generative deep neural network models, including VAEs and GANs, add a distribution matching constraint that has been shown to imbue intuitive mathematical features to the learned latent features [@arxiv:1312.6114 @arxiv:1406.2661].
For instance, a GAN learned latent features that could be manipulated with arithmetic: Subtracting out the essence of a smile from a woman and adding it to a neutral man resulted in an image vector of a smiling man [@arxiv:1511.06434].
We will evaluate the extent to which this is true for low-dimensional representations of the HCA's single-cell transcriptomes using benchmark data that contains parallel perturbations, which allows us to hold out one or more perturbations for one or more cell types to evaluate each method.
We describe two experiments using data proposed by Arjun Raj's group, but any HCA benchmark datasets with similar properties will be suitable.
We are happy to expand our evaluation to additional suitable datasets, and we are happy to evaluate low-dimensional models produced by other groups.

The Raj lab proposes to generate data related to cardiomyocte differentiation from fibroblasts.
The driving transcription factors for this process have been identified [@doi:10.1016/j.cell.2010.07.002].
The latent space vector between these two cell types should capture the key transcription factor networks (Gata4, Mef2c, and Tbx5).
To calculate this vector we will subtract the latent space projections of human fibroblasts from cardiomyocytes.
We will compare the gene composition of this differentiation vector to target calls from cistrome [@doi:10.1093/nar/gkw983], which are available for each of these TFs.

Latent space arithmetic can also be used to generate new hypothetical data, and we will test the extent to which these models can predict the results of perturbations.
Arjun Raj's lab is generating data with homogenized cell-types under various perturbations.
For each perturbation, we will hold out one or more cell types and map the rest into the latent space.
Subtracting the latent space vector of included cell types from the same cells after perturbation will produce a perturbation vector.
We will add the perturbation vector to one of the cell types for which the perturbation was withheld to generate synthetic data, and we will compare the results to the actual held out data to determine the error of the prediction.
Comparing all results for low-dimensional methods to the baseline method of performing the same subtraction operations in raw gene expression space will allow us to determine the extent to which these models capture perturbation information.

#### Rheumatic Disease Evaluation

We will input signatures from these methods into existing techniques [@doi:10.1186/s13059-016-1070-5 @doi:10.1093/bioinformatics/btv015] that decompose bulk data with cell type signatures and evaluate performance (e.g., prediction of proportions) on the simulated data from Aim 2.
By using comparing performance across multiple decomposition techniques we can measure each method's ability to define bulk-relevant signatures from single cell data in datasets with an available ground truth.
We will also use signatures to decompose the rheumatic disease compendium (Aim 2).
This compendium includes drug studies of highly-targeted therapeutics (e.g., a monoclonal antibody to IFN-gamma [anti-IFNg] in the context of systemic lupus erythematosus [@doi:10.1002/art.39248]).
This allows us to determine what cell-types change in proportion during the reduction of this cytokine or if the IFN-inducible expression is preferentially altered in one cell type.
At this stage, we will use this type of analysis to identify cases where methods disagree.
Disagreements that we identify can be directly probed in future experiments (beyond the one-year timeline) to produce periodic ground-truth benchmarks.
Real and simulated datasets will be disseminated as they are developed via periodic release on Zenodo or a similar platform.

### Statement of commitment to share

We look forward to sharing our proposals, methods, data, and code with other researchers funded by this RFA, with CZI, and with the broader research community.

* We have shared our proposal under a CC-BY license as it was being written.
* We will share our description of methods, written reports, figures, and other such outputs under a CC-BY license.
* We will share source code for our software under the BSD 3-clause license.

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
    * cgreene on github.
    * Dr. Greene has primary responsibility for managing and carrying out the proposed research, interpreting findings, disseminating findings, and preparing and submitting progress reports.

## Description of other key personnel

**Gregory Way, MS (GCB Graduate Student).**
_4.2 Cal Months_.
Gregory Way (gwaygenomics on GitHub) is a graduate student in the Genomics and Computational Biology PhD program at the University of Pennsylvania performing thesis research in Dr. Greene's laboratory.
Mr. Way has developed supervised and unsupervised methods to extract knowledge from bulk cancer transcriptomes, including the training of a VAE.
He will, with Dr. Greene, extend this VAE to learn single cell features.
He is the first author on the bulk VAE manuscript that serves as preliminary data for Aim 1.
Together with Dr. Hu, he will be responsible for developing Aim 1 – namely, training an optimal single cell VAE with an augmentation approach.
He will also lead the latent space arithmetic approach proposed in Aim 3.
Mr. Way will participate in experimental design, analysis, and manuscript preparation.

**Qiwen Hu, PhD (Postdoctoral Researcher).**
_6.0 Cal Months_.
Qiwen Hu (huqiwen0313 on GitHub) is a postdoc at the University of Pennsylvania in Dr. Greene's laboratory.
Dr.Hu has developed supervised and unsupervised methods to learn genomic translational regulatory features.
She will, with Dr. Greene, extend the methods to learn single cell features.
Together with Mr. Way, she will be responsible for developing Aim 1 – namely, training an optimal single cell VAE with an augmentation approach.
She will also work on the latent space arithmetic approach proposed in Aim 3.
Dr. Hu will participate in experimental design, analysis, and manuscript preparation.

**Jaclyn Taroni, PhD (Postdoctoral Researcher).**
_2.4 Cal Months_.
Jaclyn Taroni (jaclyn-taroni on GitHub) is a postdoctoral researcher in the Dr. Greene’s laboratory.
Dr. Taroni has expertise in rheumatic diseases and in developing computational approaches to study transcriptomic data in these disorders.
She has developed frameworks for integrating large compendia of bulk transcriptomic data and will be responsible for building the compendium of bulk rheumatic disease data in Aim 2.
Dr. Taroni will evaluate the biological relevance of the low-dimensional representations of single-cell data learned by various methods in the context of rheumatic diseases.
She will participate in experimental design, analysis, and manuscript preparation.

## Brief preliminary budget (Attached)
