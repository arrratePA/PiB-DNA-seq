# PiB-DNA-seq
DNA-seq pipeline
Background: Inactivating PTH/PTHrP signaling disorders are heterogeneous group of rare diseases in which several genes and genetic mechanisms 
(mosaicism, copy number variations or structural variants) are involved. In this context, the iPPSD-seq Next Generation DNA sequencing panel was 
designed for the study of this of patients with suspicion of this type of disorders. For the analysis of the panel sequencing data, a GATK Best Practices workflow 
based pipeline is required, capable to identified different types of variants including mosaics.
Findings: In the EpiG lab, the analysis of the iPPSD-seq panel sequencing data has been classically performed with a commercial pipeline that do not follow 
the GATK Best Practice guidelines. Through this study, the iPPSD-seq panel workflow has been implemented, and PiB-DNA-seq pipeline has been developed according 
to GATK. This workflow has been validated with data from positive and negative controls samples of the lab historical database.
Conclusion: The PiB-DNA-seq pipeline not only allowed to identify some limitations of the classic workflow but also to improve the quality metrics of 
the sequencing data. In addition, it demonstrated its utility to identified even mosaic variants.
