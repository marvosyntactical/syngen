# Talk 10.09.25

* [NOTE] use ncbi for genomic side
* [NOTE] several tools for ncbi, also new ones
* [NOTE] start with bacdive side of things
* [TODO] construct db scheme 
* [NOTE]
* [NOTE]
* [NOTE]
* [NOTE]

# TODOs

* [TODO] look at Oscar's scheduler - interacts with momentum based optimisers!
* [TODO] To pick up the key on Monday, September 8th, between 9.00 and 11.00 am, please go to the reception desk in the main DKFZ building at INF 280
* [TODO] go to B330, mkon B, and find my office
* [TODO] go through onboarding documents
* [TODO] wait for VPN access
* [TODO] once I have access (on premise/VPN), go through FIRST STEPS BEGINNER
* [TODO] whats a casino card???
* [NOTE] module & LSF
* [TODO] ocd jupyterhub?
* [NOTE] use uv? venv? conda? probably not venv; probably conda
* [NOTE] When I leave: BRING CAKE !!! THE CAKE IS NOT A LIE!
* [NOTE] Oscar ist ab Mittwoch da
* [TODO] welcher projektordner ist relevant?

## Thoughts on Pipeline

* matching - MED/Xent? FoM? 
* how to go from LLM to conditional generation?



# Thoughts before interview

* ncbi refseq? https://www.ncbi.nlm.nih.gov/refseq/
* ENA https://www.ebi.ac.uk/ena/browser/home
* compare to:
* bacdive https://bacdive.dsmz.de/

* dim reduction techniques first data science analysis lens

* matching system - something like minimum edit distance or crossentropy based? figure of merit?

* quality control - unit tests (sequence integrity etc, sanity checks)

* standardized annotation format for pheno trait labeling - do we have labels? in which database? BacDive?

* tokenisation? 1 nucleotide, 3? i think sometimes single nucleotides code for shit

* untranslated regions orchestrate gene expression, opening helix/chromatine; coding vs noncoding region?!
* single or multi nucleotide tokenisation? evo=single, aminoacids encoded by 3 nucleotides, BPE?
* where sequence starts depends on many factors, other genes etc etc etc
* dna language model are prototype
* new step: condition on output phenotype
* bacteria: use rna, single strand dna, junk rna?, bact seqs are simple, regulatory mechanism also simpler, length of seqs also much shorter, databases exist
* protein language models most advanced

# Interpretability 

* Steering!?


# Jobbeschreibung

HiWi Position – Generative AI in Genomics (DKFZ, Applied Bioinformatics Group)

The Applied Bioinformatics Group (Dr. Óscar González-Velasco) at the German Cancer Research Center (DKFZ) is offering a 6-month HiWi position for highly motivated Master’s students. The project focuses on applying generative AI methods to bacterial genomics and is available to start as soon as possible.

Preferred Background:
Applicants should have a background in computer science, mathematics, physics, or molecular biotechnology, with strong computational skills.

Project Overview: Generative AI in Genomics

1. Database framework development (50% time allocation)

* Creating an automated pipeline to extract, transform, and load bacterial genome sequences from repositories like NCBI and ENA
* Developing a matching system to connect genome sequences with phenotypic data from BacDive and other database sources
* Implementing quality control measures to ensure data integrity and harmonization
* Creating a standardized annotation format for consistent labeling of phenotypic traits

2. Baseline AI model construction (50% time allocation)

* Implementing data preprocessing workflows for genome sequence tokenization
* Developing and testing multiple attention-based architectures with conditional generation capabilities (e.g. EVO2 architecture, other LLM transformer-based architectures)
* Performing hyperparameter optimization for small-scale proof-of-concept models
* Establishing evaluation metrics to measure model performance on synthetic sequence generation tailored to phenotypic traits
