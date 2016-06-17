# Metavisitor

[Metavisitor](http://dx.doi.org/10.1101/048983) is a user-friendly and adaptable software to provide biologists, clinical researchers and possibly diagnostic clinicians with the ability to robustly detect and reconstruct viral genomes from complex deep sequence datasets. A set of modular bioinformatic tools and workflows was implemented as the Metavisitor package in the Galaxy framework. Using the graphical Galaxy workflow editor, users with minimal computational skills can use existing Metavisitor workflows or adapt them to suit specific needs by adding or modifying analysis modules.

## Availability of Metavisitor tools and workflows

Metavisitor has been developed at the [ARTbio platform](http://artbio.fr). Its tools and workflows are primarily available in [GitHub] (https://github.com/ARTbio/tools-artbio).

#### Metavisitor tools developed by ARTbio in the [ARTbio GitHub](https://github.com/ARTbio/tools-artbio)

- [`yac_clipper`](https://github.com/ARTbio/tools-artbio/tree/master/tools/yac_clipper)
- [`concatenate_multiple_datasets`](https://github.com/ARTbio/tools-artbio/tree/master/tools/concatenate_multiple_datasets)
- [`msp_sr_bowtie`](https://github.com/ARTbio/tools-artbio/tree/master/tools/msp_sr_bowtie)
- [`msp_fasta_tabular_converter`](https://github.com/ARTbio/tools-artbio/tree/master/tools/msp_fasta_tabular_converter)
- [`fetch_fasta_from_ncbi`](https://github.com/ARTbio/tools-artbio/tree/master/tools/fetch_fasta_from_ncbi)
- [`msp_cap3`](https://github.com/ARTbio/tools-artbio/tree/master/tools/msp_cap3)
- [`cherry_pick_fasta`](https://github.com/ARTbio/tools-artbio/tree/master/tools/cherry_pick_fasta)
- [`msp_blastparser_and_hits`](https://github.com/ARTbio/tools-artbio/tree/master/tools/msp_blastparser_and_hits)
- [`msp_oases`](https://github.com/ARTbio/tools-artbio/tree/master/tools/msp_oases)
- [`blast_to_scaffold`](https://github.com/ARTbio/tools-artbio/tree/master/tools/blast_to_scaffold)
- [`blastx_to_scaffold`](https://github.com/ARTbio/tools-artbio/tree/master/tools/blastx_to_scaffold)
- [`msp_sr_readmap_and_size_histograms`](https://github.com/ARTbio/tools-artbio/tree/master/tools/msp_sr_readmap_and_size_histograms).
    
Other tools from other developers are included the suite_metavisitor_1_2. Of course, these tools are not available on our GitHub repository, but are available from the [main Galaxy toolshed](https://toolshed.g2.bx.psu.edu/):
    
     name="sra_tools" owner="iuc"
     name="get_orfs_or_cdss" owner="peterjc" 
     name="trinityrnaseq" owner="anmoljh" 
     name="bowtie2" owner="devteam" 
     name="fastx_trimmer" owner="devteam" 
     name="fastq_to_fasta" owner="devteam" 
     name="fasta_filter_by_length" owner="devteam" 
     name="regex_find_replace" owner="jjohnson" 
     name="khmer_normalize_by_median" owner="iuc" 


#### Availability of Metavisitor tools and workflows for **Galaxy instance administrators**

All metavisitor tools are available in the [suite_metavisitor_1_2](https://toolshed.g2.bx.psu.edu/repository/browse_repositories?sort=name&operation=view_or_manage_repository&f-free-text-search=metavisitor&id=ca18473f5a7e691a)
Galaxy Admin can just install this suite of tools by using the `Search Tool Shed` menu in their Admin panel, searching for "metavisitor", and installing the `suite_metavisitor_1_2` tool suite.

Admins can also install the tools from the `metavisitor_workflows` repository, which will provide in addition the metavisitors workflows.

#### Availability of Metavisitors workflows for any Galaxy instance user.
We have deposited the Metavisitors workflows in the [myexperiment server](http://www.myexperiment.org/workflows), where they are searchable with "metavisitor" keyword and can be downloaded and reuploaded to the Galaxy instance.

#### Users who have already the Metavisitor suite of tools installed in their Galaxy instance, or who just want to use it on the Galaxy Mississippi Server
can skip these chapters and go directly to the chapter [Prepare input data histories](use_cases_input_data) and following.


#### In the next three chapters

We provide documentation on two methods to set up a Galaxy server instances *from scratch* with *pre-installed* Metavisitor tools and workflows.

- Based on GalaxyKickstarter: see [Metavisitor with GalaxyKickstarter (Ansible)](metavisitor_ansible.md)
- Based on Docker: see [Metavitor with Docker](metavisitor_docker.md)


