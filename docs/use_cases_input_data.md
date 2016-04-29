We are now entering into real analyses using Metavisitor !
In this section, we are going to create step by step Galaxy histories that contains the input data required to execute the analyses as described in the [metavisitor article](http://dx.doi.org/10.1101/048983). These analyses are presented as use cases. Thus let's start with the input data to run the Use Cases 1-1, 1-2, 1-3 and 1-4

# History with input data for Use Cases 1-1, 1-2, 1-3 and 1-4

1. Create a new history and rename it Input data for Use Cases 1-1, 1-2, 1-3 and 1-4
- import SRP013822 datasets
    - Use the tool `Extract reads in FASTQ/A format from NCBI SRA` and fill the SRR accession field with the first EBI SRA identifier **SRR515090**
    - repeat the exact same operation with the tool `Extract reads in FASTQ/A format from NCBI SRA` and the identifiers **SRR513993, SRR513992, SRR513990, SRR513989, SRR513981, SRR513901**
    - for the 7 datasets retrieved from EBI SRA, change the datatype **fastq** to **fastqsanger**:
    click on the pencil icon of the dataset, click the tab Datatype, and select fastqsanger in the New Type menu
- Create a dataset collection **SRP013822**
    - Click on the checked box icon at the top of history bar as indicated below
    
    ![checkbow](images/check_datasets.png)
    
    - Select the 7 datasets
    - Select **Build Dataset List** in the menu **For all selected...**
    
    ![datasetlist](images/dataset_list.png)
    
    - and type **SRP013822** in the Name field and click `Create list`
    - you can leave the checked datasets view by clicking again the check box in the history top menu
    
