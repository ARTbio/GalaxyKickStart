# Histories for Use Cases 2-1, 2-2

Now that you get more familiar with manipulations in Galaxy with the Use Cases 1-1 to 1-4 described in details in the previous chapters, we will describe the other Use Case analyses more concisely. If you experience lack of skills in basic Galaxy operations (tool usage, copy of datasets, etc), do not hesitate to go back and examine the previous chapters step by step.


## Input data for Use Cases 2-1 and 2-2

As for the previous Use Case 1, the first step is to collect all input data in an history that we will name `Input data for Use Cases 2-1 and 2-2`

1. Create a new history
2. Rename this history `Input data for Use Cases 2-1 and 2-2`
3. For the small RNA sequence datasets (ERP012577) in this study, we are going to use another tool to upload to the Galaxy Metavisitor server: the `EBI SRA ENA SRA`tool which in the "Get data" section of the left tool bar.
    - click on this tool and enter ERP012577 in the search field that shows up in the European Nucleotide Archive web page, and search. Click on the `ERP012577` link. In the column "Submitted files (galaxy)" of the table, click on the first "fastq file 1". This action should send you back to your Galaxy page automatically and you see the fastq dataset loading (yellow dataset in the history bar).
    - repeat the exact same operation, for the three other "fastq file 1".
    - at final you should upload four fastq datasets corresponding to the sequencing runs "post_infected_rep1.fastq", "post_infected_rep2.fastq", "post_non-infected_rep1.fastq" and "post_non-infected_rep2.fastq"
    - Once the 4 uploads are _completed_ (may takes minutes, depending on your network speed connection), click on the pencil icon of the 4 datasets, click on the `datatype` tab and get it to `fastqsanger`.
4. Create a dataset collection as [previously explained](use_cases_input_data/#history-with-input-data-for-use-cases-1-1-1-2-1-3-and-1-4) and name it `Small RNA reads ERP012577`
5. For the RNA sequence datasets (ERS977505) that will be used in Use Case 2-2, use again the `EBI SRA ENA SRA`tool which in the "Get data" section of the left tool bar.
    - click on this tool and enter ERS977505 in the search field that shows up in the European Nucleotide Archive web page, and search. Click on the `ERS977505` link (Sample 1 result found). In the column "Submitted files (galaxy)" of the table, click on the first "fastq file 1". This action should send you back to your Galaxy page automatically and you see the fastq dataset loading (yellow dataset in the history bar).
    - repeat the exact same operation for the other "fastq file 1" and the two other "fastq file 2"
    - at final you should upload four additional fastq datasets corresponding to the sequencing runs "IP-isoT-1_AGTCAA_L001_R_1.fastq", "IP-isoT-1_AGTCAA_L001_R_2.fastq", "IP-isoT-2_ATGTCA_L002_R_1.fastq" and "IP-isoT-2_ATGTCA_L002_R_2.fastq"
6. Create a dataset collection as explained in the previous chapter and name it `long read RNAseq datasets`
7. Using the "Upload file tool" as explained [before](metavisitor_configure_references/#b-upload-nucleotide-vir1-fasta-file), upload the Plasmodium berghei genome by pasting this URL in the `Paste/Fetch Data` tab of the tools:
```
ftp://ftp.ensemblgenomes.org/pub/release-28/protists/fasta/plasmodium_berghei/dna/Plasmodium_berghei.May_2010.28.dna_sm.genome.fa.gz
```
8. Use the `Retrieve FASTA from NCBI`, paste `phix174[title]` in the "Query to NCBI in entrez format" field and select `nucleotide` for the NCBI database. This will upload 174 fasta sequences from phix174.
9. Use the wheel icon at the top of the history bar to copy `nucleotide vir1 blast database` and `protein vir1 blast database` **from** the history `References` **to** the current history `Input data for Use Cases 2-1 and 2-2`. If you don't remember well how to copy datasets between histories, you may read again the explanation [here](use_cases_input_data/#history-with-input-data-for-use-cases-1-1-1-2-1-3-and-1-4) (step 4.)

**_Your are now ready for generating Uses Cases 2-1 and 2-2_**
    
## History for Use Case 2-1

1. Stay in the current history `Input data for Use Cases 2-1 and 2-2` !
2. In the `Workflow` menu, select the workflow `Metavisitor: Workflow for Use Case 2-1` and directly select `Run` (you may also look at the workflow using the `edit` option)
3. Be careful at selecting `Small RNA reads ERP012577` for the step 1 (Input Dataset Collection)
4. For the step 2, the option `protein vir1 blast database` is forced, because the workflow is expecting of protein blast database for this step and only one dataset with this datatype is available in the history
5. **Be careful** at selecting
```
ftp://ftp.ensemblgenomes.org/pub/release-28/protists/fasta/plasmodium_berghei/dna/Plasmodium_berghei.May_2010.28.dna_sm.genome.fa.gz
```
for **step 10** (sRbowtie)
6. Be careful at selecting
```
Retrieve FASTA from NCBI (Nucleotide) with queryString 'phix174[title]'
```
for step 11 (sRbowtie).
7. Click the `Send results to a new history` checkbox and rename the history to "History for Use Case 2-1".
8. Run Workflow !

You may follow the link to the new history when the workflow is started.

## History for Use Case 2-2

1. If you are not already in, go back to the history `Input data for Use Cases 2-1 and 2-2`
2. In the `Workflow` menu, select the workflow `Metavisitor: Workflow for Use Case 2-2` and directly select `Run` (you may also look at the workflow using the `edit` option)
3. Be careful at selecting `long read RNAseq datasets` for the step 1 (Input Dataset Collection)
4. For the step 2, the option `protein vir1 blast database` is forced, because the workflow is expecting of protein blast database for this step and only one dataset with this datatype is available in the history
5. Click the `Send results to a new history` checkbox and rename the history to "History for Use Case 2-1".
6. Run Workflow !


 















