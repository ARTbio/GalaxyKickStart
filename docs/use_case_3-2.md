## Input data for Use Case 3-2

As for the previous Use Cases 1, 2 and 3-1, the first step is to collect all input data in an history that we will name `Input data for Use Case 3-2`. 

- Create a new history
- Rename this history `Input data for Use Case 3-1`
- Using the tool `Extract reads in FASTQ/A format from NCBI SRA`, we are going to upload 42 paired end datasets. Indeed, these 42 datasets correspond to 84 fastq paired-ended sequence files. However, the `Extract reads in FASTQ/A format from NCBI SRA` directly merges two paired-end fastq datasets in a single file. In addition, some datasets derive from the same patient; in those cases we will merge those datasets using the tool `Concatenate multiple datasets tail-to-head` and delete and purge the original datasets. In all cases, we will rename the dataset with the patient id as indicated bellow, and change the datatype from fastq to fastqsanger.
Here is a table that recapitulates the actions to perform, line by line.
```
#ENA-RUN            action		                                        post-action
SRR453487           |rename	"patient 566"	                            |change datatype to "fastqsanger"
SRR453437	        |rename	"patient 438"	                            |change datatype to "fastqsanger"
SRR453443,SRR453458	|concatenate, rename merge dataset	"patient 401"	|change the merged dataset to datatype to "fastqsanger", and delete and purge the original SRR datasets
SRR453430	        |rename	"patient 382"	                            |change datatype to "fastqsanger"
SRR453491	        |rename	"patient 377"	                            |change datatype to "fastqsanger"
SRR453499	        |rename	"patient 375"	                            |change datatype to "fastqsanger"
SRR453484	        |rename	"patient 350"	                            |change datatype to "fastqsanger"
SRR453464	        |rename	"patient 349"	                            |change datatype to "fastqsanger"
SRR453506	        |rename	"patient 345"	                            |change datatype to "fastqsanger"
SRR453417	        |rename	"patient 344"	                            |change datatype to "fastqsanger"
SRR453490	        |rename	"patient 335"	                            |change datatype to "fastqsanger"
SRR453478	        |rename	"patient 331"	                            |change datatype to "fastqsanger"
SRR453465,SRR453480	|concatenate, rename merge dataset	"patient 330"	|change the merged dataset to datatype to "fastqsanger", and delete and purge the original SRR datasets
SRR453489,SRR453505	|concatenate, rename merge dataset	"patient 329"	|change the merged dataset to datatype to "fastqsanger", and delete and purge the original SRR datasets
SRR453498	        |rename	"patient 322"	                            |change datatype to "fastqsanger"
SRR453446	        |rename	"patient 321"	                            |change datatype to "fastqsanger"
SRR453427,SRR453440	|concatenate, rename merge dataset	"patient 315"	|change the merged dataset to datatype to "fastqsanger", and delete and purge the original SRR datasets
SRR453438	        |rename	"patient 282"	                            |change datatype to "fastqsanger"
SRR453450	        |rename	"patient 275"	                            |change datatype to "fastqsanger"
SRR453460	        |rename	"patient 274"	                            |change datatype to "fastqsanger"
SRR453485	        |rename	"patient 270"	                            |change datatype to "fastqsanger"
SRR453448	        |rename	"patient 266"	                            |change datatype to "fastqsanger"
SRR453424,SRR453457	|concatenate, rename merge dataset	"patient 263"	|change the merged dataset to datatype to "fastqsanger", and delete and purge the original SRR datasets
SRR453510	        |rename	"patient 193"	                            |change datatype to "fastqsanger"
SRR453456	        |rename	"patient 187"	                            |change datatype to "fastqsanger"
SRR453425,SRR453469	|concatenate, rename merge dataset	"patient 186"	|change the merged dataset to datatype to "fastqsanger", and delete and purge the original SRR datasets
SRR453481	        |rename	"patient 183"	                            |change datatype to "fastqsanger"
SRR453531	        |rename	"patient 180"	                            |change datatype to "fastqsanger"
SRR453474	        |rename	"patient 179"	                            |change datatype to "fastqsanger"
SRR453509	        |rename	"patient 171"	                            |change datatype to "fastqsanger"
SRR453451	        |rename	"patient 168"	                            |change datatype to "fastqsanger"
SRR453495,SRR453504	|concatenate, rename merge dataset	"patient 161"	|change the merged dataset to datatype to "fastqsanger", and delete and purge the original SRR datasets
SRR453500	        |rename	"patient 159"	                            |change datatype to "fastqsanger"
SRR453493	        |rename	"patient 156"	                            |change datatype to "fastqsanger"
SRR453444	        |rename	"patient 131"	                            |change datatype to "fastqsanger"
SRR453426	        |rename	"patient 78	                                |change datatype to "fastqsanger"
```
- Create a dataset collection of patient datasets: Click on the checked box icon in the history top menu, check the "patient... " datasets we have just generated (36 datasets) (you can use the `Select all` button), and `For all selected`, `Build a dataset list` that you name "Tractable Patient Datasets".
- Copy the `vir1 nucleotide BLAST database` from the `References` history to the current history `Input data for Use Case 3-2`.

## History for Use Case 3-2
1. Stay in the history `Input data for Use Case 3-2`
- pick the workflow `Metavisitor: Workflow for Use Case 3-2` in the workflows menu, and select the `run` option.
- For Step 1 (Fever Patient Sequences collection), select `Tractable Patient Datasets` (this should be already selected).
- For Step 2, select the `nucleotide vir1 blast database` (this should also be already selected)
- As usual, check the box `Send results to a new history`, edit the name of the new history to `History for Use Case 3-2`, and `Execute` the workflow ! Note, that for complex workflows with dataset collections in input, the actual warning that the workflow is started make take time to show up; you can even have a "504 Gateway Time-out" warning. This is not a serious issue: just go in your `User` -> `Saved history` menu, you will see you `History for Use Case 3-2` running and you will be able to access it.

As a last note, the workflow for Use Case 3-2 may take a long time. Be patient.
















