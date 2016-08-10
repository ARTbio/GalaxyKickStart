## Input data for Use Case 3-3

As for the previous Use Cases, the first step is to collect all input data in an history that we will name `Input data for Use Case 3-3`. 

- Create a new history
- Rename this history `Input data for Use Case 3-1`
- Using the tool `Extract reads in FASTQ/A format from NCBI SRA`, we are going to upload 63 paired end datasets. As explained in the previous chapter, the `Extract reads in FASTQ/A format from NCBI SRA` directly merges paired-end sequence data in a single file. As for the previous `Use Case 3-2`, collecting these 63 datasets is tedious. But the good news is that we are going to use repeatedly  the same `Extract reads in FASTQ/A format from NCBI SRA` tool and we only need in addition to change 63 times the datatype of these datasets from `fastq` to `fastqsanger`.

So, here are the SSR ids to enter in the `Extract reads in FASTQ/A format from NCBI SRA` tool form (leaving `fastq` as an output format):

##### For Ebola virus samples:

`SRR1613381`
`SRR1613377`
`SRR1613382`
`SRR1613378`
`SRR1613383`
`SRR1613379`
`SRR1613384`
`SRR1613380`

When you have finishid with the 8 datasets, make sure to change their datatype from `fastq`to `fastqsanger`, and create a dataset collection (as explained in the previous chapter) of these 8 datasets that you will name `Ebola virus`.

##### For Lassa virus samples:

`SRR1595772`
`SRR1595696`
`SRR1595665`
`SRR1595500`
`SRR1594619`
`SRR1595943`
`SRR1595673`
`SRR1595797`
`SRR1595763`
`SRR1595558`
`SRR1594664`
`SRR1595909`
`SRR1594651`
`SRR1595835`
`SRR1594698`
`SRR1613388`
`SRR1613389`
`SRR1613390`
`SRR1613391`
`SRR1613392`
`SRR1613393`
`SRR1613394`
`SRR1613395`
`SRR1613396`
`SRR1613397`
`SRR1613398`
`SRR1613399`
`SRR1595853`
`SRR1606288`
`SRR1613412`
`SRR1613403`
`SRR1606277`
`SRR1613386`
`SRR1613387`
`SRR1606267`
`SRR1614275`
`SRR1610580`
`SRR1595846`
`SRR1594606`
`SRR1606236`
`SRR1594723`
`SRR1594671`
`SRR1613414`
`SRR1613400`
`SRR1613401`
`SRR1613404`
`SRR1613402`
`SRR1613405`
`SRR1613407`
`SRR1613408`
`SRR1613409`
`SRR1613410`
`SRR1613406`
`SRR1613411`
`SRR1613413`

When you have finishid with the 55 datasets, make sure to change their datatype from `fastq`to `fastqsanger`, and create a dataset collection (as explained in the previous chapter) of these 8 datasets that you will name `Lassa virus`.

- Copy the `vir1 nucleotide BLAST database` from the `References` history to the current history `Input data for Use Case 3-3`.

## History for Use Case 3-3 / Ebola virus
1. Stay in the history `Input data for Use Case 3-3`
- pick the workflow `Metavisitor: Workflow for Use Case 3-3` in the workflows menu. It is possible that the workflow manager complains about settings for the Trinity tool that is used in this workflow. This is a minor issue if happens: just edit the workflow, click on the tool `Trinity` and specify the number of processors accordingly to your computing infrastructure. Save the workflow and select the `run` option.
- Before Step 1, you have to specify some parameters at run time. For Ebola virus, the field `target_virus` has to be filled with `Ebola` and the field `reference_virus` has to be filled with `NC_002549.1` (as a guide for reconstruction of the Ebola virus genome).
- For Step 1, select `Ebola virus`.
- For Step 2, select the `nucleotide vir1 blast database` (this should also be already selected)
- As usual, check the box `Send results to a new history`, edit the name of the new history to `Use Case 3-3 Ebola virus`, and `Execute` the workflow ! Note, that for complex workflows with dataset collections in input, the actual warning that the workflow is started make take time to show up; you can even have a "504 Gateway Time-out" warning. This is not a serious issue: just go in your `User` -> `Saved history` menu, you will see your `Use Case 3-3 Ebola virus` history running and you will be able to access it.

The workflow for Use Case 3-3 may take a long time. Be patient.

## History for Use Case 3-3 / Lassa virus, segment L
1. Stay in the history `Input data for Use Case 3-3`
- pick the workflow `Metavisitor: Workflow for Use Case 3-3` in the workflows menu. It is possible that the workflow manager complains about settings for the Trinity tool that is used in this workflow. This is a minor issue if happens: just edit the workflow, click on the tool `Trinity` and specify the number of processors accordingly to your computing infrastructure. Save the workflow and select the `run` option.
- Before Step 1, you have to specify some parameters at run time. For Lassa virus, the field `target_virus` has to be filled with `Lassa` and the field `reference_virus` has to be filled with `NC_004297.1` (as a guide for reconstruction of the segment L of the Lassa virus genome).
- For Step 1, select `Lassa virus`.
- For Step 2, select the `nucleotide vir1 blast database` (this should also be already selected)
- As usual, check the box `Send results to a new history`, edit the name of the new history to `Use Case 3-3 Lassa virus segment L`, and `Execute` the workflow ! Note, that for complex workflows with dataset collections in input, the actual warning that the workflow is started make take time to show up; you can even have a "504 Gateway Time-out" warning. This is not a serious issue: just go in your `User` -> `Saved history` menu, you will see your `Use Case 3-3 Lassa virus segment L` history running and you will be able to access it.

The workflow for Use Case 3-3 may take a long time. Be patient.

## History for Use Case 3-3 / Lassa virus, segment S
1. Stay in the history `Input data for Use Case 3-3`
- pick the workflow `Metavisitor: Workflow for Use Case 3-3` in the workflows menu. It is possible that the workflow manager complains about settings for the Trinity tool that is used in this workflow. This is a minor issue if happens: just edit the workflow, click on the tool `Trinity` and specify the number of processors accordingly to your computing infrastructure. Save the workflow and select the `run` option.
- Before Step 1, you have to specify some parameters at run time. For Lassa virus, the field `target_virus` has to be filled with `Lassa` and the field `reference_virus` has to be filled with `NC_004296.1` (as a guide for reconstruction of the segment S of the Lassa virus genome).
- For Step 1, select `Lassa virus` (this should be already selected).
- For Step 2, select the `nucleotide vir1 blast database` (this should also be already selected)
- As usual, check the box `Send results to a new history`, edit the name of the new history to `Use Case 3-3 Lassa virus segment S`, and `Execute` the workflow ! Note, that for complex workflows with dataset collections in input, the actual warning that the workflow is started make take time to show up; you can even have a "504 Gateway Time-out" warning. This is not a serious issue: just go in your `User` -> `Saved history` menu, you will see your `Use Case 3-3 Lassa virus segment S` history running and you will be able to access it.

The workflow for Use Case 3-3 may take a long time. Be patient.
