Now that you get more familiar with manipulations in Galaxy with the Use Cases 1-1 to 1-4 described in details in the previous chapters, we will describe the other Use Case analyses more concisely. If you experience lack of skills in basic Galaxy operations (tool usage, copy of datasets, etc), do not hesitate to go back and examine the [previous chapters](use_cases_input_data) step by step.


## Input data for Use Case 3-1

As for the previous Use Cases 1 and 2, the first step is to collect all input data in an history that we will name `Input data for Use Case 3-1`. 

1. Create a new history
2. Rename this history `Input data for Use Case 3-1`
3. We are going to upload 40 datasets for the EBI ENA SRP068722. This is a bit tedious, but you follow the instructions bellow, it is not difficult.
    - Use the tool `Extract reads in FASTQ/A format from NCBI SRA`, fill the SRR accession field with the first EBI SRA identifier **SRR3111582** and let the `output format` selected as fastq. Click the `Execute` button.
    - Now that the tool has started to run, you can click on the dataset in order to expand the information. Then click on the rerun icon (two curved arrows, when the mouse pass over the icon you can see the info "Run this job again" popping up). The only thing you have to do is to edit the SRR accession field from **SRR3111582** to **SRR3111583**, and press the  `Execute` button.
    - Just repeat this operation with all the datasets we want to upload. Here is the full list of the SRR identifiers for the 40 datasets we will upload. Note that these identifiers increment by 1 from **SRR3111582** to **SRR3111622** with *one exception*: we go directly from **SRR3111614** to **SRR3111614**.

```
SRR3111582
SRR3111583
SRR3111584
SRR3111585
SRR3111586
SRR3111587
SRR3111588
SRR3111589
SRR3111590
SRR3111591
SRR3111592
SRR3111593
SRR3111594
SRR3111595
SRR3111596
SRR3111597
SRR3111598
SRR3111599
SRR3111600
SRR3111601
SRR3111602
SRR3111603
SRR3111604
SRR3111605
SRR3111606
SRR3111607
SRR3111608
SRR3111609
SRR3111610
SRR3111611
SRR3111612
SRR3111613
SRR3111614
SRR3111616
SRR3111617
SRR3111618
SRR3111619
SRR3111620
SRR3111621
SRR3111622
```

 















