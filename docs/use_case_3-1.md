Now that you get more familiar with manipulations in Galaxy with the Use Cases 1-1 to 1-4 described in details in the previous chapters, we will describe the other Use Case analyses more concisely. If you experience lack of skills in basic Galaxy operations (tool usage, copy of datasets, etc), do not hesitate to go back and examine the [previous chapters](use_cases_input_data) step by step.


## Input data for Use Case 3-1

As for the previous Use Cases 1 and 2, the first step is to collect all input data in an history that we will name `Input data for Use Case 3-1`. 

1. Create a new history
- Rename this history `Input data for Use Case 3-1`
- We are going to upload 40 datasets for the EBI ENA SRP068722. This is a bit tedious, but you follow the instructions bellow, it is not difficult.
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

4. For each of the 40 imported datasets, click the pencil icon, and change the datatype to `fastqsanger`. Be systematic, do it for all datasets 1 to 40.
5. Click on the checked box icon in the history top menu, check all 40 datasets (`All` button), and `For all selected`, `Build a dataset list` that you name "SRP068722".
6. Copy the `vir1 nucleotide BLAST database` from the `References` history to the current history `Input data for Use Case 3-1`.
7. Now we still have to associate sequencing dataset coming from a same patient. We are going to use the tool `Concatenate multiple datasets` to merge multiple datasets in a same fastq file.
    - For patient 0450-318, use `Concatenate multiple datasets` and select the datasets SRR3111582 to SRR3111587. Run the tool and rename the dataset "patient 0450-318"
    - For patient 0387-272, use `Concatenate multiple datasets` and select the datasets SRR3111588 to SRR3111593. Run the tool and rename the dataset "patient 0387-272"
    - For patient 0629-453, use `Concatenate multiple datasets` and select the datasets SRR3111594 to SRR3111599. Run the tool and rename the dataset "patient 0629-453"
    - For patient 0444-312, use `Concatenate multiple datasets` and select the datasets SRR3111600 to SRR3111603. Run the tool and rename the dataset "patient 0444-312"
    - For patient 0500-355neg, use `Concatenate multiple datasets` and select the datasets SRR3111604 and SRR3111605. Run the tool and rename the dataset "patient 0500-355neg"
    - For patient 0292-xxxneg, use `Concatenate multiple datasets` and select the datasets SRR3111606 and SRR3111607. Run the tool and rename the dataset "patient 0292-xxxneg"
    - For patient 0394-274, use `Concatenate multiple datasets` and select the datasets SRR3111608 and SRR3111609. Run the tool and rename the dataset "patient 0394-274"
     - For patient 0218-162neg, use `Concatenate multiple datasets` and select the datasets SRR3111610 and SRR3111611. Run the tool and rename the dataset "patient 0218-162neg"
     - For patient 0311-217HIVneg, use `Concatenate multiple datasets` and select the datasets SRR3111612 and SRR3111613. Run the tool and rename the dataset "patient 0311-217HIVneg"
     - For patient 0440-307neg, use `Concatenate multiple datasets` and select the datasets SRR3111614 and SRR3111616. Run the tool and rename the dataset "patient 0440-307neg"
     - For patient 0518-370neg, use `Concatenate multiple datasets` and select the datasets SRR3111617 and SRR3111618. Run the tool and rename the dataset "patient 0518-370neg"
     - For patient 0560-420neg, use `Concatenate multiple datasets` and select the datasets SRR3111619 and SRR3111620. Run the tool and rename the dataset "patient 0560-420neg"
     - For patient 0575-419neg, use `Concatenate multiple datasets` and select the datasets SRR3111621 and SRR3111622. Run the tool and rename the dataset "patient 0575-419neg"
8. The last action to perform in this history is to create a dataset collection of patient datasets: Click on the checked box icon in the history top menu, check the "patient... " datasets we have just generated by concatenation (13 datasets), and `For all selected`, `Build a dataset list` that you name "patient collection".
9. We are done.

## History for Use Case 3-1
1. Stay in the history `Input data for Use Case 3-1`
- pick the workflow `Metavisitor: Workflow for Use Case 3-1` in the workflows menu, and select the `run` option.
- For Step 1 (Fever Patient Sequences collection), select `patient collection` (this should be already selected).
- For Step 2, select the `nucleotide vir1 blast database` (this should also be already selected)
- As usual, check the box `Send results to a new history`, edit the name of the new history to `History for Use Case 3-1`, and `Execute` the workflow ! Note, that for complex workflows with dataset collections in input, the actual warning that the workflow is started make take time to show up.
















