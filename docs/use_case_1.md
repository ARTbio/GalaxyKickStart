# Histories for Use Cases 1-1, 1-2, 1-3 and 1-4

As you will see, Histories 1-1, 1-2 and 1-3 are generated in the same way, using their corresponding workflows. These workflows are available in your Galaxy top menu. An important thing to remember is that you will always start **from** the `Input data for Use Cases 1-1, 1-2, 1-3 and 1-4` history, run the appropriate workflow, **sending the outputs of the workflow in a new history** named accordingly.


## History for Use Case 1-1.
1. As aforementioned, ensure that you are in the `Input data for Use Cases 1-1, 1-2, 1-3 and 1-4` history.

You can always control this by using the top menu **Users** --> **Saved History** and selecting the desired history. If you don't see the History right bar, just click in addition the top menu **Analyze Data**

2. Select the appropriate workflow
    - Click now on the **Workflow** top menu
    - Select the workflow **"Metavisitor: Workflow for Use Case 1-1 (imported from API)"** and to see the workflow, select the submenu **"Edit"**
    - Now that you see the workflow, you can directly execute it by clicking the top right wheel icon and selecting **"Run"**
    
    ![copydataset](images/runworkflow.png)
    
    - which will show this page whose upper part is shown here:
    
    ![copydataset](images/workflow1-1.png)
    
    - A parameter has to be provided at runtime of the workflow: the **ncbi_guide_ID**. In this Use Case as in the other 1-2 and 1-3 Use Cases, you will paste in the **ncbi_guide_ID** field the `NC_007919.3_`value. This is the NCBI identifier for the Nora virus genome sequence which will be retrieved from Genbank during the workflow and used as a guide for the final reconstruction of the Nora virus genome sequence that is "present" in the analyzed small RNA sequencing datasets.
    
    - for the 7 datasets retrieved from EBI SRA, change the datatype **fastq** to **fastqsanger**:
    click on the pencil icon of the dataset, click the tab Datatype, and select fastqsanger in the New Type menu
3. Create a dataset collection **SRP013822**
    - Click on the checked box icon at the top of history bar as indicated below
    
    ![checkbow](images/check_datasets.png)
    
    - Select the 7 datasets
    - Select **Build Dataset List** in the menu **For all selected...**
    
    ![datasetlist](images/dataset_list.png)
    
    - and type **SRP013822** in the Name field and click `Create list`
    - you can leave the checked datasets view by clicking again the check box in the history top menu
    
4. copy the vir1 blast database that we have prepared earlier in the [Reference](metavisitor_configure_references.md#3-prepare-blast-databases) history.
    - To do so, click on the little wheel icon in the history top menu (in the history right bar).
    
    ![copydataset](images/copydataset.png)
    
    - Select "Copy Datasets"
    - In the open page, select "References" in the Source History menu, check the "nucleotide vir1 blast database" dataset; select "Input data for Use Case 1_1, ..."; and click the "Copy History Items".
    - if you refresh the history, you will see the "nucleotide vir1 blast database" dataset showing up.
    
That is all for the moment. We will add latter datasets in this history `Input data for Use Cases 1-1, 1-2, 1-3 and 1-4`. However, these datasets do no exist yet: this will be produced by the Use Cases 1-1, 1-2, 1-3 workflows !