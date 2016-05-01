# Histories for Use Cases 1-1, 1-2, 1-3 and 1-4

As you will see, Histories 1-1, 1-2 and 1-3 are generated in the same way, using their corresponding workflows. These workflows are available in your Galaxy top menu. An important thing to remember is that you will always start **from** the `Input data for Use Cases 1-1, 1-2, 1-3 and 1-4` history, run the appropriate workflow, **sending the outputs of the workflow in a new history** named accordingly.


## History for Use Case 1-1.

#### 1. As aforementioned, ensure that you are in the `Input data for Use Cases 1-1, 1-2, 1-3 and 1-4` history.
You can always control this by using the top menu **Users** --> **Saved History** and selecting the desired history. If you don't see the History right bar, just click in addition the top menu **Analyze Data**

#### 2. Select the appropriate workflow
- Click now on the **Workflow** top menu
- Select the workflow **"Metavisitor: Workflow for Use Case 1-1 (imported from API)"** and to see the workflow, select the submenu **"Edit"**
- Now that you see the workflow, you can directly execute it by clicking the top right wheel icon and selecting **"Run"**
    
    ![copydataset](images/runworkflow.png)
    
    
#### 3. Select the appropriate parameters before running the workflow

You now see a page with all the workflow steps, whose top part looks like:
    
![copydataset](images/workflow1-1.png)

- As pointed by the first red arrow, a parameter has to be provided at runtime of the workflow: the **ncbi_guide_ID**. In this Use Case as in the other 1-2 and 1-3 Use Cases, you will paste in the **ncbi_guide_ID** field the `NC_007919.3_`value. This is the NCBI identifier for the Nora virus genome sequence which will be retrieved from Genbank during the workflow and used as a guide for the final reconstruction of the Nora virus genome sequence that is "present" in the analyzed small RNA sequencing datasets.

- You have to select an Input dataset collection for Step 1 (second red arrow). However, as there is only one dataset collection in the input history (the one we have prepared in the [previous chapter](use_cases_input_data/#history-with-input-data-for-use-cases-1-1-1-2-1-3-and-1-4)), there is no other option in the menu than "**SRP013822**".
    
- You have to select the viral nucleotide Blast database for Step 2. Here again there is indeed nothing else to select than the "**nucleotide vir1 blast database**", just because there is only one dataset in the input history with the "blast database" type.

- You can review the other steps of the workflow. But there is no other selection to perform before running the workflow.

#### 4. Running the workflow **_sending the outputs in a new history_**

We are almost ready, but before clicking the "**Run Workflow**" button there is an important thing to do:
- Check the "**Send results to a new history**" checkbox as shown Here

![sendtonewhistory](images/new_workflow_history.png)

- And edit the field to "History for Use Case 1-1"
- **You can now click the "Run workflow" button.**
This trigger the workflow run. After a few seconds (may be take a while for complex workflows), you will see an alert that the workflow is started, and a link to navigate to the newly created history.