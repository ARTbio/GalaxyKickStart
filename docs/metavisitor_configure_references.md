Once you know how to access to your Metavisitor Galaxy instance with a web browser and are able to perform basic start/stop/restart operations, there is still some work needed to import and configure reference data (genomes) so that they are directly available to all instance users for running tools and workflows

Here we provide the step-by-step description of what we *actually did ourselves* to prepare our Metavisitor instance before performing the analyses described [here](http://dx.doi.org/10.1101/048983).

## 1. Connect to your Metavisitor Galaxy admin account with your web browser

## 2. Import reference data in an history "References"

At first, you need to import and prepare the reference datasets you will need for most of the Metavisitor analyses. As a Galaxy admin you will make latter some of these references directly accessible to the Galaxy tools, and/or accessible to any other users by putting them in a Galaxy public library.

#### a. Preliminary actions
click on the `Analyze Data` menu
rename the `Unnamed history` to `References`

#### b. Upload nucleotide vir1 fasta file

```
Click on the small arrow icon at the top of the tool bar (left handside of the Galaxy interface)
In the open window, click on the Paste/Fetch data button
Paste the URL https://ndownloader.figshare.com/files/4949173
Click the start button.
```
A first dataset will show up in your history, first in grey (the job is starting), then yellow (the job - upload - is running), and eventually green (job is successfully done).

When finished, rename the dataset 1 `https://ndownloader.figshare.com/files/4949173` to `nucleotide vir1` for clarity (using the small pencil icon).

#### c. Upload protein vir1 fasta file
Repeat the same operation as in `b.` using the url `https://ndownloader.figshare.com/files/4949170`.
When finished, rename the dataset 2 `https://ndownloader.figshare.com/files/4949170` to `protein vir1` for clarity.


#### d. Upload the Drosophila melanogaster release 6 genome
Repeat the same operation as in `b.` using the url `ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.10_FB2016_02/fasta/dmel-all-chromosome-r6.10.fasta.gz`.
When finished, rename the dataset 3 `ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.10_FB2016_02/fasta/dmel-all-chromosome-r6.10.fasta.gz` to `dm6` for clarity.

#### e. Upload the Anopheles gambiae release P4
Repeat the same operation as in `b.` using the url `https://www.vectorbase.org/sites/default/files/ftp/downloads/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.gz`.
When finished, rename the dataset 4 `https://www.vectorbase.org/sites/default/files/ftp/downloads/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.gz` to `AgamP4` for clarity.

#### f. Upload the Plasmodium berghei genome
Repeat the same operation as in `b.` using the url `ftp://ftp.ensemblgenomes.org/pub/release-28/protists/fasta/plasmodium_berghei/dna/Plasmodium_berghei.May_2010.28.dna_sm.genome.fa.gz`.
When finished, rename the dataset 5 `ftp://ftp.ensemblgenomes.org/pub/release-28/protists/fasta/plasmodium_berghei/dna/Plasmodium_berghei.May_2010.28.dna_sm.genome.fa.gz` to `P. berghei` for clarity.

#### g. Upload the Homo sapiens release GRCh38/hg19
Repeat the same operation as in `b.` using the url `ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`.
When finished, rename the dataset 6 `ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` to `hg19` for clarity.

## 3. Prepare Blast databases
#### a. Nucleotide vir1 blast database
Use the tool `NCBI BLAST+ makeblastdb`, check the radio button "nucleotide", select the dataset 1 (nucleotide vir1), give "nucleotide vir1 blastdb" as a "Title for BLAST database", leave the rest of the tool form unchanged and click "Execute" button.

Rename the generated dataset 7 "nucleotide BLAST database from data 1" to "nucleotide vir1 blast database" for clarity

#### b. Protein vir1 blast database
Use the tool `NCBI BLAST+ makeblastdb`, check the radio button "protein", select the dataset 2 (protein vir1), give "protein vir1 blastdb" as a "Title for BLAST database", leave the rest of the tool form unchanged and click "Execute" button.

Rename the generated dataset 8 "protein BLAST database from data 2" to "protein vir1 blast database" for clarity

## 4. Creating Galaxy dbkey and fasta references accessible to tools for every user
Here we are going in the `admin` panel, click `Local data` in the left menu and select the `Create DBKey and Reference Genome` in the "**Run Data Manager Tools**" (last line of the top section).

#### a. nucleotide vir1
in the newly open browser window (from the last click on `Create DBKey and Reference Genome`)
- select "New" for `Use existing dbkey or create a new one`
- enter "vir1" in the `dbkey`field
- leave "Display name for dbkey", `Name of sequence` and `ID for sequence` empty !
- select `history` in the `Choose the source for the reference genome` menu
- select "nucleotide vir1" in the `FASTA File` menu
- let `Sort by chromosome name` selected on `As is`
- press the `execute`button !

#### b. dm6
Repeat the operation described in a., but this time

- put "dm6" for the `dbkey` field
- select "3: dm6" in the `FASTA File` menu

Be sure that the `References` history is selected in the background, otherwise the uploaded genomes will not be available in this menu.

#### c. AgamP4
Repeat the operation described in a., but this time

- put "AgamP4" for the `dbkey` field
- select "4: AgamP4" in the `FASTA File` menu

Be sure that the `References` history is selected in the background, otherwise the uploaded genomes will not be available in this menu.

#### d. hg19
Repeat the operation described in a., but this time

- put "hg19" for the `dbkey` field
- select "6: hg19" in the `FASTA File` menu

Be sure that the `References` history is selected in the background, otherwise the uploaded genomes will not be available in this menu.

## 5. Creating Galaxy bowtie indexes accessible to tools for every user
Now we are going to generate the bowtie indexes using another data manager tool.
But before doing this, we have to perform a low level Galaxy admin task: restart the Galaxy server instance so that the dbkey and the fasta genomes that we've just created for the server are registered and seen by the tools.

Depending on your skill level, you have the simple, dirty way:

- reboot the machine where the galaxy server is running !

or the clean, freaking (for non-Geek normal biologists) way:

- Connect to the server where the Galaxy instance has been installed either through an ssh connection or using your local terminal. and type

`sudo supervisorctl restart galaxy:`

If everything went fine you should see in your terminal
```
# supervisorctl restart galaxy:
galaxy_web: stopped
handler0: stopped
handler1: stopped
handler0: started
handler1: started
galaxy_web: started
```
 freaking insn't it ?

----
#### a. vir 1 bowtie index
Now, Go back to your web browser and the `admin` panel, click again `Local data` in the left menu and  select this time the `Bowtie index - builder` in the "**Run Data Manager Tools**" (top section).

select "vir1" in the `Source FASTA Sequence` menu of the Bowtie index builder tool form, leave the other options empty, and click execute.

#### b. dm6 bowtie index
Go back to your web browser and the `admin` panel, click again `Local data` in the left menu and  select this time the `Bowtie index - builder` in the "**Run Data Manager Tools**" (top section).

select "dm6" in the `Source FASTA Sequence` menu of the Bowtie index builder tool form, leave the other options empty, and click execute.

#### c. AgamP4 bowtie index
Go back to your web browser and the `admin` panel, click again `Local data` in the left menu and  select this time the `Bowtie index - builder` in the "**Run Data Manager Tools**" (top section).

select "AgamP4" in the `Source FASTA Sequence` menu of the Bowtie index builder tool form, leave the other options empty, and click execute.

#### d. hg19 bowtie index
Go back to your web browser and the `admin` panel, click again `Local data` in the left menu and  select this time the `Bowtie index - builder` in the "**Run Data Manager Tools**" (top section).

select "hg19" in the `Source FASTA Sequence` menu of the Bowtie index builder tool form, leave the other options empty, and click execute.

----
#### Note that the preparation of bowtie indexes can be long ! (several hours for the vir1 bowtie index for instance)


** When bowtie indexes are ready (green in the `Data Manager History (automatically created)`) restart the Galaxy server instance as explained above, so that the bowtie indexes that we've just created for the server are registered and seen by the tools.**

## 6. Creating Galaxy bowtie2 indexes accessible to tools for every user
Finally, we are going to generate the bowtie2 indexes using another data manager tool.
If not done before, restart the Galaxy instance as explained above.

----
#### a. vir 1 bowtie2 index
Now, Go back to your web browser and the `admin` panel, click again `Local data` in the left menu and  select this time the `Bowtie2 index - builder` in the "**Run Data Manager Tools**" (top section).

select "vir1" in the `Source FASTA Sequence` menu of the Bowtie index builder tool form, leave the other options empty, and click execute.

#### b. AgamP4 bowtie2 index
Go back to your web browser and the `admin` panel, click again `Local data` in the left menu and  select this time the `Bowtie2 index - builder` in the "**Run Data Manager Tools**" (top section).

select "AgamP4" in the `Source FASTA Sequence` menu of the Bowtie index builder tool form, leave the other options empty, and click execute.

#### c. hg19 bowtie2 index
Go back to your web browser and the `admin` panel, click again `Local data` in the left menu and  select this time the `Bowtie2 index - builder` in the "**Run Data Manager Tools**" (top section).

select "hg19" in the `Source FASTA Sequence` menu of the Bowtie index builder tool form, leave the other options empty, and click execute.

----
#### Note that the preparation of bowtie2 indexes can be long too ! (several hours for the vir1 bowtie2 index for instance)


** When bowtie2 indexes are ready (green in the `Data Manager History (automatically created)`) restart the Galaxy server instance as explained above, so that the bowtie indexes that we've just created for the server are registered and seen by the tools.**

