---
title: "Data Access & QC"
permalink: /data-access-qc/
---

![Picture1.png](Project1-DataAccessQC_Gillian_files/Picture1.png)

<div class="alert alert-block alert-info">
    <h1>BIOS 4150/BIOL 6150</h1>
    <h3>Instructor: Dr. King Jordan</h3>
    <p>Nilavrah Sensarma (nsensarma3@gatech.edu),</p> 
    <p>Lasya Pasumarthy (vpasumarthy3@gatech.edu)</p>
</p>
</div>

<div class="alert alert-block alert-warning">
    <h2>Project 1 (Data access & QC) starter notebook</h2>
    <h3>Deadline: 11:59PM, September 22nd, 2025</h3>
</div>

### Rerun the analysis as demonstrated in the lab sessions


---

# **1. Your assigned 1000 Genomes individual**
### *Total Questions: 3*
### *Total Points: 2+2+3 = 7*


---

<div class="alert alert-block alert-warning">
    <h3>1.1 Does you group have an assigned individual ID from the 1000 Genomes Project?</h3>
    <p>Find the individual on the international genomes website and paste a screenshot below</p>
</div>

![image.png](Project1-DataAccessQC_Gillian_files/image.png)

<div class="alert alert-block alert-warning">
    <h3>1.2 What population does you assigned individual belong to?</h3>
    <p>Write about where the geographical location of where this individual lived at the time of recruitment? </p>
</div>

#Write about the individual and their ancestors here.
Japanese in Tokyo, Japan
superpopulation in east asia

<div class="alert alert-block alert-warning">
    <h3>1.3 More about your individual's population?</h3>
    <p>Read more about this population and write about their recent historical  (100-500 years) ancestors.</p>
</div>


```python
#Write about the individual and their ancestors here.
```

# **2. Exome sequencing**
### *Total Questions: 2*
### *Total Points: 3.5+3.5 = 7*


---

<div class="alert alert-block alert-warning">
    <h3>2.1 Does your assigned individual has exome sequencing data? If yes, then what is the SRR ID for that sequencing run?</h3>
    <p>Does you individual have multiple SRR IDs? In case they do, which one would you pick? (Hint: The sequencing information is coming from the same individual; but larger file size might be more preferable)</p>
</div>


Yes,they have 2 SRR ID's.

Using the NHS SRA Browser, I was able to compare the file sizes and found that the bigger files is SRR715907, which is what we will be using.

<div class="alert alert-block alert-warning">
    <h3>2.2 Looking at just the exome sequencing of this individual, how many variants would you expect them to carry compared to the reference genome?</h3>
    <p>Keep in mind this is 1 individual, and learn about the expected number of variants from WGS, WES, and Genotyping</p>
</div>

# **3. Data access**
### *Total Questions: 2*
### *Total Points: 3.5+3.5 = 7*


---

<div class="alert alert-block alert-warning">
    <h3>3.1 Download the exome sequencing data for this individual and write the size of the files</h3>
    <p>3.1.1 How many files do you have? Is the sequencing paired?</p>
    <p>3.1.2 How many reads are present in each file? </p>
</div>


```python
#Show your download command and ls the directory.
!wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19010/sequence_read/SRR715907_1.filt.fastq.gz
!wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19010/sequence_read/SRR715907_2.filt.fastq.gz
    
!mkdir -p ~/Project_1/DataAccess_SRR715907
!ls $HOME/Project_1/
 
!mv SRR715907_1.filt.fastq.gz ~/Project_1/DataAccess_SRR715907/
!mv SRR715907_2.filt.fastq.gz ~/Project_1/DataAccess_SRR715907/
```

    --2025-09-16 21:33:02--  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19010/sequence_read/SRR715907_1.filt.fastq.gz
               => ‘SRR715907_1.filt.fastq.gz’
    Resolving ftp.1000genomes.ebi.ac.uk (ftp.1000genomes.ebi.ac.uk)... 193.62.193.167
    Connecting to ftp.1000genomes.ebi.ac.uk (ftp.1000genomes.ebi.ac.uk)|193.62.193.167|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /vol1/ftp/phase3/data/NA19010/sequence_read ... done.
    ==> SIZE SRR715907_1.filt.fastq.gz ... 2690326249
    ==> PASV ... done.    ==> RETR SRR715907_1.filt.fastq.gz ... done.
    Length: 2690326249 (2.5G) (unauthoritative)
    
    SRR715907_1.filt.fa 100%[===================>]   2.50G  22.4MB/s    in 4m 56s  
    
    2025-09-16 21:38:00 (8.67 MB/s) - ‘SRR715907_1.filt.fastq.gz’ saved [2690326249]
    
    --2025-09-16 21:38:00--  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA19010/sequence_read/SRR715907_2.filt.fastq.gz
               => ‘SRR715907_2.filt.fastq.gz’
    Resolving ftp.1000genomes.ebi.ac.uk (ftp.1000genomes.ebi.ac.uk)... 193.62.193.167
    Connecting to ftp.1000genomes.ebi.ac.uk (ftp.1000genomes.ebi.ac.uk)|193.62.193.167|:21... connected.
    Logging in as anonymous ... Logged in!
    ==> SYST ... done.    ==> PWD ... done.
    ==> TYPE I ... done.  ==> CWD (1) /vol1/ftp/phase3/data/NA19010/sequence_read ... done.
    ==> SIZE SRR715907_2.filt.fastq.gz ... 2697342445
    ==> PASV ... done.    ==> RETR SRR715907_2.filt.fastq.gz ... done.
    Length: 2697342445 (2.5G) (unauthoritative)
    
    SRR715907_2.filt.fa 100%[===================>]   2.51G  7.75MB/s    in 2m 56s  
    
    2025-09-16 21:40:58 (14.6 MB/s) - ‘SRR715907_2.filt.fastq.gz’ saved [2697342445]
    
     DataAccess_SRR715907		    SRR715907_1.filt.fastq.gz   exomeseq1
    'Project1-DataAccessQC (1).ipynb'   SRR715907_2.filt.fastq.gz
    


```python
#Looking at the files
!zcat ~/Project_1/DataAccess_SRR715907/SRR715907_1.filt.fastq.gz | head 
!echo "----"
!zcat ~/Project_1/DataAccess_SRR715907/SRR715907_2.filt.fastq.gz | head 
```

    @SRR715907.2 DB9DDFP1:347:D1HE2ACXX:2:1101:1826:1994/1
    GGGCTGGGCACCCGTGGAATAGGCATGAGGCCAGAAGAGAGTGACAGCGAGCTCCTTGAGGATGAGGAGGATGAAGTGGTAAGAATGAGATCGGAAGAGCA
    +
    @C@DFFFFHDFHGIHJGGHBHGGIGIGGIGDHGDHIJIIIJFAHIIGEHIIBEHGEHFDEFDCECCECDB;?CCC>5>C:>CDCCDCCC>AABB@?B@@?:
    @SRR715907.5 DB9DDFP1:347:D1HE2ACXX:2:1101:2717:1993/1
    TCATGCTCCAATATTGACATCAACATCATGCTACTTGTTGTCTTTGTGGGATTTAACTTGATGTTCACTGGGTTGGTAGTCATCTTTTCCTACATCTACAT
    +
    CCCFFFFFHHHHGJJJJJJJIJJJJJJJJJJJJJJIJHIJGIJJJJFGIJIIIJIJIIJJJJJJJJJJIIJJJJJEEHEFFFFFFEEEEEEDDD;CCDDDC
    @SRR715907.6 DB9DDFP1:347:D1HE2ACXX:2:1101:2607:2000/1
    CTCTTCTGAATTTGTTTACAGAAAATCAAAATGACCATGAGGAGGAGGAGGGGAAAGCGCCAGTGCCCCCCAGGGAACTGTGGATTTGTGGGCTGTTAGTT
    
    gzip: stdout: Broken pipe
    ----
    @SRR715907.2 DB9DDFP1:347:D1HE2ACXX:2:1101:1826:1994/2
    CATTCTTACCACTTCATCCTCCTCATCCTCAAGGAGCTCGCTGTCACTCTCTTCTGGCCTCATGCCTATTCCACGGGTGCCCAGCCCAGATCGGAAGAGCG
    +
    @CCFFFFFHAHHHJJJIIJGIGIIJJJJJGHEHGGGGGIJIIJIHGIIIJFHIJJDHGHIGGIJGIEHGGEFHHEDDA>@DDCDDDB?B<CC@D?@BBACB
    @SRR715907.5 DB9DDFP1:347:D1HE2ACXX:2:1101:2717:1993/2
    CCTGCACTAGAAGACATTTTCAGGATGGTGGCCATGATGTAGATGTAGGAAAAGATGACTACCAACCCAGTGAACATCAAGTTAAATCCCACAAAGACAAC
    +
    CCCFFFFFHHHHGJJJJJJJJJJJJJJJHIJJJJJIIJJFII4BGEIIJJEFHIIIJJJIIIIJJJJJJIIIHHHHGHFFFFFFFEDCCC@ABDDDDDDD<
    @SRR715907.6 DB9DDFP1:347:D1HE2ACXX:2:1101:2607:2000/2
    CCTTGATCCGTGGTGTCCAGATGTCACTATTGAACTAACAGCCCACAAATCCACAGTTACCTGGGGGGCACTGGCGCTTTCCCCTCCTCCTCCTCATGGTC
    
    gzip: stdout: Broken pipe
    

#FastQ file1
@SRR715907.2 DB9DDFP1
#unique instrumental name
:347
#Run ID
:D1HE2ACXX
#Flow cell ID?
:2
#number of flow cells
:1101
#Title number within flow cell lane
:1826
#x-coordinate
:1994/1
#Y-coordinate and member of pair

#Raw sequence letters
GGGCTGGGCACCCGTGGAATAGGCATGAGGCCAGAAGAGAGTGACAGCGAGCTCCTTGAGGATGAGGAGGATGAAGTGGTAAGAATGAGATCGGAAGAGCA
+
#Quality values
@C@DFFFFHDFHGIHJGGHBHGGIGIGGIGDHGDHIJIIIJFAHIIGEHIIBEHGEHFDEFDCECCECDB;?CCC>5>C:>CDCCDCCC>AABB@?B@@?:
#more identifiers and raw sequence letters
@SRR715907.5 DB9DDFP1:347:D1HE2ACXX:2:1101:2717:1993/1
TCATGCTCCAATATTGACATCAACATCATGCTACTTGTTGTCTTTGTGGGATTTAACTTGATGTTCACTGGGTTGGTAGTCATCTTTTCCTACATCTACAT
+
CCCFFFFFHHHHGJJJJJJJIJJJJJJJJJJJJJJIJHIJGIJJJJFGIJIIIJIJIIJJJJJJJJJJIIJJJJJEEHEFFFFFFEEEEEEDDD;CCDDDC
@SRR715907.6 DB9DDFP1:347:D1HE2ACXX:2:1101:2607:2000/1
CTCTTCTGAATTTGTTTACAGAAAATCAAAATGACCATGAGGAGGAGGAGGGGAAAGCGCCAGTGCCCCCCAGGGAACTGTGGATTTGTGGGCTGTTAGTT


#FastQ file2
#Sequence identifiers
@SRR715907.2 DB9DDFP1:347:D1HE2ACXX:2:1101:1826:1994/2
#Raw sequence data
CATTCTTACCACTTCATCCTCCTCATCCTCAAGGAGCTCGCTGTCACTCTCTTCTGGCCTCATGCCTATTCCACGGGTGCCCAGCCCAGATCGGAAGAGCG
+
@CCFFFFFHAHHHJJJIIJGIGIIJJJJJGHEHGGGGGIJIIJIHGIIIJFHIJJDHGHIGGIJGIEHGGEFHHEDDA>@DDCDDDB?B<CC@D?@BBACB
@SRR715907.5 DB9DDFP1:347:D1HE2ACXX:2:1101:2717:1993/2
CCTGCACTAGAAGACATTTTCAGGATGGTGGCCATGATGTAGATGTAGGAAAAGATGACTACCAACCCAGTGAACATCAAGTTAAATCCCACAAAGACAAC
+
CCCFFFFFHHHHGJJJJJJJJJJJJJJJHIJJJJJIIJJFII4BGEIIJJEFHIIIJJJIIIIJJJJJJIIIHHHHGHFFFFFFFEDCCC@ABDDDDDDD<
@SRR715907.6 DB9DDFP1:347:D1HE2ACXX:2:1101:2607:2000/2
CCTTGATCCGTGGTGTCCAGATGTCACTATTGAACTAACAGCCCACAAATCCACAGTTACCTGGGGGGCACTGGCGCTTTCCCCTCCTCCTCCTCATGGTC

# **4. Pre-QC with FastQC**
### *Total Questions: 3*
### *Total Points: 4+4+6 = 14*


---

<div class="alert alert-block alert-warning">
    <h3>4.1 Show your exact command below.</h3>
</div>


```python
#Create a directory for QC results.
!mkdir -p ~/Project_1/DataAccess_SRR715907/fastQCBeforeTrimming
!mkdir -p ~/Project_1/DataAccess_SRR715907/Trimming
!mkdir -p ~/Project_1/DataAccess_SRR715907/fastQCAfterTrimming
```


```python
#Run fastqc.
!fastqc -o ~/Project_1/DataAccess_SRR715907/fastQCBeforeTrimming/ ~/Project_1/DataAccess_SRR715907/SRR715907_1.filt.fastq.gz ~/Project_1/DataAccess_SRR715907/SRR715907_2.filt.fastq.gz

```

    application/gzip
    application/gzip
    Started analysis of SRR715907_1.filt.fastq.gz
    Approx 5% complete for SRR715907_1.filt.fastq.gz
    Approx 10% complete for SRR715907_1.filt.fastq.gz
    Approx 15% complete for SRR715907_1.filt.fastq.gz
    Approx 20% complete for SRR715907_1.filt.fastq.gz
    Approx 25% complete for SRR715907_1.filt.fastq.gz
    Approx 30% complete for SRR715907_1.filt.fastq.gz
    Approx 35% complete for SRR715907_1.filt.fastq.gz
    Approx 40% complete for SRR715907_1.filt.fastq.gz
    Approx 45% complete for SRR715907_1.filt.fastq.gz
    Approx 50% complete for SRR715907_1.filt.fastq.gz
    Approx 55% complete for SRR715907_1.filt.fastq.gz
    Approx 60% complete for SRR715907_1.filt.fastq.gz
    Approx 65% complete for SRR715907_1.filt.fastq.gz
    Approx 70% complete for SRR715907_1.filt.fastq.gz
    Approx 75% complete for SRR715907_1.filt.fastq.gz
    Approx 80% complete for SRR715907_1.filt.fastq.gz
    Approx 85% complete for SRR715907_1.filt.fastq.gz
    Approx 90% complete for SRR715907_1.filt.fastq.gz
    Approx 95% complete for SRR715907_1.filt.fastq.gz
    Analysis complete for SRR715907_1.filt.fastq.gz
    Started analysis of SRR715907_2.filt.fastq.gz
    Approx 5% complete for SRR715907_2.filt.fastq.gz
    Approx 10% complete for SRR715907_2.filt.fastq.gz
    Approx 15% complete for SRR715907_2.filt.fastq.gz
    Approx 20% complete for SRR715907_2.filt.fastq.gz
    Approx 25% complete for SRR715907_2.filt.fastq.gz
    Approx 30% complete for SRR715907_2.filt.fastq.gz
    Approx 35% complete for SRR715907_2.filt.fastq.gz
    Approx 40% complete for SRR715907_2.filt.fastq.gz
    Approx 45% complete for SRR715907_2.filt.fastq.gz
    Approx 50% complete for SRR715907_2.filt.fastq.gz
    Approx 55% complete for SRR715907_2.filt.fastq.gz
    Approx 60% complete for SRR715907_2.filt.fastq.gz
    Approx 65% complete for SRR715907_2.filt.fastq.gz
    Approx 70% complete for SRR715907_2.filt.fastq.gz
    Approx 75% complete for SRR715907_2.filt.fastq.gz
    Approx 80% complete for SRR715907_2.filt.fastq.gz
    Approx 85% complete for SRR715907_2.filt.fastq.gz
    Approx 90% complete for SRR715907_2.filt.fastq.gz
    Approx 95% complete for SRR715907_2.filt.fastq.gz
    Analysis complete for SRR715907_2.filt.fastq.gz
    

<div class="alert alert-block alert-warning">
    <h3>4.2 Show the <i>per base sequencing quality</i></h3>
    <p>Add a screenshot below</p>
</div>

file 1
![image.png](Project1-DataAccessQC_Gillian_files/image.png)

file 2
![image.png](Project1-DataAccessQC_Gillian_files/image.png)

<div class="alert alert-block alert-warning">
    <h3>4.3 Show the <i>Per base sequence content</i></h3>
    <p>Add a screenshot below</p>
</div>

file 1
![image.png](Project1-DataAccessQC_Gillian_files/image.png)

file 2
![image.png](Project1-DataAccessQC_Gillian_files/image.png)

# **5. Trimming with fastp**
### *Total Questions: 3*
### *Total Points: 7+7+7 = 21*


---

<div class="alert alert-block alert-warning">
    <h3>5.1 What was the tool and command you used?</h3>
    <p>Write the tool name and exact command below</p>
</div>


```python
!fastp -i ~/Project_1/DataAccess_SRR715907/SRR715907_1.filt.fastq.gz -I ~/Project_1/DataAccess_SRR715907/SRR715907_2.filt.fastq.gz -o ~/Project_1/DataAccess_SRR715907/Trimming/SRR715907_1.Trimmed.fastq.gz -O ~/Project_1/DataAccess_SRR715907/Trimming/SRR715907_2.Trimmed.fastq.gz -f 12 -t 15
#fastp 
```

    Read1 before filtering:
    total reads: 30296516
    total bases: 3059948116
    Q20 bases: 2960909567(96.7634%)
    Q30 bases: 2784146983(90.9867%)
    
    Read2 before filtering:
    total reads: 30296516
    total bases: 3059948116
    Q20 bases: 2911773155(95.1576%)
    Q30 bases: 2733102015(89.3186%)
    
    Read1 after filtering:
    total reads: 29200891
    total bases: 2159991968
    Q20 bases: 2126056838(98.4289%)
    Q30 bases: 2017860657(93.4198%)
    
    Read2 after filtering:
    total reads: 29200891
    total bases: 2159991968
    Q20 bases: 2118495003(98.0788%)
    Q30 bases: 2006685595(92.9025%)
    
    Filtering result:
    reads passed filter: 58401782
    reads failed due to low quality: 2186604
    reads failed due to too many N: 4646
    reads failed due to too short: 0
    reads with adapter trimmed: 1265966
    bases trimmed due to adapters: 1780482
    
    Duplication rate: 3.00044%
    
    Insert size peak (evaluated by paired-end reads): 135
    
    JSON report: fastp.json
    HTML report: fastp.html
    
    fastp -i /home/hice1/gmyers30/Project_1/DataAccess_SRR715907/SRR715907_1.filt.fastq.gz -I /home/hice1/gmyers30/Project_1/DataAccess_SRR715907/SRR715907_2.filt.fastq.gz -o /home/hice1/gmyers30/Project_1/DataAccess_SRR715907/Trimming/SRR715907_1.Trimmed.fastq.gz -O /home/hice1/gmyers30/Project_1/DataAccess_SRR715907/Trimming/SRR715907_2.Trimmed.fastq.gz -f 12 -t 15 
    fastp v0.23.4, time used: 278 seconds
    


```python
#Exact command here
#!fastqc -o ~/Project_1/DataAccess_SRR715907/fastQCAfterTrimming/ ~/Project_1/DataAccess_SRR715907/Trimming/SRR715907_1.Trimmed.fastq.gz ~/Project_1/DataAccess_SRR715907/Trimming/SRR715907_2.Trimmed.fastq.gz
#Accidentally added this earlier, this is also in the next section where fastQC is ran
```

    application/gzip
    application/gzip
    Started analysis of SRR715907_1.Trimmed.fastq.gz
    Approx 5% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 10% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 15% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 20% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 25% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 30% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 35% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 40% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 45% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 50% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 55% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 60% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 65% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 70% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 75% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 80% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 85% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 90% complete for SRR715907_1.Trimmed.fastq.gz
    Approx 95% complete for SRR715907_1.Trimmed.fastq.gz
    Analysis complete for SRR715907_1.Trimmed.fastq.gz
    Started analysis of SRR715907_2.Trimmed.fastq.gz
    Approx 5% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 10% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 15% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 20% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 25% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 30% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 35% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 40% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 45% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 50% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 55% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 60% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 65% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 70% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 75% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 80% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 85% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 90% complete for SRR715907_2.Trimmed.fastq.gz
    Approx 95% complete for SRR715907_2.Trimmed.fastq.gz
    Analysis complete for SRR715907_2.Trimmed.fastq.gz
    

<div class="alert alert-block alert-warning">
    <h3>5.2 How did you perform the trimming?</h3>
    <p>Did you use dynamic built-in features or hard coded parameters?</p>
    <p>Did you have to go back and repeat the trimming after run(s) for <b>postQC?</b></p>
</div>

#Write here

<div class="alert alert-block alert-warning">
    <h3>5.3 What percentage of reads did you lose in the trimming process?</h3>
    <p>5.3.1 Write the total number of reads before and after trimming</p>
    <p>5.3.2 If you had a paired fastq file, are there any reads after trimming present in one file and not the other? (eg: present in reverse but not in forward?)</p>
    
</div>

#Before trimming_1:30296516
#Before trimming_2:30296516


#After trimming_1:29200891
#After trimming_2:29200891

#it looks like there is the same amount of reads 

# **6. Post-QC with FastQC**
### *Total Questions: 3*
### *Total Points: 4+6+4 = 14*


---

<div class="alert alert-block alert-warning">
    <h3>6.1 Show your exact command(s) below.</h3>
</div>


```python
#Command here
!fastqc -o ~/Project_1/DataAccess_SRR715907/fastQCAfterTrimming/ ~/Project_1/DataAccess_SRR715907/Trimming/SRR715907_1.Trimmed.fastq.gz ~/Project_1/DataAccess_SRR715907/Trimming/SRR715907_2.Trimmed.fastq.gz

```

<div class="alert alert-block alert-warning">
    <h3>6.2 Show <i>per base sequencing quality</h3>
    <p>If the quality was falling below 20 for certain regions then trim it. Add a screenshot below</p>
</div>

file_1
![image.png](Project1-DataAccessQC_Gillian_files/image.png)

file_2
![image.png](Project1-DataAccessQC_Gillian_files/image.png)

<div class="alert alert-block alert-warning">
    <h3>6.3 Show the <i>Per base sequence content</i></h3>
    <p>Ideally we shouldn't see adapters at the beginning of the reads in the fastq files. Add a screenshot below</p>
</div>

File_1
![image.png](Project1-DataAccessQC_Gillian_files/image.png)

File_2
![image.png](Project1-DataAccessQC_Gillian_files/image.png)

