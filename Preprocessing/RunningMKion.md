# Instructions to run MKion on DOSBox

*Alex Brown*

*2022-09-25*

## Prepare FASTA files for MKion

Convert FASTQ to FASTA:

```
seqtk seq -a input.fastq > output.fasta
```

Reformat FASTA sequences to be in one line: 

````perl
perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' in.fasta > out.fa
````

Convert file to DOS format:

```perl
perl -pi -e 's/\n/\r\n/' input.fa
```



## Run in Boxer

On Mac Pro the MKion files live in the following location:
`/Users/alexbrown/Documents/Virtual Machines/`

Run using boxer. Click the DOS box and **drag and drop the directory you want to use**.

use `dir`

Type `mkionoc.exe`

Click the magnifying glass 🔍 icon and adjust settings to max CPU. This will allow you to run these files 2-3x faster than in DOSBox.

To run MKion, place the fasta file, and Valpha and Jalpha identifier files in the same folder. Execute the program and input the identifier files as instructed. MKion analyzes the data in several steps by menu. 

First, input the fasta file and choose **option 1** to clean up the files. If the sequence is not the positive strand, it will be converted to its complement.

In our case, the sequences are on the negative strand. Answer: **[N]**	

<u>Short and long thresholds:</u>

* 100

* 800

Change the name to something new, and give it a different name than before. Don't include .fa this will be added when the file is made. Causes a runtime error if you add it. Ex: `555C`

Wait for the program to continue (will look unresponsive for about 30sec)

Define a new name for the header. This is important. change the name of the header to something like `555JA`  **option 2**. Doing anything otherwise may screw up the analysis so that the CDR3 is not correct. It May take up to an hour to finish this step



<u>Return to the main menu</u>

select option 2. Run file from before ex: `555C` and output as `555OUT`

This time when it asks if you want to change the header name, just use the previous input for the header.

<u>Return to the main menu</u>

Export files as excel readable files (**option 4**). This outputs a .TXT file which is what you want to keep.

> Note: There is an issue with how the J column is organized. it seems that MKion is outputting these with a space or extra characters at the end. These need to be removed for the python program to run correctly. When these files are run through R the spaces are no longer present and the files can be read correctly by the python program.





## Collect Files





___

TEST with a different file that may have more sequences...

Test with a file that is set up for DOS Test with N, 



 You define the shortest and longest sequence lengths allowed to eliminate worthless sequences and those that might crash the program. Generally, pick 170bp and 600bp as the limit unless you know better. The program will then parse the sequences into 3 files `SHORTSeq.fa`,` LONGSeq.fa`and a file that you name for the cleaned-up sequences. The first two files are rewritten every time the program is run, so rename them if you want to save them.  

The fasta file with the cleaned-up sequences will have a new header for each sequence with:

1. A unique sequence identifier.

2. Sample identifier,

3. Space for the TRAV number

4. Space for the TRAJ number

5. Space for the functional state of the sequence (P=productive, O=CDR3 out of frame, T=termination codon in the CDR3 and U= some portion of the sequence not identified.

6. The CDR3 protein sequence.

Next input the file of cleaned-up sequences and pick **option 2** to analyze the V, J, and CDR3s. The header will be filled in and the sequences output to a new file that you name.

You can then input this file and pick **option 3** to output portions of the sequences to a new fasta file or **option 4** to output all sequences to an Excel-compatible text file. 

