# Instructions to run MKion

*Alex Brown*

*2023-05-16*

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

Unzip the `Virtual Machine/MKION.zip` folder and expand into the `Virtual Machine` directory.

Download and install the boxer MS-DOS emulator application:
http://boxerapp.com/



To run MKion, place the fasta file, and Valpha and Jalpha identifier files in the same folder. Execute the program and input the identifier files as instructed. MKion analyzes the data in several steps by menu. 



Run using boxer. Click the "Open a DOS prompt option".
Drag and drop the `mkionoc.exe` file into the dos prompt to get started.
Navigate to Window → Inspector, and adjust the CPU speed to max settings.



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
