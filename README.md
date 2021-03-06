# Software Engineering Course
# Week #4: C++ programming
 
Instructor: Daehwan Kim (Daehwan.Kim@UTSouthwestern.edu)

Teaching assistants: Chanhee Park (Chanhee.Park@UTSouthwestern.edu) and Micah Thornton (Micah.Thornton@UTSouthwestern.edu).

The class meets at 9am from May 23rd (Day 1) through May 27th (Day 5) in ND11.218.

Office hours will be 2 - 4pm on each of Days 1-5 in E4.350.


This course provides hands-on practice with C++ programming for implementing DNA sequence alignment. You will design and implement memory-efficient data structures and rapid algorithms for suffix arrays, Burrows-Wheeler transforms (BWT), and Ferragina-Manzini (FM) indexes. As in the previous weeks, it is important for you to write code in a modular and easy-to-maintain fashion. You are encouraged to provide documentation in the form of comments along with your code so that other future developers can understand and extend it further. You will be expected to write code so that your code will live forever.


Prerequisites: 
- Familiarity with C++ features such as array, pointer, and class
- Understanding of concepts in Object-oriented programming
- Knowledge of basic algorithm techniques such as binary search and sorting

We will provide relevant articles and slides that you can read on your own in order to carry out programming assignments as described in the daily course schedule below:

## Day 1
From 9 to 10am, we will provide a brief introduction to sequence alignment, suffix arrays, and a C++ coding environment. For the rest of the day, you will implement suffix arrays (SA) for two genomic sequences: a small 100-bp sequence and a large 51 million-bp long sequence (human chromosome 22).
 
## Day 2
From 9 to 10am, we will provide a brief introduction to sequence alignment and have a Q&A session. For the rest of the day, you will implement an alignment algorithm using the suffix array. We will provide input sequencing reads and true output alignment, which you can use to test your alignment algorithm in terms of alignment accuracy and speed. You can also implement test cases and use assertions so that you can effectively identify and fix bugs in your code.
Optional: you can implement a suffix array and an alignment algorithm for the human reference genome (about 3 billion bases) that can be used on your computer.
 
## Day 3
From 9 to 10am, we will introduce Burrows-Wheeler transforms (BWT) and have a Q&A session. For the rest of the day, you will implement BWT for two genomic sequences: a small 100-bp sequence and a large 51 million-bp long sequence (human chromosome 22).
Optional: you can implement a BWT for the human reference genome.
 
## Day 4
From 9 to 10am, we will cover sequence alignment and have a Q&A session. For the rest of the day, you will implement an alignment algorithm using BWTs for two genomic sequences: a small 100-bp sequence and a large 51 million-bp long sequence (human chromosome 22).
Optional: you can implement an alignment algorithm using a BWT for the human reference genome.
 
## Day 5
From 9 to 10am, we will provide a brief introduction to Ferragina-Manzini (FM) indexing and have a Q&A session. For the rest of the day, you will implement FM indexes and an alignment algorithm using FM indexes for two genomic sequences: a small 100-bp sequence and a large 51 million-bp long sequence (human chromosome 22).
Optional: you can make this FM index and alignment work for the human reference genome.
 
## Days 1 through 5
The performance standards for you to achieve in producing your algorithms include lean memory usage and high speed. Please document your memory usage and runtime information. The output formats for your programs are described [below](#output-formats). You will be asked to review one of your classmates??? code implementations, with the code to review assigned by the course instructors. Your review summary should not be more than 200 words.

All assignments are to be done individually, but discussion with classmates regarding general solution strategies is allowed.

## Grade
This C++ course will constitute 20% of the final grade or 20 points as follows:
 
|  | Day 1: SA implementation | Day 2: SA-based alignment implementation | Day 3: BWT implementation | Day 4: BWT-based alignment implementation | Day 5: FM index implementation |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| Correctness | 2 | 2 | 2 | 2 | 3 |
| Memory efficiency and runtime | 1 | 1 | 1 | 1 | 1 |
| Code Review | 0 | 1 | 1 | 1 | 1 |

 
The deadline for submitting a programming assignment for Day 1 is Day 2 at 9am. As such, the deadlines for submitting programming assignments for Day 2, Day 3, Day 4, and Day 5 are Day 3 at 9am, Day 4 at 9am, Day 5 at 9am, and Day 5 at 10pm.

The deadlines for submitting code review assignments for Day 2, Day 3, Day 4, and Day 5 are Day 3 at 9am, Day 4 at 9am, Day 5 at 9am, and Day 5 at 10pm.

## Output Formats
The output file formats are described as follows.

For each day we will assess the correctness of the output of the code.  Please produce codes which write 
the results of the computation to a file on disk according to each of the following: 

| | Input Format Description | Input Format Example | Output Format Description | Output Format Example | Command Example | 
| ----------- | ---------- | ----------- | ----------- | -------- | -------- | 
| Day 1: SA Impl. | Input .fa file, single sequence multi-line, > for header | [sample\_ref.fa](/Example_Files/sample_ref.fa) | A list of integers in txt format (label .sa) | [sample\_ref.sa](/Example_Files/sample_ref.sa) | <your_program_name> sample\_ref.fa sample\_ref.sa |
| Day 2: SA Search | Input .fa reference file (or input .sa from Day 1), and Multi-line query fasta file (labeled .qfa) | [abracadabra.fa](/Example_Files/abracadabra.fa) (or [abracadabra.sa](/Example_Files/abracadabra.sa)) and [abracadabra.qfa](/Example_Files/abracadabra.qfa)| A tab-delimited file where each row contains query name and alignment location, multiple rows for the same query matching multiple locations. No entries for queries that don't align | [abracadabra.aln](/Example_Files/abracadabra.aln) | <your_program_name> sample.sa sample.qfa > sample.aln |
| Day 3: BWT Impl. | Input .fa reference file, input suffix array .sa | [abracadabra.fa](/Example_Files/abracadabra.fa) and [abracadabra.sa](/Example_Files/abracadabra.sa) | text file containing burrows-wheeler transform text (with or without extra character) named .bwt | [abracadabra.bwt](/Example_Files/abracadabra.bwt) |
| Day 4: BWT Search | Input .bwt, input .sa, and queries file (.qfa) | [abracadabra.bwt](/Example_Files/abracadabra.bwt) and [abracadabra.sa](/Example_Files/abracadabra.sa) and [abracadabra.qfa](/Example_Files/abracadabra.qfa) | A tab-delimited file where each row contains query name and alignment location, multiple rows for the same query matching multiple locations. |[abracadabra.aln](/Example_Files/abracadabra.aln)  | 
| Day 5 part 1: FM Index Impl. | Input .bwt, input .sa | [abracadabra.bwt](/Example_Files/abracadabra.bwt) and [abracadabra.sa](/Example_Files/abracadabra.sa) and commandline options for (1) number of rows to skip in checkpointed tally matrix and (2) number of values to skip when subsampling suffix array | An FM Index file, with FM-Index components seperated by > and appropriate headers as shown in the example output. | [abracadabra.fm](/Example_Files/abracadabra.fm)|
| Day 5 part 2: FM Index Search | input .fm, input .qfa | [abracadabra.fm](/Example_Files/abracadabra.fm) and [abracadabra.qfa](Example_Files/abracadabra.qfa) | a tab-delimited file where each row contains query name and alignment location, multiple rows for the same query matching multiple locations | [abracadabra.aln](/Example_Files/abracadabra.aln) | 


## Datasets
Each dataset contains a reference sequence, reads, and a position of read in the reference. Student may use the Human Reference Genome to do optional assignment. When you test your program with chr22_reads.fa and genome_reads.fa, you can use a subset of reads, for example, first 1000 reads. But when you test a peformance of your program, we recommend to use whole reads.  

- Small  
   - Reference: 100bp sequence [small.fa]
   - Reads: 10 10bp reads [small_reads.fa]
   - Read positions: [small.aln]
- Human Chromosome 22
   - Reference: [chr22.fa]
   - Reads: 1Million 100bp reads [chr22_reads.fa]
   - Read positions: [chr22.aln] 
- Human Reference Genome
   - Reference: [genome.fa]
   - Reads: 1Million 100bp reads [genome_reads.fa]
   - Read positions: [genome.aln]


[small.fa]:        /Example_Files/small.fa
[small_reads.fa]:  /Example_Files/small_reads.fa
[small.aln]:       /Example_Files/small.aln
[chr22.fa]:        https://cloud.biohpc.swmed.edu/index.php/s/HBnzJPfoK46HWrX/download/chr22.fa
[chr22_reads.fa]:  https://cloud.biohpc.swmed.edu/index.php/s/RJmjzzwzi9My2Am/download/chr22_reads.fa
[chr22.aln]:       https://cloud.biohpc.swmed.edu/index.php/s/k7YjgzJ8LLDyqWa/download/chr22.aln
[genome.fa]:       https://cloud.biohpc.swmed.edu/index.php/s/SsnLjF92scW6FXH/download/genome.fa
[genome_reads.fa]: https://cloud.biohpc.swmed.edu/index.php/s/sZWTFtLyD7J7TGx/download/genome_reads.fa
[genome.aln]:      https://cloud.biohpc.swmed.edu/index.php/s/n4kgZJpcLjry7aj/download/genome.aln

## Setting up environments
We recommend using a laptop/desktop rather than the biohpc client.

### Linux

1. Install gcc/g++ compiler
   - Ubuntu
   ```sh
   sudo apt install build-essential
   ```
   - CentOS
   ```sh
   sudo yum group install "Development Tools"
   ```
   
1. Download CLion package. This package includes an evaluation license key for a free 30-day trial.  
   CLion - [Download](https://download.jetbrains.com/cpp/CLion-2022.1.1.tar.gz)  
1. Unpack the package
   ```sh
   tar zxvf CLion-2022.1.1.tar.gz
   ```
1. Run CLion.sh from the bin subdirectory
1. Open your project directory (`Day1`, `Day2`, ...)


----
### Windows
Visual Studio 2022 Community - [Download](https://c2rsetup.officeapps.live.com/c2r/downloadVS.aspx?sku=community&channel=Release&version=VS2022&source=VSLandingPage&includeRecommended=true&cid=2030)
1. Open a solution file in the project directory(`Day1.sln`, `Day2.sln`, ...)


----
### MacOS
Xcode - [Download](https://apps.apple.com/us/app/xcode/id497799835?mt=12)
1. Open a xcode project in the project directory(`Day1.xcodeproj`, `Day2.xcodeproj`, ...)
