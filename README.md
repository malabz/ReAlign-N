# ReAlign-N: an integrated realignment approach for multiple nucleic acid sequence alignment, combining global and local realignments

ReAlign-N is a tool written in C++11 for realigning the multiple nucleic acid sequence alignment. It runs on Linux.

## 🔨Installation and Usage

### 1.1 Linux/WSL(Windows Subsystem for Linux ) - from Anaconda
1.Install WSL for Windows. Instructional video [1](https://www.youtube.com/watch?v=X-DHaQLrBi8&t=5s) or [2](http://lab.malab.cn/%7Etfr/1.mp4) (Copyright belongs to the original work).

2.Download and install Anaconda. Download Anaconda for different systems [here](https://www.anaconda.com/products/distribution#Downloads). Instructional video of anaconda installation [1](https://www.youtube.com/watch?v=AshsPB3KT-E) or [2](http://lab.malab.cn/%7Etfr/Install_anaconda_in_Linux.mp4) (Copyright belongs to the original work).

3.Install ReAlign-N.
```bash
#1 Create and activate a conda environment for ReAlign-N
conda create -n realign_n_env
conda activate realign_n_env

#2 Add channels to conda
conda config --add channels malab

#3 Install ReAlign-N
conda install -c malab realign_n

#4 Test ReAlign-N
realign_n -h
```

### 1.2 Linux/WSL(Windows Subsystem for Linux ) - from the source code

1.Download and Compile the source code. (Make sure your version of gcc >= 9.4.0)
```shell
#1 Download
git clone https://github.com/malabz/ReAlign-N.git

#2 Open the folder
cd ReAlign-N

#3 Compile
make

#4 Test ReAlign-N
./realign_n -h
```

### 2 Usage
```
Usage: /.realign_n [-r] path [-a] path [-o] path [-m] distance [-e] distance [-p] pattern [-h]

  Necessary arguments:
    -r  Specify the path of raw data, a file in FASTA format.
    -a  Specify the path of initial alignment, a file in FASTA format.

  Optional arguments:
    -o  Specify the output for ReAlign-N, a file in FASTA format.
    -m  Specify the minium split distance of match (default based on the similarity).
    -e  Specify the minium split distance of entropy (default based on the similarity).
    -p  Specify the pattern of ReAlign-N (default pattern: 1).
          1 for local realignment followed by global realignment.
          2 for global realignment followed by local realignment.
    -h  Print the help message.
```
### 3 ReAlign-N Parameters Description
1. Usage: 
```
/.realign_n [-r] path [-a] path [-o] path [-m] distance [-e] distance [-p] pattern [-h]
```
2. Required Arguments:
*   -r

    Specify the path to the raw data file in FASTA format.
    
    Example: -r /path/to/raw_data.fasta

*   -a

    Specify the path to the initial alignment file in FASTA format.

    Example: -a /path/to/initial_alignment.fasta

3. Optional Arguments:

*   -o

    Specify the path for the ReAlign-N output file in FASTA format.

    If not provided, the output will be saved to the default path（/CURRENT_PATH/realign_n_result.fas）.

    Example: -o /path/to/output.fasta

*   -m

    Specify the minimum split distance for matching. 
    
    If not provided, the default value is automatically determined based on the similarity of the input data.

    Example: -m 10 sets the minimum split distance for matches to 10.

*   -e

    Specify the minimum split distance for entropy. 
    
    If not provided, the default value is automatically determined based on the similarity of the input data.

    Example: -e 5 sets the minimum split distance for entropy to 5.

*   -p

    Specify the alignment pattern for ReAlign-N (default pattern: 1).

    1: Local realignment followed by global realignment.

    2: Global realignment followed by local realignment.

    Example: -p 2 sets the pattern to global realignment followed by local realignment.

*   -h

    Print the help message displaying all available parameters and their descriptions.


4. Minimum split distances correspond to different average similarity levels and sequence counts.

*   The default value for match

    Avg similarity|Sequences Counts|Distance
    :---:|:---:|:---:
    [0.99,1]|<200|5
    [0.96,0.99)|<200|10
    [0.93,0.96)|<200|15
    [0.90,0.93)|<200|20
    [0.80,0.90)|<200|30
    [0,0.80)|<200|50
    [0.90,1]|>=200|10
    [0,0.90)|>=200|50

*   The default value for entropy

    Avg similarity|Sequences Counts|Distance
    :---:|:---:|:---:
    [0.97,1]|<200|5
    [0.90,0.97)|<200|50
    [0,0.90)|<200|5
    [0.90,1]|>=200|10
    [0,0.90)|>=200|50

*   Note
    
    If the user does not specify -e **AND** -m in the command line, their default values are set based on the similarity observed during the initial alignment. 
    
    If the user specifies only -e **OR** -m, the program will use the provided parameter for calculation, while the unspecified parameter will be determined based on the default value (i.e., the similarity of the initial alignment). 
    
    If both -e **AND** -m are specified, the program will calculate entirely based on the user’s input.


## 🔬Test dataset and the use case
### 1. Information about the test dataset

Dataset|Sequences Num|Repeats Num|Avg Length|Similarity
:---:|:---:|:---:|:---:|:---:
16S rRNA|1000|10|about 1440bp|The average similarity is about 72%
mt genome|200|10|about 16570bp|The average similarity is about 99%
23S rRNA|500|10|about 3120bp|The average similarity is about 92%
16S-like rRNA|100|9|about 1550bp|14 sets of data with different similarities (99%, 98%, 97%, 96%, 95%, 94%, 93%, 92%, 91%, 90%, 85%, 80%, 75%, 70%)
mt-like genome|100|9|about 16000bp|14 sets of data with different similarities (99%, 98%, 97%, 96%, 95%, 94%, 93%, 92%, 91%, 90%, 85%, 80%, 75%, 70%)
SARS-CoV-2-like genome|100|9|about 29000bp|14 sets of data with different similarities (99%, 98%, 97%, 96%, 95%, 94%, 93%, 92%, 91%, 90%, 85%, 80%, 75%, 70%)
CIPRES-128|255|9|about 1550bp|The average similarity is about 80%
CIPRES-256|511|9|about 1550bp|The average similarity is about 80%
CIPRES-512|1023|9|about 1550bp|The average similarity is about 80%
CIPRES-1024|2047|9|about 1550bp|The average similarity is about 80%

### 2. The use case
```shell
# Download data
wget http://lab.malab.cn/soft/ReAlign-N/data/16s_like.tar.gz

# Unzip data
tar -zxvf 16s_like.tar.gz

# Get the folder path
cd 16s_like

# Run ReAlign-N
./realign_n -r raw_data/16s_similarity_70_1.fas -a msa_results/16s_similarity_70_1_clustalo.fas -o 16s_similarity_70_1_clustalo_realign_n.fas -p 1
```
## 📍Reminder
1. Currently ReAlign-N is **ONLY** available for DNA/RNA. 
2. Please ensure that the sequence contains only AGCT (for DNA) or AGCU (for RNA).
3. Please ensure that the sequence ID entered into ReAlign is unique.
4. MAFFT installation is required for the utilization of ReAlign-N. 

## 🖥️Environment
System|GCC version
:---:|:---:
Linux|GCC 9.4.0
WSL|GCC 9.4.0

## 🔖Citation
If you use this software, please cite:

Yixiao Zhai, Tong Zhou, Yanming Wei, Quan Zou, Yansu Wang, ReAlign-N: an integrated realignment approach for multiple nucleic acid sequence alignment, combining global and local realignments, NAR Genomics and Bioinformatics, Volume 6, Issue 4, December 2024, lqae170, https://doi.org/10.1093/nargab/lqae170

## 👋Contacts
The software tools are developed and maintained by 🧑‍🏫[ZOU's lab](http://lab.malab.cn/~zq/en/index.html).

If you find any bug, welcome to contact us on the [issues page](https://github.com/malabz/ReAlign-N/issues) or email us at 👉[📩](zhai1xiao@gmail.com).

More tools and infomation can visit our [github](https://github.com/malabz).
