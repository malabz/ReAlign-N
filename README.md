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

1. Download and Compile the source code. (Make sure your version of gcc >= 9.4.0)
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
Usage: /.realign_n [-r] path [-a] path [-o] path [-m] mode

  Necessary arguments:
    -r  Specify the path of raw data, a file in FASTA format.
    -a  Specify the path of initial alignment, a file in FASTA format.

  Optional arguments:
    -o  Specify the output for ReAlign-N, a file in FASTA format.
    -m  Specify the mode of ReAlign-N (default mode: 1).
        1 for local realignment followed by global realignment.
        2 for global realignment followed by local realignment.
    -h  Print the help message.
```

## 🔬Test dataset and the use case
### 1. Information about the test dataset

Dataset|Sequences Num|Repeats Num|Avg Length|Similarity
:---:|:---:|:---:|:---:|:---:
16s simu|100|9|about 1550bp|14 sets of data with different similarities (99%, 98%, 97%, 96%, 95%, 94%, 93%, 92%, 91%, 90%, 85%, 80%, 75%, 70%)
mt simu|100|9|about 16000bp|14 sets of data with different similarities (99%, 98%, 97%, 96%, 95%, 94%, 93%, 92%, 91%, 90%, 85%, 80%, 75%, 70%)
sars2 simu|100|9|about 29000bp|14 sets of data with different similarities (99%, 98%, 97%, 96%, 95%, 94%, 93%, 92%, 91%, 90%, 85%, 80%, 75%, 70%)
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
./realign_n -r raw_data/16s_similarity_70_1.fas -a msa_results/16s_similarity_70_1_clustalo.fas -o 16s_similarity_70_1_clustalo_realign_n.fas -m 1
```
## 📍Reminder
1. Currently ReAlign-N is **ONLY** available for DNA/RNA. 
2. Ensure that the sequence ID entered into ReAlign is unique.
3. MAFFT installation is required for the utilization of ReAlign-N. 

## 🖥️Environment
System|GCC version
:---:|:---:
Linux|GCC 9.4.0
WSL|GCC 9.4.0

## 🔖Citation


## 👋Contacts
The software tools are developed and maintained by 🧑‍🏫[ZOU's lab](http://lab.malab.cn/~zq/en/index.html).

If you find any bug, welcome to contact us on the [issues page](https://github.com/malabz/ReAlign-N/issues) or email us at 👉[📩](zhai1xiao@gmail.com).

More tools and infomation can visit our [github](https://github.com/malabz).
