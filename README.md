[![Release](https://img.shields.io/github/release/saragiuliani/mum-phinder.svg)](https://github.com/saragiuliani/mum-phinder/releases)
[![Downloads](https://img.shields.io/github/downloads/saragiuliani/mum-phinder/total?logo=github)](https://github.com/saragiuliani/mum-phinder/archive/master.zip)
[![Docker Pulls](https://badgen.net/docker/pulls/maxrossi91/mum-phinder?icon=docker&label=pulls)](https://hub.docker.com/r/maxrossi91/mum-phinder/)
[![Docker Image Size](https://badgen.net/docker/size/maxrossi91/mum-phinder?icon=docker&label=image%20size)](https://hub.docker.com/r/maxrossi91/mum-phinder/)

# MUM-PHINDER 
A framework to compute MUMs on large high repetitive datasets via matching startistics computation.

MUM-PHINDER computes the set of the Maximal Unique Matches of a query pattern against large highly-repetitive texts on a commodity computer. The index of the reference text is built with [PHONI](https://github.com/koeppl/phoni) [1] beforehand. An extended Matching Statistics is used to retrieve the MUMs of the query pattern.

We require the pattern and the text to be available in form of sequences stored in the `.fa` (FASTA) format.
To use our solution, you need to have recent `cmake`, `g++`, `zsh`, and `python 3` installed.

## How to get MUM-PHINDER

### Docker

MUM-PHINDER is available on `docker`:

```console
docker pull maxrossi91/mum-phinder:latest
docker run maxrossi91/mum-phinder:latest mum-phinder -h
```
if using `singularity`:
```console
singularity pull mum-phinder_sif docker://maxrossi91/mum-phinder:latest
./mum-phinder_sif mum-phinder --help
```

### Install Packages

We provide MUM-PHINDER on a `.deb` package:
```console
wget https://github.com/saragiuliani/mum-phinder/releases/download/v0.0.1/mum-phinder-v0.0.1_amd64.deb
sudo dpkg -i mum-phinder-v0.0.1_amd64.deb
mum-phinder -h
```
We provide MUM-PHINDER on a linux `.sh` installer:
```console
wget https://github.com/saragiuliani/mum-phinder/releases/download/v0.0.1/mum-phinder-v0.0.1-Linux.sh
chmod +x mum-phinder-v0.0.1-Linux.sh
./mum-phinder-v0.0.1-Linux.sh
mum-phinder -h
```
We provide MUM-PHINDER on a pre-compiled `.tar.gz`:
```console
wget https://github.com/saragiuliani/mum-phinder/releases/download/v0.0.1/mum-phinder-v0.0.1-Linux.tar.gz
tar -xzvf mum-phinder-v0.0.1-Linux.tar.gz
mum-phinder-v0.0.1-Linux/bin/mum-phinder -h
```

### Compile and install

```console
git clone https://github.com/saragiuliani/mum-phinder
cd mum-phinder
mkdir build
cd build; cmake -DCMAKE_INSTALL_PREFIX=<path/to/install/prefix> ..
make
make install
```

Replace `<path/to/install/prefix>` with your preferred install path. If not specified the install path is `/usr/bin` by default.



## Usage
### Construction of the index:
```
usage: mum-phinder build    [-h] -r REFERENCE [-w WSIZE] [-p MOD] [-t THREADS] [-k] [-v] [-f]
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        reference file name (default: None)
  -o OUTPUT, --output OUTPUT
                        output directory path (default: same as reference)
  -w WSIZE, --wsize WSIZE
                        sliding window size (default: 10)
  -p MOD, --mod MOD     hash modulus (default: 100)
  -t THREADS, --threads THREADS
                        number of helper threads (default: 0)
  -k                    keep temporary files (default: False)
  -v                    verbose (default: False)
  -f                    read fasta (default: False)

```


### Computing the MUMs with MUM-PHINDER:
```
usage: mum-phinder build [-h] -i INDEX -p PATTERN 
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        reference index base name (default: None)
  -p PATTERN, --pattern PATTERN
                        the input query (default: None)
```


# Example

### Download MUM-PHINDER

```console
git clone https://github.com/saragiuliani/mum-phinder
```

### Compile and Install

```console
cd mum-phinder
mkdir build
cd build; cmake -DCMAKE_INSTALL_PREFIX=<path/to/install/prefix> ..
make
make install
```

Replace `<path/to/install/prefix>` with your preferred install path. If not specified the install path is `/usr/`bin` by default.

### Run

##### Download the data 

**Important:** `Java 8.1` and `unzip` are required to download the data. They can be installed through `conda` with `conda install -c conda-forge openjdk unzip`.

To download the data, run:

```console
cd experiments/sars-cov2
bash ./download.sh <dataFolder>
```

where `<dataFolder>` is the folder where you want to download the data. 


##### Build the index for a reference SARS-CoV2.2k.fa
From the `build` directory:

```console
python3 mum-phinder build -r <dataFolder>/SARS-CoV2.2k.fa -f
```
This command will produce three files `SARS-CoV2.2k.fa.slp`, `SARS-CoV2.2k.fa.phoni`, and `SARS-CoV2.2k.fa.moni.log` in the `<dataFolder>/Ref` folder. The files contain the grammar, ..., and a log file with the information about the index (such as the length of the reference sequence, the alphabet size, the number of runs of the BWT, and more).

##### Compute the MUMs of the query MZ477765.fa against the reference SARS-CoV2.2k.fa with mum-phinder 

```console
python3 mum-phinder mums -i <dataFolder>/SARS-CoV2.2k.fa -p <dataFolder>/MZ477765.fa 
```

This command will produce the output file `MZ477765.fa.mums` in the `<dataFolder>/Pattern` folder containing the MUMs of the query sequence `MZ477765.fa` against the reference `SARS-CoV2.2k.fa` in the following form: position in the reference,   position in the query,   length of the MUM.


# Authors

* [Sara Giuliani](https://github.com/saragiuliani)
* [Giuseppe Romana](https://github.com/GiuseppeRomana)
* [Massimiliano Rossi](https://github.com/maxrossi91)

# References

[1] [PHONI](https://github.com/koeppl/phoni) implementation: https://github.com/koeppl/phoni

