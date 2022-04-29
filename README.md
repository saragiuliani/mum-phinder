
# MUM-PHINDER 
A framework to compute MUMs on large high repetitive datasets via matching startistics computation.

MUM-PHINDER computes the set of the Maximal Unique Matches of a query pattern against large highly-repetitive texts on a commodity computer. The index of the reference text is built with [PHONI](https://github.com/koeppl/phoni) [1] beforehand. An extended Matching Statistics is used to retrieve the MUMs of the query pattern.

We require the pattern and the text to be available in form of sequences stored in the `.fa` (FASTA) format.
To use our solution, you need to have recent `cmake`, `g++`, `zsh`, and `python 3` installed.

# Example

### Download MUM-PHINDER

```console
git clone --branch phoni https://github.com/saragiuliani/mum-phinder
```

### Compile

```console
mkdir build
cd build; cmake ..
make
```

### Run

##### Download the data 
From the `experiments\sars-cov2` folder:

```console
bash ./download.sh <dataFolder>
```

where `<dataFolder>` is the folder where to download the data. It requires Java 8.1 and unzip.


##### Build the index for a reference SARS-CoV2.2k.fa
From the `build` directory:

```console
python3 phoni-mum <dataFolder>/Ref/SARS-CoV2.2k.fa -f
```
This command will produce three files `SARS-CoV2.2k.fa.slp`, `SARS-CoV2.2k.fa.phoni`, and `SARS-CoV2.2k.fa.moni.log` in the `<dataFolder>/Ref` folder. The files contain the grammar, ..., and a log file with the information about the index (such as the length of the reference sequence, the alphabet size, the number of runs of the BWT, and more).

##### Compute the MUMs of the query MZ477765.fa against the reference SARS-CoV2.2k.fa with mum-phinder 

```console
./test/src/phoni_mum <dataFolder>/Ref/SARS-CoV2.2k.fa -p <dataFolder>/Pattern/MZ477765.fa 
```

This command will produce the output file `MZ477765.fa.mums` in the `<dataFolder>/Pattern` folder containing the MUMs of the query sequence `MZ477765.fa` against the reference `SARS-CoV2.2k.fa` in the following form: position in the reference,   position in the query,   length of the MUM.


# Authors

* Sara Giuliani
* Giuseppe Romana
* Massimiliano Rossi

# References

[1] [PHONI](https://github.com/koeppl/phoni) implementation: https://github.com/koeppl/phoni

