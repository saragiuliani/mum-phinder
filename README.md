# MUM-PHINDER 
A framework to compute MUMs on large high repetitive datasets via matching startistics computation.

MUM-PHINDER computes the set of the Maximal Unique Matches of a query pattern against large highly-repetitive texts on a commodity computer. The index of the reference text is built with [PHONI](https://github.com/koeppl/phoni) beforehand. An extended Matching Statistics is used to retrieve the MUMs of the query pattern.

We require the pattern and the text to be available in form of sequences stored in the `.fa` (FASTA) format.
To use our solution, you need to have recent `cmake`, `g++`, `zsh`, and `python 3` installed.

#Example

### Download

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

##### Build the index for a reference SARS-CoV2.fa
From the `build` directory

```console
python3 phoni-mum ../experiments/Example/Ref/SARS-CoV2.fa -f
```
This command will produce three files `SARS-CoV2.fa.slp`, `SARS-CoV2.fa.phoni`, and `SARS-CoV2.fa.moni.log` in the `SARS` folder. The files contains the grammar, ..., and a log with information about the index (such as the length of the reference sequence, the alphabet size, the number of runs of the BWT, and more).

##### Compute the MUMs of the query MZ477765.1.fa agains the reference SARS-CoV2.fa with mum-phinder 

```console
./test/src/phoni_mum ../experiments/Example/Ref/SARS-CoV2.fa -p ../experiments/Example/Pattern/MZ477765.1.fa
```

This command will produce the output file `MZ477765.1.fa.mums` in the `Pattern` folder containing the MUMs of the query sequence `MZ477765.1.fa` against the reference `SARS-CoV2.fa` in the following form: position in the reference,   position in the query,   length of the MUM.


# Authors

* Sara Giuliani
* Giuseppe Romana
* Massimiliano Rossi

# References

[1] [PHONI](https://github.com/koeppl/phoni) implementation: https://github.com/koeppl/phoni


