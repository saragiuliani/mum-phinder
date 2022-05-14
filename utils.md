# Utils command

## MUMmer

### Install using `conda`

```console
conda create --name mummer -c bioconda -y mummer
conda activate mummer
```

### Run `MUMmer` on the example

```console
mummer -mumreference -l 1 -n ../data/yeast.fasta.plain ../data/query.fasta > mummer.out
```

Observations:
* The `-mumreference` flag makes `MUMmer` report maximal matches that are unique only in the reference. Replace it with `-mum` to unique matches.
* The `-l` flag defines the minimum match length.
* The`-n` flag forces to match only `A`, `C`, `G`, and `T` characters.
* The input fasta file should contain only one sequence, otherwise `MUMmer` finds `MUM`s on each individual subsequence, which is different from what we can do. 

## PHONI

From the `build folder`.

### Building the index

```console
python3 phoni-mum ../data/yeast.fasta -f
```

### Run `PHONI` on the example
```console
./test/src/phoni_mum ../data/yeast.fasta -p ../data/query.fasta 2>phoni_mum.out
```

## Comparison

```console
diff mummer.out phoni_mum.out
```

# Build the docker image

```console
docker build --platform linux/amd64 --no-cache -t maxrossi91/mum-phinder . 
```

docker run --platform linux/amd64 -v `pwd`/experiments/sars-cov2/data:/data  -it maxrossi91/mum-phinder bash