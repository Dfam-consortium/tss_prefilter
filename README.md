
# Tensor Slide Sketch - Prefilter For Alignment (tss_prefilter)

## Description


## Dependencies

FAIS 1.8.0  : https://github.com/facebookresearch/faiss/archive/refs/tags/v1.8.0.tar.gz

```
% tar zxvf v1.8.0.tar.gz
% cd faiss-1.8.0
% /usr/local/cmake/bin/cmake -DFAISS_ENABLE_GPU=OFF -DFAISS_ENABLE_C_API=ON -DBUILD_SHARED_LIBS=ON -B build .
% make -C build -j faiss
% make -C build install
```
### Does not appear to install C library files by default

```
% cd build/c_api
% cp lib* /usr/local/lib
```

## Build

```
% cargo build --release
```

## Example

### Build a small index on a subset of the Dfam TE consensus sequences

```
% head -n 1000 data/Dfam36-human.fa > sequences.fasta
% /bin/time target/release/tss_prefilter build
```

### Search for the nearest neighbors to query sequence(s)

First create a file called 'query.fasta' containing one or more query sequences in FASTA format.
Here we use mir.fa from the data directory.
Then run the search command.

```
% cp data/mir.fa query.fasta
% target/release/tss_prefilter search
```

### Build a larger index using a chunk of Dfam36 representing all of the human TEs (plus some more)

```
% rm faiss.index
% cp data/Dfam36-human.fa > sequences.fasta
% /bin/time target/release/tss_prefilter --build
```

### Search for the nearest neighbors to SVA
This TE family is known to contain a bit of Alu, HERVK, and LTR5 elements

```
% cp data/sva_a.fa query.fasta
% target/release/tss_prefilter --search  
```

### Align candidates and lookup TE names for identifiers
This is just an example script and will not work out-of-the-box.  It is meant to
demonstrate how to use the output of the search command to align the candidates.
It also performs a lookup of the TE names for the identifiers.  NOTE: It can be
made to work by changing the path to the 'align.pl' script if you have
installed a RepeatMasker/RepeatModeler distribution on your system.

```
/bin/time target/release/tss_prefilter --search | scripts/alignQueryToSet.pl
```

Robert Hubley 5/2024

