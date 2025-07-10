## EDSO

EDSO (Elastic Degenerate Sequence Outputter) is a little utility program that
combines a reference fasta file and its VCF file into an EDS format file.
Usage example:

`./edso reference.fasta variants.vcf.gz outfile.eds`

### Installation

First run `./pre-install.sh`, then compile _EDSO_ by running `make`.

### Troubleshooting

The first argument `reference.fasta` must be a plain text file. Gzipped files
are not supported.

The second argument `variants.vcf.gz` should have a counterpart `variants.vcf.gz.tbi`
file. This is usually available from the same source you got the vcf file from,
but if it isn't provided, you can generate it yourself:

```
sudo apt-get install tabix
tabix -f -p vcf variants.vcf.gz
```

If you get an error while running `pre-install.sh` it might be that you haven't
made it executable. Try this:

```
chmod +x pre-install.sh
./pre-install.sh
```

If you get a compilation error like this: `src/cdflib.cpp:10043:15: error: ordered comparison of pointer with integer zero (‘int*’ and ‘int’)`

This occurs with modern compilers that are more sensitive to coding bugs. Just
fix the error. Open up `vcflib/src/cdflib.cpp` and go to line 10043. Change it
to `if ( *n_data < 0 )` and save it. Now open up `pre-install.sh` and comment
out the line like so: `# tar -xvf vcflib.tar.gz`, save it. Now run: `./pre-install.sh`.
Then, `make` the _EDSO_ binary.

License: GNU GPLv3 License; Copyright (C) 2017 Ahmad Retha.
Last updated 2025 to fix Issue #2 and add troubleshooting advice.
