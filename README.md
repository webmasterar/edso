## EDSO

EDSO (Elastic Degenerate Sequence Outputter) is a little utility program that
combines a reference fasta file and its VCF file into an EDS format file.
Usage example:

`./edso reference.fasta variants.vcf.gz outfile.eds`

If you wish to generate random EDS files, please visit https://github.com/webmasterar/EDSRand

### Installation

First run `./pre-install.sh`, then compile EDSO by running `make`.

License: GNU GPLv3 License; Copyright (C) 2017 Ahmad Retha
