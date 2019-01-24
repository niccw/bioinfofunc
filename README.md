# bioinfofunc
Class and function for bioinformatics analysis

**BED**
```
./bioinfo_bed.py  -h
usage: bioinfo_bed.py [-h] [--header] [--cluster] [-d,--distance D]
                      [--bbannot] [--dictbed DBED] [--acol ACOL]
                      [--annot ANNOT] [--collapse]
                      bed

    Example:
    - Annotate bed by column in another bed (dbed)
    ./bioinfo_bed.py --bbannot main.bed --dictbed dict.bed --acol 4
    - Collapse bed (by chr,start,end), desciprtion/annotation are concatanated by ','
    ./bioinfo_bed.py --collapse sample.bed

    

positional arguments:
  bed              bed file

optional arguments:
  -h, --help       show this help message and exit
  --header         Default = False
  --cluster        Cluster featuresby distance. Work with -d/--distance
  -d,--distance D  Distance (in bp) to cluster features. Default = 2000
  --bbannot        Annotate bed based on another bed
  --dictbed DBED   BED provide information of annotation. Work with --acol.
                   This bed has no header.
  --acol ACOL      Column containing the annotation info (1-base)
  --annot ANNOT    New column name
  --collapse       Collapse bed, concat. annot by ','
  ```
  
  **GFF3**
  ```
./bioinfo_gff.py -h
usage: bioinfo_gff.py [-h] [--reformat] [--extract] [--rfmap RFMAP]
                      [--rfdata RFDATA] [--id ID] [--geneonly]
                      gff

    Example:
    - Extract gff by ID:
    ./bioinfo_gff.py --extract sample.gff3 --id id.list > result.gff3

    - Reformat attribute using data from tsv based on ID:
    ./bioinfo_gff.py --reformat sample.gff3 --rfmap sample.txt --rfdata sample.tsv  > reformated.gff3

    ##### sample.txt #####
    # To add new attribute, put '-' in the first column (old attribute)
    -   lc2
    -   condition
    -   de
    old_attr   new_attr

    ##### sample.tsv #####
    ID  lc2 condition   de  new_attr
    123   1   0.5       up    new
    ...
    ######################

    

positional arguments:
  gff              GFF3 file.

optional arguments:
  -h, --help       show this help message and exit
  --reformat       Reformat the gff attribute. Work with --rfmap and --rfdata.
  --extract        Reformat the gff attribute. Work with --id
  --rfmap RFMAP    2 columns tsv specific the rule of reformat. <attribite
                   name in gff/- (for new attribute)> <target column name in
                   rfdata>.
  --rfdata RFDATA  tsv metadata of gff. 'ID' column is used to map to gff 'ID'
                   attribute for row mapping. Make sure ID is unique!!
  --id ID          ID list (no header)
  --geneonly       Only process gene, ignore and print out mRNA and CDS and
                   other rows
 ```
 
**Repaetmasker OUT**
```
./bioinfo_repeatmasker_out.py -h
usage: bioinfo_repeatmasker_out.py [-h] [--gff3] [--saf] [--bed] out

positional arguments:
  out         RepeatMasker .out file

optional arguments:
  -h, --help  show this help message and exit
  --gff3      Convert out file to GFF3 format
  --saf       Convert out file to SAF format (for feature count)
  --bed       Convert out file to BED format , name =
              matching_repeat;repeat_class
```

**SAM/BAM**
```
./bioinfo_sam.py -h 
usage: bioinfo_sam.py [-h] [-ic] [-o,--out O] [--test] [--extract] [-f F]
                      [--flag FLAG]
                      sam

    Customised functions for SAM/BAM (sort + index)self.

    Examples:
    - Calculate insert coverage:
     ./bioinfo_sam.py -ic example.sam -o example
    - Extract reads:
     ./bioinfo_sam.py --extract example.sam -f qname.list --flag 163
    

positional arguments:
  sam          SAM or BAM file (sort + index)

optional arguments:
  -h, --help   show this help message and exit
  -ic          Calcualte the insert coverage and output .cov and .scaffoldmap
  -o,--out O   Prefix of output files
  --test       Run in test mode.
  --extract    extract read by qname. Work with -f and --flag (optional)
  -f F         list file
  --flag FLAG  extract read with specific flag
```
