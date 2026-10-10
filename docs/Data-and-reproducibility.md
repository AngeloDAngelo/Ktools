# Data and reproducibility

[Documentation](Documentation.md) · [Tutorial](Tutorial.md)

## Included examples

- **Synthetic demo:** 24 reference and 48 control sequences, generated with seed 42. The demo creates fictional RBP scores using seed 123. [Inputs and demo outputs](../examples/README.md).
- **CDR1as:** a public circBase-derived sequence mapped to GRCh38, with seed annotations and example KRS scores. The comparison FASTA is supplied by the user. [Sequence source and launch command](CDR1as-miRNA-example.md).

## Reference resources

- [PEKA matrices](../data/peka/README.md): choose one of the two included score matrices.
- [miRNA family tables](../data/mirna/README.md): complete table and human subset; the complete table is also installed with the Python package.
- [UniProt BED inputs](../data/uniprot/README.md): required file formats for K-CDS. Provide BED tracks and a matching GTF before running this module.

## Run records

Record the code commit, input checksums, annotation release, group membership, parameters and random seeds. Preserve full result tables alongside selected signatures. Save dependency versions with:

```bash
python -m pip freeze > environment.txt
```
