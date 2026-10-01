# miRNA family seeds

`miR_Family_Info.human.tsv` is the human subset (Species ID 9606) of the author-supplied table: 2,606 records, 2,064 family labels. The original checksum and filtering record are in `provenance.json`. The database release and source reuse terms were not supplied and should be recorded by the authors.

`Seed+m8` is positions 2-8 of the mature miRNA, confirmed against every included Mature sequence. It is not the target RNA sequence. `AnnotateMiRNASeeds` reverse complements it before joining target k-mers. miR-7-5p: GGAAGAC -> GUCUUCC (RNA), GTCTTCC (DNA).

See [miRNA module](../../docs/MiRNA-seeds.md) and [TargetScan seed definitions](https://www.targetscan.org/docs/seed.html).
