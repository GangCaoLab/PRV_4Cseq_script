## 1. filter reads with primer

```bash
for x in {1..33}
do
    label=primer_set_${x}
    echo $label
    sed 's/primer_set_1/'$label'/g' < primer.input.information.txt >  primer.input.information.${x}.txt
    java -cp LGL.jar LGL.chiapet.PrimerFiltering_FastQ_PET primer.input.information.${x}.txt  primer.output  primer_set_${x}
done
```

### primer.input.information.txt

fastq_file_1	chen_S1_L001_R1_001.fastq
fastq_file_2	chen_S1_L001_R2_001.fastq
primer_file	primer_set_1.txt
maximum_mismatches	2
output_data_with_ambiguous_linker_info	0

## 2. cut primer

```bash
$ cd primer.output
$ perl cutreads-fastq.pl
```

## 3. mapping

Mapping using bwa.

## Contact

Please contact XD(695711550@qq.com) or XK(xke1995@126.com).
