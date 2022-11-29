# SysISTD
A demultiplexer for SysFate Illumina spatial transcriptomic.

## Description


### How to run

```
require arguments:
  --f1 file1
                           The fastq file 1 you want analyse. (It can be compress with extension .gz)
  --f2 file2
                           The fastq file 2 you want analyse. (It can be compress with extension .gz)
  --coor coordinate
                           TSV file who contain coordinates code with nucleique barcode
  --pos position
                           TSV file who contain coordinates code with their position in the grid
  -g genome, --genome genome
                           Aligner genome index path
  --gtf gtf
                           GTF file to use with feature count
  ```
```
optional arguments:
  -h, --help               Show this help message and exit
  --version                Show program's version number and exit
  -n name, --name name     You can give a name for your analyse. By default the software use the begin date.
  -o outdir, --out outdir  Folder where send files (default: .)
  --tmp tmp                Folder who contain temporary files (default: /tmp)
  --guibson guibson        Sequences guibson (default: ACATTGAAGAACCTGTAGATAACTCGCTGT)
  -G guibsonE, --guibsonE guibsonE
                           Guibson error accepted (default: 1)
  -B barcodeE, --barcodeE barcodeE
                           Barcode error accepted (default: 1)
  --polyTE polyTE          PolyT error accepted (default: 1)
  -c core, --core core     Number of core the software can use (default: 1)
  --polyTSizeMin polyTSizeMin
                           Size min for a polyT (default: 12)
  --nopolytsearch          If you don't want validate the research with polyT (default: False)
  --trimSequence trimSequence
                           After the UMI, what do you want trim ? A fix number of base or a specific sequence. By default, trim the polyT as you config before (default: )
  --noUMIsearch            If you don't want research UMI (default: False)
  -a aligner, --aligner aligner
                           This is the aligner command. You can write what you want if the __file1__, __file2__, __genome__ is provide and the output is configured to send the results in __outputSam__. (default: bowtie2 --reorder --very-sensitive -p __core__ -x __genome__ -1 __file1__ -2 __file2__ -S __outputSam__)
  -q quantifier, --quantifier quantifier
                           This is the read quantifier command. You can write what you want if the __GTF__, __feature__ is provide and the output is configured to send the results in __out__. (default: featureCounts -T __core__ -F GTF -t __feature__ -s 0 -a __GTF__ -o __out__)
  --feature feature        Feature use during the read quantifier in your file separate by a comma (default: exon,transcript)
  --decompress decompress
                           Command use to decompress fastq files (default: pigz -dc)
  --dev                    Keep all temporar datas (default: False)
  ```

## Authors

|         |                                                                                               |
| ------- | --------------------------------------------------------------------------------------------- |
| Author  | STUDER Francois ([Github](https://github.com/studyfranco))                                    |
| Author  | MOEHLIN Julien ([Github](https://github.com/JulienMoehlin), [Gitlab](https://gitlab.com/julienmoehlin)) |
| Author  | MENDOZA PARRA Marco ([Github](https://github.com/SysFate))                                    |
| Team    | [SysFate](https://www.sysfate.org/)                                                           |
| Email   | <mmendoza@genoscope.cns.fr>                                                                   |

## How the software works
For each sequence of all two files:
  - SysISTD searches for gibson sequences with a regex query. This library allows to define the number of errors to be accepted within the regexpression.
  - Search the two barcode before and after the guibson and verify them good orientation
  - After the validation of the sequence, it align all positions
  - Launch the quantification of each gene with feature count and create a coordinate matrix with all genes
