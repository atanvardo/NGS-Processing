# NGS-Processing
 Scripts for the processing of NGS COI sequences

## Warning!

Work in progress.

Instructions incomplete/outdated/with errors.

Don't use! (yet)


## Library design

In this experiment, we are sequencing a short COI fragment from a series of individuals.

Each PCR plate contains 92 individuals, plus 2 blanks (to check for lab contamination), plus 2 repeated samples (sequencing controls).

A library is formed by 4 PCR plates, i.e. 368 individuals. These 4 PCR plates can be seen as a *superplate* with 16 rows and 24 columns.

|  |  |
| --- | --- |
| Plate 1 | Plate 2 |
| Plate 3 | Plate 4 |

The primers for each row and column are marked with tags. The forward primer marks the row, and the reverse primer the column.

For example, the sample in the F2 well of the first (upper left) plate was amplified using the primers COIBF3_**6** and COIBR2_**2**. The sample in the E4 well of the fourth (lower right) plate will be amplified using the primers COIBF3_**13** and COIBR2_**16**.

The primers have these sequences:

| Name | Direction | Sequence |
| --- | --- | --- |
| **COIBF3_*x*** | Forward | `NxxxxxxCCHGAYATRGCHTTYCCHCG` |
| **COIBR2_*x*** | Reverse | `NxxxxxxTCDGGRTGNCCRAARAAYCA` |
