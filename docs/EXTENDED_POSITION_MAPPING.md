# Extended Position Mapping for tRNAscan-SE

This document describes the extended position mapping functionality that provides comprehensive coordinate mapping for tRNA sequences, including sequence positions, alignment positions, and Sprinzl structural positions.

## Overview

The extended position mapping system allows you to:

1. **Map sequence coordinates** - Original nucleotide positions in the input sequence
2. **Map alignment coordinates** - Positions in the covariance model alignment
3. **Map Sprinzl coordinates** - Standardized structural positions in the tRNA
4. **Include genomic coordinates** - Chromosomal positions for genomic sequences

## Files Added

### Core Module
- `lib/tRNAscanSE/ExtendedSprinzlResults.pm` - Main module containing position mapping functions

### Scripts
- `scripts/extended_position_mapping.pl` - Wrapper script that calls tRNAscan-SE and generates mappings
- `scripts/generate_position_map.pl` - Direct integration script using tRNAscan-SE APIs

### Documentation
- `docs/EXTENDED_POSITION_MAPPING.md` - This file

## Functions Available

### `write_extended_position_map($global_vars, $file, $tRNA)`
Generates position mapping for a single tRNA with columns:
- `tRNA_ID` - Unique identifier
- `Seq_Pos` - Position in original sequence (1-based)
- `Nucleotide` - The nucleotide at this position
- `Align_Pos` - Position in covariance model alignment (0-based)
- `Sprinzl_Pos` - Sprinzl structural position
- `Isotype` - tRNA amino acid type
- `Anticodon` - Anticodon sequence
- `Score` - tRNAscan-SE confidence score

### `write_extended_position_db($global_vars, $file)`
Generates position mapping for all tRNAs in a dataset with additional genomic information:
- All columns from single tRNA mapping, plus:
- `Genomic_Start` - Start coordinate in input sequence
- `Genomic_End` - End coordinate in input sequence  
- `Strand` - Orientation (+ or -)
- `Sequence_Name` - Name of input sequence

### `write_genomic_coordinate_map($global_vars, $file)`
Generates mapping that includes actual genomic coordinates:
- `tRNA_ID` - Unique identifier
- `Genomic_Coord` - Actual genomic coordinate for each nucleotide
- `Seq_Pos` - Position in tRNA sequence
- `Nucleotide` - The nucleotide
- `Align_Pos` - Alignment position
- `Sprinzl_Pos` - Sprinzl position
- `Isotype` - tRNA type
- `Anticodon` - Anticodon sequence
- `Strand` - Orientation
- `Sequence_Name` - Sequence name

## Usage Examples

### Basic Usage
```bash
# Run with default eukaryotic model
perl scripts/extended_position_mapping.pl input_sequences.fa output_positions.tsv

# Use bacterial model
perl scripts/extended_position_mapping.pl -B bacterial_genome.fa positions.tsv

# Include genomic coordinate mapping
perl scripts/extended_position_mapping.pl -E --genomic eukaryotic_sequences.fa positions.tsv
```

### Using the Direct API Script
```bash
# Process sequences as mock tRNAs (for testing/development)
perl scripts/generate_position_map.pl input.fa output.tsv

# With genomic coordinates
perl scripts/generate_position_map.pl --genomic input.fa output.tsv
```

### Integration into Existing Code
```perl
use tRNAscanSE::ExtendedSprinzlResults;

# For a single tRNA
write_extended_position_map($global_vars, "single_trna_positions.tsv", $tRNA);

# For all tRNAs in collection
write_extended_position_db($global_vars, "all_trna_positions.tsv");

# With genomic coordinates
write_genomic_coordinate_map($global_vars, "genomic_positions.tsv");
```

## Output Format Example

```
tRNA_ID         Seq_Pos Nucleotide  Align_Pos   Sprinzl_Pos Isotype Anticodon   Score   Genomic_Start   Genomic_End Strand  Sequence_Name
chr1.trna1      1       G           0           1           Phe     GAA         85.2    1000            1072        +       chr1
chr1.trna1      2       C           1           2           Phe     GAA         85.2    1000            1072        +       chr1
chr1.trna1      3       G           2           3           Phe     GAA         85.2    1000            1072        +       chr1
chr2.trna5      1       G           0           1           Leu     CUN         78.9    5500            5571        -       chr2
...
```

## Understanding the Coordinate Systems

### 1. Sequence Position (`Seq_Pos`)
- 1-based numbering
- Position within the individual tRNA sequence
- Corresponds to the original nucleotide order

### 2. Alignment Position (`Align_Pos`)  
- 0-based numbering
- Position in the covariance model alignment
- May include gaps (insertions/deletions relative to the model)
- Only non-gap positions are output in the mapping

### 3. Sprinzl Position (`Sprinzl_Pos`)
- Standardized tRNA structural numbering system
- Positions like "8", "14:i1" (insertion), "20a", etc.
- Allows comparison across different tRNA sequences
- Based on the canonical cloverleaf secondary structure

### 4. Genomic Coordinate (`Genomic_Coord`)
- Actual position in the input sequence/chromosome
- Accounts for strand orientation
- For forward strand: `genomic_start + seq_pos - 1`
- For reverse strand: `genomic_end - seq_pos + 1`

## Applications

This extended mapping is useful for:

1. **Comparative analysis** - Compare homologous positions across different tRNAs
2. **Structural analysis** - Map experimental data to structural positions
3. **Evolutionary studies** - Track mutations at specific structural sites
4. **Functional analysis** - Correlate sequence variations with function
5. **Database integration** - Link tRNA positions to other genomic annotations

## Implementation Notes

### Sprinzl Position Handling
The system properly handles:
- Standard positions (1, 2, 3, ...)
- Variable positions (44, 45, 46 in V-loop)
- Insertion positions (8:i1, 14:i1, etc.)
- Modified positions (20a, 20b for long variable arms)

### Gap Handling
- Alignment gaps (insertions relative to the covariance model) are skipped
- Only actual nucleotides from the sequence are included in output
- Alignment position preserves the original CM alignment coordinates

### Memory Considerations
- For large datasets, consider processing in batches
- The genomic coordinate mapping requires additional memory for coordinate calculations
- Temporary files can be kept for debugging using the `--keep-temp` option

## Troubleshooting

### Common Issues

1. **No tRNAs found**: Check input file format and search parameters
2. **Module not found**: Ensure the ExtendedSprinzlResults.pm is in the correct lib directory
3. **Permission errors**: Check write permissions for output directory
4. **Memory issues**: Process large files in smaller batches

### Debug Mode
Use `--keep-temp` to preserve intermediate files for troubleshooting:
```bash
perl scripts/extended_position_mapping.pl --keep-temp input.fa output.tsv
```

## Integration with Existing Workflows

This system integrates seamlessly with existing tRNAscan-SE workflows. You can:

1. Run standard tRNAscan-SE analysis
2. Use the saved results with the extended mapping functions  
3. Combine with other tRNAscan-SE output formats
4. Export results for analysis in R, Python, or other tools

## Future Enhancements

Potential improvements include:
- Direct integration with main tRNAscan-SE executable
- Support for additional coordinate systems
- Export to common bioinformatics formats (BED, GFF)
- Integration with visualization tools
- Support for modified nucleotides and pseudouridines