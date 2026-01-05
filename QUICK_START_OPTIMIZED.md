# Quick Start - Optimized Implementation

## ✨ New: 12.5× Faster Direct Interval Sampling

The optimized implementation uses direct interval sampling instead of generate-then-filter, eliminating 95% of wasted computations.

---

## Generate 100,000 Exomes with 150k Variants

### Single Command (Recommended):

```bash
./target/release/synthetic-vcf-generator-rst generate-batch \
  --num-vcfs 100000 \
  --output-dir ukbb-exomes-100k \
  --num-rows 150000 \
  --bed-file S31285117_Covered.bed \
  --model-vcf ukbb-model.vcf \
  --seed 42 \
  --num-threads 8
```

### Expected Results:
- **Time**: ~17.8 hours on 8 cores (12.5× faster than before!)
- **Output**: 100,000 VCF files with ~150k variants each
- **Storage**: ~435 GB total
- **Files per sample**:
  - `.vcf.gz` - compressed VCF (~4.2 MB)
  - `.vcf.gz.tbi` - tabix index (~155 KB)
  - `.meta.csv` - metadata (sample UUID, variant count, sex, etc.)

---

## Key Parameters

### `--num-rows` (Required for 150k variants)
- Set to `150000` for approximately 150k variants per sample
- This is the TOTAL variants to generate across all chromosomes
- Variants are automatically distributed by chromosome coverage

### `--num-vcfs` 
- Number of samples to generate
- Default: 100,000

### `--output-dir`
- Directory for output files
- Will be created if doesn't exist
- Default: `ukbb-exomes`

### `--seed`
- Base random seed for reproducibility
- Each sample gets a unique seed (base_seed + index)
- Default: 42

### `--num-threads`
- Number of parallel workers
- Default: all available CPU cores
- Recommendation: Use all cores for fastest generation

---

## Monitoring Progress

### Check generation status:
```bash
# Count completed samples
ls -1 ukbb-exomes-100k/*.vcf.gz | wc -l

# View metadata for completed samples
cat ukbb-exomes-100k/*.meta.csv | tail -20

# Check batch metadata
cat ukbb-exomes-100k/batch_metadata.json
```

### Verify variant counts:
```bash
# Check a few samples
for f in ukbb-exomes-100k/*.vcf.gz | head -5; do
  echo -n "$(basename $f): "
  zcat $f | grep -v "^#" | wc -l
done
```

### Monitor disk usage:
```bash
du -sh ukbb-exomes-100k
```

---

## Running in Background

For long-running jobs, use `screen` or `tmux`:

```bash
# Start screen session
screen -S vcf-generation

# Run generation
./target/release/synthetic-vcf-generator-rst generate-batch \
  --num-vcfs 100000 \
  --output-dir ukbb-exomes-100k \
  --num-rows 150000 \
  --bed-file S31285117_Covered.bed \
  --seed 42

# Detach: Ctrl+A, then D
# Reattach later: screen -r vcf-generation
```

---

## Output Files

### Directory structure:
```
ukbb-exomes-100k/
├── batch_metadata.json                          # Batch-level metadata
├── 00000000-0000-0000-0000-000000000001.vcf.gz  # Sample 1 VCF
├── 00000000-0000-0000-0000-000000000001.vcf.gz.tbi
├── 00000000-0000-0000-0000-000000000001.meta.csv
├── 00000000-0000-0000-0000-000000000002.vcf.gz  # Sample 2 VCF
├── 00000000-0000-0000-0000-000000000002.vcf.gz.tbi
├── 00000000-0000-0000-0000-000000000002.meta.csv
└── ... (100,000 samples total)
```

### batch_metadata.json:
```json
{
  "batch_id": "uuid",
  "num_samples": 100000,
  "model_vcf_path": "ukbb-model.vcf",
  "bed_file_path": "S31285117_Covered.bed",
  "num_rows_per_chromosome": 150000,
  "base_seed": 42,
  "num_threads": 8,
  "output_directory": "ukbb-exomes-100k",
  "type": "exome"
}
```

### Sample .meta.csv:
```csv
sample_uuid,vcf_file,num_variants,sex,type,chrX_variants,chrY_variants,seed
uuid,uuid.vcf.gz,150001,male,exome,5546,536,42
```

---

## Variant Count Details

### Expected counts:
- **Males**: ~150,000 variants (includes chrY)
- **Females**: ~149,465 variants (no chrY)

The difference is biologically accurate - females don't have Y chromosomes.

### Distribution by chromosome:
Variants are distributed proportionally to each chromosome's BED coverage:
- chr1: ~11,500 variants (largest)
- chr22: ~1,800 variants (smallest autosome)
- chrX: ~5,500 variants
- chrY: ~536 variants (males only)

---

## Performance Scaling

| CPU Cores | Time for 100k samples | Time per sample |
|-----------|----------------------|-----------------|
| 8 cores   | 17.8 hours          | 0.64 seconds    |
| 16 cores  | 8.9 hours           | 0.32 seconds    |
| 32 cores  | 4.4 hours           | 0.16 seconds    |
| 64 cores  | 2.2 hours           | 0.08 seconds    |

---

## Troubleshooting

### "Not enough disk space"
- Check available space: `df -h`
- Need at least 500 GB free
- Clean up old test directories: `rm -rf test-*`

### "Out of memory"
- The optimized implementation uses minimal memory
- Each worker uses ~100-200 MB
- Should work fine with 8GB+ RAM

### "Generation seems stuck"
- Check if compression tools are available: `which bgzip tabix`
- Monitor CPU usage: `htop` or `top`
- Check for disk I/O bottleneck: `iostat -x 2`

### Variant counts slightly off
- Expected variation: ±500 variants due to:
  - Sex-specific chromosomes (males have chrY, females don't)
  - Rounding in distribution calculation
  - Deduplication (rare position collisions)

---

## Advanced Usage

### Generate specific number of samples:
```bash
# Just 1,000 samples for testing
--num-vcfs 1000

# Specific range (use different seeds)
--num-vcfs 10000 --seed 0      # Samples 0-9,999
--num-vcfs 10000 --seed 10000  # Samples 10,000-19,999
```

### Custom variant density:
```bash
# More variants (~200k per sample)
--num-rows 200000

# Fewer variants (~100k per sample)
--num-rows 100000
```

### Different seed for reproducibility:
```bash
# Each seed produces different random variants
--seed 12345
```

---

## What's Different in the Optimized Version?

### Old Approach (Generate-then-Filter):
1. Generate 400k random positions per chromosome
2. Filter to keep only positions in BED intervals
3. Keep ~6,600 positions per chromosome (5.17% hit rate)
4. Write variants
5. ❌ 95% waste

### New Approach (Direct Sampling):
1. Sample positions directly from BED intervals
2. Use weighted sampling (larger intervals = more variants)
3. 100% hit rate - all positions are valid
4. Write variants
5. ✅ Zero waste, 12.5× faster!

---

## Next Steps

1. **Test with small batch**: Start with 100-1000 samples to verify
2. **Monitor first hour**: Ensure stable performance and disk space
3. **Run full batch**: Generate all 100,000 samples
4. **Validate output**: Check variant counts and file integrity
5. **Use the data**: Ready for TileDB ingestion or other analysis!

For more details, see:
- `OPTIMIZATION_RESULTS.md` - Performance benchmarks
- `100K_SAMPLES_PLAN.md` - Full planning details
- `BED_FILTERING_EFFICIENCY_ANALYSIS.md` - Technical analysis







