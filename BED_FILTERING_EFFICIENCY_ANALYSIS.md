# BED Filtering Efficiency Analysis

## Current Approach (Post-filtering)

### How it works:
1. Generate 400,000 random positions per chromosome (range: 1 to 40,000,000)
2. Filter positions to keep only those within BED intervals
3. Write variants for valid positions

### Efficiency metrics:
- **Exome coverage**: 49,668,806 bases
- **Sampling space**: 960,000,000 bases (40M × 24 chromosomes)
- **Hit rate**: 5.17% ⚠️
- **Waste**: 94.83% of generated positions are discarded

### Cost for 100,000 samples:
- **Positions generated**: 9,600,000,000 (9.6 billion)
- **Valid positions kept**: ~497,000,000 (497 million)
- **Wasted computations**: ~9,100,000,000 (9.1 billion position checks)

### Current performance:
- Single sample: ~159,000 variants
- Generation time: Fast (Rust is efficient even with waste)
- Memory usage: Moderate (positions generated per chromosome)

---

## Recommended Approach (Direct Sampling)

### How it would work:
1. Pre-compute flattened list of all valid positions in BED intervals
2. Sample directly from valid positions only
3. Distribute samples across chromosomes based on their exome coverage
4. Write variants immediately without filtering

### Efficiency improvements:
- **Hit rate**: 100% ✓
- **Waste**: 0%
- **Positions generated**: Exactly what's needed (~15,900,000,000 for 100k samples)

### Implementation approaches:

#### Option 1: Weighted interval sampling (Best for large-scale)
```rust
// For each chromosome:
// 1. Calculate variants per chromosome based on interval coverage
let chrom_coverage = bed_intervals.get_total_bases(Some(chromosome));
let total_coverage = bed_intervals.get_total_bases(None);
let variants_for_chrom = (target_variants * chrom_coverage) / total_coverage;

// 2. Sample positions directly from intervals
let intervals = bed_intervals.get_intervals(chromosome);
for _ in 0..variants_for_chrom {
    // Pick random interval (weighted by size)
    let interval = pick_weighted_interval(intervals, rng);
    // Pick random position within interval
    let pos = rng.gen_range(interval.start..interval.end);
    generate_variant(pos);
}
```

#### Option 2: Pre-build position pool (Memory intensive)
```rust
// Build once at startup:
let mut valid_positions: HashMap<String, Vec<usize>> = HashMap::new();
for (chrom, intervals) in bed_intervals {
    for (start, end) in intervals {
        valid_positions[chrom].extend(start..end);
    }
}

// Then sample directly:
let positions = valid_positions[chrom].choose_multiple(rng, count);
```

#### Option 3: Hybrid approach (Balanced) ⭐ RECOMMENDED
```rust
// Use interval sampling with deduplication
let mut positions = HashSet::new();
let intervals = bed_intervals.get_intervals(chromosome);
let total_interval_length: usize = intervals.iter().map(|(s,e)| e-s).sum();

// Calculate cumulative weights for intervals
let cumulative_weights = build_cumulative_weights(intervals);

while positions.len() < target_count {
    // Sample interval proportional to its size
    let interval_idx = sample_weighted(&cumulative_weights, rng);
    let (start, end) = intervals[interval_idx];
    
    // Sample position within interval
    let pos = rng.gen_range(start..=end);
    positions.insert(pos);
}
```

---

## Comparison for 100,000 samples with 150k variants each

| Metric | Current (Post-filter) | Direct Sampling |
|--------|----------------------|-----------------|
| Positions generated | 9.6 billion | 15.9 billion |
| Position checks | 9.6 billion | 0 |
| Hit rate | 5.17% | 100% |
| Wasted work | 94.83% | 0% |
| Memory per sample | Low | Low-Moderate |
| Implementation complexity | Simple | Moderate |
| Total variants generated | 15.9 billion | 15.9 billion |

### Expected performance improvement:
- **CPU time**: ~15-20% faster (eliminating filtering overhead)
- **Predictability**: Exact variant counts (no randomness in filtering)
- **Scalability**: Better for larger sample counts

---

## Recommendation

For generating **100,000 samples**:

### Short term (Use current implementation):
✓ Current approach is acceptable because:
- Rust implementation is fast enough even with 95% waste
- Code is simple and proven to work
- Position filtering is cheap (just range checks)
- Main bottleneck is likely I/O (writing/compressing files)

### Long term (Optimize for scale):
⚠ Implement direct interval sampling if:
- Generating >500k samples
- Running on limited CPU resources
- Need exact variant counts (not approximate)
- BED coverage drops below 5% (even more waste)

---

## Immediate Action

**For your 100k sample batch**, I recommend:

1. **Use current Rust implementation** - it's fast enough
2. **Run in parallel** - utilize all CPU cores
3. **Monitor performance** - measure actual time for first 1000 samples
4. **Consider optimization** - if it's too slow, implement direct sampling

### Estimated time for 100k samples (current approach):
- Single sample: ~2-3 seconds (generation + compression + indexing)
- With 32 cores: ~2-3 hours total
- With 64 cores: ~1-1.5 hours total
- Storage required: ~60 GB (600 KB × 100,000)







