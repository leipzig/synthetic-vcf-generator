use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use rand::Rng;

#[derive(Debug, Clone)]
pub struct BedIntervals {
    intervals: HashMap<String, Vec<(usize, usize)>>,
}

impl BedIntervals {
    pub fn new<P: AsRef<Path>>(bed_file_path: P) -> Result<Self> {
        let mut intervals: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
        let content = fs::read_to_string(&bed_file_path)
            .with_context(|| format!("Failed to read BED file: {:?}", bed_file_path.as_ref()))?;

        for (line_num, line) in content.lines().enumerate() {
            let line = line.trim();
            // Skip empty lines and comments
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                anyhow::bail!(
                    "Invalid BED format at line {}: expected at least 3 fields, got {}",
                    line_num + 1,
                    fields.len()
                );
            }

            let chrom = fields[0].to_string();
            let start: usize = fields[1]
                .parse()
                .with_context(|| format!("Invalid start coordinate at line {}", line_num + 1))?;
            let end: usize = fields[2]
                .parse()
                .with_context(|| format!("Invalid end coordinate at line {}", line_num + 1))?;

            if start >= end {
                anyhow::bail!(
                    "Invalid interval at line {}: start={} >= end={}",
                    line_num + 1,
                    start,
                    end
                );
            }

            intervals.entry(chrom).or_insert_with(Vec::new).push((start, end));
        }

        // Sort intervals for each chromosome by start position
        for intervals_list in intervals.values_mut() {
            intervals_list.sort();
        }

        Ok(Self { intervals })
    }

    pub fn get_chromosomes(&self) -> Vec<String> {
        let mut chroms: Vec<String> = self.intervals.keys().cloned().collect();
        chroms.sort();
        chroms
    }

    pub fn get_intervals(&self, chromosome: &str) -> &[(usize, usize)] {
        self.intervals.get(chromosome).map(|v| v.as_slice()).unwrap_or(&[])
    }

    pub fn contains_position(&self, chromosome: &str, position: usize) -> bool {
        // Convert VCF 1-based position to BED 0-based coordinate
        let bed_pos = position.saturating_sub(1);

        if let Some(intervals) = self.intervals.get(chromosome) {
            for (start, end) in intervals {
                if *start <= bed_pos && bed_pos < *end {
                    return true;
                }
                if *start > bed_pos {
                    // Since sorted, no need to check further
                    break;
                }
            }
        }
        false
    }

    pub fn get_valid_positions(&self, chromosome: &str, positions: &[usize]) -> Vec<usize> {
        positions
            .iter()
            .filter(|pos| self.contains_position(chromosome, **pos))
            .copied()
            .collect()
    }

    pub fn get_total_bases(&self, chromosome: Option<&str>) -> usize {
        if let Some(chrom) = chromosome {
            self.intervals
                .get(chrom)
                .map(|intervals| {
                    intervals
                        .iter()
                        .map(|(start, end)| end - start)
                        .sum::<usize>()
                })
                .unwrap_or(0)
        } else {
            self.intervals
                .values()
                .map(|intervals| {
                    intervals
                        .iter()
                        .map(|(start, end)| end - start)
                        .sum::<usize>()
                })
                .sum()
        }
    }

    /// Sample positions directly from BED intervals using weighted sampling
    /// This is much more efficient than generate-then-filter approach
    pub fn sample_positions<R: rand::Rng>(
        &self,
        chromosome: &str,
        count: usize,
        rng: &mut R,
    ) -> Vec<usize> {
        use std::collections::HashSet;
        
        let intervals = match self.intervals.get(chromosome) {
            Some(ivs) if !ivs.is_empty() => ivs,
            _ => return Vec::new(),
        };

        // Calculate total bases in all intervals
        let total_bases: usize = intervals.iter().map(|(s, e)| e - s).sum();
        if total_bases == 0 {
            return Vec::new();
        }

        // Build cumulative weights for weighted sampling
        let mut cumulative_weights = Vec::with_capacity(intervals.len());
        let mut cumsum = 0usize;
        for (start, end) in intervals.iter() {
            cumsum += end - start;
            cumulative_weights.push(cumsum);
        }

        // Sample positions with deduplication
        let mut positions = HashSet::with_capacity(count);
        let max_attempts = count * 3; // Prevent infinite loops
        let mut attempts = 0;

        while positions.len() < count && attempts < max_attempts {
            attempts += 1;

            // Sample which interval to use (weighted by interval size)
            let weight = rng.gen_range(0..total_bases);
            let interval_idx = cumulative_weights
                .binary_search(&weight)
                .unwrap_or_else(|idx| idx);
            
            let (start, end) = intervals[interval_idx];
            
            // Sample position within the interval (convert to 1-based VCF coordinates)
            let pos = rng.gen_range(start..end) + 1;
            positions.insert(pos);
        }

        // Convert to sorted vector
        let mut result: Vec<usize> = positions.into_iter().collect();
        result.sort_unstable();
        result
    }

    /// Calculate how many variants should be generated for each chromosome
    /// based on their BED interval coverage
    pub fn distribute_variants_by_coverage(
        &self,
        total_variants: usize,
    ) -> HashMap<String, usize> {
        let total_bases = self.get_total_bases(None);
        if total_bases == 0 {
            return HashMap::new();
        }

        let mut distribution = HashMap::new();
        for chrom in self.get_chromosomes() {
            let chrom_bases = self.get_total_bases(Some(&chrom));
            let chrom_variants = ((total_variants as f64 * chrom_bases as f64) 
                / total_bases as f64).round() as usize;
            distribution.insert(chrom, chrom_variants);
        }
        
        distribution
    }
}

