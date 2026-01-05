use crate::bed_parser::BedIntervals;
use crate::sample_sheet::SampleInfo;
use crate::vcf_model::VCFModel;
use anyhow::Result;
use rand::seq::SliceRandom;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng; // For position sampling (matches Python's random.Random)
use rand_pcg::Pcg32; // For other random values (matches fastrand)
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use uuid::Uuid;

// Standard human chromosome sizes (GRCh38)
fn get_chromosome_sizes() -> HashMap<String, usize> {
    let mut sizes = HashMap::new();
    sizes.insert("chr1".to_string(), 248_956_422);
    sizes.insert("chr2".to_string(), 242_193_529);
    sizes.insert("chr3".to_string(), 198_295_559);
    sizes.insert("chr4".to_string(), 190_214_555);
    sizes.insert("chr5".to_string(), 181_538_259);
    sizes.insert("chr6".to_string(), 170_805_979);
    sizes.insert("chr7".to_string(), 159_345_973);
    sizes.insert("chr8".to_string(), 145_138_636);
    sizes.insert("chr9".to_string(), 138_394_717);
    sizes.insert("chr10".to_string(), 133_797_422);
    sizes.insert("chr11".to_string(), 135_086_622);
    sizes.insert("chr12".to_string(), 133_275_309);
    sizes.insert("chr13".to_string(), 114_364_328);
    sizes.insert("chr14".to_string(), 107_043_718);
    sizes.insert("chr15".to_string(), 101_991_189);
    sizes.insert("chr16".to_string(), 90_338_345);
    sizes.insert("chr17".to_string(), 83_257_441);
    sizes.insert("chr18".to_string(), 80_373_285);
    sizes.insert("chr19".to_string(), 58_617_616);
    sizes.insert("chr20".to_string(), 64_444_167);
    sizes.insert("chr21".to_string(), 46_709_983);
    sizes.insert("chr22".to_string(), 50_818_468);
    sizes.insert("chrX".to_string(), 156_040_895);
    sizes.insert("chrY".to_string(), 57_227_415);
    sizes
}

pub struct VCFGenerator {
    model: VCFModel,
    bed_intervals: Option<BedIntervals>,
    chromosomes: Vec<String>,
    chromosome_sizes: HashMap<String, usize>,
    num_rows: usize,
    seed: i64,
    sample_uuid: Uuid,
    sex: String,
    is_whole_genome: bool,
}

impl VCFGenerator {
    pub fn new(
        model_path: String,
        bed_path: Option<String>,
        num_rows: Option<usize>,
        seed: i64,
        sample_uuid: Option<Uuid>,
        sex: Option<String>,
    ) -> Result<Self> {
        let (bed_intervals, chromosomes, chromosome_sizes, is_whole_genome) = if let Some(bed_path_str) = bed_path {
            // Exome mode: use BED file
            let bed_intervals = BedIntervals::new(&bed_path_str)?;
            let mut chroms = bed_intervals.get_chromosomes();
            chroms.sort();
            let sizes = HashMap::new(); // Not used in exome mode
            (Some(bed_intervals), chroms, sizes, false)
        } else {
            // Whole genome mode: use standard chromosome sizes
            let sizes = get_chromosome_sizes();
            let mut chroms: Vec<String> = sizes.keys().cloned().collect();
            chroms.sort();
            (None, chroms, sizes, true)
        };

        // Calculate num_rows if not provided
        let calculated_rows = if num_rows.is_none() {
            if is_whole_genome {
                // For whole genome, use a reasonable default (e.g., 1 variant per 10kb)
                let total_bases: usize = chromosome_sizes.values().sum();
                (total_bases / 10_000).max(100)
            } else {
                // For exome, this is now the TOTAL variants desired (not per chromosome)
                // Default: ~150k variants for typical exome
                let total_bases = bed_intervals.as_ref().unwrap().get_total_bases(None);
                // Aim for 1 variant per 300-350 bases in exome regions
                (total_bases / 330).max(100)
            }
        } else {
            num_rows.unwrap()
        };

        // Load VCF model
        let model = VCFModel::new(&model_path)?;

        let uuid = sample_uuid.unwrap_or_else(Uuid::new_v4);

        // Determine sex if not provided (use UUID hash like Python)
        let sample_sex = sex.unwrap_or_else(|| {
            // Hash UUID to get deterministic sex assignment
            use std::collections::hash_map::DefaultHasher;
            use std::hash::{Hash, Hasher};
            let mut hasher = DefaultHasher::new();
            uuid.hash(&mut hasher);
            let hash = hasher.finish();
            if hash % 2 == 0 {
                "male".to_string()
            } else {
                "female".to_string()
            }
        });

        Ok(Self {
            model,
            bed_intervals,
            chromosomes,
            chromosome_sizes,
            num_rows: calculated_rows,
            seed,
            sample_uuid: uuid,
            sex: sample_sex,
            is_whole_genome,
        })
    }

    pub fn generate(&self, output_path: &Path) -> Result<SampleInfo> {
        // Use StdRng (MT19937) for position sampling to match Python's random.Random
        let mut position_rng = StdRng::seed_from_u64(self.seed as u64);
        // Use PCG32 for other random values to match fastrand
        let mut fastrand_rng = Pcg32::seed_from_u64(self.seed as u64);
        
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);

        // Write VCF header
        self.write_header(&mut writer)?;

        // Generate variants per chromosome
        let mut variant_count = 0;
        let mut variant_counts_by_chrom: HashMap<String, usize> = HashMap::new();

        for chromosome in &self.chromosomes {
            // Skip chrY for females
            if chromosome == "chrY" && self.sex == "female" {
                continue;
            }

            // Calculate number of variants for this chromosome
            let variants_per_chrom = if self.is_whole_genome {
                // Distribute variants proportionally across chromosomes
                let total_genome_size: usize = self.chromosome_sizes.values().sum();
                let chrom_size = self.chromosome_sizes.get(chromosome).unwrap_or(&100_000_000);
                ((self.num_rows as f64 * *chrom_size as f64) / total_genome_size as f64).ceil() as usize
            } else {
                // For exome mode, distribute based on BED coverage
                let bed = self.bed_intervals.as_ref().unwrap();
                let chrom_bases = bed.get_total_bases(Some(chromosome));
                let total_bases = bed.get_total_bases(None);
                if total_bases == 0 {
                    0
                } else {
                    ((self.num_rows as f64 * chrom_bases as f64) / total_bases as f64).round() as usize
                }
            };

            if variants_per_chrom == 0 {
                continue;
            }

            // Sample positions directly (new efficient approach)
            let positions = if self.is_whole_genome {
                // Whole genome: sample from entire chromosome
                let chromosome_length = *self.chromosome_sizes.get(chromosome).unwrap_or(&100_000_000);
                let sample_size = variants_per_chrom.min(chromosome_length);
                
                use rand::seq::IteratorRandom;
                let mut pos_vec: Vec<usize> = (1..=chromosome_length)
                    .choose_multiple(&mut position_rng, sample_size);
                pos_vec.sort_unstable();
                pos_vec
            } else {
                // Exome mode: sample directly from BED intervals (efficient!)
                self.bed_intervals.as_ref().unwrap()
                    .sample_positions(chromosome, variants_per_chrom, &mut position_rng)
            };

            // Generate variant for each position
            for pos in positions.iter() {
                // Skip chrY for females (double-check)
                if chromosome == "chrY" && self.sex == "female" {
                    continue;
                }

                variant_count += 1;
                *variant_counts_by_chrom
                    .entry(chromosome.clone())
                    .or_insert(0) += 1;

                // Generate variant line (matching Python logic)
                self.write_variant(&mut writer, chromosome, *pos, variant_count - 1, &mut fastrand_rng)?;
            }
        }

        writer.flush()?;

        Ok(SampleInfo {
            sample_uuid: self.sample_uuid,
            num_variants: variant_count,
            sex: self.sex.clone(),
            is_exome: !self.is_whole_genome,
            variant_counts_by_chrom,
            seed: self.seed,
        })
    }

    fn write_header(&self, writer: &mut BufWriter<File>) -> Result<()> {
        writeln!(writer, "##fileformat=VCFv4.2")?;
        writeln!(writer, "##source=VirtualVCF-UKBB-Model 0.1.0")?;

        // Write model header lines
        for line in self.model.get_header_lines() {
            writeln!(writer, "{}", line)?;
        }

        // Write contig lines
        for chrom in &self.chromosomes {
            writeln!(writer, "##contig=<ID={}>", chrom)?;
        }

        // Write column header
        writeln!(
            writer,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
            self.sample_uuid
        )?;

        Ok(())
    }

    fn write_variant(
        &self,
        writer: &mut BufWriter<File>,
        chrom: &str,
        pos: usize,
        variant_count: usize,
        rng: &mut Pcg32,
    ) -> Result<()> {
        // Generate REF/ALT - support SNVs, insertions, and deletions
        let alleles = vec!["A", "C", "G", "T"];
        
        // Decide variant type: 70% SNV, 15% insertion, 15% deletion
        let variant_type_roll = rng.gen_range(0..100);
        let (ref_allele, alt_allele) = if variant_type_roll < 70 {
            // SNV: single nucleotide substitution
            let ref_idx = rng.gen_range(0..4);
            let alt_idx = (ref_idx + rng.gen_range(1..4)) % 4;
            (alleles[ref_idx].to_string(), alleles[alt_idx].to_string())
        } else if variant_type_roll < 85 {
            // Insertion: REF is single base, ALT is REF + inserted sequence
            let ref_idx = rng.gen_range(0..4);
            let ref_base = alleles[ref_idx];
            // Insert 1-10 bases
            let insert_len = rng.gen_range(1..=10);
            let mut alt_seq = ref_base.to_string();
            for _ in 0..insert_len {
                let insert_idx = rng.gen_range(0..4);
                alt_seq.push_str(alleles[insert_idx]);
            }
            (ref_base.to_string(), alt_seq)
        } else {
            // Deletion: REF is multiple bases, ALT is single base anchor
            let ref_idx = rng.gen_range(0..4);
            let anchor_base = alleles[ref_idx];
            // Delete 1-10 bases (REF = anchor + deleted sequence)
            let del_len = rng.gen_range(1..=10);
            let mut ref_seq = anchor_base.to_string();
            for _ in 0..del_len {
                let del_idx = rng.gen_range(0..4);
                ref_seq.push_str(alleles[del_idx]);
            }
            (ref_seq, anchor_base.to_string())
        };
        
        let qual = rng.gen_range(10..=100);

        // Generate genotype (from VirtualVCF's sample generation)
        let gt = self.generate_genotype(chrom, rng);

        // Adjust genotype for sex-specific chromosomes (matching Python)
        let mut adjusted_gt = gt.clone();
        if chrom == "chrX" && self.sex == "male" {
            // Males: convert to hemizygous
            if gt == "0|0" || gt == "0/0" {
                adjusted_gt = "0".to_string();
            } else if gt == "1|1" || gt == "1/1" {
                adjusted_gt = "1".to_string();
            } else {
                // Heterozygous -> hemizygous alternate
                adjusted_gt = "1".to_string();
            }
        }
        if chrom == "chrY" && self.sex == "male" {
            if gt == "0|0" || gt == "0/0" {
                adjusted_gt = "0".to_string();
            } else if gt == "1|1" || gt == "1/1" {
                adjusted_gt = "1".to_string();
            } else {
                adjusted_gt = "1".to_string();
            }
        }

        // Build UKBB INFO fields (matching Python exactly)
        let mut ukbb_info = Vec::new();
        ukbb_info.push("AC=1".to_string());
        ukbb_info.push("AF=0.50".to_string());
        ukbb_info.push("AN=2".to_string());
        ukbb_info.push("DP=80".to_string());
        ukbb_info.push("MQ=58.2".to_string());
        ukbb_info.push("MQRankSum=0.72".to_string());
        ukbb_info.push("QD=25.3".to_string());
        ukbb_info.push("ReadPosRankSum=0.15".to_string());
        ukbb_info.push("FS=3.2".to_string());
        ukbb_info.push("SOR=1.1".to_string());
        ukbb_info.push(format!("ALLELE_ID=A{}", variant_count + 1));
        ukbb_info.push("TARGETED".to_string());

        // HML
        let hml = if adjusted_gt == "0/0" || adjusted_gt == "0|0" || adjusted_gt == "0" {
            "hom"
        } else if adjusted_gt == "1/1" || adjusted_gt == "1|1" || adjusted_gt == "1" {
            "hom_alt"
        } else {
            "het"
        };
        ukbb_info.push(format!("HML={}", hml));

        ukbb_info.push("JIDS=S1".to_string());
        ukbb_info.push("MOSAIC=0.02".to_string());
        ukbb_info.push("SoftClipRatio=0.01".to_string());
        ukbb_info.push("FractionInformativeReads=0.95".to_string());

        // Optional flags
        if variant_count % 3 == 0 {
            ukbb_info.push("Recombinant".to_string());
        }

        ukbb_info.push(format!("PS={}", 100 + variant_count * 10));
        ukbb_info.push(format!("EVENT=E{}", variant_count + 1));

        // EVENTTYPE
        let event_type = if alt_allele.len() > ref_allele.len() {
            "INS"
        } else if alt_allele.len() < ref_allele.len() {
            "DEL"
        } else {
            "SNV"
        };
        ukbb_info.push(format!("EVENTTYPE={}", event_type));

        let info_str = ukbb_info.join(";");

        // Generate FORMAT fields (matching Python exactly)
        let format_string = self.model.get_format_string();
        let format_fields: Vec<&str> = format_string.split(':').collect();
        let mut format_values = Vec::new();

        for field in &format_fields {
            if *field == "GT" {
                format_values.push(adjusted_gt.clone());
            } else if *field == "AD" {
                if adjusted_gt == "0/0" || adjusted_gt == "0|0" || adjusted_gt == "0" {
                    format_values.push("45,5".to_string());
                } else if adjusted_gt == "1/1" || adjusted_gt == "1|1" || adjusted_gt == "1" {
                    format_values.push("5,45".to_string());
                } else {
                    format_values.push("35,35".to_string());
                }
            } else if *field == "DP" {
                format_values.push("80".to_string());
            } else if *field == "AF" {
                if adjusted_gt == "0/0" || adjusted_gt == "0|0" || adjusted_gt == "0" {
                    format_values.push("0.93,0.07".to_string());
                } else if adjusted_gt == "1/1" || adjusted_gt == "1|1" || adjusted_gt == "1" {
                    format_values.push("0.07,0.93".to_string());
                } else {
                    format_values.push("0.44,0.56".to_string());
                }
            } else if *field == "GQ" {
                format_values.push("99".to_string());
            } else if *field == "PL" {
                if adjusted_gt == "0/0" || adjusted_gt == "0|0" || adjusted_gt == "0" {
                    format_values.push("0,120,890".to_string());
                } else if adjusted_gt == "1/1" || adjusted_gt == "1|1" || adjusted_gt == "1" {
                    format_values.push("1120,145,0".to_string());
                } else {
                    format_values.push("0,120,890".to_string());
                }
            } else if *field == "GP" {
                format_values.push("10,0,990".to_string());
            } else if *field == "F1R2" {
                format_values.push("32,28".to_string());
            } else if *field == "F2R1" {
                format_values.push("25,20".to_string());
            } else if *field == "SB" {
                format_values.push("0,1,39,20".to_string());
            } else if *field == "MB" {
                format_values.push("0.12".to_string());
            } else if *field == "SQ" {
                format_values.push("28.5".to_string());
            } else if *field == "PRI" {
                format_values.push("50".to_string());
            } else if *field == "PS" {
                format_values.push(format!("{}", 100 + variant_count * 10));
            } else if *field == "QL" {
                format_values.push("0.88".to_string());
            } else if field.starts_with("J") {
                if *field == "JAD" {
                    format_values.push("70,50".to_string());
                } else if *field == "JAF" {
                    format_values.push("0.15,0.32".to_string());
                } else if *field == "JDP" {
                    format_values.push("120".to_string());
                } else if *field == "JGQ" {
                    format_values.push("98".to_string());
                } else if *field == "JGT" {
                    format_values.push(adjusted_gt.clone());
                } else if *field == "JPL" {
                    format_values.push("0,145,1120".to_string());
                } else if *field == "JQL" {
                    format_values.push("0.88".to_string());
                } else if *field == "JVQL" {
                    format_values.push("0.88".to_string());
                } else {
                    format_values.push(".".to_string());
                }
            } else {
                format_values.push(".".to_string());
            }
        }

        let format_data = format_values.join(":");

        writeln!(
            writer,
            "{}\t{}\t.\t{}\t{}\t{}\tPASS\t{}\t{}\t{}",
            chrom, pos, ref_allele, alt_allele, qual, info_str, format_string, format_data
        )?;

        Ok(())
    }

    fn generate_genotype(&self, _chrom: &str, rng: &mut Pcg32) -> String {
        // Match Python's VirtualVCF sample generation logic
        // Python uses phased=True and rotates available_samples
        // Simplified: generate phased genotypes with weights
        let val = rng.gen_range(0..1000);
        if val < 700 {
            "0|0".to_string()
        } else if val < 850 {
            "1|0".to_string()
        } else if val < 950 {
            "0|1".to_string()
        } else {
            "1|1".to_string()
        }
    }
}
