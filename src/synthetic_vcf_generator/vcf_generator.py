from typing import BinaryIO, Literal

import sys
import uuid
import multiprocessing
from pathlib import Path

from synthetic_vcf_generator.virtual_vcf import VirtualVCF


def write_fileobj(
    fileobj: BinaryIO,
    virtual_vcf: VirtualVCF,
    output_type: Literal["vcf", "gzip", "bgzip"],
) -> None:
    """
    Writes VirtualVCF data to the given file object.

    Args:
        fileobj (BinaryIO): The open file object to write to.
        virtual_vcf (VirtualVCF): VirtualVCF object containing the data.
        output_type (Literal["vcf", "gzip", "bgzip"]): The type of the output data.
    """
    if output_type == "vcf":
        with virtual_vcf as v_vcf:
            for line in v_vcf:
                fileobj.write(line.encode(encoding="utf-8"))
    elif output_type == "gzip":
        import gzip

        with (
            gzip.GzipFile(fileobj=fileobj, mode="wb") as gzip_file,
            virtual_vcf as v_vcf,
        ):
            for line in v_vcf:
                gzip_file.write(line.encode(encoding="utf-8"))
    elif output_type == "bgzip":
        try:
            from Bio import bgzf

            with (
                bgzf.BgzfWriter(fileobj=fileobj, mode="wt") as bgzip_file,
                virtual_vcf as v_vcf,
            ):
                for line in v_vcf:
                    bgzip_file.write(line)
        except ImportError:
            raise RuntimeError(f"{output_type=} requires Biopython")
    else:
        raise ValueError(f"{output_type=} not supported")


def write_stdout(
    virtual_vcf: VirtualVCF, output_type: Literal["vcf", "gzip", "bgzip"]
) -> None:
    """
    Writes VirtualVCF data to the standard output.

    Args:
        virtual_vcf (VirtualVCF): VirtualVCF object containing the data.
        output_type (Literal["vcf", "gzip", "bgzip"]): The type of the output data.
    """
    stdout = sys.stdout.buffer
    write_fileobj(stdout, virtual_vcf, output_type)


def write_file(
    virtual_vcf: VirtualVCF,
    output_type: Literal["vcf", "gzip", "bgzip"],
    synthetic_vcf_path: Path,
) -> None:
    """
    Writes VirtualVCF data to a VCF file.

    Args:
        virtual_vcf (VirtualVCF): VirtualVCF object containing the data.
        output_type (Literal["vcf", "gzip", "bgzip"]): The type of the output data.
        synthetic_vcf_path (Path): Path to the output VCF file to be written.
    """
    with open(synthetic_vcf_path, "wb") as wb_file:
        write_fileobj(wb_file, virtual_vcf, output_type)


def synthetic_vcf_data(
    synthetic_vcf_path,
    output_type,
    num_rows,
    num_samples,
    chromosomes,
    seed,
    sample_prefix,
    id_type,
    phased,
    large_format,
    reference_dir_path,
    bed_intervals,
):
    """
    Generates synthetic VCF data and writes it to either a file or standard output.

    Args:
        synthetic_vcf_path (Path or None): Path to the synthetic VCF file or None to write to standard output.
        output_type (str): Type file to output ("vcp", "gzip", or "bgzip").
        num_rows (int): Number of rows (variants) to generate per chromosome.
        num_samples (int): Number of samples.
        chromosomes (List[str]): List of chromosome IDs to include in the VCF.
        seed (int): Random seed for reproducibility.
        sample_prefix (str): Prefix for sample names.
        id_type (str): Type of unique ID to use for samples ("count", "padded_count", or "uuid").
        phased (bool): Phased or unphased genotypes.
        large_format (bool): Use large format VCF.
        reference_dir_path (Path or None): Path to imported reference data.
    """
    virtual_vcf = VirtualVCF(
        num_rows=num_rows,
        num_samples=num_samples,
        chromosomes=chromosomes,
        sample_prefix=sample_prefix,
        id_type=id_type,
        random_seed=seed,
        phased=phased,
        large_format=large_format,
        reference_dir=reference_dir_path,
        bed_intervals=bed_intervals,
    )

    if synthetic_vcf_path is None:
        write_stdout(virtual_vcf, output_type)
    else:
        write_file(virtual_vcf, output_type, synthetic_vcf_path)


def batch_synthetic_vcf_data(
    synthetic_vcf_dir,
    output_type,
    num_vcfs,
    vcf_prefix,
    num_rows,
    num_samples,
    chromosomes,
    seed,
    sample_prefix,
    id_type,
    phased,
    large_format,
    reference_dir_path,
    bed_intervals,
    num_threads,
):
    """
    Generates synthetic VCF data and writes it to either a file or standard output.

    Args:
        synthetic_vcf_dir (Path): Path to the directory where output synthetic VCF files should be saved.
        output_type (str): Type file to output ("vcp", "gzip", or "bgzip").
        num_vcfs (int): Number of VCF files to generate.
        vcf_prefix (str): Prefix added to the filename of every generated VCF.
        num_rows (int): Number of rows (variants) to generate per chromosome.
        num_samples (int): Number of samples.
        chromosomes (List[str]): List of chromosome IDs to include in the VCF.
        seed (int): Random seed for reproducibility.
        sample_prefix (str): Prefix for sample names.
        id_type (str): Type of unique ID to use for samples.
        phased (bool): Phased or unphased genotypes.
        large_format (bool): Use large format VCF.
        reference_dir_path (Path or None): Path to imported reference data.
        num_threads (int): Number of threads.
    """

    if not synthetic_vcf_dir.exists():
        synthetic_vcf_dir.mkdir(parents=True)

    virtual_vcf = VirtualVCF(
        num_rows=num_rows,
        num_samples=num_samples,
        chromosomes=chromosomes,
        sample_prefix=sample_prefix,
        id_type=id_type,
        random_seed=seed,
        phased=phased,
        large_format=large_format,
        reference_dir=reference_dir_path,
        bed_intervals=bed_intervals,
    )

    ext = ".vcf.gz" if output_type in {"gzip", "bgzip"} else ".vcf"
    results = []
    with multiprocessing.Pool(num_threads) as pool:
        for i in range(num_vcfs):
            filepath = synthetic_vcf_dir / f"{vcf_prefix}{uuid.uuid4()}{ext}"
            result = pool.apply_async(write_file, (virtual_vcf, output_type, filepath))
            results.append(result)
        [r.wait() for r in results]
