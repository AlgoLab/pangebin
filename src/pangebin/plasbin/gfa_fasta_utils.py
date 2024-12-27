"""Manipulating FASTA and GFA files."""

from __future__ import annotations

import gzip
import shutil
from collections import defaultdict
from collections.abc import Iterable
from enum import StrEnum
from pathlib import Path
from typing import IO, Callable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pangebin.input_output as main_io
from pangebin.plasbin.log_errors_utils import check_num_fields, process_exception

ASSEMBLY_PENALTY_TAG = "ap"


class Assembler(StrEnum):
    """Assembler names."""

    UNICYCLER = "unicycler"
    SKESA = "skesa"
    PANGENOME = "pangenome"


class AssemblerTag(StrEnum):
    """Assembler tags."""

    UNICYCLER = "dp"
    PANGENOME = "cv"

    @classmethod
    def from_assembler(cls, assembler: Assembler) -> AssemblerTag | None:
        """Return the assembler tag for a given assembler."""
        match assembler:
            case Assembler.UNICYCLER:
                return cls.UNICYCLER
            case Assembler.PANGENOME:
                return cls.PANGENOME
            case Assembler.SKESA:
                return None


# ==================================================================================== #
#                                GENERIC FILE FUNCTIONS                                #
# ==================================================================================== #
def __open_file_read(filepath: Path) -> IO[str]:
    """Open a file for reading.

    Parameters
    ----------
    filepath : Path
        Path of file to read.

    Returns
    -------
    file object
        File to read.

    """  # REFACTOR use instead main_io
    filepath = Path(filepath)
    if main_io.is_gz_file(filepath):
        return gzip.open(filepath, "rt")
    return filepath.open()


def __open_file_write(filepath: Path) -> IO[str]:
    """Open a file for writing.

    Parameters
    ----------
    filepath : Path
        Path of file to write to.

    Returns
    -------
    file object
        File to write to.

    """  # REFACTOR use instead main_io
    filepath = Path(filepath)
    if main_io.is_gz_file(filepath):
        return gzip.open(filepath, "wt")
    return filepath.open("w")


## Reading FASTA files


def gunzip_fasta(
    input_fasta_gz: Path,
    output_fasta_gunzip: Path,
) -> None:
    """Gunzip a FASTA file.

       in_file_path (str): path to input gzipped FASTA file
       out_file_path (str): path to output FASTA file

    Returns:
       Creates FASTA file out_file_path

    """
    with gzip.open(input_fasta_gz, "rt") as handle:
        records: list[SeqIO.SeqRecord] = list(SeqIO.parse(handle, "fasta"))
        with Path(output_fasta_gunzip).open("w") as out_file:
            SeqIO.write(records, out_file, "fasta")


def read_fasta_contigs(in_file_path, record_fun, id_fun=lambda x: x):
    """Read FASTA file, processing each entry

    Args:
        - in_file_path (str): path of FASTA file to read
        - record_fun (function): processing function taking a single input of type SeqIO
        - gzipped (bool): True if file is gzipped
        - id_fun (function) processing function for record is

    Returns:
        - (Dictionary) sequence id (str) -> record_fun(sequence.record)

    """
    try:
        ctgs_dict = {
            id_fun(seq_id): record_fun(record)
            for seq_id, record in SeqIO.to_dict(
                SeqIO.parse(__open_file_read(in_file_path), "fasta"),
            ).items()
        }
    except Exception as e:
        process_exception(f"Reading {in_file_path}: {e}")
    return ctgs_dict


def read_fasta_id(in_file_path: Path, id_fun=lambda x: x):
    """Compute the list of sequences id in a FASTA file.

    Args:
        - in_file_path (str): path of FASTA file to read
        - id_fun (function) processing function for record is

    Returns:
        - (List(str)): list of sequence ids

    """
    return [
        id_fun(ctg_id)
        for ctg_id in list(
            read_fasta_contigs(
                in_file_path,
                record_fun=lambda _: None,
            ).keys(),
        )
    ]


def read_fasta_lengths(in_file_path, id_fun=lambda x: x):
    """Computes the length of entry sequences in a FASTA file

    Args:
        - in_file_path (str): path of FASTA file to read
        - id_fun (function) processing function for record is

    Returns:
        - (Dictionary): sequence id (str) -> length of sequence (int)

    """
    return read_fasta_contigs(
        in_file_path,
        record_fun=lambda x: len(x.seq),
        id_fun=id_fun,
    )


def read_fasta_sequences(in_file_path, id_fun=lambda x: x):
    """Computes entry sequences in a FASTA file

    Args:
        - in_file_path (str): path of FASTA file to read
        - id_fun (function) processing function for record is

    Returns:
        - (Dictionary): sequence id (str) -> sequence (str)

    """
    return read_fasta_contigs(
        in_file_path,
        record_fun=lambda x: str(x.seq),
        id_fun=id_fun,
    )


## Reading GFA files (contigs and links only)


def gunzip_GFA(in_file_path, out_file_path):
    """Gunzip a GFA file

    Args:
       in_file_path (str): path to input gzipped GFA file
       out_file_path (str): path to output GFA file

    Returns:
       Creates GFA file out_file_path

    """
    try:
        with gzip.open(in_file_path) as in_file, open(out_file_path, "wb") as out_file:
            shutil.copyfileobj(in_file, out_file)
    except Exception as e:
        process_exception(f"FASTA\tGunzipping {in_file_path} to {out_file_path}: {e}")


# Mandatory fields in GFA contigs and links
GFA_SEQ_KEY = "Sequence"
GFA_LEN_KEY = "Length"
GFA_FROM_KEY = "From"
GFA_FROM_ORIENT_KEY = "FromOrient"
GFA_TO_KEY = "To"
GFA_TO_ORIENT_KEY = "ToOrient"
GFA_OVERLAP_KEY = "Overlap"

# Conversionto of GFA attributes.
# Missing attributes types: B, J
GFA_ATTRIBUTE_TYPE = {
    "i": lambda x: int(float(x)),
    "f": lambda x: float(x),
    "Z": lambda x: str(x),
    "A": lambda x: str(x),
    "H": lambda x: bytes(x),
}
# 'i' converted into float at first to handle cases such as '1.20934e+06'


def __add_attributes(attributes_data: Iterable[str], attributes_list: Iterable[str]):
    """Create a dictionary of attributes for a contig/link.

    Args:
        - attributes_data (List): list of attributes in format str(key:type:value)
        - attributes_list (List(str)): list of attribute keys to read
          ['all'] for recording all attributes

    Returns:
        - (Dictionary) attribute key: attribute value (None if missing attribute)
        attributes values are converted to the type defined by GFA_ATTRIBUTE_TYPE
        if an attribute key in attributes_list is not a key in GFA_ATTRIBUTE_TYPE
        the attribute value is recorded as a string

    """
    attributes_dict: dict[str, str | None] = (
        {att_key: None for att_key in attributes_list}
        if attributes_list != ["all"]
        else {}
    )
    for att_data in attributes_data:
        att_split = att_data.split(":")
        check_num_fields(att_split, 3)
        att_key, att_type = att_split[0:2]
        att_val = ":".join(att_split[2:])
        if attributes_list == ["all"] or att_key in attributes_list:
            if att_type not in GFA_ATTRIBUTE_TYPE:
                attributes_dict[att_key] = att_val
            else:
                attributes_dict[att_key] = GFA_ATTRIBUTE_TYPE[att_type](
                    att_val,
                )  # FIXME multiple types: str to be interpred or Any already interpreted
    return attributes_dict


def __write_attributes(attributes_dict, keys_to_remove=[], sep=" "):
    """Write GFA attributes into a string

    Args:
        - attributes_dict (Dictionary): attribute key -> attribute value
        - keys_to_remove (List(str)): list of keys to not print
        - sep (str): separating string

    Returns:
        (str): list of attributes in format sep.join(key:value)
        attributes with None value are not written

    """
    return sep.join(
        [
            f"{x}:{y}"
            for x, y in attributes_dict.items()
            if x not in keys_to_remove and y is not None
        ],
    )
    # REFACTOR rename __write_attributes to gfa_segment_attribute_to_string
    # REFACTOR __write_attributes should be moved


def __assert_attributes_list(attributes_list):
    """Assert that an attributes list is either ['all'] or does not contain 'all'
    Used only for development

    Args:
        - attributes_list (List(str))

    """
    assert (
        attributes_list == ["all"] or "all" not in attributes_list
    ), f"incorrect GFA attributes list {attributes_list}"


def read_gfa_ctgs(
    in_file_path,
    attributes_list,
    ctg_fun=lambda x: x,
    id_fun=lambda x: x,
):
    """Read contigs and their attributes from a GFA files.

    Args:
        - in_file_path (str): path to GFA file to read
        - attributes_list (List(str)): list of attribute keys to read
          ['all'] for recording all attributes
        - ctg_fun: function that process a contig information
        - id_fun: function that process a contig id

    Returns:
       - (Dictionary) contig id -> ctg_fun(
           (Dictionary) attribute key: attribute value
           (None if missing attribute)
         )
         where attribute key GFA_SEQ_KEY is for the contig sequence
    Assumption:
       - every contig has an associated id and sequence (not checked)

    """
    __assert_attributes_list(attributes_list)
    result = {}
    with __open_file_read(in_file_path) as in_file:
        for gfa_line in [x for x in in_file.readlines() if x[0] == "S"]:
            line = gfa_line.rstrip()
            ctg_data = line.split("\t")
            check_num_fields(ctg_data, 2)
            ctg_id, ctg_seq = ctg_data[1], ctg_data[2]
            ctg_len = len(ctg_seq)
            att_data = [
                f"{GFA_SEQ_KEY}:Z:{ctg_seq}",
                f"{GFA_LEN_KEY}:i:{ctg_len}",
            ] + ctg_data[3:]
            result[id_fun(ctg_id)] = ctg_fun(
                __add_attributes(att_data, attributes_list),
            )
    return result


def read_gfa_id(
    gfa_filepath: Path,
    id_fun: Callable[[str], str] = lambda x: x,
) -> list[str]:
    """Compute the list of segments (contigs) id in a GFA file.

    Parameters
    ----------
    gfa_filepath : Path
        Path of GFA file to read.
    id_fun : function, optional
        Function that process a contig id, by default lambda x: x

    Returns
    -------
    list[str]
        List of sequence ids

    """
    return [
        id_fun(ctg_id)
        for ctg_id in list(
            read_gfa_ctgs(
                gfa_filepath,
                attributes_list=[],
                ctg_fun=lambda _: None,
            ).keys(),
        )
    ]


def read_gfa_attribute(
    gfa_filepath: Path,
    att_key: str,
    id_fun=lambda x: x,
):  # TODO new description probably?
    # Computes the length of segments (contigs) in a GFA file ???? # TODO CHECK
    """Return the list of contig:attribute pairs.

    Args:
        - gfa_filepath (str): path of GFA file to read
        - att_key (str): attribute key
        - id_fun: function that process a contig id

    Returns:
        - (Dictionary): sequence id (str) -> attribute value

    """
    return {
        ctg_id: ctg_attributes[att_key]
        for ctg_id, ctg_attributes in read_gfa_ctgs(
            gfa_filepath,
            [att_key],
            id_fun=id_fun,
        ).items()
    }


def read_gfa_len(gfa_filepath: Path, id_fun=lambda x: x):
    """Compute the length of segments (contigs) in a GFA file.

    Args:
        - in_file_path (str): path of GFA file to read
        - id_fun: function that process a contig id

    Returns:
        - (Dictionary): sequence id (str) -> length of sequence (int)

    """
    return read_gfa_attribute(
        gfa_filepath,
        GFA_LEN_KEY,
        id_fun=id_fun,
    )


def read_gfa_seq(gfa_filepath: Path, id_fun=lambda x: x):
    """Compute segments (contigs) sequences in a GFA file.

    Args:
        - in_file_path (str): path of GFA file to read
        - id_fun: function that process a contig id

    Returns:
        - (Dictionary): sequence id (str) -> sequence (str)

    """
    return read_gfa_attribute(
        gfa_filepath,
        GFA_SEQ_KEY,
        id_fun=id_fun,
    )


def _ctgs_normalized_coverage(attributes_dict):
    """Compute the normalized coverage for a set of contigs.

    Args:
        - attributes_dict (Dictionary): contig id -> attributes dictionary, including LN and KC keys
    Returnd:
        (Dictionary): contig id (str) -> normalized coverage (float)

    """
    # Total assembly coverage divided by k (k-mer value)
    total_coverage = sum([ctg_data["KC"] for _, ctg_data in attributes_dict.items()])
    # Total assembly length
    total_length = sum(
        [ctg_data[GFA_LEN_KEY] for _, ctg_data in attributes_dict.items()],
    )
    # Average assembly coverage: k * total_coverage / total_length
    # Contig total coverage = k * ctg.KC
    # Contig average coverage = k * ctg.KC / ctg.Length
    # Contig normalized coverage = (k * ctg.KC / ctg.Length) / (k * total_coverage / total_length)
    # = ctg.KC * total_length / ctg.LN * total_coverage
    return {
        ctg_id: (ctg_data["KC"] * total_length)
        / (ctg_data[GFA_LEN_KEY] * total_coverage)
        for ctg_id, ctg_data in attributes_dict.items()
    }


def read_gfa_ass_penalty(in_file_path, id_fun=lambda x: x):
    """Computes assemby_penalty of segments (fragments) in a Pangenome-GFA file

    Args:
        - in_file_path (str): path of FASTA file to read
        - ass_key (str or None): attribute key that records normalized coverage
           if None, coverage is based on attributes KC and LN
        - id_fun: function that process a contig id

    Returns:
        - (Dictionary): sequence id (str) -> assembly penalty (float) [0.0, 1.0]

    """
    ## se c'e', leggo dal gfa
    ## se non c'e' lo computo? come la normalized coverage

    ## If not tag ap:i:0

    return read_gfa_attribute(
        in_file_path,
        att_key=ASSEMBLY_PENALTY_TAG,
        id_fun=lambda x: x,
    )


def read_gfa_normalized_coverage(
    filepath: Path,
    assembler: Assembler,
    id_fun=lambda x: x,
):
    """Compute normalized coverage of segments (contigs) in a GFA file.

    Args:
        - in_file_path (str): path of FASTA file to read
        - cov_key (str or None): attribute key that records normalized coverage
           if None, coverage is based on attributes KC and LN
        - gzipped (bool): True if file is gzipped
        - id_fun: function that process a contig id

    Returns:
        - (Dictionary): sequence id (str) -> normalized coverage (float)

    """
    # REFACTOR use match instead of conditionnals
    cov_key = AssemblerTag.from_assembler(assembler)
    if cov_key is None:
        ctgs_data = read_gfa_ctgs(
            filepath,
            [GFA_LEN_KEY, "KC"],  # XXX what does it mean?
            id_fun=id_fun,
        )
        try:
            result = _ctgs_normalized_coverage(ctgs_data)
        except Exception as e:
            process_exception(f"GFA\tComputing normalized coverage {filepath}: {e}")
        else:
            return result
    else:
        return read_gfa_attribute(
            filepath,
            str(cov_key),
            id_fun=id_fun,
        )


def read_GFA_links(filepath):
    """Read links and their attributes from a GFA files

    Args:
        - in_file_path (str): path to GFA file to read

    Returns:
       - (Dictionary) contig id ->
         List(links from contig id
           (Dictionary) attribute key: attribute value
           graph attributes: GFA_FROM_ORIENT_KEY, GFA_TO_KEY, GFA_TO_ORIENT_KEY, GFA_OVERLAP_KEY

    """
    result = defaultdict(list)
    with __open_file_read(filepath) as in_file:
        for gfa_line in [x for x in in_file.readlines() if x[0] == "L"]:
            line = gfa_line.rstrip()
            ctg_data = line.split("\t")
            check_num_fields(ctg_data, 5)
            ctg_from = ctg_data[1]
            ctg_from_orient = ctg_data[2]
            ctg_to = ctg_data[3]
            ctg_to_orient = ctg_data[4]
            overlap = ctg_data[5]
            result[ctg_from].append(
                __add_attributes(
                    [
                        f"{GFA_TO_KEY}:Z:{ctg_to}",
                        f"{GFA_FROM_ORIENT_KEY}:A:{ctg_from_orient}",
                        f"{GFA_TO_ORIENT_KEY}:A:{ctg_to_orient}",
                        f"{GFA_OVERLAP_KEY}:Z:{overlap}",
                    ],
                    [
                        GFA_TO_KEY,
                        GFA_FROM_ORIENT_KEY,
                        GFA_TO_ORIENT_KEY,
                        GFA_OVERLAP_KEY,
                    ],
                ),
            )
    return result


def write_gfa_to_fasta(
    in_gfa_filepath: Path,
    out_fasta_filepath: Path,
    sep=" ",
):
    """Create a FASTA file from a GFA file.

    Args:
        - in_GFA_file (str): path to GFA file to read
        - out_FASTA_file (str): path to FASTA file to write
        - sep (str): string for separating GFA attributes

    Returns:
        None, creates file out_FASTA_file
        header format: <contig name> <contig name>.GFA <attributes string=sep.join(key:value)>

    """  # noqa: E501
    gfa_ctg_seqs = read_gfa_ctgs(
        in_gfa_filepath,
        attributes_list=["all"],
    )
    ctg_records = [
        SeqRecord(
            Seq(y[GFA_SEQ_KEY]),
            id=x,
            name=x,
            description=f"{x}.GFA {__write_attributes(y, keys_to_remove=[GFA_SEQ_KEY])}",  # noqa: E501
        )
        for x, y in gfa_ctg_seqs.items()
    ]  # OPTIMIZE write while iterate
    # REFACTOR iter_gfa_to_fasta function
    try:
        with __open_file_write(out_fasta_filepath) as out_file:
            SeqIO.write(ctg_records, out_file, "fasta")
    except Exception as e:
        process_exception(
            f"FASTA/GFA\tWriting {in_gfa_filepath} to {out_fasta_filepath}: {e}",
        )
