# -------------------------------------------------------------------------------------------------
# Copyright (c) 2020, DHS.
#
# This file is part of pytaxonkit (http://github.com/bioforensics/pytaxonkit)
# and is licensed under the BSD license.
#
#
# This Software was prepared for the Department of Homeland Security
# (DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
# part of contract HSHQDC-15-C-00064 to manage and operate the National
# Biodefense Analysis and Countermeasures Center (NBACC), a Federally
# Funded Research and Development Center.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# -------------------------------------------------------------------------------------------------

from builtins import list as pylist
from collections import namedtuple
from io import StringIO
import json
import os
import pandas as pd
from pandas import UInt32Dtype, StringDtype
from pytaxonkit_version import get_versions
import pytest
from subprocess import Popen, PIPE
import sys
from tempfile import NamedTemporaryFile
from warnings import warn

__version__ = get_versions()["version"]
del get_versions


class TaxonKitCLIError(RuntimeError):
    pass


class NCBITaxonomyDumpNotFoundError(FileNotFoundError):
    pass


def _get_taxonkit_version():
    proc = Popen(["taxonkit", "version"], stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise TaxonKitCLIError(err)  # pragma: no cover
    return out.strip()


def log(*args, level="debug"):  # pragma: no cover
    print(f"[pytaxonkit::{level}]", *args, file=sys.stderr)


def validate_data_dir(path):
    filepath = os.path.join(path, "nodes.dmp")
    if not os.path.isfile(filepath) or os.stat(filepath).st_size == 0:
        raise NCBITaxonomyDumpNotFoundError(path)
    return os.path.realpath(path)


def test_validate_data_dir():
    realdir = os.path.realpath(os.path.expanduser("~/.taxonkit/"))
    assert validate_data_dir(os.path.expanduser("~/.taxonkit/")) == realdir
    with pytest.raises(NCBITaxonomyDumpNotFoundError):
        assert validate_data_dir("/path/to/a/non/existent/directory/taxonkit")


def validate_threads(value):
    if value is None:
        return None
    try:
        threadcount = int(value)
        return str(threadcount)
    except ValueError:
        log(f'invalid thread count "{value}"; resetting to taxonkit default', level="warning")
        return None


def test_validate_threads(capsys):
    assert validate_threads(None) is None
    assert validate_threads(2) == "2"
    assert validate_threads("16") == "16"
    out, err = capsys.readouterr()
    assert out == err == ""
    assert validate_threads("StuffedCrust") is None
    out, err = capsys.readouterr()
    assert out == ""
    m = '[pytaxonkit::warning] invalid thread count "StuffedCrust"; resetting to taxonkit default'
    assert err.strip() == m


__taxonkitversion__ = _get_taxonkit_version()


# -------------------------------------------------------------------------------------------------
# taxonkit list
# -------------------------------------------------------------------------------------------------

BasicTaxon = namedtuple("BasicTaxon", ["taxid", "rank", "name"])


class ListResult:
    def __init__(self, jsondata):
        self._data = json.loads(jsondata)

    def __len__(self):
        return len(self._data)

    def __str__(self):
        return json.dumps(self._data, indent=4)

    def __iter__(self):
        for taxonstr, taxtree in self._data.items():
            taxid, rank, name = taxonstr.replace("[", "]").split("]")
            taxon = BasicTaxon(taxid=int(taxid.strip()), rank=rank, name=name.strip())
            if len(taxtree) > 0:
                taxtree = ListResult(json.dumps(taxtree))
            yield taxon, taxtree

    def _do_traverse(self, tree):
        for taxonstr, taxtree in tree.items():
            taxid, rank, name = taxonstr.replace("[", "]").split("]")
            taxon = BasicTaxon(taxid=int(taxid.strip()), rank=rank, name=name.strip())
            yield taxon
            if len(taxtree) > 0:
                for subtaxon in self._do_traverse(taxtree):
                    yield subtaxon

    @property
    def traverse(self):
        for taxon in self._do_traverse(self._data):
            yield taxon


def list(ids, raw=False, threads=None, data_dir=None, debug=False):
    """list taxon tree of given taxids

    Parameters
    ----------
    ids : list or iterable
        A list of taxids (ints or strings are ok)
    raw : bool, default False
        The `ListResult` object returned by default provides convenient programmatic access to the
        result, but it also introduces a significant amount of overhead; for very large results, it
        is recommended to set `raw=True` so that the raw JSON data (loaded into a dictionary) is
        returned
    threads : int, default None
        Override the default taxonkit threads setting
    data_dir : str, default None
        Specify the location of the NCBI taxonomy `.dmp` files; by default, taxonkit searches in
        `~/.taxonkit/`
    debug : bool, default False
        Print debugging output, e.g., system calls to `taxonkit`

    Returns
    -------
    ListResult or dict
        The `pytaxonkit.ListResult` class provides convenient access to the full taxonomic tree
        associated with each top level `taxonkit list` result. When `raw=True`, a dict built from
        the raw JSON is returned.

    Examples
    --------
    >>> import pytaxonkit
    >>> result = pytaxonkit.list([268197, 9903])
    >>> for taxon, tree in result:
    ...     subtaxa = [t for t in tree.traverse]
    ...     print(f'Top level result: {taxon.name} ({taxon.taxid}); {len(subtaxa)} related taxa')
    ...
    Top level result: Polistes comanchus (268197); 2 related taxa
    Top level result: Bos (9903); 33 related taxa
    >>> subtaxa[0]
    BasicTaxon(taxid=9904, rank='species', name='Bos gaurus')
    >>> pytaxonkit.list([9605], raw=True)
    {'9605 [genus] Homo': {'9606 [species] Homo sapiens': {'63221 [subspecies] Homo sapiens neanderthalensis': {}, "741158 [subspecies] Homo sapiens subsp. 'Denisova'": {}}, '1425170 [species] Homo heidelbergensis': {}, '2665952 [no rank] environmental samples': {'2665953 [species] Homo sapiens environmental sample': {}}, '2813598 [no rank] unclassified Homo': {'2813599 [species] Homo sp.': {}}}}
    """  # noqa: E501
    idlist = ",".join(map(str, ids))
    if idlist == "":
        warn("No input for pytaxonkit.list", UserWarning)
        return
    arglist = ["taxonkit", "list", "--json", "--show-name", "--show-rank", "--ids", idlist]
    if threads:
        arglist.extend(("--threads", validate_threads(threads)))
    if data_dir:
        arglist.extend(("--data-dir", validate_data_dir(data_dir)))  # pragma: no cover
    if debug:
        log(*arglist)  # pragma: no cover
    proc = Popen(arglist, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise TaxonKitCLIError(err)  # pragma: no cover
    if raw:
        return json.loads(out)
    else:
        return ListResult(out)


def test_list_leaves(capsys):
    result = list([8204, 2468], debug=True)  # Nota bene: `list` here is `pytaxonkit.list`
    assert len(result) == 2
    top_level_taxa = [taxon for taxon, tree in result]
    sub_trees = [tree for taxon, tree in result]
    assert top_level_taxa == [
        BasicTaxon(taxid=8204, rank="species", name="Anarhichas lupus"),
        BasicTaxon(taxid=2468, rank="species", name="Plasmid NR79"),
    ]
    assert sub_trees == [{}, {}]
    out, err = capsys.readouterr()
    data = "[pytaxonkit::debug] taxonkit list --json --show-name --show-rank --ids 8204,2468"
    assert err.strip() == data


@pytest.mark.parametrize(
    "taxid,taxon,subtaxon,subsubtaxon",
    [
        (
            83882,
            BasicTaxon(taxid=83882, rank="genus", name="Apogon"),
            BasicTaxon(taxid=638272, rank="subgenus", name="Apogon"),
            BasicTaxon(taxid=308069, rank="species", name="Apogon maculatus"),
        ),
        (
            44568,
            BasicTaxon(taxid=44568, rank="genus", name="Parasimulium"),
            BasicTaxon(taxid=61057, rank="subgenus", name="Parasimulium"),
            BasicTaxon(taxid=61060, rank="species", name="Parasimulium crosskeyi"),
        ),
    ],
)
def test_list_genera(taxid, taxon, subtaxon, subsubtaxon):
    result = list([taxid], threads=1)  # Nota bene: `list` here is `pytaxonkit.list`
    tax, tree = next(iter(result))
    subtax, subtree = next(iter(tree))
    subsubtax, subsubtree = next(iter(subtree))
    assert tax == taxon
    assert subtax == subtaxon
    assert subsubtax == subsubtaxon
    assert subsubtree == {}


def test_list_str():
    result = list(["20019"])
    with open("sweetleaf.json", "r") as fh:
        assert str(result) == fh.read().strip()


def test_list_empty():
    with pytest.warns(UserWarning, match="No input for pytaxonkit.list"):
        result = list([])
        assert result is None


# -------------------------------------------------------------------------------------------------
# taxonkit lineage
# -------------------------------------------------------------------------------------------------


def lineage(
    ids,
    formatstr=None,
    threads=None,
    data_dir=None,
    debug=False,
):
    """query lineage of given taxids

    Executes `taxonkit lineage` and `taxonkit reformat` to provide both a full and a "standard"
    lineage for each taxid in `ids`.

    Parameters
    ----------
    ids : list or iterable
        A list of taxids (ints or strings are ok)
    formatstr : str, default None
        See [taxonkit reformat2 documentation](https://bioinf.shenwei.me/taxonkit/usage/#reformat)
        for instructions on setting `formatstr` to override the default standard lineage format
    threads : int
        Override the default taxonkit threads setting
    data_dir : str, default None
        Specify the location of the NCBI taxonomy `.dmp` files; by default, taxonkit searches in
        `~/.taxonkit/`
    debug : bool, default False
        Print debugging output, e.g., system calls to `taxonkit`

    Returns
    -------
    DataFrame
        A two-dimensional data structure. The `Code` column refers to the status code returned by
        `taxonkit lineage`: -1 for queries not found, 0 for deleted taxids (in `delnodes.dmp`), new
        taxids for merged taxids (in `merged.dmp`), and own taxid for queries found in `nodes.dmp`.
        See [taxonkit lineage](https://bioinf.shenwei.me/taxonkit/usage/#lineage) for details.

    Examples
    --------
    >>> import pytaxonkit
    >>> result = pytaxonkit.lineage([1325911, 1649473, 1401311])
    >>> result.columns
    Index(['TaxID', 'Code', 'Name', 'Lineage', 'LineageTaxIDs', 'Rank', 'FullLineage', 'FullLineageTaxIDs', 'FullLineageRanks'], dtype='object')
    >>> result[["TaxID", "Lineage", "LineageTaxIDs"]]
         TaxID                                                                  Lineage                         LineageTaxIDs
    0  1325911      Eukaryota;Arthropoda;Insecta;Hymenoptera;Eucharitidae;Pogonocharis;  2759;6656;50557;7399;216140;1325911;
    1  1649473  Bacteria;Bacteroidota;Cytophagia;Cytophagales;Spirosomataceae;Nibrella;  2;976;768503;768507;2896860;1649473;
    2  1401311         Eukaryota;Arthropoda;Insecta;Coleoptera;Staphylinidae;Styngetus;   2759;6656;50557;7041;29026;1401311;
    >>> result = pytaxonkit.lineage(["1382510", "929505", "390333"], formatstr="{family};{genus};{species};{subspecies|strain}")
    >>> result[["TaxID", "Lineage", "LineageTaxIDs"]]
         TaxID                                                                                               Lineage           LineageTaxIDs
    0  1382510    Enterobacteriaceae;Salmonella;Salmonella bongori;Salmonella bongori serovar 48:z41:-- str. RKS3044   543;590;54736;1382510
    1   929505               Clostridiaceae;Clostridium;Clostridium botulinum;Clostridium botulinum C str. Stockholm  31979;1485;1491;929505
    2   390333  Lactobacillaceae;Lactobacillus;Lactobacillus delbrueckii;Lactobacillus delbrueckii subsp. bulgaricus    33958;1578;1584;1585
    """  # noqa: E501
    idlist = "\n".join(map(str, ids))
    if idlist == "":
        warn("No input for pytaxonkit.lineage", UserWarning)
        return
    arglist = [
        "taxonkit",
        "lineage",
        "--show-lineage-taxids",
        "--show-rank",
        "--show-status-code",
        "--show-name",
        "--show-lineage-ranks",
    ]
    if threads:
        arglist.extend(("--threads", validate_threads(threads)))
    if data_dir:
        arglist.extend(("--data-dir", validate_data_dir(data_dir)))  # pragma: no cover
    if debug:
        log(*arglist)
    with NamedTemporaryFile(suffix="-lineage.txt") as lineagefile:
        proc = Popen(arglist, stdin=PIPE, stdout=lineagefile, stderr=PIPE, universal_newlines=True)
        out, err = proc.communicate(input=idlist)
        if proc.returncode != 0:
            raise TaxonKitCLIError(err)  # pragma: no cover
        lineagefile.flush()
        os.fsync(lineagefile.fileno())
        extraargs = []
        if formatstr:
            extraargs.extend(("--format", formatstr))
        if threads:
            extraargs.extend(("--threads", validate_threads(threads)))
        if data_dir:
            extraargs.extend(("--data-dir", validate_data_dir(data_dir)))  # pragma: no cover
        arglist = [
            "taxonkit",
            "reformat2",
            *extraargs,
            "--taxid-field",
            "1",
            "--show-lineage-taxids",
            lineagefile.name,
        ]
        if debug:
            log(*arglist)
        proc = Popen(arglist, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        out, err = proc.communicate()
        if proc.returncode != 0:
            raise TaxonKitCLIError(err)  # pragma: no cover
        columnorderin = [
            "TaxID",
            "Code",
            "FullLineage",
            "FullLineageTaxIDs",
            "Name",
            "Rank",
            "FullLineageRanks",
            "Lineage",
            "LineageTaxIDs",
        ]
        columnorderout = [
            "TaxID",
            "Code",
            "Name",
            "Lineage",
            "LineageTaxIDs",
            "Rank",
            "FullLineage",
            "FullLineageTaxIDs",
            "FullLineageRanks",
        ]
        data = pd.read_csv(
            StringIO(out), sep="\t", header=None, names=columnorderin, index_col=False
        )
        data = data[columnorderout]
        return data


def name(ids, data_dir=None, debug=False):
    """rapid taxon name retrieval

    Uses the `--no-lineage` option in `taxonkit lineage` for rapid retrieval of taxon names.

    Parameters
    ----------
    ids : list or iterable
        A list of taxids (ints or strings are ok)
    data_dir : str, default None
        Specify the location of the NCBI taxonomy `.dmp` files; by default, taxonkit searches in
        `~/.taxonkit/`
    debug : bool, default False
        Print debugging output, e.g., system calls to `taxonkit`

    Returns
    -------
    DataFrame
        A two-dimensional data structure with TaxIDs and taxon names.

    Examples
    --------
    >>> import pytaxonkit
    >>> name(["151837", "2216222", "517824"])
         TaxID                                Name
    0   151837                    Hiraea smilacina
    1  2216222         Paramyia sp. BIOUG21706-A10
    2   517824  soil bacterium Cipr-S1N-M1LLLSSL-1
    """
    idlist = "\n".join(map(str, ids)) + "\n"
    if idlist == "\n":
        warn("No input for pytaxonkit.name", UserWarning)
        return
    arglist = ["taxonkit", "lineage", "--show-name", "--no-lineage"]
    if data_dir:
        arglist.extend(("--data-dir", validate_data_dir(data_dir)))  # pragma: no cover
    if debug:
        log(*arglist)
    proc = Popen(arglist, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=idlist)
    data = pd.read_csv(
        StringIO(out), sep="\t", header=None, names=["TaxID", "Name"], index_col=False
    )
    return data


def test_lineage(capsys):
    result = lineage(["1082657", "265720", "1191594", "106649", "2868953"], debug=True)
    assert result.TaxID.equals(pd.Series([1082657, 265720, 1191594, 106649, 2868953]))
    assert result.Code.equals(pd.Series([1082657, 265720, 1191594, 106649, 2868953]))
    assert result.Lineage.equals(
        pd.Series(
            [
                "Eukaryota;Discosea;;Longamoebia;Acanthamoebidae;Acanthamoeba;Acanthamoeba sp. TW95",
                "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Porphyromonadaceae;Porphyromonas;Porphyromonas genomosp. P3",
                "Eukaryota;Basidiomycota;Agaricomycetes;Russulales;Russulaceae;Russula;Russula carmesina",
                "Bacteria;Pseudomonadota;Gammaproteobacteria;Moraxellales;Moraxellaceae;Acinetobacter;Acinetobacter guillouiae",
                "Eukaryota;Arthropoda;Insecta;Hemiptera;Lygaeidae;Lygaeosoma;Lygaeosoma sardeum",
            ]
        )
    )
    print(result.LineageTaxIDs.to_list())
    assert result.LineageTaxIDs.equals(
        pd.Series(
            [
                "2759;555280;;1485168;33677;5754;1082657",
                "2;976;200643;171549;171551;836;265720",
                "2759;5204;155619;452342;5401;5402;1191593",
                "2;1224;1236;2887326;468;469;106649",
                "2759;6656;50557;7524;7533;2868952;2868953",
            ]
        )
    )
    assert result.Rank.equals(pd.Series(["species", "species", "varietas", "species", "species"]))

    out, err = capsys.readouterr()
    assert "taxonkit lineage --show-lineage-taxids --show-rank --show-status-code" in err
    assert "taxonkit reformat2 --taxid-field 1 --show-lineage-taxids" in err


def test_lineage_single_taxid():
    result = lineage([128370])
    assert result.TaxID.iloc[0] == 128370


def test_lineage_threads():
    result = lineage(["200643"], threads=1)
    assert (
        result.FullLineageRanks.iloc[0] == "cellular root;domain;kingdom;clade;clade;phylum;class"
    )
    expected = "cellular organisms;Bacteria;Pseudomonadati;FCB group;Bacteroidota/Chlorobiota group;Bacteroidota;Bacteroidia"
    assert result.FullLineage.iloc[0] == expected


def test_lineage_name():
    result = lineage(["526061"])
    assert result.Name.iloc[0] == "Henosepilachna sp. AGBA-2008"


def test_lineage_format():
    formatstr = "k__{domain|acellular root|superkingdom};p__{phylum};c__{class};o__{order};f__{family};g__{genus};s__{species}"
    result = lineage([64191], formatstr=formatstr)
    obs_out = result.Lineage.iloc[0]
    exp_out = "k__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__;f__;g__;s__magnetic proteobacterium strain rj53"
    assert exp_out == obs_out
    formatstr = "k__{domain|acellular root|superkingdom};PHYLUM:{phylum};CLASS:{class};o__{order};f__{family};g__{genus};s__{species}"
    result = lineage(["229933"], formatstr=formatstr)
    obs_out = result.Lineage.iloc[0]
    exp_out = "k__Bacteria;PHYLUM:Pseudomonadota;CLASS:Alphaproteobacteria;o__Rickettsiales;f__Anaplasmataceae;g__Wolbachia;s__Wolbachia endosymbiont of Togo hemipterus (strain 1)"
    assert exp_out == obs_out


def test_name_debug(capsys):
    result = name([207661, 1353792, 1597281], debug=True)
    assert result.Name.equals(
        pd.Series(
            [
                "Ahnfeltiopsis concinna",
                "Picobirnavirus turkey/USA-1512/2010",
                "Isopoda sp. NZAC 03013534",
            ]
        )
    )
    out, err = capsys.readouterr()
    assert "taxonkit lineage --show-name --no-lineage" in err


def test_name_regression():
    result = name([6])
    print(result)
    assert len(result) == 1
    assert result.Name[0] == "Azorhizobium"


def test_lineage_empty():
    with pytest.warns(UserWarning, match="No input for pytaxonkit.lineage"):
        result = lineage([])
        assert result is None


def test_name_empty():
    with pytest.warns(UserWarning, match="No input for pytaxonkit.name"):
        result = name([])
        assert result is None


# -------------------------------------------------------------------------------------------------
# taxonkit name2taxid
# -------------------------------------------------------------------------------------------------


def name2taxid(names, sciname=False, threads=None, data_dir=None, debug=False):
    """query taxid by taxon scientific name

    Parameters
    ----------
    names : list or iterable
        A list of species names or synonyms
    sciname: bool, default False
        By default, both scientific names and synonyms are supported; when `sciname=True`, synonyms
        are ignored
    threads : int
        Override the default taxonkit threads setting
    data_dir : str, default None
        Specify the location of the NCBI taxonomy `.dmp` files; by default, taxonkit searches in
        `~/.taxonkit/`
    debug : bool, default False
        Print debugging output, e.g., system calls to `taxonkit`

    Returns
    -------
    DataFrame
        A two-dimensional data structure.

    Examples
    --------
    >>> import pytaxonkit
    >>> names = ["Phyllobolus spinuliferus", "Alteromonas putrefaciens", "Rexia erectus"]
    >>> pytaxonkit.name2taxid(names)
                           Name   TaxID     Rank
    0  Phyllobolus spinuliferus  359607  species
    1  Alteromonas putrefaciens      24  species
    2             Rexia erectus  262902  species
    >>> pytaxonkit.name2taxid(names, sciname=True)
                           Name  TaxID  Rank
    0  Phyllobolus spinuliferus   <NA>  <NA>
    1  Alteromonas putrefaciens   <NA>  <NA>
    2             Rexia erectus   <NA>  <NA>
    """
    namelist = "\n".join(map(str, names))
    if namelist == "":
        warn("No input for pytaxonkit.name2taxid", UserWarning)
        return
    arglist = ["taxonkit", "name2taxid", "--show-rank"]
    if sciname:
        arglist.append("--sci-name")
    if threads:
        arglist.extend(("--threads", validate_threads(threads)))
    if data_dir:
        arglist.extend(("--data-dir", validate_data_dir(data_dir)))  # pragma: no cover
    if debug:
        log(*arglist)  # pragma: no cover
    proc = Popen(arglist, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=namelist)
    if proc.returncode != 0:
        raise TaxonKitCLIError(err)  # pragma: no cover
    columns = {
        "Name": StringDtype(),
        "TaxID": UInt32Dtype(),
        "Rank": StringDtype(),
    }
    data = pd.read_csv(
        StringIO(out), sep="\t", header=None, names=columns, dtype=columns, index_col=False
    )
    return data


def test_name2taxid(capsys):
    result = name2taxid(["Chaetocerotales", "Diptera", "Rickettsiales", "Hypocreales"], debug=True)
    taxids = pd.Series([265576, 7147, 766, 5125], dtype=UInt32Dtype())
    ranks = pd.Series(["order", "order", "order", "order"], dtype=StringDtype())
    assert result.TaxID.equals(taxids)
    assert result.Rank.equals(ranks)

    out, err = capsys.readouterr()
    assert "[pytaxonkit::debug] taxonkit name2taxid --show-rank" in err


def test_name2taxid_threads():
    result = name2taxid(["FCB group"], threads="1")
    assert str(result) == "        Name    TaxID   Rank\n0  FCB group  1783270  clade"


def test_name2taxid_empty():
    with pytest.warns(UserWarning, match="No input for pytaxonkit.name2taxid"):
        result = name2taxid([])
        assert result is None


# -------------------------------------------------------------------------------------------------
# taxonkit filter
# -------------------------------------------------------------------------------------------------


def filter(
    ids,
    threads=None,
    equal_to=None,
    higher_than=None,
    lower_than=None,
    discard_norank=False,
    save_predictable=False,
    discard_root=False,
    root_taxid=None,
    blacklist=None,
    rank_file=None,
    debug=False,
):
    """filter taxids by taxonomic rank (or a range of ranks)

    Executes the `taxonkit filter` command to include or exclude taxa at the specified ranks.

    Parameters
    ----------
    ids : list or iterable
        A list of taxids (ints or strings are ok)
    threads : int
        Override the default taxonkit threads setting
    equal_to : str or list, default None
        Keep only taxa at the specified rank(s); can be a string or a list of strings
    higher_than : str, default None
        Keep only taxa ranked higher than the specified rank
    lower_than : str, default None
        Keep only taxa ranked lower than the specified rank
    discard_norank : bool, default False
        Discard generic ranks without an explicit ranking order ("no rank" and "clade")
    save_predictable : bool, default False
        When `discard_norank=True`, do not discard some special ranks without order where the rank
        of the closest higher node is still lower than rank cutoff
    discard_root : bool, default False
        Discard root taxon
    root_taxid : int or str
        override taxid of the root taxon
    blacklist : list of strs
        A list of ranks to exclude
    rank_file : str, default None
        Specify the location of the rank definition and order file; by default, taxonkit uses
        `~/taxonkit/ranks.txt`
    debug : bool, default False
        Print debugging output, e.g., system calls to `taxonkit`

    Returns
    -------
    list
        A list of taxids passing the specified filters.

    >>> import pytaxonkit
    >>> taxids = [131567, 2759, 33154, 33208, 6072, 33213, 33317, 1206794, 88770, 6656, 197563, 197562, 6960, 50557, 85512, 7496, 33340, 33392, 85604, 7088]
    >>> result = pytaxonkit.filter(taxids, equal_to='phylum', higher_than='phylum')
    >>> pytaxonkit.name(result)
          TaxID                Name
    0    131567  cellular organisms
    1      2759           Eukaryota
    2     33154        Opisthokonta
    3     33208             Metazoa
    4      6072           Eumetazoa
    5     33213           Bilateria
    6     33317         Protostomia
    7   1206794           Ecdysozoa
    8     88770       Panarthropoda
    9      6656          Arthropoda
    10   197563         Mandibulata
    11   197562        Pancrustacea
    12    85512          Dicondylia
    >>> taxids = [131567, 2759, 33154, 33208, 6072, 33213, 33317, 1206794, 88770, 6656, 197563, 197562, 6960, 50557, 85512, 7496, 33340, 33342, 7524]
    >>> result = pytaxonkit.filter(taxids, lower_than='phylum', discard_norank=True)
    >>> pytaxonkit.name(result)
       TaxID          Name
    0   6960      Hexapoda
    1  50557       Insecta
    2   7496     Pterygota
    3  33340      Neoptera
    4  33342  Paraneoptera
    5   7524     Hemiptera
    """
    if higher_than is not None and lower_than is not None:
        raise ValueError('cannot specify "higher_than" and "lower_than" simultaneously')
    idlist = "\n".join(map(str, ids))
    if idlist == "":
        warn("No input for pytaxonkit.filter", UserWarning)
        return
    arglist = ["taxonkit", "filter"]
    if threads:
        arglist.extend(("--threads", validate_threads(threads)))
    if equal_to:
        if isinstance(equal_to, (pylist, tuple)):
            equal_to = ",".join(equal_to)
        arglist.extend(["--equal-to", equal_to])
    if higher_than:
        arglist.extend(["--higher-than", higher_than])
    if lower_than:
        arglist.extend(["--lower-than", lower_than])
    if discard_norank:
        arglist.append("--discard-noranks")
    if save_predictable:
        arglist.append("--save-predictable-norank")
    if discard_root:  # pragma: no cover
        arglist.append("--discard-root")
    if blacklist:
        arglist.extend(["--black-list", ",".join(blacklist)])
    if root_taxid:  # pragma: no cover
        arglist.extend(["--root-taxid", str(root_taxid)])
    if rank_file:  # pragma: no cover
        arglist.extend(["--rank-file", rank_file])
    if debug:
        log(*arglist)
    proc = Popen(arglist, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=idlist)
    data = pd.read_csv(StringIO(out), header=None, names=["TaxID"], index_col=False)
    return pylist(data.TaxID)


def list_ranks(rank_file=None, debug=False):
    """list all supported taxonomic ranks

    Parameters
    ----------
    rank_file : str, default None
        Specify the location of the rank definition and order file; by default, taxonkit uses
        `~/.taxonkit/ranks.txt`
    debug : bool, default False
        Print debugging output, e.g., system calls to `taxonkit`

    Returns
    -------
    list
        A list of taxonomic ranks.

    >>> import pytaxonkit
    >>> ranks = pytaxonkit.list_ranks()
    >>> ranks[:5]
    ['life', 'acellular', ['root', 'cellular'], 'root', ['domain', 'empire', 'realm', 'superkingdom']]
    """  # noqa: E501
    arglist = ["taxonkit", "filter", "--list-order"]
    if rank_file:  # pragma: no cover
        arglist.extend(["--rank-file", rank_file])
    if debug:
        log(*arglist)
    proc = Popen(arglist, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input="")
    ranks = pylist()
    for line in out.strip().split():
        rankvalue = line.split(",") if "," in line else line
        ranks.append(rankvalue)
    return ranks


def list_ranks_db(rank_file=None, debug=False):
    """list all taxonomic ranks present in the NCBI taxonomy database

    Parameters
    ----------
    rank_file : str, default None
        Specify the location of the rank definition and order file; by default, taxonkit uses
        `~/taxonkit/ranks.txt`
    debug : bool, default False
        Print debugging output, e.g., system calls to `taxonkit`

    Returns
    -------
    list
        A list of taxonomic ranks.

    >>> import pytaxonkit
    >>> ranks = pytaxonkit.list_ranks_db()
    >>> ranks[:5]
    ['acellular root', 'cellular root', 'domain', 'realm', 'kingdom']
    """
    arglist = ["taxonkit", "filter", "--list-ranks"]
    if rank_file:  # pragma: no cover
        arglist.extend(["--rank-file", rank_file])
    if debug:
        log(*arglist)
    proc = Popen(arglist, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input="")
    data = pd.read_csv(StringIO(out), header=None, names=["Rank"], index_col=False)
    return pylist(data.Rank)


@pytest.mark.parametrize(
    "discard_norank,exp_result",
    [
        (True, [131567, 2759, 33208, 6656, 6960, 50557, 7496, 33340, 33392, 7399]),
        (
            False,
            [
                131567,
                2759,
                33154,
                33208,
                6072,
                33213,
                33317,
                1206794,
                88770,
                6656,
                197563,
                197562,
                6960,
                50557,
                85512,
                7496,
                33340,
                33392,
                7399,
            ],
        ),
    ],
)
def test_filter_higher_than(discard_norank, exp_result):
    taxids = [
        131567,
        2759,
        33154,
        33208,
        6072,
        33213,
        33317,
        1206794,
        88770,
        6656,
        197563,
        197562,
        6960,
        50557,
        85512,
        7496,
        33340,
        33392,
        7399,
        7400,
        7434,
        34735,
        7458,
        70987,
        83322,
        481579,
        2056706,
        599582,
    ]
    obs_result = filter(
        taxids, threads=1, higher_than="order", equal_to="order", discard_norank=discard_norank
    )
    print(obs_result)
    assert obs_result == exp_result


def test_filter_lower_than(capsys):
    taxids = [
        131567,
        2759,
        33154,
        33208,
        6072,
        33213,
        33317,
        1206794,
        88770,
        6656,
        197563,
        197562,
        6960,
        50557,
        85512,
        7496,
        33340,
        33392,
        7041,
        41071,
        535382,
        41073,
        706613,
        586004,
        87479,
        412111,
    ]
    obs_result = filter(taxids, lower_than="family", discard_norank=True, debug=True)
    exp_result = [706613, 586004, 87479, 412111]
    assert obs_result == exp_result
    terminal = capsys.readouterr()
    assert "taxonkit filter --lower-than family" in terminal.err


def test_filter_equal_to_multi(capsys):
    taxids = [
        131567,
        2759,
        33154,
        33208,
        6072,
        33213,
        33317,
        1206794,
        88770,
        6656,
        197563,
        197562,
        6960,
        50557,
        85512,
        7496,
        33340,
        33392,
        7041,
        41071,
        535382,
        41073,
        706613,
        586004,
        87479,
        412111,
    ]
    obs_result = filter(taxids, threads=1, equal_to=["phylum", "genus"], debug=True)
    assert obs_result == [6656, 87479]
    terminal = capsys.readouterr()
    assert "--equal-to phylum,genus" in terminal.err


def test_filter_higher_lower_conflict():
    message = r'cannot specify "higher_than" and "lower_than" simultaneously'
    with pytest.raises(ValueError, match=message):
        filter([42], discard_norank=True, higher_than="genus", lower_than="genus")


def test_filter_save_predictable():
    taxids = [
        131567,
        2,
        1224,
        1236,
        91347,
        543,
        561,
        562,
        2605619,
        10239,
        2731341,
        2731360,
        2731618,
        2731619,
        2788787,
        1327037,
    ]
    obs_result = filter(
        taxids, threads=1, lower_than="species", equal_to="species", save_predictable=True
    )
    exp_result = [562, 2605619, 1327037]
    assert obs_result == exp_result


def test_list_ranks(capsys):
    ranks = list_ranks(debug=True)
    multiranks = [r for r in ranks if isinstance(r, pylist)]
    assert len(ranks) == 74
    assert len(multiranks) == 18
    terminal = capsys.readouterr()
    assert "taxonkit filter --list-order" in terminal.err


def test_list_ranks_db(capsys):
    ranks = list_ranks_db(debug=True)
    assert len(ranks) == 48
    terminal = capsys.readouterr()
    assert "taxonkit filter --list-ranks" in terminal.err


def test_filter_empty():
    with pytest.warns(UserWarning, match="No input for pytaxonkit.filter"):
        result = filter([])
        assert result is None


# -------------------------------------------------------------------------------------------------
# taxonkit lca
# -------------------------------------------------------------------------------------------------


def lca(
    ids,
    multi=False,
    skip_deleted=False,
    skip_unfound=False,
    keep_invalid=False,
    threads=None,
    data_dir=None,
    debug=False,
):
    """compute lowest common ancestor (LCA) for taxids

    Parameters
    ----------
    ids : list or iterable
        A list of taxids (ints or strings are ok); if `multi=True`, a list of lists is expected
    multi : bool, default False
        The input consists of multiple LCA queries, stored as a list of taxid lists; by default,
        a single list corresponding to a single LCA query is expected
    skip_deleted : bool, default False
        Ignore deleted taxids and compute LCA with the remaining taxa
    skip_unfound : bool, default False
        Ignore taxids not found in the taxonomy database and compute LCA with the remaining taxa
    keep_invalid: bool, default False
        Returns 0 when all taxids have been skipped from `skip_deleted` or `skip_unfound`
    threads : int
        Override the default taxonkit threads setting
    data_dir : str, default None
        Specify the location of the NCBI taxonomy `.dmp` files; by default, taxonkit searches in
        `~/.taxonkit/`
    debug : bool, default False
        Print debugging output, e.g., system calls to `taxonkit`

    Returns
    -------
    int or list
        A taxid or, if `multi=True`, a list of taxids

    Examples
    --------
    >>> import pytaxonkit
    >>> taxids = pytaxonkit.name2taxid(['Polistes metricus', 'Nasonia vitripennis'])
    >>> taxids
                      Name  TaxID     Rank
    0    Polistes metricus  91422  species
    1  Nasonia vitripennis   7425  species
    >>> ancestor = pytaxonkit.lca(taxids.TaxID)
    >>> pytaxonkit.name([ancestor])
       TaxID      Name
    0   7400  Apocrita
    >>> pytaxonkit.lca([239934, 239935, 349741])
    239934
    >>> pytaxonkit.lca([[63221, 2665953], [63221, 741158]], multi=True)
    [9605, 9606]
    """
    if multi:
        idstring = "\n".join([" ".join(map(str, sublist)) for sublist in ids])
    else:
        idstring = " ".join(map(str, ids))
    if idstring == "":
        warn("No input for pytaxonkit.lca", UserWarning)
        return
    arglist = ["taxonkit", "lca"]
    if skip_deleted:
        arglist.append("--skip-deleted")
    if skip_unfound:
        arglist.append("--skip-unfound")
    if keep_invalid:
        arglist.append("--keep-invalid")
    if threads:
        arglist.extend(("--threads", validate_threads(threads)))
    if data_dir:
        arglist.extend(("--data-dir", validate_data_dir(data_dir)))  # pragma: no cover
    if debug:
        log(*arglist)
    proc = Popen(arglist, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=idstring)
    if proc.returncode != 0:
        raise TaxonKitCLIError(err)  # pragma: no cover
    if out.strip() == "":
        if multi:
            return [None] * len(ids)
        else:
            return None
    results = []
    for line in out.strip().split("\n"):
        queries, lcataxid = line.split("\t")
        results.append(int(lcataxid))
    if multi:
        return results
    else:
        assert len(results) == 1
        return results[0]


def test_lca_deleted():
    assert lca([1, 2, 3]) == 0
    assert lca([1, 2, 3], skip_deleted=True, threads=1) == 1


def test_lca_unfound(capsys):
    assert lca([61021, 61022, 11111111]) == 0
    assert lca([61021, 61022, 11111111], skip_unfound=True, debug=True) == 2628496
    terminal = capsys.readouterr()
    assert "taxonkit lca --skip-unfound" in terminal.err


def test_lca_keep_invalid_single(capsys):
    assert lca([11111111], skip_deleted=True, skip_unfound=True) is None
    assert lca([22222222], skip_deleted=True, skip_unfound=True) is None
    assert lca([11111111], skip_deleted=True, skip_unfound=True, keep_invalid=True) == 0
    result = lca(
        [11111111, 22222222],
        skip_deleted=True,
        skip_unfound=True,
        keep_invalid=True,
        debug=True,
    )
    assert result == 0
    terminal = capsys.readouterr()
    assert "taxonkit lca --skip-deleted --skip-unfound --keep-invalid" in terminal.err


def test_lca_keep_invalid_multi():
    query = [
        [743375],
        [123456789],
        [987654321],
        [743375, 123456789],
        [743375, 987654321],
        [123456789, 987654321],
    ]
    observed = lca(query, skip_deleted=True, skip_unfound=True, keep_invalid=True, multi=True)
    expected = [743375, 0, 0, 743375, 743375, 0]
    assert expected == observed


@pytest.mark.parametrize(
    "domulti, ids,result",
    [
        (False, [775536, 2238728, 1121123211234321], None),
        (True, [[1766280, 406491, 2568889], [11111111111, 20487, 760325]], [None, None]),
    ],
)
def test_lca_all_missing(ids, domulti, result):
    assert lca(ids, multi=domulti, skip_deleted=True, skip_unfound=True) == result


@pytest.mark.parametrize("multi", [True, False])
def test_lca_empty(multi):
    with pytest.warns(UserWarning, match="No input for pytaxonkit.lca"):
        result = lca([], multi=multi)
        assert result is None
