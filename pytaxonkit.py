# -------------------------------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
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

from collections import namedtuple
from io import StringIO
import json
import os
import pandas
from pandas import UInt32Dtype, StringDtype
import pytest
from subprocess import Popen, PIPE
import sys
from tempfile import NamedTemporaryFile

from _version import get_versions
__version__ = get_versions()['version']
del get_versions


class TaxonKitCLIError(RuntimeError):
    pass


class NCBITaxonomyDumpNotFoundError(FileNotFoundError):
    pass


def _get_taxonkit_version():
    proc = Popen(['taxonkit', 'version'], stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise TaxonKitCLIError(err)  # pragma: no cover
    return out.strip()


def log(*args, level='debug'):  # pragma: no cover
    print(f'[pytaxonkit::{level}]', *args, file=sys.stderr)


def validate_data_dir(path):
    filepath = os.path.join(path, 'nodes.dmp')
    if not os.path.isfile(filepath) or os.stat(filepath).st_size == 0:
        raise NCBITaxonomyDumpNotFoundError(path)
    return os.path.realpath(path)


def test_validate_data_dir():
    realdir = os.path.realpath(os.path.expanduser('~/.taxonkit/'))
    assert validate_data_dir(os.path.expanduser('~/.taxonkit/')) == realdir
    with pytest.raises(NCBITaxonomyDumpNotFoundError):
        assert validate_data_dir('/path/to/a/non/existent/directory/taxonkit')


def validate_threads(value):
    if value is None:
        return None
    try:
        threadcount = int(value)
        return str(threadcount)
    except ValueError:
        log(f'invalid thread count "{value}"; resetting to taxonkit default', level='warning')
        return None


def test_validate_threads(capsys):
    assert validate_threads(None) is None
    assert validate_threads(2) == '2'
    assert validate_threads('16') == '16'
    out, err = capsys.readouterr()
    assert out == err == ''
    assert validate_threads('StuffedCrust') is None
    out, err = capsys.readouterr()
    assert out == ''
    m = '[pytaxonkit::warning] invalid thread count "StuffedCrust"; resetting to taxonkit default'
    assert err.strip() == m


__taxonkitversion__ = _get_taxonkit_version()


# -------------------------------------------------------------------------------------------------
# taxonkit list
# -------------------------------------------------------------------------------------------------

BasicTaxon = namedtuple('Taxon', ['taxid', 'rank', 'name'])


class ListResult():
    def __init__(self, jsondata):
        self._data = json.loads(jsondata)

    def __len__(self):
        return len(self._data)

    def __str__(self):
        return json.dumps(self._data, indent=4)

    def __iter__(self):
        for taxonstr, taxtree in self._data.items():
            taxid, rank, name = taxonstr.replace('[', ']').split(']')
            taxon = BasicTaxon(taxid=int(taxid.strip()), rank=rank, name=name.strip())
            if len(taxtree) > 0:
                taxtree = ListResult(json.dumps(taxtree))
            yield taxon, taxtree

    def _do_traverse(self, tree):
        for taxonstr, taxtree in tree.items():
            taxid, rank, name = taxonstr.replace('[', ']').split(']')
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
    '''list taxon tree of given taxids

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
    >>> result = pytaxonkit.list([13685, 9903])
    >>> for taxon, tree in result:
    ...     subtaxa = [t for t in tree.traverse]
    ...     print(f'Top level result: {taxon.name} ({taxon.taxid}); {len(subtaxa)} related taxa')
    ...
    Top level result: Solenopsis (13685); 198 related taxa
    Top level result: Bos (9903); 26 related taxa
    >>> subtaxa[0]
    Taxon(taxid=9904, rank='species', name='Bos gaurus')
    >>> pytaxonkit.list([9605], raw=True)
    {'9605 [genus] Homo': {'9606 [species] Homo sapiens': {'63221 [subspecies] Homo sapiens neanderthalensis': {}, "741158 [subspecies] Homo sapiens subsp. 'Denisova'": {}, '2665952 [no rank] environmental samples': {'2665953 [species] Homo sapiens environmental sample': {}}}, '1425170 [species] Homo heidelbergensis': {}}}
    '''  # noqa: E501
    idlist = ','.join(map(str, ids))
    arglist = ['taxonkit', 'list', '--json', '--show-name', '--show-rank', '--ids', idlist]
    if threads:
        arglist.extend(('--threads', validate_threads(threads)))
    if data_dir:
        arglist.extend(('--data-dir', validate_data_dir(data_dir)))  # pragma: no cover
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
        BasicTaxon(taxid=8204, rank='species', name='Anarhichas lupus'),
        BasicTaxon(taxid=2468, rank='species', name='Plasmid NR79')
    ]
    assert sub_trees == [{}, {}]
    out, err = capsys.readouterr()
    data = '[pytaxonkit::debug] taxonkit list --json --show-name --show-rank --ids 8204,2468'
    assert err.strip() == data


@pytest.mark.parametrize('taxid,taxon,subtaxon,subsubtaxon', [
    (
        83882,
        BasicTaxon(taxid=83882, rank='genus', name='Apogon'),
        BasicTaxon(taxid=638272, rank='subgenus', name='Apogon'),
        BasicTaxon(taxid=308069, rank='species', name='Apogon maculatus'),
    ),
    (
        44568,
        BasicTaxon(taxid=44568, rank='genus', name='Parasimulium'),
        BasicTaxon(taxid=61057, rank='subgenus', name='Parasimulium'),
        BasicTaxon(taxid=61060, rank='species', name='Parasimulium crosskeyi'),
    ),
])
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
    result = list(['20019'])
    with open('sweetleaf.json', 'r') as fh:
        assert str(result) == fh.read().strip()


# -------------------------------------------------------------------------------------------------
# taxonkit lineage
# -------------------------------------------------------------------------------------------------


def lineage(ids, formatstr=None, threads=None, data_dir=None, debug=False):
    '''query lineage of given taxids

    Executes `taxonkit lineage` and `taxonkit reformat` to provide both a full and a "standard"
    lineage for each taxid in `ids`.

    Parameters
    ----------
    ids : list or iterable
        A list of taxids (ints or strings are ok)
    formatstr : str, default None
        See [taxonkit reformat documentation](https://bioinf.shenwei.me/taxonkit/usage/#reformat)
        for instructions on setting `formatstr` to override the default "standard" lineage format
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
    >>> result = pytaxonkit.lineage([7399, 1973489])
    >>> result.columns
    Index(['TaxID', 'Code', 'Lineage', 'LineageTaxIDs', 'Rank', 'FullLineage', 'FullLineageTaxIDs'], dtype='object')
    >>> result[['TaxID', 'Lineage', 'LineageTaxIDs']]
         TaxID                                                                             Lineage                          LineageTaxIDs
    0     7399                                         Eukaryota;Arthropoda;Insecta;Hymenoptera;;;                2759;6656;50557;7399;;;
    1  1973489  Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus sp. ISSFR-25F  2;1239;91061;1385;186817;1386;1973489
    >>> result = pytaxonkit.lineage(['1382510', '929505', '390333'], formatstr='{f};{g};{s};{S}')
    >>> result[['TaxID', 'Lineage', 'LineageTaxIDs']]
         TaxID                                                                                               Lineage         LineageTaxIDs
    0  1382510                                                     Enterobacteriaceae;Salmonella;Salmonella bongori;        543;590;54736;
    1   929505                                                     Clostridiaceae;Clostridium;Clostridium botulinum;      31979;1485;1491;
    2   390333  Lactobacillaceae;Lactobacillus;Lactobacillus delbrueckii;Lactobacillus delbrueckii subsp. bulgaricus  33958;1578;1584;1585
    '''  # noqa: E501
    idlist = '\n'.join(map(str, ids))
    arglist = ['taxonkit', 'lineage', '--show-lineage-taxids', '--show-rank', '--show-status-code']
    if threads:
        arglist.extend(('--threads', validate_threads(threads)))
    if data_dir:
        arglist.extend(('--data-dir', validate_data_dir(data_dir)))  # pragma: no cover
    if debug:
        log(*arglist)  # pragma: no cover
    with NamedTemporaryFile(suffix='-lineage.txt') as lineagefile:
        proc = Popen(arglist, stdin=PIPE, stdout=lineagefile, stderr=PIPE, universal_newlines=True)
        out, err = proc.communicate(input=idlist)
        if proc.returncode != 0:
            raise TaxonKitCLIError(err)  # pragma: no cover
        lineagefile.flush()
        os.fsync(lineagefile.fileno())
        extraargs = []
        if formatstr:
            extraargs.extend(('--format', formatstr))
        if threads:
            extraargs.extend(('--threads', validate_threads(threads)))
        if data_dir:
            extraargs.extend(('--data-dir', validate_data_dir(data_dir)))  # pragma: no cover
        arglist = [
            'taxonkit', 'reformat', *extraargs, '--lineage-field', '3', '--show-lineage-taxids',
            lineagefile.name
        ]
        if debug:
            log(*arglist)
        proc = Popen(arglist, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        out, err = proc.communicate()
        if proc.returncode != 0:
            raise TaxonKitCLIError(err)  # pragma: no cover
        columnorderin = [
            'TaxID', 'Code', 'FullLineage', 'FullLineageTaxIDs', 'Rank', 'Lineage', 'LineageTaxIDs'
        ]
        columnorderout = [
            'TaxID', 'Code', 'Lineage', 'LineageTaxIDs', 'Rank', 'FullLineage', 'FullLineageTaxIDs'
        ]
        data = pandas.read_csv(
            StringIO(out), sep='\t', header=None, names=columnorderin, index_col=False
        )
        data = data[columnorderout]
        return data


def test_lineage(capsys):
    result = lineage(['446045', '265720', '2507530', '106649'], debug=True)
    assert result.TaxID.equals(pandas.Series([446045, 265720, 2507530, 106649]))
    assert result.Code.equals(pandas.Series([446045, 265720, 2507530, 106649]))
    assert result.Lineage.equals(pandas.Series([
        'Eukaryota;Arthropoda;Insecta;Diptera;Drosophilidae;Drosophila;',
        'Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;Porphyromonas;'
        'Porphyromonas genomosp. P3',
        'Eukaryota;Basidiomycota;Agaricomycetes;Russulales;Russulaceae;Russula;Russula species',
        'Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Moraxellaceae;Acinetobacter;'
        'Acinetobacter guillouiae',
    ]))
    assert result.LineageTaxIDs.equals(pandas.Series([
        '2759;6656;50557;7147;7214;7215;',
        '2;976;200643;171549;171551;836;265720',
        '2759;5204;155619;452342;5401;5402;2507520',
        '2;1224;1236;72274;468;469;106649',
    ]))
    assert result.Rank.equals(pandas.Series([
        'no rank', 'species', 'subspecies', 'species'
    ]))

    out, err = capsys.readouterr()
    assert 'taxonkit lineage --show-lineage-taxids --show-rank --show-status-code' in err
    assert 'taxonkit reformat --lineage-field 3 --show-lineage-taxids' in err


def test_lineage_threads():
    result = lineage(['200643'], threads=1)
    assert result.FullLineage.iloc[0] == (
        'cellular organisms;Bacteria;FCB group;Bacteroidetes/Chlorobi group;Bacteroidetes;'
        'Bacteroidia'
    )


# -------------------------------------------------------------------------------------------------
# taxonkit name2taxid
# -------------------------------------------------------------------------------------------------

def name2taxid(names, sciname=False, threads=None, data_dir=None, debug=False):
    '''query taxid by taxon scientific name

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
    >>> names = ['Phyllobolus spinuliferus', 'Alteromonas putrefaciens', 'Rexia erectus']
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
    '''
    namelist = '\n'.join(map(str, names))
    arglist = ['taxonkit', 'name2taxid', '--show-rank']
    if sciname:
        arglist.append('--sci-name')
    if threads:
        arglist.extend(('--threads', validate_threads(threads)))
    if data_dir:
        arglist.extend(('--data-dir', validate_data_dir(data_dir)))  # pragma: no cover
    if debug:
        log(*arglist)  # pragma: no cover
    proc = Popen(arglist, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=namelist)
    if proc.returncode != 0:
        raise TaxonKitCLIError(err)  # pragma: no cover
    columns = {
        'Name': StringDtype(),
        'TaxID': UInt32Dtype(),
        'Rank': StringDtype(),
    }
    data = pandas.read_csv(
        StringIO(out), sep='\t', header=None, names=columns, dtype=columns, index_col=False
    )
    return data


def test_name2taxid(capsys):
    result = name2taxid(['Chaetocerotales', 'Diptera', 'Rickettsiales', 'Hypocreales'], debug=True)
    taxids = pandas.Series([265576, 7147, 766, 5125], dtype=UInt32Dtype())
    ranks = pandas.Series(['order', 'order', 'order', 'order'], dtype=StringDtype())
    assert result.TaxID.equals(taxids)
    assert result.Rank.equals(ranks)

    out, err = capsys.readouterr()
    assert '[pytaxonkit::debug] taxonkit name2taxid --show-rank' in err


def test_name2taxid_threads():
    result = name2taxid(['FCB group'], threads='1')
    assert str(result) == '        Name    TaxID     Rank\n0  FCB group  1783270  no rank'
