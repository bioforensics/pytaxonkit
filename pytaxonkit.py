#!/usr/bin/env python
#
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

from collections import namedtuple, Mapping
from io import StringIO
import json
import os
import pandas
from subprocess import Popen, PIPE
import sys
from tempfile import NamedTemporaryFile

from _version import get_versions
__version__ = get_versions()['version']
del get_versions


def log(*args, level='debug'):  # pragma: no cover
    print(f'[pytaxonkit::{level}]', *args, file=sys.stderr)


BasicTaxon = namedtuple('Taxon', ['taxid', 'rank', 'name'])


class TaxonKitCLIError(RuntimeError):
    pass


# -------------------------------------------------------------------------------------------------
# taxonkit list
# -------------------------------------------------------------------------------------------------

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


def list(ids, raw=False, debug=False):
    '''list taxon tree of given taxids

    Returns a `pytaxonkit.ListResult` object by default. NOTE: while this data structure provides
    convenient programmatic access to the result, it introduces a significant amount of overhead.
    For very large results, it is recommended to set `raw=True` so that the raw JSON data (loaded
    into a dictionary) is returned.

    >>> import pytaxonkit
    >>> result = pytaxonkit.list([8204, 2468])
    >>> for taxon, tree in result:
    ...   print('[Result]', taxon.name, taxon.taxid, taxon.rank, sep='\t')
    ...
    [Result]    Anarhichas lupus    8204    species
    [Result]    Plasmid NR79        2468    species
    >>>
    >>>
    >>> result = pytaxonkit.list([13685, 9903])
    >>> len(result)
    2
    >>> type(result)
    <class 'pytaxonkit.ListResult'>
    >>> for taxon, tree in result:
    ...     print(f'[Top level result] {taxon.name} ({taxon.taxid})')
    ...     subtaxa = [t for t in tree.traverse]
    ...     print(f'    - {len(subtaxa)} sub taxa; first 5:')
    ...     for subtaxon in subtaxa[:5]:
    ...         print('        -', subtaxon)
    ...
    [Top level result] Solenopsis (13685)
        - 198 sub taxa; first 5:
            - Taxon(taxid=13686, rank='species', name='Solenopsis invicta')
            - Taxon(taxid=30203, rank='species', name='Solenopsis richteri')
            - Taxon(taxid=121131, rank='species', name='Solenopsis geminata')
            - Taxon(taxid=176590, rank='species', name='Solenopsis amblychila')
            - Taxon(taxid=176591, rank='species', name='Solenopsis aurea')
    [Top level result] Bos (9903)
        - 26 sub taxa; first 5:
            - Taxon(taxid=9904, rank='species', name='Bos gaurus')
            - Taxon(taxid=1383418, rank='subspecies', name='Bos gaurus gaurus')
            - Taxon(taxid=9906, rank='species', name='Bos javanicus')
            - Taxon(taxid=380177, rank='subspecies', name='Bos javanicus birmanicus')
            - Taxon(taxid=659500, rank='subspecies', name='Bos javanicus javanicus')
    >>>
    >>>
    >>> result = pytaxonkit.list([121498])
    >>> print(result)  # string representation of a `pytaxonkit.ListResult` is the raw JSON
    {
        "121498 [genus] Bothriomyrmex": {
            "121499 [species] Bothriomyrmex meridionalis": {},
            "369102 [species] Bothriomyrmex hispanicus": {},
            "604499 [species] Bothriomyrmex paradoxus": {},
            "604506 [species] Bothriomyrmex saundersi": {},
            "2625763 [no rank] unclassified Bothriomyrmex": {
                "1297180 [species] Bothriomyrmex sp. MAS001": {}
            }
        }
    }
    >>>
    >>>
    >>> result = pytaxonkit.list([9605], raw=True)  # setting `raw=True` returns only the raw JSON
    >>> print(result)
    {'9605 [genus] Homo': {'9606 [species] Homo sapiens': {'63221 [subspecies] Homo sapiens neanderthalensis': {}, "741158 [subspecies] Homo sapiens subsp. 'Denisova'": {}, '2665952 [no rank] environmental samples': {'2665953 [species] Homo sapiens environmental sample': {}}}, '1425170 [species] Homo heidelbergensis': {}}}
    >>>
    '''  # noqa: E501
    idlist = ','.join(map(str, ids))
    arglist = ['taxonkit', 'list', '--json', '--show-name', '--show-rank', '--ids', idlist]
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


def lineage(ids, formatstr=None, debug=False):
    '''query lineage of given taxids

    Runs `taxonkit lineage` and `taxonkit reformat` to provide both a full and a "standard" lineage
    for each taxid provided; by default `formatstr=None` and the standard format string is used
    with `taxonkit reformat`; see the [taxonkit reformat documentation]
    (https://bioinf.shenwei.me/taxonkit/usage/#reformat) for instructions on overriding the
    default.

    >>> import pytaxonkit
    >>> result = pytaxonkit.lineage([7399, 1973489])
    >>> type(result)
    <class 'pandas.core.frame.DataFrame'>
    >>> result
         TaxID     Code                                                                             Lineage                          LineageTaxIDs     Rank                                                                                                                                                                                                            FullLineage                                                                                                 FullLineageTaxIDs
    0     7399     7399                                         Eukaryota;Arthropoda;Insecta;Hymenoptera;;;                2759;6656;50557;7399;;;    order  cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;Ecdysozoa;Panarthropoda;Arthropoda;Mandibulata;Pancrustacea;Hexapoda;Insecta;Dicondylia;Pterygota;Neoptera;Holometabola;Hymenoptera  131567;2759;33154;33208;6072;33213;33317;1206794;88770;6656;197563;197562;6960;50557;85512;7496;33340;33392;7399
    1  1973489  1973489  Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus sp. ISSFR-25F  2;1239;91061;1385;186817;1386;1973489  species                                                                        cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus cereus group;Bacillus sp. ISSFR-25F                                                        131567;2;1783272;1239;91061;1385;186817;1386;86661;1973489
    >>>
    >>>
    >>> result = pytaxonkit.lineage(['1382510', '929505', '390333'], formatstr='{f};{g};{s};{S}', debug=True)
    >>> result
         TaxID     Code                                                                                               Lineage         LineageTaxIDs     Rank                                                                                                                                                                                                                                                FullLineage                                                 FullLineageTaxIDs
    0  1382510  1382510                                                     Enterobacteriaceae;Salmonella;Salmonella bongori;        543;590;54736;  no rank                                    cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella bongori;Salmonella bongori serovar 48:z41:--;Salmonella bongori serovar 48:z41:-- str. RKS3044              131567;2;1224;1236;91347;543;590;54736;41527;1382510
    1   929505   929505                                                     Clostridiaceae;Clostridium;Clostridium botulinum;      31979;1485;1491;  no rank                                                        cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium;Clostridium botulinum;Clostridium botulinum C;Clostridium botulinum C str. Stockholm  131567;2;1783272;1239;186801;186802;31979;1485;1491;36828;929505
    2   390333   390333  Lactobacillaceae;Lactobacillus;Lactobacillus delbrueckii;Lactobacillus delbrueckii subsp. bulgaricus  33958;1578;1584;1585  no rank  cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus delbrueckii;Lactobacillus delbrueckii subsp. bulgaricus;Lactobacillus delbrueckii subsp. bulgaricus ATCC 11842 = JCM 1002    131567;2;1783272;1239;91061;186826;33958;1578;1584;1585;390333
    >>>
    '''  # noqa: E501
    idlist = '\n'.join(map(str, ids))
    arglist = ['taxonkit', 'lineage', '--show-lineage-taxids', '--show-rank', '--show-status-code']
    if debug:
        log(*arglist)  # pragma: no cover
    with NamedTemporaryFile(suffix='-lineage.txt') as lineagefile:
        proc = Popen(arglist, stdin=PIPE, stdout=lineagefile, stderr=PIPE, universal_newlines=True)
        out, err = proc.communicate(input=idlist)
        if proc.returncode != 0:
            raise TaxonKitCLIError(err)  # pragma: no cover
        lineagefile.flush()
        os.fsync(lineagefile.fileno())
        formatargs = []
        if formatstr:
            formatargs = ['--format', formatstr]
        arglist = [
            'taxonkit', 'reformat', *formatargs, '--lineage-field', '3', '--show-lineage-taxids',
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


def name2taxid(names, sciname=False, debug=False):
    '''query taxid by taxon scientific name

    By default, both scientific names and synonyms are supported. Set `sciname=True` to ignore
    synonyms.

    >>> import pytaxonkit
    >>> names = [
    ...     'Homo sapiens', 'Akkermansia muciniphila ATCC BAA-835',
    ...     'Akkermansia muciniphila', 'Mouse Intracisternal A-particle',
    ...     'Wei Shen', 'uncultured murine large bowel bacterium BAC 54B',
    ...     'Croceibacter phage P2559Y'
    ... ]
    >>> result = pytaxonkit.name2taxid(names, debug=True)
    >>> result
                                                  Name    TaxID     Rank
    0                                     Homo sapiens     9606  species
    1             Akkermansia muciniphila ATCC BAA-835   349741  no rank
    2                          Akkermansia muciniphila   239935  species
    3                  Mouse Intracisternal A-particle    11932  species
    4                                         Wei Shen     <NA>     <NA>
    5  uncultured murine large bowel bacterium BAC 54B   314101  species
    6                        Croceibacter phage P2559Y  1327037  species
    >>>
    >>>
    >>> names = ['Phyllobolus spinuliferus', 'Alteromonas putrefaciens', 'Rexia erectus']
    >>> result = pytaxonkit.name2taxid(names)
    >>> result
                           Name   TaxID     Rank
    0  Phyllobolus spinuliferus  359607  species
    1  Alteromonas putrefaciens      24  species
    2             Rexia erectus  262902  species
    >>> result = pytaxonkit.name2taxid(names, sciname=True)
    >>> result
                           Name  TaxID  Rank
    0  Phyllobolus spinuliferus   <NA>  <NA>
    1  Alteromonas putrefaciens   <NA>  <NA>
    2             Rexia erectus   <NA>  <NA>
    >>>
    '''
    namelist = '\n'.join(map(str, names))
    arglist = ['taxonkit', 'name2taxid', '--show-rank']
    if sciname:
        arglist.append('--sci-name')
    if debug:
        log(*arglist)  # pragma: no cover
    proc = Popen(arglist, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=namelist)
    if proc.returncode != 0:
        raise TaxonKitCLIError(err)  # pragma: no cover
    columns = {
        'Name': pandas.StringDtype(),
        'TaxID': pandas.UInt32Dtype(),
        'Rank': pandas.StringDtype()
    }
    data = pandas.read_csv(
        StringIO(out), sep='\t', header=None, names=columns, dtype=columns, index_col=False
    )
    return data
