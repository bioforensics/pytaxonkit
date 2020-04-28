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
import json
from subprocess import Popen, PIPE
import sys

from _version import get_versions
__version__ = get_versions()['version']
del get_versions


def log(*args, level='debug'):  # pragma: no cover
    print(f'[pytaxonkit::{level}]', *args, file=sys.stderr)


BasicTaxon = namedtuple('Taxon', ['taxid', 'rank', 'name'])


# -------------------------------------------------------------------------------------------------
# taxonkit list
# -------------------------------------------------------------------------------------------------

class ListResult():
    def __init__(self, jsondata):
        self._data = json.loads(jsondata)

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

    By default, a `pytaxonkit.ListResult` object is returned. This provides a convenient mechanism
    for programmatic access of the results—see the doctest below. NOTE: for very large results,
    these objects introduce a significant amount of overhead—in such cases, it is recommended to
    set `raw=True` so that the raw JSON result is returned.

    >>> import pytaxonkit
    >>> result = pytaxonkit.list([8204, 2468])
    >>> for taxon, tree in result:
    ...   print('[Result]', taxon.name, taxon.taxid, taxon.rank, sep='\t')
    ...
    [Result]    Anarhichas lupus    8204    species
    [Result]    Plasmid NR79        2468    species
    >>>
    >>>
    >>> result = pytaxonkit.list(['9605', '239934'])
    >>> for taxon, tree in result:
    ...     print(f'[Top level result] {taxon.name} ({taxon.taxid})')
    ...     subtaxa = [t for t in tree.traverse]
    ...     n = len(subtaxa)
    ...     print(f'    - {n} sub taxa; first 5:')
    ...     for subtaxon in subtaxa[:5]:
    ...         print('        -', subtaxon)
    ...
    [Top level result] Homo (9605)
        - 6 sub taxa; first 5:
            - Taxon(taxid=9606, rank='species', name='Homo sapiens')
            - Taxon(taxid=63221, rank='subspecies', name='Homo sapiens neanderthalensis')
            - Taxon(taxid=741158, rank='subspecies', name="Homo sapiens subsp. 'Denisova'")
            - Taxon(taxid=2665952, rank='no rank', name='environmental samples')
            - Taxon(taxid=2665953, rank='species', name='Homo sapiens environmental sample')
    [Top level result] Akkermansia (239934)
        - 90 sub taxa; first 5:
            - Taxon(taxid=239935, rank='species', name='Akkermansia muciniphila')
            - Taxon(taxid=349741, rank='no rank', name='Akkermansia muciniphila ATCC BAA-835')
            - Taxon(taxid=512293, rank='no rank', name='environmental samples')
            - Taxon(taxid=512294, rank='species', name='uncultured Akkermansia sp.')
            - Taxon(taxid=1131822, rank='species', name='uncultured Akkermansia sp. SMG25')
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
    {'9605 [genus] Homo': {'9606 [species] Homo sapiens': {'63221 [subspecies] Homo sapiens
    neanderthalensis': {}, "741158 [subspecies] Homo sapiens subsp. 'Denisova'": {}, '2665952 [no
    rank] environmental samples': {'2665953 [species] Homo sapiens environmental sample': {}}},
    '1425170 [species] Homo heidelbergensis': {}}}
    >>>
    '''
    idlist = ','.join(map(str, ids))
    arglist = ['taxonkit', 'list', '--json', '--show-name', '--show-rank', '--ids', idlist]
    if debug:
        log(*arglist)  # pragma: no cover
    proc = Popen(arglist, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate()
    if raw:
        return json.loads(out)
    else:
        return ListResult(out)
