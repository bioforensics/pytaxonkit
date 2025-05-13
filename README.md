# pytaxonkit

Python bindings for the [TaxonKit library](https://bioinf.shenwei.me/taxonkit/).
Results of queries are returned as [pandas](https://pandas.pydata.org/) data frames, dictionaries, and other convenient data structures.


## Install

Installation with Conda is recommended.
(See `environment.yaml` for details on prerequisites if you want to try a different installation method.)

```
conda install -c bioconda pytaxonkit
```

Please follow the [taxonkit instructions](https://bioinf.shenwei.me/taxonkit/usage/#taxonkit) for downloading the NCBI Taxonomy database dump files.
By default, taxonkit and pytaxonkit search for these files in `~/.taxonkit/`, but this directory can be overridden by setting the `data_dir` appropriately for the functions described below.
Execute `help(pytaxonkit.name2taxid)` (and so on) from the Python interpreter for more details.


## Usage

`pytaxonkit` provides convenient access to TaxonKit operations from Python for use in the interactive interpreter, IPython, Jupyter notebooks, or custom Python code.

> **NOTES**
> - The `reformat2` operation is automatically run by `pytaxonkit.lineage` and cannot be invoked independently.
> - The `genautocomplete` operation is specific to the shell and is not supported.
> - Several other operations are not supported, including `cami-filter`, `create-taxdump`, `profile2cami`, and `taxid-changelog`.
> - The `pytaxonkit.__version__` variable refers to the version number of the Python bindings, while the `pytaxonkit.__taxonkitversion__` variable corresponds to the version of the installed TaxonKit program. These version numbers are not necessarily equal.

### name2taxid

```python
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
```

### lineage

```python
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
```

### name

```python
>>> import pytaxonkit
>>> name(["151837", "2216222", "517824"])
     TaxID                                Name
0   151837                    Hiraea smilacina
1  2216222         Paramyia sp. BIOUG21706-A10
2   517824  soil bacterium Cipr-S1N-M1LLLSSL-1
```

### list

```python
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
```

### filter

```python
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
```

### lca

```python
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
```

### version

```python
>>> import pytaxonkit
>>> pytaxonkit.__version__
'0.9.1'
>>> pytaxonkit.__taxonkitversion__
'taxonkit v0.20.0'
```


## License

[BSD 3-clause](LICENSE)
