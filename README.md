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
> - The `reformat` operation is automatically run by `pytaxonkit.lineage` and cannot be invoked independently.
> - The `taxid-changelog` operation is not supported.
> - The `genautocomplete` operation is specific to the shell and is not supported.

### name2taxid

```python
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
```

### lineage

```python
>>> import pytaxonkit
>>> result = pytaxonkit.lineage([7399, 1973489])
>>> result.columns
Index(['TaxID', 'Code', 'Lineage', 'LineageTaxIDs', 'Rank', 'FullLineage',
       'FullLineageTaxIDs'],
      dtype='object')
>>> result[['TaxID', 'Lineage', 'LineageTaxIDs']]
     TaxID                                            Lineage                          LineageTaxIDs
0     7399        Eukaryota;Arthropoda;Insecta;Hymenoptera;;;                2759;6656;50557;7399;;;
1  1973489  Bacteria;Firmicutes;Bacilli;Bacillales;Bacilla...  2;1239;91061;1385;186817;1386;1973489
>>> result = pytaxonkit.lineage(['1382510', '929505', '390333'], formatstr='{f};{g};{s};{S}')
>>> result['Lineage'].iloc[2]
'Lactobacillaceae;Lactobacillus;Lactobacillus delbrueckii;Lactobacillus delbrueckii subsp. bulgaricus'
```

### list

```python
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
{'9605 [genus] Homo': {'9606 [species] Homo sapiens': {'63221 [subspecies] Homo sapiens neanderthalensis': {}, "741158 [subspecies]Homo sapiens subsp. 'Denisova'": {}, '2665952 [no rank] environmental samples': {'2665953 [species] Homo sapiens environmentalsample': {}}}, '1425170 [species] Homo heidelbergensis': {}}}
```

### version

```python
>>> import pytaxonkit
>>> pytaxonkit.__version__
0.6
>>> pytaxonkit.__taxonkitversion__
0.6.0
```


## License

[BSD 3-clause](LICENSE)
