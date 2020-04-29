# pytaxonkit

Python bindings for the excellent [TaxonKit library](https://bioinf.shenwei.me/taxonkit/).


## Install

Installation with Conda is supported.

```
conda install -c bioconda pytaxonkit
```

## Usage

The TaxonKit library has an excellent CLI.
`pytaxonkit` provides convenient access to TaxonKit operations from Python for use in the interactive interpreter, IPython, Jupyter notebooks, or custom Python code.

```python
>>> import pytaxonkit
>>>
>>> # name2taxid operation
>>> names = ['Phyllobolus spinuliferus', 'Alteromonas putrefaciens', 'Rexia erectus']
>>> result = pytaxonkit.name2taxid(names)
>>> type(result)
<class 'pandas.core.frame.DataFrame'>
>>> result
                       Name   TaxID     Rank
0  Phyllobolus spinuliferus  359607  species
1  Alteromonas putrefaciens      24  species
2             Rexia erectus  262902  species
>>> result = pytaxonkit.name2taxid(names, sciname=True)  # ignore synonyms with `sciname=True`
>>> result
                       Name  TaxID  Rank
0  Phyllobolus spinuliferus   <NA>  <NA>
1  Alteromonas putrefaciens   <NA>  <NA>
2             Rexia erectus   <NA>  <NA>
>>>
>>>
>>> # lineage operation (also invokes reformat)
>>> result = pytaxonkit.lineage(['1382510', '929505', '390333'], formatstr='{f};{g};{s};{S}', debug=True)
>>> result
     TaxID     Code                                                                                               Lineage         LineageTaxIDs     Rank                                                                                                                                                                                                                                                FullLineage                                                 FullLineageTaxIDs
0  1382510  1382510                                                     Enterobacteriaceae;Salmonella;Salmonella bongori;        543;590;54736;  no rank                                    cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella bongori;Salmonella bongori serovar 48:z41:--;Salmonella bongori serovar 48:z41:-- str. RKS3044              131567;2;1224;1236;91347;543;590;54736;41527;1382510
1   929505   929505                                                     Clostridiaceae;Clostridium;Clostridium botulinum;      31979;1485;1491;  no rank                                                        cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium;Clostridium botulinum;Clostridium botulinum C;Clostridium botulinum C str. Stockholm  131567;2;1783272;1239;186801;186802;31979;1485;1491;36828;929505
2   390333   390333  Lactobacillaceae;Lactobacillus;Lactobacillus delbrueckii;Lactobacillus delbrueckii subsp. bulgaricus  33958;1578;1584;1585  no rank  cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus delbrueckii;Lactobacillus delbrueckii subsp. bulgaricus;Lactobacillus delbrueckii subsp. bulgaricus ATCC 11842 = JCM 1002    131567;2;1783272;1239;91061;186826;33958;1578;1584;1585;390333
>>>
>>>
>>> # list operation
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
>>> result = pytaxonkit.list(['9605'], raw=True)  # setting `raw=True` returns only the raw JSON
>>> print(result)
{'9605 [genus] Homo': {'9606 [species] Homo sapiens': {'63221 [subspecies] Homo sapiens neanderthalensis': {}, "741158 [subspecies]Homo sapiens subsp. 'Denisova'": {}, '2665952 [no rank] environmental samples': {'2665953 [species] Homo sapiens environmentalsample': {}}}, '1425170 [species] Homo heidelbergensis': {}}}
```

The `reformat` operation is automatically run by `pytaxonkit.lineage` and cannot be invoked independently.
The `taxid-changelog` operation is not supported.
The `genautocomplete` operation is specific to the shell and is not supported.


## License

[BSD 3-clause](LICENSE)
