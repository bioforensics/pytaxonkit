## help:     print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:     execute test suite
test:
	COLUMNS=150 pytest --cov=pytaxonkit --doctest-modules pytaxonkit.py conftest.py
test4:
	COLUMNS=150 pytest -n 4 --cov=pytaxonkit --doctest-modules pytaxonkit.py conftest.py


## style:    check code style against PEP8
style:
	pycodestyle --max-line-length=99 pytaxonkit.py


## devenv:   set up a development environment
devenv:
	conda install -y pycodestyle pytest pytest-cov
	pip install -e .
	echo 'make style' >> .git/hooks/pre-commit
	chmod 755 .git/hooks/pre-commit
