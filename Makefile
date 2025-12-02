## help:     print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:     execute test suite
test:
	COLUMNS=150 pytest --cov=pytaxonkit --doctest-modules pytaxonkit.py conftest.py
testci:
	COLUMNS=150 pytest --verbose --cov=pytaxonkit --doctest-modules pytaxonkit.py conftest.py
test4:
	COLUMNS=150 pytest -n 4 --cov=pytaxonkit --doctest-modules pytaxonkit.py conftest.py


## style:    check code style
style:
	black --line-length=99 --check pytaxonkit.py setup.py


## format:   autoformat Python code
format:
	black --line-length=99 pytaxonkit.py setup.py


## hooks:    install development hooks
hooks:
	echo "set -eo pipefail" > .git/hooks/pre-commit
	echo 'make style' >> .git/hooks/pre-commit
	chmod 755 .git/hooks/pre-commit
