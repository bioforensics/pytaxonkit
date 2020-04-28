## help:    print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:    execute test suite
test:
	pytest --cov=pytaxonkit --doctest-modules pytaxonkit.py
