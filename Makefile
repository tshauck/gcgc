# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

.PHONY: test
test:
	pytest -v -s --cov-report term-missing --cov=gcgc

.PHONY: clean
clean:
	rm -rf dist

.PHONY: build
build: clean
	poetry build

.PHONY: publish
publish:
	poetry publish

.PHONY: dev_release
dev_release: dev_version build publish

.PHONY: isort
isort:
	isort

.PHONY: docs
docs:
	mkdocs build

.PHONY: clean_docs
clean_docs:
	rm -rf site

.PHONY: docs_upload
docs_upload:
	aws s3 cp --recursive ./site s3://gcgc.trenthauck.com/

.PHONY: docs_write_good
docs_write_good:
	write-good ./docs/**/*.md

.PHONY: vulture
vulture:
	vulture gcgc

.PHONY: pydocstyle
pydocstyle:
	pydocstyle gcgc

.PHONY: mypy
mypy:
	mypy gcgc

.PHONY: test_integration
test_integration:
	pytest -m 'integration'

.PHONY: test_unit
test_unit:
	pytest --cov-report term-missing --cov=gcgc -m 'not integration'
