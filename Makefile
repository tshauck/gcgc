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
	flit build

.PHONY: publish
publish:
	python setup.py sdist bdist_wheel
	twine upload dist/*

.PHONY: dev_release
dev_release: dev_version build publish

.PHONY: isort
isort:
	isort

.PHONY: docs
docs: clean_docs
	rm -rf docs/github
	mkdir docs/github
	cp ./CHANGELOG.md ./docs/github
	sed 's/# GCGC/# README/g' README.md > ./docs/github/README.md
	mkdocs build
	rm -rf ./docs/github/

.PHONY: clean_docs
clean_docs:
	rm -rf site

.PHONY: docs_upload
docs_upload: docs
	aws s3 cp --recursive ./site s3://gcgc.trenthauck.com/

.PHONY: docs_write_good
docs_write_good:
	write-good ./docs/**/*.md

.PHONY: vulture
vulture:
	vulture gcgc

.PHONY: pydocstyle
pydocstyle:
	pydocstyle --convention=google gcgc

.PHONY: flake8
flake8:
	flake8 gcgc

.PHONY: mypy
mypy:
	mypy gcgc

.PHONY: test_integration
test_integration:
	pytest -m 'integration'

.PHONY: test_unit
test_unit:
	pytest --cov-report term-missing --cov=gcgc -m 'not integration'

.PHONY: fmt
fmt:
	black .

.PHONY: fmt-check
fmt-check:
	black --check .
