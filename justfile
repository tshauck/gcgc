# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

test:
	pytest -v -s --cov-report term-missing --cov=gcgc

clean:
	rm -rf dist

publish-python:
	poetry build
	poetry publish

publish-docker:
	docker-compose build gcgc
	docker-compose push gcgc
	TAG=latest docker-compose build gcgc
	TAG=latest docker-compose push gcgc

publish: publish-python publish-docker

isort:
	isort

docs: clean_docs
	cp ./CHANGELOG.md ./docs/
	cp ./README.md ./docs/index.md
	mkdocs build

clean_docs:
	rm -rf site

docs_upload: docs
	aws s3 cp --recursive ./site s3://gcgc.trenthauck.com/

docs_write_good:
	write-good ./docs/**/*.md

vulture:
	vulture gcgc

pydocstyle:
	pydocstyle --add-ignore=D202,D203,D301 --convention=google gcgc

pylint:
	pylint gcgc

mypy:
	mypy gcgc

test_integration:
	pytest -m 'integration'

test_unit:
	pytest --cov-report term-missing --cov=gcgc -m 'not integration'

fmt:
	black .

fmt-check:
	black --check .

build:
	docker-compose build gcgc

docker-tests: build
	docker-compose run --rm gcgc test
	docker-compose run --rm gcgc fmt-check
	docker-compose run --rm gcgc pydocstyle
	docker-compose run --rm gcgc pylint
