# (c) Copyright 2018 Trent Hauck
# All Rights Reserved

.PHONY: test
test:
	pytest -v -s --cov-report term-missing --cov=gcgc


.PHONY: build
build:
	poetry build

.PHONY: publish
publish:
	poetry publish
