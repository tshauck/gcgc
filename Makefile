.PHONY: test
test:
	pytest -v -s --cov-report term-missing --cov=gcgc
