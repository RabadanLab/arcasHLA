all: format test

.PHONY: format lint test

# Requires dev dependencies to be installed.
format:
	# Markdown files.
	fd ".*\.md$$" -H --exec mdformat {}

	# Python scripts.
	black scripts/*.py test/*.py

	# Yaml files.
	fd ".*\.ya?ml$$" -H --exclude "results" --exec yamlfmt {}

# Requires dev dependencies to be installed.
lint:
	# Markdown files.
	fd ".*\.md$$" -H --exec mdformat --check {}

	# Python scripts.
	black --check scripts/*.py test/*.py

	# Yaml files.
	fd ".*\.ya?ml$$" -H --exclude "results" --exec yamlfmt --lint {}

# Requires conda-lock to be installed.
lock:
	# Plain environment.
	conda lock \
		--file environment.yml \
		--lockfile environment.conda-lock.yml \
		--platform linux-64

	# Full dev environment.
	conda lock \
		--file environment.yml \
		--file dev_environment.yml \
		--lockfile dev_environment.conda-lock.yml \
		--platform linux-64

test:
	pytest
