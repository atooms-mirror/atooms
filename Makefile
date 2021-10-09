PROJECT = atooms
COMMIT = $$(git describe --abbrev=6 --always 2>/dev/null || echo 0)
COMMIT_DIRTY = $$(git describe --abbrev=6 --dirty --always 2>/dev/null || echo 0)
DATE=$$(git show -s --format=%ci $(COMMIT) | cut -d ' ' -f 1)

.PHONY: all test todo install docs version pep8 clean

all: version

install: version
	python setup.py install

test:	version
	coverage run --source atooms -m unittest discover -s tests
	coverage report --omit=atooms/backends/*py

docs:
	pdoc -o docs/api --force --html --skip-errors $(PROJECT)
	sed -i '/^$$/d' docs/index.html
	orgnb.py docs/index.org docs/atooms.ipynb

version:
	@echo __commit__ = \'$(COMMIT_DIRTY)\' > atooms/core/_commit.py
	@echo __date__ = \'$(DATE)\' >> atooms/core/_commit.py

pep8:
	autopep8 -r -i $(PROJECT)
	autopep8 -r -i tests
	flake8 $(PROJECT)

clean:
	rm -rf atooms/*pyc atooms/*/*pyc tests/*pyc atooms/*/*pyo atooms/*/__pycache__
