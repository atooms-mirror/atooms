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

docs: clean
        # pdoc does play nice with namespace packages -> blank __init__.py
        # go into atooms to prevent pdoc from populating docs/ with modules from the namespace package
	cd atooms; mv __init__.py __init__.py.bak; echo \"\"\"A framework for simulations of interacting particles.\"\"\" > __init__.py
	cd atooms; pdoc --overwrite --html-dir ../docs/api --html --template-dir ~/usr/pdoc_tpl/pdoc_tpl ../atooms 
	cd atooms; mv __init__.py.bak __init__.py
	rsync -uva docs/api zaphod:public_html

version:
	@echo __commit__ = \'$(COMMIT_DIRTY)\' > atooms/core/_commit.py
	@echo __date__ = \'$(DATE)\' >> atooms/core/_commit.py

pep8:
	autopep8 -r -i $(PROJECT)
	autopep8 -r -i tests
	flake8 $(PROJECT)

clean:
	rm -rf atooms/*pyc atooms/*/*pyc tests/*pyc atooms/*/*pyo atooms/*/__pycache__
