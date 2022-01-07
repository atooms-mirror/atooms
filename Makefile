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
	rm -rf docs/api
	pdoc -o docs/api --force --html --skip-errors $(PROJECT)
	emacs docs/basics.org --batch -l ~/.emacs -l ~/.emacs.d/org-mode.el -f org-rst-export-to-rst --kill
	emacs docs/simulations.org --batch -l ~/.emacs -l ~/.emacs.d/org-mode.el -f org-rst-export-to-rst --kill
	emacs docs/trajectories.org --batch -l ~/.emacs -l ~/.emacs.d/org-mode.el -f org-rst-export-to-rst --kill
	orgnb.py docs/*.org
	make -C docs/ html

version:
	@echo __commit__ = \'$(COMMIT_DIRTY)\' > atooms/core/_commit.py
	@echo __date__ = \'$(DATE)\' >> atooms/core/_commit.py

pep8:
	autopep8 -r -i $(PROJECT)
	autopep8 -r -i tests
	flake8 $(PROJECT)

clean:
	rm -rf atooms/*pyc atooms/*/*pyc tests/*pyc atooms/*/*pyo atooms/*/__pycache__
