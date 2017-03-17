COMMIT = $$(git describe --abbrev=6 --always 2>/dev/null || echo 0)
COMMIT_DIRTY = $$(git describe --abbrev=6 --dirty --always 2>/dev/null || echo 0)
DATE=$$(git show -s --format=%ci $(COMMIT))

.PHONY: all dist test todo todo_critical dist_rumd install clean

all: version

pull:
	git pull

doc:
        # pdoc does play nice with namespace packages -> blank __init__.py
        # go into atooms to prevent pdoc from populating docs/ with modules from the namespace package
	cd atooms; rm -f __init__.pyc; mv __init__.py __init__.py.bak; touch __init__.py
	cd atooms; pdoc --overwrite --html-dir ../docs --html --template-dir ~/usr/pdoc_templates/pdoc_templates ../atooms 
	cd atooms; mv __init__.py.bak __init__.py

dist:
	python setup.py sdist

test:
	python -m unittest discover -s tests

todo:
	@todo.py -S|grep '^Open'

todo_critical:
	@todo.py -S|grep '!'

dist_rumd:
	python -m unittest discover -s tests -p '*adapter*'
	tar cvf adapter_rumd.tar tests/test_adapters.py atooms/adapters/rumd.py

version:
	@echo __commit__ = \'$(COMMIT_DIRTY)\' > atooms/core/_commit.py
	@echo __date__ = \'$(DATE)\' >> atooms/core/_commit.py

install: version
	python setup.py install --user

develop: version
	python setup.py develop --user

clean:
	rm -f atooms/*pyc  atooms/*/*pyc tests/*pyc
