COMMIT = $$(git describe --abbrev=6 --always 2>/dev/null || echo 0)
VERSION = $$(git describe --abbrev=6 --dirty --always 2>/dev/null || echo 0)
DATE=$$(git show -s --format=%ci $(COMMIT))

all: version

test:
	cd ..;	make test
version:
	@echo __version__ = \'$(VERSION)\' > _version.py
	@echo __date__ = \'$(DATE)\' >> _version.py
install:
	cd ..;	python setup.py install --home=~
todo:
	todo.py .
clean:
	rm -f *.pyc
