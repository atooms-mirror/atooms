VERSION = $$(git describe --abbrev=4 --dirty --always 2>/dev/null || echo 0)

all: version

test:
	cd ..;	make test
version:
	@echo __version__ = \'$(VERSION)\' > _version.py
install:
	cd ..;	python setup.py install --home=~
todo:
	todo.py .
clean:
	rm -f *.pyc