COMMIT = $$(git describe --abbrev=6 --always 2>/dev/null || echo 0)
COMMIT_DIRTY = $$(git describe --abbrev=6 --dirty --always 2>/dev/null || echo 0)
DATE=$$(git show -s --format=%ci $(COMMIT))

all: version

test:
	cd ..;	make test

version:
	@echo __commit__ = \'$(COMMIT_DIRTY)\' > core/_commit.py
	@echo __date__ = \'$(DATE)\' >> core/_commit.py

clean:
	rm -f *.pyc */*.pyc
