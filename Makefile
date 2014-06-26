all: install

dist:
	python setup.py sdist

test:
	python -m unittest discover -s tests

todo:
	todo.py pyatooms

dist_rumd:
	python -m unittest discover -s tests -p '*adapter*'
	tar cvf adapter_rumd.tar tests/test_adapters.py pyatooms/adapters/rumd.py

install:
	python setup.py install --home=~
