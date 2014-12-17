.PHONY: all dist test todo todo_critical dist_rumd install clean

all: install

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

install:
	python setup.py install --home=~

clean:
	rm -f atooms/*pyc  atooms/*/*pyc tests/*pyc
