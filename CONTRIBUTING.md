Contributions
=============

Web platform
------------

The development of atooms takes place on [framagit](https://framagit.org/atooms) - a GitLab instance. Contributions from other platforms are welcome.

If you wish to contribute to the core library or its related packages, you can fork the repository directly on framagit: just [sign in](https://framagit.org/users/sign_in) using your GitHub or GitLab account - no registration is needed. Otherwise, point me to your repository on your platform of choice and I'll grab the code from there.

Workflow
--------

To contribute with new features, bug fixes or documentation, please follow this basic workflow (very close to [git flow](https://guides.github.com/introduction/flow/))

1. make sure your `master` branch is updated
2. create a branch from `master`
3. *hack, hack, hack...*
4. make sure the code passes all the tests
5. open a pull request or point me to your branch

Tests are automatically performed once you push your branch to [framagit](https://framagit.org) (you will see a green check close to the commit hash), but I actually recommend testing the code before pushing. You can run the tests from the command line with `make test`.

Style
-----

atooms seeks to provide an expressive and consistent interface. If a method doesn't clearly express its goal or an attribute doesn't feel natural, then I consider it a bug. Therefore, when adding new interface to the code, always use expressive, non abbreviated names for public methods and variables - even if they are long! Ex: `system.total_energy` is nice, `system.etot` is not.

The general [PEP8 style guidelines](https://www.python.org/dev/peps/pep-0008/) apply - with a grain of salt. For instance, lines up to 120 characters are fine with me. Imports in functions too. It is typically safe to execute `make pep8` (if `autopep8` and `flake8` are installed) to fix some cosmetic, whitespace issues and get a report on the outstanding PEP8 issues. The PEP8 issues I ignore are listed in `setup.cfg`.

Please follow [these simple rules](https://chris.beams.io/posts/git-commit/) to write commit messages. Try to use the verb "add" to start a commit describing a new feature (ex: "Add new trajectory class") and "fix" to describe a bug fix (ex. "Fix writing cell side in NPT simulations"). Do not leave empty commit message. Consider rebasing your branch to reword or squash commits, before pushing.
