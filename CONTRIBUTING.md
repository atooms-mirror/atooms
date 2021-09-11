Community contributions
=======================

Web platform
------------

The development of atooms takes place on [framagit](https://framagit.org/atooms) - a GitLab instance. Contributions from other platforms are welcome.

If you wish to contribute to the core library or its related packages, you can fork the repository directly on framagit: just [sign in](https://framagit.org/users/sign_in) using your GitHub or GitLab account - no registration is needed. Otherwise, point me to your repository on your platform of choice and I'll grab your code from there.

Workflow
--------

To contribute with new features and bug fixes, please follow this basic workflow (very close to [git flow](https://guides.github.com/introduction/flow/))

1. make sure your `master` branch is updated
2. create a branch from `master`
3. /hack, hack, hack.../
4. make sure the code passes of the tests
5. open a pull request or point me to your branch

Tests are automatically performed once you push your branch to [framagit](https://framagit.org) (you will see a green check close to the commit hash), but I actually recommend testing the code before pushing the changes. You can run the tests from the command line with `make test`.
