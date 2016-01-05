echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi


if [[ "2.7 3.4" =~ "$python" ]]; then
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user iModels --package demos-dev $HOME/miniconda/conda-bld/linux-64/demos-*
    conda convert $HOME/miniconda/conda-bld/linux-64/demos-* -p all
    ls
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user iModels --package demos-dev linux-32/demos-*
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user iModels --package demos-dev win-32/demos-*
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user iModels --package demos-dev win-64/demos-*
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user iModels --package demos-dev osx-64/demos-*
fi

if [[ "$python" != "2.7" ]]; then
    echo "No deploy on PYTHON_VERSION=${python}"; exit 0
fi
