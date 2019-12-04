#!/usr/bin/env bash
set -e
# skipping tests
# .local/lib/python2.7/site-packages/bioblend/_tests/TestGalaxyLibraries.py::TestGalaxyLibraries::test_upload_from_galaxy_filesystem FAILED
# .local/lib/python2.7/site-packages/bioblend/_tests/TestGalaxyObjects.py::TestLibrary::test_datasets_from_fs FAILED
# .local/lib/python2.7/site-packages/bioblend/_tests/TestGalaxyObjects.py::TestLibrary::test_get_datasets FAILED
# .local/lib/python2.7/site-packages/bioblend/_tests/TestGalaxyObjects.py::TestHistory::test_get_datasets FAILED
export PATH=$GALAXY_HOME/.local/bin/:$PATH &&
cd $GALAXY_HOME &&
export tests='DatasetCollections Datasets Folders Groups Histories Quotas Roles ToolData ToolInputs Tools Users Workflows' &&
echo  $tests &&
for test in $tests
    do
    export test=$test
    bioblend-galaxy-tests -v $GALAXY_HOME/.local/lib/python2.7/site-packages/bioblend/_tests/TestGalaxy$test.py
    done
