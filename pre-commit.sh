#!/bin/sh

# Check if this is the initial commit
if git rev-parse --verify HEAD >/dev/null 2>&1
then
    echo "pre-commit: About to create a new commit..."
    against=HEAD
else
    echo "pre-commit: About to create the first commit..."
    against=4b825dc642cb6eb9a060e54bf8d69288fbee4904
fi

# Use git diff-index to check for whitespace errors
echo "pre-commit: Testing for syntax errors..."
if ! ansible-playbook -i "localhost," galaxy.yml --syntax-check
then
    echo "pre-commit: Aborting commit due playbook syntax error"
    exit 1
else
    echo "pre-commit: No syntax errors :)"
    exit 0
fi
