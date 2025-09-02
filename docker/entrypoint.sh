#!/bin/sh -l

echo "Gene-id: $1"
version=$(bedtools --version)
echo "version=$version" >> $GITHUB_OUTPUT
