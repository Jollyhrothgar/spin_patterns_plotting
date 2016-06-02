#!/bin/sh
srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

(cd $srcdir; aclocal -I /usr/share;\
libtoolize --force; automake -a --add-missing; autoconf)

$srcdir/configure  "$@"
