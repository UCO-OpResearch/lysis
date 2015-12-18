# LIBTBX_DISPATCHER_HEAD DO NOT EDIT
#
# This file is intended to be sourced from other scripts.
# It is like the dispatcher scripts in the bin directory,
# but only sets up the environment without calling a
# command at the end.
#
# THIS IS AN AUTOMATICALLY GENERATED FILE.
# DO NOT EDIT! CHANGES WILL BE LOST.
# To customize this auto-generated script create
#
#   dispatcher_include*.sh
#
# files in "/Users/bpaynter/ownCloud/_School/_Research/Computation/BloodClotting/code/fable_build" and run
#
#   libtbx.refresh
#
# to re-generate the dispatchers (libtbx.refresh is a subset
# of the functionality of the libtbx/configure.py command).
#
# See also:
#   "/Users/bpaynter/ownCloud/_School/_Research/Computation/BloodClotting/code/fable_build/dispatcher_include_template.sh"
#

# ----------------------------------------------------------------------------
# The shellrealpath function resolves an absolute physical path of its
# first argument and stores it in a global shell variable RESULT.
# The function returns nonzero for unreadable or invalid symlinks
# and resets the RESULT to an empty string.

shellrealpath() {
    local ORGDIR="$PWD"
    local TARGET="$1"
    RESULT=""
    # This test fails for a symlink loop.  We can do without resolution
    # of symlinks that point to existing unreadable files.
    [ -r "$TARGET" ] || return $?
    # Check if the readlink command exists.
    type readlink >/dev/null || return $?
    while true; do
        cd "$(dirname "$TARGET")"
        TARGET="$(basename "$TARGET")"
        if [ -L "$TARGET" ]; then
            TARGET="$(readlink "$TARGET")"
            continue
        fi
        RESULT="$(pwd -P)/$TARGET"
        break
    done
    cd "$ORGDIR"
}
# ----------------------------------------------------------------------------

unset PYTHONHOME
LC_ALL=C
export LC_ALL
LIBTBX_BUILD="$(shellrealpath "$0" && cd "$(dirname "$RESULT")/.." && pwd)"
export LIBTBX_BUILD
LIBTBX_PYEXE_BASENAME="python"
export LIBTBX_PYEXE_BASENAME
if [ -n "$PYTHONPATH" ]; then
  PYTHONPATH="$LIBTBX_BUILD/../fable_sources:$LIBTBX_BUILD/../fable_sources/libtbx/pythonpath:$LIBTBX_BUILD/lib:$PYTHONPATH"
  export PYTHONPATH
else
  PYTHONPATH="$LIBTBX_BUILD/../fable_sources:$LIBTBX_BUILD/../fable_sources/libtbx/pythonpath:$LIBTBX_BUILD/lib"
  export PYTHONPATH
fi
if [ -n "$DYLD_LIBRARY_PATH" ]; then
  DYLD_LIBRARY_PATH="$LIBTBX_BUILD/lib:/System/Library/Frameworks/Python.framework/Versions/2.7/lib:$DYLD_LIBRARY_PATH"
  export DYLD_LIBRARY_PATH
else
  DYLD_LIBRARY_PATH="$LIBTBX_BUILD/lib:/System/Library/Frameworks/Python.framework/Versions/2.7/lib"
  export DYLD_LIBRARY_PATH
fi
if [ -n "$PATH" ]; then
  PATH="$LIBTBX_BUILD/bin:$PATH"
  export PATH
else
  PATH="$LIBTBX_BUILD/bin"
  export PATH
fi
LIBTBX_PYEXE="/usr/bin/$LIBTBX_PYEXE_BASENAME"
export LIBTBX_PYEXE
