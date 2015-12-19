# THIS IS AN AUTOMATICALLY GENERATED FILE.
# DO NOT EDIT! CHANGES WILL BE LOST.
set ocwd="$cwd"
if ($?LIBTBX_BUILD_RELOCATION_HINT) then
  cd "$LIBTBX_BUILD_RELOCATION_HINT"
  unsetenv LIBTBX_BUILD_RELOCATION_HINT
else
  cd "/Users/bpaynter/ownCloud/_School/_Research/Computation/BloodClotting/code/fable_build"
endif
setenv LIBTBX_BUILD "`/bin/sh -c 'pwd -P'`"
setenv LIBTBX_OPATH "$PATH"
setenv PATH "$LIBTBX_BUILD/bin:$PATH"
cd "$ocwd"
unset ocwd
alias libtbx.setpaths_all "source '$LIBTBX_BUILD/setpaths_all.csh'"
alias libtbx.unsetpaths "source '$LIBTBX_BUILD/unsetpaths.csh'"
if ($?PYTHONPATH) then
  setenv LIBTBX_TMPVAL "$PYTHONPATH"
else
  unsetenv LIBTBX_TMPVAL
endif
setenv PYTHONPATH "`libtbx.path_utility prepend LIBTBX_TMPVAL '$LIBTBX_BUILD/../fable_sources:$LIBTBX_BUILD/../fable_sources/libtbx/pythonpath:$LIBTBX_BUILD/lib' < /dev/null`"
if ("$PYTHONPATH" == "L_I_B_T_B_X_E_M_P_T_Y") unsetenv PYTHONPATH
if ($?DYLD_LIBRARY_PATH) then
  setenv LIBTBX_TMPVAL "$DYLD_LIBRARY_PATH"
else
  unsetenv LIBTBX_TMPVAL
endif
setenv DYLD_LIBRARY_PATH "`libtbx.path_utility prepend LIBTBX_TMPVAL '$LIBTBX_BUILD/lib:/System/Library/Frameworks/Python.framework/Versions/2.7/lib' < /dev/null`"
if ("$DYLD_LIBRARY_PATH" == "L_I_B_T_B_X_E_M_P_T_Y") unsetenv DYLD_LIBRARY_PATH
setenv LIBTBX_DIST "$LIBTBX_BUILD/../fable_sources/libtbx"
setenv FABLE_DIST "$LIBTBX_BUILD/../fable_sources/fable"
if ($?LIBTBX_OPATH) then
  setenv LIBTBX_TMPVAL "$LIBTBX_OPATH"
else
  unsetenv LIBTBX_TMPVAL
endif
setenv PATH "`libtbx.path_utility prepend LIBTBX_TMPVAL '$LIBTBX_BUILD/bin' < /dev/null`"
if ("$PATH" == "L_I_B_T_B_X_E_M_P_T_Y") unsetenv PATH
unsetenv LIBTBX_TMPVAL
unsetenv LIBTBX_OPATH
