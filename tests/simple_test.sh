#! /bin/bash

# look for TakeABreak binary. In devel mode, it's in ../build/bin directory.
# In production mode, it's in ../bin directory.
if [ -f "../bin/TakeABreak" ]
then
 bindir="../bin"
elif [ -f "../build/bin/TakeABreak" ]
then
 bindir="../build/bin"
else
 echo "could not find a compiled TakeABreak binary"
 exit 1
fi

################################################################################
# we launch the tool
################################################################################
${bindir}/TakeABreak -in gold/test4.fasta -out test4.takeabreak

################################################################################
# we check the result when md5sum is available (Linux: OK ; OSX: KO)
################################################################################
RETVAL=0
MD5_PATH=`which md5sum`
if [ ! -z "$MD5_PATH" ] ; then
  echo "ec396d93ef253a7c93c4f08edf2a82a7  test4.takeabreak.fasta" > ref
  md5sum  test4.takeabreak.fasta > check

  diff ./ref ./check
  if [ $? -eq 0 ]; then
     echo "TEST OK"
  else
     echo "TEST KO"
     RETVAL=1
  fi
fi

################################################################################
# clean up
################################################################################
rm -f  test4* ref check

# for Jenkins CI platform, we need an exit code: PASS (0) vs. FAILED (1)
exit $RETVAL
