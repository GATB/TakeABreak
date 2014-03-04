rm -f TakeABreak_$1.zip
rm -rf ../../../thirdparty/gatb-core/gatb-core/build
zip -r TakeABreak_$1.zip TakeABreak -x *svn*
echo TakeABreak_$1.zip was created

