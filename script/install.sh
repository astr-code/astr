# this is a bash script to setup a ASTR case

ASTR_DIR="$(readlink -f "$1")"
if [ ${#ASTR_DIR} -lt 1 ]; then
    echo "** input the path to astr **"
    exit 1
fi
echo "ASTR directory: $ASTR_DIR"
LOCAL_DIR="$(dirname "$(readlink -f "pwd")")"
echo "LOCAL directory: $LOCAL_DIR"
cp -v $ASTR_DIR'/Makefile' $LOCAL_DIR
cp -v $ASTR_DIR/Makefile.astr $LOCAL_DIR
cp -v $ASTR_DIR/Makefile.pastr $LOCAL_DIR
cp -rv $ASTR_DIR'/user_define_module' $LOCAL_DIR
cp -rv $ASTR_DIR'/pastr/user_define_module' $LOCAL_DIR
OLD="SRCDIR = ./astr/src ./astr/user_define_module"
NEW='SRCDIR = '$ASTR_DIR'/src '$LOCAL_DIR'/user_define_module'
echo "src path for astr: $NEW"
sed -i "s|$OLD|$NEW|g" Makefile.astr
OLD="SRCDIR = ./astr/pastr/src ./astr/pastr/user_define_module"
NEW='SRCDIR = '$ASTR_DIR'/pastr/src '$LOCAL_DIR'/user_define_module'
echo "src path for pastr: $NEW"
sed -i "s|$OLD|$NEW|g" Makefile.pastr

# cp -vr $ASTR_DIR//examples/Taylor_Green_Vortex/datin ./

# mkdir "./user_define_module"
# cp -v $ASTR_DIR/user_define_module/userdefine.F90 ./user_define_module/