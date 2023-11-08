# this is a bash script to setup a ASTR case

ASTR_DIR="$(dirname "$(readlink -f "$0")")"

echo "ASTR directory: $ASTR_DIR"

cp -v $ASTR_DIR/Makefile  ./
cp -vr $ASTR_DIR//examples/Taylor_Green_Vortex/datin ./

mkdir "./user_define_module"
cp -v $ASTR_DIR/user_define_module/userdefine.F90 ./user_define_module/