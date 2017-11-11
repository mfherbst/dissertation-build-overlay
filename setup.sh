# Patches for building the dissertation without molsturm installed.
#
# 1. Run cmake inside your build directory
# 2. Source this file from the build directory
# 3. Start the build with make

if [ "$BASH" == "" ]; then
	echo "The source script only works from the bash shell." >&1
	return 1
fi

THIS_DIR=`dirname ${BASH_SOURCE[0]}`
cp -r ${THIS_DIR}/python $PWD/python_extra
export PYTHONPATH="$PYTHONPATH:$PWD/python_extra"

for file in ${THIS_DIR}/build/*/*; do
	RELPATH=${file##${THIS_DIR}/build/}
	cp -r $file $PWD/$RELPATH || return 1
done

return 0 &> /dev/null
echo
echo "Now set PYTHONPATH to $PYTHONPATH to start your build."
