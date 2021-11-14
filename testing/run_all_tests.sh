#! /bin/sh

source ./ref_tool_glob.sh

TEST_OUT_DIR="test_output"

if [ ! -d "$TEST_OUT_DIR" ]; then
    mkdir "$TEST_OUT_DIR"
else
    rm "$TEST_OUT_DIR"/*
fi

for f in testdata/*.json; do
	TESTNAME=$(basename -s .json $f)
	echo "Running test $TESTNAME"
	LD_LIBRARY_PATH=${ECL_BIN} ./ref_tool $f 0 0 0 &> "test_output/$TESTNAME_noIF.out"
	if [ $? -ne 0 ]; then
		echo "Test $TESTNAME (no ionic effects) FAILED !!!"
	else
		echo "Test $TESTNAME (no ionic effects) PASSED"
	fi

	LD_LIBRARY_PATH=${ECL_BIN} ./ref_tool $f 1 1 0 &> "test_output/$TESTNAME_withIF.out"
	if [ $? -ne 0 ]; then
		echo "Test $TESTNAME (ionic effects) FAILED !!!"
	else
		echo "Test $TESTNAME (ionic effects) PASSED"
	fi

	echo ""
done
