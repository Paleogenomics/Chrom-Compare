#! /bin/bash


COUNTER=0
TOTAL=0

# Test 1
echo 'Test 1: Input files all identical, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/first.fa -3 t/first.fa -4 t/first.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 130 0 0 0 0 0 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 2
echo 'Test 2: 1st input file different from others, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/second.fa -2 t/first.fa -3 t/first.fa -4 t/first.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 122 8 0 0 0 0 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 3
echo 'Test 3: 2nd input file different from others, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/second.fa -3 t/first.fa -4 t/first.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 122 0 8 0 0 0 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 4
echo 'Test 4: 3rd input file different from others, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/first.fa -3 t/second.fa -4 t/first.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 122 0 0 8 0 0 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 5
echo 'Test 5: 4th input file different from others, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/first.fa -3 t/first.fa -4 t/second.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 122 0 0 0 8 0 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 6
echo 'Test 6: 1,2 input file different from others, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/second.fa -2 t/second.fa -3 t/first.fa -4 t/first.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 122 0 0 0 0 8 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 7
echo 'Test 7: 1,3 input file different from others, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/second.fa -2 t/first.fa -3 t/second.fa -4 t/first.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 122 0 0 0 0 0 8 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 8
echo 'Test 8: 1,4 input file different from others, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/second.fa -2 t/first.fa -3 t/first.fa -4 t/second.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 122 0 0 0 0 0 0 8'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 9
echo 'Test 9: 1,2,3,3 input files, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/second.fa -3 t/third.fa -4 t/third.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 118 4 4 0 0 4 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 10
echo 'Test 10: 1,2,2,3 input files, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/second.fa -3 t/second.fa -4 t/third.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 118 4 0 0 4 0 0 4'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 11
echo 'Test 11: 1,1,2,3 input files, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/first.fa -3 t/second.fa -4 t/third.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 118 0 0 4 4 4 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 12
echo 'Test 12: 1,2,1,3 input files, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/second.fa -3 t/first.fa -4 t/third.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 118 0 4 0 4 0 4 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 13
echo 'Test 13: 1,2,3,1 input files, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/second.fa -3 t/third.fa -4 t/first.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 118 0 4 4 0 0 0 4'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 14
echo 'Test 14: 1,2,3,2 input files, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/second.fa -3 t/third.fa -4 t/second.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 118 4 0 4 0 0 4 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 15
echo 'Test 15: 1,2,3,4 input files, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/second.fa -3 t/third.fa -4 t/fourth.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 118 3 3 4 0 0 1 1'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 16
echo 'Test 16: 1st file is wrapped differently, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first_wrap.fa -2 t/second.fa -3 t/third.fa -4 t/fourth.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 118 3 3 4 0 0 1 1'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi


# Test 17
echo 'Test 17: 1st file is shorter than others, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first_short.fa -2 t/second.fa -3 t/third.fa -4 t/fourth.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 120 109 3 3 3 0 0 1 1'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 18
echo 'Test 18: 1st file is all lowercase, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first_lower.fa -2 t/second.fa -3 t/third.fa -4 t/fourth.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 118 3 3 4 0 0 1 1'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 19
echo 'Test 19: 2nd file has a string of 8 Ns, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first_lower.fa -2 t/second_Ns.fa -3 t/third.fa -4 t/fourth.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 110 3 3 4 0 0 1 1'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 20
echo 'Test 20: 2nd file has a string of 8 ns, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first_lower.fa -2 t/second_ns.fa -3 t/third.fa -4 t/fourth.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 110 3 3 4 0 0 1 1'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 21
echo 'Test 21: 3rd file has non DNA letters, Whole sequence (no window)'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first_lower.fa -2 t/second.fa -3 t/third_other_char.fa -4 t/fourth.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR='0 130 110 3 3 4 0 0 1 1'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

N_PASS=$[ $TOTAL - $COUNTER ]

echo 
echo "quad-aln-report tests complete:"
echo "$N_PASS out of $TOTAL passed"
exit $COUNTER
