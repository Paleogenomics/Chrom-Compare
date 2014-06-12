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


# Test 22
echo 'Test 22: Input files all identical, -W window_size = 65'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/first.fa -3 t/first.fa -4 t/first.fa -W 65'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR=$'0 65 65 0 0 0 0 0 0 0\n65 130 65 0 0 0 0 0 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi


# Test 23
echo 'Test 23: Input files all identical, -W window_size = 60'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/first.fa -2 t/first.fa -3 t/first.fa -4 t/first.fa -W 60'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR=$'0 60 60 0 0 0 0 0 0 0\n60 120 60 0 0 0 0 0 0 0\n120 130 10 0 0 0 0 0 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi


# Test 24
echo 'Test 24: win tests win1 win1 win2 win2, -W window_size = 80'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/win1.fa -2 t/win1.fa -3 t/win2.fa -4 t/win2.fa -W 80'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR=$'0 80 40 0 0 0 0 40 0 0\n80 160 40 0 0 0 0 40 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 25
echo 'Test 25: win tests win1 win1 win2 win2, -W window_size = 40'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/win1.fa -2 t/win1.fa -3 t/win2.fa -4 t/win2.fa -W 40'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR=$'0 40 20 0 0 0 0 20 0 0\n40 80 20 0 0 0 0 20 0 0\n80 120 20 0 0 0 0 20 0 0\n120 160 20 0 0 0 0 20 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 26
echo 'Test 26: win tests win1 win1 win2 win2, -W window_size = 20'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/win1.fa -2 t/win1.fa -3 t/win2.fa -4 t/win2.fa -W 20'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR=$'0 20 20 0 0 0 0 0 0 0\n20 40 0 0 0 0 0 20 0 0\n40 60 20 0 0 0 0 0 0 0\n60 80 0 0 0 0 0 20 0 0\n80 100 20 0 0 0 0 0 0 0\n100 120 0 0 0 0 0 20 0 0\n120 140 20 0 0 0 0 0 0 0\n140 160 0 0 0 0 0 20 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 27
echo 'Test 27: win tests win1 win1 win2 win2, -W window_size = 99'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/win1.fa -2 t/win1.fa -3 t/win2.fa -4 t/win2.fa -W 99'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR=$'0 99 59 0 0 0 0 40 0 0\n99 160 21 0 0 0 0 40 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 28
echo 'Test 28: win tests win1 win1 win2 win2, -W window_size = 200'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/win1.fa -2 t/win1.fa -3 t/win2.fa -4 t/win2.fa -W 200'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR=$'0 160 80 0 0 0 0 80 0 0'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi

# Test 29
echo 'Test 28: three and four allelic sites.'
TOTAL=$[$TOTAL + 1]
TEST_CMD='./quad-aln-report -1 t/all_configs1.fa -2 t/all_configs2.fa -3 t/all_configs3.fa -4 t/all_configs4.fa'
echo "  $TEST_CMD"
TEST_STR=`${TEST_CMD}`
PASS_STR=$'0 14 1 1 1 1 1 1 1 1'
if [ "$TEST_STR" != "$PASS_STR" ];
then
    echo "  * FAIL: expected '$PASS_STR', got '$TEST_STR'"
    COUNTER=$[$COUNTER + 1]

else
    echo "  * PASS"
fi



# Final Results
N_PASS=$[ $TOTAL - $COUNTER ]
echo 
echo "quad-aln-report tests complete:"
echo "$N_PASS out of $TOTAL passed"
exit $COUNTER
