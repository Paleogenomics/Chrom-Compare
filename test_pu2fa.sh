dfs=`./pu2fa -c Scaffold_1 < t/Scaffold_1.pu | diff t/Scaffold_1.consensus.fa -`

if [ -z "$dfs" ]; then
    echo "pu2fa unit test passed!"
    exit 0
else
    echo "pu2fa unit test failed!"
fi
