function setup_file {
    >&2 echo "Making test binaries"
    make -j3 bin/test/test_state_io bin/test/test_cvec3 bin/test/test_cvec3_half \
        bin/test/test_maths_core bin/test/test_expression \
        bin/test/test_floor_nn
}

@test "test_state_io" {
    bin/test/test_state_io
}

@test "test_cvec3" {
    bin/test/test_cvec3
}

@test "test_cvec3_halg" {
    bin/test/test_cvec3_half
}

@test "test_floor_nn" {
    bin/test/test_floor_nn
}

@test "test_hash" {
    run make bin/test/test_hash
    if echo "$output" | grep "fatal error: testu01/bbattery.h" ; then
        skip;
    else
        bin/test/test_hash
    fi
}

@test "test_expression" {
    bin/test/test_expression
}