#!/usr/bin/env bats

BINARIES=bin/create_state_stationary_water
BINARIES+=bin/create_state_moving_water
BINARIES+=bin/create_state_dimers
BINARIES+=bin/test_dump_v5_graph_instance
BINARIES+=bin/extract_state_from_orch_log

function setup_file {
    >&2 echo "Making binaries"
    make -j3 ${BINARIES}
}

@test "dimersSmall" {
    >&2 echo "Test= ${BATS_TEST_NAME}"
    local WORKING
    WORKING=.test/$BATS_TEST_NAME
    mkdir -p ${WORKING}
    >&2 echo "Creating state"
    bin/create_state_dimers 4 4 4 3 0.001 > ${WORKING}/in.state
    scripts/test_state_under_xml_orch.sh ${WORKING}/in.state 800 ${WORKING}
}

@test "dimersMedium" {
    >&2 echo "Test= ${BATS_TEST_NAME}"
    local WORKING
    WORKING=.test/$BATS_TEST_NAME
    mkdir -p ${WORKING}
    >&2 echo "Creating state"
    bin/create_state_dimers 5 6 7 3 0.003 > ${WORKING}/in.state
    scripts/test_state_under_xml_orch.sh ${WORKING}/in.state 300 ${WORKING}
}

@test "dimersBiggles" {
    >&2 echo "Test= ${BATS_TEST_NAME}"
    local WORKING
    WORKING=.test/$BATS_TEST_NAME
    mkdir -p ${WORKING}
    >&2 echo "Creating state"
    bin/create_state_dimers 9 10 11 4 0.004 > ${WORKING}/in.state 2> ${WORKING}/create_state.log
    scripts/test_state_under_xml_orch.sh ${WORKING}/in.state 300 ${WORKING} 2> ${WORKING}/run_state.log
}

@test "movingWaterSmall" {
    >&2 echo "Test= ${BATS_TEST_NAME}"
    local WORKING
    WORKING=.test/$BATS_TEST_NAME
    mkdir -p ${WORKING}
    >&2 echo "Creating state"
    bin/create_state_moving_water 4 4 4 3 0.001 > ${WORKING}/in.state
    scripts/test_state_under_xml_orch.sh ${WORKING}/in.state 800 ${WORKING}
}


@test "movingWaterMedium" {
    >&2 echo "Test= ${BATS_TEST_NAME}"
    local WORKING
    WORKING=.test/$BATS_TEST_NAME
    mkdir -p ${WORKING}
    >&2 echo "Creating state"
    bin/create_state_moving_water 6 7 8  3 0.004 > ${WORKING}/in.state
    scripts/test_state_under_xml_orch.sh ${WORKING}/in.state 300 ${WORKING}
}

@test "movingWaterBiggles" {
    >&2 echo "Test= ${BATS_TEST_NAME}"
    local WORKING
    WORKING=.test/$BATS_TEST_NAME
    mkdir -p ${WORKING}
    >&2 echo "Creating state"
    bin/create_state_moving_water 11 12 13 3 0.004 > ${WORKING}/in.state
    scripts/test_state_under_xml_orch.sh ${WORKING}/in.state 300 ${WORKING}
}

@test "staticWaterSmall" {
    >&2 echo "Test= ${BATS_TEST_NAME}"
    local WORKING
    WORKING=.test/$BATS_TEST_NAME
    mkdir -p ${WORKING}
    >&2 echo "Creating state"
    bin/create_state_stationary_water 4 4 4 3 0.001 > ${WORKING}/in.state
    scripts/test_state_under_xml_orch.sh ${WORKING}/in.state 800 ${WORKING}
}

@test "staticWaterMedium" {
    >&2 echo "Test= ${BATS_TEST_NAME}"
    local WORKING
    WORKING=.test/$BATS_TEST_NAME
    mkdir -p ${WORKING}
    >&2 echo "Creating state"
    bin/create_state_stationary_water 7 8 9 3 0.004 > ${WORKING}/in.state
    scripts/test_state_under_xml_orch.sh ${WORKING}/in.state 300 ${WORKING}
}

@test "staticWaterBiggles" {
    >&2 echo "Test= ${BATS_TEST_NAME}"
    local WORKING
    WORKING=.test/$BATS_TEST_NAME
    mkdir -p ${WORKING}
    >&2 echo "Creating state"
    bin/create_state_stationary_water 11 12 13 3 0.004 > ${WORKING}/in.state
    scripts/test_state_under_xml_orch.sh ${WORKING}/in.state 300 ${WORKING}
}