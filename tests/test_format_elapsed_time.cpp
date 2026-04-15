#define TESTING_MAIN
#ifndef __BINARY_NAME__
#define __BINARY_NAME__ "test_format_elapsed_time"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Unit tests for FormatElapsedTime helper (issue #351)"
#endif

#include "testing.hpp"

#include "../dynadjust/include/functions/dnatimer.hpp"

using dynadjust::FormatElapsedTime;

TEST_CASE("Sub-second durations report milliseconds", "[time][issue351]") {
    REQUIRE(FormatElapsedTime(0.000078) == "0.078ms");
    REQUIRE(FormatElapsedTime(0.0) == "0.000ms");
    REQUIRE(FormatElapsedTime(0.5) == "500.000ms");
    REQUIRE(FormatElapsedTime(0.999) == "999.000ms");
}

TEST_CASE("Durations between 1 and 60 seconds report seconds", "[time][issue351]") {
    REQUIRE(FormatElapsedTime(1.0) == "1.000s");
    REQUIRE(FormatElapsedTime(12.580) == "12.580s");
    REQUIRE(FormatElapsedTime(59.999) == "59.999s");
}

TEST_CASE("Durations of 60 seconds or more report hh:mm:ss", "[time][issue351]") {
    REQUIRE(FormatElapsedTime(60.0) == "00:01:00");
    REQUIRE(FormatElapsedTime(68.0) == "00:01:08");
    REQUIRE(FormatElapsedTime(139.419) == "00:02:19");
    REQUIRE(FormatElapsedTime(3600.0) == "01:00:00");
    REQUIRE(FormatElapsedTime(5841.0) == "01:37:21");
    REQUIRE(FormatElapsedTime(26697.777) == "07:24:57");
    REQUIRE(FormatElapsedTime(36000.0) == "10:00:00");
}
