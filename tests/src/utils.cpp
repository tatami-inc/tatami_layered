#include <gtest/gtest.h>
#include "tatami_test/tatami_test.hpp"
#include "tatami_layered/utils.hpp"

TEST(Utils, AtLeastOne) {
    EXPECT_EQ(tatami_layered::atleastone(10), 10);
    EXPECT_EQ(tatami_layered::atleastone(1), 1);
    EXPECT_EQ(tatami_layered::atleastone(0), 1);
}

TEST(Utils, Categorize) {
    EXPECT_EQ(tatami_layered::categorize(100), tatami_layered::Category::U8);
    EXPECT_EQ(tatami_layered::categorize(1000), tatami_layered::Category::U16);
    EXPECT_EQ(tatami_layered::categorize(100000), tatami_layered::Category::U32);

    tatami_test::throws_error([]() -> void {
        tatami_layered::categorize(-1);
    }, "negative");

    tatami_test::throws_error([]() -> void {
        tatami_layered::categorize(10000000000ull);
    }, "outside of the range");

    EXPECT_EQ(tatami_layered::categorize(10.0), tatami_layered::Category::U8);

    tatami_test::throws_error([]() -> void {
        tatami_layered::categorize(1e10);
    }, "outside of the range");
}

TEST(Utils, CheckChunkSize) {
    {
        auto out = tatami_layered::check_chunk_size<int, std::uint8_t>(10);
        EXPECT_EQ(out, 10);
    }

    {
        auto out = tatami_layered::check_chunk_size<int, std::uint8_t>(1000);
        EXPECT_EQ(out, 256);
    }

    {
        auto out = tatami_layered::check_chunk_size<int, std::uint8_t>(256);
        EXPECT_EQ(out, 256);
    }

    {
        auto out = tatami_layered::check_chunk_size<int, std::uint16_t>(1000);
        EXPECT_EQ(out, 1000);
    }

    {
        auto out = tatami_layered::check_chunk_size<int, std::uint16_t>(100000);
        EXPECT_EQ(out, 65536);
    }

    {
        auto out = tatami_layered::check_chunk_size<int, std::uint16_t>(65536);
        EXPECT_EQ(out, 65536);
    }

    tatami_test::throws_error([]() -> void {
        tatami_layered::check_chunk_size<int, std::uint16_t>(0);
    }, "should be positive");
}
